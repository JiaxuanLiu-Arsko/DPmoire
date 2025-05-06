import numpy as np
from ase import Atoms
from ase.neighborlist import NeighborList
from ase.io.vasp import read_vasp, write_vasp
import torch
import re
from typing import Union, Optional, Callable, Dict
import warnings
import os

import ase.data
from ase.calculators.calculator import Calculator, all_changes
from ase.stress import full_3x3_to_voigt_6_stress

from nequip.data import AtomicData, AtomicDataDict
from nequip.data.transforms import TypeMapper
import nequip.scripts.deploy
from nequip.train.trainer import Trainer

def split_lattice(atoms, nx, ny, nz, rcut):
    # 将原子位置转换为分数坐标
    cell = atoms.get_cell()
    positions = atoms.get_positions()
    fractional_positions = atoms.get_scaled_positions() % 1.0  # 分数坐标，取模以考虑周期性

    # 获取晶格矢量的范数
    a_norms = np.linalg.norm(cell, axis=1)
    delta_f = rcut / a_norms

    # 创建一个字典来存储每个block的原子和ghost atoms
    blocks = {}

    # 计算每个原子所属的block
    atom_indices = np.arange(len(atoms))
    # 对于每个block，处理原子和ghost atoms
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                # 定义block区域的分数坐标范围
                fx_min = i / nx
                fx_max = (i + 1) / nx
                fy_min = j / ny
                fy_max = (j + 1) / ny
                fz_min = k / nz
                fz_max = (k + 1) / nz

                # 定义扩展区域（包括buffer layer）
                fx_min_ext = fx_min - delta_f[0]
                fx_max_ext = fx_max + delta_f[0]
                fy_min_ext = fy_min - delta_f[1]
                fy_max_ext = fy_max + delta_f[1]
                fz_min_ext = fz_min - delta_f[2]
                fz_max_ext = fz_max + delta_f[2]

                # 考虑周期性边界条件
                # 创建一个布尔索引，选择在扩展区域内的原子
                in_ext_region = (
                    ((fractional_positions[:, 0] >= fx_min_ext) & (fractional_positions[:, 0] < fx_max_ext)) |
                    ((fractional_positions[:, 0] + 1 >= fx_min_ext) & (fractional_positions[:, 0] + 1 < fx_max_ext)) |
                    ((fractional_positions[:, 0] - 1 >= fx_min_ext) & (fractional_positions[:, 0] - 1 < fx_max_ext))
                ) & (
                    ((fractional_positions[:, 1] >= fy_min_ext) & (fractional_positions[:, 1] < fy_max_ext)) |
                    ((fractional_positions[:, 1] + 1 >= fy_min_ext) & (fractional_positions[:, 1] + 1 < fy_max_ext)) |
                    ((fractional_positions[:, 1] - 1 >= fy_min_ext) & (fractional_positions[:, 1] - 1 < fy_max_ext))
                ) & (
                    ((fractional_positions[:, 2] >= fz_min_ext) & (fractional_positions[:, 2] < fz_max_ext)) |
                    ((fractional_positions[:, 2] + 1 >= fz_min_ext) & (fractional_positions[:, 2] + 1 < fz_max_ext)) |
                    ((fractional_positions[:, 2] - 1 >= fz_min_ext) & (fractional_positions[:, 2] - 1 < fz_max_ext))
                )

                # 获取扩展区域内的原子索引
                ext_indices = atom_indices[in_ext_region]
                ext_indices_in_ext = np.arange(len(ext_indices))
                # 获取block区域内的原子索引
                in_block = (
                    ((fractional_positions[ext_indices, 0] >= fx_min) & (fractional_positions[ext_indices, 0] < fx_max)) |
                    ((fractional_positions[ext_indices, 0] + 1 >= fx_min) & (fractional_positions[ext_indices, 0] + 1 < fx_max)) |
                    ((fractional_positions[ext_indices, 0] - 1 >= fx_min) & (fractional_positions[ext_indices, 0] - 1 < fx_max))
                ) & (
                    ((fractional_positions[ext_indices, 1] >= fy_min) & (fractional_positions[ext_indices, 1] < fy_max)) |
                    ((fractional_positions[ext_indices, 1] + 1 >= fy_min) & (fractional_positions[ext_indices, 1] + 1 < fy_max)) |
                    ((fractional_positions[ext_indices, 1] - 1 >= fy_min) & (fractional_positions[ext_indices, 1] - 1 < fy_max))
                ) & (
                    ((fractional_positions[ext_indices, 2] >= fz_min) & (fractional_positions[ext_indices, 2] < fz_max)) |
                    ((fractional_positions[ext_indices, 2] + 1 >= fz_min) & (fractional_positions[ext_indices, 2] + 1 < fz_max)) |
                    ((fractional_positions[ext_indices, 2] - 1 >= fz_min) & (fractional_positions[ext_indices, 2] - 1 < fz_max))
                )

                # 获取block区域内的原子和ghost atoms
                block_indices_in_ext = ext_indices_in_ext[in_block]

                # 存储结果
                block_key = (i, j, k)
                blocks[block_key] = (ext_indices, block_indices_in_ext)

    return blocks

def jacobian(x:torch.Tensor, y:torch.Tensor):
    jacobian = torch.zeros((y.shape[0], x.shape[0],y.shape[1], x.shape[1]))
    for i in range(y.shape[0]):
        for j in range(y.shape[1]):
            # 对 y[i, j] 求导
            gradients = torch.autograd.grad(
                outputs=y[i, j],
                inputs=x,
                create_graph=False,
                retain_graph=True,  # 保留计算图
                allow_unused=True
            )[0]
            jacobian[i, :, j, :] = gradients
    return jacobian.detach().cpu().numpy()

def generate_sc(unitcell:Atoms, sc):
    cell = unitcell.get_cell().array
    pos = unitcell.get_positions()
    atomic_numbers = unitcell.get_atomic_numbers()
    cell_sc = np.dot(np.diag(sc), cell)
    pos_sc = np.zeros((sc[0]*sc[1]*sc[2]*len(pos), 3))
    atomic_numbers_sc = np.zeros((sc[0]*sc[1]*sc[2]*len(atomic_numbers)))
    i_atoms = 0
    unitcell_idx = []
    for idx in range(len(unitcell)):
        unitcell_idx.append(i_atoms)
        for k in range(sc[2]):
            for j in range(sc[1]):
                for i in range(sc[0]):
                    pos_sc[i_atoms] = pos[idx] + np.dot([i, j, k], cell)
                    atomic_numbers_sc[i_atoms] = atomic_numbers[idx]
                    i_atoms += 1
    supercell = Atoms(numbers=atomic_numbers_sc, positions=pos_sc, pbc=[True, True, True], cell=cell_sc)
    return supercell, unitcell_idx
    


class MoirePhono():
    """NequIP ASE Calculator.

    .. warning::

        If you are running MD with custom species, please make sure to set the correct masses for ASE.

    """

    def __init__(
        self,
        model: torch.jit.ScriptModule,
        r_max: float,
        device: Union[str, torch.device],
        energy_units_to_eV: float = 1.0,
        length_units_to_A: float = 1.0,
        transform: Callable = lambda x: x,
        nx: int = 1,
        ny: int = 1, 
        nz: int = 1,
        **kwargs
    ):
        self.results = {}
        self.model = model
        assert isinstance(
            model, torch.nn.Module
        ), "To build a NequIPCalculator from a deployed model, use NequIPCalculator.from_deployed_model"
        self.r_max = r_max
        self.device = device
        self.energy_units_to_eV = energy_units_to_eV
        self.length_units_to_A = length_units_to_A
        self.transform = transform
        self.nx = nx
        self.ny = ny
        self.nz = nz

    @classmethod
    def from_training_dir(
        cls,
        train_dir,
        device: Union[str, torch.device] = "cpu",
        species_to_type_name: Optional[Dict[str, str]] = None,
        set_global_options: Union[str, bool] = "warn",
        nx: int = 1,
        ny: int = 1, 
        nz: int = 1,
        **kwargs,
    ):
        # load model
        model, config = Trainer.load_model_from_training_session(train_dir, device=device)
        r_max = float(config[nequip.scripts.deploy.R_MAX_KEY])

        # build typemapper
        type_names = config[nequip.scripts.deploy.TYPE_NAMES_KEY]
        if species_to_type_name is None:
            # Default to species names
            warnings.warn(
                "Trying to use chemical symbols as NequIP type names; this may not be correct for your model! To avoid this warning, please provide `species_to_type_name` explicitly."
            )
            species_to_type_name = {s: s for s in ase.data.chemical_symbols}
        type_name_to_index = {n: i for i, n in enumerate(type_names)}
        chemical_symbol_to_type = {
            sym: type_name_to_index[species_to_type_name[sym]]
            for sym in ase.data.chemical_symbols
            if sym in type_name_to_index
        }
        if len(chemical_symbol_to_type) != len(type_names):
            raise ValueError(
                "The default mapping of chemical symbols as type names didn't make sense; please provide an explicit mapping in `species_to_type_name`"
            )
        transform = TypeMapper(chemical_symbol_to_type=chemical_symbol_to_type)
        model.eval()
        for module in model.modules():
            if len(re.findall("StressOutput", module.__class__.__name__))>0:
                module.train()
        #self.model.eval()
        # build nequip calculator
        if "transform" in kwargs:
            raise TypeError("`transform` not allowed here")
        return cls(
            model=model, r_max=r_max, device=device, transform=transform, nx=nx, ny=ny, nz=nz, **kwargs
        )

    def calculate_force_constant(self, unitcell:Atoms, sc):
        """
        Calculate properties.

        :param atoms: ase.Atoms object
        :param properties: [str], properties to be computed, used by ASE internally
        :param system_changes: [str], system changes since last calculation, used by ASE internally
        :return:
        """
        supercell, unitcell_idx = generate_sc(unitcell, sc)
        unit_in_sc = np.ones(len(supercell), dtype=int)
        unit_in_sc *= -1
        for i, idx in enumerate(unitcell_idx):
            unit_in_sc[idx] = i
        blocks = split_lattice(supercell, self.nx, self.ny, self.nz, 2*self.r_max)
        result = np.zeros((len(unitcell), len(supercell), 3, 3))
        for key, items in blocks.items():
            ext_indices = items[0]
            block_indices_in_ext = items[1]
            required_idxs = []
            for idx in block_indices_in_ext:
                if unit_in_sc[ext_indices[idx]]>-1:
                    required_idxs.append(idx)
            # prepare data
            splitted_atoms = supercell[ext_indices]
            data = AtomicData.from_ase(atoms=splitted_atoms, r_max=self.r_max)
            for k in AtomicDataDict.ALL_ENERGY_KEYS:
                if k in data:
                    del data[k]
            data = self.transform(data)
            data = data.to(self.device)
            data = AtomicData.to_AtomicDataDict(data)

            data[AtomicDataDict.POSITIONS_KEY].requires_grad_(True)
            
            #data[AtomicDataDict.FORCE_KEY].requires_grad_(True)
            # predict + extract data
            out = self.model(data)
            force_constants = -jacobian(
                y=out[AtomicDataDict.FORCE_KEY][required_idxs],
                x=data[AtomicDataDict.POSITIONS_KEY],
            )
            # def forces(x):
            #     data[AtomicDataDict.POSITIONS_KEY] = x
            #     data[AtomicDataDict.POSITIONS_KEY].requires_grad_(True)
            #     out = self.model(data)
            #     return out[AtomicDataDict.FORCE_KEY]
            # force_constants = torch.autograd.functional.jacobian(forces, data[AtomicDataDict.POSITIONS_KEY], vectorize=False).permute(0, 2, 1, 3).detach().cpu().numpy()

            #if not AtomicDataDict.FORCE_KEY in out:
            #    raise NotImplementedError("the model don't have forces as output!!!")
            for i, idx in enumerate(required_idxs): 
                result[unit_in_sc[ext_indices[idx]], ext_indices] = (
                    self.energy_units_to_eV / self.length_units_to_A / self.length_units_to_A
                ) * force_constants[i]
            del out
            del data
            print(f"{key} force constant finished.")

        return result
    
    def write_force_constant(self, unitcell, sc, path):
        force_constants = self.calculate_force_constant(unitcell, sc)
        with open(path, "w") as outfile:
            outfile.write(f"{force_constants.shape[0]} {force_constants.shape[1]}\n")
            for i in range(force_constants.shape[0]):
                for j in range(force_constants.shape[1]):
                    outfile.write(f"{i*sc[0]*sc[1]*sc[2]+1} {j} \n")
                    for dim1 in range(3):
                        for dim2 in range(3):
                            outfile.write(f"{force_constants[i, j, dim1, dim2]}")
                            outfile.write(" ")
                        outfile.write("\n")
                