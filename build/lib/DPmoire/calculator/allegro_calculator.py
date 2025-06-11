import numpy as np
from typing import Union, Tuple, Dict, Optional, List, Set, Sequence
from ase import Atoms
from ase.neighborlist import NeighborList
from ase.io.vasp import read_vasp, write_vasp

from typing import Union, Optional, Callable, Dict
import warnings
import torch

import ase.data
from ase.calculators.calculator import Calculator, all_changes
from ase.stress import full_3x3_to_voigt_6_stress
from nequip.data.AtomicData import neighbor_list_and_relative_vec
from nequip.data import AtomicData, AtomicDataDict
from nequip.data.transforms import TypeMapper
import nequip.scripts.deploy
import nequip.ase.nequip_calculator

import os



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


class AllegroLargeCellCalculator(Calculator):
    """NequIP ASE Calculator.

    .. warning::

        If you are running MD with custom species, please make sure to set the correct masses for ASE.

    """

    implemented_properties = ["energy", "energies", "forces", "stress", "free_energy"]

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
        Calculator.__init__(self, **kwargs)
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
    def from_deployed_model(
        cls,
        model_path,
        device: Union[str, torch.device] = "cpu",
        species_to_type_name: Optional[Dict[str, str]] = None,
        set_global_options: Union[str, bool] = "warn",
        nx: int = 1,
        ny: int = 1, 
        nz: int = 1,
        **kwargs,
    ):
        # load model
        model, metadata = nequip.scripts.deploy.load_deployed_model(
            model_path=model_path,
            device=device,
            set_global_options=set_global_options,
        )
        r_max = float(metadata[nequip.scripts.deploy.R_MAX_KEY])

        # build typemapper
        type_names = metadata[nequip.scripts.deploy.TYPE_NAMES_KEY].split(" ")
        print(type_names)
        print(r_max)
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
        print(chemical_symbol_to_type)
        # build nequip calculator
        if "transform" in kwargs:
            raise TypeError("`transform` not allowed here")
        return cls(
            model=model, r_max=r_max, device=device, transform=transform, nx=nx, ny=ny, nz=nz, **kwargs
        )

    def calculate(self, atoms:ase.Atoms=None, properties=["energy"], system_changes=all_changes):
        """
        Calculate properties.

        :param atoms: ase.Atoms object
        :param properties: [str], properties to be computed, used by ASE internally
        :param system_changes: [str], system changes since last calculation, used by ASE internally
        :return:
        """
        # call to base-class to set atoms attribute
        Calculator.calculate(self, atoms)
        blocks = split_lattice(atoms, self.nx, self.ny, self.nz, self.r_max)
        self.results = {}
        self.results["energy"] = 0
        self.results["forces"] = np.zeros((len(atoms), 3))
        self.results["free_energy"] = 0
        self.results["energies"] = np.zeros((len(atoms)))
        #self.results["stress"] = np.zeros(6)
        stress = np.zeros((3, 3))
        for key, items in blocks.items():
            ext_indices = items[0]
            block_indices_in_ext = items[1]
            block_flag = np.zeros(len(ext_indices), dtype=bool)
            block_flag[block_indices_in_ext] = True
            # prepare data
            splitted_atoms = atoms[ext_indices]
            data = AtomicData.from_ase(atoms=splitted_atoms, r_max=self.r_max)
            for k in AtomicDataDict.ALL_ENERGY_KEYS:
                if k in data:
                    del data[k]
            data = self.transform(data)
            data = data.to(self.device)
            data = AtomicData.to_AtomicDataDict(data)
            #data = AtomicDataDict.with_edge_vectors(data)
            edge_indx_list = []
            for i, idx in enumerate(data[AtomicDataDict.EDGE_INDEX_KEY][0]):
                if block_flag[idx.item()]:
                    edge_indx_list.append(i)
            #data[AtomicDataDict.POSITIONS_KEY] = data[AtomicDataDict.POSITIONS_KEY][block_indices_in_ext]
            #data[AtomicDataDict.EDGE_LENGTH_KEY] = data[AtomicDataDict.EDGE_LENGTH_KEY][edge_indx_list]
            #data[AtomicDataDict.EDGE_VECTORS_KEY] = data[AtomicDataDict.EDGE_VECTORS_KEY][edge_indx_list]
            #data[AtomicDataDict.ATOM_TYPE_KEY] = data[AtomicDataDict.ATOM_TYPE_KEY][block_indices_in_ext]
            #data[AtomicDataDict.ATOMIC_NUMBERS_KEY] = data[AtomicDataDict.ATOMIC_NUMBERS_KEY][block_indices_in_ext]
            data[AtomicDataDict.EDGE_CELL_SHIFT_KEY] = data[AtomicDataDict.EDGE_CELL_SHIFT_KEY][edge_indx_list]
            data[AtomicDataDict.EDGE_INDEX_KEY] = data[AtomicDataDict.EDGE_INDEX_KEY][:, edge_indx_list]
            # predict + extract data
            out = self.model(data)
            # only store results the model actually computed to avoid KeyErrors
            for i, idx in enumerate(block_indices_in_ext):
                if AtomicDataDict.PER_ATOM_ENERGY_KEY in out:
                    self.results["energies"][ext_indices[idx]] = self.energy_units_to_eV * (
                        out[AtomicDataDict.PER_ATOM_ENERGY_KEY]
                        .detach()
                        .squeeze(-1)
                        .cpu()
                        .numpy()
                    )[idx]

            if AtomicDataDict.FORCE_KEY in out:
                # force has units eng / len:
                self.results["forces"][ext_indices] += (
                    self.energy_units_to_eV / self.length_units_to_A
                ) * out[AtomicDataDict.FORCE_KEY].detach().cpu().numpy()

            if AtomicDataDict.STRESS_KEY in out:
                partial_stress = out[AtomicDataDict.STRESS_KEY].detach().cpu().numpy()
                stress += partial_stress.reshape(3, 3) * (
                    self.energy_units_to_eV / self.length_units_to_A**3
                )

        if AtomicDataDict.STRESS_KEY in out:
            # ase wants voigt format
            stress_voigt = full_3x3_to_voigt_6_stress(stress)
            self.results["stress"] = stress_voigt

        if AtomicDataDict.TOTAL_ENERGY_KEY in out:
            self.results["energy"] = np.sum(self.results["energies"])
            # "force consistant" energy
            self.results["free_energy"] = self.results["energy"]
        
