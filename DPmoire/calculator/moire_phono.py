import torch
import re
from typing import Union, Optional, Callable, Dict
import warnings
import time
import yaml
from datetime import datetime

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh
from ase import Atoms
from ase.neighborlist import NeighborList
import ase.data
import ase.units as units
from ase.io.vasp import read_vasp

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
    cell_sc = np.dot(np.diag(2*sc-1), cell)
    pos_sc = np.zeros(((sc[0]*2-1)*(sc[1]*2-1)*(sc[2]*2-1)*len(pos), 3))
    atomic_numbers_sc = np.zeros(((sc[0]*2-1)*(sc[1]*2-1)*(sc[2]*2-1)*len(atomic_numbers)))
    i_atoms = 0
    unitcell_idx = []
    for k in range(2*sc[2]-1):
        for j in range(2*sc[1]-1):
            for i in range(2*sc[0]-1):
                for idx in range(len(unitcell)):
                    if i==sc[0]-1 and j==sc[1]-1 and k==sc[2]-1:
                        unitcell_idx.append(i_atoms)
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
        sc_to_uc = np.ones(len(supercell), dtype=int)
        sc_to_uc *= -1
        for i, idx in enumerate(unitcell_idx):
            sc_to_uc[idx] = i
        blocks = split_lattice(supercell, self.nx, self.ny, self.nz, self.r_max)
        result = {}
        result_init_flag = np.zeros((len(unitcell), len(supercell)), dtype=bool)
        for key, items in blocks.items():
            ext_indices = items[0]
            block_indices_in_ext = items[1]
            block_flag = np.zeros(len(ext_indices), dtype=bool)
            block_flag[block_indices_in_ext] = True
            required_idxs_in_ext = []
            for i, idx in enumerate(ext_indices):
                if sc_to_uc[idx]>-1:
                    required_idxs_in_ext.append(i)
            # prepare data
            splitted_atoms = supercell[ext_indices]
            data = AtomicData.from_ase(atoms=splitted_atoms, r_max=self.r_max)
            for k in AtomicDataDict.ALL_ENERGY_KEYS:
                if k in data:
                    del data[k]
            data = self.transform(data)
            data = data.to(self.device)
            data = AtomicData.to_AtomicDataDict(data)
            edge_indx_list = []
            for i, idx in enumerate(data[AtomicDataDict.EDGE_INDEX_KEY][0]):
                if block_flag[idx.item()]:
                    edge_indx_list.append(i)
            data[AtomicDataDict.POSITIONS_KEY].requires_grad_(True)
            data[AtomicDataDict.EDGE_CELL_SHIFT_KEY] = data[AtomicDataDict.EDGE_CELL_SHIFT_KEY][edge_indx_list]
            data[AtomicDataDict.EDGE_INDEX_KEY] = data[AtomicDataDict.EDGE_INDEX_KEY][:, edge_indx_list]
            #data[AtomicDataDict.FORCE_KEY].requires_grad_(True)
            # predict + extract data
            out = self.model(data)
            force_constants = -jacobian(
                y=out[AtomicDataDict.FORCE_KEY][required_idxs_in_ext],
                x=data[AtomicDataDict.POSITIONS_KEY],
            )
            for i, idx in enumerate(required_idxs_in_ext): 
                for j, idx2 in enumerate(ext_indices):
                    if not result_init_flag[sc_to_uc[ext_indices[idx]], idx2]:
                        result[(sc_to_uc[ext_indices[idx]], idx2)] = np.zeros((3, 3))
                        result_init_flag[sc_to_uc[ext_indices[idx]], idx2] = True
                    result[(sc_to_uc[ext_indices[idx]], idx2)] += (
                        self.energy_units_to_eV / self.length_units_to_A / self.length_units_to_A
                    ) * force_constants[i, j]
            del out
            del data
            print(f"{key} force constant finished.")

        return unitcell_idx, result
    
    def write_force_constant(self, unitcell, sc, path):
        unitcell_idx, force_constants = self.calculate_force_constant(unitcell, sc)
        with open(path, "w") as outfile:
            for key, items in force_constants.items():
                outfile.write(f"{key[0]} {key[1]} \n")
                for dim1 in range(3):
                    for dim2 in range(3):
                        outfile.write(f"{items[dim1, dim2]}")
                        outfile.write(" ")
                    outfile.write("\n")

def read_force_constants(filename, n_atoms_unit, n_atoms_super):
    """读取力常数矩阵文件"""
    start_time = time.time()
    print("开始读取力常数矩阵...")
    
    with open(filename, 'r') as f:
        
        # 初始化力常数矩阵
        force_constants = np.zeros((n_atoms_unit, n_atoms_super, 3, 3))
        
        # 读取力常数数据
        while True:
            # 读取原子对索引
            line = f.readline()
            if not line:
                break
            i, j = map(int, line.split())
            
            # 读取3x3力常数子矩阵
            for a in range(3):
                data = list(map(float, f.readline().split()))
                for b in range(3):
                    force_constants[i, j, a, b] = data[b]
    
    end_time = time.time()
    print(f"力常数矩阵读取完成，耗时: {end_time - start_time:.2f} 秒")
    return force_constants

def get_dynamical_matrix(force_constants, structure, q_point, supercell_matrix=np.array([2, 2, 1])):
    """计算给定q点的动力学矩阵"""
    start_time = time.time()
    masses = structure.get_masses()
    scaled_positions = structure.get_scaled_positions()
    n_atoms_unit = force_constants.shape[0]
    n_sc = np.prod(2*supercell_matrix-1)
    
    # 初始化动力学矩阵
    dyn_matrix = np.zeros((3 * n_atoms_unit, 3 * n_atoms_unit), dtype=complex)
    
    # 预计算质量因子矩阵
    mass_factors = 1.0 / np.sqrt(masses[:, None] * masses[None, :])  # shape: (n_atoms_unit, n_atoms_unit)
    
    # 生成所有超胞位移的网格
    k1, k2, k3 = np.meshgrid(
        np.arange(-supercell_matrix[0]+1, supercell_matrix[0]),
        np.arange(-supercell_matrix[1]+1, supercell_matrix[1]),
        np.arange(-supercell_matrix[2]+1, supercell_matrix[2]),
        indexing='ij'
    )
    k_vectors = np.stack([k1.flatten(), k2.flatten(), k3.flatten()], axis=1)  # shape: (n_sc, 3)
    
    # 计算所有相位因子
    lat_phases = np.exp(2j * np.pi * np.dot(k_vectors, q_point))  # shape: (n_sc,)
    
    # 预计算超胞原子索引
    sc_indices = ((k_vectors[:, 0]+supercell_matrix[0]-1) * n_atoms_unit + 
                 (k_vectors[:, 1]+supercell_matrix[1]-1) * (2*supercell_matrix[0]-1) * n_atoms_unit + 
                 (k_vectors[:, 2]+supercell_matrix[2]-1) * (2*supercell_matrix[0]-1) * (2*supercell_matrix[1]-1) * n_atoms_unit)
    
    # 对每个单胞原子对进行计算
    for i in range(n_atoms_unit):
        for j in range(n_atoms_unit):
            # 获取这对原子的所有超胞相互作用
            super_j_indices = j + sc_indices
            fc_ij = force_constants[i, super_j_indices]  # shape: (n_sc, 3, 3)
            
            # 应用质量因子和相位因子
            mass_factor = mass_factors[i, j]
            atom_phases = np.exp(2j * np.pi * np.dot(scaled_positions[j]-scaled_positions[i], q_point))
            # 计算动力学矩阵元素
            for a in range(3):
                for b in range(3):
                    # 正向项
                    dyn_matrix[3*i + a, 3*j + b] = np.sum(fc_ij[:, a, b] * mass_factor * lat_phases * atom_phases)
    
    end_time = time.time()
    print(f"q点 {q_point} 的动力学矩阵构建完成，耗时: {end_time - start_time:.2f} 秒")
    return dyn_matrix

def calculate_phonon_spectrum(force_constants, masses, q_points, supercell_matrix=np.array([2, 2, 1]), sparse_solver=True, n_lowest_bands=19):
    """计算声子谱"""
    start_time = time.time()
    print("开始计算声子谱...")
    
    n_atoms_unit = force_constants.shape[0]
    n_qpoints = len(q_points)
    if sparse_solver:
        frequencies = np.zeros((n_qpoints, n_lowest_bands))
        eigenvectors = np.zeros((n_qpoints, 3 * n_atoms_unit, n_lowest_bands), dtype=complex)
    else:
        frequencies = np.zeros((n_qpoints, 3 * n_atoms_unit))
        eigenvectors = np.zeros((n_qpoints, 3 * n_atoms_unit, 3 * n_atoms_unit), dtype=complex)
    
    for i, q in enumerate(q_points):
        # 获取动力学矩阵
        dyn_matrix = get_dynamical_matrix(force_constants, masses, q, supercell_matrix)
        
        diag_start_time = time.time()
        if sparse_solver:
            # 转换为稀疏矩阵格式
            sparse_dyn_matrix = csr_matrix(dyn_matrix)
            
            # 使用稀疏矩阵求解器
            eigenvalues, eigenvecs = eigsh(sparse_dyn_matrix, k=n_lowest_bands, 
                                         which='LM', sigma=1e-10)
            sorted_indices = np.argsort(eigenvalues)
            eigenvalues = eigenvalues[sorted_indices]
            eigenvecs = eigenvecs[:, sorted_indices]
        else:
            # 使用标准numpy求解器
            eigenvalues, eigenvecs = np.linalg.eigh(dyn_matrix)
        
        diag_end_time = time.time()
        print(f"q点 {q} 的对角化完成，耗时: {diag_end_time - diag_start_time:.2f} 秒")
        # 转换为频率（THz）
        frequencies[i] = np.sign(eigenvalues) * np.sqrt(np.abs(eigenvalues)) / (2 * np.pi * THz)
        eigenvectors[i] = eigenvecs
    
    end_time = time.time()
    print(f"声子谱计算完成，总耗时: {end_time - start_time:.2f} 秒")
    return frequencies, eigenvectors

def write_band_yaml(q_points, frequencies, structure, filename='band.yaml'):
    """将计算结果以band.yaml格式输出"""
    start_time = time.time()
    print("开始写入band.yaml文件...")
    
    # 准备数据
    nqpoint = len(q_points)
    natom = len(structure)
    nbands = frequencies.shape[1]
    reciprocal_lattice = structure.cell.reciprocal().array.tolist()
    data = {
        'nqpoint': nqpoint,
        'npath': 1,  # 假设只有一条路径
        'natom': natom,
        'reciprocal_lattice': reciprocal_lattice,
        'phonon': []
    }
    q_distance = np.linalg.norm(q_points[0] @ reciprocal_lattice)
    # 添加每个q点的数据
    for i, q in enumerate(q_points):
        q_data = {
            'q-position': q.tolist(),
            'distance': float(q_distance),
            'band': []
        }
        
        # 添加每个能带的频率
        for freq in frequencies[i]:
            q_data['band'].append({'frequency': float(freq)})
        
        data['phonon'].append(q_data)
    
    # 写入yaml文件
    with open(filename, 'w') as f:
        f.write('# Generated by phonon_analysis.py\n')
        f.write(f'# Date: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n')
        yaml.dump(data, f, default_flow_style=False)
    
    end_time = time.time()
    print(f"band.yaml文件写入完成，耗时: {end_time - start_time:.2f} 秒")
