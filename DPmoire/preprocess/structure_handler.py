import numpy as np
import os,sys
from .config import Config
from ase.io.vasp import read_vasp, write_vasp
from ase.constraints import FixedLine 
from ase.build import make_supercell, sort, stack
from ase import Atoms
import copy
from ._find_homo_twist import search_twist, adjust_atoms_d

class StructureHandler:

    '''
    Handle crystal structures during the process.
    '''

    in_dir = None  
    out_dir = None      #Directory to store shifted structures.
    n_secs = None          #number of sectors to shift
    top_atoms = None
    bot_atoms = None
    top_indexes = None
    bot_indexes = None
    new_struct = None
    d = None

    def __init__(self, config:Config|dict=None):
        if isinstance(config, (Config, dict)):
            self.n_secs = config["n_sectors"]
            self.in_dir = config["input_dir"]
            self.d = config["d"]
        else:
            raise Exception(f"Unknown type of Conifg:{type(config)}")
        self.read_all_layers(self.in_dir)
        self.new_struct, self.top_indexes, self.bot_indexes = self.build_new_struct(d=self.d)

    def read_atoms(self, in_file:str):
        atoms = read_vasp(in_file)
        return atoms
    
    def read_all_layers(self, in_dir:str):
        self.top_atoms = self.read_atoms(f"{in_dir}/top_layer.poscar")
        self.bot_atoms = self.read_atoms(f"{in_dir}/bot_layer.poscar")

    def find_layer_idx(self, atoms:Atoms):
        cell_mat = atoms.get_cell().array#.transpose()
        frac_mat = np.linalg.inv(cell_mat)
        frac_pos = np.dot(atoms.get_positions(), frac_mat)
        top_idx = []
        bot_idx = []
        for i, pos in enumerate(frac_pos):
            if pos[2]>0.5:
                top_idx.append(i)
            else:
                bot_idx.append(i)
        return top_idx, bot_idx

    
    def build_new_struct(self, d:float):
        '''
        combine new structures combining top/bot layer atoms
        '''
        top_cell_mat = self.top_atoms.get_cell().array#.transpose()
        top_cell_len = self.top_atoms.get_cell().lengths()
        bot_cell_mat = self.bot_atoms.get_cell().array#.transpose()
        bot_cell_len = self.bot_atoms.get_cell().lengths()
        fractional_pos_top = np.dot(self.top_atoms.get_positions(), np.linalg.inv(top_cell_mat))
        fractional_pos_bot = np.dot(self.bot_atoms.get_positions(), np.linalg.inv(bot_cell_mat))
        new_cell_mat = np.array([top_cell_mat[k]*(1 + bot_cell_len[k]/top_cell_len[k])/2 for k in range(3)])
        c_top = np.mean([k[2] for k in fractional_pos_top])
        c_bot = np.mean([k[2] for k in fractional_pos_bot])
        for i, _ in enumerate(fractional_pos_top):
            fractional_pos_top[i][2] += 0.5 - c_top + d/(self.top_atoms.get_cell().lengths()[2]+self.bot_atoms.get_cell().lengths()[2])
        for i, _ in enumerate(fractional_pos_bot):
            fractional_pos_bot[i][2] += 0.5 - c_bot - d/(self.top_atoms.get_cell().lengths()[2]+self.bot_atoms.get_cell().lengths()[2])
        new_pos_top = [[item[0], item[1], item[2]] for item in np.dot(fractional_pos_top, new_cell_mat)]
        new_pos_bot = [[item[0], item[1], item[2]] for item in np.dot(fractional_pos_bot, new_cell_mat)]
        new_pos = []
        new_pos.extend(new_pos_top)
        new_pos.extend(new_pos_bot)
        new_symbols = []
        new_symbols.extend(self.top_atoms.get_chemical_symbols())
        new_symbols.extend(self.bot_atoms.get_chemical_symbols())
        atoms = sort(Atoms(positions=new_pos, symbols=new_symbols, cell=new_cell_mat#.transpose()
                      , pbc = [True, True, True]))
        top_idx, bot_idx = self.find_layer_idx(atoms)
        return atoms, top_idx, bot_idx
    
    def shift_atoms(self, i:int, j:int, c_constrain:bool=True, sc:int=2):

        '''
        Returning atoms with top layer atoms shifted by (i/n_secs*lat_vec[0] + j/n_secs*lat_vec[1]).
        Set c_constrain=True to add constrain in c direction.
        '''

        atoms = copy.deepcopy(self.new_struct)
        delta = i/self.n_secs * atoms.get_cell().array[0] + j/self.n_secs * atoms.get_cell().array[1]
        pos = atoms.get_positions()
        for idx in self.top_indexes:
            pos[idx] += delta
        atoms.set_positions(pos)
        atoms_sc = sort(make_supercell(prim=atoms, P=[[sc, 0, 0], [0, sc, 0], [0, 0, 1]]))
        top_idx, bot_idx = self.find_layer_idx(atoms_sc)
        if c_constrain:
            cons = FixedLine([top_idx[0], bot_idx[0]], direction=atoms_sc.cell.array[2]/atoms_sc.cell.lengths()[2])
            atoms_sc.set_constraint(cons)
        return atoms_sc
    
    def shift(self, i:int, j:int, out_dir:str, c_constrain:bool=True, sc:int = 2):

        """
        Write POSCAR at out_dir/i_j/ of shifted structures.
        Shift vector equals (i/n_secs*lat_vec[0] + j/n_secs*lat_vec[1]).
        Set c_constrain=True to add constrain in c direction.
        """
        if not os.path.exists(f"{out_dir}/{i}_{j}/"):
            os.mkdir(f"{out_dir}/{i}_{j}/")
        atoms_sc = self.shift_atoms(i, j, c_constrain, sc)
        write_vasp(f"{out_dir}/{i}_{j}/POSCAR", atoms=atoms_sc)

    def shift_all(self, out_dir:str, c_constrain:bool=True, sc:int = 2):
        for i in range(self.n_secs):
            for j in range(self.n_secs):
                self.shift(i, j, out_dir, c_constrain=c_constrain, sc=sc)

    def make_twist_struct(self, N_min, N_max, out_dir:str):
        top_atoms = copy.deepcopy(self.top_atoms)
        bot_atoms = copy.deepcopy(self.bot_atoms)
        top_atoms, bot_atoms = adjust_atoms_d(top_atoms, bot_atoms, self.d)
        angle_list, mat_list = search_twist(N_min, N_max)
        out_atoms_list = []
        for idx, mat in enumerate(mat_list):
            top_sc = make_supercell(top_atoms, P=mat[0])
            bot_sc = make_supercell(bot_atoms, P=mat[1])
            out_atoms = stack(bot_sc, top_sc, maxstrain=None, reorder=True)
            out_atoms_list.append(out_atoms)
            if not os.path.exists(f"{out_dir}/{angle_list[idx]}/"):
                os.mkdir(f"{out_dir}/{angle_list[idx]}/")
            write_vasp(f"{out_dir}/{angle_list[idx]}/POSCAR", out_atoms)
        return angle_list, out_atoms_list

