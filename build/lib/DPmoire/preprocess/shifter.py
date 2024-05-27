import numpy as np
import os,sys
from .parameters import Parameters
from ase.io.vasp import read_vasp, write_vasp
from ase.constraints import FixedLine 
from ase import Atoms
import copy

import numpy as np
import os,sys
from .parameters import Parameters

class Shifter:

    '''
    Shift the original POSCAR to desired stackings.
    '''

    in_dir = ""  
    out_dir = ""      #Directory to store shifted structures.
    lat_vec = []
    n_secs = 0           #number of sectors to shift
    delta_vec = []
    top_atoms = None
    bot_atoms = None
    top_indexes = []
    bot_indexes = []
    new_struct = None
    def update(self, n_secs:int=None, in_dir:str=None):

        '''
        Update parameters of shifter.
        '''

        if in_dir is not None:
            self.in_dir = in_dir
        if n_secs is not None :
            self.n_secs = n_secs

    def __init__(self, param:Parameters=None, n_secs:int=None, in_dir:str=None, d:float = 6.0):
        if param is not None:
            self.update(param.nSecs, param.originDir)
        elif (n_secs is not None) and (in_dir is not None):
            self.update(n_secs=n_secs, in_dir=in_dir)
        self.read_all_layers(self.in_dir)
        self.new_struct, self.top_indexes, self.bot_indexes = self.build_new_struct(d=d)

    def read_atoms(self, in_file:str):
        atoms = read_vasp(in_file)
        return atoms
    
    def read_all_layers(self, in_dir:str):
        self.top_atoms = self.read_atoms(f"{in_dir}/top_layer.poscar")
        self.bot_atoms = self.read_atoms(f"{in_dir}/bot_layer.poscar")

    
    def build_new_struct(self, d:float):
        '''
        combine new structures combining top/bot layer atoms
        '''
        top_cell_mat = self.top_atoms.get_cell().array.transpose()
        bot_cell_mat = self.bot_atoms.get_cell().array.transpose()
        fractional_pos_top = np.dot(self.top_atoms.get_positions(), np.linalg.inv(top_cell_mat))
        fractional_pos_bot = np.dot(self.bot_atoms.get_positions(), np.linalg.inv(bot_cell_mat))
        new_cell_mat = np.array([(top_cell_mat[k] + bot_cell_mat[k])/2 for k, _ in enumerate(top_cell_mat)])
        c_top = np.mean([k[2] for k in fractional_pos_top])
        c_bot = np.mean([k[2] for k in fractional_pos_bot])
        for i, _ in enumerate(fractional_pos_top):
            fractional_pos_top[i][2] += 0.5 - c_top + d/2/self.top_atoms.get_cell().lengths()[2]
        for i, _ in enumerate(fractional_pos_bot):
            fractional_pos_bot[i][2] += 0.5 - c_bot - d/2/self.bot_atoms.get_cell().lengths()[2]
        new_pos_top = [[item[0], item[1], item[2]] for item in np.dot(fractional_pos_top, new_cell_mat)]
        new_pos_bot = [[item[0], item[1], item[2]] for item in np.dot(fractional_pos_bot, new_cell_mat)]
        new_pos = []
        new_pos.extend(new_pos_top)
        new_pos.extend(new_pos_bot)
        new_symbols = []
        new_symbols.extend(self.top_atoms.get_chemical_symbols())
        new_symbols.extend(self.bot_atoms.get_chemical_symbols())
        atoms = Atoms(positions=new_pos, symbols=new_symbols, cell=new_cell_mat.transpose()
                      , pbc = [True, True, True])
        return atoms, range(len(new_pos_top)), range(len(new_pos_top), len(new_pos_bot)+len(new_pos_top))
    
    def shift_atoms(self, i:int, j:int, c_constrain:bool=True):

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
        if c_constrain:
            cons = FixedLine([self.top_indexes[0], self.bot_indexes[0]], direction=atoms.cell.array[2]/atoms.cell.lengths()[2])
            atoms.set_constraint(cons)
        return atoms
    
    def shift(self, i:int, j:int, out_dir:str, c_constrain:bool=True):

        """
        Write POSCAR at out_dir/i_j/ of shifted structures.
        Shift vector equals (i/n_secs*lat_vec[0] + j/n_secs*lat_vec[1]).
        Set c_constrain=True to add constrain in c direction.
        """
        if not os.path.exists(f"{out_dir}/{i}_{j}/"):
            os.mkdir(f"{out_dir}/{i}_{j}/")
        write_vasp(f"{out_dir}/{i}_{j}/POSCAR", atoms=self.shift_atoms(i, j, c_constrain))



