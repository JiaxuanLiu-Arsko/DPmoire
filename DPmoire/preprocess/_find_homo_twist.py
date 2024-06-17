import numpy as np
from ase import Atoms

def search_twist(N_min:int=2, N_max:int=10):
    mat_list = []
    angle_list = []
    for i in range(N_min, N_max+1):
        mat_1 = np.array([[i, i+1, 0], [-(i+1), 2*i+1, 0], [0, 0, 1]]).T
        mat_2 = np.array([[i+1, i, 0], [-i, 2*i+1, 0], [0, 0, 1]]).T
        mat_list.append((mat_1, mat_2))
        angle_list.append("{0:.5f}".format(np.arccos((3*i**2+3*i+0.5)/(3*i**2+3*i+1))/np.pi*180))
    return angle_list, mat_list

def adjust_atoms_d(top_atoms:Atoms, bot_atoms:Atoms, d:float):
    cell = np.dot(top_atoms.get_cell().array, np.array([[1, 0, 0], [0, 1, 0], [0, 0, 0.5]]))
    top_atoms.set_cell(cell)
    bot_atoms.set_cell(cell)
    top_z = np.mean(top_atoms.get_scaled_positions()[:, 2])
    bot_z = np.mean(bot_atoms.get_scaled_positions()[:, 2])
    top_shift_c = d/top_atoms.get_cell().lengths()[2]/2-top_z
    bot_shift_c = 1-d/bot_atoms.get_cell().lengths()[2]/2-bot_z
    top_spos = top_atoms.get_scaled_positions()
    top_spos[:, 2] += top_shift_c
    top_atoms.set_scaled_positions(top_spos)
    bot_spos = bot_atoms.get_scaled_positions()
    bot_spos[:, 2] += bot_shift_c
    bot_atoms.set_scaled_positions(bot_spos)
    return top_atoms, bot_atoms