import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import LinearNDInterpolator, griddata
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from ase.io.vasp import *
from ase.build import sort, rotate, make_supercell
from ase import Atoms

def plot_distance_z(atoms:Atoms, sc:int, element:str, s=None, colormap='Spectral_r', interpolator=LinearNDInterpolator, vrange=None):

    rotate(atoms, atoms.get_cell()[0], [1, 0, 0], atoms.get_cell()[2], [0, 0, 1], rotate_cell=True)
    atoms = make_supercell(atoms, P=np.array([[sc, 0, 0], [0, sc, 0], [0, 0, 1]]))
    cell = atoms.get_cell().array
    shift_vec = cell[0]*0.5 + cell[1]*0.5
    pos = atoms.get_positions() - shift_vec
    avg_z = np.mean(pos[:, 2])
    symbols = atoms.get_chemical_symbols()

    top_idx = []
    bot_idx = []
    for i, symbol in enumerate(symbols):
        if symbol==element and pos[i, 2]>avg_z:
            top_idx.append(i)
        elif symbol==element and pos[i, 2]< avg_z:
            bot_idx.append(i)
    pos_top = pos[top_idx]
    pos_bot = pos[bot_idx]
    
    top_fn = interpolator((pos_top[:, 0], pos_top[:, 1]), pos_top[:, 2])
    bot_fn = interpolator((pos_bot[:, 0], pos_bot[:, 1]), pos_bot[:, 2])
    points_d = []
    drop_idx = []
    for i, point in enumerate(pos_top):
        z = top_fn((point[0], point[1])) - bot_fn((point[0], point[1]))
        points_d.append([point[0], point[1], z])
        if np.isnan(z):
            drop_idx.append(i)
    drop_idx.reverse()
    for idx in drop_idx:
        points_d.pop(idx)
    points_d=np.array(points_d)
    if vrange is None:
        vmax = np.max(points_d[:, 2])
        vmin = np.min(points_d[:, 2])
    else:
        vmax = vrange[1]
        vmin = vrange[0]
    print(vmax, vmin)
    cmap = mpl.colormaps[colormap]
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    if s is None:
        s = 150000/len(top_idx)
        #print("s: ", s)

    ax.scatter(points_d[:, 0], points_d[:, 1], color=cmap(norm(points_d[:, 2])), s=s)
    ax.set_aspect(1)
#    fcb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
#    fcb.set_label('Interlayer distance(Å)')
#    fcb.set_ticks(ticks=[vmin, vmax], labels=["{0:.2f}".format(vmin), "{0:.2f}".format(vmax)], fontsize=15)
    lim = np.linalg.norm(cell[0])/4
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    return ax

def plot_disp_in_plane(start_config:Atoms, end_config:Atoms, sc:int, element:str, scale=None, width=None, colormap='Spectral_r', vrange=None):

    rotate(end_config, end_config.get_cell()[0], [1, 0, 0], end_config.get_cell()[2], [0, 0, 1], rotate_cell=True)
    cell = end_config.get_cell().array
    shift_vec = cell[0]*0.5 + cell[1]*0.5
    end_pos = end_config.get_positions()
    
    rotate(start_config, start_config.get_cell()[0], [1, 0, 0], start_config.get_cell()[2], [0, 0, 1], rotate_cell=True)
    start_config.set_cell(cell, scale_atoms=True)
    shift_vec = cell[0]*0.5 + cell[1]*0.5
    start_pos = start_config.get_positions()

    avg_z = np.mean(start_pos[:, 2])
    top_idx = []
    symbols = end_config.get_chemical_symbols()
    for i, symbol in enumerate(symbols):
        if symbol==element and start_pos[i, 2]>avg_z:
            top_idx.append(i)

    intra_disp = []
    positions = []
    for i in top_idx:
        disp = end_pos[i] - start_pos[i]
        mod_vec = np.dot(disp, np.linalg.inv(cell))
        for dim in range(3):
            if mod_vec[dim] > 0.5:
                mod_vec[dim] -= 1
            elif mod_vec[dim] < -0.5:
                mod_vec[dim] += 1
        disp = np.dot(mod_vec, cell)
        disp[2] = np.sqrt(disp[0] * disp[0] + disp[1] * disp[1])
        for j in range(sc):
            for k in range(sc):
                intra_disp.append(disp)
                positions.append(start_pos[i] + j*cell[0] + k*cell[1] - shift_vec*sc)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    intra_disp_arr = np.array(intra_disp)
    positions_arr = np.array(positions)
    colors = []
    cmap = mpl.colormaps[colormap]
    if vrange is None:
        vmin = 0
        vmax = np.max(intra_disp_arr[:, 2])
    else:
        vmax = vrange[1]
        vmin = vrange[0]
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    print(vmax, vmin)
    for length in intra_disp_arr[:, 2]:
        colors.append(cmap(norm(length)))

    if scale is None:
        scale = vmax*np.sqrt(len(top_idx))*sc/3
        #print(scale)

    if width is None:
        width = 1/(np.sqrt(len(top_idx))*sc)/3.3 #0.006
        #print(width)

    ax.quiver(positions_arr[:, 0], positions_arr[:, 1],  intra_disp_arr[:, 0], intra_disp_arr[:, 1], 
            pivot="mid", color=colors, scale=scale, width=width)

    ax.set_aspect(1)
    #fcb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,)
    #fcb.set_label('Intralayer displacement (Å)')
    #fcb.set_ticks(ticks=[vmin, vmax], labels=[f"0.000", "{0:.3f}".format(vmax)], fontsize=15)

    lim = np.linalg.norm(cell[0])/4*sc
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    return ax