from DPmoire.utils.graph_plot import plot_disp_in_plane, plot_distance_z
from ase.io.vasp import read_vasp
from matplotlib import pyplot as plt

start_config = read_vasp("POSCAR/7.34.vasp")
end_config = read_vasp("POSCAR/7.34_relaxed.vasp")
ax = plot_disp_in_plane(start_config, end_config, 6, "W", )
ax.set_xticks([])
ax.set_yticks([])
plt.savefig("./WSe2_7.34_inplane.pdf")

end_config = read_vasp("POSCAR/7.34_relaxed.vasp")
ax = plot_distance_z(end_config, 6, "W")
ax.set_xticks([])
ax.set_yticks([])
plt.savefig("./WSe2_7.34_z.pdf")

end_config = read_vasp("POSCAR/3.89_relaxed.vasp")
ax = plot_distance_z(end_config, 6, "W")
ax.set_xticks([])
ax.set_yticks([])
plt.savefig("./WSe2_3.89_z.pdf")

start_config = read_vasp("POSCAR/3.89.vasp")
end_config = read_vasp("POSCAR/3.89_relaxed.vasp")
ax = plot_disp_in_plane(start_config, end_config, 6, "W", )
ax.set_xticks([])
ax.set_yticks([])
plt.savefig("./WSe2_3.89_inplane.pdf")