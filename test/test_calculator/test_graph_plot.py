from DPmoire.utils.graph_plot import plot_disp_in_plane, plot_distance_z
from ase.io.vasp import read_vasp
from matplotlib import pyplot as plt

end_config = read_vasp("3.89_relaxed.vasp")
ax = plot_distance_z(end_config, 6, "W")
ax.set_xticks([])
ax.set_yticks([])
plt.savefig("./z.pdf")