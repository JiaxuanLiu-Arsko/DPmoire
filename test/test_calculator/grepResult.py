import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from ase.io.vasp import read_vasp_out

wkdir = "."

def getAseForces(workdir):
    allForces = [[], [], []]
    with open(f"{workdir}/f_modified.txt","r") as fin:
        for lines in fin:
            words = lines.split()
            allForces[0].append(float(words[0]))
            allForces[1].append(float(words[1]))
            allForces[2].append(float(words[2]))
    return allForces

def getVaspForces(workdir, nTerm, cmap, norm):
    allForces = [[], [], []]
    with open(f"{workdir}/f_origin.txt","r") as fin:
        for lines in fin:
            words = lines.split()
            allForces[0].append(float(words[0]))
            allForces[1].append(float(words[1]))
            allForces[2].append(float(words[2]))
    return allForces

cmap = mpl.colormaps["viridis"]
norm = Normalize(vmin=0, vmax=5-1)
aseForces = getAseForces(wkdir)
dftForces = getVaspForces(wkdir, nTerm=5, cmap=cmap, norm=norm)
fe = 0
for dim in range(3):
    for i in range(len(dftForces[0])):
        fe += np.square(dftForces[dim][i]-aseForces[dim][i])
fe /= len(dftForces[0])
fe = np.sqrt(fe)
with open(f"{wkdir}/error.dat", "w") as f:
    f.write(f"rmse of forces = {fe}\n")
    f.write("fxDFT      fyDFT       fzDFT       fxMLFF      fyMLFF      fzMLFF  \n")
    for i in range(len(dftForces[0])):
        for dim in range(3):
            f.write(f"{dftForces[dim][i]} ")
        for dim in range(3):
            f.write(f"{aseForces[dim][i]} ")
        f.write("\n")

for dim in range(3):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.scatter(dftForces[dim], aseForces[dim], s=0.5)
    refLine = np.arange(np.min(dftForces[dim]), np.max(dftForces[dim]), 0.01)
    #fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax)
    plt.xlim(-0.2, 0.2)
    plt.ylim(-0.2, 0.2)
    plt.xlabel("DFT FORCE(ev/Å)")
    plt.ylabel("MACHINE LEARNED FORCE(ev/Å)")
    #plt.plot(refLine, refLine, color='black')
    ax.set_aspect(1)
    plt.savefig(f"{wkdir}/Ferror{dim}.svg")
    plt.clf()