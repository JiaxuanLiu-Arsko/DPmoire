import os, time, re, torch
import numpy as np
import multiprocessing
from ase.io.vasp import read_vasp, write_vasp
from nequip.ase import nequip_calculator as calc
from ase.md.langevin import Langevin
from ase import Atoms
from ase import units
from ..preprocess import Config, EnvironmentHandler
import matplotlib.pyplot as plt

def findError(config:Config, nTerm, termSize, mode):
    if mode=="run":
        aseForces = runAse(config, nTerm, termSize)
        dftForces = runVasp(config, nTerm)
    elif mode=="get":
        aseForces = getAseForces(params=config)
        dftForces = getVaspForces(params=config, nTerm=nTerm)
    fe = 0
    for dim in range(3):
        for i in range(len(dftForces[0])):
            fe += np.square(dftForces[dim][i]-aseForces[dim][i])
    fe /= len(dftForces[0])
    fe = np.sqrt(fe)
    with open(f"{config['test_dir']}/error.dat", "w") as f:
        f.write(f"rmse of forces = {fe}\n")
        f.write("fxDFT      fyDFT       fzDFT       fxMLFF      fyMLFF      fzMLFF  \n")
        for i in range(len(dftForces[0])):
            for dim in range(3):
                f.write(f"{dftForces[dim][i]} ")
            for dim in range(3):
                f.write(f"{aseForces[dim][i]} ")
            f.write("\n")

    for dim in range(3):
        plt.scatter(dftForces[dim], aseForces[dim])
        refLine = np.arange(np.min(dftForces[dim]), np.max(dftForces[dim]), 0.1)
        plt.xlim(-1,1)
        plt.ylim(-1,1)
        plt.xlabel("DFT FORCE(ev/Å)")
        plt.ylabel("MACHINE LEARNED FORCE(ev/Å)")
        plt.plot(refLine, refLine)
        plt.savefig(f"{config['test_dir']}/Ferror{dim}.jpg")
        plt.clf()

def runAse(config:Config, nTerm, termSize):
    allForces = [[],[],[]]
    tbg = read_vasp(config['test_dir'] + "/POSCAR")
    tbg.calc = calc.nequip_calculator(config['work_dir'] + "/main/mlff.pth")
    dyn = Langevin(tbg, timestep=units.fs, temperature_K=40, friction=0.01)
    for i in range(nTerm):
        dyn.run(termSize)
        forces = tbg.get_forces()
        for force in forces:
            allForces[0].append(force[0])
            allForces[1].append(force[1])
            allForces[2].append(force[2])
        writePath = f"{config['test_dir']}/{i}"
        if not os.path.exists(writePath):
            os.makedirs(writePath)             #create write path if not exist
        write_vasp(f"{writePath}/POSCAR", tbg)
    with open(f"{config['test_dir']}/FMLFF.dat", "w") as f:
        for i in range(len(allForces[0])):
            for dim in range(3):
                f.write(f"{allForces[dim][i]} ")
            f.write("\n")
    return allForces

def getAseForces(config:Config):
    allForces = [[], [], []]
    with open(f"{config['test_dir']}/FMLFF.dat","r") as fin:
        for lines in fin:
            words = lines.split()
            allForces[0].append(float(words[0]))
            allForces[1].append(float(words[1]))
            allForces[2].append(float(words[2]))
    return allForces

def k_return_2(lat_vec):
    return [2, 2, 1]

def runVasp(config:Config, nTerm):
    env_handler = EnvironmentHandler(config, k_generator=k_return_2)
    for i in range(nTerm):
        curDir = f"{config['test_dir']}/{i}"
        os.system(f"cp {config['script_dir']}/* {curDir}")
        env_handler.gen_KPOINTS(outDir=curDir)
        env_handler.gen_POTCAR(elements=env_handler.elements, outDir=curDir)
        env_handler.gen_INCAR(input_file=f"{config['test_dir']}/INCAR", out_dir=curDir)
        os.chdir(curDir)
        np = multiprocessing.cpu_count()
        if torch.cuda.is_available():
            np = torch.cuda.device_count()
        returncode = os.system(f"mpirun -np {np} vasp_std")
        if returncode != 0:
            print("mpirun error!!")
    return getVaspForces(config, nTerm=nTerm)


def getVaspForces(config:Config, nTerm):
    allForces = [[],[],[]]
    for i in range(nTerm):
        curDir = f"{config['test_dir']}/{i}"
        with open(curDir + "/OUTCAR", 'r') as f:
            outcar_content = f.read()
        structures = re.findall(r'POSITION\s+TOTAL-FORCE.*?\n\s+(-{0,1}\d+\.\d+\s+-{0,1}\d+\.\d+\s+-{0,1}\d+\.\d+.*?)\n\s+(?=-{3,})', outcar_content, re.DOTALL)
        for structure in structures:
            lines = structure.split("\n")
            for line in lines:
                words = line.split()
                allForces[0].append(float(words[3]))
                allForces[1].append(float(words[4]))
                allForces[2].append(float(words[5]))
    with open(f"{config['test_dir']}/FDFT.dat", "w") as f:
        for i in range(len(allForces[0])):
            for dim in range(3):
                f.write(f"{allForces[dim][i]} ")
            f.write("\n")
    return allForces

