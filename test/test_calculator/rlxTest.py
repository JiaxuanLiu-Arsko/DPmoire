from ase.optimize.fire import FIRE               #弛豫用的优化器，模拟退火
#from ase.filters import FrechetCellFilter
import DPmoire.calculator.allegro_calculator as calc
import nequip.ase.nequip_calculator as calc_n
from ase.io.vasp import write_vasp, read_vasp    #读取/写入POSCAR
from ase.io import Trajectory
import os
import numpy as np

tbg = read_vasp(f'./7.34.vasp')     
#设置机器学习势
tbg.calc = calc.AllegroLargeCellCalculator.from_deployed_model(f"./WSe2_0830.pth", device="cuda", nx=2, ny=2)

optimizer = FIRE(atoms=tbg, trajectory=f"./data.traj", logfile=f"./rlx.log")
optimizer.run(fmax=0.005)
write_vasp(f"./7.34_relaxed.vasp", tbg)    #写入弛豫后的结构

f_orig = []
f_mod = []

stress_orig = []
stress_mod = []

traj = Trajectory("./data.traj")
for atoms in traj:
    atoms.calc = calc.AllegroLargeCellCalculator.from_deployed_model(f"./WSe2_0830.pth", device="cuda", nx=2, ny=2)
    forces = atoms.get_forces()
    stress_mod.append(atoms.get_stress())
    for force in forces:
        f_mod.append(force)
    atoms.calc = calc_n.NequIPCalculator.from_deployed_model(f"./WSe2_0830.pth", device="cuda")
    forces = atoms.get_forces()
    for force in forces:
        f_orig.append(force)
    stress_orig.append(atoms.get_stress())

np.savetxt("stress_origin.txt", np.array(stress_orig))
np.savetxt("stress_modified.txt", np.array(stress_mod))
print(f"error of stress: {np.mean(np.linalg.norm(np.array(stress_orig)-np.array(stress_mod), axis=-1))}")

np.savetxt("f_origin.txt", np.array(f_orig))
np.savetxt("f_modified.txt", np.array(f_mod))