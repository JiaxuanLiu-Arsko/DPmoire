from ase.optimize.fire import FIRE               #弛豫用的优化器，模拟退火
from ase.optimize.fire2 import FIRE2
from ase.filters import FrechetCellFilter
import DPmoire.calculator.allegro_calculator as calc
import nequip.ase.nequip_calculator as calc_n
from ase.io.vasp import write_vasp, read_vasp    #读取/写入POSCAR
from ase.io import Trajectory
import os
import numpy as np

tbg = read_vasp(f'./0.88211.vasp')     
#设置机器学习势
tbg.calc = calc.AllegroLargeCellCalculator.from_deployed_model(f"./WSe2_ivdw11_0528.pth", device="cuda", nx=int(np.sqrt(len(tbg)/1300)), ny=int(np.sqrt(len(tbg)/1300)))

optimizer = FIRE(atoms=tbg, trajectory=f"./data.traj", logfile=f"./rlx.log")
optimizer.run(fmax=0.005)
write_vasp(f"./0.88211.vasp_relaxed.vasp", tbg)    #写入弛豫后的结构