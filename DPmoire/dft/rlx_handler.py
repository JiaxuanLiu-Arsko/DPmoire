from ..preprocess import Config
from .dft_handler import DFTHandler
from ..data import Dataset
from ase.io.vasp import read_vasp, write_vasp
import os
import numpy as np
class RelaxationHandler(DFTHandler):

    n_secs_layer2 = None
    n_secs_layer3 = None
    work_dir = None
    script_dir = None
    collect_freq = None
    def __init__(self, config:Config|dict, existing_job:list=None, stackings=None):
        super().__init__(script_name=config["DFT_script"], n_nodes=config["n_nodes"], existing_job=existing_job, auto_resub=config["auto_resub"])
        self.n_secs_layer2 = config["n_sectors_layer2"]
        self.n_secs_layer3 = config["n_sectors_layer3"]
        self.work_dir = config["work_dir"]
        self.script_dir = config["script_dir"]
        self.collect_freq = config["OUTCAR_collect_freq"]
        if stackings is not None:
            self.stackings = stackings
        else:
            combos = []
            seen = set()
            for i2 in range(self.n_secs_layer2):
                for j2 in range(self.n_secs_layer2):
                    vec = (i2, j2, 0, 0)
                    if vec not in seen:
                        combos.append(list(vec))
                        seen.add(vec)
            for i3 in range(self.n_secs_layer3):
                for j3 in range(self.n_secs_layer3):
                    vec = (0, 0, i3, j3)
                    if vec not in seen:
                        combos.append(list(vec))
                        seen.add(vec)
            self.stackings = np.array(combos)

    def run_calculation(self):
        for stck in self.stackings:
            i2, j2, i3, j3 = map(int, stck)
            target_dir = f"{self.work_dir}/{i2}_{j2}_{i3}_{j3}"
            os.system(f"cp {self.script_dir}/{self.script_name} {target_dir}/")
            self.submit_job(work_dir=target_dir)

    def postprocess(self):
        self.CONT_to_POS()
        rlx_dataset = self.make_dataset()
        rlx_dataset.save_extxyz(f"{self.work_dir}/rlx_data.extxyz")
        return rlx_dataset
                
    def make_dataset(self) -> Dataset:
        rlx_dataset = Dataset()
        for stck in self.stackings:
            i2, j2, i3, j3 = map(int, stck)
            rlx_dataset.load_dataset_OUTCAR(
                infile_str=f"{self.work_dir}/{i2}_{j2}_{i3}_{j3}/OUTCAR",
                freq=self.collect_freq,
            )
        return rlx_dataset
    
    def save_rlx_results(self) -> None:
        for stck in self.stackings:
            i2, j2, i3, j3 = map(int, stck)
            target_dir = f"{self.work_dir}/{i2}_{j2}_{i3}_{j3}"
            os.system(f"cp {target_dir}/OUTCAR {target_dir}/OUTCAR-relax")
            os.system(f"cp {target_dir}/OSZICAR {target_dir}/OSZICAR-relax")
            os.system(f"cp {target_dir}/XDATCAR {target_dir}/XDATCAR-relax")

    def CONT_to_POS(self):
        for stck in self.stackings:
            i2, j2, i3, j3 = map(int, stck)
            target_dir = f"{self.work_dir}/{i2}_{j2}_{i3}_{j3}"
            outfile = open(f"{target_dir}/POSCAR", "w")
            with open(f"{target_dir}/CONTCAR", "r") as infile:
                for lines in infile:
                    if len(lines.split())==0:
                        break
                    outfile.write(lines)
            outfile.close()
                    


                
        
