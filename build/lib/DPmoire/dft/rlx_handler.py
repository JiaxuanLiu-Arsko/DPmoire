from ..preprocess import Config
from .dft_handler import DFTHandler
from ..data import Dataset
from ase.io.vasp import read_vasp, write_vasp
import os

class RelaxationHandler(DFTHandler):

    n_secs = None
    work_dir = None
    script_dir = None
    collect_freq = None
    def __init__(self, config:Config|dict, existing_job:list=None):
        super().__init__(script_name=config["DFT_script"], n_nodes=config["n_nodes"], existing_job=existing_job)
        self.n_secs = config["n_sectors"]
        self.work_dir = config["work_dir"]
        self.script_dir = config["script_dir"]
        self.collect_freq = config["OUTCAR_collect_freq"]

    def run_calculation(self):
        for i in range(self.n_secs):
            for j in range(self.n_secs):
                os.system(f"cp {self.script_dir}/{self.script_name} {self.work_dir}/{i}_{j}/")
                self.submit_job(work_dir=f"{self.work_dir}/{i}_{j}")
        self.wait_until_finished()

    def postprocess(self):
        self.CONT_to_POS()
        rlx_dataset = self.make_dataset()
        rlx_dataset.save_extxyz(f"{self.work_dir}/rlx_data.extxyz")
        return rlx_dataset
                
    def make_dataset(self) -> Dataset:
        rlx_dataset = Dataset()
        for i in range(self.n_secs):
            for j in range(self.n_secs):
                rlx_dataset.load_dataset_OUTCAR(infile_str=f"{self.work_dir}/{i}_{j}/OUTCAR", freq=self.collect_freq)
        return rlx_dataset
    
    def CONT_to_POS(self):
        for i in range(self.n_secs):
            for j in range(self.n_secs):
                outfile = open(f"{self.work_dir}/{i}_{j}/POSCAR", "w")
                with open(f"{self.work_dir}/{i}_{j}/CONTCAR", "r") as infile:
                    for lines in infile:
                        if len(lines.split())==0:
                            break
                        outfile.write(lines)
                outfile.close()
                    


                
        

