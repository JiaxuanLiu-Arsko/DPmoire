from ..preprocess import Config
from .dft_handler import DFTHandler
from ..data import Dataset
from ase.io.vasp import read_vasp, write_vasp
import os

class ValidationHandler(DFTHandler):

    angles = None
    script_dir = None
    collect_freq = None
    val_dir = None
    def __init__(self, config:Config|dict, val_dir, angles, existing_job:list=None, collect_freq:int=1):
        super().__init__(script_name="val_script.sh", n_nodes=config["n_nodes"], existing_job=existing_job, auto_resub=False)
        self.script_dir = config["script_dir"]
        self.collect_freq = collect_freq
        self.angles = angles
        self.val_dir = val_dir

    def run_calculation(self):
        for angle in self.angles:
            os.system(f"cp {self.script_dir}/{self.script_name} {self.val_dir}/{angle}/")
            self.submit_job(work_dir=f"{self.val_dir}/{angle}")

    def postprocess(self):
        val_dataset = self.make_dataset()
        val_dataset.save_extxyz(f"{self.val_dir}/valid.extxyz")
        return val_dataset
                
    def make_dataset(self) -> Dataset:
        val_dataset = Dataset()
        for angle in self.angles:
            val_dataset.load_dataset_OUTCAR(infile_str=f"{self.val_dir}/{angle}/OUTCAR", freq=self.collect_freq)
        return val_dataset
                    


                
        

