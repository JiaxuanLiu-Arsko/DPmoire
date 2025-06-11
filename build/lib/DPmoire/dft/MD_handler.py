import numpy as np
import os
import time
from ..preprocess import Config
from .dft_handler import DFTHandler
from ..data import Dataset

class MDHandler(DFTHandler):
    n_secs = None
    work_dir = None
    script_dir = None
    VASP_ML = None
    collect_freq = None
    input_dir = None
    stackings = None
    def __init__(self, config:Config|dict, existing_job:list=None, stackings = None):
        super().__init__(script_name=config["DFT_script"], n_nodes=config["n_nodes"], existing_job=existing_job, auto_resub=config["auto_resub"])
        self.n_secs = config["n_sectors"]
        self.work_dir = config["work_dir"]
        self.script_dir = config["script_dir"]
        self.VASP_ML = config["VASP_ML"]
        self.collect_freq = config["OUTCAR_collect_freq"]
        self.input_dir = config['input_dir']
        if stackings is not None:
            self.stackings = stackings
        else:
            self.stackings = np.array([[i, j] for i in range(self.n_secs) for j in range(self.n_secs)])

    def run_calculation(self):
        for stck in self.stackings:
            i = stck[0]
            j = stck[1]
            os.system(f"cp {self.script_dir}/{self.script_name} {self.work_dir}/{i}_{j}/")
            self.submit_job(work_dir=f"{self.work_dir}/{i}_{j}")
        os.system(f"cp {self.script_dir}/{self.script_name} {self.work_dir}/top_layer/")
        self.submit_job(work_dir=f"{self.work_dir}/top_layer")
        os.system(f"cp {self.script_dir}/{self.script_name} {self.work_dir}/bot_layer/")
        self.submit_job(work_dir=f"{self.work_dir}/bot_layer")
        self.wait_until_finished()

    def postprocess(self):
        md_dataset = self.make_dataset()
        md_dataset.save_extxyz(f"{self.work_dir}/MD_data.extxyz")
        return self.make_dataset()

    def make_dataset(self):
        MD_dataset = Dataset()
        MD_dataset.load_dataset_AB(f"{self.input_dir}/ML_AB")
        skip_config = MD_dataset.n_configs
        for stck in self.stackings:
            i = stck[0]
            j = stck[1]
            if self.VASP_ML:
                MD_dataset.load_dataset_AB(f"{self.work_dir}/{i}_{j}/ML_ABN", skip_configs=skip_config)
            else:
                MD_dataset.load_dataset_OUTCAR(f"{self.work_dir}/{i}_{j}/OUTCAR", freq=self.collect_freq)
        if self.VASP_ML:
            MD_dataset.load_dataset_AB(f"{self.work_dir}/top_layer/ML_ABN")
            MD_dataset.load_dataset_AB(f"{self.work_dir}/bot_layer/ML_ABN")
        else:
            MD_dataset.load_dataset_OUTCAR(f"{self.work_dir}/top_layer/OUTCAR", freq=self.collect_freq)
            MD_dataset.load_dataset_OUTCAR(f"{self.work_dir}/bot_layer/OUTCAR", freq=self.collect_freq)
        MD_dataset.save_extxyz(f"{self.work_dir}/MD_data.extxyz")
        return MD_dataset