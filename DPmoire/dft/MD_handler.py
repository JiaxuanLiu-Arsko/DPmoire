import numpy as np
import os
import time
from ..preprocess import Config
from .dft_handler import DFTHandler
from ..data import Dataset

class MDHandler(DFTHandler):
    n_secs_layer2 = None
    n_secs_layer3 = None
    work_dir = None
    script_dir = None
    VASP_ML = None
    collect_freq = None
    input_dir = None
    stackings = None
    n_layers = 3
    def __init__(self, config:Config|dict, existing_job:list=None, stackings = None):
        super().__init__(script_name=config["DFT_script"], n_nodes=config["n_nodes"], existing_job=existing_job, auto_resub=config["auto_resub"])
        self.n_secs_layer2 = config["n_sectors_layer2"]
        self.n_secs_layer3 = config["n_sectors_layer3"]
        self.work_dir = config["work_dir"]
        self.script_dir = config["script_dir"]
        self.VASP_ML = config["VASP_ML"]
        self.collect_freq = config["OUTCAR_collect_freq"]
        self.input_dir = config['input_dir']
        self.n_layers = 3
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
        if len(self.stackings) > 0:
            self.n_layers = int(len(self.stackings[0]) / 2) + 1
        else:
            self.n_layers = 3

    def run_calculation(self):
        for stck in self.stackings:
            i2, j2, i3, j3 = map(int, stck)
            target_dir = f"{self.work_dir}/{i2}_{j2}_{i3}_{j3}"
            os.system(f"cp {self.script_dir}/{self.script_name} {target_dir}/")
            self.submit_job(work_dir=target_dir)
        for layer_idx in range(self.n_layers):
            layer_dir = f"{self.work_dir}/layer_{layer_idx+1}"
            os.system(f"cp {self.script_dir}/{self.script_name} {layer_dir}/")
            self.submit_job(work_dir=layer_dir)
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
            i2, j2, i3, j3 = map(int, stck)
            target_dir = f"{self.work_dir}/{i2}_{j2}_{i3}_{j3}"
            if self.VASP_ML:
                MD_dataset.load_dataset_AB(f"{target_dir}/ML_ABN", skip_configs=skip_config)
            else:
                MD_dataset.load_dataset_OUTCAR(f"{target_dir}/OUTCAR", freq=self.collect_freq)
        for layer_idx in range(self.n_layers):
            layer_dir = f"{self.work_dir}/layer_{layer_idx+1}"
            if self.VASP_ML:
                MD_dataset.load_dataset_AB(f"{layer_dir}/ML_ABN")
            else:
                MD_dataset.load_dataset_OUTCAR(f"{layer_dir}/OUTCAR", freq=self.collect_freq)
        MD_dataset.save_extxyz(f"{self.work_dir}/MD_data.extxyz")
        return MD_dataset
