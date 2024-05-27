import numpy
import os, time
import yaml

from ..data import Dataset
from ..preprocess import Config
from .mlff_trainer import MLFFTrainer

class NequIPTrainer(MLFFTrainer):

    def __init__(self, config:Config, dataset:Dataset):
        super().__init__(config=config, dataset=dataset)
        self.config_file_template = "nequIP.yaml"
        

    def make_dataset_file(self):
        if not os.path.exists(f"{self.work_dir}/data/"):
            os.mkdir(f"{self.work_dir}/data/")
        self.dataset.save_extxyz(f"{self.work_dir}/data/data.extxyz")

    def make_mlff_config(self, RCUT):
        with open(f"{self.input_dir}/{self.config_file_template}", "r") as f:
            config_tmp = yaml.load(f, Loader=yaml.FullLoader)
        config_tmp["r_max"] = float(RCUT)
        config_tmp["chemical_symbols"] = self.elements
        config_tmp["n_train"] = int(self.dataset.n_configs*0.8)
        config_tmp["n_val"] = int(self.dataset.n_configs*0.2)
        with open(f"{self.work_dir}/nequIP.yaml", "w") as f:
            yaml.dump(config_tmp, f)
    
    def preprocess(self, RCUT):
        self.get_params()
        self.make_dataset_file()
        self.make_mlff_config(RCUT)
        os.chdir(self.work_dir)
        os.system(f"cp {self.script_dir}/{self.learn_script} ./")
    
    def postprocess(self):
        os.chdir(self.work_dir)
        os.system(f"nequip-deploy build --train-dir {self.work_dir}/result/model {self.work_dir}/mlff.pth")
        time.sleep(30)
