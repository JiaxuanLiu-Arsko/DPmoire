import numpy as np
import os, time
import yaml
from ase.atoms import Atoms
from ..data import Dataset
from ..preprocess import Config

class MLFFTrainer:
    work_dir = None
    config_file_template = None
    script_dir = None
    learn_script = None
    dataset = None
    input_dir = None
    val_dataset = None
    RCUT = None
    elements = None
    sc=None

    def __init__(self, config:Config, dataset:Dataset, val_dataset:Dataset=None):
        self.work_dir = config["work_dir"] + "/main"
        if not os.path.exists(self.work_dir):
            os.mkdir(f"{self.work_dir}")
        self.learn_script = config["learn_script"]
        self.input_dir = config["input_dir"]
        self.script_dir = config["script_dir"]
        self.sc = config["sc"]
        self.d = config["d"]
        self.dataset = dataset
        self.val_dataset = val_dataset

    def get_params(self):
        step = int(self.dataset.n_configs/100)-1
        elements = []
        max_d = 0
        for atoms in self.dataset[0:step:self.dataset.n_configs]:
            for i, pos_i in enumerate(atoms.get_positions(wrap=True)):
                for j, pos_j in enumerate(atoms.get_positions(wrap=True)):
                    dis_ij = np.sqrt(np.sum([(pos_i[dim] - pos_j[dim])**2 for dim in range(2)]))/self.sc
                    dis_ij = np.sqrt(dis_ij**2 + self.d**2)
                    if dis_ij > max_d:
                        max_d = dis_ij
                if atoms.get_chemical_symbols()[i] not in elements:
                    elements.append(atoms.get_chemical_symbols()[i])
        self.RCUT = max_d * 1.1
        self.elements = elements

    def make_dataset_file(self):
        pass

    def make_mlff_config(self):
        pass

    def preprocess(self):
        self.get_params()
        self.make_dataset_file()
        self.make_mlff_config()

    def submit_training(self):
        os.chdir(self.work_dir)
        jobid = os.popen(f"sbatch {self.learn_script}").readlines()[0].split()[-1]
        time.sleep(3)
        flag = len(os.popen(f"squeue | grep {jobid}").readlines()) == 0
        while flag:
            time.sleep(30)
            flag = len(os.popen(f"squeue | grep {jobid}").readlines())
    
    def postprocess(self):
        pass
