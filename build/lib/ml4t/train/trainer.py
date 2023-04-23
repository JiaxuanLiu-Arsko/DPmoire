import numpy
import os, time
import yaml

from MLFF4TBM.data import Dataset
from MLFF4TBM.preprocess import Parameters

class Trainer:
    workDir = None
    originfile = None
    attachDir = None
    scriptLearn = None
    userName = None
    RCUT = 0
    nConfigs = 0

    def __init__(self, params:Parameters):
        self.nSecs = params.nSecs
        self.includeDir = params.workDir
        self.workDir = params.workDir + "/main"
        self.attachDir = params.attachDir
        self.scriptLearn = params.scriptLearn
        self.RCUT = params.RCUT
        self.elements = params.elements

    def makeDataset(self, dataset:Dataset):
        if not os.path.exists(self.workDir + "/data"):
            os.makedirs(self.workDir + "/data")
        dataset.save_extxyz(saveSTR=self.workDir + "/data/data.extxyz")
        self.nConfigs = dataset.nConfigs
        
    def genNequEnv(self):
        with open(self.attachDir + "/nequIP.yaml", "r") as f:
            nequConfig = yaml.load(f, Loader=yaml.Loader)
        nequConfig["r_max"] = self.RCUT
        nequConfig["chemical_symbols"] = self.elements
        nequConfig["n_train"] = int(self.nConfigs*0.8)
        nequConfig["n_val"] = int(self.nConfigs*0.2)
        with open(self.attachDir + "/nequIP.yaml", "w") as f:
            yaml.dump(nequConfig, f)
        os.system("cp -r " + self.attachDir + "/* " + self.workDir)

    def startTrain(self):
        os.chdir(self.workDir)
        jobstr = os.popen("sbatch " + self.scriptLearn)
        time.sleep(3)
        jobid = jobstr.readlines()[0].split()[-1]
        missions = os.popen("squeue | grep " + jobid)
        while len(missions.readlines())>0:
            time.sleep(30)
            missions = os.popen("squeue | grep " + jobid)
    
    def deploy(self):
        os.system(f"nequip-deploy build --train-dir {self.workDir}/result/model {self.workDir}/mlff.pth")
        time.sleep(30)
