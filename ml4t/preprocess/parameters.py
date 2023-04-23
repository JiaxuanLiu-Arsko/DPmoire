import os, yaml, re
import numpy as np

class Parameters:
    testDir = ""
    nSecs = 0
    workDir = ""
    originDir = ""
    minimizeDir = ""
    potDir = ""
    nThreads = ""
    attachDir = ""
    scriptMD = ""
    userName = ""
    scriptLearn = ""
    RCUT = 0
    elements = []
    latVec = [[], []]
    ENCUT = 0
    collectMinDat = False
    collectFreq = 0
    collectMode = ""
    def __init__(self, configFile:str):
        with open(configFile, "r") as f:
            config = yaml.load(f, Loader=yaml.Loader)
        self.nSecs = config["Sectors"]
        self.workDir = os.path.abspath(config["WorkDir"])
        self.originDir = os.path.abspath(config["originDir"])
        self.potDir = os.path.abspath(config["POTCARDir"])
        self.nThreads = config["nThreads"]
        self.attachDir = os.path.abspath(config["attachDir"])
        self.scriptMD = config["MDscript"]
        self.userName = config["username"]
        self.scriptMD = config["MDscript"]
        self.userName = config["username"]
        self.scriptLearn = config["scriptLearn"]
        self.testDir = os.path.abspath(config["testDir"])
        self.minimizeDir = os.path.abspath(config["minimizeDir"])
        if config["collectMinDat"] is None:
            self.collectMinDat = False
        else:
            self.collectMinDat = config["collectMinDat"]
        if config["collectFreq"] is None:
            self.collectFreq = 5
        else:
            self.collectFreq = config["collectFreq"]
        if config["collectMode"] is None:
            self.collectMode = "AB"
        else:
            self.collectMode = config["collectMode"]
        self.getParams()

    def getParams(self):
        with open(self.originDir + "/POSCAR",'r') as f:
            poscar = f.read()
        dstr = re.search("d\s*=\s*\d+\.*\d*", poscar.split("\n")[0]).group()
        self.RCUT = float(re.search("\d+\.*\d*", dstr).group())+2.5
        for idx in range(2):
            for dim in range(2):
                self.latVec[idx].append(float(poscar.split("\n")[idx+2].split()[dim]))
        self.elements = poscar.split("\n")[5].split()
        ens = []
        for element in self.elements:
            ens.append(os.popen("cat " + self.potDir + "/" + element + "/POTCAR|grep ENMAX").readlines()[0].split()[2].split(";")[0]) 
        ensf = [float(ens[k]) for k in range(len(ens))]
        self.ENCUT = np.max(ensf)
