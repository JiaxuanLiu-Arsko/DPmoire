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
    do_rlx = False
    collectFreq = 0
    collectMode = ""
    def __init__(self, configFile:str):
        with open(configFile, "r") as f:
            config = yaml.load(f, Loader=yaml.Loader)
        self.nSecs = config["sectors"]
        self.workDir = os.path.abspath(config["work_dir"])
        self.originDir = os.path.abspath(config["origin_dir"])
        self.potDir = os.path.abspath(config["POTCAR_dir"])
        self.nThreads = config["n_threads"]
        self.attachDir = os.path.abspath(config["attach_dir"])
        self.scriptMD = config["MD_script"]
        self.userName = config["username"]
        self.scriptMD = config["MD_script"]
        self.userName = config["username"]
        self.scriptLearn = config["learn_script"]
        self.testDir = os.path.abspath(config["test_dir"])
        self.minimizeDir = os.path.abspath(config["relaxation_dir"])
        if config["do_relaxation"] is None:
            self.do_rlx = False
        else:
            self.do_rlx = config["do_relaxation"]
        if config["collect_freq"] is None:
            self.collectFreq = 5
        else:
            self.collectFreq = config["collect_freq"]
        if config["collect_mode"] is None:
            self.collectMode = "AB"
        else:
            self.collectMode = config["collect_mode"]
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
