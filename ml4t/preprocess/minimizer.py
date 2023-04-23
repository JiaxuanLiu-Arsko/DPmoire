import numpy as np
import os
import time
from .parameters import Parameters
from .envgen import Envgen
from .shifter import Shifter

class Minimizer:

    nSecs = 0
    workDir = None
    originDir = None
    potDir = None
    nThreads = 0
    attachDir = None
    scriptMD = None
    userName = None
    params = None
    def __init__(self, params:Parameters):
        self.params = params
        self.nSecs = params.nSecs
        self.workDir = params.workDir
        self.originDir = params.originDir
        self.potDir = params.attachDir
        self.nThreads = params.nThreads
        self.attachDir = params.attachDir
        self.scriptMD = params.scriptMD
        self.userName = params.userName
    
    def runMinimize(self):
        jobs = []
        for i in range(self.nSecs):
            for j in range(self.nSecs):
                os.chdir(self.workDir+"/"+str(i)+"_"+str(j))
                os.system("cp -r " + self.attachDir + "/*" + " ." )
                jobid = os.popen("sbatch " + self.scriptMD).readlines()[0].split()[-1]
                jobs.append(jobid)
                missions = os.popen("squeue | grep " + self.userName)
                while len(missions.readlines())>=self.nThreads:
                    time.sleep(30)
                    missions = os.popen("squeue | grep " + self.userName)

        missions = 0
        for job in jobs:
            missions += len(os.popen(f"squeue | grep {job}").readlines())
        while missions>0:
            time.sleep(30)
            missions = 0
            for job in jobs:
                missions += len(os.popen(f"squeue | grep {job}").readlines())
        
        self.contToPos()
        
    def loadEnv(self):
        params = self.params
        shft = Shifter(param=params)
        shft.shiftAll(outDir=params.workDir)
        envgen = Envgen(params=self.params)
        for i in range(self.nSecs):
            for j in range(self.nSecs):
                kpoints = envgen.kpgenerator(params.latVec)
                envgen.genPOTCAR(elements=params.elements, outDir=f"{params.workDir}/{i}_{j}")
                envgen.genKPOINTS(kpoints, outDir=f"{params.workDir}/{i}_{j}")
                envgen.genINCAR(outDir=f"{params.workDir}/{i}_{j}", params=params, originDir=params.minimizeDir)
    
    def contToPos(self):
        params = self.params
        for i in range(self.nSecs):
            for j in range(self.nSecs):
                outfile = open(f"{params.workDir}/{i}_{j}/POSCAR", "w")
                with open(f"{params.workDir}/{i}_{j}/CONTCAR", "r") as infile:
                    lines = infile.readlines()
                nLine = int(len(lines)/2) + 4
                for line in lines[0:nLine]:
                    outfile.write(line)
                outfile.close()
                    


                
        

