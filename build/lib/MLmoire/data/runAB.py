import numpy as np
import os
import time
from MLmoire.preprocess import Parameters, Shifter, Envgen

class ABRunner:

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
    
    def runAB(self):
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
        
    def loadEnv(self):
        shft = Shifter(param=self.params)
        shft.shiftAll(outDir=self.workDir)
        envgen = Envgen(params=self.params)
        envgen.genAll(genDir=self.workDir, params=self.params)


                
        

