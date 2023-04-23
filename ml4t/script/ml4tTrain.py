import argparse, os, yaml, time

from ml4t.preprocess import Shifter, Envgen, Parameters, Minimizer
from ml4t.data import ABRunner, Dataset
from ml4t.train import Trainer
from ml4t.test import getAseForces, getVaspForces, findError
def main(args=None):

    parser = argparse.ArgumentParser(
        description="A module to automatically train a Machine Learning forcefield for a Twisted Bilayer Material"
    )
    parser.add_argument(
        "config", help="YAML file configuring the module"
    )

    parser.add_argument(
        "--runMode", default="all", type=str, help="runMode: all, ab, or train."
    )

    args = parser.parse_args(args=args)

    if not os.path.isfile(args.config):
        print("YAML not found!")
    params = Parameters(args.config)
    startTime = time.time()
    dataset = Dataset()
    if args.runMode == "all" or args.runMode == "ab":
        minimizer = Minimizer(params=params)
        minimizer.loadEnv()
        minimizer.runMinimize()
        if params.collectMinDat:
            dataset.loadAll_OUTCAR(nSecs=params.nSecs, freq=params.collectFreq, includeDir=params.workDir)
        print(f"Minimize finished. Time cost = {time.time()-startTime} secs.")
        dataset.save_npz(f"{params.workDir}/data.npz")
        startTime = time.time()
        abrunner = ABRunner(params=params)
        abrunner.loadEnv()
        abrunner.runAB()
        print(f"MD finished. Time cost = {time.time()-startTime} secs.")
        if params.collectMode == "AB":
            dataset.loadAll_AB(nSecs=params.nSecs, includeDir=params.workDir)
        elif params.collectMode == "OUTCAR":
            dataset.loadAll_OUTCAR(nSecs=params.nSecs, freq=params.collectFreq, includeDir=params.workDir)
        else:
            print("Warning!No Dataset in MD Loaded!")
        dataset.save_npz(f"{params.workDir}/data.npz")
        startTime = time.time()
    if args.runMode == "all" or args.runMode == "train":
        if args.runMode == "train":
            if os.path.exists(f"{params.workDir}/data.npz"):
                dataset.startFromNPZ(f"{params.workDir}/data.npz")
                dataset.nConfigs = len(dataset.data["energies"])
                print("restart from NPZ.")
            else:
                dataset.loadAll_AB(params.nSecs, params.workDir)
                dataset.nConfigs = len(dataset.data["energies"])
                print("restart from AB.")
        trainer = Trainer(params=params)
        trainer.makeDataset(dataset)
        trainer.genNequEnv()
        trainer.startTrain()
        print(f"Training finished. Time cost = {time.time()-startTime} secs.")
        trainer.deploy()

def loadEnv(params:Parameters):
    shft = Shifter(param=params)
    shft.shiftAll(outDir=params.workDir)
    envgen = Envgen(params=params)
    envgen.genAll(genDir=params.workDir, params=params)

if __name__ == "__main__":
    main()