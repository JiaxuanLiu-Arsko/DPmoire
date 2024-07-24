from DPmoire.data import Dataset
dataset = Dataset()
dataset.load_dataset_AB("./ML_AB")
#dataset.load_dataset_OUTCAR("./OUTCAR", freq=7)
dataset.save_extxyz("./dataset.extxyz")