# DPmoire

## Installation

`DPmoire` requires the `nequip` package and its dependencies; please see the [NequIP installation instructions](https://github.com/mir-group/nequip#installation) for details.

[Allegro](https://github.com/mir-group/allegro) is optional to be installed. It is a kind of updated version of `NequIP`, which allows you to use the MLFF through MPICH. 

Once `nequip` is installed, you can install `DPmoire` from source by running:
```bash
git clone https://github.com/JiaxuanLiu-Arsko/DPmoire.git
cd DPmoire
pip install .
```

## Configuration

An ready-to-run example can be found in **examples/MoS2**

### Configuration file

A *.yaml format file are required for controll the DPmoire package.
An example with full comments can be found in **examples/MoS2/config.yaml**

### Files
Several directories should be prepared before you start to train MLFF.

  - input: A directory contains essential input files for VASP and NequIP training.
    - init_INCAR*: INCAR for training **initial** VASP MLFF.

    - rlx_INCAR: INCAR for relaxation before MD, only required when `do_relaxation` was set to True.

    - MD_INCAR: INCAR for MD simulation of different stackings.

    - MD_monolayer_INCAR: INCAR for MD simulation of monolayer structures.

    - ML_AB and ML_FF: Initial VASP MLFF files. Only required when `VASP_ML` was set to `True` and `init_mlff` was set to `False`.

    - vdw_kernel.bindat: VASP vdw Kernel files. See **Important technical remarks** section in [Nonlocal vdW-DF functionals](https://www.vasp.at/wiki/index.php/Nonlocal_vdW-DF_functionals)

    - bot_layer.poscar and top_layer.poscar: Unit cell POSCAR of bottom layer and bottom layer atoms. **The c axis of poscar should be large enough!**
  - scripts: A directory contains bash scripts to submit VASP and NequIP job in slurm system.

## Running DPmoire

**Before you run the example, remeber to change two things:**
  - `POTCAR_dir` tag in config file
  - slurm scripts in `scripts` directory, and accordingly `learn_script` and `DFT_script` in config file.

After that, it should be OK to run DPmoire by following commands:

```bash
cd ./example/MoS2
nohup DPmoireTrain ./config.yaml &
```

or

```bash
nohup DPmoireTrain ./config.yaml --mode all &
```

Then, it should be running as background. The calculation and learning job will be submitted to the calculation nodes. 

After all calculation (Making dataset & Training MLFF) was done, the NequIP type MLFF file should be in [`work_dir`/main/mlff.pth](./example/MoS2/main/mlff.pth).

Also, you can use the following command to only get a dataset:

```bash
nohup DPmoireTrain ./config.yaml --mode run &
```

The dataset will be in [`work_dir`/rlx_data.extxyz](./example/MoS2/) and [`work_dir`/MD_data.extxyz](./example/MoS2/), corresponds to data collected during relaxation and MD process.

The whole dataset are merged in ./main/data/dataset.extxyz, and the nequIP.yaml was set properly. You can copy the "main" directory to the gpu machine and train the MLFF seperately:

```
  cd ./main
  nequip-train ./nequIP.yaml
```
or maybe use a script to submit a job.

You can use the following command to train a NequIP MLFF based in a directories where you had run the VASP calculation before:

```bash
nohup DPmoireTrain ./config.yaml --mode train &
```

DPmoire will try to find rlx_data.extxyz and MD_data.extxyz. If there is no MD_data.extxyz, it will try to reconstruct dataset from ML_AB or OUTCAR files generated in VASP calculation.

### Tips for new training
Before you try to train another material, these files **must** be taken care of :
  - *_INCAR: remeber to change the vdw functional settings to appropriate one. 
  - top_layer.poscar, bot_layer.poscar
  - config.yaml for DPmoireTrain. Tag "d" is the approximate value of interlayer distance. The shifted structure will be generated according to this, and the cutoff of MLFF will also be selected by this.

There's something the scripts will automatically configure:
  - ENCUT of INCAR will be automatically taken care of by ENCUT=1.2\*ENMAX.
  - POTCAR will be automatically generated
  - KPOINTS will be generated according to $n_{k_i}\cdot a_{i}>40$


Then it's ready to run a new train.

After the calculation is over, you can run a python script to check the temperature during the MD:

```python
def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def findTemp(i,j,outfile):
    with open("./"+str(i)+"_"+str(j)+'/OUTCAR','r') as f:
        max_temp = None
        for line in f:
            if 'temperature' in line:
                tempstr = line.split()[-2]
                if is_float(tempstr):
                    temp = float(tempstr)
                    if max_temp is None or temp > max_temp:
                        max_temp = temp
        if max_temp is not None:
            outfile.write(str(i)+"_"+str(j)+" ")
            outfile.write(str(max_temp))
            outfile.write("\n")
            if max_temp > 1000:
                print(f"{i}_{j}: {max_temp}")
        else:
            print('No temperature information found in the file')

n=9
outfile = open("findTemp.dat","w")

for i in range(n):
    for j in range(n):
        findTemp(i,j,outfile)

```

If the temperature is fine (nothing printed), then the data should be good to use.

To obtain the best performance, one should take the dataset as training set (for NequIP) and run another VASP calculation for some big angle twisted materials (with similar settings) and collect those data as validation set. Following script convert OUTCAR to extxyz format dataset:

```python
from DPmoire.data import Dataset
dataset = Dataset()
outcarStructures = dataset.load_dataset_OUTCAR(f"./OUTCAR", 1)
dataset.save_extxyz("./valid.extxyz")
print(dataset.n_configs)
```

After that, the you should add 
```yaml
validation_dataset: ase
validation_dataset_file_name: ./data/valid_set.extxyz
n_val: 20
```
according to your validation set. The n_train should be set to the number of data in dataset.extyz.