# DPmoire

DPmoire is an open-source tool designed to construct accurate machine learning force fields (MLFFs) for moiré systems. 

An initial MLFF will firstly be generated on the monolayer structure to ensure stability in MD simulations for the bilayer system. DPmoire will build non-twisted bilayer structures with different stacking, relaxed these structures and ran MD simulation with VASP MLFF module to build training set. An initial MLFF was generated in advance to ensure the stability. During the relaxation, the x and y coordinate of an atom from each layer are kept unchanged. Then we relaxed the twisted structures with DFT to generate the test set. Finally, the MLFF is trained on the collected datasets.

![Fig_1](./images/Fig_1.pdf)

DPmoire is structured into four functional modules. Firstly, as provided the unit cell structures of each layer, DP- moire.preprocess module will automatically combine two layers and generate shifted structures of a 2 $\times$ 2 super- cell. The twisted structure for building test set will also be prepared. The preprocess module will take care of the input files for VASP according to the provided templates. After that, the DPmoire.dft module will submit VASP calculation jobs through slurm system. When all the cal- culation is done, the DFT-calculated data in ML_ABN and OUTCAR files will be collected by DPmoire.data module. Then, DPmoire.data will generate the train- ing set file (data.extxyz) and test set file (valid.extxyz). This format can be directly read by Allegro and NequIP packages. DPmoire.train module will modify the system- dependent settings in configuration file according to given template for training Allegro or NequIP MLFF, and sub- mit the training job. After the training is done, the trained MLFF can be used in ASE or LAMMPS to perform structural relaxation.


![Fig_2](./images/Fig_2.pdf)

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Configuration](#configuration)
  - [Example Setup](#example-setup)
  - [Configuration File](#configuration-file)
  - [Directory Structure and Required Files](#directory-structure-and-required-files)
- [Running DPmoire](#running-dpmoire)
  - [Initial Setup](#initial-setup)
  - [Building Datasets & Training](#building-datasets--training)
  - [Generating Datasets Only](#generating-datasets-only)
  - [Training on Existing Datasets](#training-on-existing-datasets)
  - [Tips for a Fresh Start](#tips-for-a-fresh-start)
- [Possible Issues](#possible-issues)
  - [Loss of the Training Set Does Not Decrease](#loss-of-the-training-set-does-not-decrease)
    - [Case I](#i-check-incar-templates-settings)
    - [Case II](#ii-check-the-temperature-during-md-simulation)
  - [Loss of the Test Set Does Not Decrease](#loss-of-the-test-set-does-not-decrease)
    - [Case I]()
- [Additional Resources](#additional-resources)

## Prerequisites

Before installing DPmoire, ensure that the following software is installed on your system:

- [Slurm](https://slurm.schedmd.com/documentation.html): A workload manager for submitting and managing jobs.
- [VASP](https://www.vasp.at/): A package for performing ab-initio quantum-mechanical molecular dynamics simulations.
- [NequIP](https://github.com/mir-group/nequip): An E(3)-equivariant MLFF package with high data efficiency.
- [Allegro](https://github.com/mir-group/allegro): An E(3)-equivariant MLFF package optimized with large structures.

*note: DPmoire relies on the `NequIP` package to train machine learning force fields (MLFF). If you only need DPmoire to generate the datasets, there is no need to install `NequIP` or `Allegro`*
## Installation

Once `NequIP` is installed, proceed to install DPmoire from source:

```bash
git clone https://github.com/JiaxuanLiu-Arsko/DPmoire.git
cd DPmoire
pip install .
```

## Configuration

### Example Setup

A ready-to-run example is available in the `examples/MoS2` directory. This example provides a comprehensive template to guide you through configuring and running DPmoire for MoS₂ simulations.

### Configuration File

DPmoire requires a configuration file in `.yaml` format to control its operations. An example configuration with detailed comments can be found at `examples/MoS2/config.yaml`.

Here list some important tags:
  - **VASP_ML**: `True/False` Typically it should be True, it will greatly improve the quality of collected Dataset.
  - **do_relaxation**: `True/False` Do relaxation before performing MD simulation. If you are running in this directory for the first time and have not done relaxation before, this tag is recommanded to be True. If you are restarting from relaxed structures, you could set this tag to False.
  - **init_mlff**: `True/False` Build an initial MLFF before runing MD simulation. This tag **should be True** if you are running in this directory for the first time and **did not have ML_AB and ML_FF in `input_dir`**. If you had run this materials with same configuration here or else where, you can copy the ML_AB and ML_FF to `input_dir` and set this tag to False.
  - **symm_reduce**: `True/False` Reduce the shift vectors according to symmetry of the crystal structure. Recommand to be true since this tag can greatly reduce the computational cost.
  - **auto_resub**: `True/False` Automatically re-submit slurm jobs when they are failed/cancelled. You should carefully check your Slurm scripts and make sure everything is fine about the Slurm systems. Otherwise it might repeatedly resubmit a failed job.
  - **n_sectors**:  `int value` Number of grid points to shift. If you set this to 9, DPmoire will calculate $9\times9$ structures (before reducing symmetry) with differet in-plane shift. Typically, the larger value of this tag should improve the quality of dataset.
  - **sc**:  `int value` The size of supercells used in calculation. If you set sc to 2, DPmoire will use $2 \times 2$ supercell to calculate.
  - **d**:  `float value` Initial interlayer distance measured by averaged position of each layers in c direction. This tag will affect generated rigid structure and the **cutoff radius of VASP_MLFF** (ML_RCUT1 & ML_RCUT2 in INCAR) set by DPmoire. It will also affect the rigid interlayer distance of test set and further have impact on the quality of test set. Please carefully set **d** to a reasonable value.


In normal situation, preset values for "CALCULATION SETTINGS" in example should work fine (you should still adjust "ENVIRONMENT SETTINGS" to fit your configuration.). If you are trying to get a better performance, you could try setting larger number of **n_sectors**.

### Directory Structure and Required Files

Before training the MLFF, prepare the following directories and files:

- **Directories:**
  - `input_dir`: Contains essential input files for VASP and NequIP training.
  - `script_dir`: Contains bash scripts for submitting VASP and NequIP jobs in a Slurm environment.

- **Files in `input_dir`:**
  - `init_INCAR`: INCAR template for training the **initial** VASP MLFF.
  - `rlx_INCAR`: INCAR template for relaxation before molecular dynamics (MD).
  - `MD_INCAR`: INCAR template for MD simulations of different stackings.
  - `MD_monolayer_INCAR`: INCAR template for MD simulations of monolayer structures.
  - `val_INCAR`: INCAR template for running test set.
  - `ML_AB` and `ML_FF`: Initial VASP MLFF files. Only required if **init_mlff** is set to *False*, as explained in [Configuration File](#configuration-file)
  - `vdw_kernel.bindat`: VASP van der Waals (vdW) kernel files, necessary when enabling [Nonlocal vdW-DF functionals](https://www.vasp.at/wiki/index.php/Nonlocal_vdW-DF_functionals) in INCAR templates. Refer to the **Important Technical Remarks** section for more details.
  - `bot_layer.poscar` and `top_layer.poscar`: POSCAR files for the unit cells of the bottom and top layers, respectively. **Ensure the c-axis of POSCAR is sufficiently large!**

## Running DPmoire

### Initial Setup

**Before running the example, update the following:**

1. **POTCAR Directory:**
   - Modify the `POTCAR_dir` tag in the configuration file to point to the directory containing your VASP POTCAR files.

2. **Slurm Scripts:**
   - Update the Slurm scripts located in the `scripts` directory to match your computing environment.
   - Adjust the `learn_script` and `DFT_script` paths in the configuration file accordingly.

### Building Datasets & Training

To build datasets and train the MLFF, navigate to the example directory and execute the training command:

```bash
cd ./example/MoS2
nohup DPmoireTrain ./config.yaml &
```

Alternatively, to run all processes in one command:

```bash
nohup DPmoireTrain ./config.yaml --mode all &
```

This will initiate DPmoire as a background job, submitting calculation and learning tasks to the compute nodes. Upon completion, the trained NequIP MLFF will be located at `work_dir/main/mlff.pth`.

### Generating Datasets Only

If you prefer to generate datasets without submitting training jobs, use the following command:

```bash
nohup DPmoireTrain ./config.yaml --mode run &
```

This will create datasets in `work_dir/rlx_data.extxyz` and `work_dir/MD_data.extxyz`, corresponding to data collected during relaxation and MD processes, respectively. The complete dataset will be merged in `./main/data/`, and the `nequIP.yaml` configuration will be appropriately set. You can then transfer the `main` directory to a GPU-enabled machine to train the MLFF separately:

```bash
cd ./main
nequip-train ./nequIP.yaml
```

Alternatively, you can submit a training job using a Slurm script.

### Training on Existing Datasets

To train a NequIP MLFF based on directories where VASP calculations have already been performed, use:

```bash
nohup DPmoireTrain ./config.yaml --mode train &
```

DPmoire will search for `rlx_data.extxyz` and `MD_data.extxyz`. If `MD_data.extxyz` is absent, it will attempt to reconstruct the dataset from `ML_AB` or `OUTCAR` files generated during VASP calculations.



### Optimizing MLFF Performance

For optimal performance, train the MLFF using untwisted data and reserve some twisted structures as a validation set. Currently, there are two approaches:

1. **Automated Validation Set:**

   - Set `twist_val: True` in `config.yaml`.
   - Specify appropriate values for `max_val_n` and `min_val_n`.
   - DPmoire will automatically handle the validation set generation.
   - *Note: This method only works when both `top_layer` and `bot_layer` are the same material.*

2. **Manual Validation Set:**

   - Set `twist_val: False` in `config.yaml`.
   - Perform an additional VASP calculation for some twisted structures with similar settings.
   - Collect the data from these calculations as the validation set.
   - Convert the `OUTCAR` files to `extxyz` format using the following script:

     ```python
     from DPmoire.data import Dataset

     dataset = Dataset()
     dataset.load_dataset_OUTCAR("./OUTCAR", 1)
     dataset.save_extxyz("./valid.extxyz")
     print(f"Number of configurations: {dataset.n_configs}")
     ```

   - Update `work_dir/main/nequip.yaml` with the validation set details:

     ```yaml
     validation_dataset: ase
     validation_dataset_file_name: ./data/valid_set.extxyz
     n_val: 20
     ```

### Validation Set Configuration

Ensure that the `nequip.yaml` file accurately reflects the location and properties of your validation dataset. Adjust `validation_dataset_file_name` and `n_val` according to your specific validation data.

## Tips for a Fresh Start

When initiating a new training session for a different material, ensure the following:


1. **POSCAR Files:**
   - Modify `top_layer.poscar` and `bot_layer.poscar` to desired material
   - **Ensure the c-axis in POSCAR is sufficiently large**

2. **INCAR Files:**
   - Update all `*_INCAR` templates with appropriate van der Waals (vdW) correction settings.

3. **Configuration File:**
   - Update `config.yaml` with the correct interlayer distance (`d`) to appropriate value.

DPmoire will automatically handle the following configurations:

- **ENCUT Adjustment:**
  - ENCUT in INCAR will be set to `1.6 × ENMAX`, where `ENMAX` is the maximum plane-wave cutoff energy.

- **POTCAR Generation:**
  - POTCAR files will be automatically generated based on the provided single-layer structures.

- **KPOINTS Generation:**
  - KPOINTS will be generated to satisfy the condition `n_k_i × a_i ≈ 40`, ensuring adequate sampling in reciprocal space.

## Tips for a Restart job

When you are restarting the job in an old directory, there is something you should notice: 

1. If you had changed some DFT parameters, remeber to update `d` in configuration value.

2. If you wanna restart from MD simulation, you could set the `init_mlff` and `do_relaxation` to False, and then start DPmoire as forementioned.

3. If you mannually re-submitted some calculation for some reason and want to grub the new calculation results, you should:
   - Remove/move/rename *.extxyz  in `work_dir`, otherwise DPmoire will grub the data from thesefiles.
   - Use [mode==train](#training-on-existing-datasets) to restart DPmoire

## Additional Resources

- [NequIP GitHub Repository](https://github.com/mir-group/nequip)
- [Allegro GitHub Repository](https://github.com/mir-group/allegro)
- [VASP Documentation](https://www.vasp.at/wiki/index.php/Main_Page)
- [Slurm Documentation](https://slurm.schedmd.com/documentation.html)

## Possible Issues

### Loss of the Training Set Does Not Decrease

If you notice that the training loss decreases slowly or that $F_{rmse}$ is difficult to reduce below $1\times 10^{-2} eV/Å$ during the MLFF training process, it's typically due to incompatible structures/data occurs in training set.

#### I: Check INCAR templates settings

The first thing to do in this situation is to check if your VASP settings in INCAR templates (e.g. example/MoS2/input/*_INCAR) are compatible. If you used IVDW=10 to relax and SCAN+RVV10 vdW corrections to perform MD simulation, the calculated data will be incompatible.

#### II: Check the temperature during MD simulation

If case I does not work for you, then you should check the temperature during the MD simulation. VASP MLFF module can be unstable sometimes and will run into strange structures. Though generating initial MLFF in advance (as DPmoire do) can greatly reduce the risk, there might still be a chance to have this problem.

You can check the temperature during the MD simulation to find this issue. In normal situation, the highest temperature can be slightly higher then preset (e.g. 330K preset but 580K highest during the run). If in some trajectories, the temperature exceeds 1000K, there must be a problem.

You can use the following script to find the

```python
def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def findTemp(i, j, outfile):
    with open(f"./{i}_{j}/OUTCAR", 'r') as f:
        max_temp = None
        for line in f:
            if 'temperature' in line:
                tempstr = line.split()[-2]
                if is_float(tempstr):
                    temp = float(tempstr)
                    if max_temp is None or temp > max_temp:
                        max_temp = temp
    if max_temp is not None:
        outfile.write(f"{i}_{j} {max_temp}\n")
        if max_temp > 1000:
            print(f"{i}_{j}: {max_temp}")
    else:
        print('No temperature information found in the file')

n = 9
with open("findTemp.dat", "w") as outfile:
    for i in range(n):
        for j in range(n):
            try: 
                findTemp(i, j, outfile)
            except:
                pass
```

If temperatures exceed 1000 K, they will be printed to the console, indicating potential issues with the data.

You can go to the directories that with problem, re-submit DFT calculation manually. 

After the calculation is done, you should **remove the MD_data.extxyz in your `work_dir`** and [train again](#training-on-existing-datasets) as described in [here](#tips-for-a-restart-job)


