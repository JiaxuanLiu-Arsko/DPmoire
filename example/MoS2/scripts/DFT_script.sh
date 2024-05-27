#!/bin/bash
#
#SBATCH -N 1 ###
#SBATCH -n 52 ###
#SBATCH --job-name=DPmoire
#SBATCH --partition=short
#SBATCH --output=./log
#SBATCH --error=./err
NP=$((52))
### 执行任务所需要加载的模块
module load oneapi2022/mpi/latest
module load oneapi2022/mkl/latest
module load oneapi2022/compiler/latest
#module load fftw3/3.3.8
#module load ohpc
#module load vasp6/6.4    #VASP_MLFF模块仅在vasp6以上版本才有。建议使用6.3.0以上版本
### 一些提示性输出
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MV2_ENABLE_AFFINITY=0
echo ”The current job ID is $SLURM_JOB_ID”
echo ”Running on $SLURM_JOB_NUM_NODES nodes:”
echo $SLURM_JOB_NODELIST
echo ”Using $SLURM_NTASKS_PER_NODE tasks per node”
echo ”A total of $SLURM_NTASKS tasks is used”
### 对任务执行的内存不做限制
ulimit -s unlimited
ulimit -c unlimited
### 加载任务所需要的库
#export LD_LIBRARY_PATH=/usr/local/lib64:$LD_LIBRARY_PATH
export CUDA_LAUNCH_BLOCKING=1
echo $LD_LIBRARY_PATH
### 执行任务

mpirun -np $NP /home/users/jxliu/tools/VASP/vasp.6.3.0/bin/vasp_std