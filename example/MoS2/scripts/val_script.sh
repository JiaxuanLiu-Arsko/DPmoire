#!/bin/bash
#
#SBATCH -N 2 ###
#SBATCH -n 104 ###
#SBATCH -A hmt03
#SBATCH --partition=short
#SBATCH --job-name=compare_vdw
#SBATCH --output=./log
#SBATCH --error=./err
### 使用的 gpu 卡数，与--gres=gpu:4 的数字一致
NP=$((104))
### 执行任务所需要加载的模块
module load oneapi22.3
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
export LD_LIBRARY_PATH=/usr/local/lib64:$LD_LIBRARY_PATH
export CUDA_LAUNCH_BLOCKING=1
echo $LD_LIBRARY_PATH
### 执行任务

if [ -f "/home/jxliu/anaconda3/etc/profile.d/conda.sh" ]; then
    . "/home/jxliu/anaconda3/etc/profile.d/conda.sh"
 else
     export PATH="/home/jxliu/anaconda3/bin/:$PATH"
 fi

mpirun -np $NP /home/jxliu/tools/vasp.6.3.0_ljx/vasp.6.3.0/bin/vasp_std