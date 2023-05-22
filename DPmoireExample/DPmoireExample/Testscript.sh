#!/bin/bash
#
#SBATCH -N 1 ### 使用的节点数目
#SBATCH -n 32 ### 一个节点内部的 32 个核
#SBATCH --gres=gpu:4 ### 使用 4 张 gpu 卡
#SBATCH --partition=amdgpu40g
#SBATCH --job-name=ml4tTest
#SBATCH --output=./log
#SBATCH --error=./err
NP=$((4))       #如果是GPU节点，NP应当等于使用的卡数
### 执行任务所需要加载的模块
module load oneapi22.3
module load nvhpc/22.11
module load cuda11.8    
module load vasp/vasp6.3.0_wannier90-1.2_amd

# 不要设置OMP_NUM_THREADS=1。ASE环境的并行依赖OMP，设置OMP_NUM_THREADS=1会让MLFF的计算变得缓慢。
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

###防止conda activate environment失败，网上抄的
if [ -f "/home/jxliu/anaconda3/etc/profile.d/conda.sh" ]; then
    . "/home/jxliu/anaconda3/etc/profile.d/conda.sh"
else
     export PATH="/home/jxliu/anaconda3/bin/:$PATH"
fi

### 执行任务
conda activate nequIP
DPmoireTest config.yaml --runMode run