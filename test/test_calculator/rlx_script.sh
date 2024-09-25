#!/bin/bash
#
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --gres=gpu:1 ### 使用 1 张 gpu 卡
#SBATCH --partition=h800
#SBATCH -A hmt03
#SBATCH --job-name=rlx_WSe2
#SBATCH --output=./log
#SBATCH --error=./err ### 错误日志文件. nequIP的正常输出会在错误日志中。

NP=$((1)) ### 使用的 gpu 卡数，与--gres=gpu:1 的数字一致
### 执行任务所需要加载的模块
module load oneapi22.3
module load nvhpc/22.11
module load cuda11.8
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
#ulimit -c unlimited
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
python rlxTest.py