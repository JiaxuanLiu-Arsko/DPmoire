#!/bin/bash
#
#SBATCH -N 1 ### 使用的节点数目
#SBATCH -n 1 ### 一个节点内部的 32 个核
#SBATCH --partition=h800
#SBATCH --job-name=test
#SBATCH --output=./log
#SBATCH --error=./err
#SBATCH -A hmt03
### 执行任务所需要加载的模块

echo job_finished