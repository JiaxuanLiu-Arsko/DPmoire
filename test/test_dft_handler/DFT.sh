#!/bin/bash
#
#SBATCH -N 1 ###
#SBATCH -n 1 ###
#SBATCH --job-name=DPmoire
#SBATCH --partition=short
#SBATCH --output=./log
#SBATCH --error=./err
#SBATCH --exclude hpcc037
echo "finished"