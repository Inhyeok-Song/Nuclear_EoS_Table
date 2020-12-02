#!/bin/bash
#SBATCH --partition gpu
#SBATCH --nodelist g03
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
##SBATCH --mem 64G
#SBATCH --job-name GPU_EoS_Driver_test
#SBATCH --output ./log
#SBATCH --time 3-0

srun ./driver
