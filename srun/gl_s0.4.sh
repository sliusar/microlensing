#!/bin/env bash
 
#SBATCH --job-name gravlens_s0.4
#SBATCH --error /home/sliusar/microlensing/logs/%x_%j_%a.err
#SBATCH --output /home/sliusar/microlensing/logs/%x_%j_%a.out
#SBATCH --partition=shared-gpu
#SBATCH --time=10:00:00
#SBATCH --gpus=1

module load GCC/8.3.0  CUDA/10.1.243  OpenMPI/3.1.4 CMake/3.15.3

# if you need to know the allocated CUDA device, you can obtain it here:
echo $CUDA_VISIBLE_DEVICES

cd /home/sliusar/microlensing
srun ./bin/gravlens input/configuration_s0.4.yaml 

