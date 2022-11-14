#!/bin/bash

#SBATCH -n 1
#SBATCH -o out/cuda_chonk.out
#SBATCH -e out/cuda_chonk.err
#SBATCH --job-name=particle_sim_cuda
#SBATCH --mail-user=ethan.lo@uqconnect.edu.au
#SBATCH --mail-type=ALL

#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu
#SBATCH --constraint=A100

date

hostname

module load gnu cuda

nvidia-smi

build/./particleSim

date