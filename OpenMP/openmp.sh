#!/bin/bash

#SBATCH -n 1
#SBATCH -o out/omp1_.out
#SBATCH -e out/omp1_.err
#SBATCH --job-name=particle_sim_omp
#SBATCH --mail-user=ethan.lo@uqconnect.edu.au
#SBATCH --mail-type=ALL
#SBATCH -t 0-02:00:00

#SBATCH --mem=32G
#SBATCH --cpus-per-task=32
#SBATCH --constraint=R640

export OMP_NUM_THREADS=1

date

hostname
build/./particleSim
# valgrind --leak-check=yes build/./particleSim

date
