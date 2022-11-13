#!/bin/bash

#SBATCH -n 1
#SBATCH -o out/omp64_.out
#SBATCH -e out/omp64_.err
#SBATCH --job-name=particle_sim_omp
#SBATCH --mail-user=ethan.lo@uqconnect.edu.au
#SBATCH --mail-type=ALL
#SBATCH -t 0-02:00:00

#SBATCH --mem-per-cpu=32G

date

hostname
build/./particleSim
# valgrind --leak-check=yes build/./particleSim

date
