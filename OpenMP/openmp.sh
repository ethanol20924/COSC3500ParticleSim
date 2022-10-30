#!/bin/bash

#SBATCH -n 1
#SBATCH -o out/omp64mem.out
#SBATCH -e out/omp64mem.err
#SBATCH --job-name=particle_sim_omp
#SBATCH --mail-user=ethan.lo@uqconnect.edu.au
#SBATCH --mail-type=ALL

date

hostname
valgrind --leak-check=yes build/./particleSim

date
