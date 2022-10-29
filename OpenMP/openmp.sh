#!/bin/bash

#SBATCH -n 1
#SBATCH -o out/omp64.out
#SBATCH -e out/omp64.err
#SBATCH --job-name=particle_sim_omp
#SBATCH --mail-user=ethan.lo@uqconnect.edu.au
#SBATCH --mail-type=ALL

date

hostname
build/./particleSim

date