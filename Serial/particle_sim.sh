#!/bin/bash

#SBATCH -n 1
#SBATCH -o out/serial_.out
#SBATCH -e out/serial_.err
#SBATCH --job-name=particleSim
#SBATCH --mail-user=ethan.lo@uqconnect.edu.au
#SBATCH --mail-type=ALL
#SBATCH -t 0-02:00:00

#SBATCH --mem-per-cpu=32G

date

hostname
build/./particleSim

date
