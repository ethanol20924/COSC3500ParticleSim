#!/bin/bash

#SBATCH -n 1
#SBATCH -o out/getafix_serial_O3.out
#SBATCH -e out/getafix_serial.err
#SBATCH --job-name=particleSim
#SBATCH --mail-user=ethan.lo@uqconnect.edu.au
#SBATCH --mail-type=ALL

date

hostname
build/./particleSim

date
