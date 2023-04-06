# Particle Collision Simulation

It sucks lmao. Runs on the Getafix cluster from UQ.

Each folder contains it's own implementation. To build, go into the implementation folder and run ./build.sh. Type y/n depending on whether you want to clean or not.

To run, use the ./run.sh script. 
NOTE: IF YOU WANT TO RUN THE CUDA SCRIPT YOU MUST FIRST ASK SLURM FOR A BASH SESSION IN A GPU ENABLED NODE.

To queue a job, run:
- `sbatch particle_sim.sh` for serial
- `sbatch openmp.sh` for omp
- `sbatch cuda.sh` for cuda 
