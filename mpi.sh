#!/bin/bash
#SBATCH --job-name=your_job_name # Job name
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=20                   # number of processes = 20
#SBATCH --cpus-per-task=1      # Number of CPU cores allocated to each process (please use 1 here, in comparison with pthread)
#SBATCH --partition=Project            # Partition name: Project or Debug (Debug is default)

# cd /nfsmnt/119010355/CSC4005_2022Fall_Demo/project3_template/
clear
make clean
make mpig
mpirun -np 4 ./mpig 12 5000

# clear
# make clean
# make mpi
# mpirun -np 4 ./mpi 12 10

# mpirun -np 20 ./mpi 1000 100
# mpirun -np 40 ./mpi 1000 100