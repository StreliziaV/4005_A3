#!/bin/bash
#SBATCH --job-name=mpi_openMP # Job name
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=20                   # number of processes = 20
#SBATCH --cpus-per-task=1      # Number of CPU cores allocated to each process (please use 1 here, in comparison with pthread)
#SBATCH --partition=Project            # Partition name: Project or Debug (Debug is default)

# cd cd /nfsmnt/119010369/4005_A3/
clear
make clean
make mpip
mpirun -np 4 ./mpip 4 12 5000