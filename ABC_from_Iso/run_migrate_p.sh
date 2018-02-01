#!/bin/bash -l
#SBATCH -t 14400
#SBATCH -p bigmemm
#SBATCH -n 16
#SBATCH -N 4

# Ismail, the above config does the following:
# 1) Requests time
# 2) Requests a node or nodes with at least 120G of RAM
# 3) Requests the "bigmemm" partition. That is where the nodes with large memory are. You need this because the default is "high"
# 4) Request 16 cores
# 5) Request 4 nodes
# To run this script simply use this command, while in the same directory.
# sbatch cecjr-srun-ompi2.sh

# You need to load modules to use the software
module load migrate
module load slurm
module load openmpi

# Please use srun on all binaries
# Also, migrate-n-mpi looks in the current
# directory for params file
# You do not need to tell it to load params file
# if params file is in the current directory
srun migrate-n-mpi

