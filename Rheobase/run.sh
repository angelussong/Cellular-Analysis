#!/bin/sh
#SBATCH -A msm110
#SBATCH --job-name=""
#SBATCH --output="%j.%N.out"
#SBATCH --partition=compute
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=16
#SBATCH -t 02:30:00
#SBATCH --no-requeue

nrnivmodl optmz model

mpirun -np 256 nrniv -mpi init.hoc

