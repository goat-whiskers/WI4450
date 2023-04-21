#!/bin/bash

#SBATCH --exclusive
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 48
#SBATCH --mem=1G
#SBATCH --account=Education-EEMCS-Courses-WI4450

# run CG for a n^3 Poisson problem
n= 12
export OMP_NUM_THREADS=48
# srun ./main_cg_poisson.x $n
srun ./main_benchmarks.x
