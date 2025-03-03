#!/bin/bash

## Give the Job a descriptive name
#PBS -N run_kmeans

## Output and error files
#PBS -o run_kmeansseq_gpu.out
#PBS -e run_kmeansseq_gpu.err


##How long should the job run for?
#PBS -l walltime=00:10:00

## Start 
## Run make in the src folder (modify properly)

module load openmp
cd ~/2024-2025/a2/kmeans/
export OMP_NUM_THREADS=8
./kmeans_seq -s 1024 -n 2 -c 64 -l 10
