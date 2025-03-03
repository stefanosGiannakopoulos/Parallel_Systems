#!/bin/bash

## Give the Job a descriptive name
#PBS -N run_kmeans

## Output and error files
#PBS -o run_kmeansnaive_thread_binded.out
#PBS -e run_kmeansnaive_thread_binded.err


##How long should the job run for?
#PBS -l walltime=00:10:00

## Start
## Run make in the src folder (modify properly)

module load openmp
cd ~/2024-2025/a2/kmeans/
num_threads=(1 2 4 8 16 32 64)
for nthrds in ${num_threads[@]};
do
    export OMP_NUM_THREADS=$nthrds
    export GOMP_CPU_AFFINITY="0-$(($nthrds - 1))" 

    ./kmeans_omp_naive -s 256 -n 16 -c 32 -l 10
done

