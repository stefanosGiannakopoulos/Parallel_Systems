#!/bin/bash

## Give the Job a descriptive name
#PBS -N run_kmeans

## Output and error files
#PBS -o run_kmeans_without_locks.out
#PBS -e run_kmeans_without_locks.err

## How long should the job run for?
#PBS -l walltime=00:59:00

## Start
## Run make in the src folder (modify properly)

module load openmp
cd ~/2024-2025/a2/kmeans_locks/
num_threads=(1 2 4 8 16 32 64)

executables=("kmeans_omp_naive" "kmeans_omp_critical")

for exe in "${executables[@]}"; do
  for nthrds in ${num_threads[@]}; do
    export OMP_NUM_THREADS=$nthrds
    export GOMP_CPU_AFFINITY="0-$(($nthrds - 1))"

    ./$exe -s 32 -n 16 -c 32 -l 10
  done
done

