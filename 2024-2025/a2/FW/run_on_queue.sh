#!/bin/bash

## Give the Job a descriptive name
#PBS -N run_fw


#PBS -o run_fw_tiled_mod.out
#PBS -e run_fw.err

## How many machines should we get? 

##How long should the job run for?
#PBS -l walltime=01:00:00

## Start 
## Run make in the src folder (modify properly)

module load openmp
cd /home/parallel/parlab08/2024-2025/a2/FW

## Define matrix sizes, thread counts, and block sizes
SIZES=("1024" "2048" "4096")
THREAD_COUNTS=(1 2 4 8 16 32 64)
BSIZES=(64 128 512)

## Loop over sizes, thread counts, and block sizes
for SIZE in "${SIZES[@]}"; do
  for THREADS in "${THREAD_COUNTS[@]}"; do
    for BSIZE in "${BSIZES[@]}"; do
      export OMP_NUM_THREADS=$THREADS
      ## Set CPU affinity to bind threads to specific CPUs
      export GOMP_CPU_AFFINITY="0-$(($THREADS - 1))"     
 echo "Running fw_sr with SIZE=$SIZE, BSIZE=$BSIZE, THREADS=$THREADS"
      ./fw_tiled_mod "$SIZE" "$BSIZE"
    done
  done
done


