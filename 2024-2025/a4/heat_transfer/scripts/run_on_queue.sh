#!/bin/bash

## Give the Job a descriptive name
#PBS -N run_mpi_helloworld

## Output and error files
#PBS -o run_mpi_helloworld.out
#PBS -e run_mpi_helloworld.err

## How many machines should we get? 
#PBS -l nodes=2:ppn=8

## Start 
## Run make in the src folder (modify properly)

module load openmpi/1.8.3
cd /home/parallel/pps/2024-2025/lab_guide/mpi
mpirun -np 16 ./mpi_helloworld

