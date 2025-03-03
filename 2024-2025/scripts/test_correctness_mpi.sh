#!/bin/bash

## Give the Job a descriptive name
#PBS -N testjob

## Output and error files
#PBS -o test_correctness_mpi.out
#PBS -e test_correctness_mpi.err

## Limit memory, runtime etc.

## How many nodes:processors_per_node should we get?
## Run on parlab
#PBS -l nodes=4:ppn=8

## Start 
##echo "PBS_NODEFILE = $PBS_NODEFILE"
##cat $PBS_NODEFILE

## Run the job (use full paths to make sure we execute the correct thing)

 
## NOTE: Fix the path to show to your serial executables 
module load openmpi/1.8.3
cd /path/to/serial/programs/

for execfile in jacobi seidelsor 
do
	./${execfile} 256 256
done


## NOTE: Fix the path to show to your MPI executables
cd /path/to/mpi/programs/

for execfile in jacobi seidelsor  
do
	mpirun -np 32 --map-by node --mca btl self,tcp ${execfile} 256 256 8 4
done

## Make sure you enable convergence testing and printing

