#!/bin/bash

## Give the Job a descriptive name
#PBS -N testjob

## Output and error files
#PBS -o test_correctness_mpi_no_conv.out
#PBS -e test_correctness_mpi_no_conv.err

## Limit memory, runtime etc.
#PBS -l walltime=00:15:00

## How many nodes:processors_per_node should we get?
## Run on parlab
#PBS -l nodes=4:ppn=8

## Start 
##echo "PBS_NODEFILE = $PBS_NODEFILE"
##cat $PBS_NODEFILE

## Run the job (use full paths to make sure we execute the correct thing)

 
# NOTE: Fix the path to show to your serial executables 
module load openmpi/1.8.3
cd /home/parallel/parlab08/2024-2025/a4/heat_transfer/serial

for execfile in redblacksor
do
	./${execfile} 256 256
done


## NOTE: Fix the path to show to your MPI executables
cd ~/2024-2025/a4/heat_transfer/RedBlack

for execfile in mpi_redblack_conv  
do
	mpirun -np 32 --map-by node --mca btl self,tcp ${execfile} 256 256 8 4
done

## Make sure you enable convergence testing and printing

