#!/bin/bash

## Give the Job a descriptive name
#PBS -N testjob

## Output and error files
#PBS -o testjob.out
#PBS -e testjob.err

## Limit memory, runtime etc.
#PBS -l walltime=01:00:00

## How many nodes:processors_per_node should we get?
#PBS -l nodes=2:ppn=4

## Start 
##echo "PBS_NODEFILE = $PBS_NODEFILE"
##cat $PBS_NODEFILE

## Run the job (use full paths to make sure we execute the correct thing) 
## NOTE: Fix the path to show to your executable! 

mkdir -p ${HOME}/tmp
export TMPDIR=${HOME}/tmp
module load openmpi/1.8.3
mpirun -np 4 --map-by node --mca btl self,tcp \
/home/users/goumas/benchmarks/MPI_code/fw/fw_MPI 32
rm -r ${HOME}/tmp
