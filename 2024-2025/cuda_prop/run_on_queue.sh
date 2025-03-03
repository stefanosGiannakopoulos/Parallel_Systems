#!/bin/bash

## Give the Job a descriptive name
#PBS -N run_cuda_props

## Output and error files
#PBS -o run_cuda_props.out
#PBS -e run_cuda_props.err

## How many machines should we get? 
#PBS -l nodes=1:ppn=8

##How long should the job run for?
#PBS -l walltime=00:10:00

## Start 
## Run make in the src folder (modify properly)

cd /home/parallel/parlab08/2024-2025/cuda_prop
export CUDA_VISIBLE_DEVICES=2
./cuda_properties
