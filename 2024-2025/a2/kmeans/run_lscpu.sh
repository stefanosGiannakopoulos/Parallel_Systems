#!/bin/bash

## Give the Job a descriptive name
#PBS -N run_lscpu

## Output and error files
#PBS -o run_lscpu.out
#PBS -e run_lscpu.err


##How long should the job run for?
#PBS -l walltime=00:05:00

## Start

lscpu

