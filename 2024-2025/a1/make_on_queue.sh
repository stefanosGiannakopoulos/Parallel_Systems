#!/bin/bash

## Give the Job a descriptive name
#PBS -N make_Game_Of_Life

## Output and error files
#PBS -o make_Game_Of_Life.out
#PBS -e make_Game_Of_Life.err

## How many machines should we get?
#PBS -l nodes=1

## Start 

module load openmpi/1.8.3
cd /home/parallel/parlab08/2024-2025/a1 
make



