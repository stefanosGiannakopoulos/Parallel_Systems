#!/bin/bash
 
## Give the Job a descriptive name
#PBS -N run_game_of_life
 
## Output and error files
#PBS -o run_game_of_life.out
#PBS -e run_game_of_life.err
 
## How many machines should we get?
#PBS -l nodes=1:ppn=8
 
##How long should the job run for?
#PBS -l walltime=00:10:00
 
## Start
## Run make in the src folder (modify properly)
 
time_steps=1000
array_size=(64 1024 4096)
num_threads=(1 2 4 6 8)
 
module load openmp
cd /home/parallel/parlab08/2024-2025/a1
for size in ${array_size[@]};
do
    for nthrds in ${num_threads[@]};
    do
        export OMP_NUM_THREADS=$nthrds
        printf "Number of threads: %d - " "${nthrds}"
        ./Game_Of_Life $size $time_steps
    done
done
