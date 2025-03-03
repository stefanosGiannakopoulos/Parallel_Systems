#!/bin/bash

## Give the Job a descriptive namec
#PBS -N run_conc_ll_parallel

## Output and error files
#PBS -o conc_ll_parallel.out
#PBS -e conc_ll_parallel.err

##How long should the job run for?
#PBS -l walltime=00:15:00

## Start
## Run make in the src folder (modify properly)

cd ~/2024-2025/a2/conc_ll/

num_threads=(1 2 4 8 16 32 64 128)
executables="opt.x" #("x.cgl" "x.fgl" "x.lazy" "x.nb" "x.opt")
list_size=(1024 8192)
percentages=("100 0 0" "80 10 10" "20 40 40" "0 50 50")

for exec in "${executables[@]}"; do
    for lsize in ${list_size[@]}; do
        for perc in "${percentages[@]}"; do
            for nthrds in ${num_threads[@]};
            do
                if [ $((nthrds)) == 128 ]; then
                    aux=$(seq -s, 0 63)
                    aux+=','
                    aux+=$aux
                    export MT_CONF=$aux
                else
                    export MT_CONF=$(seq -s, 0 $((nthrds - 1)))
                fi
                ./$exec $lsize $perc
            done
        done
    done
done
