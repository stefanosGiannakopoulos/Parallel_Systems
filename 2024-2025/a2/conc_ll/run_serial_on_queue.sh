#!/bin/bash

## Give the Job a descriptive name
#PBS -N run_concurrent_ll

## Output and error files
#PBS -o run_concurrent_ll.out
#PBS -e run_concurrent_ll.err

## How long should the job run for?
#PBS -l walltime=01:05:00

## Load required modules
cd ~/2024-2025/a2/conc_ll/

## Configuration
EXECUTABLES=("x.serial")
THREADS=(1)
LIST_SIZES=(1024 8192)
WORKLOADS=("100-0-0" "80-10-10" "20-40-40" "0-50-50")

## Function to parse workload percentages
parse_workload() {
    IFS='-' read -r contains_pct add_pct remove_pct <<< "$1"
    echo "$contains_pct $add_pct $remove_pct"
}

## Run experiments
for executable in "${EXECUTABLES[@]}"; do
    OUTPUT_FILE="${executable}.out"

    echo "Results for $executable" > $OUTPUT_FILE
    echo "Threads, List Size, Workload, Throughput" >> $OUTPUT_FILE

    for size in "${LIST_SIZES[@]}"; do
        for workload in "${WORKLOADS[@]}"; do
            # Parse workload percentages
            read -r contains_pct add_pct remove_pct <<< "$(parse_workload "$workload")"

            for threads in "${THREADS[@]}"; do
                MT_CONF=$(seq -s, 0 $((threads - 1)))

                export MT_CONF=$MT_CONF

                echo "Running $executable with Threads=$threads, List Size=$size, Workload=$workload, MT_CONF=$MT_CONF"

                # Execute and capture throughput
                THROUGHPUT=$("./$executable" "$size" "$contains_pct" "$add_pct" "$remove_pct")

                # Log results
                echo "$threads, $size, $workload, $THROUGHPUT" >> $OUTPUT_FILE
            done
        done
    done

    echo "Results for $executable saved in $OUTPUT_FILE"
done

echo "All experiments completed."
