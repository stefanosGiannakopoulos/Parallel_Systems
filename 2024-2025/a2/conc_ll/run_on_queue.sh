#!/bin/bash

## Give the Job a descriptive name
#PBS -N run_conc_ll_parallel

## Output and error files
#PBS -o run_conc_ll_parallel.out
#PBS -e run_conc_ll_parallel.err

## How long should the job run for?
#PBS -l walltime=01:05:00

## Start
## Run make in the src folder (modify properly)

cd ~/2024-2025/a2/conc_ll/

num_threads=(1 2 4 8 16 32 64 128)
max_cores=32
executables=("x.cgl" "x.fgl" "x.lazy" "x.nb" "x.opt")
list_sizes=(1024 8192)
percentages=("100 0 0" "80 10 10" "20 40 40" "0 50 50")

output_file="results/all_results.out"
mkdir -p "$(dirname "$output_file")"

# Write a header to the output file
if [ ! -f "$output_file" ]; then
    echo "Experiment Results" > "$output_file"
    echo "===================" >> "$output_file"
fi

generate_mt_conf() {
    local nthrds=$1

    if [ $nthrds -le 64 ]; then
        # Sequential assignment for up to 64 threads
        echo $(seq -s, 0 $((nthrds - 1)))
    elif [ $nthrds -eq 128 ]; then
        # Repeat 0 to 63 for 128 threads
        local aux=""
        for ((i=0; i<64; i++)); do
            aux+="$i,$i,"
        done
        aux=${aux%,} # Remove trailing comma
        echo "$aux"
    else
        # Handle other cases (if necessary)
        echo "Unsupported thread count: $nthrds" >&2
        exit 1
    fi
}

for exec in "${executables[@]}"; do
    for lsize in ${list_sizes[@]}; do
        for perc in "${percentages[@]}"; do
            read -r contains_pct add_pct remove_pct <<< "$perc"

            for nthrds in ${num_threads[@]}; do
                MT_CONF=$(generate_mt_conf "$nthrds")
                export MT_CONF

                workload_name="${contains_pct}-${add_pct}-${remove_pct}"

                echo "Running $exec with Threads=$nthrds, List Size=$lsize, Workload=$perc, MT_CONF=$MT_CONF" >> "$output_file"

                throughput=$("./$exec" "$lsize" "$contains_pct" "$add_pct" "$remove_pct")
                if [ $? -ne 0 ]; then
                    echo "Error running $exec with Threads=$nthrds, List Size=$lsize, Workload=$perc" >> "$output_file"
                    echo "-------------------" >> "$output_file"
                    continue
                fi

                # Append results to the single output file
                echo "Executable: $exec" >> "$output_file"
                echo "Threads: $nthrds" >> "$output_file"
                echo "List Size: $lsize" >> "$output_file"
                echo "Workload: $workload_name" >> "$output_file"
                echo "Throughput: $throughput" >> "$output_file"
                echo "-------------------" >> "$output_file"
            done
        done
    done
done

echo "All experiments completed. Results written to $output_file."

