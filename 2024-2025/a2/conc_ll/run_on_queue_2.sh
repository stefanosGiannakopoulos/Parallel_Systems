#!/bin/bash

## Give the Job a descriptive namec
#PBS -N run_conc_ll_parallel

## Output and error files
#PBS -o run_conc_ll_serial.out
#PBS -e run_conc_ll_serial.err

##How long should the job run for?
#PBS -l walltime=01:05:00

## Start
## Run make in the src folder (modify properly)

cd ~/2024-2025/a2/conc_ll/

num_threads=(1)
max_cores=32
executables=("x.serial")
list_sizes=(1024 8192)
percentages=("100 0 0" "80 10 10" "20 40 40" "0 50 50")

output_dir="results"
mkdir -p "$output_dir"

generate_mt_conf() {
    local nthrds=$1
    local max_cores=$2

    if [ $nthrds -le $max_cores ]; then
        echo $(seq -s, 0 $((nthrds - 1)))
    else
        local threads_per_core=$((nthrds / max_cores))
        local aux=""
        
        for ((i=0; i<max_cores; i++)); do
            for ((j=0; j<threads_per_core; j++)); do
                aux+="$i,"
            done
        done
        
        local remaining=$((nthrds % max_cores))
        for ((k=0; k<remaining; k++)); do
            aux+="$k,"
        done
        
        aux=${aux%,}
        echo $aux
    fi
}

# Initialize the summary file
summary_file="${output_dir}/summary.out"
echo "Executable, List Size, Threads, Workload, Throughput" > "$summary_file"

for exec in "${executables[@]}"; do
    for lsize in ${list_sizes[@]}; do
        for perc in "${percentages[@]}"; do
            read -r contains_pct add_pct remove_pct <<< "$perc"

            for nthrds in ${num_threads[@]}; do
                MT_CONF=$(generate_mt_conf "$nthrds" "$max_cores")
                export MT_CONF

                workload_name="${contains_pct}-${add_pct}-${remove_pct}"
                output_file="${output_dir}/${exec}_size${lsize}_threads${nthrds}_workload${workload_name}.out"

                echo "Results for $exec, List Size=$lsize, Threads=$nthrds, Workload=$workload_name" > "$output_file"
                echo "Threads, List Size, Workload, Throughput" >> "$output_file"

                echo "Running $exec with Threads=$nthrds, List Size=$lsize, Workload=$perc, MT_CONF=$MT_CONF"

                throughput=$("./$exec" "$lsize" "$contains_pct" "$add_pct" "$remove_pct")
                if [ $? -ne 0 ]; then
                    echo "Error running $exec with Threads=$nthrds, List Size=$lsize, Workload=$perc" >> "$output_file"
                    continue
                fi

                echo "$nthrds, $lsize, $workload_name, $throughput" >> "$output_file"
                
                # Append throughput to the summary file
                echo "$exec, $lsize, $nthrds, $workload_name, $throughput" >> "$summary_file"
            done
        done
    done
done

echo "All experiments completed. Summary saved to $summary_file."

