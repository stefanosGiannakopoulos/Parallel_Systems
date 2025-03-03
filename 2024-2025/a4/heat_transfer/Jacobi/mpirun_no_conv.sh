#!/bin/bash
#PBS -N jacobi_mpi_no_conv
#PBS -o jacobi_mpi_no_conv.out
#PBS -e jacobi_mpi_no_conv.err
#PBS -l walltime=01:00:00
#PBS -l nodes=8:ppn=8

export TMPDIR=${HOME}/mpi_sessions
mkdir -p $TMPDIR

module load openmpi/1.8.3

cd /home/parallel/parlab08/2024-2025/a4/heat_transfer/Jacobi/

iterations=(1 2 3)
sizes=(2048 4096 6144)
threads=(1 2 4 8 16 32 64)
threads1=(1 1 2 2 4 4 8)
threads2=(1 2 2 4 4 8 8)

#for it in "${iterations[@]}"; do
    for size in "${sizes[@]}"; do
        for i in "${!threads[@]}"; do
            echo "Iteration: -, Size: $size, Threads: ${threads[$i]}"
            mpirun --mca btl self,tcp -np "${threads[$i]}" -map-by node \
            ./mpi_jacobi_no_conv "$size" "$size" "${threads1[$i]}" "${threads2[$i]}"
        done
    done
#done

rm -rf $TMPDIR


