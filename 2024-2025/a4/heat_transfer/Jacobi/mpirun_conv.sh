#!/bin/bash
#PBS -N jacobi_mpi_conv
#PBS -o jacobi_mpi_conv.out
#PBS -e jacobi_mpi_conv.err
#PBS -l walltime=01:00:00
#PBS -l nodes=8:ppn=8

export TMPDIR=${HOME}/mpi_sessions
mkdir -p $TMPDIR

module load openmpi/1.8.3

iterations=(1 2 3)

for it in "${iterations[@]}"; do
    echo "Iteration ${it} - 64 MPI processes"
    mpirun -np 64 --map-by node --mca btl self,tcp \
    /home/parallel/parlab08/2024-2025/a4/heat_transfer/Jacobi/mpi_jacobi 512 512 8 8
done

rm -rf $TMPDIR


