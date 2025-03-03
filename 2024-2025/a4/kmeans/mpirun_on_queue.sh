#!/bin/bash
#PBS -N kmeans_mpi_paralleled
#PBS -o kmeans_mpi_paralleled.out
#PBS -e kmeans_mpi_paralleled.err
#PBS -l walltime=01:00:00
#PBS -l nodes=8:ppn=8

export TMPDIR=${HOME}/mpi_sessions
mkdir -p $TMPDIR

module load openmpi/1.8.3

SIZE=256
COORDS=16
CLUSTERS=32
LOOPS=10

NP_LIST="1 2 4 8 16 32 64"

for NP in $NP_LIST; do
    echo "Running with ${NP} MPI processes"
    mpirun -np ${NP} --map-by node --mca btl self,tcp \
    ~/2024-2025/a4/kmeans/kmeans_mpi -s ${SIZE} -n ${COORDS} -c ${CLUSTERS} -l ${LOOPS}
done

rm -rf $TMPDIR

