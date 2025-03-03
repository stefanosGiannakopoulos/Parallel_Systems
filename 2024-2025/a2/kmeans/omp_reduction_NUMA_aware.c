#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "kmeans.h"

// Square of Euclidean distance between two multi-dimensional points
inline static double euclid_dist_2(int numdims, double *coord1, double *coord2) {
    int i;
    double ans = 0.0;
    for(i = 0; i < numdims; i++)
        ans += (coord1[i] - coord2[i]) * (coord1[i] - coord2[i]);
    return ans;
}

// Finds the nearest cluster for a given object
inline static int find_nearest_cluster(int numClusters, int numCoords, double *object, double *clusters) {
    int index = 0, i;
    double dist, min_dist = euclid_dist_2(numCoords, object, clusters);

    for (i = 1; i < numClusters; i++) {
        dist = euclid_dist_2(numCoords, object, &clusters[i * numCoords]);
        if (dist < min_dist) {
            min_dist = dist;
            index = i;
        }
    }
    return index;
}

void kmeans(double *objects, int numCoords, int numObjs, int numClusters, double threshold, long loop_threshold, int *membership, double *clusters) {
    int i, j, k, index, loop = 0;
    double delta, timing = 0;
    int *newClusterSize;
    double *newClusters;
    int nthreads = omp_get_max_threads();
    printf("OpenMP Kmeans - NUMA-Aware\t(number of threads: %d)\n", nthreads);

    // Initialize membership array
    #pragma omp parallel for
    for (i = 0; i < numObjs; i++)
        membership[i] = -1;

    // Allocate global cluster data
    newClusterSize = (int *) calloc(numClusters, sizeof(int));
    newClusters = (double *) calloc(numClusters * numCoords, sizeof(double));

    // Allocate thread-local arrays with "first-touch" policy
    int **local_newClusterSize = (int **) malloc(nthreads * sizeof(int *));
    double **local_newClusters = (double **) malloc(nthreads * sizeof(double *));

    #pragma omp parallel private(i, j)
    {
        int tid = omp_get_thread_num();
        local_newClusterSize[tid] = (int *) calloc(numClusters, sizeof(int));
        local_newClusters[tid] = (double *) calloc(numClusters * numCoords, sizeof(double));

        // Initialize each thread's local clusters with "first-touch"
        #pragma omp for schedule(static)
        for (i = 0; i < numClusters; i++) {
            for (j = 0; j < numCoords; j++)
                local_newClusters[tid][i * numCoords + j] = 0.0;
            local_newClusterSize[tid][i] = 0;
        }
    }

    timing = wtime();
    do {
        delta = 0.0;

        // Reset local cluster data in a NUMA-aware way
        #pragma omp parallel private(i, j, k)
        {
            int tid = omp_get_thread_num();
            for (i = 0; i < numClusters; i++) {
                local_newClusterSize[tid][i] = 0;
                for (j = 0; j < numCoords; j++)
                    local_newClusters[tid][i * numCoords + j] = 0.0;
            }
        }

        // Assign objects to clusters
        #pragma omp parallel for reduction(+:delta) private(i, j, index)
        for (i = 0; i < numObjs; i++) {
            int tid = omp_get_thread_num();
            index = find_nearest_cluster(numClusters, numCoords, &objects[i * numCoords], clusters);

            if (membership[i] != index)
                delta += 1.0;

            membership[i] = index;
            local_newClusterSize[tid][index]++;
            for (j = 0; j < numCoords; j++)
                local_newClusters[tid][index * numCoords + j] += objects[i * numCoords + j];
        }

        // Aggregate local clusters into global clusters
        for (k = 0; k < nthreads; k++) {
            for (i = 0; i < numClusters; i++) {
                newClusterSize[i] += local_newClusterSize[k][i];
                for (j = 0; j < numCoords; j++)
                    newClusters[i * numCoords + j] += local_newClusters[k][i * numCoords + j];
            }
        }

        // Update the cluster centers by averaging
        for (i = 0; i < numClusters; i++) {
            if (newClusterSize[i] > 0) {
                for (j = 0; j < numCoords; j++)
                    clusters[i * numCoords + j] = newClusters[i * numCoords + j] / newClusterSize[i];
            }
        }

        delta /= numObjs;
        loop++;
        printf("\r\tcompleted loop %d", loop);
        fflush(stdout);
    } while (delta > threshold && loop < loop_threshold);

    timing = wtime() - timing;
    printf("\n        nloops = %3d   (total = %7.4fs)  (per loop = %7.4fs)\n", loop, timing, timing / loop);

    // Free memory for thread-local arrays
    for (k = 0; k < nthreads; k++) {
        free(local_newClusterSize[k]);
        free(local_newClusters[k]);
    }
    free(local_newClusterSize);
    free(local_newClusters);
    free(newClusters);
    free(newClusterSize);
}

