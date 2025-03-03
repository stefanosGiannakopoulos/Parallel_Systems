# Parallel Systems - University Course Exercises

## Overview
This document details my work on solving all the exercises assigned for the *Parallel Systems* course during the academic year 2024-25 at the National Technical University of Athens (NTUA). The exercises covered different aspects of parallel computing, including shared memory architectures, GPU programming, and distributed memory systems.

---

## 1. Familiarization with the Programming Environment
### *Conway's Game of Life*
- Implemented a parallel version using OpenMP.
- Conducted performance measurements on different core counts and grid sizes.
- Analyzed and compared execution times.
- Explored optimizations and visualized results.

---

## 2. Parallelization and Optimization on Shared Memory Architectures
### *K-means Clustering Algorithm*
- Developed two parallel implementations:
  1. *Shared clusters* (using OpenMP synchronization primitives).
  2. *Copied clusters and reduction* (avoiding synchronization overhead).
- Analyzed performance scaling with thread binding (GOMP_CPU_AFFINITY).
- Addressed false-sharing issues and explored NUMA-aware allocation.

### *Lock Implementations*
- Implemented and benchmarked different locking mechanisms, including:
  - `pthread_mutex_lock`, `pthread_spin_lock`, `tas_lock`, `ttas_lock`, `array_lock`, `clh_lock`.
- Compared performance against OpenMP `critical` and `atomic` directives.

### *Floyd-Warshall Algorithm*
- Parallelized the recursive version using OpenMP tasks.
- Analyzed execution time scaling with matrix sizes.
- (Bonus) Implemented and compared tiled vs. recursive versions.

### *Concurrent Data Structures*
- Evaluated different synchronization mechanisms for linked lists:
  - Coarse-grain locking, fine-grain locking, optimistic synchronization, lazy synchronization, non-blocking synchronization.
- Conducted extensive benchmarking with different thread counts and operation ratios.

---

## 3. Parallelization on GPU Architectures
### *K-means Clustering with CUDA*
- Implemented multiple GPU-accelerated versions:
  1. *Naive version* (basic GPU kernel for nearest-cluster computation).
  2. *Transpose version* (optimized memory layout for better coalescing).
  3. *Shared-memory version* (utilizing GPU shared memory for improved performance).
  4. *Full-offload version* (executing entire algorithm on GPU).
  5. *Delta-reduction version* (using efficient reduction operations instead of atomic updates).
- Conducted detailed performance analysis and profiling.

---

## 4. Parallelization on Distributed Memory Architectures
### *K-means Clustering with MPI*
- Implemented a distributed-memory version using MPI.
- Compared execution time and speedup against the OpenMP version.

### *2D Heat Diffusion Simulation*
- Parallelized three numerical methods using MPI:
  1. *Jacobi method*
  2. *Gauss-Seidel SOR method*
  3. *Red-Black SOR method*
- Measured performance with and without convergence checking.
- Analyzed scaling behavior across different MPI process counts.

---

## Conclusion
This coursework involved hands-on experience with parallel programming across different architectures, including shared memory (OpenMP), GPU (CUDA), and distributed memory (MPI). The optimizations and performance evaluations provided deep insights into parallel computing challenges and best practices.