.KEEP_STATE:

CC = gcc

CFLAGS = -Wall -Wextra -Wno-unused -O3
OMPFLAGS = -fopenmp $(CFLAGS)
LDFLAGS =

H_FILES = kmeans.h
COMM_SRC = file_io.c util.c

all: kmeans_omp_reduction kmeans_seq kmeans_omp_naive kmeans_omp_reduction_with_false_sharing_handled_NUMA_optimised

kmeans_seq: main.o file_io.o util.o seq_kmeans.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

kmeans_omp_naive: main.o file_io.o util.o omp_naive_kmeans.o
	$(CC) $(OMPFLAGS) $^ -o $@ $(LDFLAGS)

kmeans_omp_reduction_with_false_sharing_handled_NUMA_optimised: main.o file_io.o util.o omp_reduction_NUMA_aware.o
	$(CC) $(OMPFLAGS) $^ -o $@ $(LDFLAGS)

kmeans_omp_reduction: main.o file_io.o util.o omp_reduction.o
	$(CC) $(OMPFLAGS) $^ -o $@ $(LDFLAGS)

main.o: main.c $(H_FILES)
	$(CC) $(CFLAGS) -c $< -o $@

seq_kmeans.o: seq_kmeans.c $(COMM_SRC) $(H_FILES)
	$(CC) $(CFLAGS) -c $< -o $@

omp_naive_kmeans.o: omp_naive_kmeans.c $(COMM_SRC) $(H_FILES)
	$(CC) $(OMPFLAGS) -c $< -o $@

omp_reduction_NUMA_aware.o: omp_reduction_NUMA_aware.c $(COMM_SRC) $(H_FILES)
	$(CC) $(OMPFLAGS) -c $< -o $@

omp_reduction.o: omp_reduction_kmeans.c $(COMM_SRC) $(H_FILES)
	$(CC) $(OMPFLAGS) -c $< -o $@

file_io.o: file_io.c
	$(CC) $(OMPFLAGS) -c $< -o $@

util.o: util.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf *.o kmeans_seq kmeans_omp_naive kmeans_omp_reduction_with_false_sharing_handled_NUMA_optimised

