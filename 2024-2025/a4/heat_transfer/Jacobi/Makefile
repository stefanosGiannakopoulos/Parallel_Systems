.KEEP_STATE:

MPICC = mpicc

CFLAGS = -Wall -Wextra -Wno-unused -O3

LDFLAGS = 

H_FILES = utils.h

COMM_SRC = utils.c

all: mpi_jacobi_no_conv

mpi_jacobi_no_conv: mpi_jacobi_no_conv.o utils.o
	$(MPICC) $(CFLAGS) $^ -o $@ $(LDFLAGS)  # <---- TAB before $(MPICC)

mpi_jacobi_no_conv.o: Jacobi_MPI_no_conv.c
	$(MPICC) $(CFLAGS) -c $< -o $@  # <---- TAB before $(MPICC)

utils.o: utils.c
	$(MPICC) $(CFLAGS) -c $< -o $@  # <---- TAB before $(MPICC)

clean: 
	rm -rf *.o mpi_jacobi_no_conv  # <---- TAB before rm

