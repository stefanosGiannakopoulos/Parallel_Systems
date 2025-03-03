.KEEP_STATE:

MPICC = mpicc

CFLAGS = -Wall -Wextra -Wno-unused -O3

LDFLAGS = -lm

H_FILES = utils.h

COMM_SRC = utils.c

all: mpi_seidelsor_conv

mpi_seidelsor_conv: mpi_seidelsor_conv.o utils.o
	$(MPICC) $(CFLAGS) $^ -o $@ $(LDFLAGS)  # <---- TAB before $(MPICC)

mpi_seidelsor_conv.o: mpi_gauss_seidel_conv.c
	$(MPICC) $(CFLAGS) -c $< -o $@  # <---- TAB before $(MPICC)

utils.o: utils.c
	$(MPICC) $(CFLAGS) -c $< -o $@  # <---- TAB before $(MPICC)

clean: 
	rm -rf *.o mpi_seidelsor_conv  # <---- TAB before rm

