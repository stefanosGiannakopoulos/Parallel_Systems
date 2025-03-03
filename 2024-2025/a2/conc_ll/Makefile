CC = gcc
CFLAGS = -Wall -Wextra -pthread -O3

all: x.serial x.cgl x.fgl x.opt x.lazy x.nb

CFILES = main.c lib/aff.c

x.serial: $(CFILES) ll/ll_serial.c
	$(CC) $(CFLAGS) $^ -o $@
x.cgl: $(CFILES) ll/ll_cgl.c
	$(CC) $(CFLAGS) $^ -o $@
x.fgl: $(CFILES) ll/ll_fgl.c
	$(CC) $(CFLAGS) $^ -o $@
x.opt: $(CFILES) ll/ll_opt.c
	$(CC) $(CFLAGS) $^ -o $@
x.lazy: $(CFILES) ll/ll_lazy.c
	$(CC) $(CFLAGS) $^ -o $@
x.nb: $(CFILES) ll/ll_nb.c
	$(CC) $(CFLAGS) $^ -o $@

clean:
	rm -f x.*
