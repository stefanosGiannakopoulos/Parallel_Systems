.PHONY: all clean

all: fw fw_sr fw_tiled fw_recursive fw_sr_no_changes fw_tiled_mod

CC = gcc
CFLAGS = -Wall -O3 -Wno-unused-variable -fopenmp

HDEPS += %.h

OBJS = util.o

fw: $(OBJS) fw.c
	$(CC) $(OBJS) fw.c -o fw $(CFLAGS)

fw_sr: fw_sr.c
	$(CC) $(OBJS) fw_sr.c -o fw_sr $(CFLAGS)

fw_sr_no_changes: fw_sr_no_changes.c
	$(CC) $(OBJS) fw_sr_no_changes.c -o fw_sr_no_changes $(CFLAGS)

fw_tiled: fw_tiled.c
	$(CC) $(OBJS) fw_tiled.c -o fw_tiled $(CFLAGS)

fw_tiled_mod: fw_tiled_mod.c
	$(CC) $(OBJS) fw_tiled_mod.c -o fw_tiled_mod $(CFLAGS)

fw_recursive: $(OBJS) fw_recursive.c
	$(CC) $(OBJS) fw_recursive.c -o fw_recursive $(CFLAGS)

%.o: %.c $(HDEPS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o fw fw_sr fw_tiled fw_recursive fw_sr_no_changes fw_tiled_mod

