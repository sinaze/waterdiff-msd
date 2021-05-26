CC=gcc
CFLAGS=-Wall
DEPS = xdrfile_xtc.h
OBJ = main.o
LIBS=-lm -lxdrfile

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

locmsd: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)
