CC = gcc

CFLAGS = -O2 -o test.exe

all: build

build:
	$(CC) $(CFLAGS) test.c ../xoshiro256plusplus.c ../random_utils.c
