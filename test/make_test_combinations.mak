CC = gcc

CFLAGS = -O2 -o test.exe

all: build

build:
	$(CC) $(CFLAGS) test_combinatorics.c ../src/combinatorics.c
