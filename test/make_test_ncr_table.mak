UNAME := $(shell uname)
CC = gcc

ifeq ($(UNAME), Linux)
EXT = 
LFLAGS = -lm
else
EXT = .exe
LFLAGS =
endif

CFLAGS = -std=c99 -O2 -o test_ncr_table$(EXT)

all: build

build:
	$(CC) $(CFLAGS) test_tables.c ../src/tables.c $(LFLAGS)
