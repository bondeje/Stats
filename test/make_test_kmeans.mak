UNAME := $(shell uname)
CC = gcc
DEBUG ?= false
DEBUG_FILE = 

IFLAGS = -I../src -I../../mallocs
CFLAGS = 
LFLAGS =
EXT = 
ifeq ($(UNAME), Linux)
LFLAGS += -lm
else
EXT += .exe
ifeq ($(DEBUG), true)
	DEBUG_FILE += ../../mallocs/debug_malloc.c
	CFLAGS += -DDEBUG_MALLOC -DDEBUG_MALLOC_VERBOSE
endif

endif

CFLAGS += -std=c99 -O3 -o test_kmeans$(EXT)

all: build

build:
	$(CC) $(CFLAGS) ../src/kmeans.c ../src/utils.c ../src/random_gen.c ../src/ext/xoshiro256plusplus.c ../src/ext/random_utils.c $(DEBUG_FILE) test_kmeans.c $(LFLAGS)
