#include <stdbool.h>
#include "utils.h"

#ifndef COMBINATORICS_H
#define COMBINATORICS_H

/*
TODO: probably a better idea to implement a void * to bit string and generalize LexBitCombo_c_str
*/

typedef unsigned char UBYTE;
#define BYTE_SIZE (sizeof(UBYTE)*8)
// set MAX_BIT appropriately in combinatorics.c if bits per element changes for UBYTE

// for internal use. Should not be directly instantiated
typedef struct LexBitCombo {
    UBYTE * combo; // char array long enough to store n bits
    size_t size; // size in bytes or number of elements in combo
    size_t n;
    size_t k;
    size_t _rightmost1_bit;
    bool stop_iteration; // next call to LexBitCombo_next will fail
} LexBitCombo;

typedef void Combo; // this is needed for the way I create Iterators on the fly by name

LexBitCombo * LexBitCombo_new(size_t n, size_t k);
void LexBitCombo_init(LexBitCombo * lbt, size_t n, size_t k);
void LexBitCombo_del(LexBitCombo * lbt);
int LexBitCombo_c_str(char * buf, LexBitCombo * lbt);
stats_status LexBitCombo_next(LexBitCombo * lbt);
stats_status LexBitCombo_get(void * dest, LexBitCombo * lbt, void * src, size_t size);

#endif // COMBINATORICS_H