#include <string.h> // memcpy
#include <stdio.h> // delete when not printing
#include <stddef.h>
#include <stdbool.h>
#include "utils.h"
#include "combinatorics.h"

/*
TODO: probably a better idea to implement a void * to bit string and generalize LexBitCombo_c_str
*/

#define MAX_BIT 0b10000000

LexBitCombo * LexBitCombo_new(size_t n, size_t k) {
    LexBitCombo * lbt;
    if (!(lbt = (LexBitCombo *) STATS_MALLOC(sizeof(LexBitCombo)))) {
        return STATS_FAIL_PTR;
    }

    lbt->size = (n+1)/BYTE_SIZE;
    if (n % BYTE_SIZE) {
        lbt->size++;
    }

    if (!(lbt->combo = (unsigned char *) STATS_MALLOC(sizeof(UBYTE)*(lbt->size)))) {
        STATS_FREE(lbt);
        lbt = NULL;
        return STATS_FAIL_PTR;
    }

    LexBitCombo_init(lbt, n, k);

    return lbt;
}

void LexBitCombo_init(LexBitCombo * lbt, size_t n, size_t k) {
    for (size_t i = 0; i < lbt->size; i++) {
        lbt->combo[i] = 0;
    }
    
    lbt->n = n;
    lbt->k = k;
    lbt->_rightmost1_bit = 0;
    // set lbt->n'th bit to 1 to mark the end of the combination
    size_t byte = n / BYTE_SIZE;
    UBYTE bit = 1 << (n % BYTE_SIZE);
    lbt->combo[byte] |= bit;

    if (k==0) { // all elements are 0 as alread initialized
        lbt->stop_iteration = true;
        return;
    }
    lbt->stop_iteration = false;

    // initialize to onces in the least significant bits
    byte = 0;
    bit = 1;
    lbt->combo[byte] |= bit;
    
    lbt->_rightmost1_bit++;
    while (lbt->_rightmost1_bit < k) {
        if (!(lbt->_rightmost1_bit % BYTE_SIZE)) {
            bit = 1; // reset to 1 and increase byte
            byte++;
        } else {
            bit <<= 1; // shift left 1 but stay in same byte
        }
        lbt->combo[byte] |= bit;
        
        lbt->_rightmost1_bit++;
    }
    lbt->_rightmost1_bit--;
}

void LexBitCombo_del(LexBitCombo * lbt) {
    STATS_FREE(lbt->combo);
    lbt->combo = NULL;
    STATS_FREE(lbt);
}

static inline void LexBitCombo_left_shift_(size_t * byte, UBYTE * bit, size_t dist) {
    for (size_t i = 0; i < dist; i++) {
        if (*bit & MAX_BIT) {
            *bit = 1;
            (*byte)++;
        } else {
            (*bit) <<= 1;
        }
    }
}

static inline void LexBitCombo_right_shift_(size_t * byte, UBYTE * bit, size_t dist) {
    size_t i = 0;
    while ((*byte || (*bit > 1)) && (i < dist)) {
        i++;
        if (*bit & 1) {
            *bit = MAX_BIT;
            (*byte)--;
        } else {
            (*bit) >>= 1;
        }
    }
    if (i < dist) {
        *bit = 0;
    }
}

// buf must be allocated to >lbt->n + 1 in length. output is number of bytes written 
int LexBitCombo_c_str(char * buf, LexBitCombo * lbt) {
    size_t byte = 0;
    UBYTE bit = 1;
    size_t j = lbt->n;
    buf[j--] = '\0';

    while (j > 0) {
        buf[j] = (lbt->combo[byte] & bit) ? '1' : '0';
        LexBitCombo_left_shift_(&byte, &bit, 1);
        j--;
    }
    buf[j] = (lbt->combo[byte] & bit) ? '1' : '0';
    return lbt->n+1;
}

inline unsigned short LexBitCombo_get_bit(LexBitCombo * lbt, size_t bitindex) {
    return lbt->combo[bitindex / BYTE_SIZE] & (1 << (bitindex % BYTE_SIZE));
}

// CONSIDER: change output to an iterator status variable?
stats_status LexBitCombo_next(LexBitCombo * lbt) {   
    size_t el = lbt->k;
    if (!el) {
        lbt->stop_iteration = true;
        return STATS_FAILURE;
    }
    size_t byte = lbt->_rightmost1_bit / BYTE_SIZE;
    UBYTE bit = (1 << (lbt->_rightmost1_bit % BYTE_SIZE));
    size_t next_byte = ((lbt->_rightmost1_bit + 1) / BYTE_SIZE);
    UBYTE next_bit = (1 << ((lbt->_rightmost1_bit + 1) % BYTE_SIZE));

    while ((lbt->combo[next_byte] & next_bit)) { // if i am on a valid element that cannot bitshift left
        el--;
        
        if (!el) { // cannot move last element
            break;
        }
        LexBitCombo_right_shift_(&byte, &bit, 1);
        lbt->_rightmost1_bit--;

        while (!(lbt->combo[byte] & bit)) { // should also check that the condition byte == 0 and bit & 1 does not fail with el > 0
            LexBitCombo_right_shift_(&byte, &bit, 1);
            lbt->_rightmost1_bit--;
        }
        // calculate next_byte and next_bit
        next_bit = bit;
        next_byte = byte;
        LexBitCombo_left_shift_(&next_byte, &next_bit, 1);
    }
    // need to handle the case of exiting on el == 0, which should be a failure
    if (!el) {
        lbt->stop_iteration = true;
        return STATS_FAILURE;
    }
    lbt->combo[byte] ^= bit; // set bit to 0
    lbt->combo[next_byte] ^= next_bit; // set next bit to 1
    lbt->_rightmost1_bit++;

    // CONSIDER: there has to be a much more efficient way to do the remainder of this

    // need to reset all bits above this bit to 1. The number of bits to add above this one is k-el
    while (el < lbt->k) {
        LexBitCombo_left_shift_(&next_byte, &next_bit, 1);
        lbt->combo[next_byte] |= next_bit;
        lbt->_rightmost1_bit++;
        el++;
    }
    
    el = lbt->_rightmost1_bit+1;
    
    while (el++ < lbt->n) {
        LexBitCombo_left_shift_(&next_byte, &next_bit, 1);
        lbt->combo[next_byte] &= ~next_bit; // turn off the bit
    }    
       
    return STATS_SUCCESS;
}

// size is the size in bytes of each element in src
// src must be allocated to lbt->n*size bytes
// dest must be allocated to lbt->k*size bytes
stats_status LexBitCombo_get(void * dest, LexBitCombo * lbt, void * src, size_t size) {
    if (!lbt->k || lbt->stop_iteration) {
        return STATS_FAILURE;
    }

    unsigned char * dest_arr = (unsigned char *) dest;
    unsigned char * src_arr = (unsigned char *) src;
    size_t byte = 0;
    UBYTE bit = 1;
    for (size_t i = 0; i < lbt->k; i++) {
        // find the next 1
        while (!(lbt->combo[byte] & bit)) { // there must be k 1s in the array so that this should only have to find the first 1 and not have to check if I've reached the end
            LexBitCombo_left_shift_(&byte, &bit, 1);
            src_arr += size;
        }

        // copy the data to the destination
        memcpy(dest_arr, src_arr, size);
        dest_arr += size;

        // advance the bit
        LexBitCombo_left_shift_(&byte, &bit, 1);
        src_arr += size;
    }

    if (!lbt->stop_iteration) {
        LexBitCombo_next(lbt);
    }
    return STATS_SUCCESS;
}

