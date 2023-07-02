#include <stdlib.h>
#include <stdint.h>
#include "random_utils.h"

unsigned int n_bits_uint64(uint64_t val) {
    unsigned int n_bits = 0;
    while (val) {
        val >>= 1;
        n_bits++;
    }
    return n_bits;
}

unsigned int n_bits_rand() {
    unsigned int bits_int = (unsigned int)sizeof(int) * 8;
    unsigned int min_bits_rand = n_bits_uint64(RAND_MAX);/*0;
    unsigned int v = RAND_MAX;
    while (v) {
        v >>= 1;
        min_bits_rand++;
    }*/
    return (min_bits_rand < bits_int) ? min_bits_rand : bits_int;
}