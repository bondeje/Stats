/*
Adaptor for the xoshiro256plusplus.c source code and my random_gen.h
*/

#include <stdint.h>

#define xoshiro256plusplus_rand xoshiro256plusplus_next

void xoshiro256plusplus_srand(unsigned int seed);

uint64_t xoshiro256plusplus_next(void);

void xoshiro256plusplus_jump(void);

void xoshiro256plusplus_long_jump(void);