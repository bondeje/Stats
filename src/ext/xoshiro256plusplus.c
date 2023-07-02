#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include "random_utils.h"
#include "xoshiro256plusplus.h"

/* 
below code retrieved 2023/02/10 from https://prng.di.unimi.it/xoshiro256plusplus.c by Jeffrey Bonde with modifications per the following list 

all function names and the seed name prepended with "xoshiro256plusplus_"
Since the seed is static and not in the header, did not necessarily have to prepend a namespace

*/

/*  Written in 2019 by David Blackman and Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

#include <stdint.h>

/* This is xoshiro256++ 1.0, one of our all-purpose, rock-solid generators.
   It has excellent (sub-ns) speed, a state (256 bits) that is large
   enough for any parallel application, and it passes all tests we are
   aware of.

   For generating just floating-point numbers, xoshiro256+ is even faster.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s. */

static inline uint64_t xoshiro256plusplus_rotl(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}


static uint64_t xoshiro256plusplus_s[4];

uint64_t xoshiro256plusplus_next(void) {
	const uint64_t result = xoshiro256plusplus_rotl(xoshiro256plusplus_s[0] + xoshiro256plusplus_s[3], 23) + xoshiro256plusplus_s[0];

	const uint64_t t = xoshiro256plusplus_s[1] << 17;

	xoshiro256plusplus_s[2] ^= xoshiro256plusplus_s[0];
	xoshiro256plusplus_s[3] ^= xoshiro256plusplus_s[1];
	xoshiro256plusplus_s[1] ^= xoshiro256plusplus_s[2];
	xoshiro256plusplus_s[0] ^= xoshiro256plusplus_s[3];

	xoshiro256plusplus_s[2] ^= t;

	xoshiro256plusplus_s[3] = xoshiro256plusplus_rotl(xoshiro256plusplus_s[3], 45);

	return result;
}


/* This is the jump function for the generator. It is equivalent
   to 2^128 calls to next(); it can be used to generate 2^128
   non-overlapping subsequences for parallel computations. */

void xoshiro256plusplus_jump(void) {
	static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };

	uint64_t s0 = 0;
	uint64_t s1 = 0;
	uint64_t s2 = 0;
	uint64_t s3 = 0;
	for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
		for(int b = 0; b < 64; b++) {
			if (JUMP[i] & UINT64_C(1) << b) {
				s0 ^= xoshiro256plusplus_s[0];
				s1 ^= xoshiro256plusplus_s[1];
				s2 ^= xoshiro256plusplus_s[2];
				s3 ^= xoshiro256plusplus_s[3];
			}
			xoshiro256plusplus_next();	
		}
		
	xoshiro256plusplus_s[0] = s0;
	xoshiro256plusplus_s[1] = s1;
	xoshiro256plusplus_s[2] = s2;
	xoshiro256plusplus_s[3] = s3;
}



/* This is the long-jump function for the generator. It is equivalent to
   2^192 calls to next(); it can be used to generate 2^64 starting points,
   from each of which jump() will generate 2^64 non-overlapping
   subsequences for parallel distributed computations. */

void xoshiro256plusplus_long_jump(void) {
	static const uint64_t LONG_JUMP[] = { 0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241, 0x39109bb02acbe635 };

	uint64_t s0 = 0;
	uint64_t s1 = 0;
	uint64_t s2 = 0;
	uint64_t s3 = 0;
	for(int i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
		for(int b = 0; b < 64; b++) {
			if (LONG_JUMP[i] & UINT64_C(1) << b) {
				s0 ^= xoshiro256plusplus_s[0];
				s1 ^= xoshiro256plusplus_s[1];
				s2 ^= xoshiro256plusplus_s[2];
				s3 ^= xoshiro256plusplus_s[3];
			}
			xoshiro256plusplus_next();	
		}
		
	xoshiro256plusplus_s[0] = s0;
	xoshiro256plusplus_s[1] = s1;
	xoshiro256plusplus_s[2] = s2;
	xoshiro256plusplus_s[3] = s3;
}

/********************************** seeding code *****************************/

// Added by Jeffrey Bonde to seed the xoshiro256plusplus prng

// since random_gen sets what SRAND will be called, presumably, either this srand xor stdlib's srand will be called and only once
// also, this way xoshiro256plusplus_srand inherits the reset behavior of stdlib's srand
void xoshiro256plusplus_srand(unsigned int seed) {
    srand(seed);

    uint64_t val = 0;
    unsigned int max_ext_bits = n_bits_rand(), bits_rem, ext_bits;

    for (int j = 0; j < 4; j++) {
        for (int i = 0; i < 8; i++) {
            val = 0;
            bits_rem = sizeof(uint64_t) * 8;
            while (bits_rem) {
                ext_bits = (max_ext_bits > bits_rem) ? bits_rem : max_ext_bits;
                val |= (((uint64_t) rand()) << (bits_rem - ext_bits));
                bits_rem -= ext_bits;
            }
            xoshiro256plusplus_s[j] = val;
        }
    }

    //printf("xoshiro seeded with:\n%llu\n%llu\n%llu\n%llu\n", s[0], s[1], s[2], s[3]);
}