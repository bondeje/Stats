// need to rename random in stdlib.h which appears in some versions of gcc
#include <stdlib.h>
//#include "random_utils.h"

#ifndef RANDOM_GEN_H
#define RANDOM_GEN_H

// list of built-in pseudo-random number generators
#define XOSHIRO256PLUSPLUS 0

// default random number generator seed
#ifndef RANDOM_GEN_SEED
#include <time.h>
#define RANDOM_GEN_SEED time(NULL);
#endif // RANDOM_GEN_SEED

// default number generator selection
#ifndef GENERATOR
// default to the xoshiro256plusplus pseudo random number generator
#define GENERATOR XOSHIRO256PLUSPLUS
#endif // GENERATOR

// load generator and set basic functions: SPRNG = seeding function, PRNG = integer generating function, PRNG_MAX_INC = max integer output from PRNG, inclusive
#if GENERATOR == XOSHIRO256PLUSPLUS
#include "./ext/xoshiro256plusplus.h"
#define SPRNG xoshiro256plusplus_srand
#define PRNG xoshiro256plusplus_rand
#define PRNG_MAX_INC UINT64_MAX
#else
#define SPRNG srand
#define PRNG rand
#define PRNG_MAX_INC ((uint64_t)RAND_MAX)
#endif // GENERATOR

// variadice selector for randrange
#define GET_RANDRANGE_MACRO(f1, f2, f3, NAME,...) NAME
#define randrange(...) GET_RANDRANGE_MACRO(__VA_ARGS__, randrange3, randrange2, randrange1, UNUSED)(__VA_ARGS__)

/*
TODO: build a pseudo-random number generator module. Preferably use what Numpy uses: PCG XSL RR 128/64. Paper from original author: https://www.cs.hmc.edu/tr/hmc-cs-2014-0905.pdf
     Numpy's implementation: https://github.com/numpy/numpy/blob/main/numpy/random/src/pcg64/pcg64.c
     C implementation from owner. https://github.com/imneme/pcg-c
     Modify the pcg to build a pcg64 similar to the API in Numpy and use MIT license.
     Per discussion in Numpy, need to get around missing __uint128_t for Microsoft compilers: https://github.com/numpy/numpy/pull/13163

     hmmm...this analysis makes it seem like there is a problem with PCG and decorrelation or lack thereof: https://pcg.di.unimi.it/pcg.php

TODO: implement a randrange that takes up to 3 parameters
    randrange(stop)
    randrange(start, stop)
    randrange(start, stop, step)

*/

/*
random 64-bit integer in the range [0, stop). Unbiased with rejection. 
Taken from the java method based on this blog post: https://www.pcg-random.org/posts/bounded-rands.html
*/
uint64_t randrange1(uint64_t stop);
uint64_t randrange2(uint64_t start, uint64_t stop);
uint64_t randrange3(uint64_t start, uint64_t stop, uint64_t step);

/*
generate a random double in the range [0, 1)
*/
double random_();

#endif // RANDOM_GEN_H