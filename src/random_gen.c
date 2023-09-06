// add a shuffle feature like python

#include <stdint.h>
#include <stdio.h> // for debugging only
#include "random_gen.h"

/*
random 64-bit integer in the range [0, stop). Unbiased with rejection. 
Taken from the java method based on this blog post: https://www.pcg-random.org/posts/bounded-rands.html
*/
uint64_t randrange1(uint64_t stop) {
    if (!stop) {
        return 0;
    }
    uint64_t x, r;
    do {
        x = PRNG();
        r = x % stop;
    } while (x - r > (-stop));
    return r;
}

uint64_t randrange2(uint64_t start, uint64_t stop) {
    return start + randrange1(stop-start);
}

uint64_t randrange3(uint64_t start, uint64_t stop, uint64_t step) {
    return start + step * randrange1((stop - start + step - 1) / step);
}

/*
generate a random double in the range [0, 1)
*/
double random() {
    return randrange(PRNG_MAX_INC) / (((double) PRNG_MAX_INC));
}

/*
generate a random double in the range [start, end)
unlike Python's random.uniform, this guarantees end is exclusive by rejection
*/
double uniform(double start, double end) {
    double ret = start + (end-start) * random();
    while (ret >= end) {
        ret = start + (end-start) * random();
    }
    return ret;
}

// probably should be inlined
void * random_select(void * parr, size_t narr, size_t size) {
    printf("randomly selecting from [0, %zu): ", narr);
    size_t val = randrange1((uint64_t) narr);
    printf("%zu\noffsetting %p by %zu (size = %zu)\n", val, parr, val * size, size);
    return (void *) (((unsigned char *)parr) + val * size);
}

// this modifies the input array in place
// need to buffer this
// TODO: need to fix dependencies on utils.h, which depends on random_gen.h
int random_selectn(void * parr, size_t nselect, size_t narr, size_t size) {
    if (!parr || nselect > narr) {
        return 1;
    }
    if (nselect == narr) {
        return 0;
    }
    unsigned char * parrc = (unsigned char *) parrc;
    uint64_t start = 0, stop = (uint64_t) narr;
    while (nselect - start) {
        size_t ind = (size_t) randrange2(start, stop);
        if (ind-start) { // swap the values at ind and start
            //swap(parrc + ind * size, parrc + start * size, size);
        }
        start++;
    }
    return 0;
}

// TODO: redo the remainder of the function signatures to be consistent with other random selects and implements
int random_selectn_unique(void * parr, size_t nselect, size_t narr, int (*comp)(void *, void *)) {
    return 0;
}

// select nselect random values without repeats from [0, stop)
int random_selectu(uint64_t * selected, size_t nselect, uint64_t stop) {
    return 0;
}

// algorithm for sorted parr where the order must be maintained
int random_select_sorted(void ** selected, size_t nselect, void ** parr, size_t narr) {
    return 0;
}