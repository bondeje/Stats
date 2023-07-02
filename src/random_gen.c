#include <stdint.h>
#include "random_gen.h"

/*
random 64-bit integer in the range [0, stop). Unbiased with rejection. 
Taken from the java method based on this blog post: https://www.pcg-random.org/posts/bounded-rands.html
*/
uint64_t _randrange1(uint64_t stop) {
    uint64_t x, r;
    do {
        x = PRNG();
        r = x % stop;
    } while (x - r > (-stop));
    return r;
}

uint64_t _randrange2(uint64_t start, uint64_t stop) {
    return start + _randrange1(stop-start);
}

uint64_t _randrange3(uint64_t start, uint64_t stop, uint64_t step) {
    return start + step * _randrange1((stop - start + step - 1) / step);
}

/*
generate a random double in the range [0, 1)
*/
double _random() {
    return randrange(PRNG_MAX_INC) / (((double) PRNG_MAX_INC));
}

/*
generate a random double in the range [start, end)
unlike Python's random.uniform, this guarantees end is exclusive by rejection
*/
double uniform(double start, double end) {
    double ret = start + (end-start) * _random();
    while (ret >= end) {
        ret = start + (end-start) * _random();
    }
    return ret;
}