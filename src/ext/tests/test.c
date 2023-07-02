#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "../xoshiro256plusplus.h"

int main() {
    xoshiro256plusplus_srand(time(NULL));

    printf("random variables:\n");
    for (int i = 0; i < 50; i++) {
        printf("%llu\n", xoshiro256plusplus_rand());
    }

    return 0;
}