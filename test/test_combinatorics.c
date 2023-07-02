#include <stdio.h>
#include "../src/combinatorics.h"

void print_double_arr(double * arr, size_t n) {
    for (size_t i = 0; i < n; i++) {
        printf("%lf ", arr[i]);
    }
    printf("\n");
}

int main() {

    LexBitCombo * lbt;

    char buf[128];
    int n, k;

    printf("shift unsigned char 1 right 8: %hhu, %hhu\n", (unsigned char) 1, ((unsigned char) 1) << 8); // should be undefined by results in 1, 256 on GCC

    n = 9;
    k = 0;
    printf("building a combination (%d %d)\n", n, k);
    lbt = LexBitCombo_new(n, k);
    printf("(%d %d) n bytes %llu\n", n, k, lbt->size);
    printf("combo:\n");
    LexBitCombo_c_str(buf, lbt);
    printf("bit c_str: %s\n", buf);
    LexBitCombo_del(lbt);

    n = 6;
    k = 2;
    printf("\nbuilding a combination (%d %d)\n", n, k);
    lbt = LexBitCombo_new(n, k);
    printf("(%d %d) n bytes %llu\n", n, k, lbt->size);
    printf("combo:\n");
    LexBitCombo_c_str(buf, lbt);
    printf("bit c_str: %s\n", buf);
    LexBitCombo_del(lbt);

    n = 13;
    k = 5;
    printf("\nbuilding a combination (%d %d)\n", n, k);
    lbt = LexBitCombo_new(n, k);
    printf("(%d %d) n bytes %llu\n", n, k, lbt->size);
    printf("combo:\n");
    LexBitCombo_c_str(buf, lbt);
    printf("bit c_str: %s\n", buf);
    LexBitCombo_del(lbt);

    n = 13;
    k = 9;
    printf("\nbuilding a combination (%d %d)\n", n, k);
    lbt = LexBitCombo_new(n, k);
    printf("(%d %d) n bytes %llu\n", n, k, lbt->size);
    printf("combo:\n");
    LexBitCombo_c_str(buf, lbt);
    printf("bit c_str: %s\n", buf);
    LexBitCombo_del(lbt);

    n = 33;
    k = 17;
    printf("\nbuilding a combination (%d %d)\n", n, k);
    lbt = LexBitCombo_new(n, k);
    printf("(%d %d) n bytes %llu\n", n, k, lbt->size);
    printf("combo:\n");
    LexBitCombo_c_str(buf, lbt);
    printf("bit c_str: %s\n", buf);

    for (int i = 0; i < 40; i++) {
        LexBitCombo_next(lbt);
        LexBitCombo_c_str(buf, lbt);
        printf("bit c_str: %s\n", buf);
    }

    LexBitCombo_del(lbt);


    // test LexBitCombo_get()
    double arr[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    double comb[4] = {0.0};

    n = 6;
    k = 4;
    printf("\nbuilding a combination (%d %d)\n", n, k);
    lbt = LexBitCombo_new(n, k);
    for (int i = 0; i < 16; i++) {
        printf("%i: ", LexBitCombo_get(comb, lbt, arr, sizeof(double)));
        print_double_arr(comb, k);
    }

    LexBitCombo_del(lbt);

    printf("done.\n");

    return 0;
}