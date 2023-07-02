#include <stdio.h>
#include "../src/tables.h"

void print_row(unsigned int n) {
    for (int k = 0; k < n+1; k++) {
        printf("%zu ", nCr(n, k));
    }
    printf("\n");
}


int main() {

    pascals_triangle_init(0);

    for (unsigned int n = 0; n < 12; n++) {
        print_row(n);
    }

    pascals_triangle_del(); // should delete

    return 0;
}