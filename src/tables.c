#include <stddef.h> // NULL
#include <math.h>
#include "tables.h"
#include "utils.h"

/*
CONSIDER: implement it as a hash table global variable 
once I have the hash table working on arrays, I can hash on (n,k), which will 
be significantly easier, but will require a significant modification

CONSIDER: storing only (n, k) for k < n/2.0 since Pascal's triangle is symmetric about k = n/2.0

CONSIDER: moving Pascal's triangle to a separate header/.c file
*/

static size_t * pascals_triangle = NULL;
static size_t pascals_triangle_size = 0;

#define PASCALS_TRIANGLE_INIT_SIZE 36

// Forward declarations of internal functions
static size_t pascals_triangle_map_nk_(size_t n, size_t k); // TODO: remove from header
static void pascals_triangle_map_index_(size_t * n, size_t * k, size_t index); // TODO: remove from header
static void pascals_triangle_populate_(size_t start, size_t end); // TODO: remove from header

void pascals_triangle_init(size_t size) {
    if (!size) {
        size = PASCALS_TRIANGLE_INIT_SIZE;
    }
    if (size <= pascals_triangle_size) {
        return;
    }
    pascals_triangle = (size_t *) STATS_MALLOC(sizeof(size_t) * size);
    pascals_triangle_size = size;
    pascals_triangle_populate_(0, size-1);
}

// get index from n,k
/*
n,k -> index
0,0 -> 0
1,0 -> 1
1,1 -> 2
2,0 -> 3
2,1 -> 4
2,2 -> 5
3,0 -> 6
3,1 -> 7
3,2 -> 8
3,3 -> 9
4,0 -> 10
.
.
.
7,7 -> 35
*/
size_t pascals_triangle_map_nk_(size_t n, size_t k) {
    size_t index;
    if (n & 1) {
        index = (n+1)/2*n + k;
    } else {
        index = n/2*(n+1) + k;
    }
    
    return index;
}

// get n,k from index
void pascals_triangle_map_index_(size_t * n, size_t * k, size_t index) {

    /* index satisfies:
    index = n^2 / 2 + n / 2 + k
    the correct value of n satisfies
    n >= n_low : n_low^2 + n_low - 2*index = 0
    n_low = -1 / 2 + sqrt(1 + 8*index) / 2

    n <= n_high : n_high^2 + 3 * n_high - 2*index = 0
    n_high = -3 / 2 + sqrt(9 + 8*index) / 2

    (-1 + sqrt(1 + 8*index))/2 <= n <= (-3 + sqrt(9 + 8*index))/2
    */

    *n = (size_t) floor((sqrt(1 + 8*index) - 1)/2);
    size_t test_index = pascals_triangle_map_nk_(*n, 0);
    if (test_index > index) {
        (*n)--;
        test_index = pascals_triangle_map_nk_(*n, 0);
    } else if (index - test_index > *n) {
        (*n)++;
        test_index = pascals_triangle_map_nk_(*n, 0);
    }

    *k = index - test_index;
}

static void * pascals_triangle_extend_(size_t size) {
    size_t cur_size = pascals_triangle_size;
    if (cur_size >= size) {
        return pascals_triangle;
    } else if (size < 2*cur_size) {
        size = 2*cur_size;
    }

    size_t * new_triangle = (size_t *) STATS_REALLOC(pascals_triangle, sizeof(size_t) * size);
    if (new_triangle) {
        pascals_triangle = new_triangle;
        pascals_triangle_size = size;
        pascals_triangle_populate_(cur_size, size-1);
    }
    return new_triangle;
}

// 0 means failure since output is size_t
size_t pascals_triangle_get(size_t n, size_t k) {
    size_t index = pascals_triangle_map_nk_(n,k);
    if (index >= pascals_triangle_size) {
        if (!pascals_triangle_extend_(index+1)) {
            return 0;
        }
    }
    return pascals_triangle[index];
}

size_t (*nCr) (size_t, size_t) = pascals_triangle_get;

// start and end are indices, inclusive
void pascals_triangle_populate_(size_t start, size_t end) {
    if (!start) { // 0 is a special value
        pascals_triangle[start++] = 1;
    }
    size_t n;
    size_t k;
    pascals_triangle_map_index_(&n, &k, start);

    while (start <= end) {
        size_t val = 0;
        if (k < n) {
            val += pascals_triangle_get(n-1, k); // this actually can fail if pascal's triangle is not big enough and cannot allocate enough memory
        }
        if (k) {
            val += pascals_triangle_get(n-1, k-1);
        }

        pascals_triangle[start] = val;

        start++;
        if (k == n) {
            n++;
            k = 0;
        } else {
            k++;
        }
    }
}

void pascals_triangle_del() {
    STATS_FREE(pascals_triangle);
    pascals_triangle = NULL;
    pascals_triangle_size = 0;
}

void tables_init() {
    pascals_triangle_init(0);
}

void tables_del() {
    pascals_triangle_del();
}