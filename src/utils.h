#define random random_orig
#include <stdlib.h>
#undef random
#include "random_gen.h"

#ifndef STATS_UTILS_H
#define STATS_UTILS_H

#ifndef STATS_MALLOC
#define STATS_MALLOC malloc
#endif // STATS_MALLOC

// may need to handle the case of STATS_REALLOC not being defined but STATS_MALLOC is
#ifndef STATS_REALLOC
#define STATS_REALLOC realloc
#endif // STATS_MALLOC

// if user sets their own STATS_MALLOC, they should set the corresponding STATS_FREE
#ifndef STATS_FREE
#define STATS_FREE free
#endif // STATS_FREE

// quickselect partition scheme
#ifndef QUICKSELECT_PARTITION_SCHEME
#define QUICKSELECT_PARTITION_SCHEME quick_partition_lomuto_
#endif // QUICKSELECT_PARTITION_SCHEME

// quickSORT partition scheme
#ifndef QUICKSORT_PARTITION_SCHEME
#define QUICKSORT_PARTITION_SCHEME quick_partition_lomuto_
#endif // QUICKSORT_PARTITION_SCHEME

#define M_PI 3.14159265358979323846
#define IS_NEG(X) (!((X) > 0) && ((X) != 0))

enum stats_status {
    STATS_SUCCESS = 0,
    STATS_FAILURE = -1
};

typedef enum stats_status stats_status;

#define STATS_FAIL_PTR NULL;

int compare_double(const double * a, const double * b);

inline unsigned int kronecker_delta(unsigned int i, unsigned int j) {
	return (i == j) ? 1 : 0;
}

void swap(void * dest, void * src, size_t size);

// Lomuto partition scheme
size_t quick_partition_lomuto_(void * arr, size_t left, size_t right, size_t pivot_index, size_t size, int (*compar)(const void *, const void *));

// selection for kth smallest
// need to test
void * quickselect(void * arr, size_t num, size_t k, size_t size, int (*compar)(const void *, const void*));

void quicksort(void * arr, size_t num, size_t size, int (*compar)(const void *, const void *));

size_t bisect_left(void * arr, void * val, size_t left, size_t right, size_t size, int (*compar)(const void *, const void *));

size_t bisect_right(void * arr, void * val, size_t left, size_t right, size_t size, int (*compar)(const void *, const void *));

#endif // STATS_UTILS_H