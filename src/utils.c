#include <string.h> // memcpy
#define random random_orig
#include <stdlib.h>
#undef random
#include "utils.h"
#include "random_gen.h"


int compare_double(const double * a, const double * b) {
    if (a < b) return -1;
    if (a > b) return 1;
    return 0;
}

// unbuffered swap
void swap(void * dest, void * src, size_t size)
{
	void * buf = STATS_MALLOC(size);
	memcpy(buf, src, size);
	memcpy(src, dest, size);
	memcpy(dest, buf, size);
	STATS_FREE(buf);
}

// Lomuto partition scheme
size_t quick_partition_lomuto_(void * arr, size_t left, size_t right, size_t pivot_index, size_t size, int (*compar)(const void *, const void *)) {
	
	unsigned char * carr = (unsigned char *) arr;
	unsigned char * pivot = carr + pivot_index * size;
	
	swap((void *) (carr + right * size), (void *)pivot, size); // arr[right] is now pivot value
	pivot = carr + right * size;
	right--;
	
	while (left < right) {
		if (compar((void *) (carr + left * size), pivot) < 0) {
			left++;
		} else {
			swap((void *) (carr + left * size), (void *) (carr + right * size), size);
			right--;
		}
	}
	
	if (compar((void *) (carr + left * size), pivot) < 0) {
		left++;
	}
	
	swap((void *) (carr + left * size), (void *) pivot, size);
		
	return left;
}

// selection for kth smallest
// need to test
void * quickselect(void * arr, size_t num, size_t k, size_t size, int (*compar)(const void *, const void*)) {
	if (num < k || !arr || !size || !compar) {
		return NULL;
	}

	size_t ind, left = 0, right = num-1;
	while ((ind = QUICKSELECT_PARTITION_SCHEME(arr, left, right, randrange(left, right+1), size, compar)) != k-1) {
		if (ind < k-1) {
			left = ind+1;
		} else {
			right = ind-1;
		}
		DEBUG_PRINT(("ind = %zu, left = %zu, right = %zu\n", ind, left, right));
	}
	
	return (void *) (((unsigned char *) arr) + ind * size);
}

void quicksort_(void * arr, size_t left, size_t right, size_t size, int (*compar)(const void *, const void *)) {
	if (left < right) {
		size_t ind = QUICKSORT_PARTITION_SCHEME(arr, left, right, randrange(left, right+1), size, compar);
		if (ind > left) {
			quicksort_(arr, left, ind - 1, size, compar);
		}
		if (ind < right) {
			quicksort_(arr, ind + 1, right, size, compar);
		}
	}
	return;
}

void quicksort(void * arr, size_t num, size_t size, int (*compar)(const void *, const void *)) {
	quicksort_(arr, 0, num-1, size, compar);
	return;
}

size_t bisect_left(void * arr, void * val, size_t left, size_t right, size_t size, int (*compar)(const void *, const void *)) {
	if (!arr || !val || !size || !compar || left > right) {
		return - 1; // TODO: does this make sense for an error
	}
	unsigned char * carr = (unsigned char *) arr;
	if (compar((void *) (carr + left * size), val) >= 0) {
		return left;
	} else if (compar((void *) (carr + right * size), val) < 0) {
		return right + 1;
	}
	
	// goal is to have arr[right] >= val & arr[left] < val at all times...arr[left]==val is dealt with at the beginning as an edge case
	// because of the condition above, left + 1 == right will result in infinite loop if allowed
	while (left + 1 < right) {
		size_t mid = left + (right - left)/2;
		if (compar((void *) (carr + mid * size), val) < 0) { // arr[mid] < val
			left = mid;
		} else { // arr[mid] >= val
			right = mid;
		}
	}
	
	// since arr[left] < val for all iterations and arr[right] >= val, the final right is the correct insertion point
	return right;
}

size_t bisect_right(void * arr, void * val, size_t left, size_t right, size_t size, int (*compar)(const void *, const void *)) {
	if (!arr || !val || !size || !compar || left > right) {
		return - 1; // TODO: does this make sense for an error
	}
	unsigned char * carr = (unsigned char *) arr;
	if (compar((void *) (carr + left * size), val) > 0) {
		return left;
	} else if (compar((void *) (carr + right * size), val) <= 0) {
		return right + 1;
	}
	
	// goal is to have arr[right] > val & arr[left] <= val at all times...arr[right]==val is dealt with at the beginning as an edge case
	// because of the condition above, left + 1 == right will result in infinite loop if allowed
	while (left + 1 < right) {
		size_t mid = left + (right - left)/2;
		if (compar((void *) (carr + mid * size), val) <= 0) { // arr[mid] <= val
			left = mid;
		} else { // arr[mid] > val
			right = mid;
		}
	}
	
	// since arr[left] < val for all iterations and arr[right] >= val, the final right is the correct insertion point
	return right;
}