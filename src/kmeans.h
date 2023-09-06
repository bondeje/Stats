#ifndef KMEANS_H
#define KMEANS_H

#ifdef DEBUG_MALLOC
    #include "../../mallocs/debug_malloc.h"
#endif
#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "utils.h" // includes random_gen.h
#include "basic_statistics.h"
#include "kmeans.h"

#define DEFAULT_KMEANS_TOLERANCE 1e-5
#define DEFAULT_KMEANS_MAX_ITER 300
#define DEFAULT_KMEANS_FLAGS 0
#define DEFAULT_KMEANS_DCOMPARE_TOLERANCE 1e-14

#ifndef KM_ID_TYPE 
    #define KM_ID_TYPE unsigned short
#endif // KM_ID_TYPE

// external typedefs TODO: move to header
typedef struct KMeans KMeans;
typedef void (*fkmmean_t)(void * mean, void * el, size_t n); // this update the mean by adding el as the nth element
typedef double (*fkmdist_t)(void *, void *);
typedef int (*fkminit_t)(void *, unsigned int, void *, size_t, size_t, fkmdist_t);

// a default function for kmeans distance calculations. Arguments are pointers to double values
double kmeans_default_dist(double * a, double * b);

void kmeans_default_mean(double * out, double * el, size_t n);

// to use defaults (data are just double values), set size = 0 or size = 8 and leave fdist and fmean NULL
int KMeans_init(KMeans * km, size_t * locs, size_t n, void * data, size_t ndata, size_t size, fkmdist_t fdist, fkmmean_t fmean);

void KMeans_clean(KMeans * km);

void KMeans_reset(KMeans * km);

// delete a KMeans object created by KMeans_new
void KMeans_del(KMeans * km);

KMeans * KMeans_new_from_buffer(void * buf, size_t buf_size, size_t * locs, size_t n, void * data, size_t ndata, size_t size, fkmdist_t fdist, fkmmean_t fmean);

KMeans * KMeans_new(size_t * locs, size_t n, void * data, size_t ndata, size_t size, fkmdist_t fdist, fkmmean_t fmean);

// if user provides a null vector, returns allocation requirement for a new KMeans object
// if a KMeans object is supplied, returns the currently required additional buffer space that will be allocated in the call to KMeans
size_t KMeans_allocation_size(KMeans * km);

// passing any of these buffers in means you take responsibility for cleaning up all of them and that if they go out of scope or are free'd, no further actions are taken on the KMeans object
int KMeans_diagnostics(KMeans * km, void * means, void * last_means, KM_ID_TYPE * ids, unsigned char * workspace, size_t nbytes);

// 0s in any element other than diagnostic step mean "no change"
int KMeans_config(KMeans * km, unsigned int diagnostic_step, unsigned int max_iter, double tolerance, fkminit_t initialize, unsigned int flags);

// dist is a metric that obeys Euclidean geometry, specifically if two points are considered the same, i.e. dist(A, B) = 0 implies dist(A, C) == dist(B, C)
int kmeans_random_init(void * means, unsigned int n, void * data, size_t ndata, size_t size, fkmdist_t dist);

size_t kmeans_lloyd_assign(KMeans * km, size_t idata);

// assign cluster ids
int kmeans_set_ids(KMeans * km);

// sets new means and returns the largest change in the means. Upon initialization, output not meaningful
double kmeans_set_means(KMeans * km);

// swaps data in kmeans at locations i and j. generally do not use after kmeans() is initialized but should be ok
void kmeans_swap(KMeans * km, size_t i, size_t j);

// partitions the data into n clusters with clusters located at [0, locs[0]), [locs[0], locs[1]), etc.
// if means is null calls a generic kmeans "kmeansg" which takes a lot longer since means will be selected from the dataset
// if data is not of floating type, means MUST be NULL. The values that are nearest the true means are located at the beginning of each cluster partition
// returning -1 indicates suspended
// returning 0 means converged
// returning a positive value means failed
int kmeans(KMeans * km);

#endif // KMEANS_H