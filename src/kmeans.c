// TODO:
// add a prediction feature. if ids are buffered, this is trivial, if they are not, have to recalculate each mean for each cluster every time 
// add kmeans++ initialization option
// KMeans->cluster_ids should be unsigned int...if not unsigned short, not size_t
// remove km->el element. It should be unnecessary since I made changed from utils\swapbuf to utils\swap and made it bufferless

#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include <stdio.h> // for debugging only
#include <assert.h>
#include <string.h>
#include "utils.h" // includes random_gen.h
#include "basic_statistics.h"
#include "kmeans.h"

// flags for KMeans
#define KMEANS_MEANS_BUFFERED 0b1
#define KMEANS_LAST_MEANS_BUFFERED 0b10
#define KMEANS_EL_BUFFERED 0b100
#define KMEANS_IDS_BUFFERED 0b1000
#define KMEANS_LOCS_BUFFERED 0b10000
#define KMEANS_BUFFERED 0b100000
#define KMEANS_BUFFERED_MASK 0b111111 // should cover all above cases

#define KMEANS_ALGO_MASK 0b1000000 // right now 0 is Lloyd's algorithm, do not have others

#define OFFSET_VOID_PTR(ptr, index, size) (void *) (((unsigned char *)ptr) + index * size)
#define INC_VOID_PTR(ptr, index, size) ptr = OFFSET_VOID_PTR(ptr, index, size)

// one configuration for each data type, which includes 
// function for averaging of the datatype, 
// function for the distance between elements of the data type
// the size of the data type in bytes, defaults to 
// opaque data type
// I = Input
// O = Output
// D = diagnostic output if set to pause
// int = internal only
// with () = optional, without () = mandatory
struct KMeans {
    unsigned char * means; // (I)/(O), means output buffer
    unsigned char * last_means; // storage for previous means
    unsigned char * data;  // I data buffer, should be size * ndata bytes
    unsigned char * el; // (I) element buffer, should be at least 'size' bytes
    unsigned char * mbuf; // pointer to the beginning of the internally malloc'd memory, if any
    unsigned char * message; // user-provided workspace address. NULL if not provided, will not allocate separately
    KM_ID_TYPE * cluster_ids;  // (I) id buffer, should be ndata long if input
    size_t * locs; // I // locations
    size_t * locs_buf; // internal buffer for locations
    fkmdist_t fdist; // I
    fkmmean_t fmean; // I
    fkminit_t initialize; // (I)
    double tolerance; // (I)
    size_t ndata; // I
    size_t size; // I
    size_t mbuf_size; // number of bytes that will or are allocated to mbuf
    size_t message_size; // size of user-provided workspace
    KM_ID_TYPE n; // I
    unsigned int max_iter; // (I)
    unsigned int iter; // int
    unsigned int iter_step; // (DI), internal flag if diagnostics are set that certain parameter cannot be updated
    unsigned int flags; // int
};

// a default function for kmeans distance calculations. Arguments are pointers to double values
double kmeans_default_dist(double * a, double * b) {
    return fabs(*a - *b);
}

void kmeans_default_mean(double * out, double * el, size_t n) {
    *out = (*out * (n-1) + *el)/n;
}

/* TODO: this is supposed to be an efficient method to initialize unique means, but it is probably only necessary for km->n very large since the current method of kmeans_check_unique is O(n^2)
// return number of points overlapping, puts them at the end of vals if sort is true
size_t kmeans_check_overlap(void ** vals, unsigned int n, double (*dist)(void *, void *)) {
    // algorithm for overlap of points. 
    // An INefficient method is to check all n * (n - 1) / 2 pairs
    // A slightly more efficient method is to 
    // (1) calculate the distance from the first point to all other points: O(n)
    //   (1a) if any distance is 0.0, increment output by 1, move it to the end of vals if sort is true
    // as long as the dist metric satisfies Euclidean geometry, then a necessary requirement for two points to be the same is they are the same distance to the first point
    // (2) so, take the distances from (1) in a hash-set, only check the pairs that have matched distances to the first point.
    // (3) for each pair found in (2), check distance between them
    //   (3a) if distance is 0.0, increment output by 1, move it to the end of vals if sort is true
    return 0;
}
*/

// slower method to check unique
// while building means, currently at size n, checks if el is significantly different from all of the current means
int kmeans_check_unique(unsigned char * means, size_t n, unsigned char * el, size_t size, fkmdist_t dist) {
    if (!means || !n || !size) { // current no elements to compare against, so el is uniuqe
        return 1;
    }
    if (!el || !dist) { // no el and no distance metric, so el is not necessarily unique
        return 0;
    }
    int result = 1;
    while (result && n--) {
        result = (dist((void *) (means + n * size), el) > DEFAULT_KMEANS_DCOMPARE_TOLERANCE) ? 1 : 0;
    }
    return result;
}

// dist is a metric that obeys Euclidean geometry, specifically if two points are considered the same, i.e. dist(A, B) = 0 implies dist(A, C) == dist(B, C)
int kmeans_random_init(void * means, unsigned int n, void * data, size_t ndata, size_t size, fkmdist_t dist) {
    printf("in kmeans_random_init\n");
    //printf("size = %zu\n", size);
    if (!data || !size || !dist || !ndata) {
        return 1;
    }
    if (!n) {
        return 0;
    }

    // slower method to ensure the initial means are different. Basically builds them one at a
    size_t start = 0;
    int status = 0;
    while (!status && start < n) {
        // TODO: switch to random_selectn when the dependency clash is fixed
        //status = random_selectn(data, 1, ndata - start, size);
        //if (status) {
        //    break;
        //}
        void * dest = OFFSET_VOID_PTR(means, start, size);
        //printf("attemting to write a piece of data into means (%p) at location (%p), offset %zu\n", means, dest, ((unsigned char *)(dest) - (unsigned char *)(means)) / sizeof(void *));
        // TODO: remove selection variable. unnecessary
        void * selection = random_select(data, ndata, size);
        //printf("attemting to write a piece of data from data (%p) at location (%p), offset %zu\n", data, selection, ((unsigned char *)(selection) - (unsigned char *)(data)) / sizeof(void *));
        memcpy(dest, selection, size);
        // if 
        //printf("writing complete, checking if means are unique\n");
        if (kmeans_check_unique((unsigned char *)means, start, (unsigned char *) dest, size, dist)) {
            start += 1;
        }
    }
    return status;
}

int KMeans_alloc_buffers(KMeans * km) {
    size_t nbytes = km->mbuf_size;
    int status = 0;
    unsigned char * buffer;
    if (nbytes) {
        buffer = (unsigned char *) STATS_MALLOC(nbytes);
        printf("allocated %zu bytes for buffers at %p\n", nbytes, (void *) buffer);
        if (!buffer) {
            return 1;
        }
        km->mbuf = buffer;
        
        size_t nbytes_needed = km->n * sizeof(size_t);
        if (!km->locs_buf) {
            if (nbytes < nbytes_needed) {
                printf("insufficient bytes allocated for km->locs_buf: have (%zu) need (%zu)\n", nbytes, nbytes_needed);
                status = 2;
                goto failalloc;
            }
            printf("set km->locs_buf from allocation %p\n", (void *) buffer);
            km->locs_buf = (size_t *) buffer;
            buffer += nbytes_needed;
            nbytes -= nbytes_needed;
        }

        nbytes_needed = km->ndata * sizeof(KM_ID_TYPE);
        if (!km->cluster_ids) {
            if (nbytes < nbytes_needed) {
                printf("insufficient bytes allocated for km->cluster_ids: have (%zu) need (%zu)\n", nbytes, nbytes_needed);
                status = 3;
                goto failalloc;
            }
            printf("set km->cluster_ids from allocation %p\n", (void *) buffer);
            km->cluster_ids = (KM_ID_TYPE *) buffer;
            buffer += nbytes_needed;
            nbytes -= nbytes_needed;
        }
        
        nbytes_needed = km->size;
        if (!km->el) {
            if (nbytes < nbytes_needed) {
                printf("insufficient bytes allocated for km->el: have (%zu) need (%zu)\n", nbytes, nbytes_needed);
                status = 4;
                goto failalloc;
            }
            printf("set km->el from allocation %p\n", (void *) buffer);
            km->el = buffer;
            buffer += nbytes_needed;
            nbytes -= nbytes_needed;
        }

        nbytes_needed = km->n * km->size;
        if (!km->means) {
            if (nbytes < nbytes_needed) {
                printf("insufficient bytes allocated for km->means: have (%zu) need (%zu)\n", nbytes, nbytes_needed);
                status = 5;
                goto failalloc;
            }
            printf("set km->means from allocation %p\n", (void *) buffer);
            km->means = buffer;
            buffer += nbytes_needed;
            nbytes -= nbytes_needed;
        }

        if (!km->last_means) {
            if (nbytes < nbytes_needed) {
                printf("insufficient bytes allocated for km->last_means: have (%zu) need (%zu)\n", nbytes, nbytes_needed);
                status = 6;
                goto failalloc;
            }
            printf("set km->last_means from allocation %p\n", (void *) buffer);
            km->last_means = buffer;
            buffer += nbytes_needed;
            nbytes -= nbytes_needed;
        }
    }
    return status;
failalloc:
    STATS_FREE(buffer);
    return status;
}

// to use defaults (data are just double values), set size = 0 or size = 8 and leave fdist and fmean NULL
// these are the mostly required data to be able to run while the last 3 do accept simple defaults
// setting n or ndata to anything other than what was set with KMeans_new causes a hard reset and all previously set buffers from user will be unset (internally allocated) without an additional call to KMeans_diagnostics
int KMeans_init(KMeans * km, size_t * locs, size_t n, void * data, size_t ndata, size_t size, fkmdist_t fdist, fkmmean_t fmean) {
    if (!km || !locs || !data || !n || !ndata) {
        return 1;
    }
    if (size == 0) {
        size = sizeof(double);
    }
    if (size == sizeof(double)) {
        if (!fdist) {
            fdist = (fkmdist_t) kmeans_default_dist;
        }
        if (!fmean) {
            fmean = (fkmmean_t) kmeans_default_mean;
        }
    } else {
        if (!fdist) {
            return 2;
        }
        if (!fmean) {
            return 3;
        }
    }
    bool hard_reset = false;
    if (n != km->n || ndata != km->ndata) { // buffer size requirements have changed. KMeans_new must set these before call KMeans_init
        hard_reset = true;
        km->ndata = ndata;
        km->n = n;
    }
    km->data = (unsigned char *) data;
    km->locs = locs; // user provided
    km->fdist = fdist;
    km->fmean = fmean;
    km->initialize = (fkminit_t) kmeans_random_init; // can be overriden with KMeans_initialize
    km->tolerance = DEFAULT_KMEANS_TOLERANCE; // tolerance for idenifying if means converged
    km->size = size;
    km->max_iter = DEFAULT_KMEANS_MAX_ITER;
    km->iter = 0;
    km->iter_step = 0; // can be set with KMeans_enable_diagnostics(specifying a step size to pause) kmeans does not pause if set to 0
    //km->flags = 0; // identify algorithms and cleanup when kmeans converges. 0 indicates means, el, cluster_ids must all free'd. retain user set flags
    km->mbuf_size = 0;
    km->mbuf = NULL;
    if (hard_reset || !(km->flags & KMEANS_MEANS_BUFFERED)) { // retain buffering from user if previously set
        km->means = NULL; // will be allocated at start of kmeans
        km->mbuf_size += n * size; // ****
    }
    if (hard_reset || !(km->flags & KMEANS_LAST_MEANS_BUFFERED)) { // retain buffering from user if previously set
        km->last_means = NULL;
        km->mbuf_size += n * size; // ****
    }
    if (hard_reset || !(km->flags & KMEANS_EL_BUFFERED)) { // retain buffering from user if previously set
        km->el = NULL; // will be allocated at start of kmeans KMeans_el
        km->mbuf_size += size; // ****
    }
    if (hard_reset || !(km->flags & KMEANS_IDS_BUFFERED)) { // retain buffering from user if previously set
        km->cluster_ids = NULL; // will be allocated at start of kmeans unless provided by KMeans_ids
        km->mbuf_size += ndata * sizeof(KM_ID_TYPE); // ****
    }
    if (hard_reset || !(km->flags & KMEANS_EL_BUFFERED)) { // retain buffering from user if previously set
        km->locs_buf = NULL;
        km->mbuf_size += n * sizeof(size_t); // ****
    }
    printf("completed initiailization, need %zu of additional buffer space\n", km->mbuf_size);
    return 0;
}

void KMeans_clean(KMeans * km) {
    if (km->mbuf) {
        STATS_FREE(km->mbuf);
        km->mbuf = NULL;
    }
    if (!(km->flags & KMEANS_MEANS_BUFFERED)) { // km->means was provided by km->mbuf
        km->means = NULL;
    }
    if (!(km->flags & KMEANS_LAST_MEANS_BUFFERED)) { // km->last_means was provided by km->mbuf
        km->last_means = NULL;
    }
    if (!(km->flags & KMEANS_IDS_BUFFERED)) { // km->cluster_ids was provided by km->mbuf
        km->cluster_ids = NULL;
    }
    if (!(km->flags & KMEANS_EL_BUFFERED)) { // km->el was provided by km->mbuf
        km->el = NULL;
    }
    if (!(km->flags & KMEANS_LOCS_BUFFERED)) { // km->locs was provided by km->mbuf
        km->locs_buf = NULL;
    }
    //km->iter = 0;
    //km->iter_step = 0;
}

// if user provides a null vector, returns allocation requirement for a new KMeans object
// if a KMeans object is supplied, returns the currently required additional buffer space that will be allocated in the call to KMeans
size_t KMeans_allocation_size(KMeans * km) {
    if (!km) {
        return sizeof(KMeans);
    }
    return km->mbuf_size;
}

// passing any of these buffers in means you take responsibility for cleaning up all of them and that if they go out of scope or are free'd, no further actions are taken on the KMeans object
int KMeans_diagnostics(KMeans * km, void * means, void * last_means, KM_ID_TYPE * ids, unsigned char * workspace, size_t nbytes) {
    if (!km || km->iter || km->mbuf) {
        return 1; // cannot modify buffers after kmeans begins or once km->mbuf is allocated
    }
    int status = 0;
    if (means) {
        if (!km->means || km->flags & KMEANS_MEANS_BUFFERED) { // km->means not already set or previously set by user
            if (!km->means) {
                km->mbuf_size -= km->n * km->size;
            }
            km->means = (unsigned char *) means;
            km->flags |= KMEANS_MEANS_BUFFERED;
            printf("buffering km->means from user (%zu) bytes at (%p)\n", km->n * km->size, means);
        } else {
            status |= KMEANS_MEANS_BUFFERED; // buffering means failed
        }
    }
    if (last_means) { 
        if (!km->last_means || km->flags & KMEANS_LAST_MEANS_BUFFERED) { // km->last_means not already set or previously set by user
            if (!km->last_means) {
                km->mbuf_size -= km->n * km->size;
            }
            km->last_means = (unsigned char *) last_means;
            km->flags |= KMEANS_LAST_MEANS_BUFFERED;
            printf("buffering km->last_means from user (%zu) bytes at (%p)\n", km->n * km->size, last_means);
        } else {
            status |= KMEANS_LAST_MEANS_BUFFERED;
        }
    }
    if (ids) {
        if (!km->cluster_ids || km->flags & KMEANS_IDS_BUFFERED) { // km->cluster_ids not already set or previously set by user
            if (!km->cluster_ids) {
                km->mbuf_size -= km->ndata * sizeof(KM_ID_TYPE);
            }
            km->cluster_ids = ids;
            km->flags |= KMEANS_IDS_BUFFERED;
            printf("buffering km->cluster_ids from user (%zu) bytes at (%p)\n", km->ndata * sizeof(KM_ID_TYPE), ids);
        } else {
            status |= KMEANS_IDS_BUFFERED;
        }
    }
    if (workspace) {
        // need to use workspace for 'el', 'locs_buf', and any of the above that hasn't already been allocated
        km->message = workspace;
        km->message_size = nbytes;
        printf("buffering km->message from user workspace (%zu) bytes at %p\n", nbytes, (void *) workspace);
        // allocate for km->locs_buf if possible
        size_t nbytes_needed = km->n * sizeof(size_t);
        if (nbytes >= (nbytes_needed) && (!km->locs_buf || km->flags & KMEANS_LOCS_BUFFERED)) {
            if (!km->locs_buf) {

                km->mbuf_size -= nbytes_needed;
            }
            printf("buffering km->locs_buf from user workspace (%zu) bytes at %p\n", nbytes_needed, (void *) workspace);
            km->locs_buf = (size_t *) workspace;
            km->flags |= KMEANS_LOCS_BUFFERED;
            workspace += nbytes_needed;
            nbytes -= nbytes_needed;
        }

        // allocate for km->el if possible
        nbytes_needed = km->size;
        if (nbytes >= km->size && (!km->el || km->flags & KMEANS_EL_BUFFERED)) {
            if (!km->el) {
                km->mbuf_size -= nbytes_needed;
            }
            printf("buffering km->el from user workspace (%zu) bytes at %p\n", nbytes_needed, (void *) workspace);
            km->el  = workspace;
            km->flags |= KMEANS_EL_BUFFERED;
            workspace += nbytes_needed;
            nbytes -= nbytes_needed; // size fewer bytes are available for other types
        }

        // for the remainder, since we have dedicated channels to set these buffers, only use from the user-provided workspace is not already set by other means
        // so they will never overwrite other user-provided buffers and never not need to decrease mbuf_size if set

        // allocate for km->means if possible
        nbytes_needed = km->n * km->size;
        if (nbytes >= nbytes_needed && !km->means) {
            printf("buffering km->means from user workspace (%zu) bytes at %p\n", nbytes_needed, (void *) workspace);
            km->mbuf_size -= nbytes_needed;
            km->means = workspace;
            km->flags |= KMEANS_MEANS_BUFFERED;
            workspace += nbytes_needed;
            nbytes -= nbytes_needed;
        }

        if (nbytes >= nbytes_needed && !km->last_means) {
            printf("buffering km->last_means from user workspace (%zu) bytes at %p\n", nbytes_needed, (void *) workspace);
            km->mbuf_size -= nbytes_needed;
            km->last_means = workspace;
            km->flags |= KMEANS_LAST_MEANS_BUFFERED;
            workspace += nbytes_needed;
            nbytes -= nbytes_needed;
        }

        nbytes_needed = km->ndata * sizeof(KM_ID_TYPE);
        if (nbytes >= nbytes_needed && !km->cluster_ids) {
            printf("buffering km->cluster_ids from user workspace (%zu) bytes at %p\n", nbytes_needed, (void *) workspace);
            km->mbuf_size -= nbytes_needed;
            km->cluster_ids = (KM_ID_TYPE *) workspace;
            km->flags |= KMEANS_IDS_BUFFERED;
            workspace += nbytes_needed;
            nbytes -= nbytes_needed;
        }
        // TODO: might have to be careful about memory alignments here
    }
    return 0;
}

// resets a kmeans run so that it can be run again. this means clear internally allocated buffers
void KMeans_reset(KMeans * km) {
    KMeans_clean(km);
    km->iter = 0;
}

// delete a KMeans object created by KMeans_new
void KMeans_del(KMeans * km) {
    KMeans_clean(km);
    if (!(km->flags & KMEANS_BUFFERED)) {
        STATS_FREE(km);
    }
}

// does not create a new object. DO NOT CALL delete on this
KMeans * KMeans_new_from_buffer(void * buf, size_t buf_size, size_t * locs, size_t n, void * data, size_t ndata, size_t size, fkmdist_t fdist, fkmmean_t fmean) {
    if (!buf || buf_size < KMeans_allocation_size(NULL)) {
        return NULL;
    }
    KMeans * km = (KMeans *) buf;
    km->flags = 0;
    km->n = n;
    km->ndata = ndata;
    int status = KMeans_init(km, locs, n, data, ndata, size, fdist, fmean);
    if (status) {
        return NULL;
    }
    km->flags |= KMEANS_BUFFERED;

    return km;
}

KMeans * KMeans_new(size_t * locs, size_t n, void * data, size_t ndata, size_t size, fkmdist_t fdist, fkmmean_t fmean) {
    printf("allocating new KMeans struct\n");
    KMeans * km = (KMeans *) STATS_MALLOC(sizeof(KMeans));
    if (!km) {
        return NULL;
    }
    km->flags = 0; // needs to be set here due to buffer checking
    km->n = n; // set here so that I can track changes in KMeans_init
    km->ndata = ndata; // set here so that I can track changes in KMeans_init
    int status = KMeans_init(km, locs, n, data, ndata, size, fdist, fmean);
    if (status) {
        STATS_FREE(km);
        return NULL;
    }
    return km;
}

// 0s in any element other than diagnostic step mean "no change"
// these are completely optional changes to the behavior of kmeans
int KMeans_config(KMeans * km, unsigned int diagnostic_step, unsigned int max_iter, double tolerance, fkminit_t initialize, unsigned int flags) {
    if (km->iter) {
        return 1;
    }
    km->iter_step = diagnostic_step;
    if (max_iter) {
        km->max_iter = max_iter;
    }
    if (tolerance > 0.0) {
        km->tolerance = tolerance;
    }
    if (initialize) {
        km->initialize = initialize;
    }
    flags &= ~KMEANS_BUFFERED_MASK; // turn off attempts to change buffering
    if (flags) {
        km->flags |= flags;
    }
    if (km->iter_step >= km->max_iter) {
        return -1;
    }
    return 0;
}

size_t kmeans_lloyd_assign(KMeans * km, size_t idata) {
    size_t closest = 0;
    size_t size = km->size;
    void * el = (void *) (km->data + idata * size);
    double min_ = km->fdist((void *) km->means, el);

    for (size_t i = 1; i < km->n; i++) {
        double d = km->fdist((void *) (km->means + i * size), el);
        if (d < min_) {
            closest = i;
            min_ = d;
        }
    }
    return closest;
}

// generalized kmeans for not-necessarily-numeric data. the means must be members of the clusters they represent and are selected to be the member that minimizes the average distance to the mean
// the means are in the data and the data are sorted so that the clusters are contiguous, i.e. [mean_of_cluster0, member_cluster0,...,member_cluster0, mean_of_cluster1, member_cluster1,...]
// locations are the locations of the means in the final data array
// this is actually going to be a lot harder than I thought since the mean must be a member of the cluster
/*
int kmeansg(size_t * locs, size_t n, void ** data, size_t ndata, double (*dist)(void *, void *), unsigned int max_iter, double tolerance, unsigned int flags, \
    int (*initialize)(void ** means, size_t n, void ** data, size_t ndata, double (*dist)(void *, void *))) {
    // TODO: flag for Lloyd's algorithm
    if (!initialize) {
        initialize = kmeans_random_init;
    }
    bool changed = false;
    asdasd

    return 0;
}
*/

int kmeans_set_ids(KMeans * km) {
    //printf("setting cluster ids\n");
    if (!km || !km->cluster_ids) {
        return 1;
    }
    for (size_t i = 0; i < km->ndata; i++) {
        km->cluster_ids[i] = kmeans_lloyd_assign(km, i);
    }
    return 0;
}

// sets new means and returns the largest change in the means
double kmeans_set_means(KMeans * km) {//unsigned char * means, size_t * locs, unsigned int n, void * data, size_t ndata, KMeansCFG * cfg) {
    //printf("setting means\n");
    memmove(km->last_means, km->means, km->size * km->n);
    for (size_t i = 0; i < km->n; i++) {
        km->locs[i] = 0;
    }
    unsigned char * el = km->data;
    size_t i = 0; 
    while (i < km->ndata) {
        km->fmean((void *) (km->means + km->cluster_ids[i] * km->size), (void *) el, ++km->locs[km->cluster_ids[i]]);
        i++;
        el += km->size;
    }
    if (km->iter) {
        double max_change = 0.0;
        for (size_t i = 0; i < km->n; i++) {
            max_change = fmax(max_change, km->fdist((void *) (km->last_means + i * km->size), (void *)(km->means + i * km->size)));
        }
        return max_change;
    } else {
        return 2 * km->tolerance;
    }
}

void kmeans_swap(KMeans * km, size_t i, size_t j) {
    /*
    swapbuf((void *) (km->data + i * km->size), (void *) (km->data + j * km->size), km->size, (void *)km->el);
    KM_ID_TYPE buf = 0;
    swapbuf((void *) (km->cluster_ids + i), (void *) (km->cluster_ids + j), sizeof(KM_ID_TYPE), (void *) &buf);
    */
    swap((void *) (km->data + i * km->size), (void *) (km->data + j * km->size), km->size);
    swap((void *) (km->cluster_ids + i), (void *) (km->cluster_ids + j), sizeof(KM_ID_TYPE));
}

int kmeans_sort(KMeans * km) {
    printf("sorting output\n");
    //printf("before:\n");
    //print_array(km->locs, km->n, "%zu");
    // km->locs must be the count of the number of elements in each partition (after kmeans_set_means call, this is OK)
    // TODO: need to eliminate this by fixing the buffer situation
    //size_t insertion_points[100] = {0};
    size_t * locs_buf = km->locs_buf;
    size_t n = 0;
    for (size_t i = 0; i < km->n; i++) {
        locs_buf[i] = n;
        //insertion_points[i] = n;
        n += km->locs[i];
        km->locs[i] = n;
    }
    //printf("after locs_buf init:\n");
    //print_array(km->locs, km->n, "%zu");
    //print_array(locs_buf, km->n, "%zu");
    //print_array(km->cluster_ids, km->ndata, "%hu");
    // km->locs now points to the ends of the partitions
    assert(n == km->ndata);
    size_t i = 0;
    KM_ID_TYPE cur_cluster = 0;
    while (i < n) {
        if (i >= km->locs[cur_cluster]) {
            cur_cluster++;
        }
        //printf("in cluster %hu, checking element %zu with id %hu\n", cur_cluster, i, km->cluster_ids[i]);
        if (km->cluster_ids[i] == cur_cluster) { // in correct cluster, move on
            //printf("element at i = %zu is in correct cluster %hu\n", i, km->cluster_ids[i]);
            //insertion_points[cur_cluster]++;
            locs_buf[cur_cluster]++;
            i++;
        } else {
            //size_t next_loc = insertion_points[km->cluster_ids[i]];
            size_t next_loc = locs_buf[km->cluster_ids[i]];
            //printf("element at i = %zu belongs to cluster %hu (current = %zu). swapping with element at %zu\n", i, km->cluster_ids[i], cur_cluster, next_loc);
            //insertion_points[km->cluster_ids[i]]++;
            locs_buf[km->cluster_ids[i]]++;
            kmeans_swap(km, i, next_loc); // will keep repeating at iuntil it finds a value that is a valid element at this point in the cluster
        }
    }

    // now the km->locs should be point at the boundary locations appropriate for the output to kmeans()
    return 0;
}

// partitions the data into n clusters with clusters located at [0, locs[0]), [locs[0], locs[1]), etc.
// if means is null calls a generic kmeans "kmeansg" which takes a lot longer since means will be selected from the dataset
// if data is not of floating type, means MUST be NULL. The values that are nearest the true means are located at the beginning of each cluster partition
// returning -1 indicates suspended
// returning 0 means converged
// returning a positive value means failed
int kmeans(KMeans * km) {
    //printf("in kmeans\n");
    unsigned int next_stop = km->max_iter;
    if (km->iter_step) {
        next_stop = km->iter + km->iter_step;
        next_stop = (next_stop < km->max_iter) ? next_stop : km->max_iter;
    }
    int status = 0;
    double max_change = 2 * km->tolerance;
    if (!km->iter) {
        // allocate buffers needed to computer kmeans
        if (status = KMeans_alloc_buffers(km)) {
            goto kminit_fail;
        }
        //printf("initiating kmeans algorithm\n");
        //printf("data size = %zu\n", km->size);
        //printf("initializing\n");
        if (status = km->initialize(km->means, km->n, (void *) km->data, km->ndata, km->size, km->fdist)) {
            //printf("initialization failed");
            goto kminit_fail;
        }
        print_array(km->locs, km->n, "%zu");
        print_array(((double *)km->means), km->n, "%.8f");
    }

    // continue/resume the iteration

    if (!(km->flags & KMEANS_ALGO_MASK)) { // Lloyd's algorithm
        //printf("continuing with Lloyd's algorithm\n%zu < %zu && %.8f, %.8f\n", km->iter, next_stop, max_change, km->tolerance);
        while (km->iter < next_stop && max_change > km->tolerance) {
            //printf("starting iteration %zu\n", km->iter);
            // reassign elements cluster by cluster
            if ((status = kmeans_set_ids(km))) {
                //printf("kmeans_set_ids failed. code %u\n", status);
                break;
            }
            // recalculate means
            max_change = kmeans_set_means(km);
            km->iter++;
        }
        //printf("exited kmeans iteration loop\n");
        //print_array(km->locs, km->n, "%zu");
        //print_array(((double *)km->means), km->n, "%.8f");
        //print_array(km->cluster_ids, km->ndata, "%hu");
    } else {
        status = 3; // algorithm not understood
    }    

    if (!status) {
        if (km->iter >= km->max_iter || max_change < km->tolerance) {
            return kmeans_sort(km);
        } else {
            return -1;
        }
    }
kminit_fail:
    KMeans_clean(km);
    return status;
}