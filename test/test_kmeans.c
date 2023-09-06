#ifdef DEBUG_MALLOC
    #include "../../mallocs/debug_malloc.h"
#endif
#include <stddef.h>
#include <stdio.h>
#include <time.h>
#include "../src/utils.h"
#include "../src/kmeans.h"

#define DATA_SIZE 21
#define BIG_DATA_SIZE (size_t)10000000

void test_big_data_double(void) {
    printf("in test_big_data_double\n");
    double * data = (double *) STATS_MALLOC(BIG_DATA_SIZE * sizeof(double));
    if (!data) {
        printf("failed to allocate enough data\n");
        return;
    }
    KM_ID_TYPE * ids = (KM_ID_TYPE *) STATS_MALLOC(BIG_DATA_SIZE * sizeof(KM_ID_TYPE));
    if (!ids) {
        printf("failed to allocate enough ids\n");
        STATS_FREE(ids);
        return;
    }
    for (size_t i = 0; i < BIG_DATA_SIZE; i++) {
        data[i] = 2 *(random() - 0.5);
        ids[i] = 0;
    }
    KM_ID_TYPE nclusters = 10;
    size_t ndata = BIG_DATA_SIZE;
    size_t partitions[10];
    double means[10] = {0.0};
    double last_means[10] = {0.0};
    unsigned char kmeans_buf[160] = {'\0'};
    unsigned char workspace[10 * sizeof(size_t) + sizeof(data[0])] = {'\0'};

    KMeans * km = KMeans_new_from_buffer(kmeans_buf, 160, partitions, nclusters, data, ndata, sizeof(double), NULL, NULL);
    KMeans_diagnostics(km, means, last_means, ids, workspace, nclusters * sizeof(size_t) + sizeof(data[0]));

    clock_t t = clock();
    int status = kmeans(km);
    t = clock() - t;
    printf ("It took me %f seconds to cluster %zu double elements into %hu clusters\n",((float)t)/CLOCKS_PER_SEC, BIG_DATA_SIZE, nclusters);

    printf("status of kmeans: %u\n", status);

    print_array(partitions, nclusters, "%zu");
    print_array(means, nclusters, "%.8f");
    //print_array(data, ndata, "%.8f");
    //print_array(ids, ndata, "%hu");
    STATS_FREE(data);
    STATS_FREE(ids);
    KMeans_del(km);
}

void basic_double_no_alloc(void) {
    // inputs
    double data[DATA_SIZE];
    size_t nclusters = 0;
    size_t ndata = DATA_SIZE;

    // outpus
    size_t partitions[DATA_SIZE] = {0};

    // buffers/outputs
    unsigned char kmeans_buf[160] = {'\0'};
    double means[DATA_SIZE] = {0.0};
    double last_means[DATA_SIZE] = {0.0};
    KM_ID_TYPE ids[DATA_SIZE] = {0};
    unsigned char workspace[DATA_SIZE * sizeof(KM_ID_TYPE) + sizeof(data[0])] = {'\0'};

    for (size_t i = 0; i < DATA_SIZE; i++) {
        data[i] = 2 *(random() - 0.5);
    }

    // two clusters
    nclusters = 2;
    KMeans * km = KMeans_new_from_buffer(kmeans_buf, 160, partitions, nclusters, data, ndata, sizeof(double), NULL, NULL);
    KMeans_diagnostics(km, means, last_means, ids, workspace, DATA_SIZE * sizeof(KM_ID_TYPE) + sizeof(data[0]));
    printf("status of kmeans: %u\n", kmeans(km));

    print_array(partitions, nclusters, "%zu");
    print_array(means, nclusters, "%.8f");
    print_array(data, ndata, "%.8f");
    print_array(ids, ndata, "%hu");

    nclusters = 5;
    for (size_t i = 0; i < DATA_SIZE; i++) {
        data[i] = 2 *(random() - 0.5);
    }
    KMeans_init(km, partitions, nclusters, data, ndata, sizeof(double), NULL, NULL);
    KMeans_diagnostics(km, means, last_means, ids, workspace, DATA_SIZE * sizeof(KM_ID_TYPE) + sizeof(data[0]));

    printf("status of kmeans: %u\n", kmeans(km));

    print_array(partitions, nclusters, "%zu");
    print_array(means, nclusters, "%.8f");
    print_array(data, ndata, "%.8f");
    print_array(ids, ndata, "%hu");
    
    KMeans_del(km);
    //debug_print_snapshot();
}

void basic_double(void) {
    // inputs
    double data[DATA_SIZE];
    size_t nclusters = 0;
    size_t ndata = DATA_SIZE;

    // outpus
    size_t partitions[DATA_SIZE] = {0};

    // buffers/outputs
    double means[DATA_SIZE] = {0.0};
    double last_means[DATA_SIZE] = {0.0};
    KM_ID_TYPE ids[DATA_SIZE] = {0};
    unsigned char workspace[DATA_SIZE * sizeof(KM_ID_TYPE) + sizeof(data[0])] = {'\0'};

    for (size_t i = 0; i < DATA_SIZE; i++) {
        data[i] = 2 *(random() - 0.5);
    }

    // two clusters
    nclusters = 2;
    KMeans * km = KMeans_new(partitions, nclusters, data, ndata, sizeof(double), NULL, NULL);
    KMeans_diagnostics(km, means, last_means, ids, workspace, DATA_SIZE * sizeof(KM_ID_TYPE) + sizeof(data[0]));
    printf("status of kmeans: %u\n", kmeans(km));

    print_array(partitions, nclusters, "%zu");
    print_array(means, nclusters, "%.8f");
    print_array(data, ndata, "%.8f");
    print_array(ids, ndata, "%hu");

    nclusters = 5;
    for (size_t i = 0; i < DATA_SIZE; i++) {
        data[i] = 2 *(random() - 0.5);
    }
    KMeans_init(km, partitions, nclusters, data, ndata, sizeof(double), NULL, NULL);
    KMeans_diagnostics(km, means, last_means, ids, workspace, DATA_SIZE * sizeof(KM_ID_TYPE) + sizeof(data[0]));

    printf("status of kmeans: %u\n", kmeans(km));

    print_array(partitions, nclusters, "%zu");
    print_array(means, nclusters, "%.8f");
    print_array(data, ndata, "%.8f");
    print_array(ids, ndata, "%hu");
    
    KMeans_del(km);
}

int main(int narg, char ** args) {
    SPRNG(0);
    //basic_double_no_alloc();
    basic_double();
    //test_big_data_double();
    //printf("size of KMeans object: %zu\n", KMeans_allocation_size(NULL));
    return 0;
}