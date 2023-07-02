#include <stdlib.h>
#include "utils.h"
#include "basic_statistics.h"
#include "histogram.h"

struct Histogram {
    double * bin_bounds;    // there are n_bins+1 bounds. Exact values are stored in the bin to the right
    double * bin_centers;   // there are n_bins center. These are the values that the bins "represent"
    size_t n_bins;          // number of frequencies/bins to track
    size_t * freqs;         // frequencies of bins
    size_t n_data;          // number of data points in the freqs
}; // really a 1-D histogram

// a copy of bin_bounds is made
Histogram * Histogram_new_from_bounds(double * data, size_t n, double * bin_bounds, size_t n_bin_bounds) {
    Histogram * hist = (Histogram *) STATS_MALLOC(sizeof(Histogram));
    Histogram_init(hist);

    // copy data
    hist->bin_bounds = (double *) STATS_MALLOC(sizeof(double)*n_bin_bounds);
    memcpy(hist->bin_bounds, bin_bounds, sizeof(double)*n_bin_bounds);
    hist->n_bins = n_bin_bounds - 1;
    hist->bin_centers = (double *) STATS_MALLOC(sizeof(double)*hist->n_bins);
    hist->freqs = (size_t *) STATS_MALLOC(sizeof(size_t)*hist->n_bins);
    for (int i = 0; i < hist->n_bins; i++) {
        hist->bin_centers[i] = 0.5 * (bin_bounds[i] + bin_bounds[i+1]);
        hist->freqs[i] = 0;
    }

    Histogram_bin_(hist, data, n);

    return hist;
}

// will attempt to center the bins on the range
Histogram * Histogram_new_from_bin_size(double * data, size_t n, double bin_size) {
    double min;
    double max;
    minmax(data, n, &min, &max);

    Histogram * hist = (Histogram *) STATS_MALLOC(sizeof(Histogram));
    Histogram_init(hist);

    hist->n_bins = (max - min) / bin_size + 1; // + 1 to ensure min/max are within a bin
    hist->bin_bounds = (double *) STATS_MALLOC(sizeof(double)*(hist->n_bins+1));
    // so that bin_bounds[hist->n_bins] =  bin_bounds[0] + hist->n_bins * bins_size = min - bin_size / 2 + (max - min) + bins_size = max + bin_size / 2 > max
    hist->bin_bounds[0] = min - bin_size / 2; 
    hist->bin_centers = (double *) STATS_MALLOC(sizeof(double)*hist->n_bins);
    hist->freqs = (size_t *) STATS_MALLOC(sizeof(size_t)*hist->n_bins);

    for (int i = 0; i < hist->n_bins;  i++) {
        hist->bin_bounds[i+1] = hist->bin_bounds[i] + bin_size;
        hist->bin_centers[i] = hist->bin_bounds[i] + bin_size;
        hist->freqs[i] = 0;
    }

    Histogram_bin_(hist, data, n);

    return hist;
}

Histogram * Histogram_new_from_n_bins(double * data, size_t n, size_t n_bins) {   
    double min;
    double max;
    minmax(data, n, &min, &max);
    // cannot make a histogram with 0 bin widths
    if (max == min) {
        return NULL; // TODO: add failure ptr to utils.h
    }

    Histogram * hist = (Histogram *) STATS_MALLOC(sizeof(Histogram));
    Histogram_init(hist);

    hist->n_bins = n_bins;
    double bin_size = (max - min) / (hist->n_bins - 1);
    hist->bin_bounds = (double *) STATS_MALLOC(sizeof(double)*(n_bins+1));
    // so that bin_bounds[hist->n_bins] =  bin_bounds[0] + hist->n_bins * bins_size = min - bin_size / 2 + (max - min) + bins_size = max + bin_size / 2 > max
    hist->bin_bounds[0] = min - bin_size / 2;
    hist->bin_centers = (double *) STATS_MALLOC(sizeof(double)*hist->n_bins);
    hist->freqs = (size_t *) STATS_MALLOC(sizeof(size_t)*hist->n_bins);
    for (int i = 0; i < hist->n_bins;  i++) {
        hist->bin_bounds[i+1] = hist->bin_bounds[i] + bin_size;
        hist->bin_centers[i] = hist->bin_bounds[i] + bin_size/2;
        hist->freqs[i] = 0;
    }

    Histogram_bin_(hist, data, n);

    return hist;
}

// a copy of bin_centers is made
Histogram * Histogram_new_from_centers(double * data, size_t n, double * bin_centers, size_t n_bin_centers) {
    double min;
    double max;
    minmax(data, n, &min, &max);

    Histogram * hist = (Histogram *) STATS_MALLOC(sizeof(Histogram));
    Histogram_init(hist);
    
    hist->n_bins = n_bin_centers;

    hist->bin_bounds = (double *) STATS_MALLOC(sizeof(double)*(hist->n_bins+1));
    hist->bin_bounds[0] = min;
    hist->bin_centers = (double *) STATS_MALLOC(sizeof(double)*hist->n_bins);
    memcpy(hist->bin_centers, bin_centers, sizeof(double) * n_bin_centers);
    size_t * freqs = (size_t *) STATS_MALLOC(sizeof(size_t)*hist->n_bins);
    for (size_t i = 0; i < hist->n_bins - 1;  i++) {
        hist->bin_bounds[i+1] = 0.5 * (bin_centers[i] + bin_centers[i+1]);
        hist->freqs[i] = 0;
    }

    hist->freqs[hist->n_bins - 1] = 0;
    hist->bin_bounds[hist->n_bins] = 1.01 * (max - bin_centers[hist->n_bins - 1]);

    Histogram_bin_(hist, data, n);

    return hist;
}

static void Histogram_bin_(Histogram * hist, double * data, size_t n) {
    for (size_t i = 0; i < n; i++) {
        size_t ind = bisect_left(hist->bin_bounds, data + i, 0, hist->n_bins - 1, sizeof(double), compare_double);
        // TODO: do i have to check if ind > hist->n_bins - 1?
        hist->freqs[ind] += 1;
        hist->n_data++;
    }
}

// TODO: consider making this the master initialization which takes all inputs
void Histogram_init(Histogram * hist) {
    hist->bin_bounds = NULL;
    hist->bin_centers = NULL;
    hist->n_bins = 0;
    hist->freqs = NULL;
    hist->n_data = 0;
}

void Histogram_del(Histogram * hist) {
    free(hist->bin_bounds);
    free(hist->bin_centers);
    free(hist->freqs);
}

/*
double hist_mean(Histogram * hist);

double hist_variance(Histogram * hist);

double hist_stddev(Histogram * hist);

double hist_quantile(Histogram * hist, double quantile); // might have to add a method

// calls hist_quantile(hist, 0.5);
double hist_median(Histogram * hist);

double hist_mode(Histogram * hist);
*/
