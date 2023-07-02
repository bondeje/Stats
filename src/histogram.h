#ifndef HISTOGRAM_H
#define HISTOGRAM_H

typedef struct Histogram Histogram;

// a copy of bin_bounds is made
Histogram * Histogram_new_from_bounds(double * data, double * bin_bounds);

Histogram * Histogram_new_from_bin_size(double * data, size_t n, double bin_size);

Histogram * Histogram_new_from_n_bins(double * data, size_t n_bins);

// a copy of bin_centers is made
Histogram * Histogram_new_from_centers(double * data, double * bin_centers);

void Histogram_init(Histogram * hist);

void Histogram_del(Histogram * hist);

// returns a copy of the frequencies normalized so that the sum of all frequencies is 1
double * Histogram_get_density(Histogram * hist);

double hist_mean(Histogram * hist);

double hist_variance(Histogram * hist);

double hist_stddev(Histogram * hist);

double hist_quantile(Histogram * hist, double quantile); // might have to add a method

// calls hist_quantile(hist, 0.5);
double hist_median(Histogram * hist);

double hist_mode(Histogram * hist);

#endif // HISTOGRAM_H