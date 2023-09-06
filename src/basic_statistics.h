#ifndef BASIC_STATISTICS_H
#define BASIC_STATISTICS_H

/*
any function or variable that has "select_method" in it is used for quantile determination
For a summary of what is going on, see the wikipedia article below or the associated paper from Hyndman & Fan
https://en.wikipedia.org/wiki/Quantile#Estimating_quantiles_from_a_sample
*/
enum select_method {
    SM_INVERTED_CDF = 1,
    SM_AVERAGED_INVERTED_CDF,
    SM_CLOSEST_OBSERVATION,
    SM_INTERPOLATED_INVERTED_CDF,
    SM_HAZEN,
    SM_WEIBULL,
    SM_LINEAR,             // This is what Python uses in most common cases
    SM_MEDIAN_UNBIASED,    // 8, should be default
    SM_NORMAL_UNBIASED
};

#define DEFAULT_SELECT_METHOD SM_MEDIAN_UNBIASED

/*
estimator for the sample mean
*/
double mean(double * data, size_t n, double * weights);

/*
estimator for the standard deviation
ddof not yet implemented, so this is a biased, but Bessel-corrected standard deviation
use dof = 0.0 for population standard deviation
dof = 1.0 for Bessel's correction estimation from the variance
dof = 1.5 for an approximation to the unbiased estimation of the standard deviation for normally distributed iid random variables
*/
double stddev_biased(double * data, size_t n, double * weights);

// unbiased estimator for the sample variance
double variance(double * data, size_t n, double * weights);

// unbiased estimator for the sample covariance
double covariance(double * data, double * ref, size_t n, double * weights);

// estimator for the sample correlation coefficient
double correlation(double * data, double * ref, size_t n, double * weights);

double quantile(double * data, size_t n, double quantile, enum select_method method);

double median(double * data, size_t n, enum select_method method);

void minmax(double * data, size_t n, double * min, double * max);

#undef min // because windef.h has the min/max macros defined
double min(double * data, size_t n);

#undef max // because windef.h has the min/max macros defined
double max(double * data, size_t n);

double range(double * data, size_t n);

double mode(double * data, size_t n);

#endif // BASIC_STATISTICS_H