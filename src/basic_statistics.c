#include <math.h>
#include <stdbool.h>
#include "utils.h"
#include "basic_statistics.h"
#include "tables.h" // if Pascal's triangle is re-factored into a separate header, include that

/*

*/

typedef struct SelectMethod_params SelectMethod_params;

/* Helper methods for the quantile selection methods */

static inline double select_method_h_factor_(size_t n, double p, double alpha, double beta) {
    return ((n + alpha) * p + beta);
}

static size_t select_method_lower_floor_(double h_factor) {
    return (size_t) floor(h_factor);
}

static size_t select_method_lower_ceil_(double h_factor) {
    return (size_t) ceil(h_factor);
}

static size_t select_method_lower_ceil_half_(double h_factor) {
    return (size_t) ceil(h_factor - 0.5);
}

static size_t select_method_lower_round_(double h_factor) {
    return (size_t) round(h_factor);
}

// lower should be NAN if n_lower(h_factor) = 0;
static double select_method_qp_sided_(double lower, double upper, double h_factor) {
    return lower;
}

static double select_method_qp_mid_(double lower, double upper, double h_factor) {
    if (isnan(lower) || ceil(h_factor - 0.5) == floor(h_factor + 0.5)) {
        return lower;
    }
    if (isnan(upper)) {
        return upper;
    }
    return 0.5 * (lower + upper);
}

static double select_method_qp_interp_(double lower, double upper, double h_factor) {
    if (isnan(lower) || isnan(upper)) {
        return NAN;
    }
    return lower + (h_factor - floor(h_factor)) * (upper - lower);
}

struct SelectMethod_params {
    double alpha;
    double beta;
    // qp returns NULL if it fails
    double (*qp) (double lower, double upper, double h_factor);
    size_t (*n_lower) (double h_factor);
};

// 
static SelectMethod_params select_method_params[] = {
    {.alpha =  0.0, .beta =  0.0, .qp = select_method_qp_sided_ , .n_lower = select_method_lower_ceil_},
    {.alpha =  0.0, .beta =  0.5, .qp = select_method_qp_mid_,    .n_lower = select_method_lower_ceil_half_},
    {.alpha =  0.0, .beta = -0.5, .qp = select_method_qp_sided_,  .n_lower = select_method_lower_round_},
    {.alpha =  0.0, .beta =  0.0, .qp = select_method_qp_interp_, .n_lower = select_method_lower_floor_},
    {.alpha =  0.0, .beta =  0.5, .qp = select_method_qp_interp_, .n_lower = select_method_lower_floor_},
    {.alpha =  1.0, .beta =  0.0, .qp = select_method_qp_interp_, .n_lower = select_method_lower_floor_}, 
    {.alpha = -1.0, .beta =  1.0, .qp = select_method_qp_interp_, .n_lower = select_method_lower_floor_}, 
    {.alpha = (double)1.0/3, .beta = (double)1.0/3, .qp = select_method_qp_interp_, .n_lower = select_method_lower_floor_}, 
    {.alpha = (double)1.0/4, .beta = (double)3.0/8, .qp = select_method_qp_interp_, .n_lower = select_method_lower_floor_}
};

// calculates the Cramer sample moments and optionally weight moments for a_p and omega_p := sum_i(w_i ^p) for p in [min_moment, max_moment]
// sample_moments and weights_moments (if provided) must be pre-allocated to at least as long as max_moment-min_moment+1
// if weight_moments is NULL, the weight_moments will not be calculated
// if weights is NULL, will be assumed to be 1.
// max_moment >= 1
// might be worthwhile to have specific failure outputs
static int accumulate_sample_moments_weights_(double * sample_moments, double * samples, size_t n, unsigned short min_moment, unsigned short max_moment, double * weights, double * weight_moments) {
    // min_moment must be at least 1
    if (!min_moment) {
        min_moment = 1;
    }
    // min_moment must be >= min_moment
    if (max_moment < min_moment) {
        return STATS_FAILURE;
    }
    
    for (unsigned short j = 0; j < max_moment-min_moment+1; j++) {
        sample_moments[j] = 0.0;
    }
    if (weight_moments) {
        for (unsigned short j = 0; j < max_moment-min_moment+1; j++) {
            weight_moments[j] = 0.0;
        }
    }
    
    double poly_;
    double weight_sum = 0.0; // needed because weight_moments of order 1 (or equivalently sample_moments of order 0)
    for (size_t i = 0; i < n; i++) {        
        double weight = (weights) ? weights[i] : 1.0;
        weight_sum += weight;
        poly_ = weight * pow(samples[i], min_moment-1);
        for (unsigned short j = 0; j < max_moment-min_moment+1; j++) {
            poly_ *= samples[i];
            sample_moments[j] += poly_;
        }
        if (weight_moments) {
            poly_ = pow(weight, min_moment-1);
            for (unsigned short j = 0; j < max_moment-min_moment+1; j++) {
                poly_ *= weight; // polynomials all start at order 1
                weight_moments[j] += poly_;
            }
        }
    }

    for (unsigned short j = 0; j < max_moment-min_moment+1; j++) {
        sample_moments[j] /= weight_sum;
    }

    return STATS_SUCCESS;
}

double sample_moment(double * samples, size_t n, unsigned short order, double * weights) {
    double sample_moment = 0.0;

    if (accumulate_sample_moments_weights_(&sample_moment, samples, n, order, order, weights, NULL) != STATS_SUCCESS) {
        return NAN;
    }

    return sample_moment;
}

// TODO: fix to be general on the order
double central_sample_moment_biased(double * samples, size_t n, unsigned short order, double * weights) {
    double out = NAN;
    double * s_acc = (double *) STATS_MALLOC(sizeof(double)*order);
    if (!s_acc) {
        return out;
    }

    if (accumulate_sample_moments_weights_(s_acc, samples, n, 1, order, weights, NULL) == STATS_SUCCESS) {
        out = 0.0;
        for (unsigned short j = 0; j < order; j++) {
            out += (((order-j+1) & 1) ? -1 : 1) * (nCr(order, j+1) - kronecker_delta(j,0)) * s_acc[j] * pow(s_acc[0], order-j-1);
        }
    }

    // clean up
    STATS_FREE(s_acc);

    return out;
}

// public API for a standard interface
// unbiased versions will only be handled up to order 4 for now. if order is above 4, will calculate biased
/* 
It appears that the only moment that might matter for having an unbiased estimate is the
variance/2nd order. The variance in the estimates of order > 2 themselves are considerable for any 
relatively small n, which is the only time the bias is significant anyway
*/
double central_sample_moment(double * samples, size_t n, unsigned short order, double * weights, bool biased) {
    double result = 0.0;
    if (biased) {
        result = central_moment_biased_(samples, n, order, weights);
    } else {
        switch (order) {
            case (1): {
                result = central_moment_biased_(samples, n, order, weights); // central_moment_biased_ of order 1 is in fact unbiased
                break;
            }
            case (2): {
                result = unbiased_second_central_sample_moment_(samples, n, weights);
            }
            case (3): {
                result = unbiased_third_central_sample_moment_(samples, n, weights);
            }
            case (4): {
                result = unbiased_fourth_central_sample_moment_(samples, n, weights);
            }
            default: {
                result = central_moment_biased_(samples, n, order, weights); // higher order unbiased not implemented. return biased
                break;
            }
        }
    }
    
    return result;
}

static double unbiased_second_central_sample_moment_(double * samples, size_t n, double * weights) {
    unsigned short order = 2; // hard-coded for 2
    double out = NAN;
    if (n < order) {
        return out;
    }
    double * s_acc = (double *) STATS_MALLOC(sizeof(double)*order);
    if (!s_acc) {
        return out;
    }
    double * w_acc = (double *) STATS_MALLOC(sizeof(double)*order);
    if (!w_acc) {
        STATS_FREE(s_acc);
        return out;
    }

    out = 0.0;
    if (accumulate_sample_moments_weights_(s_acc, samples, n, 1, order, weights, w_acc) == STATS_SUCCESS) {
        out += (s_acc[1] - pow(s_acc[0], 2.0)); // moment factor
        out *= pow(w_acc[0], 2)/(pow(w_acc[0], 2) - w_acc[1]); // weight factor
    } else {
        out = NAN;
    }

    // clean up
    STATS_FREE(w_acc);
    STATS_FREE(s_acc);

    return out;
}

static double unbiased_third_central_sample_moment_(double * samples, size_t n, double * weights) {
    unsigned short order = 3; // hard-coded for 3
    double out = NAN;
    if (n < order) {
        return out;
    }
    double * s_acc = (double *) STATS_MALLOC(sizeof(double)*order);
    if (!s_acc) {
        return out;
    }
    double * w_acc = (double *) STATS_MALLOC(sizeof(double)*order);
    if (!w_acc) {
        STATS_FREE(s_acc);
        return out;
    }

    out = 0.0;
    if (accumulate_sample_moments_weights_(s_acc, samples, n, 1, order, weights, w_acc) == STATS_SUCCESS) {
        out += (s_acc[2] - 3 * s_acc[1] * s_acc[0] + 2 * pow(s_acc[0], 3)); // moment factor
        out *= pow(w_acc[0], 3)/(pow(w_acc[0], 3) - 3 * w_acc[0] * w_acc[1] + 2 * w_acc[2]);  // weight factor
    } else {
        out = NAN;
    }

    // clean up
    STATS_FREE(w_acc);
    STATS_FREE(s_acc);

    return out;
}

double unbiased_fourth_central_sample_moment_(double * samples, size_t n, double * weights) {
    unsigned short order = 4; // hard-coded for 3
    double out = NAN;
    if (n < order) {
        return out;
    }
    double * s_acc = (double *) STATS_MALLOC(sizeof(double)*order);
    if (!s_acc) {
        return out;
    }
    double * w_acc = (double *) STATS_MALLOC(sizeof(double)*order);
    if (!w_acc) {
        STATS_FREE(s_acc);
        return out;
    }
    out = 0.0;

    if (accumulate_sample_moments_weights_(s_acc, samples, n, 1, order, weights, w_acc) == STATS_SUCCESS) {
        // TODO: this can be simplified by pre-computation and I'm pretty sure there are canceling terms in the numerator (the next two lines)
        out += (pow(w_acc[0],4) - 3*w_acc[1]*pow(w_acc[0],2) + 2*w_acc[2]*w_acc[0] + 3*pow(w_acc[1],2) - 3*w_acc[3])*(s_acc[3] - 4*s_acc[2]*s_acc[0] + 6*s_acc[1]*pow(s_acc[0],2) - 3*pow(s_acc[0],4));
        out -= 3*(2*w_acc[1]*pow(w_acc[0],2) - 2*w_acc[2]*w_acc[0] - 3*pow(w_acc[1],2) + 3*w_acc[3])*pow((s_acc[1]-pow(s_acc[0],2)),2);
        out *= pow(w_acc[0],2)/(pow(w_acc[0],6) - 7*w_acc[1]*pow(w_acc[0],4) + 8*w_acc[2]*pow(w_acc[0],3) + 9*pow(w_acc[1],2)*pow(w_acc[0],2) - 6*w_acc[3]*pow(w_acc[0],2) - 8*w_acc[2]*w_acc[1]*w_acc[0] - 3*pow(w_acc[1],3) + 6*w_acc[3]*w_acc[1]);

        /*
        There's something wrong in my algebraic attempt to simplify the above expression
        out += (pow(w_acc[0],4) - 3*w_acc[1]*pow(w_acc[0],2) + 2*w_acc[2]*w_acc[0]+2*pow(w_acc[1],2) - 3*w_acc[3]) * (s_acc[3]-4*s_acc[2]*s_acc[0]);
        out -= 3*(2*w_acc[2]*pow(w_acc[0],2) - 2*w_acc[2]*w_acc[0] - 3*pow(w_acc[1],2) + 3*w_acc[3])*pow(s_acc[1],2);
        out += 3*pow(s_acc[0],2)*pow(w_acc[0],2)*(pow(w_acc[0],2) - w_acc[1]) * (2*s_acc[1] - pow(s_acc[0],2));
        out *= pow(w_acc[0],2) / (pow(w_acc[0],6) - 7*w_acc[1]*pow(w_acc[0],4) + 8*w_acc[2]*pow(w_acc[0],3) + 9*pow(w_acc[1],2)*pow(w_acc[0],2) - 6*w_acc[3]*pow(w_acc[0],2) - 8*w_acc[2]*w_acc[1]*w_acc[0] - 3*pow(w_acc[1],3) + 6*w_acc[3]*w_acc[1]);  // weight factor
        */
    } else {
        out = NAN;
    }

    // clean up
    STATS_FREE(w_acc);
    STATS_FREE(s_acc);

    return out;
}

double mean(double * samples, size_t n, double * weights) {
    return sample_moment(samples, n, 1, weights);
}

double stddev_biased(double * samples, size_t n, double * weights) {
    return sqrt(variance(samples, n, weights));
}

/*
estimator for the standard deviation
use dof = 0.0 for population standard deviation
dof = 1.0 for Bessel's correction estimation from the variance
dof = 1.5 for an approximation to the unbiased estimation of the standard deviation for normally distributed iid random variables
*/
double stddev(double * samples, size_t n, double * weights, float dof) {
    return stddev_biased(samples, n, weights)*sqrt((n-1.0)/(n - dof));
}

double variance(double * samples, size_t n, double * weights) {
    //return covariance(data, data, n, weights);
    return unbiased_second_central_sample_moment_(samples, n, weights);
}

// TODO: figure out how to work this into the new moment calculations
double correlation(double * samples, double * ref, size_t n, double * weights) {
    double yd = 0.0, ysqd = 0.0, xd = 0.0, xsqd = 0.0, xyd = 0.0, cor = 0.0;
    if (weights) {
        double w_acc = 0.0;
        for (size_t i = 0; i < n; i++) {
            w_acc += weights[i];
            xd += weights[i] * samples[i];
            yd += weights[i] * ref[i];
            ysqd += weights[i] * ref[i] * ref[i];
            xsqd += weights[i] * samples[i] * samples[i];
            xyd += weights[i] * samples[i] * ref[i];
        }
        cor = (w_acc * xyd - xd * yd) / sqrt((w_acc * xsqd - xd * xd) * (w_acc * ysqd - yd * yd));
    } else {
        for (size_t i = 0; i < n; i++) {
            xd += samples[i];
            yd += ref[i];
            ysqd += ref[i] * ref[i];
            xsqd += samples[i] * samples[i];
            xyd += samples[i] * ref[i];
        }
        cor = (n*xyd - xd*yd)/n/(n-1);
    }
    
    return cor;
}

// TODO: figure out how to work this into the new moment calculations
double covariance(double * samples, double * ref, size_t n, double * weights) {
    double yd = 0.0, xd = 0.0, xyd = 0.0, cov = 0.0;
    if (weights) {
        double w_acc = 0.0, wsq_acc = 0.0;
        for (size_t i = 0; i < n; i++) {
            w_acc += weights[i];
            wsq_acc += weights[i] * weights[i];
            xd += weights[i] * samples[i];
            yd += weights[i] * ref[i];
            xyd += weights[i] * samples[i] * ref[i];
        }
        cov = (w_acc*xyd - xd*yd) / (w_acc * w_acc - wsq_acc);
    } else {
        for (size_t i = 0; i < n; i++) {
            xd += samples[i];
            yd += ref[i];
            xyd += samples[i] * ref[i];
        }
        cov = (n*xyd - xd*yd)/n/(n-1);
    }
    
    return cov;
}

// returns NAN if the quantile cannot be calculated
// WARNING: due to dependence on quickselect, this function will modify the data array in place by partially sorting
double quantile(double * samples, size_t n, double quantile, enum select_method method) {
    if (!method) {
        method = DEFAULT_SELECT_METHOD;
    }
    SelectMethod_params smp = select_method_params[method - 1];
    double h_factor = select_method_h_factor_(n, quantile, smp.alpha, smp.beta);
    size_t _n_lower = smp.n_lower(h_factor);
    double lower, upper;
    if (_n_lower) {
        lower = *(double *)quickselect(samples, n, _n_lower-1, sizeof(double), compare_double);
    } else {
        lower = NAN;
    }
    if (_n_lower < n) {
        upper = *(double *)quickselect(samples, n, _n_lower, sizeof(double), compare_double);
    } else {
        upper = NAN;
    }
    return smp.qp(lower, upper, h_factor);
}

double median(double * samples, size_t n, enum select_method method) {
    return quantile(samples, n, 0.5, method);
}

void minmax(double * samples, size_t n, double * min, double * max) {
    if (!n) {
        return;
    }
    double _min = samples[0];
    double _max = samples[0];
    for (size_t i = 1; i < n; i++) {
        if (samples[i] > _max) {
            _max = samples[i];
        }
        if (samples[i] < _min) {
            _min = samples[i];
        }
    }

    if (min) {
        *min = _min;
    }

    if (max) {
        *max = _max;
    }

}

// not the fastest, but convenient
double min(double * samples, size_t n) {
    double min;
    minmax(samples, n, NULL, &min);
    return min;
}

// not the fastest, but convenient
double max(double * samples, size_t n) {
    double max;
    minmax(samples, n, NULL, &max);
    return max;
}

double range(double * samples, size_t n) {
    double max, min;
    minmax(samples, n, &min, &max);
    return (max - min);
}

// TODO: make this O(n) by implementing a hash table on the array. For now, sort and then bin, which is O(n*log(n))
// currently modifies the data array in place by sorting it
double mode(double * samples, size_t n) {
    // sort O(n*log(n))
    qsort(samples, n, sizeof(double), compare_double); // could also use quicksort in utils.h

    // bin the values by counting number of repititions
    size_t max_rep = 0, cur_rep = 0;
    double cur = samples[0], _mode = samples[0];
    for (size_t i = 0; i < n; i++) {
        if (samples[i] == cur) {
            cur_rep++;
        } else {
            if (cur_rep > max_rep) {
                _mode = cur;
                max_rep = cur_rep;
            }
            cur_rep = 1;
            cur = samples[i];
        }
    }
    return _mode;
}