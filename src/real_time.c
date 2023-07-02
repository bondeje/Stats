#include <math.h>
#include "utils.h"
#include "basic_statistics.h"
#include "tables.h" // nCr
#include "real_time.h"

/*
TODO: build monitor for tracking streaming updates
*/

/*
updates the average as a rolling average using O(1) algorithm.
assumes window_size points have already been accounted for in the avg
does not handle weighted averages
avg must not be NULL
WARNING: initialization must be handled separately/results may not make sense until at least window_size data have been added
*/ 
int rt_mean(double * avg, double new_data, size_t window_size, double old_data) {
    if (avg) {
		*avg += (new_data - old_data)/window_size;
	} else { // any of the required input pointers are null
		return STATS_FAILURE;
	}
    return STATS_SUCCESS;
}

/*
updates the variance as a rolling standard deviation using O(1) algorithm. 
also updates the mean
assumes window_size points have already been accounted for in the avg
does not handle weighted standard deviations
var must not be NULL
WARNING: initialization must be handled separately/results may not make sense until at least window_size data have been added
*/
int rt_variance(double * var, double * old_mean, double new_data, size_t window_size, double old_data) {
    return rt_covariance(var, old_mean, old_mean, new_data, new_data, window_size, old_data, old_data);
}


// TODO: update for self-starting.
int rt_covariance(double * cov_ab, double * old_mean_a, double * old_mean_b, double new_data_a, double new_data_b, size_t window_size, double old_data_a, double old_data_b) {
    if (cov_ab) {
        double new_mean_a = *old_mean_a;
        double new_mean_b = *old_mean_b;
        if ((rt_average(&new_mean_a, new_data_a, window_size, old_data_a) == STATS_SUCCESS) && (rt_average(&new_mean_b, new_data_b, window_size, old_data_b) == STATS_SUCCESS)) {
            *cov_ab += ( (new_data_a - new_mean_a) * (new_data_b - new_mean_b) - (old_data_a - new_mean_a) * (old_data_b - new_mean_b) + window_size * (new_mean_b - *old_mean_b) * (new_mean_a - *old_mean_a)) / (window_size - 1);
#ifndef RT_NO_UPDATE
            if (old_mean_a != old_mean_b) { // called from variance, only have to update one of the means
                *old_mean_a = new_mean_a;
            }
            *old_mean_b = new_mean_b;
#endif // RT_NO_UPDATE
        }
    } else {
        return STATS_FAILURE;
    }
    return STATS_SUCCESS;
}

// if old_datum is NULL, window_size is treated as incrementing by 1 to include the new data point, else window_size is const
// moments must be at least max_order in length
int rt_moments(double * moments, unsigned short max_order, size_t window_size, double new_datum, double old_datum) {
    double x = 1.0;
    if (isnan(old_datum)) {
        for (unsigned short order = 0; order < max_order; order++) {
            x *= new_datum;
            moments[order] = (x + window_size * moments[order]) / (window_size+1);
        }
    } else {
        double y = 1.0;
        for (unsigned short order = 0; order < max_order; order++) {
            x *= new_datum; // x = new_datum ^ (order + 1)
            y *= old_datum; // y = old_datum ^ (order + 1)
            moments[order] += (x - y)/window_size;
        }
    }

    return STATS_SUCCESS;
}

// if old_datum is NULL, window_size is treated as incrementing by 1 to include the new data point, else window_size is const
// moments must be at least max_order-CM_ORDER_OFFSET in length
// this could be made more efficient by tracking the delta means
int rt_central_moments(double * central_moments, unsigned short max_order, size_t window_size, double new_mean, double old_mean, double new_datum, double old_datum) {
    double new_cm;
    double new_dev = new_datum - new_mean;
    double new_dev_pow = pow(new_dev, max_order + 1);
    unsigned short p;
    double mean_dev;
    double mean_dev_pow;
    if (isnan(old_datum)) { // add to window
        for (int order = max_order-CM_ORDER_OFFSET-1; order > -1; order--) {
            mean_dev = new_mean - old_mean;
            mean_dev_pow = 1.0;
            p = order + CM_ORDER_OFFSET+1;
            new_dev_pow /= new_dev; // = (new_datum - new_mean) ^ p
            new_cm = 0.0;


            for (int k = 1; k < p - 1; k++) {
                mean_dev_pow *= mean_dev; // = (new_mean - old_mean) ^ k
                new_cm += nCr(p, k) * ((k & 1) ? -1.0 : 1.0) * mean_dev_pow * central_moments[order - k];
            }

            mean_dev_pow *= mean_dev * mean_dev; // = (new_min - old_mean) ^ p
            
            new_cm = new_dev_pow + window_size * (new_cm + central_moments[order] + ((p & 1) ? -1.0 : 1.0) * mean_dev_pow);

            central_moments[order] = new_cm / (window_size + 1.0);
        }
    } else {
        double old_dev = old_datum - old_mean;
        double old_dev_pow = pow(old_dev, max_order + 1);
        for (int order = max_order-CM_ORDER_OFFSET-1; order > -1; order--) {
            mean_dev = new_mean - old_mean;
            mean_dev_pow = 1.0;
            p = order + CM_ORDER_OFFSET+1;
            new_dev_pow /= new_dev; // = (new_datum - new_mean) ^ p
            old_dev_pow /= old_dev; // = (old_datum - new_mean) ^ p
            new_cm = 0.0;


            for (int k = 1; k < p - 1; k++) {
                mean_dev_pow *= mean_dev; // = (new_mean - old_mean) ^ k
                new_cm += nCr(p, k) * ((k & 1) ? -1.0 : 1.0) * mean_dev_pow * central_moments[order - k];
            }

            mean_dev_pow *= mean_dev * mean_dev; // = (new_min - old_mean) ^ p
            
            new_cm = (new_dev_pow - old_dev_pow)/window_size + (new_cm + central_moments[order] + ((p & 1) ? -1.0 : 1.0) * mean_dev_pow);

            central_moments[order] = new_cm;
        }
    }

    return STATS_SUCCESS;
}

static double rt_unbiased_central_moment_4_(double * central_moments, size_t n) {
    if (n < 4) {
        return 0.0;
    }
    double m4 = central_moments[4 - CM_ORDER_OFFSET - 1];
    double m2 = central_moments[2 - CM_ORDER_OFFSET - 1];
    return (n*(n*n - 2*n + 3) * m4 - 2*n*(2*n-3) * m2 * m2)/((n-1) * (n-2) * (n-3));
}

static double rt_unbiased_central_moment_3_(double * central_moments, size_t n) {
    if (n < 3) {
        return 0.0;
    }
    return n*n/((n-1) * (n-2))*central_moments[3 - CM_ORDER_OFFSET - 1];
}

static double rt_unbiased_central_moment_2_(double * central_moments, size_t n) {
    if (n < 2) {
        return 0.0;
    }
    return n/(n-1) * central_moments[2 - CM_ORDER_OFFSET - 1];
}

static double (*rt_unbiased_central_moments_[]) (double *, size_t) = {rt_unbiased_central_moment_2_, rt_unbiased_central_moment_3_, rt_unbiased_central_moment_4_};


// updating unbiased_central_moments only requires the already updated central_moments, which must be supplied in the inputs of length at least max_order
int rt_unbiased_central_moments(double * unbiased_central_moments, double * central_moments, unsigned short max_order, size_t window_size) {
    max_order = (max_order > MAX_UNBIASED_CENTRAL_MOMENT) ? MAX_UNBIASED_CENTRAL_MOMENT : max_order;
    unsigned short p;
    for (unsigned short order = max_order - CM_ORDER_OFFSET - 1) {
        unbiased_central_moments[order] = rt_unbiased_central_moments_(central_moments, window_size);
    }
    return STATS_SUCCESS;
}

RTStatsMonitor * RTStatsMonitor_new(size_t n_channels, size_t window_size, size_t order_moments, size_t order_central_moments, size_t order_unbiased_central_moments, size_t n_quantiles, double * ps) {
    if (!n_channels) {
        goto failed_rtsm_alloc;
    }
    RTStatsMonitor * rtsm = (RTStatsMonitor *) STATS_MALLOC(sizeof(RTStatsMonitor));

    if (!rtsm) {
        goto failed_rtsm_alloc;
    }

    order_moments = (!order_moments) ? 1 : order_moments; // must collect at least the mean
    
    rtsm->data = NULL;
    if (window_size) {
        if (!(rtsm->data = (double *) STATS_MALLOC(window_size * n_channels * sizeof(double)))) {
            goto failed_data_alloc;
        }       
    } else {
        rtsm->data = NULL;
    }

    rtsm->moments = NULL;
    if (order_moments && !(rtsm->moments = (double *) STATS_MALLOC(order_moments * n_channels * sizeof(double)))) {
        goto failed_moments_alloc;
    }

    rtsm->central_moments = NULL;
    order_central_moments = (order_central_moments > order_unbiased_central_moments) ? order_central_moments : order_unbiased_central_moments;
    if ((order_central_moments >= CM_ORDER_OFFSET) && !(rtsm->central_moments = (double *) STATS_MALLOC((order_central_moments - CM_ORDER_OFFSET) * n_channels * sizeof(double)))) {
        goto failed_central_moments_alloc;
    }

    rtsm->unbiased_central_moments = NULL;
    if (order_unbiased_central_moments > MAX_UNBIASED_CENTRAL_MOMENT) {
        order_unbiased_central_moments = MAX_UNBIASED_CENTRAL_MOMENT;
    }
    if (order_unbiased_central_moments && !(rtsm->unbiased_central_moments = (double *) STATS_MALLOC(order_unbiased_central_moments * n_channels * sizeof(double)))) {
        goto failed_unbiased_central_moments_alloc;
    }

    rtsm->quantiles = NULL;
    rtsm->ps = NULL;
    if (n_quantiles) {
        if (!(rtsm->quantiles = (double *) STATS_MALLOC(n_quantiles * n_channels * sizeof(double)))) {
            goto failed_quantiles_alloc;
        }
        if (!(rtsm->ps = (double *) STATS_MALLOC(n_quantiles * sizeof(double)))) {
            goto failed_ps_alloc;
        }
        memcpy(rtsm->ps, ps, sizeof(double)*n_quantiles);
    }}

    rtsm->_n_buf = (order_moments > order_unbiased_central_moments - CM_ORDER_OFFSET) ? order_moments : order_unbiased_central_moments - CM_ORDER_OFFSET;
    if (!(rtsm->_buf = (double *) STATS_MALLOC(rtsm->_n_buf * sizeof(double)))) {
        goto failed_buffer_alloc;
    }

    RTStatsMonitor_init(rtsm, n_channels, window_size, order_moments, order_central_moments, order_unbiased_central_moments, n_quantiles, ps);

    return rtsm;

failed_buffer_alloc:
    STATS_FREE(rtsm->ps);

failed_ps_alloc:
    STATS_FREE(rtsm->quantiles);

failed_quantiles_alloc:
    STATS_FREE(rtsm->unbiased_central_moments);

failed_unbiased_central_moments_alloc:
    STATS_FREE(rtsm->central_moments);

failed_central_moments_alloc:
    STATS_FREE(rtsm->moments);

failed_moments_alloc:
    STATS_FREE(rtsm->data);

failed_data_alloc:
    STATS_FREE(rtsm);

failed_rtsm_alloc:

    return NULL;
}

void RTStatsMonitor_init(RTStatsMonitor * rtsm, size_t n_channels, size_t window_size, size_t order_moments, size_t order_central_moments, size_t order_unbiased_central_moments, size_t n_quantiles, double * ps) {
    if (!n_channels) {
        return;
    }
    rtsm->n_channels = n_channels;
    rtsm->_n_samples = 0;
    if (window_size) {
        for (size_t i = 0; i < n_channels*window_size; i++) {
            rtsm->data[i] = 0.0;
        }
    }
    rtsm->window_size = window_size;
    if (order_moments) {
        for (size_t i = 0; i < n_channels*order_moments; i++) {
            rtsm->moments[i] = 0.0;
        }
    }
    rtsm->order_moments = order_moments;
    if (order_central_moments) {
        for (size_t i = 0; i < n_channels*order_central_moments; i++) {
            rtsm->central_moments[i] = 0.0;
        }
    }
    rtsm->order_central_moments = order_central_moments;
    if (order_unbiased_central_moments) {
        for (size_t i = 0; i < n_channels*order_unbiased_central_moments; i++) {
            rtsm->unbiased_central_moments[i] = 0.0;
        }
    }
    rtsm->order_unbiased_central_moments = order_unbiased_central_moments;
    if (n_quantiles) {
        for (size_t i = 0; i < n_channels*n_quantiles; i++) {
            rtsm->quantiles[i] = 0.0;
        }
        rtsm->ps = ps;
    }

    for (size_t i = 0; i < rtsm->_n_buf; i++) {
        rtsm->_buf[i] = 0.0;
    }
    rtsm->n_quantiles = n_quantiles;
    
}

void RTStatsMonitor_del(RTStatsMonitor * rtsm) {
    if (rtsm->window_size) {
        STATS_FREE(rtsm->data);
    }

    if (rtsm->order_moments) {
        STATS_FREE(rtsm->moments);
    }

    if (rtsm->order_central_moments) {
        STATS_FREE(rtsm->central_moments);
    }

    if (rtsm->order_unbiased_central_moments) {
        STATS_FREE(rtsm->unbiased_central_moments);
    }

    if (rtsm->n_quantiles) {
        STATS_FREE(rtsm->quantiles);
        STATS_FREE(rtsm->ps);
    }

    STATS_FREE(rtsm->_buf);

    STATS_FREE(rtsm);
}

void RTStatsMonitor_process(RTStatsMonitor * rtsm, double * new_data) {
    double new_mean, old_mean;
    size_t offset;
    if (!rtsm->window_size) {
        for (size_t channel = 0; channel < rtsm->n_channels; channel++) {
            offset = channel * rtsm->order_moments;
            old_mean = rtsm->moments[offset];
            rt_moments(rtsm->moments + offset, rtsm->order_moments, rtsm->window_size, new_data[channel], NAN);
            new_mean = rtsm->moments[offset];
            offset = channel * (rtsm->order_central_moments - CM_ORDER_OFFSET);
            if (rtsm->order_central_moments) {
                rt_central_moments(rtsm->moments + offset, rtsm->order_central_moments, rtsm->window_size, new_mean, old_mean, new_data[channel], NAN);
            }
        }
    } else {
        size_t data_offset = (rtsm->_n_samples % rtsm->window_size) * rtsm->n_channels;
        
        for (size_t channel = 0; channel < rtsm->n_channels; channel++) {
            offset = channel * rtsm->order_moments;
            old_mean = rtsm->moments[offset];
            rt_moments(rtsm->moments + offset, rtsm->order_moments, rtsm->window_size, new_data[channel], (rtsm->_n_samples < rtsm->_window_size) ? NAN : rtsm->data[data_offset + channel]);
            new_mean = rtsm->moments[offset];
            offset = channel * (rtsm->order_central_moments - CM_ORDER_OFFSET);
            if (rtsm->order_central_moments) {
                rt_central_moments(rtsm->central_moments + offset, rtsm->order_central_moments, rtsm->window_size, new_mean, old_mean, new_data[channel], (rtsm->_n_samples < rtsm->_window_size) ? NAN : rtsm->data[data_offset + channel]);
            }
            if (rtsm->order_unbiased_central_moments) {
                rt_unbiased_central_moments(rtsm->unbiased_central_moments + channel * (rtsm->order_unbiased_central_moments - CM_ORDER_OFFSET), rtsm->central_moments + offset, rtsm->order_unbiased_central_moments, rtsm->window_size) {
            }
            rtsm->data[data_offset + channel] = new_data[channel];
        }
    }
    for (size_t channel = 0; channel < n_channel)
    // update all moments
    // process all central_moments
    // process all unbiased_central_moments
}

double RTStatsMonitor_get(RTStatsMonitor * rtsm, size_t channel, char * stat) {
    // parse 'stat' to find which statistic is requested and return the associated value for that (stat, channel) combination
    return 0.0;
}

// returns number of values written to dest. Dest must be allocated at least as much as rtsm->n_channels
size_t  RTStatsMonitor_getall(double * dest, RTStatsMonitor * rtsm, char * stat) {
    // parse 'stat' to find which statistic is requested and copy that memory into dest
    // return number of elements written to dest
    return 0;
}