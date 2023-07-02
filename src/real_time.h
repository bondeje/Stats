#include <math.h>
#include "utils.h"

#ifndef REAL_TIME_H
#define REAL_TIME_H

/*
TODO: build monitor for tracking streaming updates
TODO: add functionality to track quantiles in real-time
TODO: add functionality to track min/max in real-time (requires tracking quantiles)
TODO: allow user to configure additional stats
*/

#define CM_ORDER_OFFSET 1
#define MAX_UNBIASED_CENTRAL_MOMENT 4

typedef struct RTStatsMonitor {
    double * data;
    double * moments;
    double * central_moments;
    double * unbiased_central_moments;
    double * ps;                // not fully supported yet
    double * quantiles;         // not fully supported yet
    double * _buf;
    size_t n_channels;
    size_t window_size;
    size_t n_quantiles;         // not fully supported yet
    size_t _n_samples;          // keep track of how many samples have been added
    size_t _n_buf;      // number of elements in buf. Limited to 
    unsigned short order_moments;
    unsigned short order_central_moments;
    unsigned short order_unbiased_central_moments;
} RTStatsMonitor;

/* 
to disable update of means in rt_variance and rt_covariance methods. If tracking
 more than 2 variables, this requires more maintenance but should be marginally 
 faster and consume less memory. If tracking <= 2 variables, this will slow things down
*/
// #define RT_NO_MEAN_UPDATE

/*
updates the average as a rolling average using O(1) algorithm.
assumes window_size points have already been accounted for in the avg
does not handle weighted averages
avg must not be NULL
WARNING: initialization must be handled separately/results may not make sense until at least window_size data have been added
*/ 
int rt_mean(double * avg, double new_data, size_t window_size, double old_data);

/*
updates the standard deviation as a rolling standard deviation using O(1) algorithm. 
also update the average
assumes window_size points have already been accounted for in the avg
does not handle weighted standard deviations
stddev must not be NULL
WARNING: initialization must be handled separately/results may not make sense until at least window_size data have been added
*/
int rt_variance(double * var, double * old_mean, double new_data, size_t window_size, double old_data);

int rt_covariance(double * cov_ab, double * old_mean_a, double * old_mean_b, double new_data_a, double new_data_b, size_t window_size, double old_data_a, double old_data_b);

// in lieu of rt_correlation, just use rt_covariance(a,b) / sqrt( rt_variance(a) * rt_variance(b) ). Since each of these calls is O(1), this is also O(1).

/* use weight-balanced tree with a deque for real-time quantile selection */
/*

typedef struct {
    WeightBalancedTree * wbt;   // self-balancing binary search tree for fast selection
    Deque * deq;                // deque of nodes in the wbt
    size_t * indexes;           // indexes to select
    size_t * n_indexes;         // number of indexes to select
    size_t window;              // size of window.
    select_method method;       // method to use to select the value. Use an enum value
} RT_Selector;

RT_Selector * RT_Selector_new(double * quantiles, size_t n_quantiles, size_t window, select_method method);

void RT_Select_init(RT_Selector * rts);

void RT_Selector_del(RT_Selector * rts);

// 1a) if rts is not "full": add value to wbt, push new node to deq, update indexes
// 1b) else: pop value from deq, reinsert node with new value back into wbt
int RT_Selector_push(RT_Selector * rts, double value);

int RT_Selector_select(double * quantile_values, RT_Selector * rts); // if multiple calls are ever need, might be more efficient to call this internally after push and just store the result

// wrapper to push data and select. quantile_values must be at least as large as RT_Selector.n_indexes
int rt_quantile(double * quantile_values, RT_Selector * rts, double new_data);

with an order statistics tree, calculating min, max, range in real time become trivial
*/

#endif // REAL_TIME_H