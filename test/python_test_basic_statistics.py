import numpy as np
import matplotlib.pyplot as plt

from math import comb

fignum = 0

nCr = comb

def _accumulate_sample_moments_weights(sample_moments, samples, n, min_moment, max_moment, weights = None, weight_moments = None):
    if min_moment < 1:
        min_moment = 1
    
    for j in range(max_moment-min_moment+1):
        sample_moments[j] = 0.0
    if weight_moments:
        for j in range(max_moment-min_moment+1):
            weight_moments[j] = 0.0
    
    weight_sum = 0.0
    for i in range(n):
        weight = weights[i] if weights is not None else 1.0
        weight_sum += weight
        poly_ = weight * samples[i]**(min_moment-1)
        for j in range(max_moment-min_moment+1):
            poly_ *= samples[i]
            sample_moments[j] += poly_
        if (weight_moments):
            poly_ = pow(weight, min_moment-1)
            for j in range(max_moment-min_moment+1):
                poly_ *= weight
                weight_moments[j] += poly_

    for j in range(max_moment-min_moment+1):
        sample_moments[j] /= weight_sum

def sample_moment(samples, n, order, weights = None):
    sample_moment = [0.0]

    _accumulate_sample_moments_weights(sample_moment, samples, n, order, order, weights, None)

    return sample_moment[0]

def _unbiased_second_central_sample_moment(samples, n, weights = None):
    order = 2
    s_acc = [0 for i in range(order)]
    w_acc = [0 for i in range(order)]

    out = 0.0
    _accumulate_sample_moments_weights(s_acc, samples, n, 1, order, weights, w_acc)
    out += s_acc[1] - s_acc[0]**2
    out *= w_acc[0]**2/(w_acc[0]**2 - w_acc[1])

    return out

def _unbiased_third_central_sample_moment(samples, n, weights = None):
    order = 3
    s_acc = [0 for i in range(order)]
    w_acc = [0 for i in range(order)]

    out = 0.0
    _accumulate_sample_moments_weights(s_acc, samples, n, 1, order, weights, w_acc)
    out += s_acc[2] - 3 * s_acc[1] * s_acc[0] + 2 * s_acc[0]**3
    out *= w_acc[0]**3/(w_acc[0]**3 - 3 * w_acc[0] * w_acc[1] + 2 * w_acc[2])

    # out = n^3/(n^3 - 3n^2 + 2n) = n^2/((n-1)(n-2))

    return out

def _unbiased_fourth_central_sample_moment(samples, n, weights = None):
    order = 4
    s_acc = [0 for i in range(order)]
    w_acc = [0 for i in range(order)]

    out = 0.0

    _accumulate_sample_moments_weights(s_acc, samples, n, 1, order, weights, w_acc)
    """ DOES NOT WORK
    out += (w_acc[0]**4 - 3*w_acc[1]*w_acc[0]**2 + 2*w_acc[2]*w_acc[0]+2*w_acc[1]**2 - 3*w_acc[3]) * (s_acc[3]-4*s_acc[2]*s_acc[0])
    out -= 3*(2*w_acc[2]*w_acc[0]**2 - 2*w_acc[2]*w_acc[0] - 3*w_acc[1]**2 + 3*w_acc[3])*s_acc[1]**2
    out += 3*s_acc[0]**2*w_acc[0]**2*(w_acc[0]**2 - w_acc[1]) * (2*s_acc[1] - s_acc[0]**2)
    
    out *= w_acc[0]**2 / (w_acc[0]**6 - 7*w_acc[1]*w_acc[0]**4 + 8*w_acc[2]*w_acc[0]**3 + 9*w_acc[1]**2*w_acc[0]**2 - 6*w_acc[3]*w_acc[0]**2 - 8*w_acc[2]*w_acc[1]*w_acc[0] - 3*w_acc[1]**3 + 6*w_acc[3]*w_acc[1])
    """

    """ WORKS
    #out += (n/(n-1)) * (n**2 - 2*n + 3)/((n-2)*(n-3)) * (s_acc[3] - 4*s_acc[2]*s_acc[0] + 6*s_acc[1]*s_acc[0]**2 - 3*s_acc[0]**4) - 3*n/(n-1)*(2*n-3)/((n-2)*(n-3))*(s_acc[1]-s_acc[0]**2)**2
    out += (w_acc[0]**4 - 3*w_acc[1]*w_acc[0]**2 + 2*w_acc[2]*w_acc[0] + 3*w_acc[1]**2 - 3*w_acc[3])*(s_acc[3] - 4*s_acc[2]*s_acc[0] + 6*s_acc[1]*s_acc[0]**2 - 3*s_acc[0]**4)
    out -= 3*(2*w_acc[1]*w_acc[0]**2 - 2*w_acc[2]*w_acc[0] - 3*w_acc[1]**2 + 3*w_acc[3])*(s_acc[1]-s_acc[0]**2)**2
    out *= w_acc[0]**2/(w_acc[0]**6 - 7*w_acc[1]*w_acc[0]**4 + 8*w_acc[2]*w_acc[0]**3 + 9*w_acc[1]**2*w_acc[0]**2 - 6*w_acc[3]*w_acc[0]**2 - 8*w_acc[2]*w_acc[1]*w_acc[0] - 3*w_acc[1]**3 + 6*w_acc[3]*w_acc[1])
    """

    #"""
    # testing
    out += (w_acc[0]**4 - 3*w_acc[1]*w_acc[0]**2 + 2*w_acc[2]*w_acc[0] + 3*w_acc[1]**2 - 3*w_acc[3])*(s_acc[3] - 4*s_acc[2]*s_acc[0])
    out -= 3*(2*w_acc[1]*w_acc[0]**2 - 2*w_acc[2]*w_acc[0] - 3*w_acc[1]**2 + 3*w_acc[3])*s_acc[1]**2
    out += 3*(w_acc[0]**2 - w_acc[1])*w_acc[0]**2*s_acc[0]**2*(2*s_acc[1] - s_acc[0]**2)
    out *= w_acc[0]**2/(w_acc[0]**6 - 7*w_acc[1]*w_acc[0]**4 + 8*w_acc[2]*w_acc[0]**3 + 9*w_acc[1]**2*w_acc[0]**2 - 6*w_acc[3]*w_acc[0]**2 - 8*w_acc[2]*w_acc[1]*w_acc[0] - 3*w_acc[1]**3 + 6*w_acc[3]*w_acc[1])
    #"""
    return out

def _biased_second_central_sample_moment(samples, n, weights = None):
    order = 2
    s_acc = [0 for i in range(order)]
    w_acc = [0 for i in range(order)]

    out = 0.0
    _accumulate_sample_moments_weights(s_acc, samples, n, 1, order, weights, w_acc)
    out += s_acc[1] - s_acc[0]**2

    return out

def central_sample_moment_biased(samples, n, order, weights = None):
    s_acc = [0 for it in range(order)]

    _accumulate_sample_moments_weights(s_acc, samples, n, 1, order, weights, None)
    out = 0.0
    for j in range(order):
        out += (-1)**(order-j+1) * (nCr(order, j+1) - (1 if j==0 else 0)) * s_acc[j] * s_acc[0]**(order-j-1)

    return out

#%% test parameter estimations

# Normal distribution
loc = 3
scale = 1.2
N = 10
Nrepeat = 1000

unbiased1 = []

biased2 = []
unbiased2 = []

biased3 = []
unbiased3 = []

biased4 = []
unbiased4 = []

for i in range(Nrepeat):
    samples = np.random.normal(loc, scale, N)
    #order = 4
    #sample_weights = [0 for i in range(order)]

    unbiased1.append(sample_moment(samples, N, 1))

    #_accumulate_sample_moments_weights(sample_weights, samples, N, 0, order)
    #print("sample_weights:\n",sample_weights)
    biased2.append(central_sample_moment_biased(samples, N, 2))
    unbiased2.append(_unbiased_second_central_sample_moment(samples, N))

    biased3.append(central_sample_moment_biased(samples, N, 3))
    unbiased3.append(_unbiased_third_central_sample_moment(samples, N))

    biased4.append(central_sample_moment_biased(samples, N, 4))
    unbiased4.append(_unbiased_fourth_central_sample_moment(samples, N))

print(f"average estimated mean: {np.mean(unbiased1)}, std dev: {np.std(unbiased1, ddof=1)}")
    
print(f"average biased estimated variance: {np.mean(biased2)}, std dev: {np.std(biased2, ddof=1)}")
print(f"average unbiased estimated variance: {np.mean(unbiased2)}, std dev: {np.std(unbiased2, ddof=1)}")

print(f"average biased estimated third moment: {np.mean(biased3)}, std dev: {np.std(biased3, ddof=1)}")
print(f"average unbiased estimated third moment: {np.mean(unbiased3)}, std dev: {np.std(unbiased3, ddof=1)}")

print(f"average biased estimated fourth moment: {np.mean(biased4)}, std dev: {np.std(biased4, ddof=1)}")
print(f"average unbiased estimated fourth moment: {np.mean(unbiased4)}, std dev: {np.std(unbiased4, ddof=1)}")

"""
samples = [3,7,9,3,3,8,9,8,5,9,4,5,4,4,8,5,1,3,1,4,6,9,2,2,3,9,9,5,7,5,1,4,2,4,2,1,2,1,3,2,8,1,4,8]
for N in range(len(samples)):
    if N < 4:
        print(0.0)
    else:
        print(_unbiased_fourth_central_sample_moment(samples, N))
"""

plt.figure(fignum)
fignum += 1
plt.hist(samples, bins=[it*.5-10 for it in range(41)])
plt.title(f"normal distribution - N({loc},{scale})\nVariance={scale**2}")
plt.show()
