import os
import random # using Mersenne twister

seed = 19551105
start = 1
stop = 10
N = 1000
output_file_dir = os.path.split(__file__)[0]
output_file = 'samples.csv'

random.seed(19551105)

def a(samples, n, order, weights=None):
    if weights:
        norm = sum(weights[:n])
    else:
        norm = n
        weights = [1 for i in range(norm)]
    out = 0.0
    for i in range(n):
        out += weights[i] * samples[i]**order
    return out/norm

def m(samples, n, order, weights=None):
    smean = a(samples, n, 1, weights)
    if weights:
        norm = sum(weights[:n])
    else:
        norm = n
        weights = [1 for i in range(norm)]
    out = 0.0
    for i in range(n):
        out += weights[i] * (samples[i] - smean)**order
    return out/norm

with open(os.path.join(output_file_dir, output_file), 'w') as file_out:
    file_out.write('test\n')
    samples = [random.randrange(start, stop) for i in range(N)]
    file_out.write(','.join(['samples'] + [str(sample) for sample in samples]) + '\n')
    weights = [random.randrange(start, stop) for i in range(N)]
    file_out.write(','.join(['weights'] + [str(weight) for weight in weights]) + '\n')

    for order in range(1,5):
        ap = []
        mp = []
        for i in range(N):
            ap.append(a(samples, i+1, order, weights))
            mp.append(m(samples, i+1, order, weights))
            
        file_out.write(','.join(['a' + str(order)] + [str(f) for f in ap]) + '\n')
        file_out.write(','.join(['m' + str(order)] + [str(f) for f in mp]) + '\n')

    # test Mp unweighted
    order = 2
    Mp = []
    for i in range(N):
        n = i+1
        if (n < order):
            Mp.append(0.0)
        else:
            mp = m(samples, n, order)
            Mp.append(mp*n/(n-1))
    file_out.write(','.join(['M' + str(order)] + [str(f) for f in Mp]) + '\n')

    order = 3
    Mp = []
    for i in range(N):
        n = i+1
        if (n < order):
            Mp.append(0.0)
        else:
            mp = m(samples, n, order)
            Mp.append(mp*(n/(n-1))*(n/(n-2)))
    file_out.write(','.join(['M' + str(order)] + [str(f) for f in Mp]) + '\n')

    order = 4
    Mp = []
    for i in range(N):
        n = i+1
        if (n < order):
            Mp.append(0.0)
        else:
            m4 = m(samples, n, 4)
            m2 = m(samples, n, 2)
            Mp.append((n/(n-1)) * (n**2 - 2*n + 3)/((n-2)*(n-3)) * m4 - 3*n/(n-1)*(2*n-3)/((n-2)*(n-3))*m2**2)
    file_out.write(','.join(['M' + str(order)] + [str(f) for f in Mp]) + '\n')

#with open('./data')