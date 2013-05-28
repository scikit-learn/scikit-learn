# Author: Mathieu Blondel <mathieu@mblondel.org>
# License: BSD 3 clause
import time

import pylab as pl

from sklearn.utils import check_random_state
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.metrics.pairwise import pairwise_kernels

def plot(func):
    random_state = check_random_state(0)
    one_core = []
    multi_core = []
    sample_sizes = range(1000, 6000, 1000)

    for n_samples in sample_sizes:
        X = random_state.rand(n_samples, 300)

        start = time.time()
        func(X, n_jobs=1)
        one_core.append(time.time() - start)

        start = time.time()
        func(X, n_jobs=-1)
        multi_core.append(time.time() - start)

    pl.figure('scikit-learn parallel %s benchmark results' % func.__name__)
    pl.plot(sample_sizes, one_core, label="one core")
    pl.plot(sample_sizes, multi_core, label="multi core")
    pl.xlabel('n_samples')
    pl.ylabel('time')
    pl.title('Parallel %s' % func.__name__)
    pl.legend()

def euclidean_distances(X, n_jobs):
    return pairwise_distances(X, metric="euclidean", n_jobs=n_jobs)

def rbf_kernels(X, n_jobs):
    return pairwise_kernels(X, metric="rbf", n_jobs=n_jobs, gamma=0.1)

plot(euclidean_distances)
plot(rbf_kernels)
pl.show()
