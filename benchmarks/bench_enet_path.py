"""
Benchmark of enet_path against enet_path with strong-rules.

First, we fix a training set and increase the number of
samples. Then we plot the computation time as function of
the number of samples.

In the second benchmark, we increase the number of dimensions of the
training set. Then we plot the computation time as function of
the number of dimensions.

In both cases, only 10% of the features are informative.
"""
import gc
from time import time
import numpy as np

from sklearn.datasets.samples_generator import make_regression
from numpy.testing import assert_almost_equal


def compute_bench(n_samples, n_features, precompute):

    enet_path_sr_results = []
    enet_path_results = []
    enet_path_gram_results = []

    n_test_samples = 0
    it = 0

    for ns in n_samples:
        for nf in n_features:
            it += 1
            print '=================='
            print 'Iteration %s of %s' % (it, max(len(n_samples),
                                          len(n_features)))
            print '=================='
            n_informative = nf // 10
            X, y = make_regression(n_samples=ns, n_features=nf,
                                          n_informative=n_informative,
                                          noise=0.1)

            X /= np.sqrt(np.sum(X ** 2, axis=0))  # Normalize data

            gc.collect()
            print "enet_path with strong-rules"
            stime = time()
            enet_path(X, y, n_alphas=10, precompute=False,
                            fit_intercept=False, use_strong_rule=True)

            enet_path_sr_results.append(time() - stime)

            gc.collect()
            print "enet_path"
            stime = time()
            enet_path(X, y, n_alphas=10, precompute=True,
                                             fit_intercept=False)
            enet_path_results.append(time() - stime)

            gc.collect()
            print "enet_path (with Gram)"
            stime = time()
            enet_path(X, y, n_alphas=10, precompute=True,
                                             fit_intercept=False)
            enet_path_gram_results.append(time() - stime)

    return enet_path_sr_results, enet_path_results, enet_path_gram_results


if __name__ == '__main__':
    from sklearn.linear_model.coordinate_descent import enet_path
    import pylab as pl

    n_features = 10
#    list_n_samples = np.linspace(100, 1000000, 5).astype(np.int)
    list_n_samples = np.linspace(100, 10000, 5).astype(np.int)
    enet_path_sr_results, enet_path_results, enet_path_gram_results=  \
            compute_bench(list_n_samples, [n_features], precompute=True)

    pl.clf()
    pl.subplot(211)
    pl.plot(list_n_samples, enet_path_sr_results, 'b-',
                            label='enet_path_sr_results')
    pl.plot(list_n_samples, enet_path_results, 'r-',
                            label='enet_path_results')
    pl.plot(list_n_samples, enet_path_gram_results, 'r-',
                                 label='enet_path_gram_results')

    pl.title('Enet benchmark (%d features )' % (n_features))
    pl.legend(loc='upper left')
    pl.xlabel('number of samples')
    pl.ylabel('time (in seconds)')
    pl.axis('tight')

    n_samples = 2000
#    list_n_features = np.linspace(500, 3000, 5).astype(np.int)
    list_n_features = np.linspace(500, 1000, 5).astype(np.int)
    enet_path_sr_results, enet_path_results, enet_path_gram_results = \
            compute_bench([n_samples], list_n_features, precompute=False)
    pl.subplot(212)
    pl.plot(list_n_features, enet_path_sr_results, 'b-',
                            label='enet_path_sr_results')
    pl.plot(list_n_features, enet_path_results, 'r-',
                                 label='enet_path_results')
    pl.plot(list_n_features, enet_path_gram_results, 'r-',
                                 label='enet_path_gram_results')

    pl.title('Enet benchmark (%d samples )' % (n_samples))
    pl.legend(loc='upper left')
    pl.xlabel('number of features')
    pl.ylabel('time (in seconds)')

    pl.axis('tight')
    pl.show()
