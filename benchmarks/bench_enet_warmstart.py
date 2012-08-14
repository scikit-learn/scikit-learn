"""
Benchmarks of enet_coordinate_descent vs. enet_coordinate_descent
using the true solution as warm-start

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


def compute_bench(alpha, rho, n_samples, n_features, precompute):

    cold_start_results = []
    warm_start_results = []
    warm_start_residual_results = []

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
            X, y, coef_ = make_regression(n_samples=ns, n_features=nf,
                                          n_informative=n_informative,
                                          noise=0.1, coef=True)

            X /= np.sqrt(np.sum(X ** 2, axis=0))  # Normalize data

            l1_reg = alpha * rho * ns
            l2_reg = alpha * (1.0 - rho) * ns
            X = np.asfortranarray(X)
            w = np.zeros(nf)
            Xy = np.dot(X.T, y)

            gc.collect()
            print "enet fit"
            stime = time()
            R = y - np.dot(X, w)
            enet_coordinate_descent(w, l1_reg, l2_reg,
                        X, y, max_iter=10000, tol=1e-9, positive=False, R=R)
            cold_start_results.append(time() - stime)

            gc.collect()
            warm_start = w.copy()
            print "warmstart with solution, enet fit"
            stime = time()

            enet_coordinate_descent(warm_start, l1_reg, l2_reg,
                    X, y, max_iter=10000, tol=1e-9, positive=False)
            warm_start_results.append(time() - stime)

            gc.collect()
            warm_start = w.copy()
            print "warmstart with solution & residual, enet fit"
            stime = time()

            enet_coordinate_descent(warm_start, l1_reg, l2_reg,
                    X, y, max_iter=10000, tol=1e-9, positive=False, R=R)
            warm_start_residual_results.append(time() - stime)

            assert_almost_equal(w, warm_start)

    return cold_start_results, warm_start_results, warm_start_residual_results


if __name__ == '__main__':
    from sklearn.linear_model.coordinate_descent import ElasticNet
    from sklearn.linear_model.cd_fast import enet_coordinate_descent
    import pylab as pl

    alpha = 0.5  # regularization parameter
    rho = 0.95

    n_features = 10
    list_n_samples = np.linspace(100, 1000000, 5).astype(np.int)
    cold_start_results, warm_start_results, warm_start_residual_results =  \
            compute_bench(alpha, rho,
                         list_n_samples, [n_features], precompute=True)

    pl.clf()
    pl.subplot(211)
    pl.plot(list_n_samples, cold_start_results, 'b-',
                            label='cold_start_results')
    pl.plot(list_n_samples, warm_start_results, 'r-',
                            label='warm_start_results')
    pl.plot(list_n_samples, warm_start_residual_results, 'g-',
                                 label='warm_start_residual_results')
    pl.title('Enet benchmark (%d features - alpha=%s)' % (n_features, alpha))
    pl.legend(loc='upper left')
    pl.xlabel('number of samples')
    pl.ylabel('time (in seconds)')
    pl.axis('tight')

    n_samples = 2000
    list_n_features = np.linspace(500, 3000, 5).astype(np.int)
    cold_start_results, warm_start_results, warm_start_residual_results = \
            compute_bench(alpha, rho,
                            [n_samples], list_n_features, precompute=False)
    pl.subplot(212)
    pl.plot(list_n_features, cold_start_results, 'b-', label='cold_start_results')
    pl.plot(list_n_features, warm_start_results, 'r-',
                                 label='warm_start_results')
    pl.plot(list_n_features, warm_start_residual_results, 'g-',
                                 label='warm_start_residual_results')

    pl.title('Enet benchmark (%d samples - alpha=%s)' % (n_samples, alpha))
    pl.legend(loc='upper left')
    pl.xlabel('number of features')
    pl.ylabel('time (in seconds)')

    pl.axis('tight')
    pl.show()
