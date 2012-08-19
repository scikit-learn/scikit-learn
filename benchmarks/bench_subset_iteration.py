"""
Benchmarks of ElasticNet vs ElasticNet with strong rules

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


def compute_bench(alpha, rho, n_samples, n_features, precompute):

    full_set_results = []
    active_set_results = []
    reduced_results = []

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

            l1_reg = alpha * rho * ns
            l2_reg = alpha * (1.0 - rho) * ns
            X = np.asfortranarray(X)
            w = np.zeros(nf)
            Xy = np.dot(X.T, y)
            norm_cols_X = (X**2).sum(axis=0)

            gc.collect()
            print "iterate on full set of features"
            stime = time()

            enet_coordinate_descent(w, l1_reg, l2_reg,
                            X, y, max_iter=10000, tol=1e-9, positive=False)
            full_set_results.append(time() - stime)

            print "# nnz: " + str(len(w.nonzero()[0]))
            if len(w.nonzero()[0]) > 0:
                active_set = w.nonzero()[0]
                X_red = X[:, active_set]
                X_red = np.asfortranarray(X_red)
                w_red = w[active_set]
            else:
                active_set = None
                X_red = X
                w_red = w.copy()

            w = np.zeros(nf)
            R = y - np.dot(X, w)
            gc.collect()
            print "iterate only on active features"
            stime = time()

            enet_coordinate_descent(w, l1_reg, l2_reg,
                    X, y, max_iter=10000, tol=1e-9, positive=False,
                    iter_set=active_set, norm_cols_X=norm_cols_X, R=R)
            active_set_results.append(time() - stime)

            gc.collect()
            print "iterate on reduced problem"
            stime = time()
            enet_coordinate_descent(w_red, l1_reg, l2_reg,
                            X_red, y, max_iter=10000, tol=1e-9, positive=False)
            reduced_results.append(time() - stime)

    return full_set_results, active_set_results, reduced_results


if __name__ == '__main__':
    from sklearn.linear_model.coordinate_descent import ElasticNet
    from sklearn.linear_model.cd_fast import enet_coordinate_descent
    import pylab as pl

    alpha = 0.7  # regularization parameter
    rho = 0.95

    n_features = 10
    list_n_samples = np.linspace(100, 1000000, 5).astype(np.int)
    full_set_results, active_set_results, reduced_results = compute_bench(alpha, rho,
                         list_n_samples, [n_features], precompute=True)

    pl.clf()
    pl.subplot(211)
    pl.plot(list_n_samples, full_set_results, 'b-',
                            label='full_set_results')
    pl.plot(list_n_samples, active_set_results, 'r-',
                            label='active_set_results')
    pl.plot(list_n_samples, reduced_results, 'g-',
                            label='reduced_results')
    pl.title('Enet benchmark (%d features - alpha=%s)' % (n_features, alpha))
    pl.legend(loc='upper left')
    pl.xlabel('number of samples')
    pl.ylabel('time (in seconds)')
    pl.axis('tight')

    n_samples = 2000
    list_n_features = np.linspace(500, 3000, 5).astype(np.int)
    full_set_results, active_set_results, reduced_results = compute_bench(alpha, rho, 
                            [n_samples], list_n_features, precompute=False)
    pl.subplot(212)
    pl.plot(list_n_features, full_set_results, 'b-', label='full_set_results')
    pl.plot(list_n_features, active_set_results, 'r-',
                                 label='active_set_results')
    pl.plot(list_n_features, reduced_results, 'g-',
                            label='reduced_results')
    pl.title('Enet benchmark (%d samples - alpha=%s)' % (n_samples, alpha))
    pl.legend(loc='upper left')
    pl.xlabel('number of features')
    pl.ylabel('time (in seconds)')

    pl.axis('tight')
    pl.show()
