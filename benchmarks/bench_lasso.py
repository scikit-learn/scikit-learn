"""
Benchmarks of Lasso vs LassoLARS

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

from scikits.learn.datasets.samples_generator import make_regression_dataset


def compute_bench(alpha, n_samples, n_features, precompute):

    lasso_results = []
    larslasso_results = []

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
            X, Y, _, _, coef = make_regression_dataset(
                n_train_samples=ns, n_test_samples=n_test_samples,
                n_features=nf, noise=0.1, n_informative = n_informative)

            X /= np.sqrt(np.sum(X**2, axis=0)) # Normalize data

            gc.collect()
            print "- benching Lasso"
            clf = Lasso(alpha=alpha, fit_intercept=False)
            tstart = time()
            clf.fit(X, Y, precompute=precompute)
            lasso_results.append(time() - tstart)

            gc.collect()
            print "- benching LassoLARS"
            clf = LassoLARS(alpha=alpha, fit_intercept=False)
            tstart = time()
            clf.fit(X, Y, normalize=False, precompute=precompute)
            larslasso_results.append(time() - tstart)

    return lasso_results, larslasso_results


if __name__ == '__main__':
    from scikits.learn.linear_model import Lasso, LassoLARS
    import pylab as pl

    alpha = 0.01 # regularization parameter

    n_features = 10
    list_n_samples = np.linspace(100, 1000000, 5).astype(np.int)
    lasso_results, larslasso_results = compute_bench(alpha, list_n_samples,
                                            [n_features], precompute=True)

    pl.clf()
    pl.subplot(211)
    pl.plot(list_n_samples, lasso_results, 'b-',
                            label='Lasso (with precomputed Gram matrix)')
    pl.plot(list_n_samples, larslasso_results, 'r-',
                            label='LassoLARS (with precomputed Gram matrix)')
    pl.title('Lasso benchmark (%d features - alpha=%s)' % (n_features, alpha))
    pl.legend(loc='upper left')
    pl.xlabel('number of samples')
    pl.ylabel('time (in seconds)')
    pl.axis('tight')

    n_samples = 2000
    list_n_features = np.linspace(500, 3000, 5).astype(np.int)
    lasso_results, larslasso_results = compute_bench(alpha, [n_samples],
                                           list_n_features, precompute=False)
    pl.subplot(212)
    pl.plot(list_n_features, lasso_results, 'b-', label='Lasso')
    pl.plot(list_n_features, larslasso_results, 'r-', label='LassoLARS')
    pl.title('Lasso benchmark (%d samples - alpha=%s)' % (n_samples, alpha))
    pl.legend(loc='upper left')
    pl.xlabel('number of features')
    pl.ylabel('time (in seconds)')
    pl.axis('tight')
    pl.show()
