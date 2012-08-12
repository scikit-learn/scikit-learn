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

    enet_results = []
    enet_strong_rules_results = []

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
            X, Y, coef_ = make_regression(n_samples=ns, n_features=nf,
                                          n_informative=n_informative,
                                          noise=0.1, coef=True)

            X /= np.sqrt(np.sum(X ** 2, axis=0))  # Normalize data

            gc.collect()
            print "- benching ElasticNet"
            clf = ElasticNet(alpha=alpha, rho=rho, precompute=False, 
                             fit_intercept=False)
            tstart = time()
            clf.fit(X, Y)
            enet_results.append(time() - tstart)

            gc.collect()
            print "- benching ElasticNet with strong rules"
            clf = ElasticNet(alpha=alpha, rho=rho, precompute=False,
                          fit_intercept=False, use_strong_rule=True)
            tstart = time()
            clf.fit(X, Y)
            enet_strong_rules_results.append(time() - tstart)

    return enet_results, enet_strong_rules_results


if __name__ == '__main__':
    from sklearn.linear_model.coordinate_descent import ElasticNet
    import pylab as pl

    alpha = 100.0  # regularization parameter
    rho = 0.90

    n_features = 10
    list_n_samples = np.linspace(100, 1000000, 5).astype(np.int)
    enet_results, enet_strong_rules_results = compute_bench(alpha, rho,
                         list_n_samples, [n_features], precompute=True)

    pl.clf()
    pl.subplot(211)
    pl.plot(list_n_samples, enet_results, 'b-',
                            label='ElasticNet')
    pl.plot(list_n_samples, enet_strong_rules_results, 'r-',
                            label='ElasticNet with strong rules')
    pl.title('Enet benchmark (%d features - alpha=%s)' % (n_features, alpha))
    pl.legend(loc='upper left')
    pl.xlabel('number of samples')
    pl.ylabel('time (in seconds)')
    pl.axis('tight')

    n_samples = 2000
    list_n_features = np.linspace(500, 3000, 5).astype(np.int)
    enet_results, enet_strong_rules_results = compute_bench(alpha, rho, 
                            [n_samples], list_n_features, precompute=False)
    pl.subplot(212)
    pl.plot(list_n_features, enet_results, 'b-', label='ElasticNet')
    pl.plot(list_n_features, enet_strong_rules_results, 'r-',
                                 label='ElasticNet with strong rules')
    pl.title('Enet benchmark (%d samples - alpha=%s)' % (n_samples, alpha))
    pl.legend(loc='upper left')
    pl.xlabel('number of features')
    pl.ylabel('time (in seconds)')
    pl.axis('tight')
    pl.show()
