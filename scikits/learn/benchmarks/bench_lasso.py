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

from bench_glm import make_data

def bench(clf, X_train, Y_train, X_test, Y_test):
    gc.collect()

    # start time
    tstart = time()
    clf = clf.fit(X_train, Y_train)
    delta = (time() - tstart)
    # stop time

    as_size = np.sum(np.abs(clf.coef_) > 0)
    print "active set size: %s (%s %%)" % (as_size, float(as_size) /
                                                            X_train.shape[1])
    return delta

def compute_bench(alpha, n_samples, n_features):

    def LassoFactory(alpha):
        return Lasso(alpha=alpha, fit_intercept=False)

    def LassoLARSFactory(alpha):
        return LassoLARS(alpha=alpha, normalize=False)
        # return LassoLARS(alpha=alpha, fit_intercept=False, normalize=False)

    lasso_results = []
    larslasso_results = []

    n_tests = 1000
    it = 0

    for ns in n_samples:
        for nf in n_features:
            it += 1
            print '============'
            print 'Iteration %s' % it
            print '============'
            k = nf // 10
            X, Y, X_test, Y_test, coef_ = make_data(
                n_samples=ns, n_tests=n_tests, n_features=nf,
                noise=0.1, k=k)

            X /= np.sqrt(np.sum(X**2, axis=0)) # Normalize data

            print "benching Lasso: "
            lasso_results.append(bench(LassoFactory(alpha),
                                                X, Y, X_test, Y_test))
            print "benching LassoLARS: "
            larslasso_results.append(bench(LassoLARSFactory(alpha),
                                                X, Y, X_test, Y_test))

    return lasso_results, larslasso_results

if __name__ == '__main__':
    from scikits.learn.glm import Lasso, LassoLARS
    import pylab as pl

    alpha = 0.01 # regularization parameter

    n_features = 500
    list_n_samples = range(500, 10001, 500);
    lasso_results, larslasso_results = compute_bench(alpha, list_n_samples,
                                                                [n_features])

    pl.close('all')
    pl.title('Lasso benchmark (%d features - alpha=%s)' % (n_features, alpha))
    pl.plot(list_n_samples, lasso_results, 'b-', label='Lasso')
    pl.plot(list_n_samples, larslasso_results,'r-', label='LassoLARS')
    pl.legend()
    pl.xlabel('number of samples')
    pl.ylabel('time (in seconds)')
    pl.show()

    n_samples = 500
    list_n_features = range(500, 3001, 500);
    lasso_results, larslasso_results = compute_bench(alpha, [n_samples],
                                                                list_n_features)

    pl.figure()
    pl.title('Lasso benchmark (%d samples - alpha=%s)' % (n_samples, alpha))
    pl.plot(list_n_features, lasso_results, 'b-', label='Lasso')
    pl.plot(list_n_features, larslasso_results,'r-', label='LassoLARS')
    pl.legend()
    pl.xlabel('number of features')
    pl.ylabel('time (in seconds)')
    pl.show()

