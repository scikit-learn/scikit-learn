"""
To run this, you'll need to have installed.

  * glmnet-python
  * scikit-learn (of course)

Does two benchmarks

First, we fix a training set and increase the number of
samples. Then we plot the computation time as function of
the number of samples.

In the second benchmark, we increase the number of dimensions of the
training set. Then we plot the computation time as function of
the number of dimensions.

In both cases, only 10% of the features are informative.
"""
import numpy as np
import gc
from time import time

# alpha = 1.0
alpha = 0.1
# alpha = 0.01


def rmse(a, b):
    return np.sqrt(np.mean((a - b) ** 2))


def make_data(n_samples=100, n_tests=100, n_features=100, k=10,
              noise=0.1, seed=0):

    # deterministic test
    np.random.seed(seed)

    # generate random input set
    X = np.random.randn(n_samples, n_features)
    X_test = np.random.randn(n_tests, n_features)

    # generate a ground truth model with only the first 10 features being non
    # zeros (the other features are not correlated to Y and should be ignored by
    # the L1 regularizer)
    coef_ = np.random.randn(n_features)
    coef_[k:] = 0.0

    # generate the ground truth Y from the reference model and X
    Y = np.dot(X, coef_)
    if noise > 0.0:
        Y += np.random.normal(scale=noise, size=Y.shape)

    Y_test = np.dot(X_test, coef_)
    if noise > 0.0:
        Y_test += np.random.normal(scale=noise, size=Y_test.shape)

    return X, Y, X_test, Y_test, coef_


def bench(factory, X, Y, X_test, Y_test, ref_coef):
    gc.collect()

    # start time
    tstart = time()
    clf = factory(alpha=alpha).fit(X, Y)
    delta = (time() - tstart)
    # stop time

    print "duration: %fms" % (delta * 1000)
    print "rmse: %f" % rmse(Y_test, clf.predict(X_test))
    print "mean coef abs diff: %f" % abs(ref_coef - clf.coef_.ravel()).mean()
    return delta


if __name__ == '__main__':
    from glmnet.elastic_net import Lasso as GlmnetLasso
    from scikits.learn.glm import Lasso as ScikitLasso
    # Delayed import of pylab
    import pylab as pl

    scikit_results = []
    glmnet_results = []
    n = 20
    step = 500
    n_features = 1000
    k = n_features / 10
    n_tests = 1000
    for i in range(1, n + 1):
        print '=================='
        print 'Iteration %s of %s' % (i, n)
        print '=================='
        X, Y, X_test, Y_test, coef_ = make_data(
            n_samples=(i * step), n_tests=n_tests, n_features=n_features,
            noise=0.1, k=k)

        print "benching scikit: "
        scikit_results.append(bench(ScikitLasso, X, Y, X_test, Y_test, coef_))
        print "benching glmnet: "
        glmnet_results.append(bench(GlmnetLasso, X, Y, X_test, Y_test, coef_))

    pl.clf()
    xx = range(0, n*step, step)
    pl.title('Lasso regression on sample dataset (%d features)' % n_features)
    pl.plot(xx, scikit_results, 'b-', label='scikit-learn')
    pl.plot(xx, glmnet_results,'r-', label='glmnet')
    pl.legend()
    pl.xlabel('number of samples to classify')
    pl.ylabel('time (in seconds)')
    pl.show()

    # now do a bench where the number of points is fixed
    # and the variable is the number of features

    scikit_results = []
    glmnet_results = []
    n = 20
    step = 100
    n_samples = 500

    for i in range(1, n + 1):
        print '=================='
        print 'Iteration %02d of %02d' % (i, n)
        print '=================='
        n_features = i * step
        k = n_features / 10
        X, Y, X_test, Y_test, coef_ = make_data(
            n_samples=n_samples, n_tests=n_tests, n_features=n_features,
            noise=0.1, k=k)

        print "benching scikit: "
        scikit_results.append(bench(ScikitLasso, X, Y, X_test, Y_test, coef_))
        print "benching glmnet: "
        glmnet_results.append(bench(GlmnetLasso, X, Y, X_test, Y_test, coef_))

    xx = np.arange(100, 100 + n * step, step)
    pl.figure()
    pl.title('Regression in high dimensional spaces (%d samples)' % n_samples)
    pl.plot(xx, scikit_results, 'b-', label='scikit-learn')
    pl.plot(xx, glmnet_results,'r-', label='glmnet')
    pl.legend()
    pl.xlabel('number of features')
    pl.ylabel('time (in seconds)')
    pl.axis('tight')
    pl.show()

