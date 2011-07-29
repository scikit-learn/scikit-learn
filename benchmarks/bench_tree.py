"""
To run this, you'll need to have installed.

  * scikit-learn 

Does two benchmarks

First, we fix a training set, increase the number of
samples to classify and plot number of classified samples as a
function of time.

In the second benchmark, we increase the number of dimensions of the
training set, classify a sample and plot the time taken as a function
of the number of dimensions.
"""
import numpy as np
import pylab as pl
import gc
from datetime import datetime

# to store the results
scikit_classifier_results = []
scikit_regressor_results = []

mu_second = 0.0 + 10**6 # number of microseconds in a second


def bench_scikit_tree_classifier(X, Y):
    """
    bench with scikit-learn decision tree classifier
    """
    import scikits.learn
    from scikits.learn.tree_model import DecisionTreeClassifier

    gc.collect()

    # start time
    tstart = datetime.now()
    clf = DecisionTreeClassifier()
    clf.fit(X, Y).predict(X)
    delta = (datetime.now() - tstart)
    # stop time

    scikit_classifier_results.append(delta.seconds + delta.microseconds/mu_second)

def bench_scikit_tree_regressor(X, Y):
    """
    bench with scikit-learn decision tree regressor
    """
    import scikits.learn
    from scikits.learn.tree_model import DecisionTreeRegressor

    gc.collect()

    # start time
    tstart = datetime.now()
    clf = DecisionTreeRegressor()
    clf.fit(X, Y).predict(X)
    delta = (datetime.now() - tstart)
    # stop time

    scikit_regressor_results.append(delta.seconds + delta.microseconds/mu_second)

if __name__ == '__main__':

    n = 10
    step = 50
    n_samples = 50
    dim = 10
    K = 10 
    for i in range(n):
        print '============================================'
        print 'Entering iteration %s of %s' % (i, n)
        print '============================================'
        n_samples += step
        X = np.random.randn(n_samples, dim)
        Y = np.random.random_integers(0, K, (n_samples,))
        bench_scikit_tree_classifier(X, Y)
        Y = np.random.randn(n_samples)
        bench_scikit_tree_regressor(X, Y)


    import pylab as pl
    xx = range(0, n*step, step)
    pl.figure(1)
    pl.subplot(211)
    pl.title('Learning with varying number of samples')
    pl.plot(xx, scikit_classifier_results, 'g-', label='classification')
    pl.plot(xx, scikit_regressor_results, 'r-', label='regression')
    pl.legend()
    pl.xlabel('number of samples')
    pl.ylabel('time (in microseconds)')


    # now do a bench where the number of points is fixed
    # and the variable is the number of dimensions
    from scikits.learn.datasets.samples_generator import friedman, \
                                                         sparse_uncorrelated

    scikit_classifier_results = []
    scikit_regressor_results = []
    n = 10
    step = 50
    start_dim = 50
    K = 10 



    print '============================================'
    print 'Warning: this is going to take a looong time'
    print '============================================'

    dim = start_dim
    for i in range(0, n):
        print '============================================'
        print 'Entering iteration %s of %s' % (i, n)
        print '============================================'
        dim += step
        X = np.random.randn(100, dim)
        Y = np.random.random_integers(0, K, (100,))
        bench_scikit_tree_classifier(X, Y)
        Y = np.random.randn(100)
        bench_scikit_tree_regressor(X, Y)

    xx = np.arange(start_dim, start_dim+n*step, step)
    pl.subplot(212)
    pl.title('Learning in high dimensional spaces')
    pl.plot(xx, scikit_classifier_results, 'g-', label='classification')
    pl.plot(xx, scikit_regressor_results, 'r-', label='regression')
    pl.legend()
    pl.xlabel('number of dimensions')
    pl.ylabel('time (in seconds)')
    pl.axis('tight')
    pl.show()
