"""
To run this, you'll need to have installed.

  * pymvpa
  * libsvm and it's python bindings
  * scikit-learn (of course)

Does two benchmarks

First, we fix a training set, increase the number of
samples to classify and plot number of classified samples as a
function of time.

In the second benchmark, we increase the number of dimensions of the
training set, classify a sample and plot the time taken as a function of the number of dimensions.
"""
import numpy as np
import pylab as pl
import gc
from datetime import datetime

# to store the results
scikit_results = []
svm_results = []
mvpa_results = []

mu_second = 0.0 + 10**6 # number of microseconds in a second

def bench_scikit(X, Y):
    """
    bench with scikit-learn bindings on libsvm
    """
    import scikits.learn
    from scikits.learn.svm import SVC

    gc.collect()

    # start time
    tstart = datetime.now()
    clf = SVC(kernel='rbf')
    clf.fit(X, Y).predict(X)
    delta = (datetime.now() - tstart)
    # stop time

    scikit_results.append(delta.seconds + delta.microseconds/mu_second)

def bench_svm(X, Y):
    """
    bench with swig-generated wrappers that come with libsvm
    """

    import svm

    X1 = X.tolist()
    Y1 = Y.tolist()

    gc.collect()

    # start time
    tstart = datetime.now()
    problem = svm.svm_problem(Y1, X1)
    param = svm.svm_parameter(svm_type=0, kernel_type=2)
    model = svm.svm_model(problem, param)
    for i in X1:
        model.predict(i)
    delta = (datetime.now() - tstart)
    # stop time
    svm_results.append(delta.seconds + delta.microseconds/mu_second)

def bench_pymvpa(X, Y):
    """
    bench with pymvpa (by default uses a custom swig-generated wrapper
    around libsvm)
    """
    from mvpa.datasets import Dataset
    from mvpa.clfs import svm

    gc.collect()

    # start time
    tstart = datetime.now()
    data = Dataset(samples=X, labels=Y)
    clf = svm.RbfCSVMC(C=1.)
    clf.train(data)
    Z = clf.predict(X)
    delta = (datetime.now() - tstart)

    # stop time
    mvpa_results.append(delta.seconds + delta.microseconds/mu_second)

if __name__ == '__main__':

    n = 5
    step = 100
    n_samples = 200
    dim = 200
    for i in range(n):
        print '============================================'
        print 'Entering iteration %s of %s' % (i, n)
        print '============================================'
        n_samples += step
        X = np.random.randn(n_samples, dim)
        Y = np.random.randn(n_samples)
        bench_scikit(X, Y)
        bench_pymvpa(X, Y)
        bench_svm(X, Y)

    import pylab as pl
    xx = range(0, n*step, step)
    pl.figure(1)
    pl.subplot(211)
    pl.title('SVM with varying number of samples')
    pl.plot(xx, mvpa_results, 'g-', label='pymvpa')
    pl.plot(xx, svm_results,'r-', label='libsvm-swig')
    pl.plot(xx, scikit_results, 'b-', label='scikit-learn')
    pl.legend()
    pl.xlabel('number of samples to classify')
    pl.ylabel('time (in microseconds)')


    # now do a bench where the number of points is fixed
    # and the variable is the number of dimensions
    from scikits.learn.datasets.samples_generator import friedman, sparse_uncorrelated

    scikit_results = []
    svm_results = []
    mvpa_results = []
    n = 10
    step = 500
    start_dim = 100

    print '============================================'
    print 'Warning: this is going to take a looong time'
    print '============================================'

    dim = start_dim
    for i in range(0, n):
        print '============================================'
        print 'Entering iteration %s of %s' % (i, n)
        print '============================================'
        dim += step
        X, Y = np.random.randn(100, dim), np.random.randn(100)
        Y = Y.astype(np.int)
        bench_scikit(X, Y)
        bench_svm(X, Y)
        bench_pymvpa(X, Y)

    xx = np.arange(start_dim, start_dim+n*step, step)
    pl.subplot(212)
    pl.title('Classification in high dimensional spaces')
    pl.plot(xx, mvpa_results, 'g-', label='pymvpa')
    pl.plot(xx, svm_results,'r-', label='libsvm-swig')
    pl.plot(xx, scikit_results, 'b-', label='scikit-learn')
    pl.legend()
    pl.xlabel('number of dimensions')
    pl.ylabel('time (in seconds)')
    pl.axis('tight')
    pl.show()
