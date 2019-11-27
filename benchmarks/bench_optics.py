"""
To run this, you'll need to have installed.
  * scikit-learn
Does two benchmarks
First, we fix a training set, increase the number of
samples to cluster and plot number of clustered samples as a
function of time.
In the second benchmark, we increase the number of dimensions of the
training set, cluster a sample and plot the time taken as a function
of the number of dimensions.
"""
import numpy as np
import matplotlib.pyplot as plt
import gc
from datetime import datetime

# to store the results
scikit_cluster_results = []

mu_second = 0.0 + 10 ** 6  # number of microseconds in a second


def bench_scikit_OPTICS(X, Y):
    """Benchmark with scikit-learn OPTICS"""

    from sklearn.cluster import OPTICS
    gc.collect()

    # start time
    tstart = datetime.now()
    clf = OPTICS()
    clf.fit(X, Y)
    delta = (datetime.now() - tstart)
    # stop time

    scikit_cluster_results.append(
        delta.seconds + delta.microseconds / mu_second)


if __name__ == '__main__':

    print('============================================')
    print('Warning: this is going to take a looong time')
    print('============================================')

    n = 10
    step = 1000
    n_samples = 20000
    dim = 10
    n_classes = 10
    for i in range(n):
        print('============================================')
        print('Entering iteration %s of %s' % (i, n))
        print('============================================')
        n_samples += step
        X = np.random.randn(n_samples, dim)
        Y = np.random.randint(0, n_classes, (n_samples,))
        bench_scikit_OPTICS(X, Y)

    xx = range(0, n * step, step)
    plt.figure('scikit-learn OPTICS benchmark results')
    plt.subplot(211)
    plt.title('Learning with varying number of samples')
    plt.plot(xx, scikit_cluster_results, 'g-', label='OPTICS')
    plt.legend(loc='upper left')
    plt.xlabel('number of samples')
    plt.ylabel('Time (s)')

    scikit_cluster_results = []
    n = 10
    step = 500
    start_dim = 500
    n_classes = 10

    dim = start_dim
    for i in range(0, n):
        print('============================================')
        print('Entering iteration %s of %s' % (i, n))
        print('============================================')
        dim += step
        X = np.random.randn(100, dim)
        Y = np.random.randint(0, n_classes, (100,))
        bench_scikit_OPTICS(X, Y)

    xx = np.arange(start_dim, start_dim + n * step, step)
    plt.subplot(212)
    plt.title('OPTICS in high dimensional spaces')
    plt.plot(xx, scikit_cluster_results, 'g-', label='classification')
    plt.legend(loc='upper left')
    plt.xlabel('number of dimensions')
    plt.ylabel('Time (s)')
    plt.axis('tight')
    plt.show()
