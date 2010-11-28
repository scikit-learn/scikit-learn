"""Benchmarks of Singular Values Decomposition (Exact and Approximate)"""
import gc
from time import time
import numpy as np
from collections import defaultdict

from scipy.linalg import svd
from scikits.learn.utils.extmath import fast_svd


def make_data(n_samples, n_features, rank=None):
    """Build a fixed rank matrix with logspaced singular values"""
    np.random.seed(42)
    if rank is None:
        rank = min(n_samples, n_samples) / 10 + 1

    s = np.identity(rank) * np.logspace(10, 10e-5, rank)
    X = np.dot(np.random.randn(n_samples, rank), s)
    X = np.dot(X, np.random.randn(rank, n_features))
    return X


def compute_bench(data_gen, samples_range, features_range, q=3):

    it = 0

    results = defaultdict(lambda: [])

    max_it = len(samples_range) * len(features_range)
    for n_samples in samples_range:
        for n_features in features_range:
            it += 1
            print '===================='
            print 'Iteration %03d of %03d' % (it, max_it)
            print '===================='
            X = make_data(n_samples, n_features)
            rank = min(n_samples, n_samples) / 10 + 1

            gc.collect()
            print "benching scipy svd: "
            tstart = time()
            svd(X, full_matrices=False)
            results['scipy svd'].append(time() - tstart)

            gc.collect()
            print "benching scikit-learn fast_svd: q=0"
            tstart = time()
            fast_svd(X, rank, q=0)
            results['scikit-learn fast_svd (q=0)'].append(time() - tstart)

            gc.collect()
            print "benching scikit-learn fast_svd: q=%d " % q
            tstart = time()
            fast_svd(X, rank, q=q)
            results['scikit-learn fast_svd (q=%d)' % q].append(time() - tstart)

    return results


if __name__ == '__main__':
    from mpl_toolkits.mplot3d import axes3d
    import matplotlib.pyplot as plt

    samples_range = np.linspace(2, 1000, 4)
    features_range = np.linspace(2, 1000, 4)
    results = compute_bench(make_data, samples_range, features_range)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for c, (label, timings) in zip('rbg', sorted(results.iteritems())):
        X, Y = np.meshgrid(samples_range, features_range)
        Z = np.asarray(timings).reshape(samples_range.shape[0],
                                        features_range.shape[0])
        # plot the actual surface
        ax.plot_surface(X, Y, Z, rstride=8, cstride=8, alpha=0.3,
                        color=c)
        # dummy point plot to stick the legend to since surface plot do not
        # support legends (yet?)
        ax.plot([1], [1], [1], color=c, label=label)

    ax.set_xlabel('n_samples')
    ax.set_ylabel('n_features')
    ax.set_zlabel('time (s)')
    ax.legend()
    plt.show()

