"""
Benchmarks of Non-Negative Matrix Factorization
"""

from __future__ import print_function

from collections import defaultdict
import gc
from time import time

import numpy as np
from scipy.linalg import norm

from sklearn.decomposition.nmf import NMF, _initialize_nmf
from sklearn.datasets.samples_generator import make_low_rank_matrix
from sklearn.externals.six.moves import xrange


def alt_nnmf(V, r, max_iter=1000, tol=1e-3, init='random'):
    '''
    A, S = nnmf(X, r, tol=1e-3, R=None)

    Implement Lee & Seung's algorithm

    Parameters
    ----------
    V : 2-ndarray, [n_samples, n_features]
        input matrix
    r : integer
        number of latent features
    max_iter : integer, optional
        maximum number of iterations (default: 1000)
    tol : double
        tolerance threshold for early exit (when the update factor is within
        tol of 1., the function exits)
    init : string
        Method used to initialize the procedure.

    Returns
    -------
    A : 2-ndarray, [n_samples, r]
        Component part of the factorization

    S : 2-ndarray, [r, n_features]
        Data part of the factorization
    Reference
    ---------
    "Algorithms for Non-negative Matrix Factorization"
    by Daniel D Lee, Sebastian H Seung
    (available at http://citeseer.ist.psu.edu/lee01algorithms.html)
    '''
    # Nomenclature in the function follows Lee & Seung
    eps = 1e-5
    n, m = V.shape
    W, H = _initialize_nmf(V, r, init, random_state=0)

    for i in xrange(max_iter):
        updateH = np.dot(W.T, V) / (np.dot(np.dot(W.T, W), H) + eps)
        H *= updateH
        updateW = np.dot(V, H.T) / (np.dot(W, np.dot(H, H.T)) + eps)
        W *= updateW
        if i % 10 == 0:
            max_update = max(updateW.max(), updateH.max())
            if abs(1. - max_update) < tol:
                break
    return W, H


def report(error, time):
    print("Frobenius loss: %.5f" % error)
    print("Took: %.2fs" % time)
    print()


def benchmark(samples_range, features_range, rank=50, tolerance=1e-5):
    timeset = defaultdict(lambda: [])
    err = defaultdict(lambda: [])

    for n_samples in samples_range:
        for n_features in features_range:
            print("%2d samples, %2d features" % (n_samples, n_features))
            print('=======================')
            X = np.abs(make_low_rank_matrix(n_samples, n_features,
                       effective_rank=rank, tail_strength=0.2))

            gc.collect()
            print("benchmarking nndsvd-nmf: ")
            tstart = time()
            m = NMF(n_components=30, tol=tolerance, init='nndsvd').fit(X)
            tend = time() - tstart
            timeset['nndsvd-nmf'].append(tend)
            err['nndsvd-nmf'].append(m.reconstruction_err_)
            report(m.reconstruction_err_, tend)

            gc.collect()
            print("benchmarking nndsvda-nmf: ")
            tstart = time()
            m = NMF(n_components=30, init='nndsvda',
                    tol=tolerance).fit(X)
            tend = time() - tstart
            timeset['nndsvda-nmf'].append(tend)
            err['nndsvda-nmf'].append(m.reconstruction_err_)
            report(m.reconstruction_err_, tend)

            gc.collect()
            print("benchmarking nndsvdar-nmf: ")
            tstart = time()
            m = NMF(n_components=30, init='nndsvdar',
                    tol=tolerance).fit(X)
            tend = time() - tstart
            timeset['nndsvdar-nmf'].append(tend)
            err['nndsvdar-nmf'].append(m.reconstruction_err_)
            report(m.reconstruction_err_, tend)

            gc.collect()
            print("benchmarking random-nmf")
            tstart = time()
            m = NMF(n_components=30, init='random', max_iter=1000,
                    tol=tolerance).fit(X)
            tend = time() - tstart
            timeset['random-nmf'].append(tend)
            err['random-nmf'].append(m.reconstruction_err_)
            report(m.reconstruction_err_, tend)

            gc.collect()
            print("benchmarking alt-random-nmf")
            tstart = time()
            W, H = alt_nnmf(X, r=30, init='random', tol=tolerance)
            tend = time() - tstart
            timeset['alt-random-nmf'].append(tend)
            err['alt-random-nmf'].append(np.linalg.norm(X - np.dot(W, H)))
            report(norm(X - np.dot(W, H)), tend)

    return timeset, err


if __name__ == '__main__':
    from mpl_toolkits.mplot3d import axes3d  # register the 3d projection
    axes3d
    import matplotlib.pyplot as plt

    samples_range = np.linspace(50, 500, 3).astype(np.int)
    features_range = np.linspace(50, 500, 3).astype(np.int)
    timeset, err = benchmark(samples_range, features_range)

    for i, results in enumerate((timeset, err)):
        fig = plt.figure('scikit-learn Non-Negative Matrix Factorization'
                         'benchmark results')
        ax = fig.gca(projection='3d')
        for c, (label, timings) in zip('rbgcm', sorted(results.iteritems())):
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
        zlabel = 'Time (s)' if i == 0 else 'reconstruction error'
        ax.set_zlabel(zlabel)
        ax.legend()
        plt.show()
