"""
Benchmarks of Non-Negative Matrix Factorization
"""

import gc
from time import time
import numpy as np
from collections import defaultdict

from scikits.learn.nmf import NMF, _initialize_nmf_
from scikits.learn.datasets.samples_generator import low_rank_fat_tail


def alt_nnmf(V, r, max_iter=1000, tol=1e-3, R=None):
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
        maximum number of iterations (default: 10000)
    tol : double
        tolerance threshold for early exit (when the update factor is within
        tol of 1., the function exits)
    R : integer, optional
        random seed

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
    if R == "svd":
        W, H = _initialize_nmf_(V, r)
    elif R == None:
        R = np.random.mtrand._rand
        W = np.abs(R.standard_normal((n, r)))
        H = np.abs(R.standard_normal((r, m)))

    for i in xrange(max_iter):
        updateH = np.dot(W.T, V) / (np.dot(np.dot(W.T, W), H) + eps)
        H *= updateH
        updateW = np.dot(V, H.T) / (np.dot(W, np.dot(H, H.T)) + eps)
        W *= updateW
        if True or (i % 10) == 0:
            max_update = max(updateW.max(), updateH.max())
            if abs(1. - max_update) < tol:
                break
    return W, H


def compute_bench(samples_range, features_range, rank=50, tolerance=1e-7):
    it = 0
    timeset = defaultdict(lambda: [])
    err = defaultdict(lambda: [])

    max_it = len(samples_range) * len(features_range)
    for n_samples in samples_range:
        for n_features in features_range:
            it += 1
            print '===================='
            print 'Iteration %03d of %03d' % (it, max_it)
            print '===================='
            X = np.abs(low_rank_fat_tail(n_samples, n_features,
                       effective_rank=rank,  tail_strength=0.2))

            gc.collect()
            print "benching nndsvd-nmf: "
            tstart = time()
            m = NMF(n_components=30, tol=tolerance, init='nndsvd').fit(X)
            tend = time() - tstart
            timeset['nndsvd-nmf'].append(tend)
            err['nndsvd-nmf'].append(m.reconstruction_err_)
            print m.reconstruction_err_, tend

            gc.collect()
            print "benching nndsvda-nmf: "
            tstart = time()
            m = NMF(n_components=30, init='nndsvda',
                    tol=tolerance).fit(X)
            tend = time() - tstart
            timeset['nndsvda-nmf'].append(tend)
            err['nndsvda-nmf'].append(m.reconstruction_err_)
            print m.reconstruction_err_, tend

            gc.collect()
            print "benching nndsvdar-nmf: "
            tstart = time()
            m = NMF(n_components=30, init='nndsvdar',
                    tol=tolerance).fit(X)
            tend = time() - tstart
            timeset['nndsvdar-nmf'].append(tend)
            err['nndsvdar-nmf'].append(m.reconstruction_err_)
            print m.reconstruction_err_, tend

            gc.collect()
            print "benching random-nmf"
            tstart = time()
            m = NMF(n_components=30, init=None, max_iter=1000,
                    tol=tolerance).fit(X)
            tend = time() - tstart
            timeset['random-nmf'].append(tend)
            err['random-nmf'].append(m.reconstruction_err_)
            print m.reconstruction_err_, tend

            gc.collect()
            print "benching alt-random-nmf"
            tstart = time()
            W, H = alt_nnmf(X, r=30, R=None, tol=tolerance)
            tend = time() - tstart
            timeset['alt-random-nmf'].append(tend)
            err['alt-random-nmf'].append(np.linalg.norm(X - np.dot(W, H)))
            print np.linalg.norm(X - np.dot(W, H)), tend

    return timeset, err


if __name__ == '__main__':
    from mpl_toolkits.mplot3d import axes3d  # register the 3d projection
    import matplotlib.pyplot as plt

    samples_range = np.linspace(50, 500, 3).astype(np.int)
    features_range = np.linspace(50, 500, 3).astype(np.int)
    timeset, err = compute_bench(samples_range, features_range)

    for i, results in enumerate((timeset, err)):
        fig = plt.figure()
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
        zlabel = 'time (s)' if i == 0 else 'reconstruction error'
        ax.set_zlabel(zlabel)
        ax.legend()
        plt.show()
