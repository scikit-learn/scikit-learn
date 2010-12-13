"""Benchmarks of Singular Values Decomposition (Exact and Approximate)

The data is mostly low rank but is a fat infinite tail.
"""
import gc
from time import time
import numpy as np
from collections import defaultdict

from scikits.learn.nmf import NMF, _initialize_nmf_
from scikits.learn.datasets.samples_generator import low_rank_fat_tail


def alt_nnmf(V, r, cost='norm2', max_iter=1000, tol=1e-3, R=None):
    '''
    A, S = nnmf(X, r, cost='norm2', tol=1e-5, R=None)

    Implement Lee & Seung's algorithm

    Parameters
    ----------
    V : 2-ndarray
        input matrix
    r : integer
        nr of latent features
    cost : one of:
        'norm2' : minimise || X - AS ||_2 (default)
        'i-div' : minimise D(X||AS), where D is I-divergence (generalisation of K-L divergence)
    max_iter : integer, optional
        maximum number of iterations (default: 10000)
    tol : double
        tolerance threshold for early exit (when the update factor is with tol
        of 1., the function exits)
    R : integer, optional
        random seed

    Returns
    -------
    A : 2-ndarray
    S : 2-ndarray

    Reference
    ---------
    "Algorithms for Non-negative Matrix Factorization"
    by Daniel D Lee, Sebastian H Seung
    (available at http://citeseer.ist.psu.edu/lee01algorithms.html)
    '''
    # Nomenclature in the function follows lee & seung, while outside nomenclature follows
    eps = 1e-5
    n,m = V.shape
    if R == "svd":
        W, H = _initialize_nmf_(V, r)
    elif R == None:
        R = np.random.mtrand._rand
        W = np.abs(R.standard_normal((n,r)))
        H = np.abs(R.standard_normal((r,m)))

    for i in xrange(max_iter):
        if cost == 'norm2':
            updateH = np.dot(W.T, V) / (np.dot(np.dot(W.T, W), H) + eps)
            H *= updateH
            updateW = np.dot(V, H.T) / (np.dot(W, np.dot(H, H.T)) + eps)
            W *= updateW
        elif cost == 'i-div':
            raise NotImplementedError,'I-Div not implemented in lee_seung.nnmf'
        if True or (i % 10) == 0:
            max_update = max(updateW.max(),updateH.max())
            if abs(1.-max_update) < tol:
                break
    return W, H


def compute_bench(samples_range, features_range, rank=50, tolerance=1e-3):

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
            X = np.abs(low_rank_fat_tail(n_samples, n_features, effective_rank=rank,
                                  tail_strength=0.2))
            
            gc.collect()
            print "benching svd-nmf: "
            tstart = time()
            m = NMF(50, max_iter=1000, tolerance=tolerance).fit(X)
            timeset['svd-nmf'].append(time() - tstart)
            err['svd-nmf'].append(m.reconstruction_err_)

        #    gc.collect()
        #    print "benching alt-svd-nmf: "
        #    tstart = time()
        #    W, H = alt_nnmf(X, r=50, R="svd")
        #    timeset['alt-svd-nmf'].append(time() - tstart)
        #    err['alt-svd-nmf'].append(np.linalg.norm(X - np.dot(W,H)))

            gc.collect()
            print "benching random-nmf"
            tstart = time()
            m = NMF(50, initial=None, max_iter=1000, tolerance=tolerance).fit(X)
            timeset['random-nmf'].append(time() - tstart)
            err['random-nmf'].append(m.reconstruction_err_)

            gc.collect()
            print "benching alt-random-nmf"
            tstart = time()
            W, H = alt_nnmf(X, r=50, R=None, tol=tolerance)
            timeset['alt-random-nmf'].append(time() - tstart)
            err['alt-random-nmf'].append(np.linalg.norm(X - np.dot(W,H)))

    return timeset, err


if __name__ == '__main__':
    from mpl_toolkits.mplot3d import axes3d # register the 3d projection
    import matplotlib.pyplot as plt

    samples_range = np.linspace(50, 1000, 4).astype(np.int)
    features_range = np.linspace(50, 1000, 4).astype(np.int)
    timeset, err = compute_bench(samples_range, features_range)

    for results in (timeset,err):
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
