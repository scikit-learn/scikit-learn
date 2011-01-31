"""
This script compares the performance of the Ball Tree code
with the cKDTree from scipy.spatial

"""

from scikits.learn.ball_tree import BallTree
import numpy as np
from time import time

from scipy.spatial import cKDTree
import pylab as pl


def compare_nbrs(nbrs1, nbrs2):
    assert nbrs1.shape == nbrs2.shape
    if(nbrs1.ndim == 2):
        n_samples, k = nbrs1.shape
        for i in range(n_samples):
            for j in range(k):
                if nbrs1[i, j] == i:
                    continue
                elif nbrs1[i, j] not in nbrs2[i]:
                    return False
        return True
    elif(nbrs1.ndim == 1):
        return np.all(nbrs1 == nbrs2)

if __name__ == '__main__':
    n_samples = 1000
    leaf_size = 1 # leaf size
    k = 20
    BT_results = []
    KDT_results = []

    for i in range(1, 10):
        print 'Iteration %s' %i
        n_features = i*100
        X = np.random.random([n_samples, n_features])

        t0 = time()
        BT = BallTree(X, leaf_size)
        d, nbrs1 = BT.query(X, k)
        delta = time() - t0
        BT_results.append(delta)

        t0 = time()
        KDT = cKDTree(X, leaf_size)
        d, nbrs2 = KDT.query(X, k)
        delta = time() - t0
        KDT_results.append(delta)

        # this checks we get the correct result
        assert compare_nbrs(nbrs1, nbrs2)

    xx = 100 * np.arange(1, 10)
    pl.plot(xx, BT_results, label='scikits.learn (BallTree)')
    pl.plot(xx, KDT_results, label='scipy (cKDTree)')
    pl.xlabel('number of dimensions')
    pl.ylabel('time (seconds)')
    pl.legend()
    pl.show()
