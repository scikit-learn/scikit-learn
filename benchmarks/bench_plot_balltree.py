
from scikits.learn.ball_tree import BallTree, knn_brute
import numpy as np
from time import time

from scipy.spatial import cKDTree
import sys
import pylab as pl

def compare_nbrs(nbrs1,nbrs2):
    assert nbrs1.shape == nbrs2.shape
    if(nbrs1.ndim == 2):
        N,k = nbrs1.shape
        for i in range(N):
            for j in range(k):
                if nbrs1[i,j]==i:
                    continue
                elif nbrs1[i,j] not in nbrs2[i]:
                    return False
        return True
    elif(nbrs1.ndim == 1):
        N = len(nbrs1)
        return numpy.all(nbrs1 == nbrs2)

N = 1000
ls = 1 # leaf size
k = 20
BT_results = []
KDT_results = []

for i in range(1, 10):
    print 'Iteration %s' %i
    D = i*100
    M = np.random.random([N, D])

    t0 = time()
    BT = BallTree(M, ls)
    d, nbrs1 = BT.query(M, k)
    delta = time() - t0
    BT_results.append(delta)

    t0 = time()
    KDT = cKDTree(M, ls)
    d, nbrs2 = KDT.query(M, k)
    delta = time() - t0
    KDT_results.append(delta)

    # this checks we get the correct result
    assert compare_nbrs(nbrs1,nbrs2)

xx = 100*np.arange(1, 10)
pl.plot(xx, BT_results, label='scikits.learn (BallTree)')
pl.plot(xx, KDT_results, label='scipy (cKDTree)')
pl.xlabel('number of dimensions')
pl.ylabel('time (seconds)')
pl.legend()
pl.show()

