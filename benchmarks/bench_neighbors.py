"""
This script compares the performance of the various nearest neighbors
algorithms available in NearestNeighbors: ball_tree, kd_tree, and brute

Then run the simple timings script:
 python bench_neighbors.py 1000 100
"""

import sys
from time import time
import numpy as np
from sklearn.neighbors import NearestNeighbors


def compare_nbrs(nbrs1, nbrs2):
    assert nbrs1.shape == nbrs2.shape
    if(nbrs1.ndim == 2):
        N, k = nbrs1.shape
        for i in range(N):
            for j in range(k):
                if nbrs1[i, j] == i:
                    continue
                elif nbrs1[i, j] not in nbrs2[i]:
                    return False
        return True
    elif(nbrs1.ndim == 1):
        N = len(nbrs1)
        return np.all(nbrs1 == nbrs2)


def test_time(n_samples=1000, n_features=100, leaf_size=20, k=20):
    X = np.random.random([n_samples, n_features])

    print "---------------------------------------------------"
    print "%i neighbors of %i points in %i dimensions:" \
        % (k, n_samples, n_features)
    print "   (leaf size = %i)" % leaf_size
    print "  -------------"
    BT = NearestNeighbors(algorithm='ball_tree',
                          leaf_size=leaf_size)
    KDT = NearestNeighbors(algorithm='kd_tree',
                           leaf_size=leaf_size)
    Brute = NearestNeighbors(algorithm='brute')

    t0 = time()
    BT.fit(X)
    print "  Ball Tree construction     : %.3g sec" % (time() - t0)
    d, nbrs1 = BT.kneighbors(X, k)
    print "  total (construction+query) : %.3g sec" % (time() - t0)
    print "  -------------"

    t0 = time()
    KDT.fit(X)
    print "  KD tree construction       : %.3g sec" % (time() - t0)
    d, nbrs2 = KDT.kneighbors(X, k)
    print "  total (construction+query) : %.3g sec" % (time() - t0)
    print "  -------------"

    t0 = time()
    Brute.fit(X)
    print "  Brute Force construction   : %.3g sec" % (time() - t0)
    d, nbrs3 = Brute.kneighbors(X, k)
    print "  total (construction+query) : %.3g sec" % (time() - t0)
    print "  -------------"

    print "  neighbors match: ",
    print compare_nbrs(nbrs1, nbrs2) and compare_nbrs(nbrs2, nbrs3)
    print "  -------------"

if __name__ == '__main__':
    if len(sys.argv) == 3:
        n_samples, n_features = map(int, sys.argv[1:])
        leaf_size = 20
        k = min(20, n_samples)

    elif len(sys.argv) == 4:
        n_samples, n_features, leaf_size = map(int, sys.argv[1:])
        k = min(20, n_samples)

    elif len(sys.argv) == 5:
        n_samples, n_features, leaf_size, k = map(int, sys.argv[1:])

    else:
        print "usage: bench_neighbors.py n_samples n_features " + \
              "[leafsize=20], [k=20]"
        sys.exit()

    test_time(n_samples, n_features, leaf_size, k)
