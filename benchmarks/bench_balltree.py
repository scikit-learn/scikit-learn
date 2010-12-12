"""
This script compares the performance of the Ball Tree code with
scipy.spatial.cKDTree.

Then run the simple timings script:
 python bench_kdtree.py 1000 100
"""

from scikits.learn.ball_tree import BallTree
import numpy
from time import time

from scipy.spatial import cKDTree
import sys


def compare_nbrs(nbrs1, nbrs2):
    assert nbrs1.shape == nbrs2.shape
    if(nbrs1.ndim == 2):
        N, k = nbrs1.shape
        for i in range(N):
            for j in range(k):
                if nbrs1[i, j]==i:
                    continue
                elif nbrs1[i, j] not in nbrs2[i]:
                    return False
        return True
    elif(nbrs1.ndim == 1):
        N = len(nbrs1)
        return numpy.all(nbrs1 == nbrs2)


def test_time(n_samples=1000, n_features=100, leaf_size=1, k=20):
    X = numpy.random.random([n_samples, n_features])

    print "---------------------------------------------------"
    print "%i neighbors of %i points in %i dimensions:" % (k, n_samples, n_features)
    print "   (leaf size = %i)" % leaf_size
    print "  -------------"

    t0 = time()
    BT = BallTree(X, leaf_size)
    print "  Ball Tree construction     : %.3g sec" % (time() - t0)
    d, nbrs1 = BT.query(X, k)
    print "  total (construction+query) : %.3g sec" % (time() - t0)
    print "  -------------"


    t0 = time()
    KDT = cKDTree(X, leaf_size)
    print "  KD tree construction       : %.3g sec" % (time() - t0)
    d, nbrs2 = KDT.query(X, k)
    print "  total (construction+query) : %.3g sec" % (time() - t0)
    print "  -------------"

    print "  neighbors match: ",
    print compare_nbrs(nbrs1, nbrs2)
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
        print "usage: bench_balltree.py n_samples n_features " + \
              "[leafsize=20], [k=20]"
        exit()

    test_time(n_samples, n_features, leaf_size, k)
