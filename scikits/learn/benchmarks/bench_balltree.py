"""
This script compares the performance of the Ball Tree code with
scipy.spatial.cKDTree.

Then run the simple timings script:
 python bench_kdtree.py 1000 100
"""

from scikits.learn.ball_tree import BallTree, knn_brute
import numpy
from time import time

from scipy.spatial import cKDTree
import sys


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


def test_time(N=1000, D=100, ls=1, k=20):
    M = numpy.random.random([N,D])

    print "---------------------------------------------------"
    print "%i neighbors of %i points in %i dimensions:" % (k,N,D)
    print "   (leaf size = %i)" % ls
    print "  -------------"

    t0 = time()
    BT = BallTree(M,ls)
    print "  Ball Tree construction     : %.3g sec" % ( time()-t0 )
    d,nbrs1 = BT.query(M,k)
    print "  total (construction+query) : %.3g sec" % ( time()-t0 )
    print "  -------------"


    t0 = time()
    KDT = cKDTree(M,ls)
    print "  KD tree construction       : %.3g sec" % ( time()-t0 )
    d,nbrs2 = KDT.query(M,k)
    print "  total (construction+query) : %.3g sec" % ( time()-t0 )
    print "  -------------"

    print "  neighbors match: ",
    print ( compare_nbrs(nbrs1,nbrs2) )
    print "  -------------"

if __name__ == '__main__':
    if len(sys.argv)==3:
        N,D = map(int,sys.argv[1:])
        ls = 20
        k = min(20,N)

    elif len(sys.argv)==4:
        N,D,ls = map(int,sys.argv[1:])
        k = min(20,N)

    elif len(sys.argv)==5:
        N,D,ls,k = map(int,sys.argv[1:])

    else:
        print "usage: bench_balltree.py N D [leafsize=20], [k=20]"
        exit()


    test_time(N,D,ls,k)



