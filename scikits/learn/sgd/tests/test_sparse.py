import numpy as np
from scipy import sparse
from scikits.learn import datasets, sgd
from numpy.testing import assert_array_almost_equal, \
     assert_array_equal, assert_equal

from nose.tools import assert_raises

import bolt

# test sample 1
X = np.array([[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]])
Y = [1, 1, 1, 2, 2, 2]
T = np.array([[-1, -1], [2, 2], [3, 2]])
true_result = [1, 2, 2]

# test sample 2
X2 = np.array([[0, 0, 0], [1, 1, 1], [2, 0, 0, ],
               [0, 0, 2], [3, 3, 3]])
Y2 = [1, 2, 2, 2, 3]
T2 = np.array([[-1, -1, -1], [1, 1, 1], [2, 2, 2]])
true_result2 = [1, 2, 3]


def test_sgd():
    """Check that sparse SGD gives any results :-)"""
    
    clf = sgd.sparse.SGD(penalty='l2', alpha = 0.001, fit_intercept = True)
    clf.fit(X, Y)
    assert_array_equal(clf.coef_, np.zeros((len(clf.coef_),)))

def test_rcv1():
    ds = bolt.io.MemoryDataset.load("/home/peter/corpora/rcv1-ccat/test.npy")
    m, n = ds.n, ds.dim
    X = sparse.lil_matrix((m, n), dtype = np.float32)
    print "Build sparse matrix... ", 
    for i, x in enumerate(ds.iterinstances()):
	X[i, x['f0']] = x['f1']
    print "[done]"

    X = X.tocsr()
    Y = ds.labels
    print "Fitting model... "
    clf = sgd.sparse.SGD(penalty='l2', alpha = 0.0001, fit_intercept = True)
    clf.fit(X, Y)
    score = clf.score(X,Y)
    print "training score: ", score
    print "training error: ", ((1.0 - score) * 100.0)

if __name__ == "__main__":
    test_rcv1()
    
