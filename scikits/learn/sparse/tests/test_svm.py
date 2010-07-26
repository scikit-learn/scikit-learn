
import numpy as np
from numpy.testing import *
from scikits.learn import sparse, svm

# test sample 1
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
Y = [1, 1, 1, 2, 2, 2]
T = [[-1, -1], [2, 2], [3, 2]]
true_result = [1, 2, 2]

# test sample 2
X2 = [[0, 0, 0], [1, 1, 1], [2, 0, 0, ],
      [0, 0, 2], [3, 3, 3]]
Y2 = [1, 2, 2, 2, 3]
T2 = [[-1, -1, -1], [1, 1, 1], [2, 2, 2]]
true_result2 = [1, 2, 3]


def test_SVC():
    """
    Check that sparse SVC gives the same result as SVC
    """

    clf = svm.SVC(kernel='linear').fit(X, Y)
    sp_clf = sparse.svm.SVC(kernel='linear').fit(X, Y)
    assert_array_equal (clf.support_, sp_clf.support_.todense())
    assert_array_equal(clf.predict(T), sp_clf.predict(T))
    
    clf.fit(X2, Y2)
    sp_clf.fit(X2, Y2)
    assert_array_equal (clf.support_, sp_clf.support_.todense())
    assert_array_equal(clf.predict(T2), sp_clf.predict(T2))

