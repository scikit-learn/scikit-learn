import numpy as np
from scikits.learn import sgd
from numpy.testing import assert_array_equal

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

# test sample 3
X3 = np.array([[1,1,0,0,0,0], [1,1,0,0,0,0],
               [0,0,1,0,0,0], [0,0,1,0,0,0],
               [0,0,0,0,1,1], [0,0,0,0,1,1],
               [0,0,0,1,0,0], [0,0,0,1,0,0]])
Y3 = np.array([1, 1, 1, 1, 2, 2, 2, 2])

X4 = np.array([[1,0.9,0.8,0,0,0], [1,.84,.98,0,0,0],
               [1,.96,.88,0,0,0], [1,.91,.99,0,0,0],
               [0,0,0,.89,.91,1], [0,0,0,.79,.84,1],
               [0,0,0,.91,.95,1], [0,0,0,.93,1,1]])
Y4 = np.array([1, 1, 1, 1, 2, 2, 2, 2])

X5 = np.array([[1,1,1,0,0,0], [1,1,1,0,0,0],
               [1,1,1,0,0,0], [1,1,1,0,0,0],
               [0,0,0,1,1,1], [0,0,0,1,1,1],
               [0,0,0,1,1,1], [0,0,0,1,1,1]])
Y5 = np.array([1, 1, 1, 1, 2, 2, 2, 2])


def test_sgd():
    """Check that sparse SGD gives any results :-)"""

    clf = sgd.sparse.SGD(penalty='l2', alpha=0.01,
                         fit_intercept=True,
                         n_iter=10, shuffle=True)
    clf.fit(X, Y)
    #assert_almost_equal(clf.coef_[0], clf.coef_[1], decimal=7)
    assert_array_equal(clf.predict(T), true_result)


def test_sgd_penalties():
    """Check whether penalties and hyperparameters are set properly"""
    clf = sgd.sparse.SGD(penalty='l2')
    assert clf.rho == 1.0
    clf = sgd.sparse.SGD(penalty='l1')
    assert clf.rho == 0.0
    clf = sgd.sparse.SGD(penalty='elasticnet', rho=0.85)
    assert clf.rho == 0.85


def test_sgd_params():
    """Test parameter validity check"""
    try:
        clf = sgd.sparse.SGD(n_iter=0)
        clf = sgd.sparse.SGD(n_iter=-10000)
    except ValueError:
        pass
    else:
        assert False

    try:
        clf = sgd.sparse.SGD(shuffle="false")
    except ValueError:
        pass
    else:
        assert False


def test_sgd_multiclass():
    """SGD is not able to handle multi class problems"""
    clf = sgd.sparse.SGD()
    try:
        clf.fit(X2, Y2)
    except ValueError:
        pass
    else:
        assert False


def test_sgd_l1():
    n = len(X4)
    np.random.seed(13)
    idx = np.arange(n)
    np.random.shuffle(idx)
    X = X4[idx, :]
    Y = Y4[idx, :]
    clf = sgd.sparse.SGD(penalty='l1', alpha=.2, fit_intercept=False,
                         n_iter=1000)
    clf.fit(X, Y)
    assert_array_equal(clf.coef_[1:-1], np.zeros((4,)))

