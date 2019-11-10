import numpy as np
from numpy.testing import assert_array_almost_equal
import pytest

from sklearn.manifold import _mds as mds


def test_smacof():
    # test metric smacof using the data of "Modern Multidimensional Scaling",
    # Borg & Groenen, p 154
    sim = np.array([[0, 5, 3, 4],
                    [5, 0, 2, 2],
                    [3, 2, 0, 1],
                    [4, 2, 1, 0]])
    Z = np.array([[-.266, -.539],
                  [.451, .252],
                  [.016, -.238],
                  [-.200, .524]])
    X, _ = mds.smacof(sim, init=Z, n_components=2, max_iter=1, n_init=1)
    X_true = np.array([[-1.415, -2.471],
                       [1.633, 1.107],
                       [.249, -.067],
                       [-.468, 1.431]])
    assert_array_almost_equal(X, X_true, decimal=3)


def test_smacof_error():
    # Not symmetric similarity matrix:
    sim = np.array([[0, 5, 9, 4],
                    [5, 0, 2, 2],
                    [3, 2, 0, 1],
                    [4, 2, 1, 0]])

    with pytest.raises(ValueError):
        mds.smacof(sim)

    # Not squared similarity matrix:
    sim = np.array([[0, 5, 9, 4],
                    [5, 0, 2, 2],
                    [4, 2, 1, 0]])

    with pytest.raises(ValueError):
        mds.smacof(sim)

    # init not None and not correct format:
    sim = np.array([[0, 5, 3, 4],
                    [5, 0, 2, 2],
                    [3, 2, 0, 1],
                    [4, 2, 1, 0]])

    Z = np.array([[-.266, -.539],
                  [.016, -.238],
                  [-.200, .524]])
    with pytest.raises(ValueError):
        mds.smacof(sim, init=Z, n_init=1)


def test_svd():
    # Test svd using example data from "An Introduction to MDS"
    # Florian Wickelmaier, p 11
    sim = np.array([[0, 93, 82, 133],
                    [93, 0, 52, 60],
                    [82, 52, 0, 111],
                    [133, 60, 111, 0]])

    X, stress = mds.svd(sim, n_components=2)
    X_true_1 = np.array([[-62.831, -32.97448],
                         [18.403, 12.02697],
                         [-24.960, 39.71091],
                         [69.388, -18.76340]])
    X_true_2 = np.copy(X_true_1)
    X_true_2[:, 0] *= -1

    # Signs of columns are dependent on signs of computed eigenvectors
    # which are arbitrary and meaningless
    assert(np.allclose(X, X_true_1)
           or np.allclose(X, -X_true_1)
           or np.allclose(X, X_true_2)
           or np.allclose(X, -X_true_2))


def test_svd_error():
    # Non symmetric (dis)similarity matrix:
    sim = np.array([[0, 5, 9, 4],
                    [5, 0, 2, 2],
                    [3, 2, 0, 1],
                    [4, 2, 1, 0]])

    with pytest.raises(ValueError):
        mds.svd(sim)

    # Non squared (dis)similarity matrix:
    sim = np.array([[0, 5, 9, 4],
                    [5, 0, 2, 2],
                    [4, 2, 1, 0]])

    with pytest.raises(ValueError):
        mds.svd(sim)

    # Non Euclidean (dis)similarity matrix:
    sim = np.array([[0, 12, 3, 4],
                    [12, 0, 2, 2],
                    [3, 2, 0, 1],
                    [4, 2, 1, 0]])

    with pytest.raises(ValueError):
        mds.svd(sim)


def test_MDS_error():
    # Bad method name
    sim = np.ones((2, 2))
    mdc_clf = mds.MDS(method='bad name')
    with pytest.raises(ValueError):
        mdc_clf.fit(sim)

    # SVD with metric=False
    sim = np.ones((2, 2))
    mdc_clf = mds.MDS(metric=False, method='svd')
    with pytest.raises(ValueError):
        mdc_clf.fit(sim)


def test_MDS():
    sim = np.array([[0, 5, 3, 4],
                    [5, 0, 2, 2],
                    [3, 2, 0, 1],
                    [4, 2, 1, 0]])
    mds_clf = mds.MDS(metric=False, n_jobs=3, dissimilarity="precomputed")
    mds_clf.fit(sim)


def test_MDS_svd():
    sim = np.array([[0, 5, 3, 4],
                    [5, 0, 2, 2],
                    [3, 2, 0, 1],
                    [4, 2, 1, 0]])
    mds_clf = mds.MDS(metric=True, method="svd")
    mds_clf.fit(sim)
