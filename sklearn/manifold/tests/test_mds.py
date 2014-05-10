import numpy as np
from numpy.testing import assert_array_almost_equal

from nose.tools import assert_raises
from sklearn.manifold import mds
from sklearn.metrics import euclidean_distances
from sklearn.utils.testing import assert_less


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

    assert_raises(ValueError, mds.smacof, sim)

    # Not squared similarity matrix:
    sim = np.array([[0, 5, 9, 4],
                    [5, 0, 2, 2],
                    [4, 2, 1, 0]])

    assert_raises(ValueError, mds.smacof, sim)

    # init not None and not correct format:
    sim = np.array([[0, 5, 3, 4],
                    [5, 0, 2, 2],
                    [3, 2, 0, 1],
                    [4, 2, 1, 0]])

    Z = np.array([[-.266, -.539],
                  [.016, -.238],
                  [-.200, .524]])
    assert_raises(ValueError, mds.smacof, sim, init=Z, n_init=1)


def test_MDS():
    sim = np.array([[0, 5, 3, 4],
                    [5, 0, 2, 2],
                    [3, 2, 0, 1],
                    [4, 2, 1, 0]])
    mds_clf = mds.MDS(metric=False, n_jobs=3, dissimilarity="precomputed")
    mds_clf.fit(sim)


def test_svd_mds():
    # Generate 4 randomly chosen points
    Y = np.array([[1,  0,  1],
                  [-1, 3,  2],
                  [1, -2,  3],
                  [2, -1, -3]])
    sim = euclidean_distances(Y)
    # calculate error or smacof-based solution
    X_smacof, smacof_stress = mds.smacof(sim, n_components=2, random_state=42)
    X_svd, svd_stress = mds.svd_mds(sim, n_components=2)
    assert_less(svd_stress, smacof_stress)


def test_svd_mds_non_square_matrix():
    sim = np.ones((3, 4))
    assert_raises(ValueError, mds.svd_mds, sim)


def test_svd_mds_non_symmetric_matrix():
    sim = np.ones((3, 3))
    sim[0, 1] = 0
    assert_raises(ValueError, mds.svd_mds, sim)


def test_svd_mds_non_euclidean():
    # non euclidean similarities (no triangular inequality)
    sim = np.array([[0, 12, 3, 4],
                    [12, 0, 2, 2],
                    [3,  2, 0, 1],
                    [4,  2, 1, 0]])
    assert_raises(ValueError, mds.svd_mds, sim)


def test_MDS_svd():
    Y = np.array([[1, 0, 1],
                  [-1, 3, 2],
                  [1, -2, 3],
                  [2, -1, -3]])
    mds_clf = mds.MDS(metric=True, method="svd")
    mds_clf.fit(Y)


def test_MDS_svd_non_metric():
    Y = np.array([[1, 0, 1],
                  [-1, 3, 2],
                  [1, -2, 3],
                  [2, -1, -3]])
    mds_clf = mds.MDS(metric=False, method="svd")
    assert_raises(ValueError, mds_clf.fit, Y)


def test_MDS_unknown_method():
    Y = np.array([[1, 0, 1],
                  [-1, 3, 2],
                  [1, -2, 3],
                  [2, -1, -3]])
    mds_clf = mds.MDS(metric=True, method="unknown_method")
    assert_raises(ValueError, mds_clf.fit, Y)
