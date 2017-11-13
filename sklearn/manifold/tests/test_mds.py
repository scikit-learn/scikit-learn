import numpy as np
from numpy.testing import assert_array_almost_equal

from sklearn.manifold import mds
from sklearn.utils.testing import assert_raises


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


def test_fit_error():
    # fit() should throw an error if the dissimilarity metric is unknown
    sim = np.array([[0, 5, 3, 4],
                    [5, 0, 2, 2],
                    [3, 2, 0, 1],
                    [4, 2, 1, 0]])

    mds_clf = mds.MDS(dissimilarity="unknown")
    assert_raises(ValueError, mds_clf.fit, sim)

    mds_clf = mds.MDS(dissimilarity="unknown", method="inductive")
    assert_raises(ValueError, mds_clf.fit, sim)


def test_fit_transform_error():
    sim = np.array([[0, 5, 3, 4],
                    [5, 0, 2, 2],
                    [3, 2, 0, 1],
                    [4, 2, 1, 0]])

    Z = np.array([[-.266, -.539],
                  [.451, .252],
                  [.016, -.238],
                  [-.200, .524]])

    # fit_transform should throw an error if init is not None
    # and dissimilarity="euclidean"
    mds_clf = mds.MDS(dissimilarity="euclidean", method="inductive")
    assert_raises(ValueError, mds_clf.fit_transform, sim, init=Z)


def test_MDS():
    sim = np.array([[0, 5, 3, 4],
                    [5, 0, 2, 2],
                    [3, 2, 0, 1],
                    [4, 2, 1, 0]])
    mds_clf = mds.MDS(metric=False, n_jobs=3, dissimilarity="precomputed")
    mds_clf.fit(sim)


def test_MDS_fit_transform_inductive():
    sim = np.array([[0, 5, 3, 4],
                    [5, 0, 2, 2],
                    [3, 2, 0, 1],
                    [4, 2, 1, 0]])

    # test fit_transform under the inductive case
    mds_clf = mds.MDS(dissimilarity="euclidean", method="inductive")
    result1 = mds_clf.fit_transform(sim)

    # test fit and transform for precomputed case
    mds_clf = mds.MDS(dissimilarity="euclidean", method="inductive")
    mds_clf.fit(sim)
    result2 = mds_clf.transform(sim)

    # test that fit_transform is the same as fit + transform
    assert_array_almost_equal(result1, result2, decimal=3)


def test_MDS_transform():
    sim = np.array([[0, 5, 3, 4],
                    [5, 0, 2, 2],
                    [3, 2, 0, 1],
                    [4, 2, 1, 0]])

    mds_clf = mds.MDS(dissimilarity="euclidean", method="inductive")
    mds_clf.fit(sim)

    # Testing for extending MDS to new points
    sim2 = np.array([[3, 1, 1, 2],
                     [4, 1, 2, 2]])
    result = mds_clf.transform(sim2)
    expected = np.array([[-.705, -.452],
                         [-1.051, -1.329]])
    assert_array_almost_equal(result, expected, decimal=3)


def test_MDS_transform_error():
    sim = np.array([[0, 5, 3, 4],
                    [5, 0, 2, 2],
                    [3, 2, 0, 1],
                    [4, 2, 1, 0]])

    # calling transform with inductive=False causes an error
    mds_clf = mds.MDS(dissimilarity="euclidean")
    assert_raises(ValueError, mds_clf.transform, sim)

    # it is necessary to fit the projection before transforming new points
    mds_clf = mds.MDS(dissimilarity="euclidean", method="inductive")
    assert_raises(ValueError, mds_clf.transform, sim)

    # transform should throw an error if the dissimilarity metric is unknown
    mds_clf = mds.MDS(dissimilarity="unknown", method="inductive")
    assert_raises(ValueError, mds_clf.transform, sim)
