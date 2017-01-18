from __future__ import division

from collections import defaultdict
from functools import partial

import numpy as np
import scipy.sparse as sp
from sklearn.externals.six.moves import zip

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_raises

from sklearn.datasets import make_classification
from sklearn.datasets import make_multilabel_classification
from sklearn.datasets import make_hastie_10_2
from sklearn.datasets import make_regression
from sklearn.datasets import make_blobs
from sklearn.datasets import make_friedman1
from sklearn.datasets import make_friedman2
from sklearn.datasets import make_friedman3
from sklearn.datasets import make_low_rank_matrix
from sklearn.datasets import make_moons
from sklearn.datasets import make_sparse_coded_signal
from sklearn.datasets import make_sparse_uncorrelated
from sklearn.datasets import make_spd_matrix
from sklearn.datasets import make_swiss_roll
from sklearn.datasets import make_s_curve
from sklearn.datasets import make_biclusters
from sklearn.datasets import make_checkerboard

from sklearn.utils.validation import assert_all_finite


def test_make_classification():
    X, y = make_classification(n_samples=100, n_features=20, n_informative=5,
                               n_redundant=1, n_repeated=1, n_classes=3,
                               n_clusters_per_class=1, hypercube=False,
                               shift=None, scale=None, weights=[0.1, 0.25],
                               random_state=0)

    assert_equal(X.shape, (100, 20), "X shape mismatch")
    assert_equal(y.shape, (100,), "y shape mismatch")
    assert_equal(np.unique(y).shape, (3,), "Unexpected number of classes")
    assert_equal(sum(y == 0), 10, "Unexpected number of samples in class #0")
    assert_equal(sum(y == 1), 25, "Unexpected number of samples in class #1")
    assert_equal(sum(y == 2), 65, "Unexpected number of samples in class #2")


def test_make_classification_informative_features():
    """Test the construction of informative features in make_classification

    Also tests `n_clusters_per_class`, `n_classes`, `hypercube` and
    fully-specified `weights`.
    """
    # Create very separate clusters; check that vertices are unique and
    # correspond to classes
    class_sep = 1e6
    make = partial(make_classification, class_sep=class_sep, n_redundant=0,
                   n_repeated=0, flip_y=0, shift=0, scale=1, shuffle=False)

    for n_informative, weights, n_clusters_per_class in [(2, [1], 1),
                                                         (2, [1/3] * 3, 1),
                                                         (2, [1/4] * 4, 1),
                                                         (2, [1/2] * 2, 2),
                                                         (2, [3/4, 1/4], 2),
                                                         (10, [1/3] * 3, 10)
                                                         ]:
        n_classes = len(weights)
        n_clusters = n_classes * n_clusters_per_class
        n_samples = n_clusters * 50

        for hypercube in (False, True):
            X, y = make(n_samples=n_samples, n_classes=n_classes,
                        weights=weights, n_features=n_informative,
                        n_informative=n_informative,
                        n_clusters_per_class=n_clusters_per_class,
                        hypercube=hypercube, random_state=0)

            assert_equal(X.shape, (n_samples, n_informative))
            assert_equal(y.shape, (n_samples,))

            # Cluster by sign, viewed as strings to allow uniquing
            signs = np.sign(X)
            signs = signs.view(dtype='|S{0}'.format(signs.strides[0]))
            unique_signs, cluster_index = np.unique(signs,
                                                    return_inverse=True)

            assert_equal(len(unique_signs), n_clusters,
                         "Wrong number of clusters, or not in distinct "
                         "quadrants")

            clusters_by_class = defaultdict(set)
            for cluster, cls in zip(cluster_index, y):
                clusters_by_class[cls].add(cluster)
            for clusters in clusters_by_class.values():
                assert_equal(len(clusters), n_clusters_per_class,
                             "Wrong number of clusters per class")
            assert_equal(len(clusters_by_class), n_classes,
                         "Wrong number of classes")

            assert_array_almost_equal(np.bincount(y) / len(y) // weights,
                                      [1] * n_classes,
                                      err_msg="Wrong number of samples "
                                              "per class")

            # Ensure on vertices of hypercube
            for cluster in range(len(unique_signs)):
                centroid = X[cluster_index == cluster].mean(axis=0)
                if hypercube:
                    assert_array_almost_equal(np.abs(centroid),
                                              [class_sep] * n_informative,
                                              decimal=0,
                                              err_msg="Clusters are not "
                                                      "centered on hypercube "
                                                      "vertices")
                else:
                    assert_raises(AssertionError,
                                  assert_array_almost_equal,
                                  np.abs(centroid),
                                  [class_sep] * n_informative,
                                  decimal=0,
                                  err_msg="Clusters should not be cenetered "
                                          "on hypercube vertices")

    assert_raises(ValueError, make, n_features=2, n_informative=2, n_classes=5,
                  n_clusters_per_class=1)
    assert_raises(ValueError, make, n_features=2, n_informative=2, n_classes=3,
                  n_clusters_per_class=2)


def test_make_multilabel_classification_return_sequences():
    for allow_unlabeled, min_length in zip((True, False), (0, 1)):
        X, Y = make_multilabel_classification(n_samples=100, n_features=20,
                                              n_classes=3, random_state=0,
                                              return_indicator=False,
                                              allow_unlabeled=allow_unlabeled)
        assert_equal(X.shape, (100, 20), "X shape mismatch")
        if not allow_unlabeled:
            assert_equal(max([max(y) for y in Y]), 2)
        assert_equal(min([len(y) for y in Y]), min_length)
        assert_true(max([len(y) for y in Y]) <= 3)


def test_make_multilabel_classification_return_indicator():
    for allow_unlabeled, min_length in zip((True, False), (0, 1)):
        X, Y = make_multilabel_classification(n_samples=25, n_features=20,
                                              n_classes=3, random_state=0,
                                              allow_unlabeled=allow_unlabeled)
        assert_equal(X.shape, (25, 20), "X shape mismatch")
        assert_equal(Y.shape, (25, 3), "Y shape mismatch")
        assert_true(np.all(np.sum(Y, axis=0) > min_length))

    # Also test return_distributions and return_indicator with True
    X2, Y2, p_c, p_w_c = make_multilabel_classification(
        n_samples=25, n_features=20, n_classes=3, random_state=0,
        allow_unlabeled=allow_unlabeled, return_distributions=True)

    assert_array_equal(X, X2)
    assert_array_equal(Y, Y2)
    assert_equal(p_c.shape, (3,))
    assert_almost_equal(p_c.sum(), 1)
    assert_equal(p_w_c.shape, (20, 3))
    assert_almost_equal(p_w_c.sum(axis=0), [1] * 3)

def test_make_multilabel_classification_return_indicator_sparse():
    for allow_unlabeled, min_length in zip((True, False), (0, 1)):
        X, Y = make_multilabel_classification(n_samples=25, n_features=20,
                                              n_classes=3, random_state=0,
                                              return_indicator='sparse',
                                              allow_unlabeled=allow_unlabeled)
        assert_equal(X.shape, (25, 20), "X shape mismatch")
        assert_equal(Y.shape, (25, 3), "Y shape mismatch")
        assert_true(sp.issparse(Y))

def test_make_hastie_10_2():
    X, y = make_hastie_10_2(n_samples=100, random_state=0)
    assert_equal(X.shape, (100, 10), "X shape mismatch")
    assert_equal(y.shape, (100,), "y shape mismatch")
    assert_equal(np.unique(y).shape, (2,), "Unexpected number of classes")


def test_make_regression():
    X, y, c = make_regression(n_samples=100, n_features=10, n_informative=3,
                              effective_rank=5, coef=True, bias=0.0,
                              noise=1.0, random_state=0)

    assert_equal(X.shape, (100, 10), "X shape mismatch")
    assert_equal(y.shape, (100,), "y shape mismatch")
    assert_equal(c.shape, (10,), "coef shape mismatch")
    assert_equal(sum(c != 0.0), 3, "Unexpected number of informative features")

    # Test that y ~= np.dot(X, c) + bias + N(0, 1.0).
    assert_almost_equal(np.std(y - np.dot(X, c)), 1.0, decimal=1)

    # Test with small number of features.
    X, y = make_regression(n_samples=100, n_features=1)  # n_informative=3
    assert_equal(X.shape, (100, 1))


def test_make_regression_multitarget():
    X, y, c = make_regression(n_samples=100, n_features=10, n_informative=3,
                              n_targets=3, coef=True, noise=1., random_state=0)

    assert_equal(X.shape, (100, 10), "X shape mismatch")
    assert_equal(y.shape, (100, 3), "y shape mismatch")
    assert_equal(c.shape, (10, 3), "coef shape mismatch")
    assert_array_equal(sum(c != 0.0), 3,
                       "Unexpected number of informative features")

    # Test that y ~= np.dot(X, c) + bias + N(0, 1.0)
    assert_almost_equal(np.std(y - np.dot(X, c)), 1.0, decimal=1)


def test_make_blobs():
    cluster_stds = np.array([0.05, 0.2, 0.4])
    cluster_centers = np.array([[0.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    X, y = make_blobs(random_state=0, n_samples=50, n_features=2,
                      centers=cluster_centers, cluster_std=cluster_stds)

    assert_equal(X.shape, (50, 2), "X shape mismatch")
    assert_equal(y.shape, (50,), "y shape mismatch")
    assert_equal(np.unique(y).shape, (3,), "Unexpected number of blobs")
    for i, (ctr, std) in enumerate(zip(cluster_centers, cluster_stds)):
        assert_almost_equal((X[y == i] - ctr).std(), std, 1, "Unexpected std")


def test_make_friedman1():
    X, y = make_friedman1(n_samples=5, n_features=10, noise=0.0,
                          random_state=0)

    assert_equal(X.shape, (5, 10), "X shape mismatch")
    assert_equal(y.shape, (5,), "y shape mismatch")

    assert_array_almost_equal(y,
                              10 * np.sin(np.pi * X[:, 0] * X[:, 1])
                              + 20 * (X[:, 2] - 0.5) ** 2
                              + 10 * X[:, 3] + 5 * X[:, 4])


def test_make_friedman2():
    X, y = make_friedman2(n_samples=5, noise=0.0, random_state=0)

    assert_equal(X.shape, (5, 4), "X shape mismatch")
    assert_equal(y.shape, (5,), "y shape mismatch")

    assert_array_almost_equal(y,
                              (X[:, 0] ** 2
                               + (X[:, 1] * X[:, 2] - 1
                                  / (X[:, 1] * X[:, 3])) ** 2) ** 0.5)


def test_make_friedman3():
    X, y = make_friedman3(n_samples=5, noise=0.0, random_state=0)

    assert_equal(X.shape, (5, 4), "X shape mismatch")
    assert_equal(y.shape, (5,), "y shape mismatch")

    assert_array_almost_equal(y, np.arctan((X[:, 1] * X[:, 2]
                                            - 1 / (X[:, 1] * X[:, 3]))
                                           / X[:, 0]))


def test_make_low_rank_matrix():
    X = make_low_rank_matrix(n_samples=50, n_features=25, effective_rank=5,
                             tail_strength=0.01, random_state=0)

    assert_equal(X.shape, (50, 25), "X shape mismatch")

    from numpy.linalg import svd
    u, s, v = svd(X)
    assert_less(sum(s) - 5, 0.1, "X rank is not approximately 5")


def test_make_sparse_coded_signal():
    Y, D, X = make_sparse_coded_signal(n_samples=5, n_components=8,
                                       n_features=10, n_nonzero_coefs=3,
                                       random_state=0)
    assert_equal(Y.shape, (10, 5), "Y shape mismatch")
    assert_equal(D.shape, (10, 8), "D shape mismatch")
    assert_equal(X.shape, (8, 5), "X shape mismatch")
    for col in X.T:
        assert_equal(len(np.flatnonzero(col)), 3, 'Non-zero coefs mismatch')
    assert_array_almost_equal(np.dot(D, X), Y)
    assert_array_almost_equal(np.sqrt((D ** 2).sum(axis=0)),
                              np.ones(D.shape[1]))


def test_make_sparse_uncorrelated():
    X, y = make_sparse_uncorrelated(n_samples=5, n_features=10, random_state=0)

    assert_equal(X.shape, (5, 10), "X shape mismatch")
    assert_equal(y.shape, (5,), "y shape mismatch")


def test_make_spd_matrix():
    X = make_spd_matrix(n_dim=5, random_state=0)

    assert_equal(X.shape, (5, 5), "X shape mismatch")
    assert_array_almost_equal(X, X.T)

    from numpy.linalg import eig
    eigenvalues, _ = eig(X)
    assert_array_equal(eigenvalues > 0, np.array([True] * 5),
                       "X is not positive-definite")


def test_make_swiss_roll():
    X, t = make_swiss_roll(n_samples=5, noise=0.0, random_state=0)

    assert_equal(X.shape, (5, 3), "X shape mismatch")
    assert_equal(t.shape, (5,), "t shape mismatch")
    assert_array_almost_equal(X[:, 0], t * np.cos(t))
    assert_array_almost_equal(X[:, 2], t * np.sin(t))


def test_make_s_curve():
    X, t = make_s_curve(n_samples=5, noise=0.0, random_state=0)

    assert_equal(X.shape, (5, 3), "X shape mismatch")
    assert_equal(t.shape, (5,), "t shape mismatch")
    assert_array_almost_equal(X[:, 0], np.sin(t))
    assert_array_almost_equal(X[:, 2], np.sign(t) * (np.cos(t) - 1))


def test_make_biclusters():
    X, rows, cols = make_biclusters(
        shape=(100, 100), n_clusters=4, shuffle=True, random_state=0)
    assert_equal(X.shape, (100, 100), "X shape mismatch")
    assert_equal(rows.shape, (4, 100), "rows shape mismatch")
    assert_equal(cols.shape, (4, 100,), "columns shape mismatch")
    assert_all_finite(X)
    assert_all_finite(rows)
    assert_all_finite(cols)

    X2, _, _ = make_biclusters(shape=(100, 100), n_clusters=4,
                               shuffle=True, random_state=0)
    assert_array_almost_equal(X, X2)


def test_make_checkerboard():
    X, rows, cols = make_checkerboard(
        shape=(100, 100), n_clusters=(20, 5),
        shuffle=True, random_state=0)
    assert_equal(X.shape, (100, 100), "X shape mismatch")
    assert_equal(rows.shape, (100, 100), "rows shape mismatch")
    assert_equal(cols.shape, (100, 100,), "columns shape mismatch")

    X, rows, cols = make_checkerboard(
        shape=(100, 100), n_clusters=2, shuffle=True, random_state=0)
    assert_all_finite(X)
    assert_all_finite(rows)
    assert_all_finite(cols)

    X1, _, _ = make_checkerboard(shape=(100, 100), n_clusters=2,
                                 shuffle=True, random_state=0)
    X2, _, _ = make_checkerboard(shape=(100, 100), n_clusters=2,
                                 shuffle=True, random_state=0)
    assert_array_equal(X1, X2)


def test_make_moons():
    X, y = make_moons(3, shuffle=False)
    for x, label in zip(X, y):
        center = [0.0, 0.0] if label == 0 else [1.0, 0.5]
        dist_sqr = ((x - center) ** 2).sum()
        assert_almost_equal(dist_sqr, 1.0,
                            err_msg="Point is not on expected unit circle")
