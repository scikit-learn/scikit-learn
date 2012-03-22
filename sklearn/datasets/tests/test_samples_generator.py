import numpy as np
from numpy.testing import assert_equal, assert_approx_equal, \
                          assert_array_almost_equal
from nose.tools import assert_true

from .. import make_classification
from .. import make_multilabel_classification
from .. import make_regression
from .. import make_blobs
from .. import make_friedman1
from .. import make_friedman2
from .. import make_friedman3
from .. import make_low_rank_matrix
from .. import make_sparse_coded_signal
from .. import make_sparse_uncorrelated
from .. import make_spd_matrix
from .. import make_swiss_roll
from .. import make_s_curve


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


def test_make_multilabel_classification():
    for allow_unlabeled, min_length in zip((True, False), (0, 1)):
        X, Y = make_multilabel_classification(n_samples=100, n_features=20,
                                              n_classes=3, random_state=0,
                                              allow_unlabeled=allow_unlabeled)
        assert_equal(X.shape, (100, 20), "X shape mismatch")
        if not allow_unlabeled:
            assert_equal(max([max(y) for y in Y]), 2)
        assert_equal(min([len(y) for y in Y]), min_length)
        assert_true(max([len(y) for y in Y]) <= 3)


def test_make_regression():
    X, y, c = make_regression(n_samples=100, n_features=10, n_informative=3,
                              effective_rank=5, coef=True, bias=0.0,
                              noise=1.0, random_state=0)

    assert_equal(X.shape, (100, 10), "X shape mismatch")
    assert_equal(y.shape, (100,), "y shape mismatch")
    assert_equal(c.shape, (10,), "coef shape mismatch")
    assert_equal(sum(c != 0.0), 3, "Unexpected number of informative features")

    # Test that y ~= np.dot(X, c) + bias + N(0, 1.0)
    assert_approx_equal(np.std(y - np.dot(X, c)), 1.0, significant=2)


def test_make_blobs():
    X, y = make_blobs(n_samples=50, n_features=2,
                      centers=[[0.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
                      random_state=0)

    assert_equal(X.shape, (50, 2), "X shape mismatch")
    assert_equal(y.shape, (50,), "y shape mismatch")
    assert_equal(np.unique(y).shape, (3,), "Unexpected number of blobs")


def test_make_friedman1():
    X, y = make_friedman1(n_samples=5, n_features=10, noise=0.0,
                          random_state=0)

    assert_equal(X.shape, (5, 10), "X shape mismatch")
    assert_equal(y.shape, (5,), "y shape mismatch")

    assert_array_almost_equal(y, 10 * np.sin(np.pi * X[:, 0] * X[:, 1])
                                 + 20 * (X[:, 2] - 0.5) ** 2 \
                                 + 10 * X[:, 3] + 5 * X[:, 4])


def test_make_friedman2():
    X, y = make_friedman2(n_samples=5, noise=0.0, random_state=0)

    assert_equal(X.shape, (5, 4), "X shape mismatch")
    assert_equal(y.shape, (5,), "y shape mismatch")

    assert_array_almost_equal(y, (X[:, 0] ** 2
                                 + (X[:, 1] * X[:, 2]
                                    - 1 / (X[:, 1] * X[:, 3])) ** 2) ** 0.5)


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
    assert_true(sum(s) - 5 < 0.1, "X rank is not approximately 5")


def test_make_sparse_coded_signal():
    Y, D, X = make_sparse_coded_signal(n_samples=5, n_components=8,
                                           n_features=10, n_nonzero_coefs=3,
                                           random_state=0)
    assert_equal(Y.shape, (10, 5), "Y shape mismatch")
    assert_equal(D.shape, (10, 8), "D shape mismatch")
    assert_equal(X.shape, (8, 5), "X shape mismatch")
    for col in X.T:
        assert_equal(len(np.flatnonzero(col)), 3, 'Non-zero coefs mismatch')
    assert_equal(np.dot(D, X), Y)
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
    assert_equal(eigenvalues > 0, np.array([True] * 5),
                 "X is not positive-definite")


def test_make_swiss_roll():
    X, t = make_swiss_roll(n_samples=5, noise=0.0, random_state=0)

    assert_equal(X.shape, (5, 3), "X shape mismatch")
    assert_equal(t.shape, (5,), "t shape mismatch")
    assert_equal(X[:, 0], t * np.cos(t))
    assert_equal(X[:, 2], t * np.sin(t))


def test_make_s_curve():
    X, t = make_s_curve(n_samples=5, noise=0.0, random_state=0)

    assert_equal(X.shape, (5, 3), "X shape mismatch")
    assert_equal(t.shape, (5,), "t shape mismatch")
    assert_equal(X[:, 0], np.sin(t))
    assert_equal(X[:, 2], np.sign(t) * (np.cos(t) - 1))
