from re import search
from warnings import simplefilter

import numpy as np
import scipy.sparse as sp

from numpy.testing import assert_array_equal, assert_array_almost_equal

from sklearn.cluster import BisectKMeans

import pytest


@pytest.mark.parametrize("bisect_strategy", ["biggest_sse", "largest_cluster"])
def test_three_clusters(bisect_strategy):
    """Tries to perform bisect k-means for three clusters to check
    if splitting data is performed correctly
    """

    # X = np.array([[1, 2], [1, 4], [1, 0],
    #               [10, 2], [10, 4], [10, 0],
    #               [10, 6], [10, 8], [10, 10]])

    # X[0][1] swapped with X[1][1] intentionally for checking labeling
    X = np.array(
        [[1, 2], [10, 4], [1, 0], [10, 2], [1, 4], [10, 0], [10, 6], [10, 8], [10, 10]]
    )
    bisect_means = BisectKMeans(
        n_clusters=3, random_state=0, bisect_strategy=bisect_strategy
    )
    bisect_means.fit(X)

    expected_centers = [[1, 2], [10, 2], [10, 8]]
    expected_predict = [0, 1]
    expected_labels = [0, 1, 0, 1, 0, 1, 2, 2, 2]

    assert_array_equal(expected_centers, bisect_means.cluster_centers_)
    assert_array_equal(expected_predict, bisect_means.predict([[0, 0], [12, 3]]))
    assert_array_equal(expected_labels, bisect_means.labels_)


def test_sparse():
    """Test Bisecting K-Means with sparse data
    Also test if results obtained from fit(X) are the same as fit(X).predict(X)
    """

    rng = np.random.RandomState(0)

    X = rng.rand(20, 2)
    X[X < 0.8] = 0
    X_csr = sp.csr_matrix(X)

    bisect_means = BisectKMeans(n_clusters=3, random_state=0)
    bisect_means.fit(X_csr)

    # Check if labels obtained from fit(X) are equal to fit(X).predict(X)
    # X_csr is passed here to check also predict for sparse
    assert_array_almost_equal(bisect_means.labels_, bisect_means.predict(X_csr))

    sparse_centers = bisect_means.cluster_centers_

    bisect_means.fit(X)
    normal_centers = bisect_means.cluster_centers_

    # Check if results is the same for dense and sparse data
    assert_array_almost_equal(normal_centers, sparse_centers)


@pytest.mark.parametrize("n_clusters", [4, 5])
def test_n_clusters(n_clusters):
    """Test if resulting labels are in range [0, n_clusters - 1]"""

    rng = np.random.RandomState(0)
    X = rng.rand(10, 2)

    bisect_means = BisectKMeans(n_clusters=n_clusters, random_state=0)
    bisect_means.fit(X)

    assert_array_equal(np.unique(bisect_means.labels_), np.arange(n_clusters))


def test_one_cluster():
    """Test warnings and performance for n_cluster = 1"""

    X = np.array([[1, 2], [10, 2], [10, 8]])

    with pytest.warns(None) as w:
        bisect_means = BisectKMeans(n_clusters=1, random_state=0)
        bisect_means.fit(X)
        msg = (
            "BisectKMeans might be inefficient for n_cluster smaller than 3 "
            + " - Use Normal KMeans from sklearn.cluster instead."
        )
        assert str(w[0].message) == msg

        msg = "Bisection won't be performed - needs at least two clusters to run."
        assert str(w[1].message) == msg

        # All labels from fit or predict should be equal 0
        assert all(bisect_means.predict(X) == 0)


@pytest.mark.parametrize(
    "param, match, single_value",
    [
        # Test bisect_strategy param
        (
            {"bisect_strategy": "None"},
            r"Bisect Strategy must be 'biggest_sse', or 'largest_cluster'",
            False,
        ),
        # Test init array
        (
            {"init": np.ones((5, 2))},
            "BisectKMeans does not support init as array.",
            False,
        ),
        # Test single X value
        (
            {"n_clusters": 1},
            "BisectKMeans needs more than one sample to perform bisection.",
            True,
        ),
    ],
)
def test_wrong_params(param, match, single_value):
    """Test Exceptions at check_params function"""

    simplefilter("ignore")

    if single_value:
        X = np.ones((1, 1))
    else:
        rng = np.random.RandomState(0)
        X = rng.rand(5, 2)

    with pytest.raises(ValueError, match=match):
        bisect_means = BisectKMeans(n_clusters=3, n_init=1)
        bisect_means.set_params(**param)
        bisect_means.fit(X)


def test_verbose(capsys):
    """Test Verbose mode"""
    rng = np.random.RandomState(0)
    X = rng.rand(5, 2)

    bisect_means = BisectKMeans(n_clusters=3, verbose=1)
    bisect_means.fit(X)

    captured = capsys.readouterr()

    assert search(r"Running Bisecting K-Means", captured.out)
