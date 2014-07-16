"""
Testing for the Locality Sensitive Hashing Forest
module (sklearn.neighbors.LSHForest).
"""

# Author: Maheshakya Wijewardena

import numpy as np

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_warns

from sklearn.metrics import euclidean_distances
from sklearn.neighbors import LSHForest


def test_neighbors_accuracy_with_c():
    """Accuracy increases as `c` increases."""
    c_values = np.array([10, 50, 250])
    samples = 1000
    dim = 50
    n_iter = 10
    n_points = 20
    accuracies = np.zeros(c_values.shape[0], dtype=float)
    X = np.random.rand(samples, dim)

    for i, c in enumerate(c_values):
        lshf = LSHForest(c=c)
        lshf.fit(X)
        for j in range(n_iter):
            point = X[np.random.randint(0, samples)]
            neighbors = lshf.kneighbors(point, n_neighbors=n_points,
                                        return_distance=False)
            distances = euclidean_distances(point, X)
            ranks = np.argsort(distances)[0, :n_points]

            intersection = np.intersect1d(ranks, neighbors).shape[0]
            ratio = intersection/float(n_points)
            accuracies[i] = accuracies[i] + ratio

        accuracies[i] = accuracies[i]/float(n_iter)

    # Sorted accuracies should be equal to original accuracies
    assert_array_equal(accuracies, np.sort(accuracies),
                       err_msg="Accuracies are not non-decreasing.")


def test_neighbors_accuracy_with_n_trees():
    """Accuracy increases as `n_trees` increases."""
    n_trees = np.array([1, 10, 100])
    samples = 1000
    dim = 50
    n_iter = 10
    n_points = 20
    accuracies = np.zeros(n_trees.shape[0], dtype=float)
    X = np.random.rand(samples, dim)

    for i, t in enumerate(n_trees):
        lshf = LSHForest(c=500, n_trees=t)
        lshf.fit(X)
        for j in range(n_iter):
            point = X[np.random.randint(0, samples)]
            neighbors = lshf.kneighbors(point, n_neighbors=n_points,
                                        return_distance=False)
            distances = euclidean_distances(point, X)
            ranks = np.argsort(distances)[0, :n_points]

            intersection = np.intersect1d(ranks, neighbors).shape[0]
            ratio = intersection/float(n_points)
            accuracies[i] = accuracies[i] + ratio

        accuracies[i] = accuracies[i]/float(n_iter)

    # Sorted accuracies should be equal to original accuracies
    assert_array_equal(accuracies, np.sort(accuracies),
                       err_msg="Accuracies are not non-decreasing.")


def test_kneighbors():
    samples = 1000
    dim = 50
    n_iter = 100
    X = np.random.rand(samples, dim)

    lshf = LSHForest(lower_bound=0)
    lshf.fit(X)

    for i in range(n_iter):
        n_neighbors = np.random.randint(0, samples)
        point = X[np.random.randint(0, samples)]
        neighbors = lshf.kneighbors(point, n_neighbors=n_neighbors, 
                                    return_distance=False)
        # Desired number of neighbors should be returned.
        assert_equal(neighbors.shape[1], n_neighbors)

    # Test whether a value error is raised when X=None
    assert_raises(ValueError, lshf.kneighbors, None)


def test_distances():
    samples = 1000
    dim = 50
    n_iter = 100
    X = np.random.rand(samples, dim)

    lshf = LSHForest()
    lshf.fit(X)

    for i in range(n_iter):
        n_neighbors = np.random.randint(0, samples)
        point = X[np.random.randint(0, samples)]
        neighbors = lshf.kneighbors(point, n_neighbors=n_neighbors,
                                    return_distance=True)
        # Returned distances should be in sorted order.
        assert_array_equal(neighbors[1][0], np.sort(neighbors[1][0]))


def test_fit():
    samples = 1000
    dim = 50
    n_trees = 5
    X = np.random.rand(samples, dim)

    lshf = LSHForest(n_trees=n_trees)

    # Test whether a value error is raised when X=None
    assert_raises(ValueError, lshf.fit, None)

    lshf.fit(X)

    # _input_array = X
    assert_array_equal(X, lshf._input_array)

    # A hash function g(p) for each tree
    assert_equal(n_trees, lshf.hash_functions_.shape[0])

    # Hash length = 32
    assert_equal(32, lshf.hash_functions_.shape[1])

    # Number of trees in the forest
    assert_equal(n_trees, lshf._trees.shape[0])

    # Each tree has entries for every data point
    assert_equal(samples, lshf._trees.shape[1])


if __name__ == "__main__":
    import nose
    nose.runmodule()
