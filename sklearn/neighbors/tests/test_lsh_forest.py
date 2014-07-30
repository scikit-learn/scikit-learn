"""
Testing for the Locality Sensitive Hashing Forest
module (sklearn.neighbors.LSHForest).
"""

# Author: Maheshakya Wijewardena

import numpy as np

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_array_less
from sklearn.utils.testing import assert_greater

from sklearn.metrics import euclidean_distances
from sklearn.neighbors import LSHForest


def test_neighbors_accuracy_with_c():
    """Accuracy increases as `c` increases."""
    c_values = np.array([10, 50, 250])
    samples = 100
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
    samples = 100
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
    samples = 100
    dim = 50
    n_iter = 10
    X = np.random.rand(samples, dim)

    lshf = LSHForest(lower_bound=0)
    # Test unfitted estimator
    assert_raises(ValueError, lshf.kneighbors, X[0])

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

    # Multiple points
    n_points = 10
    points = X[np.random.randint(0, samples, n_points)]
    neighbors, distances = lshf.kneighbors(points,
                                           n_neighbors=1,
                                           return_distance=True)
    assert_equal(neighbors.shape[0], n_points)
    assert_equal(distances.shape[0], n_points)
    # Test only neighbors
    neighbors = lshf.kneighbors(points, n_neighbors=1)
    assert_equal(neighbors.shape[0], n_points)
    # Test random point(not in the data set)
    point = np.random.randn(dim)
    lshf.kneighbors(point, n_neighbors=1,
                    return_distance=False)


def test_radius_neighbors():
    samples = 100
    dim = 50
    n_iter = 10
    X = np.random.rand(samples, dim)

    lshf = LSHForest()
    # Test unfitted estimator
    assert_raises(ValueError, lshf.radius_neighbors, X[0])

    lshf.fit(X)

    for i in range(n_iter):
        point = X[np.random.randint(0, samples)]
        mean_dist = np.mean(euclidean_distances(point, X))
        neighbors = lshf.radius_neighbors(point, radius=mean_dist)
        # At least one neighbor should be returned.
        assert_greater(neighbors.shape[1], 0)
        # All distances should be less than mean_dist
        neighbors, distances = lshf.radius_neighbors(point,
                                                     radius=mean_dist,
                                                     return_distance=True)
        assert_array_less(distances, mean_dist)

    # Test whether a value error is raised when X=None
    assert_raises(ValueError, lshf.radius_neighbors, None)

    # Multiple points
    n_points = 10
    points = X[np.random.randint(0, samples, n_points)]
    neighbors, distances = lshf.radius_neighbors(points,
                                                 return_distance=True)
    assert_equal(neighbors.shape[0], n_points)
    assert_equal(distances.shape[0], n_points)


def test_distances():
    samples = 100
    dim = 50
    n_iter = 10
    X = np.random.rand(samples, dim)

    lshf = LSHForest()
    lshf.fit(X)

    for i in range(n_iter):
        n_neighbors = np.random.randint(0, samples)
        point = X[np.random.randint(0, samples)]
        neighbors, distances = lshf.kneighbors(point,
                                               n_neighbors=n_neighbors,
                                               return_distance=True)
        # Returned distances should be in sorted order.
        assert_array_equal(distances[0], np.sort(distances[0]))

        mean_dist = np.mean(euclidean_distances(point, X))
        neighbors, distances = lshf.radius_neighbors(point,
                                                     radius=mean_dist,
                                                     return_distance=True)
        assert_array_less(distances, mean_dist)


def test_fit():
    samples = 100
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
    assert_equal(n_trees, len(lshf._trees))
    # Each tree has entries for every data point
    assert_equal(samples, len(lshf._trees[0]))
    # Original indices after sorting the hashes
    assert_equal(n_trees, len(lshf._original_indices))
    # Each set of original indices in a tree has entries for every data point
    assert_equal(samples, len(lshf._original_indices[0]))


def test_insert():
    samples = 100
    samples_insert = 10
    dim = 50
    X = np.random.rand(samples, dim)
    X_insert = np.random.rand(samples_insert, dim)

    lshf = LSHForest()
    # Test unfitted estimator
    assert_raises(ValueError, lshf.insert, X[0])

    lshf.fit(X)

    # Insert wrong dimension
    assert_raises(ValueError, lshf.insert,
                  np.random.randn(samples_insert, dim-1))

    lshf.insert(X_insert)

    # size of _input_array = samples + 1 after insertion
    assert_equal(lshf._input_array.shape[0],
                 samples+samples_insert)
    # size of _original_indices[1] = samples + 1
    assert_equal(len(lshf._original_indices[0]),
                 samples+samples_insert)
    # size of _trees[1] = samples + 1
    assert_equal(len(lshf._trees[1]),
                 samples+samples_insert)


if __name__ == "__main__":
    import nose
    nose.runmodule()
