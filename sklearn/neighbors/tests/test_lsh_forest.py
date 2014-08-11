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


def test_neighbors_accuracy_with_n_candidates():
    """Checks whether accuracy increases as `n_candidates` increases."""
    n_candidates_values = np.array([10, 50, 250])
    samples = 12
    dim = 2
    n_iter = 10
    n_points = 5
    accuracies = np.zeros(n_candidates_values.shape[0], dtype=float)
    X = np.random.rand(samples, dim)

    for i, n_candidates in enumerate(n_candidates_values):
        lshf = LSHForest(n_candidates=n_candidates)
        lshf.fit(X)
        for j in range(n_iter):
            point = X[np.random.randint(0, samples)]
            neighbors = lshf.kneighbors(point, n_neighbors=n_points)
            distances = euclidean_distances(point, X)
            ranks = np.argsort(distances)[0, :n_points]

            intersection = np.intersect1d(ranks, neighbors).shape[0]
            ratio = intersection/float(n_points)
            accuracies[i] = accuracies[i] + ratio

        accuracies[i] = accuracies[i]/float(n_iter)
    # Sorted accuracies should be equal to original accuracies
    assert_array_equal(accuracies, np.sort(accuracies),
                       err_msg="Accuracies are not non-decreasing.")


def test_neighbors_accuracy_with_n_estimators():
    """Checks whether accuracy increases as `n_estimators` increases."""
    n_estimators = np.array([1, 10, 100])
    samples = 12
    dim = 2
    n_iter = 10
    n_points = 20
    accuracies = np.zeros(n_estimators.shape[0], dtype=float)
    X = np.random.rand(samples, dim)

    for i, t in enumerate(n_estimators):
        lshf = LSHForest(n_candidates=500, n_estimators=t)
        lshf.fit(X)
        for j in range(n_iter):
            point = X[np.random.randint(0, samples)]
            neighbors = lshf.kneighbors(point, n_neighbors=n_points)
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
    """Checks whether desired number of neighbors are returned.

    It is guaranteed to return the requested number of neighbors
    if `min_hash_length` is set to 0. Returned distances should be
    in ascending order.
    """
    samples = 12
    dim = 2
    n_iter = 10
    X = np.random.rand(samples, dim)

    lshf = LSHForest(min_hash_length=0)
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
    """Checks whether Returned distances are less than `radius`

    At least one point should be returned when the `radius` is set
    to mean distance from the considering point to other points in
    the database.
    """
    samples = 12
    dim = 2
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

    # Multiple points
    n_points = 10
    points = X[np.random.randint(0, samples, n_points)]
    neighbors, distances = lshf.radius_neighbors(points,
                                                 return_distance=True)
    assert_equal(neighbors.shape[0], n_points)
    assert_equal(distances.shape[0], n_points)


def test_distances():
    """Checks whether returned distances are in ascending order."""
    samples = 12
    dim = 2
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
    """Checks whether `fit` method sets all attribute values correctly."""
    samples = 12
    dim = 2
    n_estimators = 5
    X = np.random.rand(samples, dim)

    lshf = LSHForest(n_estimators=n_estimators)

    lshf.fit(X)

    # _input_array = X
    assert_array_equal(X, lshf._fit_X)
    # A hash function g(p) for each tree
    assert_equal(n_estimators, lshf.hash_functions_.shape[0])
    # Hash length = 32
    assert_equal(32, lshf.hash_functions_.shape[1])
    # Number of trees in the forest
    assert_equal(n_estimators, len(lshf._trees))
    # Each tree has entries for every data point
    assert_equal(samples, len(lshf._trees[0]))
    # Original indices after sorting the hashes
    assert_equal(n_estimators, len(lshf._original_indices))
    # Each set of original indices in a tree has entries for every data point
    assert_equal(samples, len(lshf._original_indices[0]))


def test_partial_fit():
    """Checks whether inserting array is consitent with fitted data.

    `partial_fit` method should set all attribute values correctly.
    """
    samples = 12
    samples_partial_fit = 3
    dim = 2
    X = np.random.rand(samples, dim)
    X_partial_fit = np.random.rand(samples_partial_fit, dim)

    lshf = LSHForest()
    # Test unfitted estimator
    assert_raises(ValueError, lshf.partial_fit, X[0])

    lshf.fit(X)

    # Insert wrong dimension
    assert_raises(ValueError, lshf.partial_fit,
                  np.random.randn(samples_partial_fit, dim-1))

    lshf.partial_fit(X_partial_fit)

    # size of _input_array = samples + 1 after insertion
    assert_equal(lshf._fit_X.shape[0],
                 samples+samples_partial_fit)
    # size of _original_indices[1] = samples + 1
    assert_equal(len(lshf._original_indices[0]),
                 samples+samples_partial_fit)
    # size of _trees[1] = samples + 1
    assert_equal(len(lshf._trees[1]),
                 samples+samples_partial_fit)


if __name__ == "__main__":
    import nose
    nose.runmodule()
