"""
Testing for the approximate neighbor search using
Locality Sensitive Hashing Forest module
(sklearn.neighbors.LSHForest).
"""

# Author: Maheshakya Wijewardena, Joel Nothman

import numpy as np
import scipy.sparse as sp

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_array_less
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_not_equal
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import ignore_warnings

from sklearn.metrics.pairwise import pairwise_distances
from sklearn.neighbors import LSHForest
from sklearn.neighbors import NearestNeighbors


def test_neighbors_accuracy_with_n_candidates():
    # Checks whether accuracy increases as `n_candidates` increases.
    n_candidates_values = np.array([.1, 50, 500])
    n_samples = 100
    n_features = 10
    n_iter = 10
    n_points = 5
    rng = np.random.RandomState(42)
    accuracies = np.zeros(n_candidates_values.shape[0], dtype=float)
    X = rng.rand(n_samples, n_features)

    for i, n_candidates in enumerate(n_candidates_values):
        lshf = LSHForest(n_candidates=n_candidates)
        ignore_warnings(lshf.fit)(X)
        for j in range(n_iter):
            query = X[rng.randint(0, n_samples)].reshape(1, -1)

            neighbors = lshf.kneighbors(query, n_neighbors=n_points,
                                        return_distance=False)
            distances = pairwise_distances(query, X, metric='cosine')
            ranks = np.argsort(distances)[0, :n_points]

            intersection = np.intersect1d(ranks, neighbors).shape[0]
            ratio = intersection / float(n_points)
            accuracies[i] = accuracies[i] + ratio

        accuracies[i] = accuracies[i] / float(n_iter)
    # Sorted accuracies should be equal to original accuracies
    assert_true(np.all(np.diff(accuracies) >= 0),
                msg="Accuracies are not non-decreasing.")
    # Highest accuracy should be strictly greater than the lowest
    assert_true(np.ptp(accuracies) > 0,
                msg="Highest accuracy is not strictly greater than lowest.")


def test_neighbors_accuracy_with_n_estimators():
    # Checks whether accuracy increases as `n_estimators` increases.
    n_estimators = np.array([1, 10, 100])
    n_samples = 100
    n_features = 10
    n_iter = 10
    n_points = 5
    rng = np.random.RandomState(42)
    accuracies = np.zeros(n_estimators.shape[0], dtype=float)
    X = rng.rand(n_samples, n_features)

    for i, t in enumerate(n_estimators):
        lshf = LSHForest(n_candidates=500, n_estimators=t)
        ignore_warnings(lshf.fit)(X)
        for j in range(n_iter):
            query = X[rng.randint(0, n_samples)].reshape(1, -1)
            neighbors = lshf.kneighbors(query, n_neighbors=n_points,
                                        return_distance=False)
            distances = pairwise_distances(query, X, metric='cosine')
            ranks = np.argsort(distances)[0, :n_points]

            intersection = np.intersect1d(ranks, neighbors).shape[0]
            ratio = intersection / float(n_points)
            accuracies[i] = accuracies[i] + ratio

        accuracies[i] = accuracies[i] / float(n_iter)
    # Sorted accuracies should be equal to original accuracies
    assert_true(np.all(np.diff(accuracies) >= 0),
                msg="Accuracies are not non-decreasing.")
    # Highest accuracy should be strictly greater than the lowest
    assert_true(np.ptp(accuracies) > 0,
                msg="Highest accuracy is not strictly greater than lowest.")


@ignore_warnings
def test_kneighbors():
    # Checks whether desired number of neighbors are returned.
    # It is guaranteed to return the requested number of neighbors
    # if `min_hash_match` is set to 0. Returned distances should be
    # in ascending order.
    n_samples = 12
    n_features = 2
    n_iter = 10
    rng = np.random.RandomState(42)
    X = rng.rand(n_samples, n_features)

    lshf = LSHForest(min_hash_match=0)
    # Test unfitted estimator
    assert_raises(ValueError, lshf.kneighbors, X[0])

    ignore_warnings(lshf.fit)(X)

    for i in range(n_iter):
        n_neighbors = rng.randint(0, n_samples)
        query = X[rng.randint(0, n_samples)].reshape(1, -1)
        neighbors = lshf.kneighbors(query, n_neighbors=n_neighbors,
                                    return_distance=False)
        # Desired number of neighbors should be returned.
        assert_equal(neighbors.shape[1], n_neighbors)

    # Multiple points
    n_queries = 5
    queries = X[rng.randint(0, n_samples, n_queries)]
    distances, neighbors = lshf.kneighbors(queries,
                                           n_neighbors=1,
                                           return_distance=True)
    assert_equal(neighbors.shape[0], n_queries)
    assert_equal(distances.shape[0], n_queries)
    # Test only neighbors
    neighbors = lshf.kneighbors(queries, n_neighbors=1,
                                return_distance=False)
    assert_equal(neighbors.shape[0], n_queries)
    # Test random point(not in the data set)
    query = rng.randn(n_features).reshape(1, -1)
    lshf.kneighbors(query, n_neighbors=1,
                    return_distance=False)
    # Test n_neighbors at initialization
    neighbors = lshf.kneighbors(query, return_distance=False)
    assert_equal(neighbors.shape[1], 5)
    # Test `neighbors` has an integer dtype
    assert_true(neighbors.dtype.kind == 'i',
                msg="neighbors are not in integer dtype.")


def test_radius_neighbors():
    # Checks whether Returned distances are less than `radius`
    # At least one point should be returned when the `radius` is set
    # to mean distance from the considering point to other points in
    # the database.
    # Moreover, this test compares the radius neighbors of LSHForest
    # with the `sklearn.neighbors.NearestNeighbors`.
    n_samples = 12
    n_features = 2
    n_iter = 10
    rng = np.random.RandomState(42)
    X = rng.rand(n_samples, n_features)

    lshf = LSHForest()
    # Test unfitted estimator
    assert_raises(ValueError, lshf.radius_neighbors, X[0])

    ignore_warnings(lshf.fit)(X)

    for i in range(n_iter):
        # Select a random point in the dataset as the query
        query = X[rng.randint(0, n_samples)].reshape(1, -1)

        # At least one neighbor should be returned when the radius is the
        # mean distance from the query to the points of the dataset.
        mean_dist = np.mean(pairwise_distances(query, X, metric='cosine'))
        neighbors = lshf.radius_neighbors(query, radius=mean_dist,
                                          return_distance=False)

        assert_equal(neighbors.shape, (1,))
        assert_equal(neighbors.dtype, object)
        assert_greater(neighbors[0].shape[0], 0)
        # All distances to points in the results of the radius query should
        # be less than mean_dist
        distances, neighbors = lshf.radius_neighbors(query,
                                                     radius=mean_dist,
                                                     return_distance=True)
        assert_array_less(distances[0], mean_dist)

    # Multiple points
    n_queries = 5
    queries = X[rng.randint(0, n_samples, n_queries)]
    distances, neighbors = lshf.radius_neighbors(queries,
                                                 return_distance=True)

    # dists and inds should not be 1D arrays or arrays of variable lengths
    # hence the use of the object dtype.
    assert_equal(distances.shape, (n_queries,))
    assert_equal(distances.dtype, object)
    assert_equal(neighbors.shape, (n_queries,))
    assert_equal(neighbors.dtype, object)

    # Compare with exact neighbor search
    query = X[rng.randint(0, n_samples)].reshape(1, -1)
    mean_dist = np.mean(pairwise_distances(query, X, metric='cosine'))
    nbrs = NearestNeighbors(algorithm='brute', metric='cosine').fit(X)

    distances_exact, _ = nbrs.radius_neighbors(query, radius=mean_dist)
    distances_approx, _ = lshf.radius_neighbors(query, radius=mean_dist)

    # Radius-based queries do not sort the result points and the order
    # depends on the method, the random_state and the dataset order. Therefore
    # we need to sort the results ourselves before performing any comparison.
    sorted_dists_exact = np.sort(distances_exact[0])
    sorted_dists_approx = np.sort(distances_approx[0])

    # Distances to exact neighbors are less than or equal to approximate
    # counterparts as the approximate radius query might have missed some
    # closer neighbors.
    assert_true(np.all(np.less_equal(sorted_dists_exact,
                                     sorted_dists_approx)))


@ignore_warnings
def test_radius_neighbors_boundary_handling():
    X = [[0.999, 0.001], [0.5, 0.5], [0, 1.], [-1., 0.001]]
    n_points = len(X)

    # Build an exact nearest neighbors model as reference model to ensure
    # consistency between exact and approximate methods
    nnbrs = NearestNeighbors(algorithm='brute', metric='cosine').fit(X)

    # Build a LSHForest model with hyperparameter values that always guarantee
    # exact results on this toy dataset.
    lsfh = LSHForest(min_hash_match=0, n_candidates=n_points).fit(X)

    # define a query aligned with the first axis
    query = [[1., 0.]]

    # Compute the exact cosine distances of the query to the four points of
    # the dataset
    dists = pairwise_distances(query, X, metric='cosine').ravel()

    # The first point is almost aligned with the query (very small angle),
    # the cosine distance should therefore be almost null:
    assert_almost_equal(dists[0], 0, decimal=5)

    # The second point form an angle of 45 degrees to the query vector
    assert_almost_equal(dists[1], 1 - np.cos(np.pi / 4))

    # The third point is orthogonal from the query vector hence at a distance
    # exactly one:
    assert_almost_equal(dists[2], 1)

    # The last point is almost colinear but with opposite sign to the query
    # therefore it has a cosine 'distance' very close to the maximum possible
    # value of 2.
    assert_almost_equal(dists[3], 2, decimal=5)

    # If we query with a radius of one, all the samples except the last sample
    # should be included in the results. This means that the third sample
    # is lying on the boundary of the radius query:
    exact_dists, exact_idx = nnbrs.radius_neighbors(query, radius=1)
    approx_dists, approx_idx = lsfh.radius_neighbors(query, radius=1)

    assert_array_equal(np.sort(exact_idx[0]), [0, 1, 2])
    assert_array_equal(np.sort(approx_idx[0]), [0, 1, 2])
    assert_array_almost_equal(np.sort(exact_dists[0]), dists[:-1])
    assert_array_almost_equal(np.sort(approx_dists[0]), dists[:-1])

    # If we perform the same query with a slightly lower radius, the third
    # point of the dataset that lay on the boundary of the previous query
    # is now rejected:
    eps = np.finfo(np.float64).eps
    exact_dists, exact_idx = nnbrs.radius_neighbors(query, radius=1 - eps)
    approx_dists, approx_idx = lsfh.radius_neighbors(query, radius=1 - eps)

    assert_array_equal(np.sort(exact_idx[0]), [0, 1])
    assert_array_equal(np.sort(approx_idx[0]), [0, 1])
    assert_array_almost_equal(np.sort(exact_dists[0]), dists[:-2])
    assert_array_almost_equal(np.sort(approx_dists[0]), dists[:-2])


def test_distances():
    # Checks whether returned neighbors are from closest to farthest.
    n_samples = 12
    n_features = 2
    n_iter = 10
    rng = np.random.RandomState(42)
    X = rng.rand(n_samples, n_features)

    lshf = LSHForest()
    ignore_warnings(lshf.fit)(X)

    for i in range(n_iter):
        n_neighbors = rng.randint(0, n_samples)
        query = X[rng.randint(0, n_samples)].reshape(1, -1)
        distances, neighbors = lshf.kneighbors(query,
                                               n_neighbors=n_neighbors,
                                               return_distance=True)

        # Returned neighbors should be from closest to farthest, that is
        # increasing distance values.
        assert_true(np.all(np.diff(distances[0]) >= 0))

        # Note: the radius_neighbors method does not guarantee the order of
        # the results.


def test_fit():
    # Checks whether `fit` method sets all attribute values correctly.
    n_samples = 12
    n_features = 2
    n_estimators = 5
    rng = np.random.RandomState(42)
    X = rng.rand(n_samples, n_features)

    lshf = LSHForest(n_estimators=n_estimators)
    ignore_warnings(lshf.fit)(X)

    # _input_array = X
    assert_array_equal(X, lshf._fit_X)
    # A hash function g(p) for each tree
    assert_equal(n_estimators, len(lshf.hash_functions_))
    # Hash length = 32
    assert_equal(32, lshf.hash_functions_[0].components_.shape[0])
    # Number of trees_ in the forest
    assert_equal(n_estimators, len(lshf.trees_))
    # Each tree has entries for every data point
    assert_equal(n_samples, len(lshf.trees_[0]))
    # Original indices after sorting the hashes
    assert_equal(n_estimators, len(lshf.original_indices_))
    # Each set of original indices in a tree has entries for every data point
    assert_equal(n_samples, len(lshf.original_indices_[0]))


def test_partial_fit():
    # Checks whether inserting array is consistent with fitted data.
    # `partial_fit` method should set all attribute values correctly.
    n_samples = 12
    n_samples_partial_fit = 3
    n_features = 2
    rng = np.random.RandomState(42)
    X = rng.rand(n_samples, n_features)
    X_partial_fit = rng.rand(n_samples_partial_fit, n_features)

    lshf = LSHForest()

    # Test unfitted estimator
    ignore_warnings(lshf.partial_fit)(X)
    assert_array_equal(X, lshf._fit_X)

    ignore_warnings(lshf.fit)(X)

    # Insert wrong dimension
    assert_raises(ValueError, lshf.partial_fit,
                  np.random.randn(n_samples_partial_fit, n_features - 1))

    ignore_warnings(lshf.partial_fit)(X_partial_fit)

    # size of _input_array = samples + 1 after insertion
    assert_equal(lshf._fit_X.shape[0],
                 n_samples + n_samples_partial_fit)
    # size of original_indices_[1] = samples + 1
    assert_equal(len(lshf.original_indices_[0]),
                 n_samples + n_samples_partial_fit)
    # size of trees_[1] = samples + 1
    assert_equal(len(lshf.trees_[1]),
                 n_samples + n_samples_partial_fit)


def test_hash_functions():
    # Checks randomness of hash functions.
    # Variance and mean of each hash function (projection vector)
    # should be different from flattened array of hash functions.
    # If hash functions are not randomly built (seeded with
    # same value), variances and means of all functions are equal.
    n_samples = 12
    n_features = 2
    n_estimators = 5
    rng = np.random.RandomState(42)
    X = rng.rand(n_samples, n_features)

    lshf = LSHForest(n_estimators=n_estimators,
                     random_state=rng.randint(0, np.iinfo(np.int32).max))
    ignore_warnings(lshf.fit)(X)

    hash_functions = []
    for i in range(n_estimators):
        hash_functions.append(lshf.hash_functions_[i].components_)

    for i in range(n_estimators):
        assert_not_equal(np.var(hash_functions),
                         np.var(lshf.hash_functions_[i].components_))

    for i in range(n_estimators):
        assert_not_equal(np.mean(hash_functions),
                         np.mean(lshf.hash_functions_[i].components_))


def test_candidates():
    # Checks whether candidates are sufficient.
    # This should handle the cases when number of candidates is 0.
    # User should be warned when number of candidates is less than
    # requested number of neighbors.
    X_train = np.array([[5, 5, 2], [21, 5, 5], [1, 1, 1], [8, 9, 1],
                        [6, 10, 2]], dtype=np.float32)
    X_test = np.array([7, 10, 3], dtype=np.float32).reshape(1, -1)

    # For zero candidates
    lshf = LSHForest(min_hash_match=32)
    ignore_warnings(lshf.fit)(X_train)

    message = ("Number of candidates is not sufficient to retrieve"
               " %i neighbors with"
               " min_hash_match = %i. Candidates are filled up"
               " uniformly from unselected"
               " indices." % (3, 32))
    assert_warns_message(UserWarning, message, lshf.kneighbors,
                         X_test, n_neighbors=3)
    distances, neighbors = lshf.kneighbors(X_test, n_neighbors=3)
    assert_equal(distances.shape[1], 3)

    # For candidates less than n_neighbors
    lshf = LSHForest(min_hash_match=31)
    ignore_warnings(lshf.fit)(X_train)

    message = ("Number of candidates is not sufficient to retrieve"
               " %i neighbors with"
               " min_hash_match = %i. Candidates are filled up"
               " uniformly from unselected"
               " indices." % (5, 31))
    assert_warns_message(UserWarning, message, lshf.kneighbors,
                         X_test, n_neighbors=5)
    distances, neighbors = lshf.kneighbors(X_test, n_neighbors=5)
    assert_equal(distances.shape[1], 5)


def test_graphs():
    # Smoke tests for graph methods.
    n_samples_sizes = [5, 10, 20]
    n_features = 3
    rng = np.random.RandomState(42)

    for n_samples in n_samples_sizes:
        X = rng.rand(n_samples, n_features)
        lshf = LSHForest(min_hash_match=0)
        ignore_warnings(lshf.fit)(X)

        kneighbors_graph = lshf.kneighbors_graph(X)
        radius_neighbors_graph = lshf.radius_neighbors_graph(X)

        assert_equal(kneighbors_graph.shape[0], n_samples)
        assert_equal(kneighbors_graph.shape[1], n_samples)
        assert_equal(radius_neighbors_graph.shape[0], n_samples)
        assert_equal(radius_neighbors_graph.shape[1], n_samples)


def test_sparse_input():
    # note: Fixed random state in sp.rand is not supported in older scipy.
    #       The test should succeed regardless.
    X1 = sp.rand(50, 100)
    X2 = sp.rand(10, 100)
    forest_sparse = LSHForest(radius=1, random_state=0).fit(X1)
    forest_dense = LSHForest(radius=1, random_state=0).fit(X1.A)

    d_sparse, i_sparse = forest_sparse.kneighbors(X2, return_distance=True)
    d_dense, i_dense = forest_dense.kneighbors(X2.A, return_distance=True)

    assert_almost_equal(d_sparse, d_dense)
    assert_almost_equal(i_sparse, i_dense)

    d_sparse, i_sparse = forest_sparse.radius_neighbors(X2,
                                                        return_distance=True)
    d_dense, i_dense = forest_dense.radius_neighbors(X2.A,
                                                     return_distance=True)
    assert_equal(d_sparse.shape, d_dense.shape)
    for a, b in zip(d_sparse, d_dense):
        assert_almost_equal(a, b)
    for a, b in zip(i_sparse, i_dense):
        assert_almost_equal(a, b)
