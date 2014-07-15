"""
Testing for the Locality Sensitive Hashing Forest
module (sklearn.neighbors.LSHForest).
"""

# Author: Gilles Louppe

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


mask_length = 7
left_masks = np.array([int('000000', 2), int('100000', 2),
                       int('110000', 2), int('111000', 2),
                       int('111100', 2), int('111110', 2),
                       int('111111', 2)])
right_masks = np.array([int('111111', 2), int('011111', 2),
                        int('001111', 2), int('000111', 2),
                        int('000011', 2), int('000001', 2),
                        int('000000', 2)])

test_strings = np.array(['000001', '001010', '001101',
                         '001110', '001111', '100001',
                         '101001', '101010'])

test_array = np.zeros(test_strings.shape[0], dtype=int)
for i in range(test_array.shape[0]):
    test_array[i] = int(test_strings[i], 2)

test_value = int('001110', 2)


def test_neighbors_accuracy_with_c():
    """Accuracy increases as `c` increases."""
    c_values = np.array([10, 50, 250])
    samples = 1000
    dim = 50
    n_iter = 10
    n_points = 20
    accuracies = np.zeros(c_values.shape[0], dtype=float)
    X = np.random.rand(samples, dim)

    for i in range(c_values.shape[0]):
        lshf = LSHForest(c=c_values[i])
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

    for i in range(n_trees.shape[0]):
        lshf = LSHForest(c=500, n_trees=n_trees[i])
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


if __name__ == "__main__":
    import nose
    nose.runmodule()
