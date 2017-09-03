from __future__ import division
import numpy as np

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_false

from sklearn.preprocessing.imputation import KNNImputer
from sklearn.neighbors import NearestNeighbors
from sklearn.random_projection import sparse_random_matrix


def test_knn_imputation_shape():
    # Verify the shapes of the imputed matrix for different weights and
    # number of neighbors.
    n_rows = 10
    n_cols = 2
    X = np.random.rand(n_rows, n_cols)
    X[0, 0] = np.nan

    for weights in ['uniform', 'distance']:
        for n_neighbors in range(1, 6):
            imputer = KNNImputer(n_neighbors=n_neighbors, weights=weights)
            X_imputed = imputer.fit_transform(X)
            assert_equal(X_imputed.shape, (n_rows, n_cols))


def test_knn_imputation_zero():
    # Test imputation when missing_values == 0
    missing_values = 0
    n_neighbors = 2
    imputer = KNNImputer(missing_values=missing_values,
                         n_neighbors=n_neighbors,
                         weights="uniform")
    # imputer_nan = KNNImputer(missing_values=np.nan,
    #                          n_neighbors=n_neighbors,
    #                          weights="uniform")

    # Test with missing_values=0 when NaN present
    X = np.array([
        [np.nan, 0, 0, 0, 5],
        [np.nan, 1, 0, np.nan, 3],
        [np.nan, 2, 0, 0, 0],
        [np.nan, 6, 0, 5, 13],
    ])
    assert_raises(ValueError, imputer.fit, X)

    # Test with % zeros in column > col_max_missing
    X = np.array([
        [1, 0, 0, 0, 5],
        [2, 1, 0, 2, 3],
        [3, 2, 0, 0, 0],
        [4, 6, 0, 5, 13],
    ])
    assert_raises(ValueError, imputer.fit, X)

    # Test with an imputable matrix and also compare with missing_values="NaN"
    # X = np.array([
    #     [1, 0, 1, 0, 5],
    #     [2, 1, 2, 2, 3],
    #     [3, 2, 3, 0, 0],
    #     [6, 6, 0, 5, 13],
    # ])

    X = np.array([
        [1, 0, 1, 0, 1],
        [2, 1, 2, 2, 3],
        [3, 2, 3, 0, 0],
        [6, 6, 0, 5, 17],
    ])

    # X_nan = np.array([
    #     [1, np.nan, 1, np.nan, 1],
    #     [2, 1, 2, 2, 3],
    #     [3, 2, 3, np.nan, np.nan],
    #     [6, 6, np.nan, 5, 17],
    # ])

    statistics_mean = [3, 3, 2, 3.5, 7]
    X_imputed = np.array([
        [1, 1.5, 1,   2, 1],
        [2, 1,   2,   2, 3],
        [3, 2,   3,   2, 2],
        [6, 6,   1.5, 5, 17],
    ])

    assert_array_equal(imputer.fit_transform(X), X_imputed)
    assert_array_equal(imputer.statistics_, statistics_mean)
    # The following fails at the moment as NearestNeighbors object does not
    # pass missing_values=0 to pairwise_distances()
    # assert_array_equal(imputer.fit_transform(X), imputer_nan.fit_transform(
    #     X_nan))


def test_knn_imputation_default():
    # Test imputation with default values
    # imputer = KNNImputer()

    # Test with % missing in a column > col_max_missing
    X = np.array([
        [np.nan, 0, 0, 0, 5],
        [np.nan, 1, 0, np.nan, 3],
        [np.nan, 2, 0, 0, 0],
        [np.nan, 6, 0, 5, 13],
        [np.nan, 7, 0, 7, 8],
        [np.nan, 8, 0, 8, 9],
    ])
    assert_raises(ValueError, KNNImputer().fit, X)

    # Test with insufficient number of neighbors
    imputer = KNNImputer()
    X = np.array([
        [1, 1, 1, 2, np.nan],
        [2, 1, 2, 2, 3],
        [3, 2, 3, 3, 8],
        [6, 6, 2, 5, 13],
    ])
    assert_raises(ValueError, KNNImputer().fit, X)

    # Test with inf present
    X = np.array([
        [np.inf, 1, 1, 2, np.nan],
        [2, 1, 2, 2, 3],
        [3, 2, 3, 3, 8],
        [np.nan, 6, 0, 5, 13],
        [np.nan, 7, 0, 7, 8],
        [6, 6, 2, 5, 7],
    ])
    assert_raises(ValueError, KNNImputer().fit, X)

    # Test with an imputable matrix
    X = np.array([
        [1,      0,      0,      1],
        [2,      1,      2,      np.nan],
        [3,      2,      3,      np.nan],
        [np.nan, 4,      5,      5],
        [6,      np.nan, 6,      7],
        [8,      8,      8,      8],
        [16,     15,     18,    19],
    ])

    statistics_mean = [6, 5, 6, 8]

    X_imputed = np.array([
        [1,      0,      0,      1],
        [2,      1,      2,      5.25],
        [3,      2,      3,      5.25],
        [4,      4,      5,      5],
        [6,      3,      6,      7],
        [8,      8,      8,      8],
        [16,     15,     18,    19],
    ])

    imputer = KNNImputer()
    assert_array_equal(imputer.fit_transform(X), X_imputed)
    assert_array_equal(imputer.statistics_, statistics_mean)

    # Test with % missing in row > row_max_missing
    X = np.array([
        [1,      0,      0,      1],
        [2,      1,      2,      np.nan],
        [3,      2,      3,      np.nan],
        [np.nan, 4,      5,      5],
        [6,      np.nan, 6,      7],
        [8,      8,      8,      8],
        [np.nan, np.nan, np.nan, 19],
    ])

    statistics_mean = [4, 3, 4, 8]
    X_imputed = np.array([
        [1,      0,      0,      1],
        [2,      1,      2,      5.25],
        [3,      2,      3,      5.25],
        [4,      4,      5,      5],
        [6,      3,      6,      7],
        [8,      8,      8,      8],
        [4,      3,      4,      19],
    ])

    imputer = KNNImputer()
    assert_array_equal(imputer.fit_transform(X), X_imputed)
    assert_array_equal(imputer.statistics_, statistics_mean)

    # Test with all neighboring donors also having missing feature values
    X = np.array([
        [1, 0, 0, np.nan],
        [2, 1, 2, np.nan],
        [3, 2, 3, np.nan],
        [np.nan, 4, 5, 5],
        [6, np.nan, 6, 7],
        [8, 8, 8, 8],
        [np.nan, np.nan, np.nan, 20],
    ])

    statistics_mean = [4, 3, 4, 10]

    X_imputed = np.array([
        [1, 0, 0, 10],
        [2, 1, 2, 10],
        [3, 2, 3, 5],
        [4.5, 4, 5, 5],
        [6, 6, 6, 7],
        [8, 8, 8, 8],
        [4, 3, 4, 20],
    ])

    imputer = KNNImputer(n_neighbors=2)
    assert_array_equal(imputer.fit_transform(X), X_imputed)
    assert_array_equal(imputer.statistics_, statistics_mean)

    # Test when data in fit() and transform() are different
    X = np.array([
        [0,      0],
        [np.nan, 2],
        [4,      3],
        [5,      6],
        [7,      7],
        [9,      8],
        [11,     16]
    ])
    statistics_mean = [6, 6]

    Y = np.array([
        [1,      0],
        [3,      2],
        [4,      np.nan]
        ])

    Y_imputed = np.array([
        [1,      0],
        [3,      2],
        [4,      4.8]
        ])

    imputer = KNNImputer()
    assert_array_equal(imputer.fit(X).transform(Y), Y_imputed)
    assert_array_equal(imputer.statistics_, statistics_mean)


def test_knn_n_neighbors():

    X = np.array([
        [0,      0],
        [np.nan, 2],
        [4,      3],
        [5,      np.nan],
        [7,      7],
        [np.nan, 8],
        [14,      13]
    ])
    statistics_mean = [6, 5.5]

    # Test with 1 neighbor
    X_imputed_1NN = np.array([
        [0,      0],
        [4,      2],
        [4,      3],
        [5,      3],
        [7,      7],
        [7,      8],
        [14,     13]
    ])

    n_neighbors = 1
    imputer = KNNImputer(n_neighbors=n_neighbors)
    imputer_plus1 = KNNImputer(n_neighbors=n_neighbors+1)

    assert_array_equal(imputer.fit_transform(X), X_imputed_1NN)
    assert_array_equal(imputer.statistics_, statistics_mean)
    assert_array_equal(imputer.fit_transform(X), imputer_plus1.fit(
        X).transform(X))

    # Test with 6 neighbors
    X = np.array([
        [0,      0],
        [np.nan, 2],
        [4,      3],
        [5,      np.nan],
        [7,      7],
        [np.nan, 8],
        [14,      13]
    ])

    X_imputed_6NN = np.array([
        [0,      0],
        [6,      2],
        [4,      3],
        [5,      5.5],
        [7,      7],
        [6,      8],
        [14,     13]
    ])

    n_neighbors = 6
    imputer = KNNImputer(n_neighbors=6)
    imputer_plus1 = KNNImputer(n_neighbors=n_neighbors + 1)

    assert_array_equal(imputer.fit_transform(X), X_imputed_6NN)
    assert_array_equal(imputer.statistics_, statistics_mean)
    assert_array_equal(imputer.fit_transform(X), imputer_plus1.fit(
        X).transform(X))


def test_weight_type():
    X = np.array([
        [0,      0],
        [np.nan, 2],
        [4,      3],
        [5,      6],
        [7,      7],
        [9,      8],
        [11,     10]
    ])

    # Test with "uniform" weight (or unweighted)
    X_imputed_uniform = np.array([
        [0,      0],
        [5,      2],
        [4,      3],
        [5,      6],
        [7,      7],
        [9,      8],
        [11,     10]
    ])

    imputer = KNNImputer(weights="uniform")
    assert_array_equal(imputer.fit_transform(X), X_imputed_uniform)

    # Test with "callable" weight
    def no_weight(dist):
        return None

    imputer = KNNImputer(weights=no_weight)
    assert_array_equal(imputer.fit_transform(X), X_imputed_uniform)

    # Test with "distance" weight
    nn = NearestNeighbors(metric="masked_euclidean")
    nn.fit(X)
    # Get distance of "n_neighbors" neighbors of row 1
    dist, index = nn.kneighbors()
    dist = dist[1, :]
    index = index[1, :]
    weights = 1 / dist
    values = X[index, 0]
    imputed = np.dot(values, weights) / np.sum(weights)

    # Manual calculation
    X_imputed_distance1 = np.array([
        [0,                 0],
        [3.850393700787402, 2],
        [4,                 3],
        [5,                 6],
        [7,                 7],
        [9,                 8],
        [11,                10]
    ])

    # NearestNeighbor calculation
    X_imputed_distance2 = np.array([
        [0,                 0],
        [imputed,           2],
        [4,                 3],
        [5,                 6],
        [7,                 7],
        [9,                 8],
        [11,                10]
    ])

    imputer = KNNImputer(weights="distance")
    assert_array_almost_equal(imputer.fit_transform(X), X_imputed_distance1,
                              decimal=6)
    assert_array_almost_equal(imputer.fit_transform(X), X_imputed_distance2,
                              decimal=6)

    # Test with weights = "distance" and n_neighbors=2
    X = np.array([
        [np.nan, 0,      0],
        [2,      1,      2],
        [3,      2,      3],
        [4,      5,      5],
    ])

    statistics_mean = [3, 2, 2.5]

    X_imputed = np.array([
        [2.3828, 0,     0],
        [2,      1,     2],
        [3,      2,     3],
        [4,      5,     5],
    ])

    imputer = KNNImputer(n_neighbors=2, weights="distance")
    assert_array_almost_equal(imputer.fit_transform(X), X_imputed,
                              decimal=4)
    assert_array_equal(imputer.statistics_, statistics_mean)


def test_metric_type():
    X = np.array([
        [0,      0],
        [np.nan, 2],
        [4,      3],
        [5,      6],
        [7,      7],
        [9,      8],
        [11,     10]
    ])

    # Test with a metric type without NaN support
    imputer = KNNImputer(metric="euclidean")
    assert_raises(ValueError, imputer.fit, X)


def test_imputation_copy():
    # Test imputation with copy
    X_orig = sparse_random_matrix(10, 10, density=0.75, random_state=0)

    # copy=True, dense => copy
    X = X_orig.copy().toarray()
    imputer = KNNImputer(missing_values=0, copy=True)
    Xt = imputer.fit_transform(X)
    Xt[0, 0] = -1
    assert_false(np.all(X == Xt))

    # copy=False, dense => no copy
    X = X_orig.copy().toarray()
    imputer = KNNImputer(missing_values=0, copy=False)
    Xt = imputer.fit_transform(X)
    Xt[0, 0] = -1
    assert_array_equal(X, Xt)
