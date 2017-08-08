from __future__ import division
import numpy as np

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_false

from sklearn.preprocessing.knn_imputation import KNNImputer
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

    # Test with an imputable matrix
    X = np.array([
        [1, 0, 1, 0, 5],
        [2, 1, 2, 2, 3],
        [3, 2, 3, 0, 0],
        [6, 6, 0, 5, 13],
    ])

    statistics_mean = [3, 3, 2, 3.5, 7]
    X_imputed = np.array([
        [1, 1.5, 1,   2, 5],
        [2, 1,   2,   2, 3],
        [3, 2,   3,   2, 4],
        [6, 6,   1.5, 5, 13],
    ])

    assert_array_equal(imputer.fit(X).transform(X), X_imputed)
    assert_array_equal(imputer.statistics_[1], statistics_mean)


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
    assert_array_equal(imputer.fit(X).transform(X), X_imputed)
    assert_array_equal(imputer.statistics_[1], statistics_mean)

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
    assert_array_equal(imputer.fit(X).transform(X), X_imputed)
    assert_array_equal(imputer.statistics_[1], statistics_mean)

    # Test with weights = "distance"
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
    assert_array_almost_equal(imputer.fit(X).transform(X), X_imputed,
                              decimal=4)
    assert_array_equal(imputer.statistics_[1], statistics_mean)


def test_imputation_pickle():
    # Test for pickling imputers.
    import pickle

    l = 100
    X = np.random.rand(l, l+1)

    imputer = KNNImputer()
    imputer.fit(X)

    imputer_pickled = pickle.loads(pickle.dumps(imputer))

    assert_array_equal(imputer.transform(X.copy()),
                       imputer_pickled.transform(X.copy()),
                       "Fail to transform the data after pickling ")


def test_imputation_copy():
    # Test imputation with copy
    X_orig = sparse_random_matrix(10, 10, density=0.75, random_state=0)

    # copy=True, dense => copy
    X = X_orig.copy().toarray()
    imputer = KNNImputer(missing_values=0, copy=True)
    Xt = imputer.fit(X).transform(X)
    Xt[0, 0] = -1
    assert_false(np.all(X == Xt))

    # copy=False, dense => no copy
    X = X_orig.copy().toarray()
    imputer = KNNImputer(missing_values=0, copy=False)
    Xt = imputer.fit(X).transform(X)
    Xt[0, 0] = -1
    assert_array_equal(X, Xt)
