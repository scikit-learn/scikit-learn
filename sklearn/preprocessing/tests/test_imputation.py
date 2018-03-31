from __future__ import division
import numpy as np
from scipy import sparse

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import ignore_warnings

from sklearn.preprocessing.imputation import Imputer
from sklearn.preprocessing.imputation import KNNImputer
from sklearn.metrics.pairwise import masked_euclidean_distances
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.neighbors import NearestNeighbors
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV
from sklearn import tree
from sklearn.random_projection import sparse_random_matrix


@ignore_warnings
def _check_statistics(X, X_true,
                      strategy, statistics, missing_values):
    """Utility function for testing imputation for a given strategy.

    Test:
        - along the two axes
        - with dense and sparse arrays

    Check that:
        - the statistics (mean, median, mode) are correct
        - the missing values are imputed correctly"""

    err_msg = "Parameters: strategy = %s, missing_values = %s, " \
              "axis = {0}, sparse = {1}" % (strategy, missing_values)

    assert_ae = assert_array_equal
    if X.dtype.kind == 'f' or X_true.dtype.kind == 'f':
        assert_ae = assert_array_almost_equal

    # Normal matrix, axis = 0
    imputer = Imputer(missing_values, strategy=strategy, axis=0)
    X_trans = imputer.fit(X).transform(X.copy())
    assert_ae(imputer.statistics_, statistics,
              err_msg=err_msg.format(0, False))
    assert_ae(X_trans, X_true, err_msg=err_msg.format(0, False))

    # Normal matrix, axis = 1
    imputer = Imputer(missing_values, strategy=strategy, axis=1)
    imputer.fit(X.transpose())
    if np.isnan(statistics).any():
        assert_raises(ValueError, imputer.transform, X.copy().transpose())
    else:
        X_trans = imputer.transform(X.copy().transpose())
        assert_ae(X_trans, X_true.transpose(),
                  err_msg=err_msg.format(1, False))

    # Sparse matrix, axis = 0
    imputer = Imputer(missing_values, strategy=strategy, axis=0)
    imputer.fit(sparse.csc_matrix(X))
    X_trans = imputer.transform(sparse.csc_matrix(X.copy()))

    if sparse.issparse(X_trans):
        X_trans = X_trans.toarray()

    assert_ae(imputer.statistics_, statistics,
              err_msg=err_msg.format(0, True))
    assert_ae(X_trans, X_true, err_msg=err_msg.format(0, True))

    # Sparse matrix, axis = 1
    imputer = Imputer(missing_values, strategy=strategy, axis=1)
    imputer.fit(sparse.csc_matrix(X.transpose()))
    if np.isnan(statistics).any():
        assert_raises(ValueError, imputer.transform,
                      sparse.csc_matrix(X.copy().transpose()))
    else:
        X_trans = imputer.transform(sparse.csc_matrix(X.copy().transpose()))

        if sparse.issparse(X_trans):
            X_trans = X_trans.toarray()

        assert_ae(X_trans, X_true.transpose(),
                  err_msg=err_msg.format(1, True))


@ignore_warnings
def test_imputation_shape():
    # Verify the shapes of the imputed matrix for different strategies.
    X = np.random.randn(10, 2)
    X[::2] = np.nan

    for strategy in ['mean', 'median', 'most_frequent']:
        imputer = Imputer(strategy=strategy)
        X_imputed = imputer.fit_transform(X)
        assert_equal(X_imputed.shape, (10, 2))
        X_imputed = imputer.fit_transform(sparse.csr_matrix(X))
        assert_equal(X_imputed.shape, (10, 2))


@ignore_warnings
def test_imputation_mean_median_only_zero():
    # Test imputation using the mean and median strategies, when
    # missing_values == 0.
    X = np.array([
        [np.nan, 0, 0, 0, 5],
        [np.nan, 1, 0, np.nan, 3],
        [np.nan, 2, 0, 0, 0],
        [np.nan, 6, 0, 5, 13],
    ])

    X_imputed_mean = np.array([
        [3, 5],
        [1, 3],
        [2, 7],
        [6, 13],
    ])
    statistics_mean = [np.nan, 3, np.nan, np.nan, 7]

    # Behaviour of median with NaN is undefined, e.g. different results in
    # np.median and np.ma.median
    X_for_median = X[:, [0, 1, 2, 4]]
    X_imputed_median = np.array([
        [2, 5],
        [1, 3],
        [2, 5],
        [6, 13],
    ])
    statistics_median = [np.nan, 2, np.nan, 5]

    _check_statistics(X, X_imputed_mean, "mean", statistics_mean, 0)
    _check_statistics(X_for_median, X_imputed_median, "median",
                      statistics_median, 0)


def safe_median(arr, *args, **kwargs):
    # np.median([]) raises a TypeError for numpy >= 1.10.1
    length = arr.size if hasattr(arr, 'size') else len(arr)
    return np.nan if length == 0 else np.median(arr, *args, **kwargs)


def safe_mean(arr, *args, **kwargs):
    # np.mean([]) raises a RuntimeWarning for numpy >= 1.10.1
    length = arr.size if hasattr(arr, 'size') else len(arr)
    return np.nan if length == 0 else np.mean(arr, *args, **kwargs)


@ignore_warnings
def test_imputation_mean_median():
    # Test imputation using the mean and median strategies, when
    # missing_values != 0.
    rng = np.random.RandomState(0)

    dim = 10
    dec = 10
    shape = (dim * dim, dim + dec)

    zeros = np.zeros(shape[0])
    values = np.arange(1, shape[0] + 1)
    values[4::2] = - values[4::2]

    tests = [("mean", "NaN", lambda z, v, p: safe_mean(np.hstack((z, v)))),
             ("mean", 0, lambda z, v, p: np.mean(v)),
             ("median", "NaN", lambda z, v, p: safe_median(np.hstack((z, v)))),
             ("median", 0, lambda z, v, p: np.median(v))]

    for strategy, test_missing_values, true_value_fun in tests:
        X = np.empty(shape)
        X_true = np.empty(shape)
        true_statistics = np.empty(shape[1])

        # Create a matrix X with columns
        #    - with only zeros,
        #    - with only missing values
        #    - with zeros, missing values and values
        # And a matrix X_true containing all true values
        for j in range(shape[1]):
            nb_zeros = (j - dec + 1 > 0) * (j - dec + 1) * (j - dec + 1)
            nb_missing_values = max(shape[0] + dec * dec
                                    - (j + dec) * (j + dec), 0)
            nb_values = shape[0] - nb_zeros - nb_missing_values

            z = zeros[:nb_zeros]
            p = np.repeat(test_missing_values, nb_missing_values)
            v = values[rng.permutation(len(values))[:nb_values]]

            true_statistics[j] = true_value_fun(z, v, p)

            # Create the columns
            X[:, j] = np.hstack((v, z, p))

            if 0 == test_missing_values:
                X_true[:, j] = np.hstack((v,
                                          np.repeat(
                                              true_statistics[j],
                                              nb_missing_values + nb_zeros)))
            else:
                X_true[:, j] = np.hstack((v,
                                          z,
                                          np.repeat(true_statistics[j],
                                                    nb_missing_values)))

            # Shuffle them the same way
            np.random.RandomState(j).shuffle(X[:, j])
            np.random.RandomState(j).shuffle(X_true[:, j])

        # Mean doesn't support columns containing NaNs, median does
        if strategy == "median":
            cols_to_keep = ~np.isnan(X_true).any(axis=0)
        else:
            cols_to_keep = ~np.isnan(X_true).all(axis=0)

        X_true = X_true[:, cols_to_keep]

        _check_statistics(X, X_true, strategy,
                          true_statistics, test_missing_values)


@ignore_warnings
def test_imputation_median_special_cases():
    # Test median imputation with sparse boundary cases
    X = np.array([
        [0, np.nan, np.nan],  # odd: implicit zero
        [5, np.nan, np.nan],  # odd: explicit nonzero
        [0, 0, np.nan],    # even: average two zeros
        [-5, 0, np.nan],   # even: avg zero and neg
        [0, 5, np.nan],    # even: avg zero and pos
        [4, 5, np.nan],    # even: avg nonzeros
        [-4, -5, np.nan],  # even: avg negatives
        [-1, 2, np.nan],   # even: crossing neg and pos
    ]).transpose()

    X_imputed_median = np.array([
        [0, 0, 0],
        [5, 5, 5],
        [0, 0, 0],
        [-5, 0, -2.5],
        [0, 5, 2.5],
        [4, 5, 4.5],
        [-4, -5, -4.5],
        [-1, 2, .5],
    ]).transpose()
    statistics_median = [0, 5, 0, -2.5, 2.5, 4.5, -4.5, .5]

    _check_statistics(X, X_imputed_median, "median",
                      statistics_median, 'NaN')


@ignore_warnings
def test_imputation_most_frequent():
    # Test imputation using the most-frequent strategy.
    X = np.array([
        [-1, -1, 0, 5],
        [-1, 2, -1, 3],
        [-1, 1, 3, -1],
        [-1, 2, 3, 7],
    ])

    X_true = np.array([
        [2, 0, 5],
        [2, 3, 3],
        [1, 3, 3],
        [2, 3, 7],
    ])

    # scipy.stats.mode, used in Imputer, doesn't return the first most
    # frequent as promised in the doc but the lowest most frequent. When this
    # test will fail after an update of scipy, Imputer will need to be updated
    # to be consistent with the new (correct) behaviour
    _check_statistics(X, X_true, "most_frequent", [np.nan, 2, 3, 3], -1)


@ignore_warnings
def test_imputation_pipeline_grid_search():
    # Test imputation within a pipeline + gridsearch.
    pipeline = Pipeline([('imputer', Imputer(missing_values=0)),
                         ('tree', tree.DecisionTreeRegressor(random_state=0))])

    parameters = {
        'imputer__strategy': ["mean", "median", "most_frequent"],
        'imputer__axis': [0, 1]
    }

    l = 100
    X = sparse_random_matrix(l, l, density=0.10)
    Y = sparse_random_matrix(l, 1, density=0.10).toarray()
    gs = GridSearchCV(pipeline, parameters)
    gs.fit(X, Y)


@ignore_warnings
def test_imputation_pickle():
    # Test for pickling imputers.
    import pickle

    l = 100
    X = sparse_random_matrix(l, l, density=0.10)

    for strategy in ["mean", "median", "most_frequent"]:
        imputer = Imputer(missing_values=0, strategy=strategy)
        imputer.fit(X)

        imputer_pickled = pickle.loads(pickle.dumps(imputer))

        assert_array_almost_equal(
            imputer.transform(X.copy()),
            imputer_pickled.transform(X.copy()),
            err_msg="Fail to transform the data after pickling "
            "(strategy = %s)" % (strategy)
        )


@ignore_warnings
def test_imputation_copy():
    # Test imputation with copy
    X_orig = sparse_random_matrix(10, 10, density=0.75, random_state=0)
    imputers = {Imputer: {"missing_values": 0, "strategy": "mean"},
                KNNImputer: {"missing_values": 0}}

    # copy=True, dense => copy
    # copy=False, dense => no copy
    for imputer_cls, params in imputers.items():
        for copy in [True, False]:
            X = X_orig.copy().toarray()
            params["copy"] = copy
            imputer = imputer_cls(**params)
            Xt = imputer.fit(X).transform(X)
            Xt[0, 0] = -1
            if copy:
                assert_false(np.all(X == Xt))
            else:
                assert_array_almost_equal(X, Xt)

    # copy=True, sparse csr => copy
    X = X_orig.copy()
    imputer = Imputer(missing_values=X.data[0], strategy="mean", copy=True)
    Xt = imputer.fit(X).transform(X)
    Xt.data[0] = -1
    assert_false(np.all(X.data == Xt.data))

    # copy=False, dense => no copy
    X = X_orig.copy().toarray()
    imputer = Imputer(missing_values=0, strategy="mean", copy=False)
    Xt = imputer.fit(X).transform(X)
    Xt[0, 0] = -1
    assert_array_almost_equal(X, Xt)

    # copy=False, sparse csr, axis=1 => no copy
    X = X_orig.copy()
    imputer = Imputer(missing_values=X.data[0], strategy="mean",
                      copy=False, axis=1)
    Xt = imputer.fit(X).transform(X)
    Xt.data[0] = -1
    assert_array_almost_equal(X.data, Xt.data)

    # copy=False, sparse csc, axis=0 => no copy
    X = X_orig.copy().tocsc()
    imputer = Imputer(missing_values=X.data[0], strategy="mean",
                      copy=False, axis=0)
    Xt = imputer.fit(X).transform(X)
    Xt.data[0] = -1
    assert_array_almost_equal(X.data, Xt.data)

    # copy=False, sparse csr, axis=0 => copy
    X = X_orig.copy()
    imputer = Imputer(missing_values=X.data[0], strategy="mean",
                      copy=False, axis=0)
    Xt = imputer.fit(X).transform(X)
    Xt.data[0] = -1
    assert_false(np.all(X.data == Xt.data))

    # copy=False, sparse csc, axis=1 => copy
    X = X_orig.copy().tocsc()
    imputer = Imputer(missing_values=X.data[0], strategy="mean",
                      copy=False, axis=1)
    Xt = imputer.fit(X).transform(X)
    Xt.data[0] = -1
    assert_false(np.all(X.data == Xt.data))

    # copy=False, sparse csr, axis=1, missing_values=0 => copy
    X = X_orig.copy()
    imputer = Imputer(missing_values=0, strategy="mean",
                      copy=False, axis=1)
    Xt = imputer.fit(X).transform(X)
    assert_false(sparse.issparse(Xt))

    # Note: If X is sparse and if missing_values=0, then a (dense) copy of X is
    # made, even if copy=False.


#############################################################################
# BEGIN KNNIMPUTER TEST


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
    imputer_nan = KNNImputer(missing_values="NaN",
                             n_neighbors=n_neighbors,
                             weights="uniform")

    # Test with missing_values=0 when NaN present
    X = np.array([
        [np.nan, 0, 0, 0, 5],
        [np.nan, 1, 0, np.nan, 3],
        [np.nan, 2, 0, 0, 0],
        [np.nan, 6, 0, 5, 13],
    ])
    msg = "Input contains NaN, infinity or a value too large for %r." % X.dtype
    assert_raise_message(ValueError, msg, imputer.fit, X)

    # Test with % zeros in column > col_max_missing
    X = np.array([
        [1, 0, 0, 0, 5],
        [2, 1, 0, 2, 3],
        [3, 2, 0, 0, 0],
        [4, 6, 0, 5, 13],
    ])
    msg = "Some column(s) have more than {}% missing values".format(
        imputer.col_max_missing * 100)
    assert_raise_message(ValueError, msg, imputer.fit, X)

    # Test with an imputable matrix and also compare with missing_values="NaN"
    X = np.array([
        [1, 0, 1, 0, 1.],
        [2, 1, 2, 2, 3],
        [3, 2, 3, 0, 0],
        [6, 6, 0, 5, 17],
    ])

    X_nan = np.array([
        [1, np.nan, 1,      np.nan, 1.],
        [2, 1,      2,      2,      3],
        [3, 2,      3,      np.nan, np.nan],
        [6, 6,      np.nan, 5,      17],
    ])
    statistics_mean = np.nanmean(X_nan, axis=0)

    X_imputed = np.array([
        [1, 1.5, 1,   2, 1.],
        [2, 1,   2,   2, 3],
        [3, 2,   3,   2, 2],
        [6, 6,   2.5, 5, 17],
    ])

    assert_array_equal(imputer.fit_transform(X), X_imputed)
    assert_array_equal(imputer.statistics_, statistics_mean)
    assert_array_equal(imputer.fit_transform(X), imputer_nan.fit_transform(
        X_nan))


def test_knn_imputation_default():
    # Test imputation with default parameter values

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
    statistics_mean = np.nanmean(X, axis=0)

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
    statistics_mean = np.nanmean(X, axis=0)

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
        [4, 4, 5, np.nan],
        [6, 7, 6, np.nan],
        [8, 8, 8, np.nan],
        [20, 20, 20, 20],
        [22, 22, 22, 22]
    ])
    statistics_mean = np.nanmean(X, axis=0)

    X_imputed = np.array([
        [1, 0, 0, 21],
        [2, 1, 2, 21],
        [3, 2, 3, 21],
        [4, 4, 5, 21],
        [6, 7, 6, 21],
        [8, 8, 8, 21],
        [20, 20, 20, 20],
        [22, 22, 22, 22]
    ])

    imputer = KNNImputer()
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
    statistics_mean = np.nanmean(X, axis=0)

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


def test_default_with_invalid_input():
    # Test imputation with default values and invalid input

    # Test with % missing in a column > col_max_missing
    X = np.array([
        [np.nan, 0, 0, 0, 5],
        [np.nan, 1, 0, np.nan, 3],
        [np.nan, 2, 0, 0, 0],
        [np.nan, 6, 0, 5, 13],
        [np.nan, 7, 0, 7, 8],
        [np.nan, 8, 0, 8, 9],
    ])
    imputer = KNNImputer()
    msg = "Some column(s) have more than {}% missing values".format(
        imputer.col_max_missing * 100)
    assert_raise_message(ValueError, msg, imputer.fit, X)

    # Test with insufficient number of neighbors
    X = np.array([
        [1, 1, 1, 2, np.nan],
        [2, 1, 2, 2, 3],
        [3, 2, 3, 3, 8],
        [6, 6, 2, 5, 13],
    ])
    msg = "There are only %d samples, but n_neighbors=%d." % \
          (X.shape[0], imputer.n_neighbors)
    assert_raise_message(ValueError, msg, imputer.fit, X)

    # Test with inf present
    X = np.array([
        [np.inf, 1, 1, 2, np.nan],
        [2, 1, 2, 2, 3],
        [3, 2, 3, 3, 8],
        [np.nan, 6, 0, 5, 13],
        [np.nan, 7, 0, 7, 8],
        [6, 6, 2, 5, 7],
    ])
    msg = "+/- inf values are not allowed."
    assert_raise_message(ValueError, msg, KNNImputer().fit, X)

    # Test with inf present in matrix passed in transform()
    X = np.array([
        [np.inf, 1, 1, 2, np.nan],
        [2, 1, 2, 2, 3],
        [3, 2, 3, 3, 8],
        [np.nan, 6, 0, 5, 13],
        [np.nan, 7, 0, 7, 8],
        [6, 6, 2, 5, 7],
    ])

    X_fit = np.array([
        [0, 1, 1, 2, np.nan],
        [2, 1, 2, 2, 3],
        [3, 2, 3, 3, 8],
        [np.nan, 6, 0, 5, 13],
        [np.nan, 7, 0, 7, 8],
        [6, 6, 2, 5, 7],
    ])
    msg = "+/- inf values are not allowed in data to be transformed."
    assert_raise_message(ValueError, msg, KNNImputer().fit(X_fit).transform, X)


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
    statistics_mean = np.nanmean(X, axis=0)

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
    imputer_plus1 = KNNImputer(n_neighbors=n_neighbors + 1)

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


def test_weight_uniform():
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
    def no_weight(dist=None):
        return None

    imputer = KNNImputer(weights=no_weight)
    assert_array_equal(imputer.fit_transform(X), X_imputed_uniform)


def test_weight_distance():
    X = np.array([
        [0, 0],
        [np.nan, 2],
        [4, 3],
        [5, 6],
        [7, 7],
        [9, 8],
        [11, 10]
    ])

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
        [3.850394,          2],
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
    statistics_mean = np.nanmean(X, axis=0)

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

    # Test with varying missingness patterns
    X = np.array([
        [1,         0,      0,  1],
        [0,         np.nan, 1,  np.nan],
        [1,         1,      1,  np.nan],
        [0,         1,      0,  0],
        [0,         np.nan, 1,  0],
        [1,         1,      1,  1],
        [10,        10,     10, 10],
    ])
    statistics_mean = np.nanmean(X, axis=0)

    # Get weights of donor neighbors
    dist = masked_euclidean_distances(X)
    row1_nbor_dists = dist[1, :6]
    row1_nbor_dists[np.array([1, 2, 4])] = np.inf  # Degenerate neighbors
    row1_nbor_wt = 1/row1_nbor_dists

    row2_nbor_dists = dist[2, :6]
    row2_nbor_dists[np.array([1, 2])] = np.inf  # Degenerate neighbors
    row2_nbor_wt = 1/row2_nbor_dists
    # A non-degenerate donor has zero distance so it's weight is 1 and
    # others have weight 0
    row2_nbor_wt[~np.isinf(row2_nbor_wt)] = 0
    row2_nbor_wt[np.isinf(row2_nbor_wt)] = 1

    row4_nbor_dists = dist[4, :6]
    row4_nbor_dists[np.array([1, 4])] = np.inf  # Degenerate neighbors
    row4_nbor_wt = 1/row4_nbor_dists

    # Collect donor values
    col1_donor_values = np.ma.masked_invalid(X[:6, 1].copy())
    col3_donor_values = np.ma.masked_invalid(X[:6, 3].copy())

    # Final imputed values
    r1c1_imp = np.ma.average(col1_donor_values, weights=row1_nbor_wt)
    r1c3_imp = np.ma.average(col3_donor_values, weights=row1_nbor_wt)
    r2c3_imp = np.ma.average(col3_donor_values, weights=row2_nbor_wt)
    r4c1_imp = np.ma.average(col1_donor_values, weights=row4_nbor_wt)

    X_imputed = np.array([
        [1,         0,          0,  1],
        [0,         r1c1_imp,   1,  r1c3_imp],
        [1,         1,          1,  r2c3_imp],
        [0,         1,          0,  0],
        [0,         r4c1_imp,   1,  0],
        [1,         1,          1,  1],
        [10,        10,         10, 10],
    ])

    imputer = KNNImputer(weights="distance")
    assert_array_almost_equal(imputer.fit_transform(X), X_imputed, decimal=6)
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


def test_callable_metric():

    # Define callable metric that returns the l1 norm:
    def custom_callable(x, y, missing_values="NaN"):
        x = np.ma.array(x, mask=np.isnan(x))
        y = np.ma.array(y, mask=np.isnan(y))
        dist = np.nansum(np.abs(x-y))
        return dist

    X = np.array([
        [4, 3, 3, np.nan],
        [6, 9, 6, 9],
        [4, 8, 6, 9],
        [np.nan, 9, 11, 10.]
    ])

    X_imputed = np.array([
        [4, 3, 3, 9],
        [6, 9, 6, 9],
        [4, 8, 6, 9],
        [5, 9, 11, 10.]
    ])

    imputer = KNNImputer(n_neighbors=2, metric=custom_callable)
    assert_array_equal(imputer.fit_transform(X), X_imputed)


def test_complete_features():

    # Test with use_complete=True
    X = np.array([
        [0,      np.nan,    0,       np.nan],
        [1,      1,         1,       np.nan],
        [2,      2,         np.nan,  2],
        [3,      3,         3,       3],
        [4,      4,         4,       4],
        [5,      5,         5,       5],
        [6,      6,         6,       6],
        [np.nan, 7,         7,       7]
    ])

    r0c1 = np.mean(X[1:6, 1])
    r0c3 = np.mean(X[2:-1, -1])
    r1c3 = np.mean(X[2:-1, -1])
    r2c2 = np.nanmean(X[:6, 2])
    r7c0 = np.mean(X[2:-1, 0])

    X_imputed = np.array([
        [0,     r0c1,   0,    r0c3],
        [1,     1,      1,    r1c3],
        [2,     2,      r2c2, 2],
        [3,     3,      3,    3],
        [4,     4,      4,    4],
        [5,     5,      5,    5],
        [6,     6,      6,    6],
        [r7c0,  7,      7,    7]
    ])

    imputer_comp = KNNImputer(use_complete=True)
    assert_array_almost_equal(imputer_comp.fit_transform(X), X_imputed)


def test_complete_features_weighted():

    # Test with use_complete=True
    X = np.array([
        [0,      0,     0,       np.nan],
        [1,      1,     1,       np.nan],
        [2,      2,     np.nan,  2],
        [3,      3,     3,       3],
        [4,      4,     4,       4],
        [5,      5,     5,       5],
        [6,      6,     6,       6],
        [np.nan, 7,     7,       7]
    ])

    dist = pairwise_distances(X,
                              metric="masked_euclidean",
                              squared=False)

    # Calculate weights
    r0c3_w = 1.0 / dist[0, 2:-1]
    r1c3_w = 1.0 / dist[1, 2:-1]
    r2c2_w = 1.0 / dist[2, (0, 1, 3, 4, 5)]
    r7c0_w = 1.0 / dist[7, 2:7]

    # Calculate weighted averages
    r0c3 = np.average(X[2:-1, -1], weights=r0c3_w)
    r1c3 = np.average(X[2:-1, -1], weights=r1c3_w)
    r2c2 = np.average(X[(0, 1, 3, 4, 5), 2], weights=r2c2_w)
    r7c0 = np.average(X[2:7, 0], weights=r7c0_w)

    X_imputed = np.array([
        [0,     0,  0,    r0c3],
        [1,     1,  1,    r1c3],
        [2,     2,  r2c2, 2],
        [3,     3,  3,    3],
        [4,     4,  4,    4],
        [5,     5,  5,    5],
        [6,     6,  6,    6],
        [r7c0,  7,  7,    7]
    ])

    imputer_comp_wt = KNNImputer(weights="distance", use_complete=True)
    assert_array_almost_equal(imputer_comp_wt.fit_transform(X), X_imputed)
