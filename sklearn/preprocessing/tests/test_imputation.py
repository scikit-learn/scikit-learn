import warnings
import numpy as np
import numpy.linalg as la
from scipy import sparse

from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_false

from sklearn.utils.sparsefuncs import mean_variance_axis0
from sklearn.preprocessing import Binarizer
from sklearn.preprocessing import KernelCenterer
from sklearn.preprocessing import LabelBinarizer
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import Normalizer
from sklearn.preprocessing import normalize
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import scale
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import add_dummy_feature

from sklearn.preprocessing import Imputer
from sklearn.pipeline import Pipeline
from sklearn import grid_search
from sklearn import tree
from sklearn.random_projection import sparse_random_matrix

from sklearn import datasets
from sklearn.linear_model.stochastic_gradient import SGDClassifier

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
              "axis = %%s, sparse = %%s".format(strategy, missing_values)

    # Normal matrix, axis = 0
    imputer = Imputer(missing_values, strategy=strategy, axis=0)
    X_trans = imputer.fit(X).transform(X.copy())
    assert_array_equal(imputer.statistics_, statistics,
                       err_msg.format(0, False))
    assert_array_equal(X_trans, X_true, err_msg.format(0, False))

    # Normal matrix, axis = 1
    imputer = Imputer(missing_values, strategy=strategy, axis=1)
    imputer.fit(X.transpose())
    if np.isnan(statistics).any():
        assert_raises(ValueError, imputer.transform, X.copy().transpose())
    else:
        X_trans = imputer.transform(X.copy().transpose())
        assert_array_equal(imputer.statistics_, statistics,
                           err_msg.format(1, False))
        assert_array_equal(X_trans, X_true.transpose(),
                           err_msg.format(1, False))

    # Sparse matrix, axis = 0
    imputer = Imputer(missing_values, strategy=strategy, axis=0)
    imputer.fit(sparse.csc_matrix(X))
    X_trans = imputer.transform(sparse.csc_matrix(X.copy()))

    if sparse.issparse(X_trans):
        X_trans = X_trans.toarray()

    assert_array_equal(imputer.statistics_, statistics,
                       err_msg.format(0, True))
    assert_array_equal(X_trans, X_true, err_msg.format(0, True))

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

        assert_array_equal(imputer.statistics_, statistics,
                           err_msg.format(1, True))
        assert_array_equal(X_trans, X_true.transpose(),
                           err_msg.format(1, True))


def test_imputation_mean_median_only_zero():
    """Test imputation using the mean and median strategies, when
       missing_values == 0."""
    X = np.array([
        [np.nan, 0, 0,  0,  5],
        [np.nan, 1, 0,  np.nan,  3],
        [np.nan, 2, 0,  0, 0],
        [np.nan, 6, 0,  5,  13],
    ])

    X_imputed_mean = np.array([
        [3,  5],
        [1,  3],
        [2,  7],
        [6, 13],
    ])
    statistics_mean = [np.nan, 3, np.nan, np.nan, 7]

    X_imputed_median = np.array([
        [2, 5,  5],
        [1, np.nan,  3],
        [2, 5, 5],
        [6, 5,  13],
    ])
    statistics_median = [np.nan, 2, np.nan, 5, 5]

    _check_statistics(X, X_imputed_mean, "mean", statistics_mean, 0)
    _check_statistics(X, X_imputed_median, "median", statistics_median, 0)


def test_imputation_mean_median():
    """Test imputation using the mean and median strategies, when
       missing_values != 0."""
    rng = np.random.RandomState(0)

    dim = 10
    dec = 10
    shape = (dim * dim, dim + dec)

    zeros = np.zeros(shape[0])
    values = np.arange(1, shape[0]+1)
    values[4::2] = - values[4::2]

    tests = [("mean", "NaN", lambda z, v, p: np.mean(np.hstack((z, v)))),
             ("mean", 0, lambda z, v, p: np.mean(v)),
             ("median", "NaN", lambda z, v, p: np.median(np.hstack((z, v)))),
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


def test_imputation_most_frequent():
    """Test imputation using the most-frequent strategy."""
    X = np.array([
        [-1, -1,  0,  5],
        [-1,  2, -1,  3],
        [-1,  1,  3, -1],
        [-1,  2,  3,  7],
    ])

    X_true = np.array([
        [2,  0,  5],
        [2,  3,  3],
        [1,  3,  3],
        [2,  3,  7],
    ])

    # scipy.stats.mode, used in Imputer, doesn't return the first most
    # frequent as promised in the doc but the lowest most frequent. When this
    # test will fail after an update of scipy, Imputer will need to be updated
    # to be consistent with the new (correct) behaviour
    _check_statistics(X, X_true, "most_frequent", [np.nan, 2, 3, 3], -1)


def test_imputation_pipeline_grid_search():
    """Test imputation within a pipeline + gridsearch."""
    pipeline = Pipeline([('imputer', Imputer(missing_values=0)),
                         ('tree', tree.DecisionTreeRegressor(random_state=0))])

    parameters = {
        'imputer__strategy': ["mean", "median", "most_frequent"],
        'imputer__axis': [0, 1]
    }

    l = 100
    X = sparse_random_matrix(l, l, density=0.10)
    Y = sparse_random_matrix(l, 1, density=0.10).todense()
    gs = grid_search.GridSearchCV(pipeline, parameters)
    gs.fit(X, Y)


def test_imputation_pickle():
    """Test for pickling imputers."""
    import pickle

    l = 100
    X = sparse_random_matrix(l, l, density=0.10)

    for strategy in ["mean", "median", "most_frequent"]:
        imputer = Imputer(missing_values=0, strategy=strategy)
        imputer.fit(X)

        imputer_pickled = pickle.loads(pickle.dumps(imputer))

        assert_array_equal(imputer.transform(X.copy()),
                           imputer_pickled.transform(X.copy()),
                           "Fail to transform the data after pickling "
                           "(strategy = %s)" % (strategy))


def test_imputation_copy():
    """Test imputation with copy=True."""
    l = 5

    # Test default behaviour and with copy=True
    for params in [{}, {'copy': True}]:
        X = sparse_random_matrix(l, l, density=0.75, random_state=0)

        # Dense
        imputer = Imputer(missing_values=0, strategy="mean", **params)
        Xt = imputer.fit(X).transform(X)
        Xt[0, 0] = np.nan
        # Check that the objects are different and that they don't use
        # the same buffer
        assert_false(np.all(X.todense() == Xt))

        # Sparse
        imputer = Imputer(missing_values=0, strategy="mean", **params)
        X = X.todense()
        Xt = imputer.fit(X).transform(X)
        Xt[0, 0] = np.nan
        # Check that the objects are different and that they don't use
        # the same buffer
        assert_false(np.all(X == Xt))
