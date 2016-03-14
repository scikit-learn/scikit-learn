import numpy as np
import scipy.sparse as sp

from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raises_regex

from sklearn import datasets
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.linear_model import Lasso
from sklearn.multioutput import MultiOutputRegressor


def test_multi_target_regression():
    X, y = datasets.make_regression(n_targets=3)
    X_train, y_train = X[:50], y[:50]
    X_test, y_test = X[50:], y[50:]

    references = np.zeros_like(y_test)
    for n in range(3):
        rgr = GradientBoostingRegressor(random_state=0)
        rgr.fit(X_train, y_train[:, n])
        references[:,n] = rgr.predict(X_test)

    rgr = MultiOutputRegressor(GradientBoostingRegressor(random_state=0))
    rgr.fit(X_train, y_train)
    y_pred = rgr.predict(X_test)

    assert_almost_equal(references, y_pred)


def test_multi_target_regression_one_target():
    # Test multi target regression raises
    X, y = datasets.make_regression(n_targets=1)
    X_train, y_train = X[:50], y[:50]
    X_test, y_test = X[50:], y[50:]

    rgr = MultiOutputRegressor(GradientBoostingRegressor(random_state=0))
    assert_raises(ValueError, rgr.fit, X_train, y_train)


def test_multi_target_sparse_regression():
    X, y = datasets.make_regression(n_targets=3)
    X_train, y_train = X[:50], y[:50]
    X_test, y_test = X[50:], y[50:]

    for sparse in [sp.csr_matrix, sp.csc_matrix, sp.coo_matrix, sp.dok_matrix,
                   sp.lil_matrix]:
        rgr = MultiOutputRegressor(Lasso(random_state=0))
        rgr_sparse = MultiOutputRegressor(Lasso(random_state=0))

        rgr.fit(X_train, y_train)
        rgr_sparse.fit(sparse(X_train), y_train)

        assert_almost_equal(rgr.predict(X_test), rgr_sparse.predict(sparse(X_test)))


def test_multi_target_sample_weights_api():
    X = [[1,2,3], [4,5,6]]
    y = [[3.141, 2.718], [2.718, 3.141]]
    w = [0.8, 0.6]

    rgr = MultiOutputRegressor(Lasso())
    assert_raises_regex(ValueError, "does not support sample weights",
                        rgr.fit, X, y, w)

    # no exception should be raised if the base estimator supports weights
    rgr = MultiOutputRegressor(GradientBoostingRegressor(random_state=0))
    rgr.fit(X, y, w)


def test_multi_target_sample_weights():
    # weighted regressor
    Xw = [[1,2,3], [4,5,6]]
    yw = [[3.141, 2.718], [2.718, 3.141]]
    w = [2., 1.]
    rgr_w = MultiOutputRegressor(GradientBoostingRegressor(random_state=0))
    rgr_w.fit(Xw, yw, w)

    # unweighted, but with repeated samples
    X = [[1,2,3], [1,2,3], [4,5,6]]
    y = [[3.141, 2.718], [3.141, 2.718], [2.718, 3.141]]
    rgr = MultiOutputRegressor(GradientBoostingRegressor(random_state=0))
    rgr.fit(X, y)

    X_test = [[1.5,2.5,3.5], [3.5,4.5,5.5]]
    assert_almost_equal(rgr.predict(X_test), rgr_w.predict(X_test))
