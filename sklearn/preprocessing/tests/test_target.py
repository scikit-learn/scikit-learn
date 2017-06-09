import numpy as np

from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises_regex
from sklearn.utils.testing import assert_allclose

from sklearn.preprocessing import TransformedTargetRegressor

from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler

from sklearn import datasets

friedman = datasets.make_friedman1(random_state=0)


def test_transformed_target_regressor_error_kwargs():
    X = friedman[0]
    y = friedman[1]
    # provide a transformer and functions at the same time
    regr = TransformedTargetRegressor(regressor=LinearRegression(),
                                      transformer=StandardScaler(),
                                      func=np.exp, inverse_func=np.log)
    assert_raises_regex(ValueError, "Both 'transformer' and functions"
                        " 'func'/'inverse_func' cannot be set at the"
                        " same time.", regr.fit, X, y)


def test_transformed_target_regressor_invertible():
    X = friedman[0]
    y = friedman[1]
    regr = TransformedTargetRegressor(regressor=LinearRegression(),
                                      func=np.sqrt, inverse_func=np.log,
                                      check_inverse=True)
    assert_raises_regex(ValueError, "The provided functions or transformer"
                        " are not strictly inverse of each other.",
                        regr.fit, X, y)
    regr = TransformedTargetRegressor(regressor=LinearRegression(),
                                      func=np.sqrt, inverse_func=np.log,
                                      check_inverse=False)
    # the transformer/functions are not checked to be invertible the fitting
    # should pass
    regr.fit(X, y)


def test_transformed_target_regressor_friedman():
    X = friedman[0]
    y = friedman[1]
    # pass some functions
    regr = TransformedTargetRegressor(regressor=LinearRegression(),
                                      func=np.log, inverse_func=np.exp)
    y_pred = regr.fit(X, y).predict(X)
    y_tran = np.ravel(regr.transformer_.transform(y))
    assert_array_almost_equal(np.log(y), y_tran)
    assert_array_almost_equal(y, np.ravel(regr.transformer_.inverse_transform(
        y_tran.reshape(-1, 1))))
    assert_equal(y.shape, y_pred.shape)
    assert_allclose(y_pred, regr.inverse_func(regr.regressor_.predict(X)))
    lr = LinearRegression().fit(X, regr.func(y))
    assert_allclose(y_pred, regr.inverse_func(lr.predict(X)))
    assert_array_equal(regr.regressor_.coef_.ravel(),
                       lr.coef_.ravel())
    # pass a transformer
    regr = TransformedTargetRegressor(regressor=LinearRegression(),
                                      transformer=StandardScaler())
    y_pred = regr.fit(X, y).predict(X)
    assert_equal(y.shape, y_pred.shape)
    y_mean = np.mean(y)
    y_std = np.std(y)
    y_tran = np.ravel(regr.transformer_.transform(y.reshape(-1, 1)))
    assert_array_almost_equal((y - y_mean) / y_std, y_tran)
    assert_array_almost_equal(y, np.ravel(regr.transformer_.inverse_transform(
        y_tran.reshape(-1, 1))))
    assert_equal(y.shape, y_pred.shape)
    lr = LinearRegression()
    ss = StandardScaler()
    lr.fit(X, ss.fit_transform(y[:, None])[:, 0])
    assert_allclose(y_pred, ss.inverse_transform(lr.predict(X)))
    assert_array_equal(regr.regressor_.coef_.ravel(),
                       lr.coef_.ravel())


def test_transformed_target_regressor_multioutput():
    X = friedman[0]
    y = friedman[1]
    y = np.vstack((y, y ** 2 + 1)).T
    # pass some functions
    regr = TransformedTargetRegressor(regressor=LinearRegression(),
                                      func=np.log, inverse_func=np.exp)
    y_pred = regr.fit(X, y).predict(X)
    y_tran = regr.transformer_.transform(y)
    assert_array_almost_equal(np.log(y), y_tran)
    assert_array_almost_equal(y, regr.transformer_.inverse_transform(y_tran))
    assert_equal(y.shape, y_pred.shape)
    assert_allclose(y_pred, regr.inverse_func(regr.regressor_.predict(X)))
    lr = LinearRegression().fit(X, regr.func(y))
    assert_allclose(y_pred, regr.inverse_func(lr.predict(X)))
    assert_array_equal(regr.regressor_.coef_.ravel(),
                       lr.coef_.ravel())
    # pass a transformer
    regr = TransformedTargetRegressor(regressor=LinearRegression(),
                                      transformer=StandardScaler())
    y_pred = regr.fit(X, y).predict(X)
    assert_equal(y.shape, y_pred.shape)
    y_mean = np.mean(y, axis=0)
    y_std = np.std(y, axis=0)
    y_tran = regr.transformer_.transform(y)
    assert_array_almost_equal((y - y_mean) / y_std, y_tran)
    assert_array_almost_equal(y, regr.transformer_.inverse_transform(
        y_tran))
    assert_equal(y.shape, y_pred.shape)
    ss = StandardScaler()
    lr.fit(X, ss.fit_transform(y))
    assert_allclose(y_pred, ss.inverse_transform(lr.predict(X)))
    assert_array_equal(regr.regressor_.coef_.ravel(),
                       lr.coef_.ravel())
