import numpy as np

from sklearn.base import clone

from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises_regex
from sklearn.utils.testing import assert_allclose

from sklearn.preprocessing import TransformTargetRegressor
from sklearn.preprocessing import MaxAbsScaler

from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler

from sklearn import datasets

friedman = datasets.make_friedman1(random_state=0)


def test_transform_target_regressor_error_kwargs():
    X = friedman[0]
    y = friedman[1]
    # provide a transformer and functions at the same time
    regr = TransformTargetRegressor(regressor=LinearRegression(),
                                    transformer=StandardScaler(),
                                    func=np.exp, inverse_func=np.log)
    assert_raises_regex(ValueError, "Both 'transformer' and functions"
                        " 'func'/'inverse_func' cannot be set at the"
                        " same time.", regr.fit, X, y)


def test_transform_target_regressor_invertible():
    X = friedman[0]
    y = friedman[1]
    regr = TransformTargetRegressor(regressor=LinearRegression(),
                                    func=np.sqrt, inverse_func=np.log,
                                    check_inverse=True)
    assert_raises_regex(ValueError, "The provided functions or transformer"
                        " are not strictly inverse of each other.",
                        regr.fit, X, y)
    regr = TransformTargetRegressor(regressor=LinearRegression(),
                                    func=np.sqrt, inverse_func=np.log,
                                    check_inverse=False)
    # the transformer/functions are not checked to be invertible the fitting
    # should pass
    regr.fit(X, y)


def _check_standard_scaler(y, y_pred):
    y_mean = np.mean(y, axis=0)
    y_std = np.std(y, axis=0)
    assert_allclose((y - y_mean) / y_std, y_pred)


def _check_max_abs_scaler(y, y_pred):
    max_abs = np.abs(y).max(axis=0)
    assert_allclose(y / max_abs, y_pred)


def test_transform_target_regressor_friedman():
    X = friedman[0]
    y = friedman[1]
    # pass some functions
    regr = TransformTargetRegressor(regressor=LinearRegression(),
                                    func=np.log, inverse_func=np.exp)
    y_pred = regr.fit(X, y).predict(X)
    y_tran = np.ravel(regr.transformer_.transform(y))
    assert_allclose(np.log(y), y_tran)
    assert_allclose(y, np.ravel(regr.transformer_.inverse_transform(
        y_tran.reshape(-1, 1))))
    assert_equal(y.shape, y_pred.shape)
    assert_allclose(y_pred, regr.inverse_func(regr.regressor_.predict(X)))
    lr = LinearRegression().fit(X, regr.func(y))
    assert_allclose(y_pred, regr.inverse_func(lr.predict(X)))
    assert_array_equal(regr.regressor_.coef_.ravel(),
                       lr.coef_.ravel())
    # StandardScaler support 1d array while MaxAbsScaler support only 2d array
    for transformer in (StandardScaler(), MaxAbsScaler()):
        regr = TransformTargetRegressor(regressor=LinearRegression(),
                                        transformer=transformer)
        y_pred = regr.fit(X, y).predict(X)
        assert_equal(y.shape, y_pred.shape)
        y_tran = np.ravel(regr.transformer_.transform(y.reshape(-1, 1)))
        if issubclass(StandardScaler, transformer.__class__):
            _check_standard_scaler(y, y_tran)
        else:
            _check_max_abs_scaler(y, y_tran)
        assert_equal(y.shape, y_pred.shape)
        assert_allclose(y, regr.transformer_.inverse_transform(
                y_tran.reshape(-1, 1)).squeeze())
        lr = LinearRegression()
        transformer2 = clone(transformer)
        lr.fit(X, transformer2.fit_transform(y.reshape(-1, 1)).squeeze())
        assert_allclose(y_pred, transformer2.inverse_transform(
            lr.predict(X).reshape(-1, 1)).squeeze())
        assert_array_equal(regr.regressor_.coef_.squeeze(),
                           lr.coef_.squeeze())


def test_transform_target_regressor_multioutput():
    X = friedman[0]
    y = friedman[1]
    y = np.vstack((y, y ** 2 + 1)).T
    # pass some functions
    regr = TransformTargetRegressor(regressor=LinearRegression(),
                                    func=np.log, inverse_func=np.exp)
    y_pred = regr.fit(X, y).predict(X)
    y_tran = regr.transformer_.transform(y)
    assert_allclose(np.log(y), y_tran)
    assert_allclose(y, regr.transformer_.inverse_transform(y_tran))
    assert_equal(y.shape, y_pred.shape)
    assert_allclose(y_pred, regr.inverse_func(regr.regressor_.predict(X)))
    lr = LinearRegression().fit(X, regr.func(y))
    assert_allclose(y_pred, regr.inverse_func(lr.predict(X)))
    assert_array_equal(regr.regressor_.coef_.ravel(),
                       lr.coef_.ravel())
    # StandardScaler support 1d array while MaxAbsScaler support only 2d array
    for transformer in (StandardScaler(), MaxAbsScaler()):
        regr = TransformTargetRegressor(regressor=LinearRegression(),
                                        transformer=transformer)
        y_pred = regr.fit(X, y).predict(X)
        assert_equal(y.shape, y_pred.shape)
        y_tran = regr.transformer_.transform(y)
        if issubclass(StandardScaler, transformer.__class__):
            _check_standard_scaler(y, y_tran)
        else:
            _check_max_abs_scaler(y, y_tran)
        assert_equal(y.shape, y_pred.shape)
        assert_allclose(y, regr.transformer_.inverse_transform(y_tran))
        transformer2 = clone(transformer)
        lr.fit(X, transformer2.fit_transform(y))
        assert_allclose(y_pred, transformer2.inverse_transform(lr.predict(X)))
        assert_array_equal(regr.regressor_.coef_.squeeze(),
                           lr.coef_.squeeze())


def test_transform_target_regressor_identity():
    X = friedman[0]
    y = friedman[1]
    regr = TransformTargetRegressor()
    y_pred = regr.fit(X, y).predict(X)
    y_pred_2 = LinearRegression().fit(X, y).predict(X)
    assert_array_equal(y_pred, y_pred_2)
