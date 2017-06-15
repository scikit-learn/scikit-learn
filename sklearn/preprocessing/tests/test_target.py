import numpy as np

from sklearn.base import clone

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises_regex
from sklearn.utils.testing import assert_allclose

from sklearn.preprocessing import TransformTargetRegressor
from sklearn.preprocessing import MaxAbsScaler

from sklearn.linear_model import LinearRegression, Lasso
from sklearn.preprocessing import StandardScaler

from sklearn import datasets

friedman = datasets.make_friedman1(random_state=0)


def test_transform_target_regressor_error():
    X = friedman[0]
    y = friedman[1]
    # provide a transformer and functions at the same time
    regr = TransformTargetRegressor(regressor=LinearRegression(),
                                    transformer=StandardScaler(),
                                    func=np.exp, inverse_func=np.log)
    assert_raises_regex(ValueError, "'transformer' and functions"
                        " 'func'/'inverse_func' cannot both be set.",
                        regr.fit, X, y)
    # fit with sample_weight with a regressor which does not support it
    sample_weight = np.ones((y.shape[0],))
    regr = TransformTargetRegressor(regressor=Lasso(),
                                    transformer=StandardScaler())
    assert_raises_regex(ValueError, "The regressor Lasso does not support"
                        " sample weight.", regr.fit, X, y,
                        sample_weight=sample_weight)


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
    # create a multioutput y
    # keep why to be 2d and it will squeezed when relevant
    Y = [y.reshape(-1, 1), np.vstack((y, y ** 2 + 1)).T]
    for y_2d in Y:
        regr = TransformTargetRegressor(regressor=LinearRegression(),
                                        func=np.log, inverse_func=np.exp)
        y_pred = regr.fit(X, y_2d.squeeze()).predict(X)
        y_tran = regr.transformer_.transform(y_2d).squeeze()
        assert_allclose(np.log(y_2d.squeeze()), y_tran)
        assert_allclose(y_2d.squeeze(), regr.transformer_.inverse_transform(
            y_tran).squeeze())
        assert_equal(y_2d.squeeze().shape, y_pred.shape)
        assert_allclose(y_pred, regr.inverse_func(regr.regressor_.predict(X)))
        lr = LinearRegression().fit(X, regr.func(y_2d.squeeze()))
        assert_allclose(y_pred, regr.inverse_func(lr.predict(X)))
        assert_allclose(regr.regressor_.coef_.ravel(), lr.coef_.ravel())
        # StandardScaler support 1d array while MaxAbsScaler support only 2d
        # array
        for transformer in (StandardScaler(), MaxAbsScaler()):
            regr = TransformTargetRegressor(regressor=LinearRegression(),
                                            transformer=transformer)
            y_pred = regr.fit(X, y_2d.squeeze()).predict(X)
            assert_equal(y_2d.squeeze().shape, y_pred.shape)
            y_tran = regr.transformer_.transform(y_2d).squeeze()
            if issubclass(StandardScaler, transformer.__class__):
                _check_standard_scaler(y_2d.squeeze(), y_tran)
            else:
                _check_max_abs_scaler(y_2d.squeeze(), y_tran)
            assert_equal(y_2d.squeeze().shape, y_pred.shape)
            if y_tran.ndim == 1:
                y_tran = y_tran.reshape(-1, 1)
            assert_allclose(y_2d.squeeze(),
                            regr.transformer_.inverse_transform(
                                y_tran).squeeze())
            lr = LinearRegression()
            transformer2 = clone(transformer)
            lr.fit(X, transformer2.fit_transform(y_2d).squeeze())
            y_lr_pred = lr.predict(X)
            if y_lr_pred.ndim == 1:
                y_lr_pred = y_lr_pred.reshape(-1, 1)
            assert_allclose(y_pred, transformer2.inverse_transform(
                y_lr_pred).squeeze())
            assert_allclose(regr.regressor_.coef_.squeeze(),
                            lr.coef_.squeeze())
