import numpy as np

from sklearn.base import clone

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises_regex
from sklearn.utils.testing import assert_allclose

from sklearn.preprocessing import FunctionTransformer
from sklearn.preprocessing import TransformTargetRegressor
from sklearn.preprocessing import StandardScaler

from sklearn.linear_model import LinearRegression, LogisticRegression, Lasso

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
    assert_raises_regex(TypeError, "fit\(\) got an unexpected keyword argument"
                        " 'sample_weight'", regr.fit, X, y,
                        sample_weight=sample_weight)
    # provide a classifier instead of a regressor
    regr = TransformTargetRegressor(regressor=LogisticRegression())
    assert_raises_regex(TypeError, "The regressor LogisticRegression is of"
                        " type classifier. Provide a regressor instead.",
                        regr.fit, X, y)


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


def _check_custom_scaler(y, y_pred):
    assert_allclose(y + 1, y_pred)


def test_transform_target_regressor_functions():
    X = friedman[0]
    y = friedman[1]
    regr = TransformTargetRegressor(regressor=LinearRegression(),
                                    func=np.log, inverse_func=np.exp)
    y_pred = regr.fit(X, y).predict(X)
    # check the transformer output
    y_tran = regr.transformer_.transform(y.reshape(-1, 1)).squeeze()
    assert_allclose(np.log(y), y_tran)
    assert_allclose(y, regr.transformer_.inverse_transform(
        y_tran.reshape(-1, 1)).squeeze())
    assert_equal(y.shape, y_pred.shape)
    assert_allclose(y_pred, regr.inverse_func(regr.regressor_.predict(X)))
    # check the regressor output
    lr = LinearRegression().fit(X, regr.func(y))
    assert_allclose(regr.regressor_.coef_.ravel(), lr.coef_.ravel())


def test_transform_target_regressor_functions_multioutput():
    X = friedman[0]
    y = np.vstack((friedman[1], friedman[1] ** 2 + 1)).T
    regr = TransformTargetRegressor(regressor=LinearRegression(),
                                    func=np.log, inverse_func=np.exp)
    y_pred = regr.fit(X, y).predict(X)
    # check the transformer output
    y_tran = regr.transformer_.transform(y)
    assert_allclose(np.log(y), y_tran)
    assert_allclose(y, regr.transformer_.inverse_transform(y_tran))
    assert_equal(y.shape, y_pred.shape)
    assert_allclose(y_pred, regr.inverse_func(regr.regressor_.predict(X)))
    # check the regressor output
    lr = LinearRegression().fit(X, regr.func(y))
    assert_allclose(regr.regressor_.coef_.ravel(), lr.coef_.ravel())


def test_transform_target_regressor_1d_transformer():
    X = friedman[0]
    y = friedman[1]
    transformer = FunctionTransformer(func=lambda x: x + 1,
                                      inverse_func=lambda x: x - 1,
                                      validate=False)
    regr = TransformTargetRegressor(regressor=LinearRegression(),
                                    transformer=transformer)
    y_pred = regr.fit(X, y).predict(X)
    assert_equal(y.shape, y_pred.shape)
    # consistency forward transform
    y_tran = regr.transformer_.transform(y)
    _check_custom_scaler(y, y_tran)
    assert_equal(y.shape, y_pred.shape)
    # consistency inverse transform
    assert_allclose(y, regr.transformer_.inverse_transform(
        y_tran).squeeze())
    # consistency of the regressor
    lr = LinearRegression()
    transformer2 = clone(transformer)
    lr.fit(X, transformer2.fit_transform(y))
    y_lr_pred = lr.predict(X)
    assert_allclose(y_pred, transformer2.inverse_transform(y_lr_pred))
    assert_allclose(regr.regressor_.coef_, lr.coef_)


def test_transform_target_regressor_1d_transformer_multioutput():
    X = friedman[0]
    y = np.vstack((friedman[1], friedman[1] ** 2 + 1)).T
    transformer = FunctionTransformer(func=lambda x: x + 1,
                                      inverse_func=lambda x: x - 1,
                                      validate=False)
    regr = TransformTargetRegressor(regressor=LinearRegression(),
                                    transformer=transformer)
    y_pred = regr.fit(X, y).predict(X)
    assert_equal(y.shape, y_pred.shape)
    # consistency forward transform
    y_tran = regr.transformer_.transform(y)
    _check_custom_scaler(y, y_tran)
    assert_equal(y.shape, y_pred.shape)
    # consistency inverse transform
    assert_allclose(y, regr.transformer_.inverse_transform(
        y_tran).squeeze())
    # consistency of the regressor
    lr = LinearRegression()
    transformer2 = clone(transformer)
    lr.fit(X, transformer2.fit_transform(y))
    y_lr_pred = lr.predict(X)
    assert_allclose(y_pred, transformer2.inverse_transform(y_lr_pred))
    assert_allclose(regr.regressor_.coef_, lr.coef_)


def test_transform_target_regressor_2d_transformer():
    X = friedman[0]
    y = friedman[1]
    transformer = StandardScaler()
    regr = TransformTargetRegressor(regressor=LinearRegression(),
                                    transformer=transformer)
    y_pred = regr.fit(X, y).predict(X)
    assert_equal(y.shape, y_pred.shape)
    # consistency forward transform
    y_tran = regr.transformer_.transform(y.reshape(-1, 1)).squeeze()
    _check_standard_scaler(y, y_tran)
    assert_equal(y.shape, y_pred.shape)
    # consistency inverse transform
    assert_allclose(y, regr.transformer_.inverse_transform(
        y_tran).squeeze())
    # consistency of the regressor
    lr = LinearRegression()
    transformer2 = clone(transformer)
    lr.fit(X, transformer2.fit_transform(y.reshape(-1, 1)).squeeze())
    y_lr_pred = lr.predict(X)
    assert_allclose(y_pred, transformer2.inverse_transform(y_lr_pred))
    assert_allclose(regr.regressor_.coef_, lr.coef_)


def test_transform_target_regressor_2d_transformer_multioutput():
    X = friedman[0]
    y = np.vstack((friedman[1], friedman[1] ** 2 + 1)).T
    transformer = StandardScaler()
    regr = TransformTargetRegressor(regressor=LinearRegression(),
                                    transformer=transformer)
    y_pred = regr.fit(X, y).predict(X)
    assert_equal(y.shape, y_pred.shape)
    # consistency forward transform
    y_tran = regr.transformer_.transform(y)
    _check_standard_scaler(y, y_tran)
    assert_equal(y.shape, y_pred.shape)
    # consistency inverse transform
    assert_allclose(y, regr.transformer_.inverse_transform(
        y_tran).squeeze())
    # consistency of the regressor
    lr = LinearRegression()
    transformer2 = clone(transformer)
    lr.fit(X, transformer2.fit_transform(y))
    y_lr_pred = lr.predict(X)
    assert_allclose(y_pred, transformer2.inverse_transform(y_lr_pred))
    assert_allclose(regr.regressor_.coef_, lr.coef_)


def test_transform_target_regressor_single_to_multi():
    X = friedman[0]
    y = friedman[1] + 1j * (friedman[1] ** 2 + 1)

    def func(y):
        """Stack the real and imaginary part."""
        return np.hstack((np.real(y), np.imag(y)))

    tt = TransformTargetRegressor(func=func, check_inverse=False)
    tt.fit(X, y)
    y_pred = tt.predict(X)
    assert_equal(y_pred.shape, (100, 2))


def test_transform_target_regressor_multi_to_single():
    X = friedman[0]
    y = np.transpose([friedman[1], (friedman[1] ** 2 + 1)])

    def func(y):
        return np.sqrt(y[:, 0] ** 2 + y[:, 1] ** 2)

    tt = TransformTargetRegressor(func=func, check_inverse=False)
    tt.fit(X, y)
    y_pred = tt.predict(X)
    assert_equal(y_pred.shape, (100, 1))
