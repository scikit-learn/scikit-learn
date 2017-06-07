import numpy as np

from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises_regex

from sklearn.preprocessing import TransformedTargetRegressor

from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler

from sklearn import datasets

friedman = datasets.make_friedman1(random_state=0)


def test_transformed_target_regressor_error_kwargs():
    X = friedman[0]
    y = friedman[1]
    # provide a transformer and functions at the same time
    clf = TransformedTargetRegressor(estimator=LinearRegression(),
                                     transformer=StandardScaler(),
                                     func=np.exp, inverse_func=np.log)
    assert_raises_regex(ValueError, "Both 'transformer' and functions"
                        " 'func'/'inverse_func' cannot be set at the"
                        " same time.", clf.fit, X, y)


def test_transformed_target_regressor_invertible():
    X = friedman[0]
    y = friedman[1]
    clf = TransformedTargetRegressor(estimator=LinearRegression(),
                                     func=np.exp, inverse_func=np.exp,
                                     check_invertible=True)
    assert_raises_regex(ValueError, "The provided functions or transformer"
                        " are not strictly inverse of each other.",
                        clf.fit, X, y)
    clf = TransformedTargetRegressor(estimator=LinearRegression(),
                                     func=np.exp, inverse_func=np.exp,
                                     check_invertible=False)
    # the transformer/functions are not checked to be invertible the fitting
    # should pass
    clf.fit(X, y)


def test_transformed_target_regressor_friedman():
    X = friedman[0]
    y = friedman[1]
    # pass some functions
    clf = TransformedTargetRegressor(estimator=LinearRegression(),
                                     func=np.log, inverse_func=np.exp)
    pred = clf.fit(X, y).predict(X)
    y_tran = np.ravel(clf.transformer_.transform(y))
    assert_array_almost_equal(np.log(y), y_tran)
    assert_array_almost_equal(y, np.ravel(clf.transformer_.inverse_transform(
        y_tran.reshape(-1, 1))))
    assert_equal(y.shape, pred.shape)
    # pass a transformer
    clf = TransformedTargetRegressor(estimator=LinearRegression(),
                                     transformer=StandardScaler())
    pred = clf.fit(X, y).predict(X)
    assert_equal(y.shape, pred.shape)
    y_mean = np.mean(y)
    y_std = np.std(y)
    y_tran = np.ravel(clf.transformer_.transform(y.reshape(-1, 1)))
    assert_array_almost_equal((y - y_mean) / y_std, y_tran)
    assert_array_almost_equal(y, np.ravel(clf.transformer_.inverse_transform(
        y_tran.reshape(-1, 1))))
    assert_equal(y.shape, pred.shape)


def test_transformed_target_regressor_multioutput():
    X = friedman[0]
    y = friedman[1]
    y = np.vstack((y, y ** 2 + 1)).T
    # pass some functions
    clf = TransformedTargetRegressor(estimator=LinearRegression(),
                                     func=np.log, inverse_func=np.exp)
    pred = clf.fit(X, y).predict(X)
    y_tran = clf.transformer_.transform(y)
    assert_array_almost_equal(np.log(y), y_tran)
    assert_array_almost_equal(y, clf.transformer_.inverse_transform(y_tran))
    assert_equal(y.shape, pred.shape)
    # pass a transformer
    clf = TransformedTargetRegressor(estimator=LinearRegression(),
                                     transformer=StandardScaler())
    pred = clf.fit(X, y).predict(X)
    assert_equal(y.shape, pred.shape)
    y_mean = np.mean(y, axis=0)
    y_std = np.std(y, axis=0)
    y_tran = clf.transformer_.transform(y)
    assert_array_almost_equal((y - y_mean) / y_std, y_tran)
    assert_array_almost_equal(y, clf.transformer_.inverse_transform(
        y_tran))
    assert_equal(y.shape, pred.shape)
