"""
Module from skgarden
"""

import numpy as np
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal

from sklearn.datasets import load_boston
from sklearn.model_selection import train_test_split

from sklearn.ensemble._qrf import ExtraTreesQuantileRegressor
from sklearn.ensemble._qrf import RandomForestQuantileRegressor
from sklearn.ensemble._forest import RandomForestRegressor

boston = load_boston()
X, y = boston.data, boston.target
X_train, X_test, y_train, y_test = train_test_split(
    X, y, train_size=0.6, test_size=0.4, random_state=0)
X_train = np.array(X_train, dtype=np.float32)
X_test = np.array(X_test, dtype=np.float32)
estimators = [
    RandomForestQuantileRegressor(random_state=0),
    ExtraTreesQuantileRegressor(random_state=0)]
approx_estimators = [
    RandomForestQuantileRegressor(random_state=0, method='sample', n_estimators=250),
    ExtraTreesQuantileRegressor(random_state=0, method='sample', n_estimators=250)
]


def test_quantile_attributes():
    for est in estimators:
        est.fit(X_train, y_train)

        # If a sample is not present in a particular tree, that
        # corresponding leaf is marked as -1.
        assert_array_equal(
            np.vstack(np.where(est.y_train_leaves_ == -1)),
            np.vstack(np.where(est.y_weights_ == 0))
        )

        # Should sum up to number of leaf nodes.
        assert_array_equal(
            np.sum(est.y_weights_, axis=1),
            [sum(tree.tree_.children_left == -1) for tree in est.estimators_]
        )

        n_est = est.n_estimators
        est.set_params(bootstrap=False)
        est.fit(X_train, y_train)
        assert_array_almost_equal(
            np.sum(est.y_weights_, axis=1),
            [sum(tree.tree_.children_left == -1) for tree in est.estimators_],
            6
        )
        assert np.all(est.y_train_leaves_ != -1)


def test_random_sample_RF_difference():
    # The QRF with method='sample' only operates on different values stored in the tree_.value array
    # So when calling model.apply the results should be the same, but the indexed values should differ
    qrf = RandomForestQuantileRegressor(random_state=0, n_estimators=10, method='sample', max_depth=2)
    rf = RandomForestRegressor(random_state=0, n_estimators=10, max_depth=2)
    qrf.fit(X_train, y_train)
    rf.fit(X_train, y_train)

    # indices from apply should be the same
    assert_array_equal(qrf.apply(X_test), rf.apply(X_test))

    # the result from indexing into tree_.value array with these indices should be different
    assert not np.array_equal(qrf.estimators_[0].tree_.value[qrf.estimators_[0].apply(X_test)],
                              rf.estimators_[0].tree_.value[rf.estimators_[0].apply(X_test)])


def test_max_depth_None_rfqr():
    # Since each leaf is pure and has just one unique value.
    # the median equals any quantile.
    rng = np.random.RandomState(0)
    X = rng.randn(10, 1)
    y = np.linspace(0.0, 100.0, 10)

    rfqr_estimators = [
        RandomForestQuantileRegressor(random_state=0, bootstrap=False, max_depth=None),
        RandomForestQuantileRegressor(random_state=0, bootstrap=False, max_depth=None, method='sample')
    ]
    for rfqr in rfqr_estimators:
        rfqr.fit(X, y)
        rfqr.quantiles = 0.5
        a = rfqr.predict(X)
        for quantile in (0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9, 1):
            rfqr.quantiles = quantile
            assert_array_almost_equal(
                a, rfqr.predict(X), 5)


def test_forest_toy_data():
    rng = np.random.RandomState(105)
    x1 = rng.randn(1, 10)
    X1 = np.tile(x1, (10000, 1))
    x2 = 20.0 * rng.randn(1, 10)
    X2 = np.tile(x2, (10000, 1))
    X = np.vstack((X1, X2))

    y1 = rng.randn(10000)
    y2 = 5.0 + rng.randn(10000)
    y = np.concatenate((y1, y2))

    for est in estimators:
        est.set_params(max_depth=1)
        est.fit(X, y)
        for quantile in (0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9, 1):
            est.quantiles = quantile
            assert_array_almost_equal(
                est.predict(x1),
                [np.quantile(y1, quantile)], 3)
            assert_array_almost_equal(
                est.predict(x2),
                [np.quantile(y2, quantile)], 3)

    # the approximate methods have a lower precision, which is to be expected
    for est in approx_estimators:
        est.set_params(max_depth=1)
        est.fit(X, y)
        for quantile in (0.2, 0.3, 0.5, 0.7):
            est.quantiles = quantile
            assert_array_almost_equal(
                est.predict(x1),
                [np.quantile(y1, quantile)], 0)
            assert_array_almost_equal(
                est.predict(x2),
                [np.quantile(y2, quantile)], 0)


if __name__ == "sklearn.ensemble.tests.test_ensemble":
    print("Test ensemble")
    test_quantile_attributes()
    test_random_sample_RF_difference()
    test_max_depth_None_rfqr()
    test_forest_toy_data()
