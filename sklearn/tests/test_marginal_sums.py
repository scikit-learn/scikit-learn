import numpy as np
import pytest
import scipy
from sklearn.marginal_sums import MarginalSumsRegression
from sklearn.utils._testing import assert_array_almost_equal

weights = np.array([300, 700, 600, 200])
y_sums = np.array([66000, 231000, 120000, 60000])
features = np.array(
    [
        [1.0, 0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0, 0.0],
        [0.0, 1.0, 0.0, 1.0],
        [0.0, 1.0, 1.0, 0.0],
    ]
)
factors = np.array([1.04907574, 0.95345681, 1.18709246, 0.79148547])
y_pred = [220.03697532, 330.01772326, 199.98151234, 299.9379686]


def test_with_given_weights():
    y = y_sums.reshape(-1, 1)
    X = np.insert(features, 0, weights, axis=1)

    msr = MarginalSumsRegression()
    msr.fit(X, y)

    assert_array_almost_equal(msr.factors, factors)
    assert_array_almost_equal(msr.predict(features), y_pred)
    assert_array_almost_equal(msr.fit_predict(X, y), y_pred)


def test_missing_weights():
    y = y_sums.reshape(-1, 1)
    X = features

    msr = MarginalSumsRegression()
    msg = r"0 detected in first column. Expected weights > 0."
    with pytest.raises(ValueError, match=msg):
        msr.fit(X, y)


def test_add_weights():
    y = np.repeat(y_sums / weights, weights, axis=0).reshape(-1, 1)
    X = np.repeat(features, weights, axis=0)

    msr = MarginalSumsRegression(add_weights=True)
    msr.fit(X, y)

    assert_array_almost_equal(msr.factors, factors)
    assert_array_almost_equal(msr.predict(features), y_pred)


def test_convergence_warning():
    y = y_sums.reshape(-1, 1)
    X = np.insert(features, 0, weights, axis=1)

    msr = MarginalSumsRegression(max_iter=3)

    msg = r"not converge"
    with pytest.warns(UserWarning, match=msg):
        msr.fit(X, y)


def test_sparse_input():
    y = y_sums.reshape(-1, 1)
    X = scipy.sparse.csr_matrix(np.insert(features, 0, weights, axis=1))

    msr = MarginalSumsRegression()
    msr.fit(X, y)

    assert_array_almost_equal(msr.factors, factors)
    assert_array_almost_equal(msr.predict(features), y_pred)
    assert_array_almost_equal(msr.fit_predict(X, y), y_pred)


def test_not_onehot_encoded_input():
    not_encoded_features = np.array(
        [
            [2, 0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0, 0.0],
            [3, 0.0, 0.0, 1.0],
            [1, 0.0, 1.0, 0.0],
        ]
    )

    y = y_sums.reshape(-1, 1)
    X = np.insert(not_encoded_features, 0, weights, axis=1)

    msr = MarginalSumsRegression()
    msg = r"Value different from 1 or 0 detected. Only onehot encoded values expected."
    with pytest.raises(ValueError, match=msg):
        msr.fit(X, y)
