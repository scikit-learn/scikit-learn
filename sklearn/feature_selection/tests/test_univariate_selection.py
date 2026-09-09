import warnings

import numpy as np
import pytest
from scipy import sparse

from sklearn.feature_selection import f_regression, r_regression


@pytest.mark.parametrize(
    "container", [lambda X: X, sparse.csr_matrix, sparse.csc_matrix]
)
def test_r_regression_constant_features_no_runtime_warning(container):
    X = np.array([[2.0, 1.0], [2.0, 0.0], [2.0, 10.0], [2.0, 4.0]])
    y = np.array([0.0, 1.0, 1.0, 0.0])

    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        correlation = r_regression(container(X), y)

    assert correlation[0] == 0.0
    assert correlation[1] == pytest.approx(0.32075015)


@pytest.mark.parametrize(
    "container", [lambda X: X, sparse.csr_matrix, sparse.csc_matrix]
)
def test_f_regression_constant_features_no_runtime_warning(container):
    X = np.array([[2.0, 1.0], [2.0, 0.0], [2.0, 10.0], [2.0, 4.0]])
    y = np.array([0.0, 1.0, 1.0, 0.0])

    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        f_statistic, p_values = f_regression(container(X), y)

    assert f_statistic[0] == 0.0
    assert p_values[0] == 1.0
    assert f_statistic[1] == pytest.approx(0.2293578)


@pytest.mark.parametrize(
    "container", [lambda X: X, sparse.csr_matrix, sparse.csc_matrix]
)
def test_r_regression_constant_features_force_finite_false(container):
    X = np.array([[2.0, 1.0], [2.0, 0.0], [2.0, 10.0], [2.0, 4.0]])
    y = np.array([0.0, 1.0, 1.0, 0.0])

    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        correlation = r_regression(container(X), y, force_finite=False)

    assert np.isnan(correlation[0])
    assert correlation[1] == pytest.approx(0.32075015)


def test_r_regression_cancellation_is_clamped():
    X = np.column_stack(
        [
            # The perturbation is below the ULP at this magnitude, so the
            # represented values are constant while the moment subtraction
            # still suffers cancellation.
            np.full(20, 1e12) + np.arange(20) * 1e-8,
            np.arange(20, dtype=float),
        ]
    )
    y = np.arange(20, dtype=float)

    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        correlation = r_regression(X, y)

    assert correlation[0] == 0.0
    assert correlation[1] == 1.0


def test_r_regression_active_features_are_unchanged():
    X = np.array([[2.0, 1.0], [2.0, 0.0], [2.0, 10.0], [2.0, 4.0]])
    y = np.array([0.0, 1.0, 1.0, 0.0])
    expected = r_regression(X[:, 1:], y)[0]

    actual = r_regression(X, y)[1]

    assert actual == expected
