"""
Regression tests for Issue #34003:
  Array API: Ridge with alpha as an array fails due to misconfigured checks

Two root causes:
  1. _BaseRidge._parameter_constraints uses np.ndarray instead of "array-like",
     rejecting non-NumPy array-like objects (e.g. torch.Tensor) at _validate_params.
  2. _ridge_regression's scalar detection uses isinstance(alpha, type(xp.asarray([0.0]))),
     which is broken when xp is derived from X (e.g. torch) but alpha is np.ndarray.
"""

import numpy as np
import pytest

from sklearn import config_context
from sklearn.datasets import make_regression
from sklearn.linear_model import Ridge
from sklearn.utils._array_api import yield_namespace_device_dtype_combinations

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_multi_target_data(n_samples=200, n_targets=5, random_state=0):
    X, y = make_regression(
        n_samples=n_samples,
        n_features=10,
        n_targets=n_targets,
        noise=0.1,
        random_state=random_state,
    )
    return X, y


# ---------------------------------------------------------------------------
# Tests for numpy-only (no array_api_dispatch) - should be unaffected
# ---------------------------------------------------------------------------


def test_ridge_alpha_scalar_numpy():
    """Scalar float alpha works with numpy (baseline, must not regress)."""
    X, y = _make_multi_target_data()
    reg = Ridge(alpha=1.0, solver="svd")
    reg.fit(X, y)
    assert reg.coef_.shape == (5, 10)


def test_ridge_alpha_numpy_array_numpy_dispatch():
    """np.ndarray alpha works in pure numpy mode (baseline, must not regress)."""
    X, y = _make_multi_target_data(n_targets=5)
    reg = Ridge(alpha=np.array([0.1, 0.5, 1.0, 2.0, 5.0]), solver="svd")
    reg.fit(X, y)
    assert reg.coef_.shape == (5, 10)


def test_ridge_alpha_list_numpy_dispatch():
    """List alpha is accepted (converted internally)."""
    X, y = _make_multi_target_data(n_targets=3)
    reg = Ridge(alpha=[0.1, 1.0, 5.0], solver="svd")
    reg.fit(X, y)
    assert reg.coef_.shape == (3, 10)


# ---------------------------------------------------------------------------
# Tests for Array API dispatch with numpy.array alpha (Bug 1)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "array_namespace, device, dtype_name",
    yield_namespace_device_dtype_combinations(),
)
def test_ridge_array_api_numpy_array_alpha(array_namespace, device, dtype_name):
    """
    Regression test for Bug 1:
      Ridge(alpha=np.ndarray).fit(xp_tensor, xp_tensor) must succeed.

    Previously failed at _ridge_regression line 703:
      isinstance(np.array([...]), type(xp.asarray([0.0])))
      = isinstance(np.ndarray, torch.Tensor)  [with torch namespace]
      = False  ->  check_scalar was wrongly called -> TypeError
    """
    from sklearn.utils._array_api import _array_api_for_tests

    xp = _array_api_for_tests(array_namespace, device)
    n_targets = 4

    X_np, y_np = _make_multi_target_data(n_targets=n_targets)
    X_np = X_np.astype(dtype_name)
    y_np = y_np.astype(dtype_name)

    X_xp = xp.asarray(X_np, device=device)
    y_xp = xp.asarray(y_np, device=device)

    # alpha is np.ndarray, X and y are backend tensors
    alpha = np.array([0.1, 0.5, 1.0, 2.0], dtype=np.float64)

    reg = Ridge(alpha=alpha, solver="svd")
    with config_context(array_api_dispatch=True):
        reg.fit(X_xp, y_xp)

    assert reg.coef_.shape == (n_targets, X_np.shape[1])


@pytest.mark.parametrize(
    "array_namespace, device, dtype_name",
    yield_namespace_device_dtype_combinations(),
)
def test_ridge_array_api_xp_array_alpha(array_namespace, device, dtype_name):
    """
    Regression test for Bug 2:
      Ridge(alpha=xp.tensor([...])).fit(xp_tensor, xp_tensor) must succeed.

    Previously failed at _validate_params:
      _parameter_constraints has np.ndarray as the only array type;
      torch.Tensor is not np.ndarray -> InvalidParameterError.
    """
    from sklearn.utils._array_api import _array_api_for_tests

    xp = _array_api_for_tests(array_namespace, device)
    n_targets = 4

    X_np, y_np = _make_multi_target_data(n_targets=n_targets)
    X_np = X_np.astype(dtype_name)
    y_np = y_np.astype(dtype_name)

    X_xp = xp.asarray(X_np, device=device)
    y_xp = xp.asarray(y_np, device=device)

    # alpha is an xp-native tensor
    alpha_np = np.array([0.1, 0.5, 1.0, 2.0], dtype=np.float64)
    alpha_xp = xp.asarray(alpha_np, device=device)

    reg = Ridge(alpha=alpha_xp, solver="svd")
    with config_context(array_api_dispatch=True):
        reg.fit(X_xp, y_xp)

    assert reg.coef_.shape == (n_targets, X_np.shape[1])


@pytest.mark.parametrize(
    "array_namespace, device, dtype_name",
    yield_namespace_device_dtype_combinations(),
)
def test_ridge_array_api_scalar_alpha(array_namespace, device, dtype_name):
    """
    Scalar alpha must continue to work with Array API dispatch (no regression).
    """
    from sklearn.utils._array_api import _array_api_for_tests

    xp = _array_api_for_tests(array_namespace, device)
    n_targets = 3

    X_np, y_np = _make_multi_target_data(n_targets=n_targets)
    X_np = X_np.astype(dtype_name)
    y_np = y_np.astype(dtype_name)

    X_xp = xp.asarray(X_np, device=device)
    y_xp = xp.asarray(y_np, device=device)

    reg = Ridge(alpha=1.0, solver="svd")
    with config_context(array_api_dispatch=True):
        reg.fit(X_xp, y_xp)

    assert reg.coef_.shape == (n_targets, X_np.shape[1])


# ---------------------------------------------------------------------------
# Tests for coefficient correctness (array alpha must match scalar-alpha results)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "array_namespace, device, dtype_name",
    yield_namespace_device_dtype_combinations(),
)
def test_ridge_array_api_alpha_array_coef_correctness(
    array_namespace, device, dtype_name
):
    """
    Uniform array alpha (all same value) must produce the same coefficients
    as the equivalent scalar alpha.
    """
    from sklearn.utils._array_api import _array_api_for_tests, _convert_to_numpy
    from sklearn.utils._testing import assert_allclose

    xp = _array_api_for_tests(array_namespace, device)
    n_targets = 3
    alpha_val = 1.0

    X_np, y_np = _make_multi_target_data(n_targets=n_targets)
    X_np = X_np.astype(dtype_name)
    y_np = y_np.astype(dtype_name)

    X_xp = xp.asarray(X_np, device=device)
    y_xp = xp.asarray(y_np, device=device)

    # Scalar alpha reference
    reg_scalar = Ridge(alpha=alpha_val, solver="svd")
    with config_context(array_api_dispatch=True):
        reg_scalar.fit(X_xp, y_xp)

    # np.ndarray alpha (all same value)
    reg_arr_np = Ridge(alpha=np.full(n_targets, alpha_val), solver="svd")
    with config_context(array_api_dispatch=True):
        reg_arr_np.fit(X_xp, y_xp)

    # xp-native array alpha
    reg_arr_xp = Ridge(
        alpha=xp.asarray(np.full(n_targets, alpha_val), device=device), solver="svd"
    )
    with config_context(array_api_dispatch=True):
        reg_arr_xp.fit(X_xp, y_xp)

    coef_scalar = _convert_to_numpy(reg_scalar.coef_, xp=xp)
    coef_arr_np = _convert_to_numpy(reg_arr_np.coef_, xp=xp)
    coef_arr_xp = _convert_to_numpy(reg_arr_xp.coef_, xp=xp)

    atol = 1e-4 if dtype_name == "float32" else 1e-8
    assert_allclose(coef_arr_np, coef_scalar, atol=atol)
    assert_allclose(coef_arr_xp, coef_scalar, atol=atol)


# ---------------------------------------------------------------------------
# Parameter validation edge cases
# ---------------------------------------------------------------------------


def test_ridge_invalid_alpha_rejected():
    """Negative scalar alpha must still be rejected."""
    from sklearn.utils._param_validation import InvalidParameterError

    with pytest.raises((InvalidParameterError, ValueError)):
        Ridge(alpha=-1.0).fit(*_make_multi_target_data())


def test_ridge_wrong_shape_alpha_rejected():
    """Array alpha with wrong number of elements must raise ValueError."""
    X, y = _make_multi_target_data(n_targets=5)
    # alpha has 3 values but y has 5 targets
    reg = Ridge(alpha=np.array([1.0, 2.0, 3.0]), solver="svd")
    with pytest.raises(ValueError, match="Number of targets"):
        reg.fit(X, y)


def test_ridge_single_element_array_alpha():
    """Single-element array alpha broadcasts to all targets."""
    X, y = _make_multi_target_data(n_targets=5)
    reg = Ridge(alpha=np.array([1.0]), solver="svd")
    reg.fit(X, y)
    assert reg.coef_.shape == (5, 10)
