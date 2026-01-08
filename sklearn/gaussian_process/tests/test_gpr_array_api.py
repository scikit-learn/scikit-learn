"""Array API tests for GaussianProcessRegressor."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from sklearn import config_context
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, Matern
from sklearn.utils._array_api import (
    _convert_to_numpy,
    _get_namespace_device_dtype_ids,
    get_namespace,
    yield_namespace_device_dtype_combinations,
)
from sklearn.utils._testing import (
    _array_api_for_tests,
    skip_if_array_api_compat_not_configured,
)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_gpr_fit_predict_array_api(array_namespace, device_, dtype_name):
    """Test GPR fit and predict work with array API inputs."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(20, 3).astype(dtype_name)
    y_np = rng.randn(20).astype(dtype_name)
    X_test_np = rng.randn(5, 3).astype(dtype_name)

    X_xp = xp.asarray(X_np, device=device_)
    y_xp = xp.asarray(y_np, device=device_)
    X_test_xp = xp.asarray(X_test_np, device=device_)

    kernel = RBF(length_scale=1.0)

    # Fit with numpy
    gpr_np = GaussianProcessRegressor(kernel=kernel, random_state=42)
    gpr_np.fit(X_np, y_np)
    y_pred_np = gpr_np.predict(X_test_np)

    # Fit with array API
    gpr_xp = GaussianProcessRegressor(kernel=kernel, random_state=42)
    with config_context(array_api_dispatch=True):
        gpr_xp.fit(X_xp, y_xp)
        y_pred_xp = gpr_xp.predict(X_test_xp)

        # Check output is in correct namespace
        result_xp, _ = get_namespace(y_pred_xp)
        input_xp, _ = get_namespace(X_xp)
        assert result_xp.__name__ == input_xp.__name__

    # Compare results
    y_pred_xp_np = _convert_to_numpy(y_pred_xp, xp)
    rtol = 1e-4 if dtype_name == "float32" else 1e-6
    assert_allclose(y_pred_xp_np, y_pred_np, rtol=rtol)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_gpr_predict_std_array_api(array_namespace, device_, dtype_name):
    """Test GPR predict with return_std works with array API inputs."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(20, 3).astype(dtype_name)
    y_np = rng.randn(20).astype(dtype_name)
    X_test_np = rng.randn(5, 3).astype(dtype_name)

    X_xp = xp.asarray(X_np, device=device_)
    y_xp = xp.asarray(y_np, device=device_)
    X_test_xp = xp.asarray(X_test_np, device=device_)

    kernel = RBF(length_scale=1.0)

    # Fit with numpy
    gpr_np = GaussianProcessRegressor(kernel=kernel, random_state=42)
    gpr_np.fit(X_np, y_np)
    y_pred_np, y_std_np = gpr_np.predict(X_test_np, return_std=True)

    # Fit with array API
    gpr_xp = GaussianProcessRegressor(kernel=kernel, random_state=42)
    with config_context(array_api_dispatch=True):
        gpr_xp.fit(X_xp, y_xp)
        y_pred_xp, y_std_xp = gpr_xp.predict(X_test_xp, return_std=True)

        # Check outputs are in correct namespace
        result_xp, _ = get_namespace(y_std_xp)
        input_xp, _ = get_namespace(X_xp)
        assert result_xp.__name__ == input_xp.__name__

    # Compare results
    y_pred_xp_np = _convert_to_numpy(y_pred_xp, xp)
    y_std_xp_np = _convert_to_numpy(y_std_xp, xp)
    rtol = 1e-4 if dtype_name == "float32" else 1e-6
    assert_allclose(y_pred_xp_np, y_pred_np, rtol=rtol)
    assert_allclose(y_std_xp_np, y_std_np, rtol=rtol)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_gpr_predict_cov_array_api(array_namespace, device_, dtype_name):
    """Test GPR predict with return_cov works with array API inputs."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(20, 3).astype(dtype_name)
    y_np = rng.randn(20).astype(dtype_name)
    X_test_np = rng.randn(5, 3).astype(dtype_name)

    X_xp = xp.asarray(X_np, device=device_)
    y_xp = xp.asarray(y_np, device=device_)
    X_test_xp = xp.asarray(X_test_np, device=device_)

    kernel = RBF(length_scale=1.0)

    # Fit with numpy
    gpr_np = GaussianProcessRegressor(kernel=kernel, random_state=42)
    gpr_np.fit(X_np, y_np)
    y_pred_np, y_cov_np = gpr_np.predict(X_test_np, return_cov=True)

    # Fit with array API
    gpr_xp = GaussianProcessRegressor(kernel=kernel, random_state=42)
    with config_context(array_api_dispatch=True):
        gpr_xp.fit(X_xp, y_xp)
        y_pred_xp, y_cov_xp = gpr_xp.predict(X_test_xp, return_cov=True)

        # Check outputs are in correct namespace
        result_xp, _ = get_namespace(y_cov_xp)
        input_xp, _ = get_namespace(X_xp)
        assert result_xp.__name__ == input_xp.__name__

    # Compare results
    y_pred_xp_np = _convert_to_numpy(y_pred_xp, xp)
    y_cov_xp_np = _convert_to_numpy(y_cov_xp, xp)
    rtol = 1e-4 if dtype_name == "float32" else 1e-6
    assert_allclose(y_pred_xp_np, y_pred_np, rtol=rtol)
    assert_allclose(y_cov_xp_np, y_cov_np, rtol=rtol)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_gpr_sample_y_array_api(array_namespace, device_, dtype_name):
    """Test GPR sample_y works with array API inputs."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(20, 3).astype(dtype_name)
    y_np = rng.randn(20).astype(dtype_name)
    X_test_np = rng.randn(5, 3).astype(dtype_name)

    X_xp = xp.asarray(X_np, device=device_)
    y_xp = xp.asarray(y_np, device=device_)
    X_test_xp = xp.asarray(X_test_np, device=device_)

    kernel = RBF(length_scale=1.0)

    # Fit with array API
    gpr_xp = GaussianProcessRegressor(kernel=kernel, random_state=42)
    with config_context(array_api_dispatch=True):
        gpr_xp.fit(X_xp, y_xp)
        y_samples = gpr_xp.sample_y(X_test_xp, n_samples=3, random_state=42)

        # Check output is in correct namespace
        result_xp, _ = get_namespace(y_samples)
        input_xp, _ = get_namespace(X_xp)
        assert result_xp.__name__ == input_xp.__name__

        # Check output shape
        assert y_samples.shape == (5, 3)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_gpr_log_marginal_likelihood_array_api(array_namespace, device_, dtype_name):
    """Test GPR log_marginal_likelihood works with array API inputs."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(20, 3).astype(dtype_name)
    y_np = rng.randn(20).astype(dtype_name)

    X_xp = xp.asarray(X_np, device=device_)
    y_xp = xp.asarray(y_np, device=device_)

    kernel = RBF(length_scale=1.0)

    # Fit with numpy
    gpr_np = GaussianProcessRegressor(kernel=kernel, random_state=42, optimizer=None)
    gpr_np.fit(X_np, y_np)
    lml_np = gpr_np.log_marginal_likelihood()

    # Fit with array API
    gpr_xp = GaussianProcessRegressor(kernel=kernel, random_state=42, optimizer=None)
    with config_context(array_api_dispatch=True):
        gpr_xp.fit(X_xp, y_xp)
        lml_xp = gpr_xp.log_marginal_likelihood()

    # Compare results
    rtol = 1e-4 if dtype_name == "float32" else 1e-6
    assert_allclose(float(lml_xp), float(lml_np), rtol=rtol)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_gpr_log_marginal_likelihood_gradient_array_api(
    array_namespace, device_, dtype_name
):
    """Test GPR log_marginal_likelihood with gradient works with array API."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(20, 3).astype(dtype_name)
    y_np = rng.randn(20).astype(dtype_name)

    X_xp = xp.asarray(X_np, device=device_)
    y_xp = xp.asarray(y_np, device=device_)

    kernel = RBF(length_scale=1.0)
    theta = np.array([0.5])

    # Fit with numpy
    gpr_np = GaussianProcessRegressor(kernel=kernel, random_state=42, optimizer=None)
    gpr_np.fit(X_np, y_np)
    lml_np, grad_np = gpr_np.log_marginal_likelihood(theta, eval_gradient=True)

    # Fit with array API
    gpr_xp = GaussianProcessRegressor(kernel=kernel, random_state=42, optimizer=None)
    with config_context(array_api_dispatch=True):
        gpr_xp.fit(X_xp, y_xp)
        lml_xp, grad_xp = gpr_xp.log_marginal_likelihood(theta, eval_gradient=True)

    # Compare results
    rtol = 1e-3 if dtype_name == "float32" else 1e-6
    assert_allclose(float(lml_xp), float(lml_np), rtol=rtol)
    grad_xp_np = _convert_to_numpy(grad_xp, xp)
    assert_allclose(grad_xp_np, grad_np, rtol=rtol)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_gpr_with_optimizer_array_api(array_namespace, device_, dtype_name):
    """Test GPR with optimizer works with array API inputs."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(20, 3).astype(dtype_name)
    y_np = rng.randn(20).astype(dtype_name)
    X_test_np = rng.randn(5, 3).astype(dtype_name)

    X_xp = xp.asarray(X_np, device=device_)
    y_xp = xp.asarray(y_np, device=device_)
    X_test_xp = xp.asarray(X_test_np, device=device_)

    kernel = ConstantKernel(1.0) * RBF(length_scale=1.0)

    # Fit with numpy (with optimization)
    gpr_np = GaussianProcessRegressor(
        kernel=kernel, random_state=42, n_restarts_optimizer=2
    )
    gpr_np.fit(X_np, y_np)
    y_pred_np = gpr_np.predict(X_test_np)

    # Fit with array API (with optimization)
    gpr_xp = GaussianProcessRegressor(
        kernel=kernel, random_state=42, n_restarts_optimizer=2
    )
    with config_context(array_api_dispatch=True):
        gpr_xp.fit(X_xp, y_xp)
        y_pred_xp = gpr_xp.predict(X_test_xp)

    # Compare results - use larger tolerance due to optimizer behavior
    y_pred_xp_np = _convert_to_numpy(y_pred_xp, xp)
    rtol = 1e-3 if dtype_name == "float32" else 1e-5
    assert_allclose(y_pred_xp_np, y_pred_np, rtol=rtol)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_gpr_normalize_y_array_api(array_namespace, device_, dtype_name):
    """Test GPR with normalize_y=True works with array API inputs."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(20, 3).astype(dtype_name)
    y_np = (rng.randn(20) * 100 + 500).astype(dtype_name)  # Large values
    X_test_np = rng.randn(5, 3).astype(dtype_name)

    X_xp = xp.asarray(X_np, device=device_)
    y_xp = xp.asarray(y_np, device=device_)
    X_test_xp = xp.asarray(X_test_np, device=device_)

    kernel = RBF(length_scale=1.0)

    # Fit with numpy
    gpr_np = GaussianProcessRegressor(
        kernel=kernel, random_state=42, normalize_y=True, optimizer=None
    )
    gpr_np.fit(X_np, y_np)
    y_pred_np = gpr_np.predict(X_test_np)

    # Fit with array API
    gpr_xp = GaussianProcessRegressor(
        kernel=kernel, random_state=42, normalize_y=True, optimizer=None
    )
    with config_context(array_api_dispatch=True):
        gpr_xp.fit(X_xp, y_xp)
        y_pred_xp = gpr_xp.predict(X_test_xp)

    # Compare results
    y_pred_xp_np = _convert_to_numpy(y_pred_xp, xp)
    rtol = 1e-4 if dtype_name == "float32" else 1e-6
    assert_allclose(y_pred_xp_np, y_pred_np, rtol=rtol)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_gpr_matern_kernel_array_api(array_namespace, device_, dtype_name):
    """Test GPR with Matern kernel works with array API inputs."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(20, 3).astype(dtype_name)
    y_np = rng.randn(20).astype(dtype_name)
    X_test_np = rng.randn(5, 3).astype(dtype_name)

    X_xp = xp.asarray(X_np, device=device_)
    y_xp = xp.asarray(y_np, device=device_)
    X_test_xp = xp.asarray(X_test_np, device=device_)

    kernel = Matern(length_scale=1.0, nu=1.5)

    # Fit with numpy
    gpr_np = GaussianProcessRegressor(kernel=kernel, random_state=42, optimizer=None)
    gpr_np.fit(X_np, y_np)
    y_pred_np = gpr_np.predict(X_test_np)

    # Fit with array API
    gpr_xp = GaussianProcessRegressor(kernel=kernel, random_state=42, optimizer=None)
    with config_context(array_api_dispatch=True):
        gpr_xp.fit(X_xp, y_xp)
        y_pred_xp = gpr_xp.predict(X_test_xp)

    # Compare results
    y_pred_xp_np = _convert_to_numpy(y_pred_xp, xp)
    rtol = 1e-4 if dtype_name == "float32" else 1e-6
    assert_allclose(y_pred_xp_np, y_pred_np, rtol=rtol)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    [
        combo
        for combo in yield_namespace_device_dtype_combinations()
        if combo[0] != "numpy"  # Skip numpy since numpy->numpy won't raise
    ],
    ids=_get_namespace_device_dtype_ids,
)
def test_gpr_namespace_mismatch_error(array_namespace, device_, dtype_name):
    """Test GPR raises error when predict namespace doesn't match fit namespace."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(20, 3).astype(dtype_name)
    y_np = rng.randn(20).astype(dtype_name)
    X_test_np = rng.randn(5, 3).astype(dtype_name)

    X_xp = xp.asarray(X_np, device=device_)
    y_xp = xp.asarray(y_np, device=device_)

    kernel = RBF(length_scale=1.0)
    gpr = GaussianProcessRegressor(kernel=kernel, random_state=42, optimizer=None)

    with config_context(array_api_dispatch=True):
        gpr.fit(X_xp, y_xp)

        # Predict with numpy should raise an error about namespace mismatch
        with pytest.raises(ValueError, match="array namespace"):
            gpr.predict(X_test_np)
