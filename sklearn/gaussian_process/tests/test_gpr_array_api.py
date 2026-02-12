"""Array API tests for GaussianProcessRegressor."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from sklearn import config_context
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel
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
from sklearn.utils.estimator_checks import check_array_api_input_and_values


@pytest.fixture
def gpr_data(request):
    """Generate train/test data for GPR array API tests.

    Returns a dict with numpy arrays and their array API equivalents.
    """
    array_namespace = request.node.funcargs["array_namespace"]
    device_ = request.node.funcargs["device_"]
    dtype_name = request.node.funcargs["dtype_name"]

    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(20, 3).astype(dtype_name)
    y_np = rng.randn(20).astype(dtype_name)
    X_test_np = rng.randn(5, 3).astype(dtype_name)

    return {
        "xp": xp,
        "device": device_,
        "dtype_name": dtype_name,
        "X_np": X_np,
        "y_np": y_np,
        "X_test_np": X_test_np,
        "X_xp": xp.asarray(X_np, device=device_),
        "y_xp": xp.asarray(y_np, device=device_),
        "X_test_xp": xp.asarray(X_test_np, device=device_),
    }


def _check_array_api_result(result_xp, expected_np, data, rtol=None):
    """Check result is in correct namespace and matches expected values."""
    result_ns, _ = get_namespace(result_xp)
    input_ns, _ = get_namespace(data["X_xp"])
    assert result_ns.__name__ == input_ns.__name__

    result_np = _convert_to_numpy(result_xp, data["xp"])
    if rtol is None:
        rtol = 1e-4 if data["dtype_name"] == "float32" else 1e-6
    assert_allclose(result_np, expected_np, rtol=rtol)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_gpr_array_api_compliance(array_namespace, device_, dtype_name):
    """Test GPR compliance with array API using standard checks."""
    # GPR's intermediate computations (Cholesky decomposition, solving linear
    # systems) can have larger numerical differences between backends for
    # float32, so we use a looser rtol for float32.
    rtol = 1e-2 if dtype_name == "float32" else None
    check_array_api_input_and_values(
        "GaussianProcessRegressor",
        GaussianProcessRegressor(kernel=RBF(length_scale=0.9), random_state=42),
        array_namespace=array_namespace,
        device=device_,
        dtype_name=dtype_name,
        rtol=rtol,
    )


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
@pytest.mark.parametrize("return_type", ["std", "cov"])
def test_gpr_predict_uncertainty_array_api(
    array_namespace, device_, dtype_name, return_type, gpr_data
):
    """Test GPR predict with return_std/return_cov works with array API inputs."""
    predict_kwargs = (
        {"return_std": True} if return_type == "std" else {"return_cov": True}
    )
    kernel = RBF(length_scale=1.0)

    # Fit with numpy
    gpr_np = GaussianProcessRegressor(kernel=kernel, random_state=42)
    gpr_np.fit(gpr_data["X_np"], gpr_data["y_np"])
    y_pred_np, uncertainty_np = gpr_np.predict(gpr_data["X_test_np"], **predict_kwargs)

    # Fit with array API
    gpr_xp = GaussianProcessRegressor(kernel=kernel, random_state=42)
    with config_context(array_api_dispatch=True):
        gpr_xp.fit(gpr_data["X_xp"], gpr_data["y_xp"])
        y_pred_xp, uncertainty_xp = gpr_xp.predict(
            gpr_data["X_test_xp"], **predict_kwargs
        )

    _check_array_api_result(y_pred_xp, y_pred_np, gpr_data)
    _check_array_api_result(uncertainty_xp, uncertainty_np, gpr_data)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_gpr_sample_y_array_api(array_namespace, device_, dtype_name, gpr_data):
    """Test GPR sample_y works with array API inputs."""
    kernel = RBF(length_scale=1.0)

    gpr_xp = GaussianProcessRegressor(kernel=kernel, random_state=42)
    with config_context(array_api_dispatch=True):
        gpr_xp.fit(gpr_data["X_xp"], gpr_data["y_xp"])
        y_samples = gpr_xp.sample_y(gpr_data["X_test_xp"], n_samples=3, random_state=42)

        # Check output is in correct namespace
        result_xp, _ = get_namespace(y_samples)
        input_xp, _ = get_namespace(gpr_data["X_xp"])
        assert result_xp.__name__ == input_xp.__name__

        # Check output shape
        assert y_samples.shape == (5, 3)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
@pytest.mark.parametrize("eval_gradient", [False, True])
def test_gpr_log_marginal_likelihood_array_api(
    array_namespace, device_, dtype_name, eval_gradient, gpr_data
):
    """Test GPR log_marginal_likelihood works with array API inputs."""
    kernel = RBF(length_scale=1.0)
    theta = np.array([0.5]) if eval_gradient else None

    # Fit with numpy
    gpr_np = GaussianProcessRegressor(kernel=kernel, random_state=42, optimizer=None)
    gpr_np.fit(gpr_data["X_np"], gpr_data["y_np"])
    result_np = gpr_np.log_marginal_likelihood(theta, eval_gradient=eval_gradient)

    # Fit with array API
    gpr_xp = GaussianProcessRegressor(kernel=kernel, random_state=42, optimizer=None)
    with config_context(array_api_dispatch=True):
        gpr_xp.fit(gpr_data["X_xp"], gpr_data["y_xp"])
        result_xp = gpr_xp.log_marginal_likelihood(theta, eval_gradient=eval_gradient)

    # Compare results
    rtol = 1e-3 if gpr_data["dtype_name"] == "float32" else 1e-6
    if eval_gradient:
        lml_np, grad_np = result_np
        lml_xp, grad_xp = result_xp
        assert_allclose(float(lml_xp), float(lml_np), rtol=rtol)
        grad_xp_np = _convert_to_numpy(grad_xp, gpr_data["xp"])
        assert_allclose(grad_xp_np, grad_np, rtol=rtol)
    else:
        assert_allclose(float(result_xp), float(result_np), rtol=rtol)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_gpr_with_optimizer_array_api(array_namespace, device_, dtype_name, gpr_data):
    """Test GPR with optimizer works with array API inputs."""
    kernel = ConstantKernel(1.0) * RBF(length_scale=1.0)

    # Fit with numpy (with optimization)
    gpr_np = GaussianProcessRegressor(
        kernel=kernel, random_state=42, n_restarts_optimizer=2
    )
    gpr_np.fit(gpr_data["X_np"], gpr_data["y_np"])
    y_pred_np = gpr_np.predict(gpr_data["X_test_np"])

    # Fit with array API (with optimization)
    gpr_xp = GaussianProcessRegressor(
        kernel=kernel, random_state=42, n_restarts_optimizer=2
    )
    with config_context(array_api_dispatch=True):
        gpr_xp.fit(gpr_data["X_xp"], gpr_data["y_xp"])
        y_pred_xp = gpr_xp.predict(gpr_data["X_test_xp"])

    # Compare results - use larger tolerance due to optimizer behavior
    rtol = 1e-3 if gpr_data["dtype_name"] == "float32" else 1e-5
    _check_array_api_result(y_pred_xp, y_pred_np, gpr_data, rtol=rtol)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_gpr_normalize_y_array_api(array_namespace, device_, dtype_name, gpr_data):
    """Test GPR with normalize_y=True works with array API inputs."""
    # Use large y values to test normalization
    rng = np.random.RandomState(42)
    y_np = (rng.randn(20) * 100 + 500).astype(dtype_name)
    y_xp = gpr_data["xp"].asarray(y_np, device=gpr_data["device"])

    kernel = RBF(length_scale=1.0)

    # Fit with numpy
    gpr_np = GaussianProcessRegressor(
        kernel=kernel, random_state=42, normalize_y=True, optimizer=None
    )
    gpr_np.fit(gpr_data["X_np"], y_np)
    y_pred_np = gpr_np.predict(gpr_data["X_test_np"])

    # Fit with array API
    gpr_xp = GaussianProcessRegressor(
        kernel=kernel, random_state=42, normalize_y=True, optimizer=None
    )
    with config_context(array_api_dispatch=True):
        gpr_xp.fit(gpr_data["X_xp"], y_xp)
        y_pred_xp = gpr_xp.predict(gpr_data["X_test_xp"])

    _check_array_api_result(y_pred_xp, y_pred_np, gpr_data)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
@pytest.mark.parametrize("n_targets", [1, 3])
@pytest.mark.parametrize("return_type", ["mean", "std", "cov"])
def test_gpr_multi_output_array_api(
    array_namespace, device_, dtype_name, n_targets, return_type
):
    """Test GPR with multi-output y works with array API inputs."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(20, 3).astype(dtype_name)
    y_np = rng.randn(20, n_targets).astype(dtype_name)
    X_test_np = rng.randn(5, 3).astype(dtype_name)

    kernel = RBF(length_scale=1.0)

    # Fit with numpy
    gpr_np = GaussianProcessRegressor(kernel=kernel, random_state=42, optimizer=None)
    gpr_np.fit(X_np, y_np)

    # Fit with array API
    gpr_xp = GaussianProcessRegressor(kernel=kernel, random_state=42, optimizer=None)
    X_xp = xp.asarray(X_np, device=device_)
    y_xp = xp.asarray(y_np, device=device_)
    X_test_xp = xp.asarray(X_test_np, device=device_)

    with config_context(array_api_dispatch=True):
        gpr_xp.fit(X_xp, y_xp)

        if return_type == "mean":
            y_pred_np = gpr_np.predict(X_test_np)
            y_pred_xp = gpr_xp.predict(X_test_xp)
        elif return_type == "std":
            y_pred_np, std_np = gpr_np.predict(X_test_np, return_std=True)
            y_pred_xp, std_xp = gpr_xp.predict(X_test_xp, return_std=True)
        else:
            y_pred_np, cov_np = gpr_np.predict(X_test_np, return_cov=True)
            y_pred_xp, cov_xp = gpr_xp.predict(X_test_xp, return_cov=True)

    # GPR intermediate computations (Cholesky, triangular solves) accumulate
    # larger numerical differences in float32
    rtol = 1e-3 if dtype_name == "float32" else 1e-6

    y_pred_xp_np = _convert_to_numpy(y_pred_xp, xp)
    assert_allclose(y_pred_xp_np, y_pred_np, rtol=rtol)
    assert y_pred_xp.shape == y_pred_np.shape

    if return_type == "std":
        std_xp_np = _convert_to_numpy(std_xp, xp)
        assert_allclose(std_xp_np, std_np, rtol=rtol)
        assert std_xp.shape == std_np.shape
    elif return_type == "cov":
        cov_xp_np = _convert_to_numpy(cov_xp, xp)
        assert_allclose(cov_xp_np, cov_np, rtol=rtol)
        assert cov_xp.shape == cov_np.shape


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
@pytest.mark.parametrize("return_type", ["std", "cov"])
def test_gpr_normalize_y_uncertainty_array_api(
    array_namespace, device_, dtype_name, return_type, gpr_data
):
    """Test GPR with normalize_y=True and return_std/return_cov."""
    rng = np.random.RandomState(42)
    y_np = (rng.randn(20) * 100 + 500).astype(dtype_name)
    y_xp = gpr_data["xp"].asarray(y_np, device=gpr_data["device"])

    kernel = RBF(length_scale=1.0)

    predict_kwargs = (
        {"return_std": True} if return_type == "std" else {"return_cov": True}
    )

    # Fit with numpy
    gpr_np = GaussianProcessRegressor(
        kernel=kernel, random_state=42, normalize_y=True, optimizer=None
    )
    gpr_np.fit(gpr_data["X_np"], y_np)
    y_pred_np, uncertainty_np = gpr_np.predict(gpr_data["X_test_np"], **predict_kwargs)

    # Fit with array API
    gpr_xp = GaussianProcessRegressor(
        kernel=kernel, random_state=42, normalize_y=True, optimizer=None
    )
    with config_context(array_api_dispatch=True):
        gpr_xp.fit(gpr_data["X_xp"], y_xp)
        y_pred_xp, uncertainty_xp = gpr_xp.predict(
            gpr_data["X_test_xp"], **predict_kwargs
        )

    # GPR covariance with large y-values and normalize_y accumulates more
    # numerical error on GPU backends (MPS), so use a looser tolerance
    rtol = 1e-3 if gpr_data["dtype_name"] == "float32" else None
    _check_array_api_result(y_pred_xp, y_pred_np, gpr_data, rtol=rtol)
    _check_array_api_result(uncertainty_xp, uncertainty_np, gpr_data, rtol=rtol)


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
def test_gpr_namespace_mismatch_error(array_namespace, device_, dtype_name, gpr_data):
    """Test GPR raises error when predict namespace doesn't match fit namespace."""
    kernel = RBF(length_scale=1.0)
    gpr = GaussianProcessRegressor(kernel=kernel, random_state=42, optimizer=None)

    with config_context(array_api_dispatch=True):
        gpr.fit(gpr_data["X_xp"], gpr_data["y_xp"])

        # Predict with numpy should raise an error about namespace mismatch
        with pytest.raises(ValueError, match="array namespace"):
            gpr.predict(gpr_data["X_test_np"])
