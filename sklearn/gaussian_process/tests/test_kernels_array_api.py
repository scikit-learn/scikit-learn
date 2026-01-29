"""Array API tests for Gaussian process kernels."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from sklearn import config_context
from sklearn.gaussian_process.kernels import (
    RBF,
    ConstantKernel,
    DotProduct,
    Exponentiation,
    ExpSineSquared,
    Matern,
    Product,
    RationalQuadratic,
    Sum,
    WhiteKernel,
)
from sklearn.utils._array_api import (
    _atol_for_type,
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
def test_rbf_kernel_array_api(array_namespace, device_, dtype_name):
    """Test RBF kernel works with array API inputs and returns correct namespace."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(10, 3).astype(dtype_name)
    Y_np = rng.randn(5, 3).astype(dtype_name)

    X_xp = xp.asarray(X_np, device=device_)
    Y_xp = xp.asarray(Y_np, device=device_)

    kernel = RBF(length_scale=1.0)

    # Test k(X, X) without Y
    K_np = kernel(X_np)
    with config_context(array_api_dispatch=True):
        K_xp = kernel(X_xp)

        # Check output is in correct namespace
        result_xp, _ = get_namespace(K_xp)
        input_xp, _ = get_namespace(X_xp)
        assert result_xp.__name__ == input_xp.__name__

    # Convert back to numpy for comparison
    K_xp_np = _convert_to_numpy(K_xp, xp)
    rtol = 1e-5 if dtype_name == "float64" else 1e-4
    assert_allclose(K_xp_np, K_np, rtol=rtol)

    # Test k(X, Y) with Y provided
    K_np_xy = kernel(X_np, Y_np)
    with config_context(array_api_dispatch=True):
        K_xp_xy = kernel(X_xp, Y_xp)

    K_xp_xy_np = _convert_to_numpy(K_xp_xy, xp)
    assert_allclose(K_xp_xy_np, K_np_xy, rtol=rtol)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_rbf_kernel_gradient_array_api(array_namespace, device_, dtype_name):
    """Test RBF kernel gradient computation works with array API."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(10, 3).astype(dtype_name)
    X_xp = xp.asarray(X_np, device=device_)

    kernel = RBF(length_scale=1.0)

    # Test gradient computation
    K_np, K_grad_np = kernel(X_np, eval_gradient=True)
    with config_context(array_api_dispatch=True):
        K_xp, K_grad_xp = kernel(X_xp, eval_gradient=True)

    K_xp_np = _convert_to_numpy(K_xp, xp)
    K_grad_xp_np = _convert_to_numpy(K_grad_xp, xp)

    rtol = 1e-5 if dtype_name == "float64" else 1e-4
    assert_allclose(K_xp_np, K_np, rtol=rtol)
    assert_allclose(K_grad_xp_np, K_grad_np, rtol=rtol, atol=_atol_for_type(dtype_name))


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_rbf_kernel_anisotropic_array_api(array_namespace, device_, dtype_name):
    """Test RBF kernel with anisotropic length scales and array API."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(10, 3).astype(dtype_name)
    X_xp = xp.asarray(X_np, device=device_)

    # Anisotropic kernel with different length scales per dimension
    kernel = RBF(length_scale=[1.0, 2.0, 0.5])

    K_np, K_grad_np = kernel(X_np, eval_gradient=True)
    with config_context(array_api_dispatch=True):
        K_xp, K_grad_xp = kernel(X_xp, eval_gradient=True)

    K_xp_np = _convert_to_numpy(K_xp, xp)
    K_grad_xp_np = _convert_to_numpy(K_grad_xp, xp)

    rtol = 1e-5 if dtype_name == "float64" else 1e-4
    assert_allclose(K_xp_np, K_np, rtol=rtol)
    assert_allclose(K_grad_xp_np, K_grad_np, rtol=rtol, atol=_atol_for_type(dtype_name))


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_rbf_kernel_diag_array_api(array_namespace, device_, dtype_name):
    """Test RBF kernel diag method works with array API."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(10, 3).astype(dtype_name)
    X_xp = xp.asarray(X_np, device=device_)

    kernel = RBF(length_scale=1.0)

    # RBF inherits diag from NormalizedKernelMixin which returns ones
    diag_np = kernel.diag(X_np)
    with config_context(array_api_dispatch=True):
        diag_xp = kernel.diag(X_xp)

    diag_xp_np = _convert_to_numpy(diag_xp, xp)
    assert_allclose(diag_xp_np, diag_np)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_constant_kernel_array_api(array_namespace, device_, dtype_name):
    """Test ConstantKernel works with array API inputs."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(10, 3).astype(dtype_name)
    Y_np = rng.randn(5, 3).astype(dtype_name)

    X_xp = xp.asarray(X_np, device=device_)
    Y_xp = xp.asarray(Y_np, device=device_)

    kernel = ConstantKernel(constant_value=2.0)

    # Test k(X, X)
    K_np = kernel(X_np)
    with config_context(array_api_dispatch=True):
        K_xp = kernel(X_xp)

    K_xp_np = _convert_to_numpy(K_xp, xp)
    rtol = 1e-5 if dtype_name == "float64" else 1e-4
    assert_allclose(K_xp_np, K_np, rtol=rtol)

    # Test k(X, Y)
    K_np_xy = kernel(X_np, Y_np)
    with config_context(array_api_dispatch=True):
        K_xp_xy = kernel(X_xp, Y_xp)

    K_xp_xy_np = _convert_to_numpy(K_xp_xy, xp)
    assert_allclose(K_xp_xy_np, K_np_xy, rtol=rtol)

    # Test diag
    diag_np = kernel.diag(X_np)
    with config_context(array_api_dispatch=True):
        diag_xp = kernel.diag(X_xp)

    diag_xp_np = _convert_to_numpy(diag_xp, xp)
    assert_allclose(diag_xp_np, diag_np, rtol=rtol)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_white_kernel_array_api(array_namespace, device_, dtype_name):
    """Test WhiteKernel works with array API inputs."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(10, 3).astype(dtype_name)
    Y_np = rng.randn(5, 3).astype(dtype_name)

    X_xp = xp.asarray(X_np, device=device_)
    Y_xp = xp.asarray(Y_np, device=device_)

    kernel = WhiteKernel(noise_level=0.5)

    # Test k(X, X) - returns diagonal matrix
    K_np = kernel(X_np)
    with config_context(array_api_dispatch=True):
        K_xp = kernel(X_xp)

    K_xp_np = _convert_to_numpy(K_xp, xp)
    rtol = 1e-5 if dtype_name == "float64" else 1e-4
    assert_allclose(K_xp_np, K_np, rtol=rtol)

    # Test k(X, Y) - returns zeros
    K_np_xy = kernel(X_np, Y_np)
    with config_context(array_api_dispatch=True):
        K_xp_xy = kernel(X_xp, Y_xp)

    K_xp_xy_np = _convert_to_numpy(K_xp_xy, xp)
    assert_allclose(K_xp_xy_np, K_np_xy, rtol=rtol)

    # Test diag
    diag_np = kernel.diag(X_np)
    with config_context(array_api_dispatch=True):
        diag_xp = kernel.diag(X_xp)

    diag_xp_np = _convert_to_numpy(diag_xp, xp)
    assert_allclose(diag_xp_np, diag_np, rtol=rtol)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize("nu", [0.5, 1.5, 2.5, np.inf, 3.0])
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_matern_kernel_array_api(nu, array_namespace, device_, dtype_name):
    """Test Matern kernel works with array API inputs for different nu values."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(10, 3).astype(dtype_name)
    Y_np = rng.randn(5, 3).astype(dtype_name)

    X_xp = xp.asarray(X_np, device=device_)
    Y_xp = xp.asarray(Y_np, device=device_)

    kernel = Matern(length_scale=1.0, nu=nu)

    # Test k(X, X)
    K_np = kernel(X_np)
    with config_context(array_api_dispatch=True):
        K_xp = kernel(X_xp)

    K_xp_np = _convert_to_numpy(K_xp, xp)
    # Higher tolerance for general nu case which uses Bessel functions
    rtol = 1e-4 if dtype_name == "float64" else 1e-3
    assert_allclose(K_xp_np, K_np, rtol=rtol)

    # Test k(X, Y)
    K_np_xy = kernel(X_np, Y_np)
    with config_context(array_api_dispatch=True):
        K_xp_xy = kernel(X_xp, Y_xp)

    K_xp_xy_np = _convert_to_numpy(K_xp_xy, xp)
    assert_allclose(K_xp_xy_np, K_np_xy, rtol=rtol)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_rational_quadratic_kernel_array_api(array_namespace, device_, dtype_name):
    """Test RationalQuadratic kernel works with array API inputs."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(10, 3).astype(dtype_name)
    Y_np = rng.randn(5, 3).astype(dtype_name)

    X_xp = xp.asarray(X_np, device=device_)
    Y_xp = xp.asarray(Y_np, device=device_)

    kernel = RationalQuadratic(length_scale=1.0, alpha=1.5)

    # Test k(X, X)
    K_np = kernel(X_np)
    with config_context(array_api_dispatch=True):
        K_xp = kernel(X_xp)

    K_xp_np = _convert_to_numpy(K_xp, xp)
    rtol = 1e-5 if dtype_name == "float64" else 1e-4
    assert_allclose(K_xp_np, K_np, rtol=rtol)

    # Test k(X, Y)
    K_np_xy = kernel(X_np, Y_np)
    with config_context(array_api_dispatch=True):
        K_xp_xy = kernel(X_xp, Y_xp)

    K_xp_xy_np = _convert_to_numpy(K_xp_xy, xp)
    assert_allclose(K_xp_xy_np, K_np_xy, rtol=rtol)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_rational_quadratic_gradient_array_api(array_namespace, device_, dtype_name):
    """Test RationalQuadratic kernel gradient works with array API inputs."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(10, 3).astype(dtype_name)
    X_xp = xp.asarray(X_np, device=device_)

    kernel = RationalQuadratic(length_scale=1.0, alpha=1.5)

    K_np, K_grad_np = kernel(X_np, eval_gradient=True)
    with config_context(array_api_dispatch=True):
        K_xp, K_grad_xp = kernel(X_xp, eval_gradient=True)

    K_xp_np = _convert_to_numpy(K_xp, xp)
    K_grad_xp_np = _convert_to_numpy(K_grad_xp, xp)

    rtol = 1e-5 if dtype_name == "float64" else 2e-4
    assert_allclose(K_xp_np, K_np, rtol=rtol)
    assert_allclose(K_grad_xp_np, K_grad_np, rtol=rtol, atol=_atol_for_type(dtype_name))


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_exp_sine_squared_kernel_array_api(array_namespace, device_, dtype_name):
    """Test ExpSineSquared kernel works with array API inputs."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(10, 3).astype(dtype_name)
    Y_np = rng.randn(5, 3).astype(dtype_name)

    X_xp = xp.asarray(X_np, device=device_)
    Y_xp = xp.asarray(Y_np, device=device_)

    kernel = ExpSineSquared(length_scale=1.0, periodicity=2.0)

    # Test k(X, X)
    K_np = kernel(X_np)
    with config_context(array_api_dispatch=True):
        K_xp = kernel(X_xp)

    K_xp_np = _convert_to_numpy(K_xp, xp)
    rtol = 1e-5 if dtype_name == "float64" else 1e-4
    assert_allclose(K_xp_np, K_np, rtol=rtol)

    # Test k(X, Y)
    K_np_xy = kernel(X_np, Y_np)
    with config_context(array_api_dispatch=True):
        K_xp_xy = kernel(X_xp, Y_xp)

    K_xp_xy_np = _convert_to_numpy(K_xp_xy, xp)
    assert_allclose(K_xp_xy_np, K_np_xy, rtol=rtol)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_exp_sine_squared_gradient_array_api(array_namespace, device_, dtype_name):
    """Test ExpSineSquared kernel gradient works with array API inputs."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(10, 3).astype(dtype_name)
    X_xp = xp.asarray(X_np, device=device_)

    kernel = ExpSineSquared(length_scale=1.0, periodicity=2.0)

    K_np, K_grad_np = kernel(X_np, eval_gradient=True)
    with config_context(array_api_dispatch=True):
        K_xp, K_grad_xp = kernel(X_xp, eval_gradient=True)

    K_xp_np = _convert_to_numpy(K_xp, xp)
    K_grad_xp_np = _convert_to_numpy(K_grad_xp, xp)

    rtol = 1e-5 if dtype_name == "float64" else 1e-4
    assert_allclose(K_xp_np, K_np, rtol=rtol)
    assert_allclose(K_grad_xp_np, K_grad_np, rtol=rtol, atol=_atol_for_type(dtype_name))


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_dot_product_kernel_array_api(array_namespace, device_, dtype_name):
    """Test DotProduct kernel works with array API inputs."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(10, 3).astype(dtype_name)
    Y_np = rng.randn(5, 3).astype(dtype_name)

    X_xp = xp.asarray(X_np, device=device_)
    Y_xp = xp.asarray(Y_np, device=device_)

    kernel = DotProduct(sigma_0=1.0)

    # Test k(X, X)
    K_np = kernel(X_np)
    with config_context(array_api_dispatch=True):
        K_xp = kernel(X_xp)

    K_xp_np = _convert_to_numpy(K_xp, xp)
    rtol = 1e-5 if dtype_name == "float64" else 1e-4
    assert_allclose(K_xp_np, K_np, rtol=rtol)

    # Test k(X, Y)
    K_np_xy = kernel(X_np, Y_np)
    with config_context(array_api_dispatch=True):
        K_xp_xy = kernel(X_xp, Y_xp)

    K_xp_xy_np = _convert_to_numpy(K_xp_xy, xp)
    assert_allclose(K_xp_xy_np, K_np_xy, rtol=rtol)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_dot_product_diag_array_api(array_namespace, device_, dtype_name):
    """Test DotProduct kernel diag method works with array API inputs."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(10, 3).astype(dtype_name)
    X_xp = xp.asarray(X_np, device=device_)

    kernel = DotProduct(sigma_0=1.0)

    diag_np = kernel.diag(X_np)
    with config_context(array_api_dispatch=True):
        diag_xp = kernel.diag(X_xp)

    diag_xp_np = _convert_to_numpy(diag_xp, xp)
    rtol = 1e-5 if dtype_name == "float64" else 1e-4
    assert_allclose(diag_xp_np, diag_np, rtol=rtol)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_dot_product_gradient_array_api(array_namespace, device_, dtype_name):
    """Test DotProduct kernel gradient works with array API inputs."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(10, 3).astype(dtype_name)
    X_xp = xp.asarray(X_np, device=device_)

    kernel = DotProduct(sigma_0=1.0)

    K_np, K_grad_np = kernel(X_np, eval_gradient=True)
    with config_context(array_api_dispatch=True):
        K_xp, K_grad_xp = kernel(X_xp, eval_gradient=True)

    K_xp_np = _convert_to_numpy(K_xp, xp)
    K_grad_xp_np = _convert_to_numpy(K_grad_xp, xp)

    rtol = 1e-5 if dtype_name == "float64" else 1e-4
    assert_allclose(K_xp_np, K_np, rtol=rtol)
    assert_allclose(K_grad_xp_np, K_grad_np, rtol=rtol, atol=_atol_for_type(dtype_name))


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_sum_kernel_array_api(array_namespace, device_, dtype_name):
    """Test Sum kernel works with array API inputs."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(10, 3).astype(dtype_name)
    Y_np = rng.randn(5, 3).astype(dtype_name)

    X_xp = xp.asarray(X_np, device=device_)
    Y_xp = xp.asarray(Y_np, device=device_)

    kernel = Sum(RBF(length_scale=1.0), ConstantKernel(constant_value=2.0))

    # Test k(X, X)
    K_np = kernel(X_np)
    with config_context(array_api_dispatch=True):
        K_xp = kernel(X_xp)

    K_xp_np = _convert_to_numpy(K_xp, xp)
    rtol = 1e-5 if dtype_name == "float64" else 1e-4
    assert_allclose(K_xp_np, K_np, rtol=rtol)

    # Test k(X, Y)
    K_np_xy = kernel(X_np, Y_np)
    with config_context(array_api_dispatch=True):
        K_xp_xy = kernel(X_xp, Y_xp)

    K_xp_xy_np = _convert_to_numpy(K_xp_xy, xp)
    assert_allclose(K_xp_xy_np, K_np_xy, rtol=rtol)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_sum_kernel_gradient_array_api(array_namespace, device_, dtype_name):
    """Test Sum kernel gradient works with array API inputs."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(10, 3).astype(dtype_name)
    X_xp = xp.asarray(X_np, device=device_)

    kernel = Sum(RBF(length_scale=1.0), ConstantKernel(constant_value=2.0))

    K_np, K_grad_np = kernel(X_np, eval_gradient=True)
    with config_context(array_api_dispatch=True):
        K_xp, K_grad_xp = kernel(X_xp, eval_gradient=True)

    K_xp_np = _convert_to_numpy(K_xp, xp)
    K_grad_xp_np = _convert_to_numpy(K_grad_xp, xp)

    rtol = 1e-5 if dtype_name == "float64" else 1e-4
    assert_allclose(K_xp_np, K_np, rtol=rtol)
    assert_allclose(K_grad_xp_np, K_grad_np, rtol=rtol, atol=_atol_for_type(dtype_name))


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_product_kernel_array_api(array_namespace, device_, dtype_name):
    """Test Product kernel works with array API inputs."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(10, 3).astype(dtype_name)
    Y_np = rng.randn(5, 3).astype(dtype_name)

    X_xp = xp.asarray(X_np, device=device_)
    Y_xp = xp.asarray(Y_np, device=device_)

    kernel = Product(RBF(length_scale=1.0), ConstantKernel(constant_value=2.0))

    # Test k(X, X)
    K_np = kernel(X_np)
    with config_context(array_api_dispatch=True):
        K_xp = kernel(X_xp)

    K_xp_np = _convert_to_numpy(K_xp, xp)
    rtol = 1e-5 if dtype_name == "float64" else 1e-4
    assert_allclose(K_xp_np, K_np, rtol=rtol)

    # Test k(X, Y)
    K_np_xy = kernel(X_np, Y_np)
    with config_context(array_api_dispatch=True):
        K_xp_xy = kernel(X_xp, Y_xp)

    K_xp_xy_np = _convert_to_numpy(K_xp_xy, xp)
    assert_allclose(K_xp_xy_np, K_np_xy, rtol=rtol)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_product_kernel_gradient_array_api(array_namespace, device_, dtype_name):
    """Test Product kernel gradient works with array API inputs."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(10, 3).astype(dtype_name)
    X_xp = xp.asarray(X_np, device=device_)

    kernel = Product(RBF(length_scale=1.0), ConstantKernel(constant_value=2.0))

    K_np, K_grad_np = kernel(X_np, eval_gradient=True)
    with config_context(array_api_dispatch=True):
        K_xp, K_grad_xp = kernel(X_xp, eval_gradient=True)

    K_xp_np = _convert_to_numpy(K_xp, xp)
    K_grad_xp_np = _convert_to_numpy(K_grad_xp, xp)

    rtol = 1e-5 if dtype_name == "float64" else 1e-4
    assert_allclose(K_xp_np, K_np, rtol=rtol)
    assert_allclose(K_grad_xp_np, K_grad_np, rtol=rtol, atol=_atol_for_type(dtype_name))


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_exponentiation_kernel_array_api(array_namespace, device_, dtype_name):
    """Test Exponentiation kernel works with array API inputs."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(10, 3).astype(dtype_name)
    Y_np = rng.randn(5, 3).astype(dtype_name)

    X_xp = xp.asarray(X_np, device=device_)
    Y_xp = xp.asarray(Y_np, device=device_)

    kernel = Exponentiation(RBF(length_scale=1.0), exponent=2)

    # Test k(X, X)
    K_np = kernel(X_np)
    with config_context(array_api_dispatch=True):
        K_xp = kernel(X_xp)

    K_xp_np = _convert_to_numpy(K_xp, xp)
    rtol = 1e-5 if dtype_name == "float64" else 1e-4
    assert_allclose(K_xp_np, K_np, rtol=rtol)

    # Test k(X, Y)
    K_np_xy = kernel(X_np, Y_np)
    with config_context(array_api_dispatch=True):
        K_xp_xy = kernel(X_xp, Y_xp)

    K_xp_xy_np = _convert_to_numpy(K_xp_xy, xp)
    assert_allclose(K_xp_xy_np, K_np_xy, rtol=rtol)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_exponentiation_kernel_gradient_array_api(array_namespace, device_, dtype_name):
    """Test Exponentiation kernel gradient works with array API inputs."""
    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(10, 3).astype(dtype_name)
    X_xp = xp.asarray(X_np, device=device_)

    kernel = Exponentiation(RBF(length_scale=1.0), exponent=2)

    K_np, K_grad_np = kernel(X_np, eval_gradient=True)
    with config_context(array_api_dispatch=True):
        K_xp, K_grad_xp = kernel(X_xp, eval_gradient=True)

    K_xp_np = _convert_to_numpy(K_xp, xp)
    K_grad_xp_np = _convert_to_numpy(K_grad_xp, xp)

    rtol = 1e-5 if dtype_name == "float64" else 1e-4
    assert_allclose(K_xp_np, K_np, rtol=rtol)
    assert_allclose(K_grad_xp_np, K_grad_np, rtol=rtol, atol=_atol_for_type(dtype_name))
