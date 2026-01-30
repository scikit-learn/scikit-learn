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


@pytest.fixture
def kernel_data(request):
    """Generate train/test data for kernel array API tests.

    Returns a dict with numpy arrays and their array API equivalents.
    """
    array_namespace = request.node.funcargs["array_namespace"]
    device_ = request.node.funcargs["device_"]
    dtype_name = request.node.funcargs["dtype_name"]

    xp = _array_api_for_tests(array_namespace, device_)

    rng = np.random.RandomState(42)
    X_np = rng.randn(10, 3).astype(dtype_name)
    Y_np = rng.randn(5, 3).astype(dtype_name)

    return {
        "xp": xp,
        "device": device_,
        "dtype_name": dtype_name,
        "X_np": X_np,
        "Y_np": Y_np,
        "X_xp": xp.asarray(X_np, device=device_),
        "Y_xp": xp.asarray(Y_np, device=device_),
    }


def _check_kernel_result(K_xp, K_np, data, rtol=None):
    """Compare kernel output with numpy baseline."""
    K_xp_np = _convert_to_numpy(K_xp, data["xp"])
    if rtol is None:
        rtol = 1e-5 if data["dtype_name"] == "float64" else 1e-4
    assert_allclose(K_xp_np, K_np, rtol=rtol)


def _get_kernel_id(kernel):
    """Generate human-readable ID for kernel parametrization."""
    return type(kernel).__name__


# Kernels for basic call tests (k(X) and k(X, Y))
BASIC_KERNELS = [
    RBF(),
    ConstantKernel(constant_value=2.0),
    WhiteKernel(noise_level=0.5),
    RationalQuadratic(alpha=1.5),
    ExpSineSquared(periodicity=2.0),
    DotProduct(),
    Sum(RBF(), ConstantKernel(constant_value=2.0)),
    Product(RBF(), ConstantKernel(constant_value=2.0)),
    Exponentiation(RBF(), exponent=2),
]

# Kernels for gradient tests (eval_gradient=True)
GRADIENT_KERNELS = [
    RBF(),
    RationalQuadratic(alpha=1.5),
    ExpSineSquared(periodicity=2.0),
    DotProduct(),
    Sum(RBF(), ConstantKernel(constant_value=2.0)),
    Product(RBF(), ConstantKernel(constant_value=2.0)),
    Exponentiation(RBF(), exponent=2),
]


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
@pytest.mark.parametrize("kernel", BASIC_KERNELS, ids=_get_kernel_id)
def test_kernel_call_array_api(
    kernel, array_namespace, device_, dtype_name, kernel_data
):
    """Test kernel __call__ works with array API inputs."""
    # Test k(X, X)
    K_np = kernel(kernel_data["X_np"])
    with config_context(array_api_dispatch=True):
        K_xp = kernel(kernel_data["X_xp"])

    _check_kernel_result(K_xp, K_np, kernel_data)

    # Test k(X, Y)
    K_np_xy = kernel(kernel_data["X_np"], kernel_data["Y_np"])
    with config_context(array_api_dispatch=True):
        K_xp_xy = kernel(kernel_data["X_xp"], kernel_data["Y_xp"])

    _check_kernel_result(K_xp_xy, K_np_xy, kernel_data)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
@pytest.mark.parametrize("kernel", GRADIENT_KERNELS, ids=_get_kernel_id)
def test_kernel_gradient_array_api(
    kernel, array_namespace, device_, dtype_name, kernel_data
):
    """Test kernel gradient computation works with array API."""
    K_np, K_grad_np = kernel(kernel_data["X_np"], eval_gradient=True)
    with config_context(array_api_dispatch=True):
        K_xp, K_grad_xp = kernel(kernel_data["X_xp"], eval_gradient=True)

    _check_kernel_result(K_xp, K_np, kernel_data)

    K_grad_xp_np = _convert_to_numpy(K_grad_xp, kernel_data["xp"])
    rtol = 1e-5 if kernel_data["dtype_name"] == "float64" else 1e-4
    assert_allclose(
        K_grad_xp_np,
        K_grad_np,
        rtol=rtol,
        atol=_atol_for_type(kernel_data["dtype_name"]),
    )


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_kernel_namespace_check_array_api(
    array_namespace, device_, dtype_name, kernel_data
):
    """Test that kernel output is in the correct namespace."""
    kernel = RBF()

    with config_context(array_api_dispatch=True):
        K_xp = kernel(kernel_data["X_xp"])

        result_xp, _ = get_namespace(K_xp)
        input_xp, _ = get_namespace(kernel_data["X_xp"])
        assert result_xp.__name__ == input_xp.__name__


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_rbf_kernel_anisotropic_gradient_array_api(
    array_namespace, device_, dtype_name, kernel_data
):
    """Test RBF kernel with anisotropic length scales and array API."""
    # Anisotropic kernel with different length scales per dimension
    kernel = RBF(length_scale=[1.0, 2.0, 0.5])

    K_np, K_grad_np = kernel(kernel_data["X_np"], eval_gradient=True)
    with config_context(array_api_dispatch=True):
        K_xp, K_grad_xp = kernel(kernel_data["X_xp"], eval_gradient=True)

    _check_kernel_result(K_xp, K_np, kernel_data)

    K_grad_xp_np = _convert_to_numpy(K_grad_xp, kernel_data["xp"])
    rtol = 1e-5 if kernel_data["dtype_name"] == "float64" else 1e-4
    assert_allclose(
        K_grad_xp_np,
        K_grad_np,
        rtol=rtol,
        atol=_atol_for_type(kernel_data["dtype_name"]),
    )


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize("nu", [0.5, 1.5, 2.5, np.inf, 3.0])
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_matern_kernel_array_api(nu, array_namespace, device_, dtype_name, kernel_data):
    """Test Matern kernel works with array API inputs for different nu values."""
    kernel = Matern(nu=nu)

    # Test k(X, X)
    K_np = kernel(kernel_data["X_np"])
    with config_context(array_api_dispatch=True):
        K_xp = kernel(kernel_data["X_xp"])

    # Higher tolerance for general nu case which uses Bessel functions
    rtol = 1e-4 if kernel_data["dtype_name"] == "float64" else 1e-3
    _check_kernel_result(K_xp, K_np, kernel_data, rtol=rtol)

    # Test k(X, Y)
    K_np_xy = kernel(kernel_data["X_np"], kernel_data["Y_np"])
    with config_context(array_api_dispatch=True):
        K_xp_xy = kernel(kernel_data["X_xp"], kernel_data["Y_xp"])

    _check_kernel_result(K_xp_xy, K_np_xy, kernel_data, rtol=rtol)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_rbf_kernel_diag_array_api(array_namespace, device_, dtype_name, kernel_data):
    """Test RBF kernel diag method works with array API."""
    kernel = RBF()

    # RBF inherits diag from NormalizedKernelMixin which returns ones
    diag_np = kernel.diag(kernel_data["X_np"])
    with config_context(array_api_dispatch=True):
        diag_xp = kernel.diag(kernel_data["X_xp"])

    diag_xp_np = _convert_to_numpy(diag_xp, kernel_data["xp"])
    assert_allclose(diag_xp_np, diag_np)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_dot_product_diag_array_api(array_namespace, device_, dtype_name, kernel_data):
    """Test DotProduct kernel diag method works with array API inputs."""
    kernel = DotProduct()

    diag_np = kernel.diag(kernel_data["X_np"])
    with config_context(array_api_dispatch=True):
        diag_xp = kernel.diag(kernel_data["X_xp"])

    diag_xp_np = _convert_to_numpy(diag_xp, kernel_data["xp"])
    rtol = 1e-5 if kernel_data["dtype_name"] == "float64" else 1e-4
    assert_allclose(diag_xp_np, diag_np, rtol=rtol)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_constant_kernel_diag_array_api(
    array_namespace, device_, dtype_name, kernel_data
):
    """Test ConstantKernel diag method works with array API inputs."""
    kernel = ConstantKernel(constant_value=2.0)

    diag_np = kernel.diag(kernel_data["X_np"])
    with config_context(array_api_dispatch=True):
        diag_xp = kernel.diag(kernel_data["X_xp"])

    diag_xp_np = _convert_to_numpy(diag_xp, kernel_data["xp"])
    rtol = 1e-5 if kernel_data["dtype_name"] == "float64" else 1e-4
    assert_allclose(diag_xp_np, diag_np, rtol=rtol)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device_, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
def test_white_kernel_diag_array_api(array_namespace, device_, dtype_name, kernel_data):
    """Test WhiteKernel diag method works with array API inputs."""
    kernel = WhiteKernel(noise_level=0.5)

    diag_np = kernel.diag(kernel_data["X_np"])
    with config_context(array_api_dispatch=True):
        diag_xp = kernel.diag(kernel_data["X_xp"])

    diag_xp_np = _convert_to_numpy(diag_xp, kernel_data["xp"])
    rtol = 1e-5 if kernel_data["dtype_name"] == "float64" else 1e-4
    assert_allclose(diag_xp_np, diag_np, rtol=rtol)
