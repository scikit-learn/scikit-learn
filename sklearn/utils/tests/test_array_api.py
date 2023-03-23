import numpy
from numpy.testing import assert_allclose, assert_array_equal
import pytest

from sklearn.base import BaseEstimator
from sklearn.utils._array_api import get_namespace
from sklearn.utils._array_api import _ArrayAPIWrapper

from sklearn.utils._array_api import _asarray_with_order
from sklearn.utils._array_api import _convert_to_numpy
from sklearn.utils._array_api import _estimator_with_converted_arrays

import sklearn.externals._array_api_compat.numpy as array_api_compat_numpy
from sklearn._config import config_context

pytestmark = pytest.mark.filterwarnings(
    "ignore:The numpy.array_api submodule:UserWarning"
)


def test_get_namespace_ndarray():
    """Test get_namespace on NumPy ndarrays."""
    pytest.importorskip("numpy.array_api")

    X_np = numpy.asarray([[1, 2, 3]])

    for array_api_dispatch in [True, False]:
        with config_context(array_api_dispatch=array_api_dispatch):
            xp_out, is_array_api = get_namespace(X_np)
            assert is_array_api == array_api_dispatch
            assert xp_out is array_api_compat_numpy


def test_get_namespace_array_api():
    """Test get_namespace for ArrayAPI arrays."""
    xp = pytest.importorskip("numpy.array_api")

    X_np = numpy.asarray([[1, 2, 3]])
    X_xp = xp.asarray(X_np)
    with config_context(array_api_dispatch=True):
        xp_out, is_array_api = get_namespace(X_xp)
        assert is_array_api
        assert isinstance(xp_out, _ArrayAPIWrapper)

        # check errors
        with pytest.raises(TypeError, match="Multiple namespaces"):
            get_namespace(X_np, X_xp)

        xp_out, is_array_api = get_namespace(1)
        assert xp_out == array_api_compat_numpy
        assert not is_array_api


class _AdjustableNameAPITestWrapper(_ArrayAPIWrapper):
    """API wrapper that has an adjustable name. Used for testing."""

    def __init__(self, array_namespace, name):
        super().__init__(array_namespace=array_namespace)
        self.__name__ = name


def test_array_api_wrapper_astype():
    """Test _ArrayAPIWrapper for ArrayAPIs that is not NumPy."""
    numpy_array_api = pytest.importorskip("numpy.array_api")
    xp_ = _AdjustableNameAPITestWrapper(numpy_array_api, "wrapped_numpy.array_api")
    xp = _ArrayAPIWrapper(xp_)

    X = xp.asarray(([[1, 2, 3], [3, 4, 5]]), dtype=xp.float64)
    X_converted = xp.astype(X, xp.float32)
    assert X_converted.dtype == xp.float32

    X_converted = xp.asarray(X, dtype=xp.float32)
    assert X_converted.dtype == xp.float32


def test_array_api_wrapper_take_for_numpy_api():
    """Test that fast path is called for numpy.array_api."""
    numpy_array_api = pytest.importorskip("numpy.array_api")
    # USe the same name as numpy.array_api
    xp_ = _AdjustableNameAPITestWrapper(numpy_array_api, "numpy.array_api")
    xp = _ArrayAPIWrapper(xp_)

    X = xp.asarray(([[1, 2, 3], [3, 4, 5]]), dtype=xp.float64)
    X_take = xp.take(X, xp.asarray([1]), axis=0)
    assert hasattr(X_take, "__array_namespace__")
    assert_array_equal(X_take, numpy.take(X, [1], axis=0))


def test_array_api_wrapper_take():
    """Test _ArrayAPIWrapper API for take."""
    numpy_array_api = pytest.importorskip("numpy.array_api")
    xp_ = _AdjustableNameAPITestWrapper(numpy_array_api, "wrapped_numpy.array_api")
    xp = _ArrayAPIWrapper(xp_)

    # Check take compared to NumPy's with axis=0
    X_1d = xp.asarray([1, 2, 3], dtype=xp.float64)
    X_take = xp.take(X_1d, xp.asarray([1]), axis=0)
    assert hasattr(X_take, "__array_namespace__")
    assert_array_equal(X_take, numpy.take(X_1d, [1], axis=0))

    X = xp.asarray(([[1, 2, 3], [3, 4, 5]]), dtype=xp.float64)
    X_take = xp.take(X, xp.asarray([0]), axis=0)
    assert hasattr(X_take, "__array_namespace__")
    assert_array_equal(X_take, numpy.take(X, [0], axis=0))

    # Check take compared to NumPy's with axis=1
    X_take = xp.take(X, xp.asarray([0, 2]), axis=1)
    assert hasattr(X_take, "__array_namespace__")
    assert_array_equal(X_take, numpy.take(X, [0, 2], axis=1))

    with pytest.raises(ValueError, match=r"Only axis in \(0, 1\) is supported"):
        xp.take(X, xp.asarray([0]), axis=2)

    with pytest.raises(ValueError, match=r"Only X.ndim in \(1, 2\) is supported"):
        xp.take(xp.asarray([[[0]]]), xp.asarray([0]), axis=0)


@pytest.mark.parametrize("is_array_api", [True, False])
def test_asarray_with_order(is_array_api):
    """Test _asarray_with_order passes along order for NumPy arrays."""
    if is_array_api:
        xp = pytest.importorskip("numpy.array_api")
    else:
        xp = numpy

    X = xp.asarray([1.2, 3.4, 5.1])
    X_new = _asarray_with_order(X, order="F")

    X_new_np = numpy.asarray(X_new)
    assert X_new_np.flags["F_CONTIGUOUS"]


def test_asarray_with_order_ignored():
    """Test _asarray_with_order ignores order for Generic ArrayAPI."""
    xp = pytest.importorskip("numpy.array_api")
    xp_ = _AdjustableNameAPITestWrapper(xp, "wrapped.array_api")

    X = numpy.asarray([[1.2, 3.4, 5.1], [3.4, 5.5, 1.2]], order="C")
    X = xp_.asarray(X)

    X_new = _asarray_with_order(X, order="F", xp=xp_)

    X_new_np = numpy.asarray(X_new)
    assert X_new_np.flags["C_CONTIGUOUS"]
    assert not X_new_np.flags["F_CONTIGUOUS"]


def test_convert_to_numpy_error():
    """Test convert to numpy errors for unsupported namespaces."""
    xp = pytest.importorskip("numpy.array_api")
    xp_ = _AdjustableNameAPITestWrapper(xp, "wrapped.array_api")

    X = xp_.asarray([1.2, 3.4])

    with pytest.raises(ValueError, match="wrapped.array_api is an unsupported"):
        _convert_to_numpy(X, xp=xp_)


@pytest.mark.parametrize("library", ["cupy", "torch", "cupy.array_api"])
def test_convert_to_numpy_gpu(library):
    """Check convert_to_numpy for GPU backed libraries."""
    xp = pytest.importorskip(library)

    if library == "torch":
        if not xp.has_cuda:
            pytest.skip("test requires cuda")
        X_gpu = xp.asarray([1.0, 2.0, 3.0], device="cuda")
    else:
        X_gpu = xp.asarray([1.0, 2.0, 3.0])

    X_cpu = _convert_to_numpy(X_gpu, xp=xp)
    expected_output = numpy.asarray([1.0, 2.0, 3.0])
    assert_allclose(X_cpu, expected_output)


def test_convert_to_numpy_cpu():
    """Check convert_to_numpy for PyTorch CPU arrays."""
    torch = pytest.importorskip("torch")
    X_torch = torch.asarray([1.0, 2.0, 3.0], device="cpu")

    X_cpu = _convert_to_numpy(X_torch, xp=torch)
    expected_output = numpy.asarray([1.0, 2.0, 3.0])
    assert_allclose(X_cpu, expected_output)


class SimpleEstimator(BaseEstimator):
    def fit(self, X, y=None):
        self.X_ = X
        self.n_features_ = X.shape[0]
        return self


@pytest.mark.parametrize(
    "array_namespace, converter",
    [
        ("torch", lambda array: array.cpu().numpy()),
        ("numpy.array_api", lambda array: numpy.asarray(array)),
        ("cupy.array_api", lambda array: array._array.get()),
    ],
)
def test_convert_estimator_to_ndarray(array_namespace, converter):
    """Convert estimator attributes to ndarray."""
    xp = pytest.importorskip(array_namespace)

    X = xp.asarray([[1.3, 4.5]])
    est = SimpleEstimator().fit(X)

    new_est = _estimator_with_converted_arrays(est, converter)
    assert isinstance(new_est.X_, numpy.ndarray)


def test_convert_estimator_to_array_api():
    """Convert estimator attributes to ArrayAPI arrays."""
    xp = pytest.importorskip("numpy.array_api")

    X_np = numpy.asarray([[1.3, 4.5]])
    est = SimpleEstimator().fit(X_np)

    new_est = _estimator_with_converted_arrays(est, lambda array: xp.asarray(array))
    assert hasattr(new_est.X_, "__array_namespace__")


@pytest.mark.parametrize("array_api_dispatch", [True, False])
def test_get_namespace_array_api_isdtype(array_api_dispatch):
    """Test isdtype implementation from _ArrayAPIWrapper and array_api_compat."""
    xp = pytest.importorskip("numpy.array_api")

    X_xp = xp.asarray([[1, 2, 3]])
    with config_context(array_api_dispatch=array_api_dispatch):
        xp_out, _ = get_namespace(X_xp)
        assert xp_out.isdtype(xp_out.float32, "real floating")
        assert xp_out.isdtype(xp_out.float64, "real floating")
        assert not xp_out.isdtype(xp_out.int32, "real floating")

        assert xp_out.isdtype(xp_out.bool, "bool")
        assert not xp_out.isdtype(xp_out.float32, "bool")

        assert xp_out.isdtype(xp_out.int16, "signed integer")
        assert not xp_out.isdtype(xp_out.uint32, "signed integer")

        assert xp_out.isdtype(xp_out.uint16, "unsigned integer")
        assert not xp_out.isdtype(xp_out.int64, "unsigned integer")

        assert xp_out.isdtype(xp_out.int64, "numeric")
        assert xp_out.isdtype(xp_out.float32, "numeric")
        assert xp_out.isdtype(xp_out.uint32, "numeric")


@pytest.mark.parametrize("array_api_dispatch", [True, False])
def test_get_namespace_list(array_api_dispatch):
    """Test get_namespace for lists."""

    X = [1, 2, 3]
    with config_context(array_api_dispatch=array_api_dispatch):
        xp_out, is_array = get_namespace(X)
        assert not is_array
        assert xp_out is array_api_compat_numpy
