import numpy
from numpy.testing import assert_array_equal
import pytest

from sklearn.base import BaseEstimator
from sklearn.utils._array_api import get_namespace
from sklearn.utils._array_api import _NumPyAPIWrapper
from sklearn.utils._array_api import _ArrayAPIWrapper
from sklearn.utils._array_api import _asarray_with_order
from sklearn.utils._array_api import _convert_to_numpy
from sklearn.utils._array_api import _estimator_with_converted_arrays
from sklearn._config import config_context

pytestmark = pytest.mark.filterwarnings(
    "ignore:The numpy.array_api submodule:UserWarning"
)


def test_get_namespace_ndarray():
    """Test get_namespace on NumPy ndarrays."""
    pytest.importorskip("numpy.array_api")

    X_np = numpy.asarray([[1, 2, 3]])

    # Dispatching on Numpy regardless or the value of array_api_dispatch.
    for array_api_dispatch in [True, False]:
        with config_context(array_api_dispatch=array_api_dispatch):
            xp_out, is_array_api = get_namespace(X_np)
            assert not is_array_api
            assert isinstance(xp_out, _NumPyAPIWrapper)


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
        with pytest.raises(ValueError, match="Multiple namespaces"):
            get_namespace(X_np, X_xp)

        with pytest.raises(ValueError, match="Unrecognized array input"):
            get_namespace(1)


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

    with pytest.raises(ValueError, match="Supported namespaces are:"):
        _convert_to_numpy(X, xp=xp_)


class SimpleEstimator(BaseEstimator):
    def fit(self, X, y=None):
        self.X_ = X
        self.n_features_ = X.shape[0]
        return self


@pytest.mark.parametrize("array_namespace", ["numpy.array_api", "cupy.array_api"])
def test_convert_estimator_to_ndarray(array_namespace):
    """Convert estimator attributes to ndarray."""
    xp = pytest.importorskip(array_namespace)

    if array_namespace == "numpy.array_api":
        converter = lambda array: numpy.asarray(array)  # noqa
    else:  # pragma: no cover
        converter = lambda array: array._array.get()  # noqa

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


@pytest.mark.parametrize("wrapper", [_ArrayAPIWrapper, _NumPyAPIWrapper])
def test_get_namespace_array_api_isdtype(wrapper):
    """Test isdtype implementation from _ArrayAPIWrapper and _NumPyAPIWrapper."""

    if wrapper == _ArrayAPIWrapper:
        xp_ = pytest.importorskip("numpy.array_api")
        xp = _ArrayAPIWrapper(xp_)
    else:
        xp = _NumPyAPIWrapper()

    assert xp.isdtype(xp.float32, xp.float32)
    assert xp.isdtype(xp.float32, "real floating")
    assert xp.isdtype(xp.float64, "real floating")
    assert not xp.isdtype(xp.int32, "real floating")

    assert xp.isdtype(xp.bool, "bool")
    assert not xp.isdtype(xp.float32, "bool")

    assert xp.isdtype(xp.int16, "signed integer")
    assert not xp.isdtype(xp.uint32, "signed integer")

    assert xp.isdtype(xp.uint16, "unsigned integer")
    assert not xp.isdtype(xp.int64, "unsigned integer")

    assert xp.isdtype(xp.int64, "numeric")
    assert xp.isdtype(xp.float32, "numeric")
    assert xp.isdtype(xp.uint32, "numeric")

    assert not xp.isdtype(xp.float32, "complex floating")

    if wrapper == _NumPyAPIWrapper:
        assert not xp.isdtype(xp.int8, "complex floating")
        assert xp.isdtype(xp.complex64, "complex floating")
        assert xp.isdtype(xp.complex128, "complex floating")

    with pytest.raises(ValueError, match="Unrecognized data type"):
        assert xp.isdtype(xp.int16, "unknown")
