import numpy
from numpy.testing import assert_array_equal
import pytest

from sklearn.utils._array_api import get_namespace
from sklearn.utils._array_api import _NumPyApiWrapper
from sklearn.utils._array_api import _ArrayAPIWrapper
from sklearn.utils._array_api import _asarray_with_order
from sklearn.utils._array_api import _convert_to_numpy
from sklearn.utils._array_api import _convert_estimator_to_ndarray
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
            assert isinstance(xp_out, _NumPyApiWrapper)


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


class _NumPyArrayAPITestWrapper(_ArrayAPIWrapper):
    """NumPy Array API wrapper that has a different name. Used for testing."""

    __name__ = "wrapped_numpy.array_api"

    def __init__(self, array_namespace):
        super().__init__(array_namespace=array_namespace)

    def __getattr__(self, name):
        return getattr(self._namespace, name)


def test_array_api_wrapper():
    """Test _ArrayAPIWrapper for ArrayAPIs that is not NumPy."""
    xp = pytest.importorskip("numpy.array_api")
    xp_ = _NumPyArrayAPITestWrapper(xp)

    X = xp.asarray(([[1, 2, 3], [3, 4, 5]]), dtype=xp.float64)
    X_converted = xp_.astype(X, xp.float32)
    assert X_converted.dtype == xp.float32

    X = xp.asarray(([[1, 2, 3], [3, 4, 5]]), dtype=xp.float64)
    X_converted = xp_.asarray(X, dtype=xp.float32)
    assert X_converted.dtype == xp.float32

    # Check take compared to NumPy's with axis=0
    X_take = xp_.take(X, xp.asarray([0]), axis=0)
    assert_array_equal(X_take, numpy.take(X, [0], axis=0))

    # Check take compared to NumPy's with axis=1
    X_take = xp_.take(X, xp.asarray([0, 2]), axis=1)
    assert_array_equal(X_take, numpy.take(X, [0, 2], axis=1))


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
    xp_ = _NumPyArrayAPITestWrapper(xp)

    X = numpy.asarray([[1.2, 3.4, 5.1], [3.4, 5.5, 1.2]], order="C")
    X = xp_.asarray(X)

    X_new = _asarray_with_order(X, order="F", xp=xp_)

    X_new_np = numpy.asarray(X_new)
    assert X_new_np.flags["C_CONTIGUOUS"]
    assert not X_new_np.flags["F_CONTIGUOUS"]


def test_convert_to_numpy_error():
    """Test convert to numpy errors for unsupported namespaces."""
    xp = pytest.importorskip("numpy.array_api")
    xp_ = _NumPyArrayAPITestWrapper(xp)

    X = xp_.asarray([1.2, 3.4])

    with pytest.raises(ValueError, match="Supported namespaces are:"):
        _convert_to_numpy(X, xp=xp_)


@pytest.mark.parametrize("array_namespace", ["numpy.array_api", "cupy.array_api"])
def test_convert_estimator_to_ndarray(array_namespace):
    xp = pytest.importorskip(array_namespace)

    class Estimator:
        def fit(self, X, y=None):
            self.X_ = X
            self.n_features_ = X.shape[0]
            return self

    X = xp.asarray([[1.3, 4.5]])
    est = Estimator().fit(X)

    _convert_estimator_to_ndarray(est)
    assert isinstance(est.X_, numpy.ndarray)
