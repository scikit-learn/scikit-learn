import numpy
from numpy.testing import assert_array_equal
import pytest

from sklearn.utils._array_api import get_namespace
from sklearn.utils._array_api import _NumPyApiWrapper
from sklearn.utils._array_api import _ArrayAPIWrapper
from sklearn._config import config_context


@pytest.mark.filterwarnings("ignore:The numpy.array_api submodule:UserWarning")
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


@pytest.mark.filterwarnings("ignore:The numpy.array_api submodule:UserWarning")
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


class _NumPyArrayAPIWrapper:
    __name__ = "wrapped_numpy.array_api"

    def __init__(self, array_namespace):
        self._namespace = array_namespace

    def __getattr__(self, name):
        return getattr(self._namespace, name)


@pytest.mark.filterwarnings("ignore:The numpy.array_api submodule:UserWarning")
def test_array_api_wrapper():
    """Test _ArrayAPIWrapper for ArrayAPIs that is not NumPy."""
    xp = pytest.importorskip("numpy.array_api")

    wrapped_np = _NumPyArrayAPIWrapper(xp)
    xp_ = _ArrayAPIWrapper(wrapped_np)

    X = xp.asarray(([[1, 2, 3], [3, 4, 5]]), dtype=xp.float64)
    X_converted = xp_.astype(X, xp.float32)
    assert X_converted.dtype == xp.float32

    X = xp.asarray(([[1, 2, 3], [3, 4, 5]]), dtype=xp.float64)
    X_converted = xp_.asarray(X, dtype=xp.float32)
    assert X_converted.dtype == xp.float32

    # Check take compared to NumPy's
    X_take = xp_.take(X, xp.asarray([0]), axis=0)
    assert_array_equal(X_take, numpy.take(X, [0], axis=0))

    # Check take compared to NumPy's
    X_take = xp_.take(X, xp.asarray([0, 2]), axis=1)
    assert_array_equal(X_take, numpy.take(X, [0, 2], axis=1))


def test_array_api_raises():
    """Check that _ArrayAPIWrapper raises for non NumPy ArrayAPI arrays."""
    xp = pytest.importorskip("numpy.array_api")

    wrapped_np = _NumPyArrayAPIWrapper(xp)
    xp_ = _ArrayAPIWrapper(wrapped_np)

    X = xp.asarray(([[1, 2, 3], [3, 4, 5]]), dtype=xp.float64)

    msg = "Only NumPy's Array API implementation supports casting != 'unsafe'"
    with pytest.raises(ValueError, match=msg):
        xp_.astype(X, dtype=xp.float32, casting="same_kind")

    msg = "Only NumPy's Array API implementation supports order"
    with pytest.raises(ValueError, match=msg):
        xp_.asarray(X, order="F")
