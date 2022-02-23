import numpy
import pytest

from sklearn.utils._array_api import get_namespace
from sklearn.utils._array_api import _NumPyApiWrapper
from sklearn.utils._array_api import _ArrayAPIWrapper
from sklearn._config import config_context


def test_get_namespace():
    pytest.importorskip("numpy", minversion="1.22", reason="Requires Array API")
    xp = pytest.importorskip("numpy.array_api")

    X_np = numpy.asarray([[1, 2, 3]])
    xp_out, is_array_api = get_namespace(X_np)
    assert not is_array_api
    assert isinstance(xp_out, _NumPyApiWrapper)

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
