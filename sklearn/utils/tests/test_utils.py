import joblib
import pytest

from sklearn.utils import parallel_backend, register_parallel_backend, tosequence


# TODO(1.7): remove
def test_is_pypy_deprecated():
    with pytest.warns(FutureWarning, match="IS_PYPY is deprecated"):
        from sklearn.utils import IS_PYPY  # noqa


# TODO(1.7): remove
def test_tosequence_deprecated():
    with pytest.warns(FutureWarning, match="tosequence was deprecated in 1.5"):
        tosequence([1, 2, 3])


# TODO(1.7): remove
def test_parallel_backend_deprecated():
    with pytest.warns(FutureWarning, match="parallel_backend is deprecated"):
        parallel_backend("loky", None)

    with pytest.warns(FutureWarning, match="register_parallel_backend is deprecated"):
        register_parallel_backend("a_backend", None)

    del joblib.parallel.BACKENDS["a_backend"]
