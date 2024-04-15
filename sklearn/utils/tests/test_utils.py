import pytest

from sklearn.utils import tosequence
from sklearn.utils._testing import assert_no_warnings


def test_deprecation_joblib_api(tmpdir):
    # Only parallel_backend and register_parallel_backend are not deprecated in
    # sklearn.utils
    from sklearn.utils import parallel_backend, register_parallel_backend

    assert_no_warnings(parallel_backend, "loky", None)
    assert_no_warnings(register_parallel_backend, "failing", None)

    from sklearn.utils._joblib import joblib

    del joblib.parallel.BACKENDS["failing"]


# TODO(1.7): remove
def test_is_pypy_deprecated():
    with pytest.warns(FutureWarning, match="IS_PYPY is deprecated"):
        from sklearn.utils import IS_PYPY  # noqa


# TODO(1.7): remove
def test_tosequence_deprecated():
    with pytest.warns(FutureWarning, match="tosequence was deprecated in 1.5"):
        tosequence([1, 2, 3])
