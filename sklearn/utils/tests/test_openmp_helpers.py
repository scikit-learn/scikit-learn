import os
import pytest

from sklearn.utils._openmp_helpers import _openmp_supported


def test_openmp_supported():
    # Check that sklearn is built with OpenMP supported.
    # This test can be disabled by setting the environment variable
    # ``SKLEARN_SKIP_OPENMP_TEST``.
    if os.getenv("SKLEARN_SKIP_OPENMP_TEST"):
        pytest.skip("test explicitly disabled (SKLEARN_SKIP_OPENMP_TEST)")

    msg = ("This test fails because scikit-learn has been built without OpenMP"
           " support. You can skip this test by setting the environment "
           "variable SKLEARN_SKIP_OPENMP_TEST")

    assert _openmp_supported(), msg
