import os
import pytest
import textwrap

from sklearn import __version__
from sklearn.utils._openmp_helpers import _openmp_supported


def test_openmp_supported():
    # Check that sklearn is built with OpenMP supported.
    # This test can be disabled by setting the environment variable
    # ``SKLEARN_SKIP_OPENMP_TEST``.
    if os.getenv("SKLEARN_SKIP_OPENMP_TEST"):
        pytest.skip("test explicitly disabled (SKLEARN_SKIP_OPENMP_TEST)")

    base_url = "dev" if __version__.endswith(".dev0") else "stable"
    err_msg = textwrap.dedent(
            """
            This test fails because scikit-learn has been built without OpenMP
            support. This is not recommended since some estimators will run in
            sequential mode instead of leveraging thread-based parallelism.

            You can find instructions to build scikit-learn with OpenMP support
            at this adress:

                https://scikit-learn.org/{}/developers/advanced_installation.html

            You can skip this test by setting the environment variable
            SKLEARN_SKIP_OPENMP_TEST to any value.
            """).format(base_url)

    assert _openmp_supported(), err_msg
