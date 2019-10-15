# Basic unittests to test functioning of module's top-level


__author__ = 'Yaroslav Halchenko'
__license__ = 'BSD'


import pytest

try:
    from sklearn import *  # noqa
    _top_import_error = None
except Exception as e:
    _top_import_error = e

from sklearn.utils.testing import assert_run_python_script
from sklearn.utils._openmp_helpers import _openmp_supported


def test_import_skl():
    # Test either above import has failed for some reason
    # "import *" is discouraged outside of the module level, hence we
    # rely on setting up the variable above
    assert _top_import_error is None


def test_init_openmp_warning():
    # Check that a warning is printed when sklearn has been built without
    # OpenMP
    if _openmp_supported():
        assert_run_python_script("import sklearn", ignore_init_warnings=False)
    else:
        with pytest.raises(AssertionError,
                           match="Scikit-learn has been built"
                                 " without OpenMP support"):
            assert_run_python_script("import sklearn",
                                     ignore_init_warnings=False)
