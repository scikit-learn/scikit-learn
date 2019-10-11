# Basic unittests to test functioning of module's top-level


__author__ = 'Yaroslav Halchenko'
__license__ = 'BSD'


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


def test_import_openmp_warning():
    # Check that a warning is printed when sklearn has been built without
    # OpenMP
    try:
        assert_run_python_script("import sklearn", ignore_warnings=False)
    except AssertionError as err:
        match = "Scikit-learn has been built without OpenMP support"
        assert not _openmp_supported() and match in str(err)
