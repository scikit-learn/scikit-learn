# Basic unittests to test functioning of module's top-level

__author__ = 'Yaroslav Halchenko'
__license__ = 'BSD'


from sklearn.utils.testing import assert_equal

try:
    from sklearn import *  # noqa
    _top_import_error = None
except Exception as e:
    _top_import_error = e


def test_import_skl():
    # Test either above import has failed for some reason
    # "import *" is discouraged outside of the module level, hence we
    # rely on setting up the variable above
    assert_equal(_top_import_error, None)
