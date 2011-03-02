#emacs: -*- mode: python-mode; py-indent-offset: 4; tab-width: 4; indent-tabs-mode: nil -*-
#ex: set sts=4 ts=4 sw=4 noet:

# Basic unittests to test functioning of module's top-level

__author__ = 'Yaroslav Halchenko'
__license__ = 'BSD'


from nose.tools import assert_true, assert_false, assert_equal, \
    assert_raises

try:
    from scikits.learn import *
    _top_import_error = None
except Exception, e:
    _top_import_error = e

def test_import_skl():
    """Test either above import has failed for some reason

    "import *" is discouraged outside of the module level, hence we
    rely on setting up the variable above
    """
    assert_equal(_top_import_error, None)
