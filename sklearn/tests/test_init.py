# Basic unittests to test functioning of module's top-level

import subprocess

from sklearn.utils.testing import assert_equal

__author__ = 'Yaroslav Halchenko'
__license__ = 'BSD'


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


def test_import_sklearn_no_warnings():
    # Test that importing scikit-learn doesn't raise any warnings.

    message = subprocess.check_output(['python', '-Wdefault',
                                       '-c', 'import sklearn'],
                                      stderr=subprocess.STDOUT)
    message = message.decode("utf-8")
    # ignore the ImportWarning
    message = '\n'.join([line for line in message.splitlines()
                         if "ImportWarning" not in line])
    assert 'Warning' not in message
    assert 'Error' not in message
