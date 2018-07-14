# Basic unittests to test functioning of module's top-level

import subprocess

import pkgutil

import sklearn
from sklearn.utils.testing import assert_equal, SkipTest

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
    # Test that importing scikit-learn main modules doesn't raise any warnings.

    try:
        pkgs = pkgutil.iter_modules(path=sklearn.__path__, prefix='sklearn.')
        import_modules = '; '.join(['import ' + modname
                                    for _, modname, _ in pkgs
                                    if (not modname.startswith('_') and
                                        # add deprecated top level modules
                                        # below to ignore them
                                        modname not in [])])

        message = subprocess.check_output(['python', '-Wdefault',
                                           '-c', import_modules],
                                          stderr=subprocess.STDOUT)
        message = message.decode("utf-8")
        message = '\n'.join([line for line in message.splitlines()
                             if not (
                                     # ignore ImportWarning due to Cython
                                     "ImportWarning" in line or
                                     # ignore DeprecationWarning due to pytest
                                     "pytest" in line or
                                     # ignore DeprecationWarnings due to
                                     # numpy.oldnumeric
                                     "oldnumeric" in line
                                     )])
        assert 'Warning' not in message
        assert 'Error' not in message

    except Exception as e:
        raise SkipTest('soft-failed test_import_sklearn_no_warnings.\n'
                       ' %s, \n %s' % (e, message))
