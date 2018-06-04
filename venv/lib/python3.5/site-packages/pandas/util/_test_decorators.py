"""
This module provides decorator functions which can be applied to test objects
in order to skip those objects when certain conditions occur. A sample use case
is to detect if the platform is missing ``matplotlib``. If so, any test objects
which require ``matplotlib`` and decorated with ``@td.skip_if_no_mpl`` will be
skipped by ``pytest`` during the execution of the test suite.

To illustrate, after importing this module:

import pandas.util._test_decorators as td

The decorators can be applied to classes:

@td.skip_if_some_reason
class Foo():
    ...

Or individual functions:

@td.skip_if_some_reason
def test_foo():
    ...

For more information, refer to the ``pytest`` documentation on ``skipif``.
"""

import pytest
import locale
from distutils.version import LooseVersion

from pandas.compat import (is_platform_windows, is_platform_32bit, PY3,
                           import_lzma)
from pandas.core.computation.expressions import (_USE_NUMEXPR,
                                                 _NUMEXPR_INSTALLED)


def safe_import(mod_name, min_version=None):
    """
    Parameters:
    -----------
    mod_name : str
        Name of the module to be imported
    min_version : str, default None
        Minimum required version of the specified mod_name

    Returns:
    --------
    object
        The imported module if successful, or False
    """
    try:
        mod = __import__(mod_name)
    except ImportError:
        return False

    if not min_version:
        return mod
    else:
        import sys
        try:
            version = getattr(sys.modules[mod_name], '__version__')
        except AttributeError:
            # xlrd uses a capitalized attribute name
            version = getattr(sys.modules[mod_name], '__VERSION__')
        if version:
            from distutils.version import LooseVersion
            if LooseVersion(version) >= LooseVersion(min_version):
                return mod

    return False


def _skip_if_no_mpl():
    mod = safe_import("matplotlib")
    if mod:
        mod.use("Agg", warn=False)
    else:
        return True


def _skip_if_mpl_1_5():
    mod = safe_import("matplotlib")

    if mod:
        v = mod.__version__
        if LooseVersion(v) > LooseVersion('1.4.3') or str(v)[0] == '0':
            return True
        else:
            mod.use("Agg", warn=False)


def _skip_if_mpl_2_2():
    mod = safe_import("matplotlib")

    if mod:
        v = mod.__version__
        if LooseVersion(v) > LooseVersion('2.1.2'):
            return True
        else:
            mod.use("Agg", warn=False)


def _skip_if_has_locale():
    lang, _ = locale.getlocale()
    if lang is not None:
        return True


def _skip_if_not_us_locale():
    lang, _ = locale.getlocale()
    if lang != 'en_US':
        return True


def _skip_if_no_scipy():
    return not (safe_import('scipy.stats') and safe_import('scipy.sparse') and
                safe_import('scipy.interpolate'))


def _skip_if_no_lzma():
    try:
        import_lzma()
    except ImportError:
        return True


def skip_if_no(package, min_version=None):
    """
    Generic function to help skip test functions when required packages are not
    present on the testing system.

    Intended for use as a decorator, this function will wrap the decorated
    function with a pytest ``skip_if`` mark. During a pytest test suite
    execution, that mark will attempt to import the specified ``package`` and
    optionally ensure it meets the ``min_version``. If the import and version
    check are unsuccessful, then the decorated function will be skipped.

    Parameters
    ----------
    package: str
        The name of the package required by the decorated function
    min_version: str or None, default None
        Optional minimum version of the package required by the decorated
        function

    Returns
    -------
    decorated_func: function
        The decorated function wrapped within a pytest ``skip_if`` mark
    """
    def decorated_func(func):
        msg = "Could not import '{}'".format(package)
        if min_version:
            msg += " satisfying a min_version of {}".format(min_version)
        return pytest.mark.skipif(
            not safe_import(package, min_version=min_version), reason=msg
        )(func)
    return decorated_func


skip_if_no_mpl = pytest.mark.skipif(_skip_if_no_mpl(),
                                    reason="Missing matplotlib dependency")
skip_if_mpl = pytest.mark.skipif(not _skip_if_no_mpl(),
                                 reason="matplotlib is present")
skip_if_mpl_1_5 = pytest.mark.skipif(_skip_if_mpl_1_5(),
                                     reason="matplotlib 1.5")
xfail_if_mpl_2_2 = pytest.mark.xfail(_skip_if_mpl_2_2(),
                                     reason="matplotlib 2.2")
skip_if_32bit = pytest.mark.skipif(is_platform_32bit(),
                                   reason="skipping for 32 bit")
skip_if_windows = pytest.mark.skipif(is_platform_windows(),
                                     reason="Running on Windows")
skip_if_windows_python_3 = pytest.mark.skipif(is_platform_windows() and PY3,
                                              reason=("not used on python3/"
                                                      "win32"))
skip_if_has_locale = pytest.mark.skipif(_skip_if_has_locale(),
                                        reason="Specific locale is set {lang}"
                                        .format(lang=locale.getlocale()[0]))
skip_if_not_us_locale = pytest.mark.skipif(_skip_if_not_us_locale(),
                                           reason="Specific locale is set "
                                           "{lang}".format(
                                               lang=locale.getlocale()[0]))
skip_if_no_scipy = pytest.mark.skipif(_skip_if_no_scipy(),
                                      reason="Missing SciPy requirement")
skip_if_no_lzma = pytest.mark.skipif(_skip_if_no_lzma(),
                                     reason="need backports.lzma to run")
skip_if_no_ne = pytest.mark.skipif(not _USE_NUMEXPR,
                                   reason="numexpr enabled->{enabled}, "
                                   "installed->{installed}".format(
                                       enabled=_USE_NUMEXPR,
                                       installed=_NUMEXPR_INSTALLED))


def parametrize_fixture_doc(*args):
    """
    Intended for use as a decorator for parametrized fixture,
    this function will wrap the decorated function with a pytest
    ``parametrize_fixture_doc`` mark. That mark will format
    initial fixture docstring by replacing placeholders {0}, {1} etc
    with parameters passed as arguments.

    Parameters:
    ----------
        args: iterable
            Positional arguments for docstring.

    Returns:
    -------
    documented_fixture: function
        The decorated function wrapped within a pytest
        ``parametrize_fixture_doc`` mark
    """
    def documented_fixture(fixture):
        fixture.__doc__ = fixture.__doc__.format(*args)
        return fixture
    return documented_fixture
