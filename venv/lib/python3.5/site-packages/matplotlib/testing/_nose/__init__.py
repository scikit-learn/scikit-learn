from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import sys


def get_extra_test_plugins():
    from .plugins.performgc import PerformGC
    from .plugins.knownfailure import KnownFailure
    from nose.plugins import attrib

    return [PerformGC, KnownFailure, attrib.Plugin]


def get_env():
    env = {'NOSE_COVER_PACKAGE': ['matplotlib', 'mpl_toolkits'],
           'NOSE_COVER_HTML': 1,
           'NOSE_COVER_NO_PRINT': 1}
    return env


def check_deps():
    try:
        import nose
        try:
            from unittest import mock
        except ImportError:
            import mock
    except ImportError:
        print("matplotlib.test requires nose and mock to run.")
        raise


def test(verbosity=None, coverage=False, switch_backend_warn=True,
         recursionlimit=0, **kwargs):
    from ... import default_test_modules, get_backend, use

    old_backend = get_backend()
    old_recursionlimit = sys.getrecursionlimit()
    try:
        use('agg')
        if recursionlimit:
            sys.setrecursionlimit(recursionlimit)
        import nose
        from nose.plugins import multiprocess

        # Nose doesn't automatically instantiate all of the plugins in the
        # child processes, so we have to provide the multiprocess plugin with
        # a list.
        extra_plugins = get_extra_test_plugins()
        multiprocess._instantiate_plugins = extra_plugins

        env = get_env()
        if coverage:
            env['NOSE_WITH_COVERAGE'] = 1

        if verbosity is not None:
            env['NOSE_VERBOSE'] = verbosity

        success = nose.run(
            addplugins=[plugin() for plugin in extra_plugins],
            env=env,
            defaultTest=default_test_modules,
            **kwargs
        )
    finally:
        if old_backend.lower() != 'agg':
            use(old_backend, warn=switch_backend_warn)
        if recursionlimit:
            sys.setrecursionlimit(old_recursionlimit)

    return success


def knownfail(msg):
    from .exceptions import KnownFailureTest
    # Keep the next ultra-long comment so it shows in console.
    raise KnownFailureTest(msg)  # An error here when running nose means that you don't have the matplotlib.testing.nose.plugins:KnownFailure plugin in use.  # noqa
