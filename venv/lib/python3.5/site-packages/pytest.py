# PYTHON_ARGCOMPLETE_OK
"""
pytest: unit and functional testing with Python.
"""


# else we are imported

from _pytest.config import (
    main, UsageError, cmdline,
    hookspec, hookimpl
)
from _pytest.fixtures import fixture, yield_fixture
from _pytest.assertion import register_assert_rewrite
from _pytest.freeze_support import freeze_includes
from _pytest import __version__
from _pytest.debugging import pytestPDB as __pytestPDB
from _pytest.recwarn import warns, deprecated_call
from _pytest.outcomes import fail, skip, importorskip, exit, xfail
from _pytest.mark import MARK_GEN as mark, param
from _pytest.main import Session
from _pytest.nodes import Item, Collector, File
from _pytest.fixtures import fillfixtures as _fillfuncargs
from _pytest.python import (
    Module, Class, Instance, Function, Generator,
)

from _pytest.python_api import approx, raises

set_trace = __pytestPDB.set_trace

__all__ = [
    'main',
    'UsageError',
    'cmdline',
    'hookspec',
    'hookimpl',
    '__version__',
    'register_assert_rewrite',
    'freeze_includes',
    'set_trace',
    'warns',
    'deprecated_call',
    'fixture',
    'yield_fixture',
    'fail',
    'skip',
    'xfail',
    'importorskip',
    'exit',
    'mark',
    'param',
    'approx',
    '_fillfuncargs',

    'Item',
    'File',
    'Collector',
    'Session',
    'Module',
    'Class',
    'Instance',
    'Function',
    'Generator',
    'raises',


]

if __name__ == '__main__':
    # if run as a script or by 'python -m pytest'
    # we trigger the below "else" condition by the following import
    import pytest
    raise SystemExit(pytest.main())
else:

    from _pytest.compat import _setup_collect_fakemodule
    _setup_collect_fakemodule()
