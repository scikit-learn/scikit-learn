from __future__ import absolute_import, division, print_function

import warnings
from contextlib import contextmanager

import pytest

from _pytest import compat


def _setoption(wmod, arg):
    """
    Copy of the warning._setoption function but does not escape arguments.
    """
    parts = arg.split(':')
    if len(parts) > 5:
        raise wmod._OptionError("too many fields (max 5): %r" % (arg,))
    while len(parts) < 5:
        parts.append('')
    action, message, category, module, lineno = [s.strip()
                                                 for s in parts]
    action = wmod._getaction(action)
    category = wmod._getcategory(category)
    if lineno:
        try:
            lineno = int(lineno)
            if lineno < 0:
                raise ValueError
        except (ValueError, OverflowError):
            raise wmod._OptionError("invalid lineno %r" % (lineno,))
    else:
        lineno = 0
    wmod.filterwarnings(action, message, category, module, lineno)


def pytest_addoption(parser):
    group = parser.getgroup("pytest-warnings")
    group.addoption(
        '-W', '--pythonwarnings', action='append',
        help="set which warnings to report, see -W option of python itself.")
    parser.addini("filterwarnings", type="linelist",
                  help="Each line specifies a pattern for "
                  "warnings.filterwarnings. "
                  "Processed after -W and --pythonwarnings.")


@contextmanager
def catch_warnings_for_item(item):
    """
    catches the warnings generated during setup/call/teardown execution
    of the given item and after it is done posts them as warnings to this
    item.
    """
    args = item.config.getoption('pythonwarnings') or []
    inifilters = item.config.getini("filterwarnings")
    with warnings.catch_warnings(record=True) as log:
        for arg in args:
            warnings._setoption(arg)

        for arg in inifilters:
            _setoption(warnings, arg)

        mark = item.get_marker('filterwarnings')
        if mark:
            for arg in mark.args:
                warnings._setoption(arg)

        yield

        for warning in log:
            warn_msg = warning.message
            unicode_warning = False

            if compat._PY2 and any(isinstance(m, compat.UNICODE_TYPES) for m in warn_msg.args):
                new_args = []
                for m in warn_msg.args:
                    new_args.append(compat.ascii_escaped(m) if isinstance(m, compat.UNICODE_TYPES) else m)
                unicode_warning = list(warn_msg.args) != new_args
                warn_msg.args = new_args

            msg = warnings.formatwarning(
                warn_msg, warning.category,
                warning.filename, warning.lineno, warning.line)
            item.warn("unused", msg)

            if unicode_warning:
                warnings.warn(
                    "Warning is using unicode non convertible to ascii, "
                    "converting to a safe representation:\n  %s" % msg,
                    UnicodeWarning)


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_protocol(item):
    with catch_warnings_for_item(item):
        yield
