"""
Testing utilities.
"""

import os
import platform
import re
import struct
import sys
import functools
import inspect
from tempfile import NamedTemporaryFile

import numpy as np
from numpy import testing
from numpy.testing import (
    TestCase,
    assert_,
    assert_warns,
    assert_no_warnings,
    assert_equal,
    assert_almost_equal,
    assert_array_equal,
    assert_allclose,
    assert_array_almost_equal,
    assert_array_almost_equal_nulp,
    assert_array_less,
)

from .. import data, io
from ..data._fetchers import _fetch
from ..util import img_as_uint, img_as_float, img_as_int, img_as_ubyte
from ._warnings import expected_warnings
from ._dependency_checks import is_wasm

import pytest


skipif = pytest.mark.skipif
xfail = pytest.mark.xfail
parametrize = pytest.mark.parametrize
raises = pytest.raises
fixture = pytest.fixture

SKIP_RE = re.compile(r"(\s*>>>.*?)(\s*)#\s*skip\s+if\s+(.*)$")

# true if python is running in 32bit mode
# Calculate the size of a void * pointer in bits
# https://docs.python.org/3/library/struct.html
arch32 = struct.calcsize("P") * 8 == 32


def assert_less(a, b, msg=None):
    message = f"{a!r} is not lower than {b!r}"
    if msg is not None:
        message += ": " + msg
    assert a < b, message


def assert_greater(a, b, msg=None):
    message = f"{a!r} is not greater than {b!r}"
    if msg is not None:
        message += ": " + msg
    assert a > b, message


def doctest_skip_parser(func):
    """Decorator replaces custom skip test markup in doctests

    Say a function has a docstring::

        >>> something, HAVE_AMODULE, HAVE_BMODULE = 0, False, False
        >>> something # skip if not HAVE_AMODULE
        0
        >>> something # skip if HAVE_BMODULE
        0

    This decorator will evaluate the expression after ``skip if``.  If this
    evaluates to True, then the comment is replaced by ``# doctest: +SKIP``. If
    False, then the comment is just removed. The expression is evaluated in the
    ``globals`` scope of `func`.

    For example, if the module global ``HAVE_AMODULE`` is False, and module
    global ``HAVE_BMODULE`` is False, the returned function will have docstring::

        >>> something # doctest: +SKIP
        >>> something + else # doctest: +SKIP
        >>> something # doctest: +SKIP

    """
    lines = func.__doc__.split('\n')
    new_lines = []
    for line in lines:
        match = SKIP_RE.match(line)
        if match is None:
            new_lines.append(line)
            continue
        code, space, expr = match.groups()

        try:
            # Works as a function decorator
            if eval(expr, func.__globals__):
                code = code + space + "# doctest: +SKIP"
        except AttributeError:
            # Works as a class decorator
            if eval(expr, func.__init__.__globals__):
                code = code + space + "# doctest: +SKIP"

        new_lines.append(code)
    func.__doc__ = "\n".join(new_lines)
    return func


def roundtrip(image, plugin, suffix):
    """Save and read an image using a specified plugin"""
    if '.' not in suffix:
        suffix = '.' + suffix
    with NamedTemporaryFile(suffix=suffix, delete=False) as temp_file:
        fname = temp_file.name
    io.imsave(fname, image, plugin=plugin)
    new = io.imread(fname, plugin=plugin)
    try:
        os.remove(fname)
    except Exception:
        pass
    return new


def color_check(plugin, fmt='png'):
    """Check roundtrip behavior for color images.

    All major input types should be handled as ubytes and read
    back correctly.
    """
    img = img_as_ubyte(data.chelsea())
    r1 = roundtrip(img, plugin, fmt)
    testing.assert_allclose(img, r1)

    img2 = img > 128
    r2 = roundtrip(img2, plugin, fmt)
    testing.assert_allclose(img2, r2.astype(bool))

    img3 = img_as_float(img)
    r3 = roundtrip(img3, plugin, fmt)
    testing.assert_allclose(r3, img)

    img4 = img_as_int(img)
    if fmt.lower() in (('tif', 'tiff')):
        img4 -= 100
        r4 = roundtrip(img4, plugin, fmt)
        testing.assert_allclose(r4, img4)
    else:
        r4 = roundtrip(img4, plugin, fmt)
        testing.assert_allclose(r4, img_as_ubyte(img4))

    img5 = img_as_uint(img)
    r5 = roundtrip(img5, plugin, fmt)
    testing.assert_allclose(r5, img)


def mono_check(plugin, fmt='png'):
    """Check the roundtrip behavior for images that support most types.

    All major input types should be handled.
    """

    img = img_as_ubyte(data.moon())
    r1 = roundtrip(img, plugin, fmt)
    testing.assert_allclose(img, r1)

    img2 = img > 128
    r2 = roundtrip(img2, plugin, fmt)
    testing.assert_allclose(img2, r2.astype(bool))

    img3 = img_as_float(img)
    r3 = roundtrip(img3, plugin, fmt)
    if r3.dtype.kind == 'f':
        testing.assert_allclose(img3, r3)
    else:
        testing.assert_allclose(r3, img_as_uint(img))

    img4 = img_as_int(img)
    if fmt.lower() in (('tif', 'tiff')):
        img4 -= 100
        r4 = roundtrip(img4, plugin, fmt)
        testing.assert_allclose(r4, img4)
    else:
        r4 = roundtrip(img4, plugin, fmt)
        testing.assert_allclose(r4, img_as_uint(img4))

    img5 = img_as_uint(img)
    r5 = roundtrip(img5, plugin, fmt)
    testing.assert_allclose(r5, img5)


def fetch(data_filename):
    """Attempt to fetch data, but if unavailable, skip the tests."""
    try:
        return _fetch(data_filename)
    except (ConnectionError, ModuleNotFoundError):
        pytest.skip(f'Unable to download {data_filename}', allow_module_level=True)


# Ref: about the lack of threading support in WASM, please see
# https://github.com/pyodide/pyodide/issues/237
def run_in_parallel(num_threads=2, warnings_matching=None):
    """Decorator to run the same function multiple times in parallel.

    This decorator is useful to ensure that separate threads execute
    concurrently and correctly while releasing the GIL.

    It is currently skipped when running on WASM-based platforms, as
    the threading module is not supported.

    Parameters
    ----------
    num_threads : int, optional
        The number of times the function is run in parallel.

    warnings_matching: list or None
        This parameter is passed on to `expected_warnings` so as not to have
        race conditions with the warnings filters. A single
        `expected_warnings` context manager is used for all threads.
        If None, then no warnings are checked.

    """

    assert num_threads > 0

    def wrapper(func):
        if is_wasm:
            # Threading isn't supported on WASM, return early
            return func

        import threading

        @functools.wraps(func)
        def inner(*args, **kwargs):
            with expected_warnings(warnings_matching):
                threads = []
                for i in range(num_threads - 1):
                    thread = threading.Thread(target=func, args=args, kwargs=kwargs)
                    threads.append(thread)
                for thread in threads:
                    thread.start()

                func(*args, **kwargs)

                for thread in threads:
                    thread.join()

        return inner

    return wrapper


def assert_stacklevel(warnings, *, offset=-1):
    """Assert correct stacklevel of captured warnings.

    When scikit-image raises warnings, the stacklevel should ideally be set
    so that the origin of the warnings will point to the public function
    that was called by the user and not necessarily the very place where the
    warnings were emitted (which may be inside some internal function).
    This utility function helps with checking that
    the stacklevel was set correctly on warnings captured by `pytest.warns`.

    Parameters
    ----------
    warnings : collections.abc.Iterable[warning.WarningMessage]
        Warnings that were captured by `pytest.warns`.
    offset : int, optional
        Offset from the line this function is called to the line were the
        warning is supposed to originate from. For multiline calls, the
        first line is relevant. Defaults to -1 which corresponds to the line
        right above the one where this function is called.

    Raises
    ------
    AssertionError
        If a warning in `warnings` does not match the expected line number or
        file name.

    Examples
    --------
    >>> def test_something():
    ...     with pytest.warns(UserWarning, match="some message") as record:
    ...         something_raising_a_warning()
    ...     assert_stacklevel(record)
    ...
    >>> def test_another_thing():
    ...     with pytest.warns(UserWarning, match="some message") as record:
    ...         iam_raising_many_warnings(
    ...             "A long argument that forces the call to wrap."
    ...         )
    ...     assert_stacklevel(record, offset=-3)
    """
    __tracebackhide__ = True  # Hide traceback for py.test

    frame = inspect.stack()[1].frame  # 0 is current frame, 1 is outer frame
    line_number = frame.f_lineno + offset
    filename = frame.f_code.co_filename
    expected = f"{filename}:{line_number}"
    for warning in warnings:
        actual = f"{warning.filename}:{warning.lineno}"
        msg = (
            "Warning with wrong stacklevel:\n"
            f"  Expected: {expected}\n"
            f"  Actual: {actual}\n"
            f"  {warning.category.__name__}: {warning.message}"
        )
        assert actual == expected, msg
