""" recording warnings during test function execution. """
from __future__ import absolute_import, division, print_function

import inspect

import _pytest._code
import py
import sys
import warnings

import re

from _pytest.fixtures import yield_fixture
from _pytest.outcomes import fail


@yield_fixture
def recwarn():
    """Return a :class:`WarningsRecorder` instance that records all warnings emitted by test functions.

    See http://docs.python.org/library/warnings.html for information
    on warning categories.
    """
    wrec = WarningsRecorder()
    with wrec:
        warnings.simplefilter('default')
        yield wrec


def deprecated_call(func=None, *args, **kwargs):
    """context manager that can be used to ensure a block of code triggers a
    ``DeprecationWarning`` or ``PendingDeprecationWarning``::

        >>> import warnings
        >>> def api_call_v2():
        ...     warnings.warn('use v3 of this api', DeprecationWarning)
        ...     return 200

        >>> with deprecated_call():
        ...    assert api_call_v2() == 200

    ``deprecated_call`` can also be used by passing a function and ``*args`` and ``*kwargs``,
    in which case it will ensure calling ``func(*args, **kwargs)`` produces one of the warnings
    types above.
    """
    if not func:
        return _DeprecatedCallContext()
    else:
        __tracebackhide__ = True
        with _DeprecatedCallContext():
            return func(*args, **kwargs)


class _DeprecatedCallContext(object):
    """Implements the logic to capture deprecation warnings as a context manager."""

    def __enter__(self):
        self._captured_categories = []
        self._old_warn = warnings.warn
        self._old_warn_explicit = warnings.warn_explicit
        warnings.warn_explicit = self._warn_explicit
        warnings.warn = self._warn

    def _warn_explicit(self, message, category, *args, **kwargs):
        self._captured_categories.append(category)

    def _warn(self, message, category=None, *args, **kwargs):
        if isinstance(message, Warning):
            self._captured_categories.append(message.__class__)
        else:
            self._captured_categories.append(category)

    def __exit__(self, exc_type, exc_val, exc_tb):
        warnings.warn_explicit = self._old_warn_explicit
        warnings.warn = self._old_warn

        if exc_type is None:
            deprecation_categories = (DeprecationWarning, PendingDeprecationWarning)
            if not any(issubclass(c, deprecation_categories) for c in self._captured_categories):
                __tracebackhide__ = True
                msg = "Did not produce DeprecationWarning or PendingDeprecationWarning"
                raise AssertionError(msg)


def warns(expected_warning, *args, **kwargs):
    """Assert that code raises a particular class of warning.

    Specifically, the parameter ``expected_warning`` can be a warning class or
    sequence of warning classes, and the inside the ``with`` block must issue a warning of that class or
    classes.

    This helper produces a list of :class:`warnings.WarningMessage` objects,
    one for each warning raised.

    This function can be used as a context manager, or any of the other ways
    ``pytest.raises`` can be used::

        >>> with warns(RuntimeWarning):
        ...    warnings.warn("my warning", RuntimeWarning)

    In the context manager form you may use the keyword argument ``match`` to assert
    that the exception matches a text or regex::

        >>> with warns(UserWarning, match='must be 0 or None'):
        ...     warnings.warn("value must be 0 or None", UserWarning)

        >>> with warns(UserWarning, match=r'must be \d+$'):
        ...     warnings.warn("value must be 42", UserWarning)

        >>> with warns(UserWarning, match=r'must be \d+$'):
        ...     warnings.warn("this is not here", UserWarning)
        Traceback (most recent call last):
          ...
        Failed: DID NOT WARN. No warnings of type ...UserWarning... was emitted...

    """
    match_expr = None
    if not args:
        if "match" in kwargs:
            match_expr = kwargs.pop("match")
        return WarningsChecker(expected_warning, match_expr=match_expr)
    elif isinstance(args[0], str):
        code, = args
        assert isinstance(code, str)
        frame = sys._getframe(1)
        loc = frame.f_locals.copy()
        loc.update(kwargs)

        with WarningsChecker(expected_warning, match_expr=match_expr):
            code = _pytest._code.Source(code).compile()
            py.builtin.exec_(code, frame.f_globals, loc)
    else:
        func = args[0]
        with WarningsChecker(expected_warning, match_expr=match_expr):
            return func(*args[1:], **kwargs)


class WarningsRecorder(warnings.catch_warnings):
    """A context manager to record raised warnings.

    Adapted from `warnings.catch_warnings`.
    """

    def __init__(self):
        super(WarningsRecorder, self).__init__(record=True)
        self._entered = False
        self._list = []

    @property
    def list(self):
        """The list of recorded warnings."""
        return self._list

    def __getitem__(self, i):
        """Get a recorded warning by index."""
        return self._list[i]

    def __iter__(self):
        """Iterate through the recorded warnings."""
        return iter(self._list)

    def __len__(self):
        """The number of recorded warnings."""
        return len(self._list)

    def pop(self, cls=Warning):
        """Pop the first recorded warning, raise exception if not exists."""
        for i, w in enumerate(self._list):
            if issubclass(w.category, cls):
                return self._list.pop(i)
        __tracebackhide__ = True
        raise AssertionError("%r not found in warning list" % cls)

    def clear(self):
        """Clear the list of recorded warnings."""
        self._list[:] = []

    def __enter__(self):
        if self._entered:
            __tracebackhide__ = True
            raise RuntimeError("Cannot enter %r twice" % self)
        self._list = super(WarningsRecorder, self).__enter__()
        warnings.simplefilter('always')
        return self

    def __exit__(self, *exc_info):
        if not self._entered:
            __tracebackhide__ = True
            raise RuntimeError("Cannot exit %r without entering first" % self)
        super(WarningsRecorder, self).__exit__(*exc_info)


class WarningsChecker(WarningsRecorder):
    def __init__(self, expected_warning=None, match_expr=None):
        super(WarningsChecker, self).__init__()

        msg = ("exceptions must be old-style classes or "
               "derived from Warning, not %s")
        if isinstance(expected_warning, tuple):
            for exc in expected_warning:
                if not inspect.isclass(exc):
                    raise TypeError(msg % type(exc))
        elif inspect.isclass(expected_warning):
            expected_warning = (expected_warning,)
        elif expected_warning is not None:
            raise TypeError(msg % type(expected_warning))

        self.expected_warning = expected_warning
        self.match_expr = match_expr

    def __exit__(self, *exc_info):
        super(WarningsChecker, self).__exit__(*exc_info)

        # only check if we're not currently handling an exception
        if all(a is None for a in exc_info):
            if self.expected_warning is not None:
                if not any(issubclass(r.category, self.expected_warning)
                           for r in self):
                    __tracebackhide__ = True
                    fail("DID NOT WARN. No warnings of type {0} was emitted. "
                         "The list of emitted warnings is: {1}.".format(
                             self.expected_warning,
                             [each.message for each in self]))
                elif self.match_expr is not None:
                    for r in self:
                        if issubclass(r.category, self.expected_warning):
                            if re.compile(self.match_expr).search(str(r.message)):
                                break
                    else:
                        fail("DID NOT WARN. No warnings of type {0} matching"
                             " ('{1}') was emitted. The list of emitted warnings"
                             " is: {2}.".format(self.expected_warning, self.match_expr,
                                                [each.message for each in self]))
