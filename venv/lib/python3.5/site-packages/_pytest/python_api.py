import math
import sys

import py
from six import binary_type, text_type
from six.moves import zip, filterfalse
from more_itertools.more import always_iterable

from _pytest.compat import isclass
from _pytest.outcomes import fail
import _pytest._code


def _cmp_raises_type_error(self, other):
    """__cmp__ implementation which raises TypeError. Used
    by Approx base classes to implement only == and != and raise a
    TypeError for other comparisons.

    Needed in Python 2 only, Python 3 all it takes is not implementing the
    other operators at all.
    """
    __tracebackhide__ = True
    raise TypeError('Comparison operators other than == and != not supported by approx objects')


# builtin pytest.approx helper


class ApproxBase(object):
    """
    Provide shared utilities for making approximate comparisons between numbers
    or sequences of numbers.
    """

    # Tell numpy to use our `__eq__` operator instead of its
    __array_ufunc__ = None
    __array_priority__ = 100

    def __init__(self, expected, rel=None, abs=None, nan_ok=False):
        self.expected = expected
        self.abs = abs
        self.rel = rel
        self.nan_ok = nan_ok

    def __repr__(self):
        raise NotImplementedError

    def __eq__(self, actual):
        return all(
            a == self._approx_scalar(x)
            for a, x in self._yield_comparisons(actual))

    __hash__ = None

    def __ne__(self, actual):
        return not (actual == self)

    if sys.version_info[0] == 2:
        __cmp__ = _cmp_raises_type_error

    def _approx_scalar(self, x):
        return ApproxScalar(x, rel=self.rel, abs=self.abs, nan_ok=self.nan_ok)

    def _yield_comparisons(self, actual):
        """
        Yield all the pairs of numbers to be compared.  This is used to
        implement the `__eq__` method.
        """
        raise NotImplementedError


class ApproxNumpy(ApproxBase):
    """
    Perform approximate comparisons for numpy arrays.
    """

    def __repr__(self):
        # It might be nice to rewrite this function to account for the
        # shape of the array...
        import numpy as np

        return "approx({0!r})".format(list(
            self._approx_scalar(x) for x in np.asarray(self.expected)))

    if sys.version_info[0] == 2:
        __cmp__ = _cmp_raises_type_error

    def __eq__(self, actual):
        import numpy as np

        # self.expected is supposed to always be an array here

        if not np.isscalar(actual):
            try:
                actual = np.asarray(actual)
            except:  # noqa
                raise TypeError("cannot compare '{0}' to numpy.ndarray".format(actual))

        if not np.isscalar(actual) and actual.shape != self.expected.shape:
            return False

        return ApproxBase.__eq__(self, actual)

    def _yield_comparisons(self, actual):
        import numpy as np

        # `actual` can either be a numpy array or a scalar, it is treated in
        # `__eq__` before being passed to `ApproxBase.__eq__`, which is the
        # only method that calls this one.

        if np.isscalar(actual):
            for i in np.ndindex(self.expected.shape):
                yield actual, np.asscalar(self.expected[i])
        else:
            for i in np.ndindex(self.expected.shape):
                yield np.asscalar(actual[i]), np.asscalar(self.expected[i])


class ApproxMapping(ApproxBase):
    """
    Perform approximate comparisons for mappings where the values are numbers
    (the keys can be anything).
    """

    def __repr__(self):
        return "approx({0!r})".format(dict(
            (k, self._approx_scalar(v))
            for k, v in self.expected.items()))

    def __eq__(self, actual):
        if set(actual.keys()) != set(self.expected.keys()):
            return False

        return ApproxBase.__eq__(self, actual)

    def _yield_comparisons(self, actual):
        for k in self.expected.keys():
            yield actual[k], self.expected[k]


class ApproxSequence(ApproxBase):
    """
    Perform approximate comparisons for sequences of numbers.
    """

    def __repr__(self):
        seq_type = type(self.expected)
        if seq_type not in (tuple, list, set):
            seq_type = list
        return "approx({0!r})".format(seq_type(
            self._approx_scalar(x) for x in self.expected))

    def __eq__(self, actual):
        if len(actual) != len(self.expected):
            return False
        return ApproxBase.__eq__(self, actual)

    def _yield_comparisons(self, actual):
        return zip(actual, self.expected)


class ApproxScalar(ApproxBase):
    """
    Perform approximate comparisons for single numbers only.
    """
    DEFAULT_ABSOLUTE_TOLERANCE = 1e-12
    DEFAULT_RELATIVE_TOLERANCE = 1e-6

    def __repr__(self):
        """
        Return a string communicating both the expected value and the tolerance
        for the comparison being made, e.g. '1.0 +- 1e-6'.  Use the unicode
        plus/minus symbol if this is python3 (it's too hard to get right for
        python2).
        """
        if isinstance(self.expected, complex):
            return str(self.expected)

        # Infinities aren't compared using tolerances, so don't show a
        # tolerance.
        if math.isinf(self.expected):
            return str(self.expected)

        # If a sensible tolerance can't be calculated, self.tolerance will
        # raise a ValueError.  In this case, display '???'.
        try:
            vetted_tolerance = '{:.1e}'.format(self.tolerance)
        except ValueError:
            vetted_tolerance = '???'

        if sys.version_info[0] == 2:
            return '{0} +- {1}'.format(self.expected, vetted_tolerance)
        else:
            return u'{0} \u00b1 {1}'.format(self.expected, vetted_tolerance)

    def __eq__(self, actual):
        """
        Return true if the given value is equal to the expected value within
        the pre-specified tolerance.
        """
        if _is_numpy_array(actual):
            return ApproxNumpy(actual, self.abs, self.rel, self.nan_ok) == self.expected

        # Short-circuit exact equality.
        if actual == self.expected:
            return True

        # Allow the user to control whether NaNs are considered equal to each
        # other or not.  The abs() calls are for compatibility with complex
        # numbers.
        if math.isnan(abs(self.expected)):
            return self.nan_ok and math.isnan(abs(actual))

        # Infinity shouldn't be approximately equal to anything but itself, but
        # if there's a relative tolerance, it will be infinite and infinity
        # will seem approximately equal to everything.  The equal-to-itself
        # case would have been short circuited above, so here we can just
        # return false if the expected value is infinite.  The abs() call is
        # for compatibility with complex numbers.
        if math.isinf(abs(self.expected)):
            return False

        # Return true if the two numbers are within the tolerance.
        return abs(self.expected - actual) <= self.tolerance

    __hash__ = None

    @property
    def tolerance(self):
        """
        Return the tolerance for the comparison.  This could be either an
        absolute tolerance or a relative tolerance, depending on what the user
        specified or which would be larger.
        """
        def set_default(x, default):
            return x if x is not None else default

        # Figure out what the absolute tolerance should be.  ``self.abs`` is
        # either None or a value specified by the user.
        absolute_tolerance = set_default(self.abs, self.DEFAULT_ABSOLUTE_TOLERANCE)

        if absolute_tolerance < 0:
            raise ValueError("absolute tolerance can't be negative: {}".format(absolute_tolerance))
        if math.isnan(absolute_tolerance):
            raise ValueError("absolute tolerance can't be NaN.")

        # If the user specified an absolute tolerance but not a relative one,
        # just return the absolute tolerance.
        if self.rel is None:
            if self.abs is not None:
                return absolute_tolerance

        # Figure out what the relative tolerance should be.  ``self.rel`` is
        # either None or a value specified by the user.  This is done after
        # we've made sure the user didn't ask for an absolute tolerance only,
        # because we don't want to raise errors about the relative tolerance if
        # we aren't even going to use it.
        relative_tolerance = set_default(self.rel, self.DEFAULT_RELATIVE_TOLERANCE) * abs(self.expected)

        if relative_tolerance < 0:
            raise ValueError("relative tolerance can't be negative: {}".format(absolute_tolerance))
        if math.isnan(relative_tolerance):
            raise ValueError("relative tolerance can't be NaN.")

        # Return the larger of the relative and absolute tolerances.
        return max(relative_tolerance, absolute_tolerance)


class ApproxDecimal(ApproxScalar):
    from decimal import Decimal

    DEFAULT_ABSOLUTE_TOLERANCE = Decimal('1e-12')
    DEFAULT_RELATIVE_TOLERANCE = Decimal('1e-6')


def approx(expected, rel=None, abs=None, nan_ok=False):
    """
    Assert that two numbers (or two sets of numbers) are equal to each other
    within some tolerance.

    Due to the `intricacies of floating-point arithmetic`__, numbers that we
    would intuitively expect to be equal are not always so::

        >>> 0.1 + 0.2 == 0.3
        False

    __ https://docs.python.org/3/tutorial/floatingpoint.html

    This problem is commonly encountered when writing tests, e.g. when making
    sure that floating-point values are what you expect them to be.  One way to
    deal with this problem is to assert that two floating-point numbers are
    equal to within some appropriate tolerance::

        >>> abs((0.1 + 0.2) - 0.3) < 1e-6
        True

    However, comparisons like this are tedious to write and difficult to
    understand.  Furthermore, absolute comparisons like the one above are
    usually discouraged because there's no tolerance that works well for all
    situations.  ``1e-6`` is good for numbers around ``1``, but too small for
    very big numbers and too big for very small ones.  It's better to express
    the tolerance as a fraction of the expected value, but relative comparisons
    like that are even more difficult to write correctly and concisely.

    The ``approx`` class performs floating-point comparisons using a syntax
    that's as intuitive as possible::

        >>> from pytest import approx
        >>> 0.1 + 0.2 == approx(0.3)
        True

    The same syntax also works for sequences of numbers::

        >>> (0.1 + 0.2, 0.2 + 0.4) == approx((0.3, 0.6))
        True

    Dictionary *values*::

        >>> {'a': 0.1 + 0.2, 'b': 0.2 + 0.4} == approx({'a': 0.3, 'b': 0.6})
        True

    ``numpy`` arrays::

        >>> import numpy as np                                                          # doctest: +SKIP
        >>> np.array([0.1, 0.2]) + np.array([0.2, 0.4]) == approx(np.array([0.3, 0.6])) # doctest: +SKIP
        True

    And for a ``numpy`` array against a scalar::

        >>> import numpy as np                                         # doctest: +SKIP
        >>> np.array([0.1, 0.2]) + np.array([0.2, 0.1]) == approx(0.3) # doctest: +SKIP
        True

    By default, ``approx`` considers numbers within a relative tolerance of
    ``1e-6`` (i.e. one part in a million) of its expected value to be equal.
    This treatment would lead to surprising results if the expected value was
    ``0.0``, because nothing but ``0.0`` itself is relatively close to ``0.0``.
    To handle this case less surprisingly, ``approx`` also considers numbers
    within an absolute tolerance of ``1e-12`` of its expected value to be
    equal.  Infinity and NaN are special cases.  Infinity is only considered
    equal to itself, regardless of the relative tolerance.  NaN is not
    considered equal to anything by default, but you can make it be equal to
    itself by setting the ``nan_ok`` argument to True.  (This is meant to
    facilitate comparing arrays that use NaN to mean "no data".)

    Both the relative and absolute tolerances can be changed by passing
    arguments to the ``approx`` constructor::

        >>> 1.0001 == approx(1)
        False
        >>> 1.0001 == approx(1, rel=1e-3)
        True
        >>> 1.0001 == approx(1, abs=1e-3)
        True

    If you specify ``abs`` but not ``rel``, the comparison will not consider
    the relative tolerance at all.  In other words, two numbers that are within
    the default relative tolerance of ``1e-6`` will still be considered unequal
    if they exceed the specified absolute tolerance.  If you specify both
    ``abs`` and ``rel``, the numbers will be considered equal if either
    tolerance is met::

        >>> 1 + 1e-8 == approx(1)
        True
        >>> 1 + 1e-8 == approx(1, abs=1e-12)
        False
        >>> 1 + 1e-8 == approx(1, rel=1e-6, abs=1e-12)
        True

    If you're thinking about using ``approx``, then you might want to know how
    it compares to other good ways of comparing floating-point numbers.  All of
    these algorithms are based on relative and absolute tolerances and should
    agree for the most part, but they do have meaningful differences:

    - ``math.isclose(a, b, rel_tol=1e-9, abs_tol=0.0)``:  True if the relative
      tolerance is met w.r.t. either ``a`` or ``b`` or if the absolute
      tolerance is met.  Because the relative tolerance is calculated w.r.t.
      both ``a`` and ``b``, this test is symmetric (i.e.  neither ``a`` nor
      ``b`` is a "reference value").  You have to specify an absolute tolerance
      if you want to compare to ``0.0`` because there is no tolerance by
      default.  Only available in python>=3.5.  `More information...`__

      __ https://docs.python.org/3/library/math.html#math.isclose

    - ``numpy.isclose(a, b, rtol=1e-5, atol=1e-8)``: True if the difference
      between ``a`` and ``b`` is less that the sum of the relative tolerance
      w.r.t. ``b`` and the absolute tolerance.  Because the relative tolerance
      is only calculated w.r.t. ``b``, this test is asymmetric and you can
      think of ``b`` as the reference value.  Support for comparing sequences
      is provided by ``numpy.allclose``.  `More information...`__

      __ http://docs.scipy.org/doc/numpy-1.10.0/reference/generated/numpy.isclose.html

    - ``unittest.TestCase.assertAlmostEqual(a, b)``: True if ``a`` and ``b``
      are within an absolute tolerance of ``1e-7``.  No relative tolerance is
      considered and the absolute tolerance cannot be changed, so this function
      is not appropriate for very large or very small numbers.  Also, it's only
      available in subclasses of ``unittest.TestCase`` and it's ugly because it
      doesn't follow PEP8.  `More information...`__

      __ https://docs.python.org/3/library/unittest.html#unittest.TestCase.assertAlmostEqual

    - ``a == pytest.approx(b, rel=1e-6, abs=1e-12)``: True if the relative
      tolerance is met w.r.t. ``b`` or if the absolute tolerance is met.
      Because the relative tolerance is only calculated w.r.t. ``b``, this test
      is asymmetric and you can think of ``b`` as the reference value.  In the
      special case that you explicitly specify an absolute tolerance but not a
      relative tolerance, only the absolute tolerance is considered.

    .. warning::

       .. versionchanged:: 3.2

       In order to avoid inconsistent behavior, ``TypeError`` is
       raised for ``>``, ``>=``, ``<`` and ``<=`` comparisons.
       The example below illustrates the problem::

           assert approx(0.1) > 0.1 + 1e-10  # calls approx(0.1).__gt__(0.1 + 1e-10)
           assert 0.1 + 1e-10 > approx(0.1)  # calls approx(0.1).__lt__(0.1 + 1e-10)

       In the second example one expects ``approx(0.1).__le__(0.1 + 1e-10)``
       to be called. But instead, ``approx(0.1).__lt__(0.1 + 1e-10)`` is used to
       comparison. This is because the call hierarchy of rich comparisons
       follows a fixed behavior. `More information...`__

       __ https://docs.python.org/3/reference/datamodel.html#object.__ge__
    """

    from collections import Mapping, Sequence
    from _pytest.compat import STRING_TYPES as String
    from decimal import Decimal

    # Delegate the comparison to a class that knows how to deal with the type
    # of the expected value (e.g. int, float, list, dict, numpy.array, etc).
    #
    # This architecture is really driven by the need to support numpy arrays.
    # The only way to override `==` for arrays without requiring that approx be
    # the left operand is to inherit the approx object from `numpy.ndarray`.
    # But that can't be a general solution, because it requires (1) numpy to be
    # installed and (2) the expected value to be a numpy array.  So the general
    # solution is to delegate each type of expected value to a different class.
    #
    # This has the advantage that it made it easy to support mapping types
    # (i.e. dict).  The old code accepted mapping types, but would only compare
    # their keys, which is probably not what most people would expect.

    if _is_numpy_array(expected):
        cls = ApproxNumpy
    elif isinstance(expected, Mapping):
        cls = ApproxMapping
    elif isinstance(expected, Sequence) and not isinstance(expected, String):
        cls = ApproxSequence
    elif isinstance(expected, Decimal):
        cls = ApproxDecimal
    else:
        cls = ApproxScalar

    return cls(expected, rel, abs, nan_ok)


def _is_numpy_array(obj):
    """
    Return true if the given object is a numpy array.  Make a special effort to
    avoid importing numpy unless it's really necessary.
    """
    import inspect

    for cls in inspect.getmro(type(obj)):
        if cls.__module__ == 'numpy':
            try:
                import numpy as np
                return isinstance(obj, np.ndarray)
            except ImportError:
                pass

    return False


# builtin pytest.raises helper

def raises(expected_exception, *args, **kwargs):
    """
    Assert that a code block/function call raises ``expected_exception``
    and raise a failure exception otherwise.

    :arg message: if specified, provides a custom failure message if the
        exception is not raised
    :arg match: if specified, asserts that the exception matches a text or regex

    This helper produces a ``ExceptionInfo()`` object (see below).

    You may use this function as a context manager::

        >>> with raises(ZeroDivisionError):
        ...    1/0

    .. versionchanged:: 2.10

    In the context manager form you may use the keyword argument
    ``message`` to specify a custom failure message::

        >>> with raises(ZeroDivisionError, message="Expecting ZeroDivisionError"):
        ...    pass
        Traceback (most recent call last):
          ...
        Failed: Expecting ZeroDivisionError

    .. note::

       When using ``pytest.raises`` as a context manager, it's worthwhile to
       note that normal context manager rules apply and that the exception
       raised *must* be the final line in the scope of the context manager.
       Lines of code after that, within the scope of the context manager will
       not be executed. For example::

           >>> value = 15
           >>> with raises(ValueError) as exc_info:
           ...     if value > 10:
           ...         raise ValueError("value must be <= 10")
           ...     assert exc_info.type == ValueError  # this will not execute

       Instead, the following approach must be taken (note the difference in
       scope)::

           >>> with raises(ValueError) as exc_info:
           ...     if value > 10:
           ...         raise ValueError("value must be <= 10")
           ...
           >>> assert exc_info.type == ValueError


    Since version ``3.1`` you can use the keyword argument ``match`` to assert that the
    exception matches a text or regex::

        >>> with raises(ValueError, match='must be 0 or None'):
        ...     raise ValueError("value must be 0 or None")

        >>> with raises(ValueError, match=r'must be \d+$'):
        ...     raise ValueError("value must be 42")

    **Legacy forms**

    The forms below are fully supported but are discouraged for new code because the
    context manager form is regarded as more readable and less error-prone.

    It is possible to specify a callable by passing a to-be-called lambda::

        >>> raises(ZeroDivisionError, lambda: 1/0)
        <ExceptionInfo ...>

    or you can specify an arbitrary callable with arguments::

        >>> def f(x): return 1/x
        ...
        >>> raises(ZeroDivisionError, f, 0)
        <ExceptionInfo ...>
        >>> raises(ZeroDivisionError, f, x=0)
        <ExceptionInfo ...>

    It is also possible to pass a string to be evaluated at runtime::

        >>> raises(ZeroDivisionError, "f(0)")
        <ExceptionInfo ...>

    The string will be evaluated using the same ``locals()`` and ``globals()``
    at the moment of the ``raises`` call.

    .. currentmodule:: _pytest._code

    Consult the API of ``excinfo`` objects: :class:`ExceptionInfo`.

    .. note::
        Similar to caught exception objects in Python, explicitly clearing
        local references to returned ``ExceptionInfo`` objects can
        help the Python interpreter speed up its garbage collection.

        Clearing those references breaks a reference cycle
        (``ExceptionInfo`` --> caught exception --> frame stack raising
        the exception --> current frame stack --> local variables -->
        ``ExceptionInfo``) which makes Python keep all objects referenced
        from that cycle (including all local variables in the current
        frame) alive until the next cyclic garbage collection run. See the
        official Python ``try`` statement documentation for more detailed
        information.

    """
    __tracebackhide__ = True
    base_type = (type, text_type, binary_type)
    for exc in filterfalse(isclass, always_iterable(expected_exception, base_type)):
        msg = ("exceptions must be old-style classes or"
               " derived from BaseException, not %s")
        raise TypeError(msg % type(exc))

    message = "DID NOT RAISE {0}".format(expected_exception)
    match_expr = None

    if not args:
        if "message" in kwargs:
            message = kwargs.pop("message")
        if "match" in kwargs:
            match_expr = kwargs.pop("match")
        if kwargs:
            msg = 'Unexpected keyword arguments passed to pytest.raises: '
            msg += ', '.join(kwargs.keys())
            raise TypeError(msg)
        return RaisesContext(expected_exception, message, match_expr)
    elif isinstance(args[0], str):
        code, = args
        assert isinstance(code, str)
        frame = sys._getframe(1)
        loc = frame.f_locals.copy()
        loc.update(kwargs)
        # print "raises frame scope: %r" % frame.f_locals
        try:
            code = _pytest._code.Source(code).compile()
            py.builtin.exec_(code, frame.f_globals, loc)
            # XXX didn'T mean f_globals == f_locals something special?
            #     this is destroyed here ...
        except expected_exception:
            return _pytest._code.ExceptionInfo()
    else:
        func = args[0]
        try:
            func(*args[1:], **kwargs)
        except expected_exception:
            return _pytest._code.ExceptionInfo()
    fail(message)


raises.Exception = fail.Exception


class RaisesContext(object):
    def __init__(self, expected_exception, message, match_expr):
        self.expected_exception = expected_exception
        self.message = message
        self.match_expr = match_expr
        self.excinfo = None

    def __enter__(self):
        self.excinfo = object.__new__(_pytest._code.ExceptionInfo)
        return self.excinfo

    def __exit__(self, *tp):
        __tracebackhide__ = True
        if tp[0] is None:
            fail(self.message)
        self.excinfo.__init__(tp)
        suppress_exception = issubclass(self.excinfo.type, self.expected_exception)
        if sys.version_info[0] == 2 and suppress_exception:
            sys.exc_clear()
        if self.match_expr and suppress_exception:
            self.excinfo.match(self.match_expr)
        return suppress_exception
