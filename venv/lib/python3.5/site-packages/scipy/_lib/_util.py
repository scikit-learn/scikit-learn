from __future__ import division, print_function, absolute_import

import functools
import operator
import sys
import warnings
import numbers
from collections import namedtuple
import inspect

import numpy as np


def _valarray(shape, value=np.nan, typecode=None):
    """Return an array of all value.
    """

    out = np.ones(shape, dtype=bool) * value
    if typecode is not None:
        out = out.astype(typecode)
    if not isinstance(out, np.ndarray):
        out = np.asarray(out)
    return out


def _lazywhere(cond, arrays, f, fillvalue=None, f2=None):
    """
    np.where(cond, x, fillvalue) always evaluates x even where cond is False.
    This one only evaluates f(arr1[cond], arr2[cond], ...).
    For example,
    >>> a, b = np.array([1, 2, 3, 4]), np.array([5, 6, 7, 8])
    >>> def f(a, b):
        return a*b
    >>> _lazywhere(a > 2, (a, b), f, np.nan)
    array([ nan,  nan,  21.,  32.])

    Notice it assumes that all `arrays` are of the same shape, or can be
    broadcasted together.

    """
    if fillvalue is None:
        if f2 is None:
            raise ValueError("One of (fillvalue, f2) must be given.")
        else:
            fillvalue = np.nan
    else:
        if f2 is not None:
            raise ValueError("Only one of (fillvalue, f2) can be given.")

    arrays = np.broadcast_arrays(*arrays)
    temp = tuple(np.extract(cond, arr) for arr in arrays)
    tcode = np.mintypecode([a.dtype.char for a in arrays])
    out = _valarray(np.shape(arrays[0]), value=fillvalue, typecode=tcode)
    np.place(out, cond, f(*temp))
    if f2 is not None:
        temp = tuple(np.extract(~cond, arr) for arr in arrays)
        np.place(out, ~cond, f2(*temp))

    return out


def _lazyselect(condlist, choicelist, arrays, default=0):
    """
    Mimic `np.select(condlist, choicelist)`.

    Notice it assumes that all `arrays` are of the same shape, or can be
    broadcasted together.

    All functions in `choicelist` must accept array arguments in the order
    given in `arrays` and must return an array of the same shape as broadcasted
    `arrays`.

    Examples
    --------
    >>> x = np.arange(6)
    >>> np.select([x <3, x > 3], [x**2, x**3], default=0)
    array([  0,   1,   4,   0,  64, 125])

    >>> _lazyselect([x < 3, x > 3], [lambda x: x**2, lambda x: x**3], (x,))
    array([   0.,    1.,    4.,   0.,   64.,  125.])

    >>> a = -np.ones_like(x)
    >>> _lazyselect([x < 3, x > 3],
    ...             [lambda x, a: x**2, lambda x, a: a * x**3],
    ...             (x, a), default=np.nan)
    array([   0.,    1.,    4.,   nan,  -64., -125.])

    """
    arrays = np.broadcast_arrays(*arrays)
    tcode = np.mintypecode([a.dtype.char for a in arrays])
    out = _valarray(np.shape(arrays[0]), value=default, typecode=tcode)
    for index in range(len(condlist)):
        func, cond = choicelist[index], condlist[index]
        if np.all(cond is False):
            continue
        cond, _ = np.broadcast_arrays(cond, arrays[0])
        temp = tuple(np.extract(cond, arr) for arr in arrays)
        np.place(out, cond, func(*temp))
    return out


def _aligned_zeros(shape, dtype=float, order="C", align=None):
    """Allocate a new ndarray with aligned memory.

    Primary use case for this currently is working around a f2py issue
    in Numpy 1.9.1, where dtype.alignment is such that np.zeros() does
    not necessarily create arrays aligned up to it.

    """
    dtype = np.dtype(dtype)
    if align is None:
        align = dtype.alignment
    if not hasattr(shape, '__len__'):
        shape = (shape,)
    size = functools.reduce(operator.mul, shape) * dtype.itemsize
    buf = np.empty(size + align + 1, np.uint8)
    offset = buf.__array_interface__['data'][0] % align
    if offset != 0:
        offset = align - offset
    # Note: slices producing 0-size arrays do not necessarily change
    # data pointer --- so we use and allocate size+1
    buf = buf[offset:offset+size+1][:-1]
    data = np.ndarray(shape, dtype, buf, order=order)
    data.fill(0)
    return data


def _prune_array(array):
    """Return an array equivalent to the input array. If the input
    array is a view of a much larger array, copy its contents to a
    newly allocated array. Otherwise, return the input unchanged.
    """
    if array.base is not None and array.size < array.base.size // 2:
        return array.copy()
    return array


class DeprecatedImport(object):
    """
    Deprecated import, with redirection + warning.

    Examples
    --------
    Suppose you previously had in some module::

        from foo import spam

    If this has to be deprecated, do::

        spam = DeprecatedImport("foo.spam", "baz")

    to redirect users to use "baz" module instead.

    """

    def __init__(self, old_module_name, new_module_name):
        self._old_name = old_module_name
        self._new_name = new_module_name
        __import__(self._new_name)
        self._mod = sys.modules[self._new_name]

    def __dir__(self):
        return dir(self._mod)

    def __getattr__(self, name):
        warnings.warn("Module %s is deprecated, use %s instead"
                      % (self._old_name, self._new_name),
                      DeprecationWarning)
        return getattr(self._mod, name)


# copy-pasted from scikit-learn utils/validation.py
def check_random_state(seed):
    """Turn seed into a np.random.RandomState instance

    If seed is None (or np.random), return the RandomState singleton used
    by np.random.
    If seed is an int, return a new RandomState instance seeded with seed.
    If seed is already a RandomState instance, return it.
    Otherwise raise ValueError.
    """
    if seed is None or seed is np.random:
        return np.random.mtrand._rand
    if isinstance(seed, (numbers.Integral, np.integer)):
        return np.random.RandomState(seed)
    if isinstance(seed, np.random.RandomState):
        return seed
    raise ValueError('%r cannot be used to seed a numpy.random.RandomState'
                     ' instance' % seed)


def _asarray_validated(a, check_finite=True,
                       sparse_ok=False, objects_ok=False, mask_ok=False,
                       as_inexact=False):
    """
    Helper function for scipy argument validation.

    Many scipy linear algebra functions do support arbitrary array-like
    input arguments.  Examples of commonly unsupported inputs include
    matrices containing inf/nan, sparse matrix representations, and
    matrices with complicated elements.

    Parameters
    ----------
    a : array_like
        The array-like input.
    check_finite : bool, optional
        Whether to check that the input matrices contain only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.
        Default: True
    sparse_ok : bool, optional
        True if scipy sparse matrices are allowed.
    objects_ok : bool, optional
        True if arrays with dype('O') are allowed.
    mask_ok : bool, optional
        True if masked arrays are allowed.
    as_inexact : bool, optional
        True to convert the input array to a np.inexact dtype.

    Returns
    -------
    ret : ndarray
        The converted validated array.

    """
    if not sparse_ok:
        import scipy.sparse
        if scipy.sparse.issparse(a):
            msg = ('Sparse matrices are not supported by this function. '
                   'Perhaps one of the scipy.sparse.linalg functions '
                   'would work instead.')
            raise ValueError(msg)
    if not mask_ok:
        if np.ma.isMaskedArray(a):
            raise ValueError('masked arrays are not supported')
    toarray = np.asarray_chkfinite if check_finite else np.asarray
    a = toarray(a)
    if not objects_ok:
        if a.dtype is np.dtype('O'):
            raise ValueError('object arrays are not supported')
    if as_inexact:
        if not np.issubdtype(a.dtype, np.inexact):
            a = toarray(a, dtype=np.float_)
    return a


# Add a replacement for inspect.getargspec() which is deprecated in python 3.5
# The version below is borrowed from Django,
# https://github.com/django/django/pull/4846

# Note an inconsistency between inspect.getargspec(func) and
# inspect.signature(func). If `func` is a bound method, the latter does *not*
# list `self` as a first argument, while the former *does*.
# Hence cook up a common ground replacement: `getargspec_no_self` which
# mimics `inspect.getargspec` but does not list `self`.
#
# This way, the caller code does not need to know whether it uses a legacy
# .getargspec or bright and shiny .signature.

try:
    # is it python 3.3 or higher?
    inspect.signature

    # Apparently, yes. Wrap inspect.signature

    ArgSpec = namedtuple('ArgSpec', ['args', 'varargs', 'keywords', 'defaults'])

    def getargspec_no_self(func):
        """inspect.getargspec replacement using inspect.signature.

        inspect.getargspec is deprecated in python 3. This is a replacement
        based on the (new in python 3.3) `inspect.signature`.

        Parameters
        ----------
        func : callable
            A callable to inspect

        Returns
        -------
        argspec : ArgSpec(args, varargs, varkw, defaults)
            This is similar to the result of inspect.getargspec(func) under
            python 2.x.
            NOTE: if the first argument of `func` is self, it is *not*, I repeat
            *not* included in argspec.args.
            This is done for consistency between inspect.getargspec() under
            python 2.x, and inspect.signature() under python 3.x.
        """
        sig = inspect.signature(func)
        args = [
            p.name for p in sig.parameters.values()
            if p.kind == inspect.Parameter.POSITIONAL_OR_KEYWORD
        ]
        varargs = [
            p.name for p in sig.parameters.values()
            if p.kind == inspect.Parameter.VAR_POSITIONAL
        ]
        varargs = varargs[0] if varargs else None
        varkw = [
            p.name for p in sig.parameters.values()
            if p.kind == inspect.Parameter.VAR_KEYWORD
        ]
        varkw = varkw[0] if varkw else None
        defaults = [
            p.default for p in sig.parameters.values()
            if (p.kind == inspect.Parameter.POSITIONAL_OR_KEYWORD and
               p.default is not p.empty)
        ] or None
        return ArgSpec(args, varargs, varkw, defaults)

except AttributeError:
    # python 2.x
    def getargspec_no_self(func):
        """inspect.getargspec replacement for compatibility with python 3.x.

        inspect.getargspec is deprecated in python 3. This wraps it, and
        *removes* `self` from the argument list of `func`, if present.
        This is done for forward compatibility with python 3.

        Parameters
        ----------
        func : callable
            A callable to inspect

        Returns
        -------
        argspec : ArgSpec(args, varargs, varkw, defaults)
            This is similar to the result of inspect.getargspec(func) under
            python 2.x.
            NOTE: if the first argument of `func` is self, it is *not*, I repeat
            *not* included in argspec.args.
            This is done for consistency between inspect.getargspec() under
            python 2.x, and inspect.signature() under python 3.x.
        """
        argspec = inspect.getargspec(func)
        if argspec.args[0] == 'self':
            argspec.args.pop(0)
        return argspec
