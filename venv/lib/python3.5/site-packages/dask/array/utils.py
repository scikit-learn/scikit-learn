from __future__ import absolute_import, division, print_function

import difflib
import functools
import math
import os

import numpy as np
from toolz import frequencies, concat

from .core import Array
from ..local import get_sync
from ..sharedict import ShareDict


def allclose(a, b, equal_nan=False, **kwargs):
    if getattr(a, 'dtype', None) != 'O':
        return np.allclose(a, b, equal_nan=equal_nan, **kwargs)
    if equal_nan:
        return (a.shape == b.shape and
                all(np.isnan(b) if np.isnan(a) else a == b
                    for (a, b) in zip(a.flat, b.flat)))
    return (a == b).all()


def same_keys(a, b):
    def key(k):
        if isinstance(k, str):
            return (k, -1, -1, -1)
        else:
            return k
    return sorted(a.dask, key=key) == sorted(b.dask, key=key)


def _not_empty(x):
    return x.shape and 0 not in x.shape


def _check_dsk(dsk):
    """ Check that graph is well named and non-overlapping """
    if not isinstance(dsk, ShareDict):
        return

    assert all(isinstance(k, (tuple, str)) for k in dsk.dicts)
    freqs = frequencies(concat(dsk.dicts.values()))
    non_one = {k: v for k, v in freqs.items() if v != 1}
    assert not non_one, non_one


def assert_eq_shape(a, b, check_nan=True):
    for aa, bb in zip(a, b):
        if math.isnan(aa) or math.isnan(bb):
            if check_nan:
                assert math.isnan(aa) == math.isnan(bb)
        else:
            assert aa == bb


def assert_eq(a, b, check_shape=True, **kwargs):
    a_original = a
    b_original = b
    if isinstance(a, Array):
        assert a.dtype is not None
        adt = a.dtype
        _check_dsk(a.dask)
        a = a.compute(get=get_sync)
        if hasattr(a, 'todense'):
            a = a.todense()
        if not hasattr(a, 'dtype'):
            a = np.array(a, dtype='O')
        if _not_empty(a):
            assert a.dtype == a_original.dtype
        if check_shape:
            assert_eq_shape(a_original.shape, a.shape, check_nan=False)
    else:
        if not hasattr(a, 'dtype'):
            a = np.array(a, dtype='O')
        adt = getattr(a, 'dtype', None)

    if isinstance(b, Array):
        assert b.dtype is not None
        bdt = b.dtype
        _check_dsk(b.dask)
        b = b.compute(get=get_sync)
        if not hasattr(b, 'dtype'):
            b = np.array(b, dtype='O')
        if hasattr(b, 'todense'):
            b = b.todense()
        if _not_empty(b):
            assert b.dtype == b_original.dtype
        if check_shape:
            assert_eq_shape(b_original.shape, b.shape, check_nan=False)
    else:
        if not hasattr(b, 'dtype'):
            b = np.array(b, dtype='O')
        bdt = getattr(b, 'dtype', None)

    if str(adt) != str(bdt):
        diff = difflib.ndiff(str(adt).splitlines(), str(bdt).splitlines())
        raise AssertionError('string repr are different' + os.linesep +
                             os.linesep.join(diff))

    try:
        assert a.shape == b.shape
        assert allclose(a, b, **kwargs)
        return True
    except TypeError:
        pass

    c = a == b

    if isinstance(c, np.ndarray):
        assert c.all()
    else:
        assert c

    return True


def safe_wraps(wrapped, assigned=functools.WRAPPER_ASSIGNMENTS):
    """Like functools.wraps, but safe to use even if wrapped is not a function.

    Only needed on Python 2.
    """
    if all(hasattr(wrapped, attr) for attr in assigned):
        return functools.wraps(wrapped, assigned=assigned)
    else:
        return lambda x: x
