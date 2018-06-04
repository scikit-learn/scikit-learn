from __future__ import absolute_import, division, print_function

from functools import wraps
from distutils.version import LooseVersion

import numpy as np

from ..base import normalize_token
from .core import (concatenate_lookup, tensordot_lookup, map_blocks,
                   asanyarray, atop)


if LooseVersion(np.__version__) < '1.11.2':
    raise ImportError("dask.array.ma requires numpy >= 1.11.2")


@normalize_token.register(np.ma.masked_array)
def normalize_masked_array(x):
    data = normalize_token(x.data)
    mask = normalize_token(x.mask)
    fill_value = normalize_token(x.fill_value)
    return (data, mask, fill_value)


@concatenate_lookup.register(np.ma.masked_array)
def _concatenate(arrays, axis=0):
    out = np.ma.concatenate(arrays, axis=axis)
    fill_values = [i.fill_value for i in arrays if hasattr(i, 'fill_value')]
    if any(isinstance(f, np.ndarray) for f in fill_values):
        raise ValueError("Dask doesn't support masked array's with "
                         "non-scalar `fill_value`s")
    if fill_values:
        # If all the fill_values are the same copy over the fill value
        fill_values = np.unique(fill_values)
        if len(fill_values) == 1:
            out.fill_value = fill_values[0]
    return out


@tensordot_lookup.register(np.ma.masked_array)
def _tensordot(a, b, axes=2):
    # Much of this is stolen from numpy/core/numeric.py::tensordot
    # Please see license at https://github.com/numpy/numpy/blob/master/LICENSE.txt
    try:
        iter(axes)
    except TypeError:
        axes_a = list(range(-axes, 0))
        axes_b = list(range(0, axes))
    else:
        axes_a, axes_b = axes
    try:
        na = len(axes_a)
        axes_a = list(axes_a)
    except TypeError:
        axes_a = [axes_a]
        na = 1
    try:
        nb = len(axes_b)
        axes_b = list(axes_b)
    except TypeError:
        axes_b = [axes_b]
        nb = 1

    # a, b = asarray(a), asarray(b)  # <--- modified
    as_ = a.shape
    nda = a.ndim
    bs = b.shape
    ndb = b.ndim
    equal = True
    if na != nb:
        equal = False
    else:
        for k in range(na):
            if as_[axes_a[k]] != bs[axes_b[k]]:
                equal = False
                break
            if axes_a[k] < 0:
                axes_a[k] += nda
            if axes_b[k] < 0:
                axes_b[k] += ndb
    if not equal:
        raise ValueError("shape-mismatch for sum")

    # Move the axes to sum over to the end of "a"
    # and to the front of "b"
    notin = [k for k in range(nda) if k not in axes_a]
    newaxes_a = notin + axes_a
    N2 = 1
    for axis in axes_a:
        N2 *= as_[axis]
    newshape_a = (-1, N2)
    olda = [as_[axis] for axis in notin]

    notin = [k for k in range(ndb) if k not in axes_b]
    newaxes_b = axes_b + notin
    N2 = 1
    for axis in axes_b:
        N2 *= bs[axis]
    newshape_b = (N2, -1)
    oldb = [bs[axis] for axis in notin]

    at = a.transpose(newaxes_a).reshape(newshape_a)
    bt = b.transpose(newaxes_b).reshape(newshape_b)
    res = np.ma.dot(at, bt)
    return res.reshape(olda + oldb)


@wraps(np.ma.filled)
def filled(a, fill_value=None):
    a = asanyarray(a)
    return a.map_blocks(np.ma.filled, fill_value=fill_value)


def _wrap_masked(f):

    @wraps(f)
    def _(a, value):
        a = asanyarray(a)
        value = asanyarray(value)
        ainds = tuple(range(a.ndim))[::-1]
        vinds = tuple(range(value.ndim))[::-1]
        oinds = max(ainds, vinds, key=len)
        return atop(f, oinds, a, ainds, value, vinds, dtype=a.dtype)

    return _


masked_greater = _wrap_masked(np.ma.masked_greater)
masked_greater_equal = _wrap_masked(np.ma.masked_greater_equal)
masked_less = _wrap_masked(np.ma.masked_less)
masked_less_equal = _wrap_masked(np.ma.masked_less_equal)
masked_not_equal = _wrap_masked(np.ma.masked_not_equal)


@wraps(np.ma.masked_equal)
def masked_equal(a, value):
    a = asanyarray(a)
    if getattr(value, 'shape', ()):
        raise ValueError("da.ma.masked_equal doesn't support array `value`s")
    inds = tuple(range(a.ndim))
    return atop(np.ma.masked_equal, inds, a, inds, value, (), dtype=a.dtype)


@wraps(np.ma.masked_invalid)
def masked_invalid(a):
    return asanyarray(a).map_blocks(np.ma.masked_invalid)


@wraps(np.ma.masked_inside)
def masked_inside(x, v1, v2):
    x = asanyarray(x)
    return x.map_blocks(np.ma.masked_inside, v1, v2)


@wraps(np.ma.masked_outside)
def masked_outside(x, v1, v2):
    x = asanyarray(x)
    return x.map_blocks(np.ma.masked_outside, v1, v2)


@wraps(np.ma.masked_where)
def masked_where(condition, a):
    cshape = getattr(condition, 'shape', ())
    if cshape and cshape != a.shape:
        raise IndexError("Inconsistant shape between the condition and the "
                         "input (got %s and %s)" % (cshape, a.shape))
    condition = asanyarray(condition)
    a = asanyarray(a)
    ainds = tuple(range(a.ndim))
    cinds = tuple(range(condition.ndim))
    return atop(np.ma.masked_where, ainds, condition, cinds, a, ainds,
                dtype=a.dtype)


@wraps(np.ma.masked_values)
def masked_values(x, value, rtol=1e-05, atol=1e-08, shrink=True):
    x = asanyarray(x)
    if getattr(value, 'shape', ()):
        raise ValueError("da.ma.masked_values doesn't support array `value`s")
    return map_blocks(np.ma.masked_values, x, value, rtol=rtol,
                      atol=atol, shrink=shrink)


@wraps(np.ma.fix_invalid)
def fix_invalid(a, fill_value=None):
    a = asanyarray(a)
    return a.map_blocks(np.ma.fix_invalid, fill_value=fill_value)


@wraps(np.ma.getdata)
def getdata(a):
    a = asanyarray(a)
    return a.map_blocks(np.ma.getdata)


@wraps(np.ma.getmaskarray)
def getmaskarray(a):
    a = asanyarray(a)
    return a.map_blocks(np.ma.getmaskarray)


def _masked_array(data, mask=np.ma.nomask, **kwargs):
    dtype = kwargs.pop('masked_dtype', None)
    return np.ma.masked_array(data, mask=mask, dtype=dtype, **kwargs)


@wraps(np.ma.masked_array)
def masked_array(data, mask=np.ma.nomask, fill_value=None,
                 **kwargs):
    data = asanyarray(data)
    inds = tuple(range(data.ndim))
    arginds = [inds, data, inds]

    if getattr(fill_value, 'shape', ()):
        raise ValueError("non-scalar fill_value not supported")
    kwargs['fill_value'] = fill_value

    if mask is not np.ma.nomask:
        mask = asanyarray(mask)
        if mask.size == 1:
            mask = mask.reshape((1,) * data.ndim)
        elif data.shape != mask.shape:
            raise np.ma.MaskError("Mask and data not compatible: data shape "
                                  "is %s, and mask shape is "
                                  "%s." % (repr(data.shape), repr(mask.shape)))
        arginds.extend([mask, inds])

    if 'dtype' in kwargs:
        kwargs['masked_dtype'] = kwargs['dtype']
    else:
        kwargs['dtype'] = data.dtype

    return atop(_masked_array, *arginds, **kwargs)


def _set_fill_value(x, fill_value):
    if isinstance(x, np.ma.masked_array):
        x = x.copy()
        np.ma.set_fill_value(x, fill_value=fill_value)
    return x


@wraps(np.ma.set_fill_value)
def set_fill_value(a, fill_value):
    a = asanyarray(a)
    if getattr(fill_value, 'shape', ()):
        raise ValueError("da.ma.set_fill_value doesn't support array `value`s")
    fill_value = np.ma.core._check_fill_value(fill_value, a.dtype)
    res = a.map_blocks(_set_fill_value, fill_value)
    a.dask = res.dask
    a.name = res.name
