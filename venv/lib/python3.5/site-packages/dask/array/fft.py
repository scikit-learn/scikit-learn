from __future__ import absolute_import, division, print_function

import collections
from functools import wraps
import inspect

import numpy as np

try:
    import scipy
    import scipy.fftpack
except ImportError:
    scipy = None

from .core import concatenate as _concatenate
from .creation import arange as _arange


chunk_error = ("Dask array only supports taking an FFT along an axis that \n"
               "has a single chunk. An FFT operation was tried on axis %s \n"
               "which has chunks %s. To change the array's chunks use "
               "dask.Array.rechunk.")

fft_preamble = """
    Wrapping of %s

    The axis along which the FFT is applied must have a one chunk. To change
    the array's chunking use dask.Array.rechunk.

    The %s docstring follows below:

    """


def _fft_out_chunks(a, s, axes):
    """ For computing the output chunks of [i]fft*"""
    if s is None:
        return a.chunks
    chunks = list(a.chunks)
    for i, axis in enumerate(axes):
        chunks[axis] = (s[i],)
    return chunks


def _rfft_out_chunks(a, s, axes):
    """ For computing the output chunks of rfft*"""
    if s is None:
        s = [a.chunks[axis][0] for axis in axes]
    s = list(s)
    s[-1] = s[-1] // 2 + 1
    chunks = list(a.chunks)
    for i, axis in enumerate(axes):
        chunks[axis] = (s[i],)
    return chunks


def _irfft_out_chunks(a, s, axes):
    """ For computing the output chunks of irfft*"""
    if s is None:
        s = [a.chunks[axis][0] for axis in axes]
        s[-1] = 2 * (s[-1] - 1)
    chunks = list(a.chunks)
    for i, axis in enumerate(axes):
        chunks[axis] = (s[i],)
    return chunks


def _hfft_out_chunks(a, s, axes):
    assert len(axes) == 1

    axis = axes[0]

    if s is None:
        s = [2 * (a.chunks[axis][0] - 1)]

    n = s[0]

    chunks = list(a.chunks)
    chunks[axis] = (n,)
    return chunks


def _ihfft_out_chunks(a, s, axes):
    assert len(axes) == 1

    axis = axes[0]

    if s is None:
        s = [a.chunks[axis][0]]
    else:
        assert len(s) == 1

    n = s[0]

    chunks = list(a.chunks)
    if n % 2 == 0:
        m = (n // 2) + 1
    else:
        m = (n + 1) // 2
    chunks[axis] = (m,)
    return chunks


_out_chunk_fns = {'fft': _fft_out_chunks,
                  'ifft': _fft_out_chunks,
                  'rfft': _rfft_out_chunks,
                  'irfft': _irfft_out_chunks,
                  'hfft': _hfft_out_chunks,
                  'ihfft': _ihfft_out_chunks}


def fft_wrap(fft_func, kind=None, dtype=None):
    """ Wrap 1D, 2D, and ND real and complex FFT functions

    Takes a function that behaves like ``numpy.fft`` functions and
    a specified kind to match it to that are named after the functions
    in the ``numpy.fft`` API.

    Supported kinds include:

        * fft
        * fft2
        * fftn
        * ifft
        * ifft2
        * ifftn
        * rfft
        * rfft2
        * rfftn
        * irfft
        * irfft2
        * irfftn
        * hfft
        * ihfft

    Examples
    --------
    >>> parallel_fft = fft_wrap(np.fft.fft)
    >>> parallel_ifft = fft_wrap(np.fft.ifft)
    """
    if scipy is not None:
        if fft_func is scipy.fftpack.rfft:
            raise ValueError("SciPy's `rfft` doesn't match the NumPy API.")
        elif fft_func is scipy.fftpack.irfft:
            raise ValueError("SciPy's `irfft` doesn't match the NumPy API.")

    if kind is None:
        kind = fft_func.__name__
    try:
        out_chunk_fn = _out_chunk_fns[kind.rstrip("2n")]
    except KeyError:
        raise ValueError("Given unknown `kind` %s." % kind)

    def func(a, s=None, axes=None):
        if axes is None:
            if kind.endswith('2'):
                axes = (-2, -1)
            elif kind.endswith('n'):
                if s is None:
                    axes = tuple(range(a.ndim))
                else:
                    axes = tuple(range(len(s)))
            else:
                axes = (-1,)
        else:
            if len(set(axes)) < len(axes):
                raise ValueError("Duplicate axes not allowed.")

        _dtype = dtype
        if _dtype is None:
            _dtype = fft_func(np.ones(len(axes) * (8,),
                                      dtype=a.dtype)).dtype

        for each_axis in axes:
            if len(a.chunks[each_axis]) != 1:
                raise ValueError(chunk_error % (each_axis, a.chunks[each_axis]))

        chunks = out_chunk_fn(a, s, axes)

        args = (s, axes)
        if kind.endswith('fft'):
            axis = None if axes is None else axes[0]
            n = None if s is None else s[0]
            args = (n, axis)

        return a.map_blocks(fft_func, *args, dtype=_dtype,
                            chunks=chunks)

    if kind.endswith('fft'):
        _func = func

        def func(a, n=None, axis=None):
            s = None
            if n is not None:
                s = (n,)

            axes = None
            if axis is not None:
                axes = (axis,)

            return _func(a, s, axes)

    func_mod = inspect.getmodule(fft_func)
    func_name = fft_func.__name__
    func_fullname = func_mod.__name__ + "." + func_name
    if fft_func.__doc__ is not None:
        func.__doc__ = (fft_preamble % (2 * (func_fullname,)))
        func.__doc__ += fft_func.__doc__
    func.__name__ = func_name
    return func


fft = fft_wrap(np.fft.fft, dtype=np.complex_)
fft2 = fft_wrap(np.fft.fft2, dtype=np.complex_)
fftn = fft_wrap(np.fft.fftn, dtype=np.complex_)
ifft = fft_wrap(np.fft.ifft, dtype=np.complex_)
ifft2 = fft_wrap(np.fft.ifft2, dtype=np.complex_)
ifftn = fft_wrap(np.fft.ifftn, dtype=np.complex_)
rfft = fft_wrap(np.fft.rfft, dtype=np.complex_)
rfft2 = fft_wrap(np.fft.rfft2, dtype=np.complex_)
rfftn = fft_wrap(np.fft.rfftn, dtype=np.complex_)
irfft = fft_wrap(np.fft.irfft, dtype=np.float_)
irfft2 = fft_wrap(np.fft.irfft2, dtype=np.float_)
irfftn = fft_wrap(np.fft.irfftn, dtype=np.float_)
hfft = fft_wrap(np.fft.hfft, dtype=np.float_)
ihfft = fft_wrap(np.fft.ihfft, dtype=np.complex_)


def _fftfreq_block(i, n, d):
    r = i.copy()
    r[i >= (n + 1) // 2] -= n
    r /= n * d
    return r


@wraps(np.fft.fftfreq)
def fftfreq(n, d=1.0, chunks=None):
    n = int(n)
    d = float(d)

    r = _arange(n, dtype=float, chunks=chunks)

    return r.map_blocks(_fftfreq_block, dtype=float, n=n, d=d)


@wraps(np.fft.rfftfreq)
def rfftfreq(n, d=1.0, chunks=None):
    n = int(n)
    d = float(d)

    r = _arange(n // 2 + 1, dtype=float, chunks=chunks)
    r /= n * d

    return r


def _fftshift_helper(x, axes=None, inverse=False):
    if axes is None:
        axes = list(range(x.ndim))
    elif not isinstance(axes, collections.Sequence):
        axes = (axes,)

    y = x
    for i in axes:
        n = y.shape[i]
        n_2 = (n + int(inverse is False)) // 2

        l = y.ndim * [slice(None)]
        l[i] = slice(None, n_2)
        l = tuple(l)

        r = y.ndim * [slice(None)]
        r[i] = slice(n_2, None)
        r = tuple(r)

        y = _concatenate([y[r], y[l]], axis=i)

        if len(x.chunks[i]) == 1:
            y = y.rechunk({i: x.chunks[i]})

    return y


@wraps(np.fft.fftshift)
def fftshift(x, axes=None):
    return _fftshift_helper(x, axes=axes, inverse=False)


@wraps(np.fft.ifftshift)
def ifftshift(x, axes=None):
    return _fftshift_helper(x, axes=axes, inverse=True)
