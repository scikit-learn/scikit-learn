from __future__ import division, print_function, absolute_import

import operator
from numpy import (arange, array, asarray, atleast_1d, intc, integer,
                   isscalar, issubdtype, take, unique, where)
from numpy.fft.helper import fftshift, ifftshift, fftfreq
from bisect import bisect_left

__all__ = ['fftshift', 'ifftshift', 'fftfreq', 'rfftfreq', 'next_fast_len']


def rfftfreq(n, d=1.0):
    """DFT sample frequencies (for usage with rfft, irfft).

    The returned float array contains the frequency bins in
    cycles/unit (with zero at the start) given a window length `n` and a
    sample spacing `d`::

      f = [0,1,1,2,2,...,n/2-1,n/2-1,n/2]/(d*n)   if n is even
      f = [0,1,1,2,2,...,n/2-1,n/2-1,n/2,n/2]/(d*n)   if n is odd

    Parameters
    ----------
    n : int
        Window length.
    d : scalar, optional
        Sample spacing. Default is 1.

    Returns
    -------
    out : ndarray
        The array of length `n`, containing the sample frequencies.

    Examples
    --------
    >>> from scipy import fftpack
    >>> sig = np.array([-2, 8, 6, 4, 1, 0, 3, 5], dtype=float)
    >>> sig_fft = fftpack.rfft(sig)
    >>> n = sig_fft.size
    >>> timestep = 0.1
    >>> freq = fftpack.rfftfreq(n, d=timestep)
    >>> freq
    array([ 0.  ,  1.25,  1.25,  2.5 ,  2.5 ,  3.75,  3.75,  5.  ])

    """
    n = operator.index(n)
    if n < 0:
        raise ValueError("n = %s is not valid. "
                         "n must be a nonnegative integer." % n)

    return (arange(1, n + 1, dtype=int) // 2) / float(n * d)


def next_fast_len(target):
    """
    Find the next fast size of input data to `fft`, for zero-padding, etc.

    SciPy's FFTPACK has efficient functions for radix {2, 3, 4, 5}, so this
    returns the next composite of the prime factors 2, 3, and 5 which is
    greater than or equal to `target`. (These are also known as 5-smooth
    numbers, regular numbers, or Hamming numbers.)

    Parameters
    ----------
    target : int
        Length to start searching from.  Must be a positive integer.

    Returns
    -------
    out : int
        The first 5-smooth number greater than or equal to `target`.

    Notes
    -----
    .. versionadded:: 0.18.0

    Examples
    --------
    On a particular machine, an FFT of prime length takes 133 ms:

    >>> from scipy import fftpack
    >>> min_len = 10007  # prime length is worst case for speed
    >>> a = np.random.randn(min_len)
    >>> b = fftpack.fft(a)

    Zero-padding to the next 5-smooth length reduces computation time to
    211 us, a speedup of 630 times:

    >>> fftpack.helper.next_fast_len(min_len)
    10125
    >>> b = fftpack.fft(a, 10125)

    Rounding up to the next power of 2 is not optimal, taking 367 us to
    compute, 1.7 times as long as the 5-smooth size:

    >>> b = fftpack.fft(a, 16384)

    """
    hams = (8, 9, 10, 12, 15, 16, 18, 20, 24, 25, 27, 30, 32, 36, 40, 45, 48,
            50, 54, 60, 64, 72, 75, 80, 81, 90, 96, 100, 108, 120, 125, 128,
            135, 144, 150, 160, 162, 180, 192, 200, 216, 225, 240, 243, 250,
            256, 270, 288, 300, 320, 324, 360, 375, 384, 400, 405, 432, 450,
            480, 486, 500, 512, 540, 576, 600, 625, 640, 648, 675, 720, 729,
            750, 768, 800, 810, 864, 900, 960, 972, 1000, 1024, 1080, 1125,
            1152, 1200, 1215, 1250, 1280, 1296, 1350, 1440, 1458, 1500, 1536,
            1600, 1620, 1728, 1800, 1875, 1920, 1944, 2000, 2025, 2048, 2160,
            2187, 2250, 2304, 2400, 2430, 2500, 2560, 2592, 2700, 2880, 2916,
            3000, 3072, 3125, 3200, 3240, 3375, 3456, 3600, 3645, 3750, 3840,
            3888, 4000, 4050, 4096, 4320, 4374, 4500, 4608, 4800, 4860, 5000,
            5120, 5184, 5400, 5625, 5760, 5832, 6000, 6075, 6144, 6250, 6400,
            6480, 6561, 6750, 6912, 7200, 7290, 7500, 7680, 7776, 8000, 8100,
            8192, 8640, 8748, 9000, 9216, 9375, 9600, 9720, 10000)

    target = int(target)

    if target <= 6:
        return target

    # Quickly check if it's already a power of 2
    if not (target & (target-1)):
        return target

    # Get result quickly for small sizes, since FFT itself is similarly fast.
    if target <= hams[-1]:
        return hams[bisect_left(hams, target)]

    match = float('inf')  # Anything found will be smaller
    p5 = 1
    while p5 < target:
        p35 = p5
        while p35 < target:
            # Ceiling integer division, avoiding conversion to float
            # (quotient = ceil(target / p35))
            quotient = -(-target // p35)

            # Quickly find next power of 2 >= quotient
            p2 = 2**((quotient - 1).bit_length())

            N = p2 * p35
            if N == target:
                return N
            elif N < match:
                match = N
            p35 *= 3
            if p35 == target:
                return p35
        if p35 < match:
            match = p35
        p5 *= 5
        if p5 == target:
            return p5
    if p5 < match:
        match = p5
    return match


def _init_nd_shape_and_axes(x, shape, axes):
    """Handle shape and axes arguments for n-dimensional transforms.

    Returns the shape and axes in a standard form, taking into account negative
    values and checking for various potential errors.

    Parameters
    ----------
    x : array_like
        The input array.
    shape : int or array_like of ints or None
        The shape of the result.  If both `shape` and `axes` (see below) are
        None, `shape` is ``x.shape``; if `shape` is None but `axes` is
        not None, then `shape` is ``scipy.take(x.shape, axes, axis=0)``.
        If `shape` is -1, the size of the corresponding dimension of `x` is
        used.
    axes : int or array_like of ints or None
        Axes along which the calculation is computed.
        The default is over all axes.
        Negative indices are automatically converted to their positive
        counterpart.

    Returns
    -------
    shape : array
        The shape of the result. It is a 1D integer array.
    axes : array
        The shape of the result. It is a 1D integer array.

    """
    x = asarray(x)
    noshape = shape is None
    noaxes = axes is None

    if noaxes:
        axes = arange(x.ndim, dtype=intc)
    else:
        axes = atleast_1d(axes)

    if axes.size == 0:
        axes = axes.astype(intc)

    if not axes.ndim == 1:
        raise ValueError("when given, axes values must be a scalar or vector")
    if not issubdtype(axes.dtype, integer):
        raise ValueError("when given, axes values must be integers")

    axes = where(axes < 0, axes + x.ndim, axes)

    if axes.size != 0 and (axes.max() >= x.ndim or axes.min() < 0):
        raise ValueError("axes exceeds dimensionality of input")
    if axes.size != 0 and unique(axes).shape != axes.shape:
        raise ValueError("all axes must be unique")

    if not noshape:
        shape = atleast_1d(shape)
    elif isscalar(x):
        shape = array([], dtype=intc)
    elif noaxes:
        shape = array(x.shape, dtype=intc)
    else:
        shape = take(x.shape, axes)

    if shape.size == 0:
        shape = shape.astype(intc)

    if shape.ndim != 1:
        raise ValueError("when given, shape values must be a scalar or vector")
    if not issubdtype(shape.dtype, integer):
        raise ValueError("when given, shape values must be integers")
    if axes.shape != shape.shape:
        raise ValueError("when given, axes and shape arguments"
                         " have to be of the same length")

    shape = where(shape == -1, array(x.shape)[axes], shape)

    if shape.size != 0 and (shape < 1).any():
        raise ValueError(
            "invalid number of data points ({0}) specified".format(shape))

    return shape, axes


def _init_nd_shape_and_axes_sorted(x, shape, axes):
    """Handle and sort shape and axes arguments for n-dimensional transforms.

    This is identical to `_init_nd_shape_and_axes`, except the axes are
    returned in sorted order and the shape is reordered to match.

    Parameters
    ----------
    x : array_like
        The input array.
    shape : int or array_like of ints or None
        The shape of the result.  If both `shape` and `axes` (see below) are
        None, `shape` is ``x.shape``; if `shape` is None but `axes` is
        not None, then `shape` is ``scipy.take(x.shape, axes, axis=0)``.
        If `shape` is -1, the size of the corresponding dimension of `x` is
        used.
    axes : int or array_like of ints or None
        Axes along which the calculation is computed.
        The default is over all axes.
        Negative indices are automatically converted to their positive
        counterpart.

    Returns
    -------
    shape : array
        The shape of the result. It is a 1D integer array.
    axes : array
        The shape of the result. It is a 1D integer array.

    """
    noaxes = axes is None
    shape, axes = _init_nd_shape_and_axes(x, shape, axes)

    if not noaxes:
        shape = shape[axes.argsort()]
        axes.sort()

    return shape, axes
