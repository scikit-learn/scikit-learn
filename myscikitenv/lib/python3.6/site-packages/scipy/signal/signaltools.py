# Author: Travis Oliphant
# 1999 -- 2002

from __future__ import division, print_function, absolute_import

import operator
import threading
import sys
import timeit

from . import sigtools, dlti
from ._upfirdn import upfirdn, _output_len
from scipy._lib.six import callable
from scipy._lib._version import NumpyVersion
from scipy import fftpack, linalg
from scipy.fftpack.helper import _init_nd_shape_and_axes_sorted
from numpy import (allclose, angle, arange, argsort, array, asarray,
                   atleast_1d, atleast_2d, cast, dot, exp, expand_dims,
                   iscomplexobj, mean, ndarray, newaxis, ones, pi,
                   poly, polyadd, polyder, polydiv, polymul, polysub, polyval,
                   prod, r_, ravel, real_if_close, reshape,
                   roots, sort, take, transpose, unique, where, zeros,
                   zeros_like)
import numpy as np
import math
from scipy.special import factorial
from .windows import get_window
from ._arraytools import axis_slice, axis_reverse, odd_ext, even_ext, const_ext
from .filter_design import cheby1, _validate_sos
from .fir_filter_design import firwin

if sys.version_info.major >= 3 and sys.version_info.minor >= 5:
    from math import gcd
else:
    from fractions import gcd


__all__ = ['correlate', 'fftconvolve', 'convolve', 'convolve2d', 'correlate2d',
           'order_filter', 'medfilt', 'medfilt2d', 'wiener', 'lfilter',
           'lfiltic', 'sosfilt', 'deconvolve', 'hilbert', 'hilbert2',
           'cmplx_sort', 'unique_roots', 'invres', 'invresz', 'residue',
           'residuez', 'resample', 'resample_poly', 'detrend',
           'lfilter_zi', 'sosfilt_zi', 'sosfiltfilt', 'choose_conv_method',
           'filtfilt', 'decimate', 'vectorstrength']


_modedict = {'valid': 0, 'same': 1, 'full': 2}

_boundarydict = {'fill': 0, 'pad': 0, 'wrap': 2, 'circular': 2, 'symm': 1,
                 'symmetric': 1, 'reflect': 4}


_rfft_mt_safe = (NumpyVersion(np.__version__) >= '1.9.0.dev-e24486e')

_rfft_lock = threading.Lock()


def _valfrommode(mode):
    try:
        return _modedict[mode]
    except KeyError:
        raise ValueError("Acceptable mode flags are 'valid',"
                         " 'same', or 'full'.")


def _bvalfromboundary(boundary):
    try:
        return _boundarydict[boundary] << 2
    except KeyError:
        raise ValueError("Acceptable boundary flags are 'fill', 'circular' "
                         "(or 'wrap'), and 'symmetric' (or 'symm').")


def _inputs_swap_needed(mode, shape1, shape2):
    """
    If in 'valid' mode, returns whether or not the input arrays need to be
    swapped depending on whether `shape1` is at least as large as `shape2` in
    every dimension.

    This is important for some of the correlation and convolution
    implementations in this module, where the larger array input needs to come
    before the smaller array input when operating in this mode.

    Note that if the mode provided is not 'valid', False is immediately
    returned.
    """
    if mode == 'valid':
        ok1, ok2 = True, True

        for d1, d2 in zip(shape1, shape2):
            if not d1 >= d2:
                ok1 = False
            if not d2 >= d1:
                ok2 = False

        if not (ok1 or ok2):
            raise ValueError("For 'valid' mode, one must be at least "
                             "as large as the other in every dimension")

        return not ok1

    return False


def correlate(in1, in2, mode='full', method='auto'):
    r"""
    Cross-correlate two N-dimensional arrays.

    Cross-correlate `in1` and `in2`, with the output size determined by the
    `mode` argument.

    Parameters
    ----------
    in1 : array_like
        First input.
    in2 : array_like
        Second input. Should have the same number of dimensions as `in1`.
    mode : str {'full', 'valid', 'same'}, optional
        A string indicating the size of the output:

        ``full``
           The output is the full discrete linear cross-correlation
           of the inputs. (Default)
        ``valid``
           The output consists only of those elements that do not
           rely on the zero-padding. In 'valid' mode, either `in1` or `in2`
           must be at least as large as the other in every dimension.
        ``same``
           The output is the same size as `in1`, centered
           with respect to the 'full' output.
    method : str {'auto', 'direct', 'fft'}, optional
        A string indicating which method to use to calculate the correlation.

        ``direct``
           The correlation is determined directly from sums, the definition of
           correlation.
        ``fft``
           The Fast Fourier Transform is used to perform the correlation more
           quickly (only available for numerical arrays.)
        ``auto``
           Automatically chooses direct or Fourier method based on an estimate
           of which is faster (default).  See `convolve` Notes for more detail.

           .. versionadded:: 0.19.0

    Returns
    -------
    correlate : array
        An N-dimensional array containing a subset of the discrete linear
        cross-correlation of `in1` with `in2`.

    See Also
    --------
    choose_conv_method : contains more documentation on `method`.

    Notes
    -----
    The correlation z of two d-dimensional arrays x and y is defined as::

        z[...,k,...] = sum[..., i_l, ...] x[..., i_l,...] * conj(y[..., i_l - k,...])

    This way, if x and y are 1-D arrays and ``z = correlate(x, y, 'full')``
    then

    .. math::

          z[k] = (x * y)(k - N + 1)
               = \sum_{l=0}^{||x||-1}x_l y_{l-k+N-1}^{*}

    for :math:`k = 0, 1, ..., ||x|| + ||y|| - 2`

    where :math:`||x||` is the length of ``x``, :math:`N = \max(||x||,||y||)`,
    and :math:`y_m` is 0 when m is outside the range of y.

    ``method='fft'`` only works for numerical arrays as it relies on
    `fftconvolve`. In certain cases (i.e., arrays of objects or when
    rounding integers can lose precision), ``method='direct'`` is always used.

    Examples
    --------
    Implement a matched filter using cross-correlation, to recover a signal
    that has passed through a noisy channel.

    >>> from scipy import signal
    >>> sig = np.repeat([0., 1., 1., 0., 1., 0., 0., 1.], 128)
    >>> sig_noise = sig + np.random.randn(len(sig))
    >>> corr = signal.correlate(sig_noise, np.ones(128), mode='same') / 128

    >>> import matplotlib.pyplot as plt
    >>> clock = np.arange(64, len(sig), 128)
    >>> fig, (ax_orig, ax_noise, ax_corr) = plt.subplots(3, 1, sharex=True)
    >>> ax_orig.plot(sig)
    >>> ax_orig.plot(clock, sig[clock], 'ro')
    >>> ax_orig.set_title('Original signal')
    >>> ax_noise.plot(sig_noise)
    >>> ax_noise.set_title('Signal with noise')
    >>> ax_corr.plot(corr)
    >>> ax_corr.plot(clock, corr[clock], 'ro')
    >>> ax_corr.axhline(0.5, ls=':')
    >>> ax_corr.set_title('Cross-correlated with rectangular pulse')
    >>> ax_orig.margins(0, 0.1)
    >>> fig.tight_layout()
    >>> fig.show()

    """
    in1 = asarray(in1)
    in2 = asarray(in2)

    if in1.ndim == in2.ndim == 0:
        return in1 * in2.conj()
    elif in1.ndim != in2.ndim:
        raise ValueError("in1 and in2 should have the same dimensionality")

    # Don't use _valfrommode, since correlate should not accept numeric modes
    try:
        val = _modedict[mode]
    except KeyError:
        raise ValueError("Acceptable mode flags are 'valid',"
                         " 'same', or 'full'.")

    # this either calls fftconvolve or this function with method=='direct'
    if method in ('fft', 'auto'):
        return convolve(in1, _reverse_and_conj(in2), mode, method)

    elif method == 'direct':
        # fastpath to faster numpy.correlate for 1d inputs when possible
        if _np_conv_ok(in1, in2, mode):
            return np.correlate(in1, in2, mode)

        # _correlateND is far slower when in2.size > in1.size, so swap them
        # and then undo the effect afterward if mode == 'full'.  Also, it fails
        # with 'valid' mode if in2 is larger than in1, so swap those, too.
        # Don't swap inputs for 'same' mode, since shape of in1 matters.
        swapped_inputs = ((mode == 'full') and (in2.size > in1.size) or
                          _inputs_swap_needed(mode, in1.shape, in2.shape))

        if swapped_inputs:
            in1, in2 = in2, in1

        if mode == 'valid':
            ps = [i - j + 1 for i, j in zip(in1.shape, in2.shape)]
            out = np.empty(ps, in1.dtype)

            z = sigtools._correlateND(in1, in2, out, val)

        else:
            ps = [i + j - 1 for i, j in zip(in1.shape, in2.shape)]

            # zero pad input
            in1zpadded = np.zeros(ps, in1.dtype)
            sc = tuple(slice(0, i) for i in in1.shape)
            in1zpadded[sc] = in1.copy()

            if mode == 'full':
                out = np.empty(ps, in1.dtype)
            elif mode == 'same':
                out = np.empty(in1.shape, in1.dtype)

            z = sigtools._correlateND(in1zpadded, in2, out, val)

        if swapped_inputs:
            # Reverse and conjugate to undo the effect of swapping inputs
            z = _reverse_and_conj(z)

        return z

    else:
        raise ValueError("Acceptable method flags are 'auto',"
                         " 'direct', or 'fft'.")


def _centered(arr, newshape):
    # Return the center newshape portion of the array.
    newshape = asarray(newshape)
    currshape = array(arr.shape)
    startind = (currshape - newshape) // 2
    endind = startind + newshape
    myslice = [slice(startind[k], endind[k]) for k in range(len(endind))]
    return arr[tuple(myslice)]


def fftconvolve(in1, in2, mode="full", axes=None):
    """Convolve two N-dimensional arrays using FFT.

    Convolve `in1` and `in2` using the fast Fourier transform method, with
    the output size determined by the `mode` argument.

    This is generally much faster than `convolve` for large arrays (n > ~500),
    but can be slower when only a few output values are needed, and can only
    output float arrays (int or object array inputs will be cast to float).

    As of v0.19, `convolve` automatically chooses this method or the direct
    method based on an estimation of which is faster.

    Parameters
    ----------
    in1 : array_like
        First input.
    in2 : array_like
        Second input. Should have the same number of dimensions as `in1`.
    mode : str {'full', 'valid', 'same'}, optional
        A string indicating the size of the output:

        ``full``
           The output is the full discrete linear convolution
           of the inputs. (Default)
        ``valid``
           The output consists only of those elements that do not
           rely on the zero-padding. In 'valid' mode, either `in1` or `in2`
           must be at least as large as the other in every dimension.
        ``same``
           The output is the same size as `in1`, centered
           with respect to the 'full' output.
           axis : tuple, optional
    axes : int or array_like of ints or None, optional
        Axes over which to compute the convolution.
        The default is over all axes.


    Returns
    -------
    out : array
        An N-dimensional array containing a subset of the discrete linear
        convolution of `in1` with `in2`.

    Examples
    --------
    Autocorrelation of white noise is an impulse.

    >>> from scipy import signal
    >>> sig = np.random.randn(1000)
    >>> autocorr = signal.fftconvolve(sig, sig[::-1], mode='full')

    >>> import matplotlib.pyplot as plt
    >>> fig, (ax_orig, ax_mag) = plt.subplots(2, 1)
    >>> ax_orig.plot(sig)
    >>> ax_orig.set_title('White noise')
    >>> ax_mag.plot(np.arange(-len(sig)+1,len(sig)), autocorr)
    >>> ax_mag.set_title('Autocorrelation')
    >>> fig.tight_layout()
    >>> fig.show()

    Gaussian blur implemented using FFT convolution.  Notice the dark borders
    around the image, due to the zero-padding beyond its boundaries.
    The `convolve2d` function allows for other types of image boundaries,
    but is far slower.

    >>> from scipy import misc
    >>> face = misc.face(gray=True)
    >>> kernel = np.outer(signal.gaussian(70, 8), signal.gaussian(70, 8))
    >>> blurred = signal.fftconvolve(face, kernel, mode='same')

    >>> fig, (ax_orig, ax_kernel, ax_blurred) = plt.subplots(3, 1,
    ...                                                      figsize=(6, 15))
    >>> ax_orig.imshow(face, cmap='gray')
    >>> ax_orig.set_title('Original')
    >>> ax_orig.set_axis_off()
    >>> ax_kernel.imshow(kernel, cmap='gray')
    >>> ax_kernel.set_title('Gaussian kernel')
    >>> ax_kernel.set_axis_off()
    >>> ax_blurred.imshow(blurred, cmap='gray')
    >>> ax_blurred.set_title('Blurred')
    >>> ax_blurred.set_axis_off()
    >>> fig.show()

    """
    in1 = asarray(in1)
    in2 = asarray(in2)
    noaxes = axes is None

    if in1.ndim == in2.ndim == 0:  # scalar inputs
        return in1 * in2
    elif in1.ndim != in2.ndim:
        raise ValueError("in1 and in2 should have the same dimensionality")
    elif in1.size == 0 or in2.size == 0:  # empty arrays
        return array([])

    _, axes = _init_nd_shape_and_axes_sorted(in1, shape=None, axes=axes)

    if not noaxes and not axes.size:
        raise ValueError("when provided, axes cannot be empty")

    if noaxes:
        other_axes = array([], dtype=np.intc)
    else:
        other_axes = np.setdiff1d(np.arange(in1.ndim), axes)

    s1 = array(in1.shape)
    s2 = array(in2.shape)

    if not np.all((s1[other_axes] == s2[other_axes])
                  | (s1[other_axes] == 1) | (s2[other_axes] == 1)):
        raise ValueError("incompatible shapes for in1 and in2:"
                         " {0} and {1}".format(in1.shape, in2.shape))

    complex_result = (np.issubdtype(in1.dtype, np.complexfloating)
                      or np.issubdtype(in2.dtype, np.complexfloating))
    shape = np.maximum(s1, s2)
    shape[axes] = s1[axes] + s2[axes] - 1

    # Check that input sizes are compatible with 'valid' mode
    if _inputs_swap_needed(mode, s1, s2):
        # Convolution is commutative; order doesn't have any effect on output
        in1, s1, in2, s2 = in2, s2, in1, s1

    # Speed up FFT by padding to optimal size for FFTPACK
    fshape = [fftpack.helper.next_fast_len(d) for d in shape[axes]]
    fslice = tuple([slice(sz) for sz in shape])
    # Pre-1.9 NumPy FFT routines are not threadsafe.  For older NumPys, make
    # sure we only call rfftn/irfftn from one thread at a time.
    if not complex_result and (_rfft_mt_safe or _rfft_lock.acquire(False)):
        try:
            sp1 = np.fft.rfftn(in1, fshape, axes=axes)
            sp2 = np.fft.rfftn(in2, fshape, axes=axes)
            ret = np.fft.irfftn(sp1 * sp2, fshape, axes=axes)[fslice].copy()
        finally:
            if not _rfft_mt_safe:
                _rfft_lock.release()
    else:
        # If we're here, it's either because we need a complex result, or we
        # failed to acquire _rfft_lock (meaning rfftn isn't threadsafe and
        # is already in use by another thread).  In either case, use the
        # (threadsafe but slower) SciPy complex-FFT routines instead.
        sp1 = fftpack.fftn(in1, fshape, axes=axes)
        sp2 = fftpack.fftn(in2, fshape, axes=axes)
        ret = fftpack.ifftn(sp1 * sp2, axes=axes)[fslice].copy()
        if not complex_result:
            ret = ret.real

    if mode == "full":
        return ret
    elif mode == "same":
        return _centered(ret, s1)
    elif mode == "valid":
        shape_valid = shape.copy()
        shape_valid[axes] = s1[axes] - s2[axes] + 1
        return _centered(ret, shape_valid)
    else:
        raise ValueError("acceptable mode flags are 'valid',"
                         " 'same', or 'full'")


def _numeric_arrays(arrays, kinds='buifc'):
    """
    See if a list of arrays are all numeric.

    Parameters
    ----------
    ndarrays : array or list of arrays
        arrays to check if numeric.
    numeric_kinds : string-like
        The dtypes of the arrays to be checked. If the dtype.kind of
        the ndarrays are not in this string the function returns False and
        otherwise returns True.
    """
    if type(arrays) == ndarray:
        return arrays.dtype.kind in kinds
    for array_ in arrays:
        if array_.dtype.kind not in kinds:
            return False
    return True


def _prod(iterable):
    """
    Product of a list of numbers.
    Faster than np.prod for short lists like array shapes.
    """
    product = 1
    for x in iterable:
        product *= x
    return product


def _fftconv_faster(x, h, mode):
    """
    See if using `fftconvolve` or `_correlateND` is faster. The boolean value
    returned depends on the sizes and shapes of the input values.

    The big O ratios were found to hold across different machines, which makes
    sense as it's the ratio that matters (the effective speed of the computer
    is found in both big O constants). Regardless, this had been tuned on an
    early 2015 MacBook Pro with 8GB RAM and an Intel i5 processor.
    """
    if mode == 'full':
        out_shape = [n + k - 1 for n, k in zip(x.shape, h.shape)]
        big_O_constant = 10963.92823819 if x.ndim == 1 else 8899.1104874
    elif mode == 'same':
        out_shape = x.shape
        if x.ndim == 1:
            if h.size <= x.size:
                big_O_constant = 7183.41306773
            else:
                big_O_constant = 856.78174111
        else:
            big_O_constant = 34519.21021589
    elif mode == 'valid':
        out_shape = [n - k + 1 for n, k in zip(x.shape, h.shape)]
        big_O_constant = 41954.28006344 if x.ndim == 1 else 66453.24316434
    else:
        raise ValueError("Acceptable mode flags are 'valid',"
                         " 'same', or 'full'.")

    # see whether the Fourier transform convolution method or the direct
    # convolution method is faster (discussed in scikit-image PR #1792)
    direct_time = (x.size * h.size * _prod(out_shape))
    fft_time = sum(n * math.log(n) for n in (x.shape + h.shape +
                                             tuple(out_shape)))
    return big_O_constant * fft_time < direct_time


def _reverse_and_conj(x):
    """
    Reverse array `x` in all dimensions and perform the complex conjugate
    """
    reverse = (slice(None, None, -1),) * x.ndim
    return x[reverse].conj()


def _np_conv_ok(volume, kernel, mode):
    """
    See if numpy supports convolution of `volume` and `kernel` (i.e. both are
    1D ndarrays and of the appropriate shape).  NumPy's 'same' mode uses the
    size of the larger input, while SciPy's uses the size of the first input.

    Invalid mode strings will return False and be caught by the calling func.
    """
    if volume.ndim == kernel.ndim == 1:
        if mode in ('full', 'valid'):
            return True
        elif mode == 'same':
            return volume.size >= kernel.size
    else:
        return False


def _timeit_fast(stmt="pass", setup="pass", repeat=3):
    """
    Returns the time the statement/function took, in seconds.

    Faster, less precise version of IPython's timeit. `stmt` can be a statement
    written as a string or a callable.

    Will do only 1 loop (like IPython's timeit) with no repetitions
    (unlike IPython) for very slow functions.  For fast functions, only does
    enough loops to take 5 ms, which seems to produce similar results (on
    Windows at least), and avoids doing an extraneous cycle that isn't
    measured.

    """
    timer = timeit.Timer(stmt, setup)

    # determine number of calls per rep so total time for 1 rep >= 5 ms
    x = 0
    for p in range(0, 10):
        number = 10**p
        x = timer.timeit(number)  # seconds
        if x >= 5e-3 / 10:  # 5 ms for final test, 1/10th that for this one
            break
    if x > 1:  # second
        # If it's macroscopic, don't bother with repetitions
        best = x
    else:
        number *= 10
        r = timer.repeat(repeat, number)
        best = min(r)

    sec = best / number
    return sec


def choose_conv_method(in1, in2, mode='full', measure=False):
    """
    Find the fastest convolution/correlation method.

    This primarily exists to be called during the ``method='auto'`` option in
    `convolve` and `correlate`, but can also be used when performing many
    convolutions of the same input shapes and dtypes, determining
    which method to use for all of them, either to avoid the overhead of the
    'auto' option or to use accurate real-world measurements.

    Parameters
    ----------
    in1 : array_like
        The first argument passed into the convolution function.
    in2 : array_like
        The second argument passed into the convolution function.
    mode : str {'full', 'valid', 'same'}, optional
        A string indicating the size of the output:

        ``full``
           The output is the full discrete linear convolution
           of the inputs. (Default)
        ``valid``
           The output consists only of those elements that do not
           rely on the zero-padding.
        ``same``
           The output is the same size as `in1`, centered
           with respect to the 'full' output.
    measure : bool, optional
        If True, run and time the convolution of `in1` and `in2` with both
        methods and return the fastest. If False (default), predict the fastest
        method using precomputed values.

    Returns
    -------
    method : str
        A string indicating which convolution method is fastest, either
        'direct' or 'fft'
    times : dict, optional
        A dictionary containing the times (in seconds) needed for each method.
        This value is only returned if ``measure=True``.

    See Also
    --------
    convolve
    correlate

    Notes
    -----
    For large n, ``measure=False`` is accurate and can quickly determine the
    fastest method to perform the convolution.  However, this is not as
    accurate for small n (when any dimension in the input or output is small).

    In practice, we found that this function estimates the faster method up to
    a multiplicative factor of 5 (i.e., the estimated method is *at most* 5
    times slower than the fastest method). The estimation values were tuned on
    an early 2015 MacBook Pro with 8GB RAM but we found that the prediction
    held *fairly* accurately across different machines.

    If ``measure=True``, time the convolutions. Because this function uses
    `fftconvolve`, an error will be thrown if it does not support the inputs.
    There are cases when `fftconvolve` supports the inputs but this function
    returns `direct` (e.g., to protect against floating point integer
    precision).

    .. versionadded:: 0.19

    Examples
    --------
    Estimate the fastest method for a given input:

    >>> from scipy import signal
    >>> a = np.random.randn(1000)
    >>> b = np.random.randn(1000000)
    >>> method = signal.choose_conv_method(a, b, mode='same')
    >>> method
    'fft'

    This can then be applied to other arrays of the same dtype and shape:

    >>> c = np.random.randn(1000)
    >>> d = np.random.randn(1000000)
    >>> # `method` works with correlate and convolve
    >>> corr1 = signal.correlate(a, b, mode='same', method=method)
    >>> corr2 = signal.correlate(c, d, mode='same', method=method)
    >>> conv1 = signal.convolve(a, b, mode='same', method=method)
    >>> conv2 = signal.convolve(c, d, mode='same', method=method)

    """
    volume = asarray(in1)
    kernel = asarray(in2)

    if measure:
        times = {}
        for method in ['fft', 'direct']:
            times[method] = _timeit_fast(lambda: convolve(volume, kernel,
                                         mode=mode, method=method))

        chosen_method = 'fft' if times['fft'] < times['direct'] else 'direct'
        return chosen_method, times

    # fftconvolve doesn't support complex256
    fftconv_unsup = "complex256" if sys.maxsize > 2**32 else "complex192"
    if hasattr(np, fftconv_unsup):
        if volume.dtype == fftconv_unsup or kernel.dtype == fftconv_unsup:
            return 'direct'

    # for integer input,
    # catch when more precision required than float provides (representing an
    # integer as float can lose precision in fftconvolve if larger than 2**52)
    if any([_numeric_arrays([x], kinds='ui') for x in [volume, kernel]]):
        max_value = int(np.abs(volume).max()) * int(np.abs(kernel).max())
        max_value *= int(min(volume.size, kernel.size))
        if max_value > 2**np.finfo('float').nmant - 1:
            return 'direct'

    if _numeric_arrays([volume, kernel], kinds='b'):
        return 'direct'

    if _numeric_arrays([volume, kernel]):
        if _fftconv_faster(volume, kernel, mode):
            return 'fft'

    return 'direct'


def convolve(in1, in2, mode='full', method='auto'):
    """
    Convolve two N-dimensional arrays.

    Convolve `in1` and `in2`, with the output size determined by the
    `mode` argument.

    Parameters
    ----------
    in1 : array_like
        First input.
    in2 : array_like
        Second input. Should have the same number of dimensions as `in1`.
    mode : str {'full', 'valid', 'same'}, optional
        A string indicating the size of the output:

        ``full``
           The output is the full discrete linear convolution
           of the inputs. (Default)
        ``valid``
           The output consists only of those elements that do not
           rely on the zero-padding. In 'valid' mode, either `in1` or `in2`
           must be at least as large as the other in every dimension.
        ``same``
           The output is the same size as `in1`, centered
           with respect to the 'full' output.
    method : str {'auto', 'direct', 'fft'}, optional
        A string indicating which method to use to calculate the convolution.

        ``direct``
           The convolution is determined directly from sums, the definition of
           convolution.
        ``fft``
           The Fourier Transform is used to perform the convolution by calling
           `fftconvolve`.
        ``auto``
           Automatically chooses direct or Fourier method based on an estimate
           of which is faster (default).  See Notes for more detail.

           .. versionadded:: 0.19.0

    Returns
    -------
    convolve : array
        An N-dimensional array containing a subset of the discrete linear
        convolution of `in1` with `in2`.

    See Also
    --------
    numpy.polymul : performs polynomial multiplication (same operation, but
                    also accepts poly1d objects)
    choose_conv_method : chooses the fastest appropriate convolution method
    fftconvolve

    Notes
    -----
    By default, `convolve` and `correlate` use ``method='auto'``, which calls
    `choose_conv_method` to choose the fastest method using pre-computed
    values (`choose_conv_method` can also measure real-world timing with a
    keyword argument). Because `fftconvolve` relies on floating point numbers,
    there are certain constraints that may force `method=direct` (more detail
    in `choose_conv_method` docstring).

    Examples
    --------
    Smooth a square pulse using a Hann window:

    >>> from scipy import signal
    >>> sig = np.repeat([0., 1., 0.], 100)
    >>> win = signal.hann(50)
    >>> filtered = signal.convolve(sig, win, mode='same') / sum(win)

    >>> import matplotlib.pyplot as plt
    >>> fig, (ax_orig, ax_win, ax_filt) = plt.subplots(3, 1, sharex=True)
    >>> ax_orig.plot(sig)
    >>> ax_orig.set_title('Original pulse')
    >>> ax_orig.margins(0, 0.1)
    >>> ax_win.plot(win)
    >>> ax_win.set_title('Filter impulse response')
    >>> ax_win.margins(0, 0.1)
    >>> ax_filt.plot(filtered)
    >>> ax_filt.set_title('Filtered signal')
    >>> ax_filt.margins(0, 0.1)
    >>> fig.tight_layout()
    >>> fig.show()

    """
    volume = asarray(in1)
    kernel = asarray(in2)

    if volume.ndim == kernel.ndim == 0:
        return volume * kernel
    elif volume.ndim != kernel.ndim:
        raise ValueError("volume and kernel should have the same "
                         "dimensionality")

    if _inputs_swap_needed(mode, volume.shape, kernel.shape):
        # Convolution is commutative; order doesn't have any effect on output
        volume, kernel = kernel, volume

    if method == 'auto':
        method = choose_conv_method(volume, kernel, mode=mode)

    if method == 'fft':
        out = fftconvolve(volume, kernel, mode=mode)
        result_type = np.result_type(volume, kernel)
        if result_type.kind in {'u', 'i'}:
            out = np.around(out)
        return out.astype(result_type)
    elif method == 'direct':
        # fastpath to faster numpy.convolve for 1d inputs when possible
        if _np_conv_ok(volume, kernel, mode):
            return np.convolve(volume, kernel, mode)

        return correlate(volume, _reverse_and_conj(kernel), mode, 'direct')
    else:
        raise ValueError("Acceptable method flags are 'auto',"
                         " 'direct', or 'fft'.")


def order_filter(a, domain, rank):
    """
    Perform an order filter on an N-dimensional array.

    Perform an order filter on the array in.  The domain argument acts as a
    mask centered over each pixel.  The non-zero elements of domain are
    used to select elements surrounding each input pixel which are placed
    in a list.   The list is sorted, and the output for that pixel is the
    element corresponding to rank in the sorted list.

    Parameters
    ----------
    a : ndarray
        The N-dimensional input array.
    domain : array_like
        A mask array with the same number of dimensions as `a`.
        Each dimension should have an odd number of elements.
    rank : int
        A non-negative integer which selects the element from the
        sorted list (0 corresponds to the smallest element, 1 is the
        next smallest element, etc.).

    Returns
    -------
    out : ndarray
        The results of the order filter in an array with the same
        shape as `a`.

    Examples
    --------
    >>> from scipy import signal
    >>> x = np.arange(25).reshape(5, 5)
    >>> domain = np.identity(3)
    >>> x
    array([[ 0,  1,  2,  3,  4],
           [ 5,  6,  7,  8,  9],
           [10, 11, 12, 13, 14],
           [15, 16, 17, 18, 19],
           [20, 21, 22, 23, 24]])
    >>> signal.order_filter(x, domain, 0)
    array([[  0.,   0.,   0.,   0.,   0.],
           [  0.,   0.,   1.,   2.,   0.],
           [  0.,   5.,   6.,   7.,   0.],
           [  0.,  10.,  11.,  12.,   0.],
           [  0.,   0.,   0.,   0.,   0.]])
    >>> signal.order_filter(x, domain, 2)
    array([[  6.,   7.,   8.,   9.,   4.],
           [ 11.,  12.,  13.,  14.,   9.],
           [ 16.,  17.,  18.,  19.,  14.],
           [ 21.,  22.,  23.,  24.,  19.],
           [ 20.,  21.,  22.,  23.,  24.]])

    """
    domain = asarray(domain)
    size = domain.shape
    for k in range(len(size)):
        if (size[k] % 2) != 1:
            raise ValueError("Each dimension of domain argument "
                             " should have an odd number of elements.")
    return sigtools._order_filterND(a, domain, rank)


def medfilt(volume, kernel_size=None):
    """
    Perform a median filter on an N-dimensional array.

    Apply a median filter to the input array using a local window-size
    given by `kernel_size`. The array will automatically be zero-padded.

    Parameters
    ----------
    volume : array_like
        An N-dimensional input array.
    kernel_size : array_like, optional
        A scalar or an N-length list giving the size of the median filter
        window in each dimension.  Elements of `kernel_size` should be odd.
        If `kernel_size` is a scalar, then this scalar is used as the size in
        each dimension. Default size is 3 for each dimension.

    Returns
    -------
    out : ndarray
        An array the same size as input containing the median filtered
        result.

    See also
    --------
    scipy.ndimage.median_filter

    Notes
    -------
    The more general function `scipy.ndimage.median_filter` has a more
    efficient implementation of a median filter and therefore runs much faster.
    """
    volume = atleast_1d(volume)
    if kernel_size is None:
        kernel_size = [3] * volume.ndim
    kernel_size = asarray(kernel_size)
    if kernel_size.shape == ():
        kernel_size = np.repeat(kernel_size.item(), volume.ndim)

    for k in range(volume.ndim):
        if (kernel_size[k] % 2) != 1:
            raise ValueError("Each element of kernel_size should be odd.")

    domain = ones(kernel_size)

    numels = prod(kernel_size, axis=0)
    order = numels // 2
    return sigtools._order_filterND(volume, domain, order)


def wiener(im, mysize=None, noise=None):
    """
    Perform a Wiener filter on an N-dimensional array.

    Apply a Wiener filter to the N-dimensional array `im`.

    Parameters
    ----------
    im : ndarray
        An N-dimensional array.
    mysize : int or array_like, optional
        A scalar or an N-length list giving the size of the Wiener filter
        window in each dimension.  Elements of mysize should be odd.
        If mysize is a scalar, then this scalar is used as the size
        in each dimension.
    noise : float, optional
        The noise-power to use. If None, then noise is estimated as the
        average of the local variance of the input.

    Returns
    -------
    out : ndarray
        Wiener filtered result with the same shape as `im`.

    """
    im = asarray(im)
    if mysize is None:
        mysize = [3] * im.ndim
    mysize = asarray(mysize)
    if mysize.shape == ():
        mysize = np.repeat(mysize.item(), im.ndim)

    # Estimate the local mean
    lMean = correlate(im, ones(mysize), 'same') / prod(mysize, axis=0)

    # Estimate the local variance
    lVar = (correlate(im ** 2, ones(mysize), 'same') /
            prod(mysize, axis=0) - lMean ** 2)

    # Estimate the noise power if needed.
    if noise is None:
        noise = mean(ravel(lVar), axis=0)

    res = (im - lMean)
    res *= (1 - noise / lVar)
    res += lMean
    out = where(lVar < noise, lMean, res)

    return out


def convolve2d(in1, in2, mode='full', boundary='fill', fillvalue=0):
    """
    Convolve two 2-dimensional arrays.

    Convolve `in1` and `in2` with output size determined by `mode`, and
    boundary conditions determined by `boundary` and `fillvalue`.

    Parameters
    ----------
    in1 : array_like
        First input.
    in2 : array_like
        Second input. Should have the same number of dimensions as `in1`.
    mode : str {'full', 'valid', 'same'}, optional
        A string indicating the size of the output:

        ``full``
           The output is the full discrete linear convolution
           of the inputs. (Default)
        ``valid``
           The output consists only of those elements that do not
           rely on the zero-padding. In 'valid' mode, either `in1` or `in2`
           must be at least as large as the other in every dimension.
        ``same``
           The output is the same size as `in1`, centered
           with respect to the 'full' output.
    boundary : str {'fill', 'wrap', 'symm'}, optional
        A flag indicating how to handle boundaries:

        ``fill``
           pad input arrays with fillvalue. (default)
        ``wrap``
           circular boundary conditions.
        ``symm``
           symmetrical boundary conditions.

    fillvalue : scalar, optional
        Value to fill pad input arrays with. Default is 0.

    Returns
    -------
    out : ndarray
        A 2-dimensional array containing a subset of the discrete linear
        convolution of `in1` with `in2`.

    Examples
    --------
    Compute the gradient of an image by 2D convolution with a complex Scharr
    operator.  (Horizontal operator is real, vertical is imaginary.)  Use
    symmetric boundary condition to avoid creating edges at the image
    boundaries.

    >>> from scipy import signal
    >>> from scipy import misc
    >>> ascent = misc.ascent()
    >>> scharr = np.array([[ -3-3j, 0-10j,  +3 -3j],
    ...                    [-10+0j, 0+ 0j, +10 +0j],
    ...                    [ -3+3j, 0+10j,  +3 +3j]]) # Gx + j*Gy
    >>> grad = signal.convolve2d(ascent, scharr, boundary='symm', mode='same')

    >>> import matplotlib.pyplot as plt
    >>> fig, (ax_orig, ax_mag, ax_ang) = plt.subplots(3, 1, figsize=(6, 15))
    >>> ax_orig.imshow(ascent, cmap='gray')
    >>> ax_orig.set_title('Original')
    >>> ax_orig.set_axis_off()
    >>> ax_mag.imshow(np.absolute(grad), cmap='gray')
    >>> ax_mag.set_title('Gradient magnitude')
    >>> ax_mag.set_axis_off()
    >>> ax_ang.imshow(np.angle(grad), cmap='hsv') # hsv is cyclic, like angles
    >>> ax_ang.set_title('Gradient orientation')
    >>> ax_ang.set_axis_off()
    >>> fig.show()

    """
    in1 = asarray(in1)
    in2 = asarray(in2)

    if not in1.ndim == in2.ndim == 2:
        raise ValueError('convolve2d inputs must both be 2D arrays')

    if _inputs_swap_needed(mode, in1.shape, in2.shape):
        in1, in2 = in2, in1

    val = _valfrommode(mode)
    bval = _bvalfromboundary(boundary)
    out = sigtools._convolve2d(in1, in2, 1, val, bval, fillvalue)
    return out


def correlate2d(in1, in2, mode='full', boundary='fill', fillvalue=0):
    """
    Cross-correlate two 2-dimensional arrays.

    Cross correlate `in1` and `in2` with output size determined by `mode`, and
    boundary conditions determined by `boundary` and `fillvalue`.

    Parameters
    ----------
    in1 : array_like
        First input.
    in2 : array_like
        Second input. Should have the same number of dimensions as `in1`.
    mode : str {'full', 'valid', 'same'}, optional
        A string indicating the size of the output:

        ``full``
           The output is the full discrete linear cross-correlation
           of the inputs. (Default)
        ``valid``
           The output consists only of those elements that do not
           rely on the zero-padding. In 'valid' mode, either `in1` or `in2`
           must be at least as large as the other in every dimension.
        ``same``
           The output is the same size as `in1`, centered
           with respect to the 'full' output.
    boundary : str {'fill', 'wrap', 'symm'}, optional
        A flag indicating how to handle boundaries:

        ``fill``
           pad input arrays with fillvalue. (default)
        ``wrap``
           circular boundary conditions.
        ``symm``
           symmetrical boundary conditions.

    fillvalue : scalar, optional
        Value to fill pad input arrays with. Default is 0.

    Returns
    -------
    correlate2d : ndarray
        A 2-dimensional array containing a subset of the discrete linear
        cross-correlation of `in1` with `in2`.

    Examples
    --------
    Use 2D cross-correlation to find the location of a template in a noisy
    image:

    >>> from scipy import signal
    >>> from scipy import misc
    >>> face = misc.face(gray=True) - misc.face(gray=True).mean()
    >>> template = np.copy(face[300:365, 670:750])  # right eye
    >>> template -= template.mean()
    >>> face = face + np.random.randn(*face.shape) * 50  # add noise
    >>> corr = signal.correlate2d(face, template, boundary='symm', mode='same')
    >>> y, x = np.unravel_index(np.argmax(corr), corr.shape)  # find the match

    >>> import matplotlib.pyplot as plt
    >>> fig, (ax_orig, ax_template, ax_corr) = plt.subplots(3, 1,
    ...                                                     figsize=(6, 15))
    >>> ax_orig.imshow(face, cmap='gray')
    >>> ax_orig.set_title('Original')
    >>> ax_orig.set_axis_off()
    >>> ax_template.imshow(template, cmap='gray')
    >>> ax_template.set_title('Template')
    >>> ax_template.set_axis_off()
    >>> ax_corr.imshow(corr, cmap='gray')
    >>> ax_corr.set_title('Cross-correlation')
    >>> ax_corr.set_axis_off()
    >>> ax_orig.plot(x, y, 'ro')
    >>> fig.show()

    """
    in1 = asarray(in1)
    in2 = asarray(in2)

    if not in1.ndim == in2.ndim == 2:
        raise ValueError('correlate2d inputs must both be 2D arrays')

    swapped_inputs = _inputs_swap_needed(mode, in1.shape, in2.shape)
    if swapped_inputs:
        in1, in2 = in2, in1

    val = _valfrommode(mode)
    bval = _bvalfromboundary(boundary)
    out = sigtools._convolve2d(in1, in2.conj(), 0, val, bval, fillvalue)

    if swapped_inputs:
        out = out[::-1, ::-1]

    return out


def medfilt2d(input, kernel_size=3):
    """
    Median filter a 2-dimensional array.

    Apply a median filter to the `input` array using a local window-size
    given by `kernel_size` (must be odd). The array is zero-padded
    automatically.

    Parameters
    ----------
    input : array_like
        A 2-dimensional input array.
    kernel_size : array_like, optional
        A scalar or a list of length 2, giving the size of the
        median filter window in each dimension.  Elements of
        `kernel_size` should be odd.  If `kernel_size` is a scalar,
        then this scalar is used as the size in each dimension.
        Default is a kernel of size (3, 3).

    Returns
    -------
    out : ndarray
        An array the same size as input containing the median filtered
        result.

    See also
    --------
    scipy.ndimage.median_filter

    Notes
    -------
    The more general function `scipy.ndimage.median_filter` has a more
    efficient implementation of a median filter and therefore runs much faster.
    """
    image = asarray(input)
    if kernel_size is None:
        kernel_size = [3] * 2
    kernel_size = asarray(kernel_size)
    if kernel_size.shape == ():
        kernel_size = np.repeat(kernel_size.item(), 2)

    for size in kernel_size:
        if (size % 2) != 1:
            raise ValueError("Each element of kernel_size should be odd.")

    return sigtools._medfilt2d(image, kernel_size)


def lfilter(b, a, x, axis=-1, zi=None):
    """
    Filter data along one-dimension with an IIR or FIR filter.

    Filter a data sequence, `x`, using a digital filter.  This works for many
    fundamental data types (including Object type).  The filter is a direct
    form II transposed implementation of the standard difference equation
    (see Notes).

    Parameters
    ----------
    b : array_like
        The numerator coefficient vector in a 1-D sequence.
    a : array_like
        The denominator coefficient vector in a 1-D sequence.  If ``a[0]``
        is not 1, then both `a` and `b` are normalized by ``a[0]``.
    x : array_like
        An N-dimensional input array.
    axis : int, optional
        The axis of the input data array along which to apply the
        linear filter. The filter is applied to each subarray along
        this axis.  Default is -1.
    zi : array_like, optional
        Initial conditions for the filter delays.  It is a vector
        (or array of vectors for an N-dimensional input) of length
        ``max(len(a), len(b)) - 1``.  If `zi` is None or is not given then
        initial rest is assumed.  See `lfiltic` for more information.

    Returns
    -------
    y : array
        The output of the digital filter.
    zf : array, optional
        If `zi` is None, this is not returned, otherwise, `zf` holds the
        final filter delay values.

    See Also
    --------
    lfiltic : Construct initial conditions for `lfilter`.
    lfilter_zi : Compute initial state (steady state of step response) for
                 `lfilter`.
    filtfilt : A forward-backward filter, to obtain a filter with linear phase.
    savgol_filter : A Savitzky-Golay filter.
    sosfilt: Filter data using cascaded second-order sections.
    sosfiltfilt: A forward-backward filter using second-order sections.

    Notes
    -----
    The filter function is implemented as a direct II transposed structure.
    This means that the filter implements::

       a[0]*y[n] = b[0]*x[n] + b[1]*x[n-1] + ... + b[M]*x[n-M]
                             - a[1]*y[n-1] - ... - a[N]*y[n-N]

    where `M` is the degree of the numerator, `N` is the degree of the
    denominator, and `n` is the sample number.  It is implemented using
    the following difference equations (assuming M = N)::

         a[0]*y[n] = b[0] * x[n]               + d[0][n-1]
           d[0][n] = b[1] * x[n] - a[1] * y[n] + d[1][n-1]
           d[1][n] = b[2] * x[n] - a[2] * y[n] + d[2][n-1]
         ...
         d[N-2][n] = b[N-1]*x[n] - a[N-1]*y[n] + d[N-1][n-1]
         d[N-1][n] = b[N] * x[n] - a[N] * y[n]

    where `d` are the state variables.

    The rational transfer function describing this filter in the
    z-transform domain is::

                             -1              -M
                 b[0] + b[1]z  + ... + b[M] z
         Y(z) = -------------------------------- X(z)
                             -1              -N
                 a[0] + a[1]z  + ... + a[N] z

    Examples
    --------
    Generate a noisy signal to be filtered:

    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt
    >>> t = np.linspace(-1, 1, 201)
    >>> x = (np.sin(2*np.pi*0.75*t*(1-t) + 2.1) +
    ...      0.1*np.sin(2*np.pi*1.25*t + 1) +
    ...      0.18*np.cos(2*np.pi*3.85*t))
    >>> xn = x + np.random.randn(len(t)) * 0.08

    Create an order 3 lowpass butterworth filter:

    >>> b, a = signal.butter(3, 0.05)

    Apply the filter to xn.  Use lfilter_zi to choose the initial condition of
    the filter:

    >>> zi = signal.lfilter_zi(b, a)
    >>> z, _ = signal.lfilter(b, a, xn, zi=zi*xn[0])

    Apply the filter again, to have a result filtered at an order the same as
    filtfilt:

    >>> z2, _ = signal.lfilter(b, a, z, zi=zi*z[0])

    Use filtfilt to apply the filter:

    >>> y = signal.filtfilt(b, a, xn)

    Plot the original signal and the various filtered versions:

    >>> plt.figure
    >>> plt.plot(t, xn, 'b', alpha=0.75)
    >>> plt.plot(t, z, 'r--', t, z2, 'r', t, y, 'k')
    >>> plt.legend(('noisy signal', 'lfilter, once', 'lfilter, twice',
    ...             'filtfilt'), loc='best')
    >>> plt.grid(True)
    >>> plt.show()

    """
    a = np.atleast_1d(a)
    if len(a) == 1:
        # This path only supports types fdgFDGO to mirror _linear_filter below.
        # Any of b, a, x, or zi can set the dtype, but there is no default
        # casting of other types; instead a NotImplementedError is raised.
        b = np.asarray(b)
        a = np.asarray(a)
        if b.ndim != 1 and a.ndim != 1:
            raise ValueError('object of too small depth for desired array')
        x = np.asarray(x)
        inputs = [b, a, x]
        if zi is not None:
            # _linear_filter does not broadcast zi, but does do expansion of
            # singleton dims.
            zi = np.asarray(zi)
            if zi.ndim != x.ndim:
                raise ValueError('object of too small depth for desired array')
            expected_shape = list(x.shape)
            expected_shape[axis] = b.shape[0] - 1
            expected_shape = tuple(expected_shape)
            # check the trivial case where zi is the right shape first
            if zi.shape != expected_shape:
                strides = zi.ndim * [None]
                if axis < 0:
                    axis += zi.ndim
                for k in range(zi.ndim):
                    if k == axis and zi.shape[k] == expected_shape[k]:
                        strides[k] = zi.strides[k]
                    elif k != axis and zi.shape[k] == expected_shape[k]:
                        strides[k] = zi.strides[k]
                    elif k != axis and zi.shape[k] == 1:
                        strides[k] = 0
                    else:
                        raise ValueError('Unexpected shape for zi: expected '
                                         '%s, found %s.' %
                                         (expected_shape, zi.shape))
                zi = np.lib.stride_tricks.as_strided(zi, expected_shape,
                                                     strides)
            inputs.append(zi)
        dtype = np.result_type(*inputs)

        if dtype.char not in 'fdgFDGO':
            raise NotImplementedError("input type '%s' not supported" % dtype)

        b = np.array(b, dtype=dtype)
        a = np.array(a, dtype=dtype, copy=False)
        b /= a[0]
        x = np.array(x, dtype=dtype, copy=False)

        out_full = np.apply_along_axis(lambda y: np.convolve(b, y), axis, x)
        ind = out_full.ndim * [slice(None)]
        if zi is not None:
            ind[axis] = slice(zi.shape[axis])
            out_full[tuple(ind)] += zi

        ind[axis] = slice(out_full.shape[axis] - len(b) + 1)
        out = out_full[tuple(ind)]

        if zi is None:
            return out
        else:
            ind[axis] = slice(out_full.shape[axis] - len(b) + 1, None)
            zf = out_full[tuple(ind)]
            return out, zf
    else:
        if zi is None:
            return sigtools._linear_filter(b, a, x, axis)
        else:
            return sigtools._linear_filter(b, a, x, axis, zi)


def lfiltic(b, a, y, x=None):
    """
    Construct initial conditions for lfilter given input and output vectors.

    Given a linear filter (b, a) and initial conditions on the output `y`
    and the input `x`, return the initial conditions on the state vector zi
    which is used by `lfilter` to generate the output given the input.

    Parameters
    ----------
    b : array_like
        Linear filter term.
    a : array_like
        Linear filter term.
    y : array_like
        Initial conditions.

        If ``N = len(a) - 1``, then ``y = {y[-1], y[-2], ..., y[-N]}``.

        If `y` is too short, it is padded with zeros.
    x : array_like, optional
        Initial conditions.

        If ``M = len(b) - 1``, then ``x = {x[-1], x[-2], ..., x[-M]}``.

        If `x` is not given, its initial conditions are assumed zero.

        If `x` is too short, it is padded with zeros.

    Returns
    -------
    zi : ndarray
        The state vector ``zi = {z_0[-1], z_1[-1], ..., z_K-1[-1]}``,
        where ``K = max(M, N)``.

    See Also
    --------
    lfilter, lfilter_zi

    """
    N = np.size(a) - 1
    M = np.size(b) - 1
    K = max(M, N)
    y = asarray(y)
    if y.dtype.kind in 'bui':
        # ensure calculations are floating point
        y = y.astype(np.float64)
    zi = zeros(K, y.dtype)
    if x is None:
        x = zeros(M, y.dtype)
    else:
        x = asarray(x)
        L = np.size(x)
        if L < M:
            x = r_[x, zeros(M - L)]
    L = np.size(y)
    if L < N:
        y = r_[y, zeros(N - L)]

    for m in range(M):
        zi[m] = np.sum(b[m + 1:] * x[:M - m], axis=0)

    for m in range(N):
        zi[m] -= np.sum(a[m + 1:] * y[:N - m], axis=0)

    return zi


def deconvolve(signal, divisor):
    """Deconvolves ``divisor`` out of ``signal`` using inverse filtering.

    Returns the quotient and remainder such that
    ``signal = convolve(divisor, quotient) + remainder``

    Parameters
    ----------
    signal : array_like
        Signal data, typically a recorded signal
    divisor : array_like
        Divisor data, typically an impulse response or filter that was
        applied to the original signal

    Returns
    -------
    quotient : ndarray
        Quotient, typically the recovered original signal
    remainder : ndarray
        Remainder

    Examples
    --------
    Deconvolve a signal that's been filtered:

    >>> from scipy import signal
    >>> original = [0, 1, 0, 0, 1, 1, 0, 0]
    >>> impulse_response = [2, 1]
    >>> recorded = signal.convolve(impulse_response, original)
    >>> recorded
    array([0, 2, 1, 0, 2, 3, 1, 0, 0])
    >>> recovered, remainder = signal.deconvolve(recorded, impulse_response)
    >>> recovered
    array([ 0.,  1.,  0.,  0.,  1.,  1.,  0.,  0.])

    See Also
    --------
    numpy.polydiv : performs polynomial division (same operation, but
                    also accepts poly1d objects)

    """
    num = atleast_1d(signal)
    den = atleast_1d(divisor)
    N = len(num)
    D = len(den)
    if D > N:
        quot = []
        rem = num
    else:
        input = zeros(N - D + 1, float)
        input[0] = 1
        quot = lfilter(num, den, input)
        rem = num - convolve(den, quot, mode='full')
    return quot, rem


def hilbert(x, N=None, axis=-1):
    """
    Compute the analytic signal, using the Hilbert transform.

    The transformation is done along the last axis by default.

    Parameters
    ----------
    x : array_like
        Signal data.  Must be real.
    N : int, optional
        Number of Fourier components.  Default: ``x.shape[axis]``
    axis : int, optional
        Axis along which to do the transformation.  Default: -1.

    Returns
    -------
    xa : ndarray
        Analytic signal of `x`, of each 1-D array along `axis`

    See Also
    --------
    scipy.fftpack.hilbert : Return Hilbert transform of a periodic sequence x.

    Notes
    -----
    The analytic signal ``x_a(t)`` of signal ``x(t)`` is:

    .. math:: x_a = F^{-1}(F(x) 2U) = x + i y

    where `F` is the Fourier transform, `U` the unit step function,
    and `y` the Hilbert transform of `x`. [1]_

    In other words, the negative half of the frequency spectrum is zeroed
    out, turning the real-valued signal into a complex signal.  The Hilbert
    transformed signal can be obtained from ``np.imag(hilbert(x))``, and the
    original signal from ``np.real(hilbert(x))``.

    Examples
    ---------
    In this example we use the Hilbert transform to determine the amplitude
    envelope and instantaneous frequency of an amplitude-modulated signal.

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.signal import hilbert, chirp

    >>> duration = 1.0
    >>> fs = 400.0
    >>> samples = int(fs*duration)
    >>> t = np.arange(samples) / fs

    We create a chirp of which the frequency increases from 20 Hz to 100 Hz and
    apply an amplitude modulation.

    >>> signal = chirp(t, 20.0, t[-1], 100.0)
    >>> signal *= (1.0 + 0.5 * np.sin(2.0*np.pi*3.0*t) )

    The amplitude envelope is given by magnitude of the analytic signal. The
    instantaneous frequency can be obtained by differentiating the
    instantaneous phase in respect to time. The instantaneous phase corresponds
    to the phase angle of the analytic signal.

    >>> analytic_signal = hilbert(signal)
    >>> amplitude_envelope = np.abs(analytic_signal)
    >>> instantaneous_phase = np.unwrap(np.angle(analytic_signal))
    >>> instantaneous_frequency = (np.diff(instantaneous_phase) /
    ...                            (2.0*np.pi) * fs)

    >>> fig = plt.figure()
    >>> ax0 = fig.add_subplot(211)
    >>> ax0.plot(t, signal, label='signal')
    >>> ax0.plot(t, amplitude_envelope, label='envelope')
    >>> ax0.set_xlabel("time in seconds")
    >>> ax0.legend()
    >>> ax1 = fig.add_subplot(212)
    >>> ax1.plot(t[1:], instantaneous_frequency)
    >>> ax1.set_xlabel("time in seconds")
    >>> ax1.set_ylim(0.0, 120.0)

    References
    ----------
    .. [1] Wikipedia, "Analytic signal".
           https://en.wikipedia.org/wiki/Analytic_signal
    .. [2] Leon Cohen, "Time-Frequency Analysis", 1995. Chapter 2.
    .. [3] Alan V. Oppenheim, Ronald W. Schafer. Discrete-Time Signal
           Processing, Third Edition, 2009. Chapter 12.
           ISBN 13: 978-1292-02572-8

    """
    x = asarray(x)
    if iscomplexobj(x):
        raise ValueError("x must be real.")
    if N is None:
        N = x.shape[axis]
    if N <= 0:
        raise ValueError("N must be positive.")

    Xf = fftpack.fft(x, N, axis=axis)
    h = zeros(N)
    if N % 2 == 0:
        h[0] = h[N // 2] = 1
        h[1:N // 2] = 2
    else:
        h[0] = 1
        h[1:(N + 1) // 2] = 2

    if x.ndim > 1:
        ind = [newaxis] * x.ndim
        ind[axis] = slice(None)
        h = h[tuple(ind)]
    x = fftpack.ifft(Xf * h, axis=axis)
    return x


def hilbert2(x, N=None):
    """
    Compute the '2-D' analytic signal of `x`

    Parameters
    ----------
    x : array_like
        2-D signal data.
    N : int or tuple of two ints, optional
        Number of Fourier components. Default is ``x.shape``

    Returns
    -------
    xa : ndarray
        Analytic signal of `x` taken along axes (0,1).

    References
    ----------
    .. [1] Wikipedia, "Analytic signal",
        https://en.wikipedia.org/wiki/Analytic_signal

    """
    x = atleast_2d(x)
    if x.ndim > 2:
        raise ValueError("x must be 2-D.")
    if iscomplexobj(x):
        raise ValueError("x must be real.")
    if N is None:
        N = x.shape
    elif isinstance(N, int):
        if N <= 0:
            raise ValueError("N must be positive.")
        N = (N, N)
    elif len(N) != 2 or np.any(np.asarray(N) <= 0):
        raise ValueError("When given as a tuple, N must hold exactly "
                         "two positive integers")

    Xf = fftpack.fft2(x, N, axes=(0, 1))
    h1 = zeros(N[0], 'd')
    h2 = zeros(N[1], 'd')
    for p in range(2):
        h = eval("h%d" % (p + 1))
        N1 = N[p]
        if N1 % 2 == 0:
            h[0] = h[N1 // 2] = 1
            h[1:N1 // 2] = 2
        else:
            h[0] = 1
            h[1:(N1 + 1) // 2] = 2
        exec("h%d = h" % (p + 1), globals(), locals())

    h = h1[:, newaxis] * h2[newaxis, :]
    k = x.ndim
    while k > 2:
        h = h[:, newaxis]
        k -= 1
    x = fftpack.ifft2(Xf * h, axes=(0, 1))
    return x


def cmplx_sort(p):
    """Sort roots based on magnitude.

    Parameters
    ----------
    p : array_like
        The roots to sort, as a 1-D array.

    Returns
    -------
    p_sorted : ndarray
        Sorted roots.
    indx : ndarray
        Array of indices needed to sort the input `p`.

    Examples
    --------
    >>> from scipy import signal
    >>> vals = [1, 4, 1+1.j, 3]
    >>> p_sorted, indx = signal.cmplx_sort(vals)
    >>> p_sorted
    array([1.+0.j, 1.+1.j, 3.+0.j, 4.+0.j])
    >>> indx
    array([0, 2, 3, 1])

    """
    p = asarray(p)
    if iscomplexobj(p):
        indx = argsort(abs(p))
    else:
        indx = argsort(p)
    return take(p, indx, 0), indx


def unique_roots(p, tol=1e-3, rtype='min'):
    """
    Determine unique roots and their multiplicities from a list of roots.

    Parameters
    ----------
    p : array_like
        The list of roots.
    tol : float, optional
        The tolerance for two roots to be considered equal. Default is 1e-3.
    rtype : {'max', 'min, 'avg'}, optional
        How to determine the returned root if multiple roots are within
        `tol` of each other.

          - 'max': pick the maximum of those roots.
          - 'min': pick the minimum of those roots.
          - 'avg': take the average of those roots.

    Returns
    -------
    pout : ndarray
        The list of unique roots, sorted from low to high.
    mult : ndarray
        The multiplicity of each root.

    Notes
    -----
    This utility function is not specific to roots but can be used for any
    sequence of values for which uniqueness and multiplicity has to be
    determined. For a more general routine, see `numpy.unique`.

    Examples
    --------
    >>> from scipy import signal
    >>> vals = [0, 1.3, 1.31, 2.8, 1.25, 2.2, 10.3]
    >>> uniq, mult = signal.unique_roots(vals, tol=2e-2, rtype='avg')

    Check which roots have multiplicity larger than 1:

    >>> uniq[mult > 1]
    array([ 1.305])

    """
    if rtype in ['max', 'maximum']:
        comproot = np.max
    elif rtype in ['min', 'minimum']:
        comproot = np.min
    elif rtype in ['avg', 'mean']:
        comproot = np.mean
    else:
        raise ValueError("`rtype` must be one of "
                         "{'max', 'maximum', 'min', 'minimum', 'avg', 'mean'}")
    p = asarray(p) * 1.0
    tol = abs(tol)
    p, indx = cmplx_sort(p)
    pout = []
    mult = []
    indx = -1
    curp = p[0] + 5 * tol
    sameroots = []
    for k in range(len(p)):
        tr = p[k]
        if abs(tr - curp) < tol:
            sameroots.append(tr)
            curp = comproot(sameroots)
            pout[indx] = curp
            mult[indx] += 1
        else:
            pout.append(tr)
            curp = tr
            sameroots = [tr]
            indx += 1
            mult.append(1)
    return array(pout), array(mult)


def invres(r, p, k, tol=1e-3, rtype='avg'):
    """
    Compute b(s) and a(s) from partial fraction expansion.

    If `M` is the degree of numerator `b` and `N` the degree of denominator
    `a`::

              b(s)     b[0] s**(M) + b[1] s**(M-1) + ... + b[M]
      H(s) = ------ = ------------------------------------------
              a(s)     a[0] s**(N) + a[1] s**(N-1) + ... + a[N]

    then the partial-fraction expansion H(s) is defined as::

               r[0]       r[1]             r[-1]
           = -------- + -------- + ... + --------- + k(s)
             (s-p[0])   (s-p[1])         (s-p[-1])

    If there are any repeated roots (closer together than `tol`), then H(s)
    has terms like::

          r[i]      r[i+1]              r[i+n-1]
        -------- + ----------- + ... + -----------
        (s-p[i])  (s-p[i])**2          (s-p[i])**n

    This function is used for polynomials in positive powers of s or z,
    such as analog filters or digital filters in controls engineering.  For
    negative powers of z (typical for digital filters in DSP), use `invresz`.

    Parameters
    ----------
    r : array_like
        Residues.
    p : array_like
        Poles.
    k : array_like
        Coefficients of the direct polynomial term.
    tol : float, optional
        The tolerance for two roots to be considered equal. Default is 1e-3.
    rtype : {'max', 'min, 'avg'}, optional
        How to determine the returned root if multiple roots are within
        `tol` of each other.

          - 'max': pick the maximum of those roots.
          - 'min': pick the minimum of those roots.
          - 'avg': take the average of those roots.

    Returns
    -------
    b : ndarray
        Numerator polynomial coefficients.
    a : ndarray
        Denominator polynomial coefficients.

    See Also
    --------
    residue, invresz, unique_roots

    """
    extra = k
    p, indx = cmplx_sort(p)
    r = take(r, indx, 0)
    pout, mult = unique_roots(p, tol=tol, rtype=rtype)
    p = []
    for k in range(len(pout)):
        p.extend([pout[k]] * mult[k])
    a = atleast_1d(poly(p))
    if len(extra) > 0:
        b = polymul(extra, a)
    else:
        b = [0]
    indx = 0
    for k in range(len(pout)):
        temp = []
        for l in range(len(pout)):
            if l != k:
                temp.extend([pout[l]] * mult[l])
        for m in range(mult[k]):
            t2 = temp[:]
            t2.extend([pout[k]] * (mult[k] - m - 1))
            b = polyadd(b, r[indx] * atleast_1d(poly(t2)))
            indx += 1
    b = real_if_close(b)
    while allclose(b[0], 0, rtol=1e-14) and (b.shape[-1] > 1):
        b = b[1:]
    return b, a


def residue(b, a, tol=1e-3, rtype='avg'):
    """
    Compute partial-fraction expansion of b(s) / a(s).

    If `M` is the degree of numerator `b` and `N` the degree of denominator
    `a`::

              b(s)     b[0] s**(M) + b[1] s**(M-1) + ... + b[M]
      H(s) = ------ = ------------------------------------------
              a(s)     a[0] s**(N) + a[1] s**(N-1) + ... + a[N]

    then the partial-fraction expansion H(s) is defined as::

               r[0]       r[1]             r[-1]
           = -------- + -------- + ... + --------- + k(s)
             (s-p[0])   (s-p[1])         (s-p[-1])

    If there are any repeated roots (closer together than `tol`), then H(s)
    has terms like::

          r[i]      r[i+1]              r[i+n-1]
        -------- + ----------- + ... + -----------
        (s-p[i])  (s-p[i])**2          (s-p[i])**n

    This function is used for polynomials in positive powers of s or z,
    such as analog filters or digital filters in controls engineering.  For
    negative powers of z (typical for digital filters in DSP), use `residuez`.

    Parameters
    ----------
    b : array_like
        Numerator polynomial coefficients.
    a : array_like
        Denominator polynomial coefficients.

    Returns
    -------
    r : ndarray
        Residues.
    p : ndarray
        Poles.
    k : ndarray
        Coefficients of the direct polynomial term.

    See Also
    --------
    invres, residuez, numpy.poly, unique_roots

    """

    b, a = map(asarray, (b, a))
    rscale = a[0]
    k, b = polydiv(b, a)
    p = roots(a)
    r = p * 0.0
    pout, mult = unique_roots(p, tol=tol, rtype=rtype)
    p = []
    for n in range(len(pout)):
        p.extend([pout[n]] * mult[n])
    p = asarray(p)
    # Compute the residue from the general formula
    indx = 0
    for n in range(len(pout)):
        bn = b.copy()
        pn = []
        for l in range(len(pout)):
            if l != n:
                pn.extend([pout[l]] * mult[l])
        an = atleast_1d(poly(pn))
        # bn(s) / an(s) is (s-po[n])**Nn * b(s) / a(s) where Nn is
        # multiplicity of pole at po[n]
        sig = mult[n]
        for m in range(sig, 0, -1):
            if sig > m:
                # compute next derivative of bn(s) / an(s)
                term1 = polymul(polyder(bn, 1), an)
                term2 = polymul(bn, polyder(an, 1))
                bn = polysub(term1, term2)
                an = polymul(an, an)
            r[indx + m - 1] = (polyval(bn, pout[n]) / polyval(an, pout[n]) /
                               factorial(sig - m))
        indx += sig
    return r / rscale, p, k


def residuez(b, a, tol=1e-3, rtype='avg'):
    """
    Compute partial-fraction expansion of b(z) / a(z).

    If `M` is the degree of numerator `b` and `N` the degree of denominator
    `a`::

                b(z)     b[0] + b[1] z**(-1) + ... + b[M] z**(-M)
        H(z) = ------ = ------------------------------------------
                a(z)     a[0] + a[1] z**(-1) + ... + a[N] z**(-N)

    then the partial-fraction expansion H(z) is defined as::

                 r[0]                   r[-1]
         = --------------- + ... + ---------------- + k[0] + k[1]z**(-1) ...
           (1-p[0]z**(-1))         (1-p[-1]z**(-1))

    If there are any repeated roots (closer than `tol`), then the partial
    fraction expansion has terms like::

             r[i]              r[i+1]                    r[i+n-1]
        -------------- + ------------------ + ... + ------------------
        (1-p[i]z**(-1))  (1-p[i]z**(-1))**2         (1-p[i]z**(-1))**n

    This function is used for polynomials in negative powers of z,
    such as digital filters in DSP.  For positive powers, use `residue`.

    Parameters
    ----------
    b : array_like
        Numerator polynomial coefficients.
    a : array_like
        Denominator polynomial coefficients.

    Returns
    -------
    r : ndarray
        Residues.
    p : ndarray
        Poles.
    k : ndarray
        Coefficients of the direct polynomial term.

    See Also
    --------
    invresz, residue, unique_roots

    """
    b, a = map(asarray, (b, a))
    gain = a[0]
    brev, arev = b[::-1], a[::-1]
    krev, brev = polydiv(brev, arev)
    if krev == []:
        k = []
    else:
        k = krev[::-1]
    b = brev[::-1]
    p = roots(a)
    r = p * 0.0
    pout, mult = unique_roots(p, tol=tol, rtype=rtype)
    p = []
    for n in range(len(pout)):
        p.extend([pout[n]] * mult[n])
    p = asarray(p)
    # Compute the residue from the general formula (for discrete-time)
    #  the polynomial is in z**(-1) and the multiplication is by terms
    #  like this (1-p[i] z**(-1))**mult[i].  After differentiation,
    #  we must divide by (-p[i])**(m-k) as well as (m-k)!
    indx = 0
    for n in range(len(pout)):
        bn = brev.copy()
        pn = []
        for l in range(len(pout)):
            if l != n:
                pn.extend([pout[l]] * mult[l])
        an = atleast_1d(poly(pn))[::-1]
        # bn(z) / an(z) is (1-po[n] z**(-1))**Nn * b(z) / a(z) where Nn is
        # multiplicity of pole at po[n] and b(z) and a(z) are polynomials.
        sig = mult[n]
        for m in range(sig, 0, -1):
            if sig > m:
                # compute next derivative of bn(s) / an(s)
                term1 = polymul(polyder(bn, 1), an)
                term2 = polymul(bn, polyder(an, 1))
                bn = polysub(term1, term2)
                an = polymul(an, an)
            r[indx + m - 1] = (polyval(bn, 1.0 / pout[n]) /
                               polyval(an, 1.0 / pout[n]) /
                               factorial(sig - m) / (-pout[n]) ** (sig - m))
        indx += sig
    return r / gain, p, k


def invresz(r, p, k, tol=1e-3, rtype='avg'):
    """
    Compute b(z) and a(z) from partial fraction expansion.

    If `M` is the degree of numerator `b` and `N` the degree of denominator
    `a`::

                b(z)     b[0] + b[1] z**(-1) + ... + b[M] z**(-M)
        H(z) = ------ = ------------------------------------------
                a(z)     a[0] + a[1] z**(-1) + ... + a[N] z**(-N)

    then the partial-fraction expansion H(z) is defined as::

                 r[0]                   r[-1]
         = --------------- + ... + ---------------- + k[0] + k[1]z**(-1) ...
           (1-p[0]z**(-1))         (1-p[-1]z**(-1))

    If there are any repeated roots (closer than `tol`), then the partial
    fraction expansion has terms like::

             r[i]              r[i+1]                    r[i+n-1]
        -------------- + ------------------ + ... + ------------------
        (1-p[i]z**(-1))  (1-p[i]z**(-1))**2         (1-p[i]z**(-1))**n

    This function is used for polynomials in negative powers of z,
    such as digital filters in DSP.  For positive powers, use `invres`.

    Parameters
    ----------
    r : array_like
        Residues.
    p : array_like
        Poles.
    k : array_like
        Coefficients of the direct polynomial term.
    tol : float, optional
        The tolerance for two roots to be considered equal. Default is 1e-3.
    rtype : {'max', 'min, 'avg'}, optional
        How to determine the returned root if multiple roots are within
        `tol` of each other.

          - 'max': pick the maximum of those roots.
          - 'min': pick the minimum of those roots.
          - 'avg': take the average of those roots.

    Returns
    -------
    b : ndarray
        Numerator polynomial coefficients.
    a : ndarray
        Denominator polynomial coefficients.

    See Also
    --------
    residuez, unique_roots, invres

    """
    extra = asarray(k)
    p, indx = cmplx_sort(p)
    r = take(r, indx, 0)
    pout, mult = unique_roots(p, tol=tol, rtype=rtype)
    p = []
    for k in range(len(pout)):
        p.extend([pout[k]] * mult[k])
    a = atleast_1d(poly(p))
    if len(extra) > 0:
        b = polymul(extra, a)
    else:
        b = [0]
    indx = 0
    brev = asarray(b)[::-1]
    for k in range(len(pout)):
        temp = []
        # Construct polynomial which does not include any of this root
        for l in range(len(pout)):
            if l != k:
                temp.extend([pout[l]] * mult[l])
        for m in range(mult[k]):
            t2 = temp[:]
            t2.extend([pout[k]] * (mult[k] - m - 1))
            brev = polyadd(brev, (r[indx] * atleast_1d(poly(t2)))[::-1])
            indx += 1
    b = real_if_close(brev[::-1])
    return b, a


def resample(x, num, t=None, axis=0, window=None):
    """
    Resample `x` to `num` samples using Fourier method along the given axis.

    The resampled signal starts at the same value as `x` but is sampled
    with a spacing of ``len(x) / num * (spacing of x)``.  Because a
    Fourier method is used, the signal is assumed to be periodic.

    Parameters
    ----------
    x : array_like
        The data to be resampled.
    num : int
        The number of samples in the resampled signal.
    t : array_like, optional
        If `t` is given, it is assumed to be the sample positions
        associated with the signal data in `x`.
    axis : int, optional
        The axis of `x` that is resampled.  Default is 0.
    window : array_like, callable, string, float, or tuple, optional
        Specifies the window applied to the signal in the Fourier
        domain.  See below for details.

    Returns
    -------
    resampled_x or (resampled_x, resampled_t)
        Either the resampled array, or, if `t` was given, a tuple
        containing the resampled array and the corresponding resampled
        positions.

    See Also
    --------
    decimate : Downsample the signal after applying an FIR or IIR filter.
    resample_poly : Resample using polyphase filtering and an FIR filter.

    Notes
    -----
    The argument `window` controls a Fourier-domain window that tapers
    the Fourier spectrum before zero-padding to alleviate ringing in
    the resampled values for sampled signals you didn't intend to be
    interpreted as band-limited.

    If `window` is a function, then it is called with a vector of inputs
    indicating the frequency bins (i.e. fftfreq(x.shape[axis]) ).

    If `window` is an array of the same length as `x.shape[axis]` it is
    assumed to be the window to be applied directly in the Fourier
    domain (with dc and low-frequency first).

    For any other type of `window`, the function `scipy.signal.get_window`
    is called to generate the window.

    The first sample of the returned vector is the same as the first
    sample of the input vector.  The spacing between samples is changed
    from ``dx`` to ``dx * len(x) / num``.

    If `t` is not None, then it represents the old sample positions,
    and the new sample positions will be returned as well as the new
    samples.

    As noted, `resample` uses FFT transformations, which can be very
    slow if the number of input or output samples is large and prime;
    see `scipy.fftpack.fft`.

    Examples
    --------
    Note that the end of the resampled data rises to meet the first
    sample of the next cycle:

    >>> from scipy import signal

    >>> x = np.linspace(0, 10, 20, endpoint=False)
    >>> y = np.cos(-x**2/6.0)
    >>> f = signal.resample(y, 100)
    >>> xnew = np.linspace(0, 10, 100, endpoint=False)

    >>> import matplotlib.pyplot as plt
    >>> plt.plot(x, y, 'go-', xnew, f, '.-', 10, y[0], 'ro')
    >>> plt.legend(['data', 'resampled'], loc='best')
    >>> plt.show()
    """
    x = asarray(x)
    X = fftpack.fft(x, axis=axis)
    Nx = x.shape[axis]
    if window is not None:
        if callable(window):
            W = window(fftpack.fftfreq(Nx))
        elif isinstance(window, ndarray):
            if window.shape != (Nx,):
                raise ValueError('window must have the same length as data')
            W = window
        else:
            W = fftpack.ifftshift(get_window(window, Nx))
        newshape = [1] * x.ndim
        newshape[axis] = len(W)
        W.shape = newshape
        X = X * W
        W.shape = (Nx,)
    sl = [slice(None)] * x.ndim
    newshape = list(x.shape)
    newshape[axis] = num
    N = int(np.minimum(num, Nx))
    Y = zeros(newshape, 'D')
    sl[axis] = slice(0, (N + 1) // 2)
    Y[tuple(sl)] = X[tuple(sl)]
    sl[axis] = slice(-(N - 1) // 2, None)
    Y[tuple(sl)] = X[tuple(sl)]

    if N % 2 == 0:  # special treatment if low number of points is even. So far we have set Y[-N/2]=X[-N/2]
        if N < Nx:  # if downsampling
            sl[axis] = slice(N//2,N//2+1,None)  # select the component at frequency N/2
            Y[tuple(sl)] += X[tuple(sl)]  # add the component of X at N/2
        elif N < num:  # if upsampling
            sl[axis] = slice(num-N//2,num-N//2+1,None)  # select the component at frequency -N/2
            Y[tuple(sl)] /= 2  # halve the component at -N/2
            temp = Y[tuple(sl)]
            sl[axis] = slice(N//2,N//2+1,None)  # select the component at +N/2
            Y[tuple(sl)] = temp  # set that equal to the component at -N/2

    y = fftpack.ifft(Y, axis=axis) * (float(num) / float(Nx))

    if x.dtype.char not in ['F', 'D']:
        y = y.real

    if t is None:
        return y
    else:
        new_t = arange(0, num) * (t[1] - t[0]) * Nx / float(num) + t[0]
        return y, new_t


def resample_poly(x, up, down, axis=0, window=('kaiser', 5.0)):
    """
    Resample `x` along the given axis using polyphase filtering.

    The signal `x` is upsampled by the factor `up`, a zero-phase low-pass
    FIR filter is applied, and then it is downsampled by the factor `down`.
    The resulting sample rate is ``up / down`` times the original sample
    rate. Values beyond the boundary of the signal are assumed to be zero
    during the filtering step.

    Parameters
    ----------
    x : array_like
        The data to be resampled.
    up : int
        The upsampling factor.
    down : int
        The downsampling factor.
    axis : int, optional
        The axis of `x` that is resampled. Default is 0.
    window : string, tuple, or array_like, optional
        Desired window to use to design the low-pass filter, or the FIR filter
        coefficients to employ. See below for details.

    Returns
    -------
    resampled_x : array
        The resampled array.

    See Also
    --------
    decimate : Downsample the signal after applying an FIR or IIR filter.
    resample : Resample up or down using the FFT method.

    Notes
    -----
    This polyphase method will likely be faster than the Fourier method
    in `scipy.signal.resample` when the number of samples is large and
    prime, or when the number of samples is large and `up` and `down`
    share a large greatest common denominator. The length of the FIR
    filter used will depend on ``max(up, down) // gcd(up, down)``, and
    the number of operations during polyphase filtering will depend on
    the filter length and `down` (see `scipy.signal.upfirdn` for details).

    The argument `window` specifies the FIR low-pass filter design.

    If `window` is an array_like it is assumed to be the FIR filter
    coefficients. Note that the FIR filter is applied after the upsampling
    step, so it should be designed to operate on a signal at a sampling
    frequency higher than the original by a factor of `up//gcd(up, down)`.
    This function's output will be centered with respect to this array, so it
    is best to pass a symmetric filter with an odd number of samples if, as
    is usually the case, a zero-phase filter is desired.

    For any other type of `window`, the functions `scipy.signal.get_window`
    and `scipy.signal.firwin` are called to generate the appropriate filter
    coefficients.

    The first sample of the returned vector is the same as the first
    sample of the input vector. The spacing between samples is changed
    from ``dx`` to ``dx * down / float(up)``.

    Examples
    --------
    Note that the end of the resampled data rises to meet the first
    sample of the next cycle for the FFT method, and gets closer to zero
    for the polyphase method:

    >>> from scipy import signal

    >>> x = np.linspace(0, 10, 20, endpoint=False)
    >>> y = np.cos(-x**2/6.0)
    >>> f_fft = signal.resample(y, 100)
    >>> f_poly = signal.resample_poly(y, 100, 20)
    >>> xnew = np.linspace(0, 10, 100, endpoint=False)

    >>> import matplotlib.pyplot as plt
    >>> plt.plot(xnew, f_fft, 'b.-', xnew, f_poly, 'r.-')
    >>> plt.plot(x, y, 'ko-')
    >>> plt.plot(10, y[0], 'bo', 10, 0., 'ro')  # boundaries
    >>> plt.legend(['resample', 'resamp_poly', 'data'], loc='best')
    >>> plt.show()
    """
    x = asarray(x)
    if up != int(up):
        raise ValueError("up must be an integer")
    if down != int(down):
        raise ValueError("down must be an integer")
    up = int(up)
    down = int(down)
    if up < 1 or down < 1:
        raise ValueError('up and down must be >= 1')

    # Determine our up and down factors
    # Use a rational approximation to save computation time on really long
    # signals
    g_ = gcd(up, down)
    up //= g_
    down //= g_
    if up == down == 1:
        return x.copy()
    n_out = x.shape[axis] * up
    n_out = n_out // down + bool(n_out % down)

    if isinstance(window, (list, np.ndarray)):
        window = array(window)  # use array to force a copy (we modify it)
        if window.ndim > 1:
            raise ValueError('window must be 1-D')
        half_len = (window.size - 1) // 2
        h = window
    else:
        # Design a linear-phase low-pass FIR filter
        max_rate = max(up, down)
        f_c = 1. / max_rate  # cutoff of FIR filter (rel. to Nyquist)
        half_len = 10 * max_rate  # reasonable cutoff for our sinc-like function
        h = firwin(2 * half_len + 1, f_c, window=window)
    h *= up

    # Zero-pad our filter to put the output samples at the center
    n_pre_pad = (down - half_len % down)
    n_post_pad = 0
    n_pre_remove = (half_len + n_pre_pad) // down
    # We should rarely need to do this given our filter lengths...
    while _output_len(len(h) + n_pre_pad + n_post_pad, x.shape[axis],
                      up, down) < n_out + n_pre_remove:
        n_post_pad += 1
    h = np.concatenate((np.zeros(n_pre_pad, dtype=h.dtype), h,
                        np.zeros(n_post_pad, dtype=h.dtype)))
    n_pre_remove_end = n_pre_remove + n_out

    # filter then remove excess
    y = upfirdn(h, x, up, down, axis=axis)
    keep = [slice(None), ]*x.ndim
    keep[axis] = slice(n_pre_remove, n_pre_remove_end)
    return y[tuple(keep)]


def vectorstrength(events, period):
    '''
    Determine the vector strength of the events corresponding to the given
    period.

    The vector strength is a measure of phase synchrony, how well the
    timing of the events is synchronized to a single period of a periodic
    signal.

    If multiple periods are used, calculate the vector strength of each.
    This is called the "resonating vector strength".

    Parameters
    ----------
    events : 1D array_like
        An array of time points containing the timing of the events.
    period : float or array_like
        The period of the signal that the events should synchronize to.
        The period is in the same units as `events`.  It can also be an array
        of periods, in which case the outputs are arrays of the same length.

    Returns
    -------
    strength : float or 1D array
        The strength of the synchronization.  1.0 is perfect synchronization
        and 0.0 is no synchronization.  If `period` is an array, this is also
        an array with each element containing the vector strength at the
        corresponding period.
    phase : float or array
        The phase that the events are most strongly synchronized to in radians.
        If `period` is an array, this is also an array with each element
        containing the phase for the corresponding period.

    References
    ----------
    van Hemmen, JL, Longtin, A, and Vollmayr, AN. Testing resonating vector
        strength: Auditory system, electric fish, and noise.
        Chaos 21, 047508 (2011);
        :doi:`10.1063/1.3670512`.
    van Hemmen, JL.  Vector strength after Goldberg, Brown, and von Mises:
        biological and mathematical perspectives.  Biol Cybern.
        2013 Aug;107(4):385-96. :doi:`10.1007/s00422-013-0561-7`.
    van Hemmen, JL and Vollmayr, AN.  Resonating vector strength: what happens
        when we vary the "probing" frequency while keeping the spike times
        fixed.  Biol Cybern. 2013 Aug;107(4):491-94.
        :doi:`10.1007/s00422-013-0560-8`.
    '''
    events = asarray(events)
    period = asarray(period)
    if events.ndim > 1:
        raise ValueError('events cannot have dimensions more than 1')
    if period.ndim > 1:
        raise ValueError('period cannot have dimensions more than 1')

    # we need to know later if period was originally a scalar
    scalarperiod = not period.ndim

    events = atleast_2d(events)
    period = atleast_2d(period)
    if (period <= 0).any():
        raise ValueError('periods must be positive')

    # this converts the times to vectors
    vectors = exp(dot(2j*pi/period.T, events))

    # the vector strength is just the magnitude of the mean of the vectors
    # the vector phase is the angle of the mean of the vectors
    vectormean = mean(vectors, axis=1)
    strength = abs(vectormean)
    phase = angle(vectormean)

    # if the original period was a scalar, return scalars
    if scalarperiod:
        strength = strength[0]
        phase = phase[0]
    return strength, phase


def detrend(data, axis=-1, type='linear', bp=0, overwrite_data=False):
    """
    Remove linear trend along axis from data.

    Parameters
    ----------
    data : array_like
        The input data.
    axis : int, optional
        The axis along which to detrend the data. By default this is the
        last axis (-1).
    type : {'linear', 'constant'}, optional
        The type of detrending. If ``type == 'linear'`` (default),
        the result of a linear least-squares fit to `data` is subtracted
        from `data`.
        If ``type == 'constant'``, only the mean of `data` is subtracted.
    bp : array_like of ints, optional
        A sequence of break points. If given, an individual linear fit is
        performed for each part of `data` between two break points.
        Break points are specified as indices into `data`.
    overwrite_data : bool, optional
        If True, perform in place detrending and avoid a copy. Default is False

    Returns
    -------
    ret : ndarray
        The detrended input data.

    Examples
    --------
    >>> from scipy import signal
    >>> randgen = np.random.RandomState(9)
    >>> npoints = 1000
    >>> noise = randgen.randn(npoints)
    >>> x = 3 + 2*np.linspace(0, 1, npoints) + noise
    >>> (signal.detrend(x) - noise).max() < 0.01
    True

    """
    if type not in ['linear', 'l', 'constant', 'c']:
        raise ValueError("Trend type must be 'linear' or 'constant'.")
    data = asarray(data)
    dtype = data.dtype.char
    if dtype not in 'dfDF':
        dtype = 'd'
    if type in ['constant', 'c']:
        ret = data - expand_dims(mean(data, axis), axis)
        return ret
    else:
        dshape = data.shape
        N = dshape[axis]
        bp = sort(unique(r_[0, bp, N]))
        if np.any(bp > N):
            raise ValueError("Breakpoints must be less than length "
                             "of data along given axis.")
        Nreg = len(bp) - 1
        # Restructure data so that axis is along first dimension and
        #  all other dimensions are collapsed into second dimension
        rnk = len(dshape)
        if axis < 0:
            axis = axis + rnk
        newdims = r_[axis, 0:axis, axis + 1:rnk]
        newdata = reshape(transpose(data, tuple(newdims)),
                          (N, _prod(dshape) // N))
        if not overwrite_data:
            newdata = newdata.copy()  # make sure we have a copy
        if newdata.dtype.char not in 'dfDF':
            newdata = newdata.astype(dtype)
        # Find leastsq fit and remove it for each piece
        for m in range(Nreg):
            Npts = bp[m + 1] - bp[m]
            A = ones((Npts, 2), dtype)
            A[:, 0] = cast[dtype](arange(1, Npts + 1) * 1.0 / Npts)
            sl = slice(bp[m], bp[m + 1])
            coef, resids, rank, s = linalg.lstsq(A, newdata[sl])
            newdata[sl] = newdata[sl] - dot(A, coef)
        # Put data back in original shape.
        tdshape = take(dshape, newdims, 0)
        ret = reshape(newdata, tuple(tdshape))
        vals = list(range(1, rnk))
        olddims = vals[:axis] + [0] + vals[axis:]
        ret = transpose(ret, tuple(olddims))
        return ret


def lfilter_zi(b, a):
    """
    Construct initial conditions for lfilter for step response steady-state.

    Compute an initial state `zi` for the `lfilter` function that corresponds
    to the steady state of the step response.

    A typical use of this function is to set the initial state so that the
    output of the filter starts at the same value as the first element of
    the signal to be filtered.

    Parameters
    ----------
    b, a : array_like (1-D)
        The IIR filter coefficients. See `lfilter` for more
        information.

    Returns
    -------
    zi : 1-D ndarray
        The initial state for the filter.

    See Also
    --------
    lfilter, lfiltic, filtfilt

    Notes
    -----
    A linear filter with order m has a state space representation (A, B, C, D),
    for which the output y of the filter can be expressed as::

        z(n+1) = A*z(n) + B*x(n)
        y(n)   = C*z(n) + D*x(n)

    where z(n) is a vector of length m, A has shape (m, m), B has shape
    (m, 1), C has shape (1, m) and D has shape (1, 1) (assuming x(n) is
    a scalar).  lfilter_zi solves::

        zi = A*zi + B

    In other words, it finds the initial condition for which the response
    to an input of all ones is a constant.

    Given the filter coefficients `a` and `b`, the state space matrices
    for the transposed direct form II implementation of the linear filter,
    which is the implementation used by scipy.signal.lfilter, are::

        A = scipy.linalg.companion(a).T
        B = b[1:] - a[1:]*b[0]

    assuming `a[0]` is 1.0; if `a[0]` is not 1, `a` and `b` are first
    divided by a[0].

    Examples
    --------
    The following code creates a lowpass Butterworth filter. Then it
    applies that filter to an array whose values are all 1.0; the
    output is also all 1.0, as expected for a lowpass filter.  If the
    `zi` argument of `lfilter` had not been given, the output would have
    shown the transient signal.

    >>> from numpy import array, ones
    >>> from scipy.signal import lfilter, lfilter_zi, butter
    >>> b, a = butter(5, 0.25)
    >>> zi = lfilter_zi(b, a)
    >>> y, zo = lfilter(b, a, ones(10), zi=zi)
    >>> y
    array([1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.])

    Another example:

    >>> x = array([0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0])
    >>> y, zf = lfilter(b, a, x, zi=zi*x[0])
    >>> y
    array([ 0.5       ,  0.5       ,  0.5       ,  0.49836039,  0.48610528,
        0.44399389,  0.35505241])

    Note that the `zi` argument to `lfilter` was computed using
    `lfilter_zi` and scaled by `x[0]`.  Then the output `y` has no
    transient until the input drops from 0.5 to 0.0.

    """

    # FIXME: Can this function be replaced with an appropriate
    # use of lfiltic?  For example, when b,a = butter(N,Wn),
    #    lfiltic(b, a, y=numpy.ones_like(a), x=numpy.ones_like(b)).
    #

    # We could use scipy.signal.normalize, but it uses warnings in
    # cases where a ValueError is more appropriate, and it allows
    # b to be 2D.
    b = np.atleast_1d(b)
    if b.ndim != 1:
        raise ValueError("Numerator b must be 1-D.")
    a = np.atleast_1d(a)
    if a.ndim != 1:
        raise ValueError("Denominator a must be 1-D.")

    while len(a) > 1 and a[0] == 0.0:
        a = a[1:]
    if a.size < 1:
        raise ValueError("There must be at least one nonzero `a` coefficient.")

    if a[0] != 1.0:
        # Normalize the coefficients so a[0] == 1.
        b = b / a[0]
        a = a / a[0]

    n = max(len(a), len(b))

    # Pad a or b with zeros so they are the same length.
    if len(a) < n:
        a = np.r_[a, np.zeros(n - len(a))]
    elif len(b) < n:
        b = np.r_[b, np.zeros(n - len(b))]

    IminusA = np.eye(n - 1) - linalg.companion(a).T
    B = b[1:] - a[1:] * b[0]
    # Solve zi = A*zi + B
    zi = np.linalg.solve(IminusA, B)

    # For future reference: we could also use the following
    # explicit formulas to solve the linear system:
    #
    # zi = np.zeros(n - 1)
    # zi[0] = B.sum() / IminusA[:,0].sum()
    # asum = 1.0
    # csum = 0.0
    # for k in range(1,n-1):
    #     asum += a[k]
    #     csum += b[k] - a[k]*b[0]
    #     zi[k] = asum*zi[0] - csum

    return zi


def sosfilt_zi(sos):
    """
    Construct initial conditions for sosfilt for step response steady-state.

    Compute an initial state `zi` for the `sosfilt` function that corresponds
    to the steady state of the step response.

    A typical use of this function is to set the initial state so that the
    output of the filter starts at the same value as the first element of
    the signal to be filtered.

    Parameters
    ----------
    sos : array_like
        Array of second-order filter coefficients, must have shape
        ``(n_sections, 6)``. See `sosfilt` for the SOS filter format
        specification.

    Returns
    -------
    zi : ndarray
        Initial conditions suitable for use with ``sosfilt``, shape
        ``(n_sections, 2)``.

    See Also
    --------
    sosfilt, zpk2sos

    Notes
    -----
    .. versionadded:: 0.16.0

    Examples
    --------
    Filter a rectangular pulse that begins at time 0, with and without
    the use of the `zi` argument of `scipy.signal.sosfilt`.

    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt

    >>> sos = signal.butter(9, 0.125, output='sos')
    >>> zi = signal.sosfilt_zi(sos)
    >>> x = (np.arange(250) < 100).astype(int)
    >>> f1 = signal.sosfilt(sos, x)
    >>> f2, zo = signal.sosfilt(sos, x, zi=zi)

    >>> plt.plot(x, 'k--', label='x')
    >>> plt.plot(f1, 'b', alpha=0.5, linewidth=2, label='filtered')
    >>> plt.plot(f2, 'g', alpha=0.25, linewidth=4, label='filtered with zi')
    >>> plt.legend(loc='best')
    >>> plt.show()

    """
    sos = np.asarray(sos)
    if sos.ndim != 2 or sos.shape[1] != 6:
        raise ValueError('sos must be shape (n_sections, 6)')

    n_sections = sos.shape[0]
    zi = np.empty((n_sections, 2))
    scale = 1.0
    for section in range(n_sections):
        b = sos[section, :3]
        a = sos[section, 3:]
        zi[section] = scale * lfilter_zi(b, a)
        # If H(z) = B(z)/A(z) is this section's transfer function, then
        # b.sum()/a.sum() is H(1), the gain at omega=0.  That's the steady
        # state value of this section's step response.
        scale *= b.sum() / a.sum()

    return zi


def _filtfilt_gust(b, a, x, axis=-1, irlen=None):
    """Forward-backward IIR filter that uses Gustafsson's method.

    Apply the IIR filter defined by `(b,a)` to `x` twice, first forward
    then backward, using Gustafsson's initial conditions [1]_.

    Let ``y_fb`` be the result of filtering first forward and then backward,
    and let ``y_bf`` be the result of filtering first backward then forward.
    Gustafsson's method is to compute initial conditions for the forward
    pass and the backward pass such that ``y_fb == y_bf``.

    Parameters
    ----------
    b : scalar or 1-D ndarray
        Numerator coefficients of the filter.
    a : scalar or 1-D ndarray
        Denominator coefficients of the filter.
    x : ndarray
        Data to be filtered.
    axis : int, optional
        Axis of `x` to be filtered.  Default is -1.
    irlen : int or None, optional
        The length of the nonnegligible part of the impulse response.
        If `irlen` is None, or if the length of the signal is less than
        ``2 * irlen``, then no part of the impulse response is ignored.

    Returns
    -------
    y : ndarray
        The filtered data.
    x0 : ndarray
        Initial condition for the forward filter.
    x1 : ndarray
        Initial condition for the backward filter.

    Notes
    -----
    Typically the return values `x0` and `x1` are not needed by the
    caller.  The intended use of these return values is in unit tests.

    References
    ----------
    .. [1] F. Gustaffson. Determining the initial states in forward-backward
           filtering. Transactions on Signal Processing, 46(4):988-992, 1996.

    """
    # In the comments, "Gustafsson's paper" and [1] refer to the
    # paper referenced in the docstring.

    b = np.atleast_1d(b)
    a = np.atleast_1d(a)

    order = max(len(b), len(a)) - 1
    if order == 0:
        # The filter is just scalar multiplication, with no state.
        scale = (b[0] / a[0])**2
        y = scale * x
        return y, np.array([]), np.array([])

    if axis != -1 or axis != x.ndim - 1:
        # Move the axis containing the data to the end.
        x = np.swapaxes(x, axis, x.ndim - 1)

    # n is the number of samples in the data to be filtered.
    n = x.shape[-1]

    if irlen is None or n <= 2*irlen:
        m = n
    else:
        m = irlen

    # Create Obs, the observability matrix (called O in the paper).
    # This matrix can be interpreted as the operator that propagates
    # an arbitrary initial state to the output, assuming the input is
    # zero.
    # In Gustafsson's paper, the forward and backward filters are not
    # necessarily the same, so he has both O_f and O_b.  We use the same
    # filter in both directions, so we only need O. The same comment
    # applies to S below.
    Obs = np.zeros((m, order))
    zi = np.zeros(order)
    zi[0] = 1
    Obs[:, 0] = lfilter(b, a, np.zeros(m), zi=zi)[0]
    for k in range(1, order):
        Obs[k:, k] = Obs[:-k, 0]

    # Obsr is O^R (Gustafsson's notation for row-reversed O)
    Obsr = Obs[::-1]

    # Create S.  S is the matrix that applies the filter to the reversed
    # propagated initial conditions.  That is,
    #     out = S.dot(zi)
    # is the same as
    #     tmp, _ = lfilter(b, a, zeros(), zi=zi)  # Propagate ICs.
    #     out = lfilter(b, a, tmp[::-1])          # Reverse and filter.

    # Equations (5) & (6) of [1]
    S = lfilter(b, a, Obs[::-1], axis=0)

    # Sr is S^R (row-reversed S)
    Sr = S[::-1]

    # M is [(S^R - O), (O^R - S)]
    if m == n:
        M = np.hstack((Sr - Obs, Obsr - S))
    else:
        # Matrix described in section IV of [1].
        M = np.zeros((2*m, 2*order))
        M[:m, :order] = Sr - Obs
        M[m:, order:] = Obsr - S

    # Naive forward-backward and backward-forward filters.
    # These have large transients because the filters use zero initial
    # conditions.
    y_f = lfilter(b, a, x)
    y_fb = lfilter(b, a, y_f[..., ::-1])[..., ::-1]

    y_b = lfilter(b, a, x[..., ::-1])[..., ::-1]
    y_bf = lfilter(b, a, y_b)

    delta_y_bf_fb = y_bf - y_fb
    if m == n:
        delta = delta_y_bf_fb
    else:
        start_m = delta_y_bf_fb[..., :m]
        end_m = delta_y_bf_fb[..., -m:]
        delta = np.concatenate((start_m, end_m), axis=-1)

    # ic_opt holds the "optimal" initial conditions.
    # The following code computes the result shown in the formula
    # of the paper between equations (6) and (7).
    if delta.ndim == 1:
        ic_opt = linalg.lstsq(M, delta)[0]
    else:
        # Reshape delta so it can be used as an array of multiple
        # right-hand-sides in linalg.lstsq.
        delta2d = delta.reshape(-1, delta.shape[-1]).T
        ic_opt0 = linalg.lstsq(M, delta2d)[0].T
        ic_opt = ic_opt0.reshape(delta.shape[:-1] + (M.shape[-1],))

    # Now compute the filtered signal using equation (7) of [1].
    # First, form [S^R, O^R] and call it W.
    if m == n:
        W = np.hstack((Sr, Obsr))
    else:
        W = np.zeros((2*m, 2*order))
        W[:m, :order] = Sr
        W[m:, order:] = Obsr

    # Equation (7) of [1] says
    #     Y_fb^opt = Y_fb^0 + W * [x_0^opt; x_{N-1}^opt]
    # `wic` is (almost) the product on the right.
    # W has shape (m, 2*order), and ic_opt has shape (..., 2*order),
    # so we can't use W.dot(ic_opt).  Instead, we dot ic_opt with W.T,
    # so wic has shape (..., m).
    wic = ic_opt.dot(W.T)

    # `wic` is "almost" the product of W and the optimal ICs in equation
    # (7)--if we're using a truncated impulse response (m < n), `wic`
    # contains only the adjustments required for the ends of the signal.
    # Here we form y_opt, taking this into account if necessary.
    y_opt = y_fb
    if m == n:
        y_opt += wic
    else:
        y_opt[..., :m] += wic[..., :m]
        y_opt[..., -m:] += wic[..., -m:]

    x0 = ic_opt[..., :order]
    x1 = ic_opt[..., -order:]
    if axis != -1 or axis != x.ndim - 1:
        # Restore the data axis to its original position.
        x0 = np.swapaxes(x0, axis, x.ndim - 1)
        x1 = np.swapaxes(x1, axis, x.ndim - 1)
        y_opt = np.swapaxes(y_opt, axis, x.ndim - 1)

    return y_opt, x0, x1


def filtfilt(b, a, x, axis=-1, padtype='odd', padlen=None, method='pad',
             irlen=None):
    """
    Apply a digital filter forward and backward to a signal.

    This function applies a linear digital filter twice, once forward and
    once backwards.  The combined filter has zero phase and a filter order
    twice that of the original.

    The function provides options for handling the edges of the signal.

    Parameters
    ----------
    b : (N,) array_like
        The numerator coefficient vector of the filter.
    a : (N,) array_like
        The denominator coefficient vector of the filter.  If ``a[0]``
        is not 1, then both `a` and `b` are normalized by ``a[0]``.
    x : array_like
        The array of data to be filtered.
    axis : int, optional
        The axis of `x` to which the filter is applied.
        Default is -1.
    padtype : str or None, optional
        Must be 'odd', 'even', 'constant', or None.  This determines the
        type of extension to use for the padded signal to which the filter
        is applied.  If `padtype` is None, no padding is used.  The default
        is 'odd'.
    padlen : int or None, optional
        The number of elements by which to extend `x` at both ends of
        `axis` before applying the filter.  This value must be less than
        ``x.shape[axis] - 1``.  ``padlen=0`` implies no padding.
        The default value is ``3 * max(len(a), len(b))``.
    method : str, optional
        Determines the method for handling the edges of the signal, either
        "pad" or "gust".  When `method` is "pad", the signal is padded; the
        type of padding is determined by `padtype` and `padlen`, and `irlen`
        is ignored.  When `method` is "gust", Gustafsson's method is used,
        and `padtype` and `padlen` are ignored.
    irlen : int or None, optional
        When `method` is "gust", `irlen` specifies the length of the
        impulse response of the filter.  If `irlen` is None, no part
        of the impulse response is ignored.  For a long signal, specifying
        `irlen` can significantly improve the performance of the filter.

    Returns
    -------
    y : ndarray
        The filtered output with the same shape as `x`.

    See Also
    --------
    sosfiltfilt, lfilter_zi, lfilter, lfiltic, savgol_filter, sosfilt

    Notes
    -----
    When `method` is "pad", the function pads the data along the given axis
    in one of three ways: odd, even or constant.  The odd and even extensions
    have the corresponding symmetry about the end point of the data.  The
    constant extension extends the data with the values at the end points. On
    both the forward and backward passes, the initial condition of the
    filter is found by using `lfilter_zi` and scaling it by the end point of
    the extended data.

    When `method` is "gust", Gustafsson's method [1]_ is used.  Initial
    conditions are chosen for the forward and backward passes so that the
    forward-backward filter gives the same result as the backward-forward
    filter.

    The option to use Gustaffson's method was added in scipy version 0.16.0.

    References
    ----------
    .. [1] F. Gustaffson, "Determining the initial states in forward-backward
           filtering", Transactions on Signal Processing, Vol. 46, pp. 988-992,
           1996.

    Examples
    --------
    The examples will use several functions from `scipy.signal`.

    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt

    First we create a one second signal that is the sum of two pure sine
    waves, with frequencies 5 Hz and 250 Hz, sampled at 2000 Hz.

    >>> t = np.linspace(0, 1.0, 2001)
    >>> xlow = np.sin(2 * np.pi * 5 * t)
    >>> xhigh = np.sin(2 * np.pi * 250 * t)
    >>> x = xlow + xhigh

    Now create a lowpass Butterworth filter with a cutoff of 0.125 times
    the Nyquist frequency, or 125 Hz, and apply it to ``x`` with `filtfilt`.
    The result should be approximately ``xlow``, with no phase shift.

    >>> b, a = signal.butter(8, 0.125)
    >>> y = signal.filtfilt(b, a, x, padlen=150)
    >>> np.abs(y - xlow).max()
    9.1086182074789912e-06

    We get a fairly clean result for this artificial example because
    the odd extension is exact, and with the moderately long padding,
    the filter's transients have dissipated by the time the actual data
    is reached.  In general, transient effects at the edges are
    unavoidable.

    The following example demonstrates the option ``method="gust"``.

    First, create a filter.

    >>> b, a = signal.ellip(4, 0.01, 120, 0.125)  # Filter to be applied.
    >>> np.random.seed(123456)

    `sig` is a random input signal to be filtered.

    >>> n = 60
    >>> sig = np.random.randn(n)**3 + 3*np.random.randn(n).cumsum()

    Apply `filtfilt` to `sig`, once using the Gustafsson method, and
    once using padding, and plot the results for comparison.

    >>> fgust = signal.filtfilt(b, a, sig, method="gust")
    >>> fpad = signal.filtfilt(b, a, sig, padlen=50)
    >>> plt.plot(sig, 'k-', label='input')
    >>> plt.plot(fgust, 'b-', linewidth=4, label='gust')
    >>> plt.plot(fpad, 'c-', linewidth=1.5, label='pad')
    >>> plt.legend(loc='best')
    >>> plt.show()

    The `irlen` argument can be used to improve the performance
    of Gustafsson's method.

    Estimate the impulse response length of the filter.

    >>> z, p, k = signal.tf2zpk(b, a)
    >>> eps = 1e-9
    >>> r = np.max(np.abs(p))
    >>> approx_impulse_len = int(np.ceil(np.log(eps) / np.log(r)))
    >>> approx_impulse_len
    137

    Apply the filter to a longer signal, with and without the `irlen`
    argument.  The difference between `y1` and `y2` is small.  For long
    signals, using `irlen` gives a significant performance improvement.

    >>> x = np.random.randn(5000)
    >>> y1 = signal.filtfilt(b, a, x, method='gust')
    >>> y2 = signal.filtfilt(b, a, x, method='gust', irlen=approx_impulse_len)
    >>> print(np.max(np.abs(y1 - y2)))
    1.80056858312e-10

    """
    b = np.atleast_1d(b)
    a = np.atleast_1d(a)
    x = np.asarray(x)

    if method not in ["pad", "gust"]:
        raise ValueError("method must be 'pad' or 'gust'.")

    if method == "gust":
        y, z1, z2 = _filtfilt_gust(b, a, x, axis=axis, irlen=irlen)
        return y

    # method == "pad"
    edge, ext = _validate_pad(padtype, padlen, x, axis,
                              ntaps=max(len(a), len(b)))

    # Get the steady state of the filter's step response.
    zi = lfilter_zi(b, a)

    # Reshape zi and create x0 so that zi*x0 broadcasts
    # to the correct value for the 'zi' keyword argument
    # to lfilter.
    zi_shape = [1] * x.ndim
    zi_shape[axis] = zi.size
    zi = np.reshape(zi, zi_shape)
    x0 = axis_slice(ext, stop=1, axis=axis)

    # Forward filter.
    (y, zf) = lfilter(b, a, ext, axis=axis, zi=zi * x0)

    # Backward filter.
    # Create y0 so zi*y0 broadcasts appropriately.
    y0 = axis_slice(y, start=-1, axis=axis)
    (y, zf) = lfilter(b, a, axis_reverse(y, axis=axis), axis=axis, zi=zi * y0)

    # Reverse y.
    y = axis_reverse(y, axis=axis)

    if edge > 0:
        # Slice the actual signal from the extended signal.
        y = axis_slice(y, start=edge, stop=-edge, axis=axis)

    return y


def _validate_pad(padtype, padlen, x, axis, ntaps):
    """Helper to validate padding for filtfilt"""
    if padtype not in ['even', 'odd', 'constant', None]:
        raise ValueError(("Unknown value '%s' given to padtype.  padtype "
                          "must be 'even', 'odd', 'constant', or None.") %
                         padtype)

    if padtype is None:
        padlen = 0

    if padlen is None:
        # Original padding; preserved for backwards compatibility.
        edge = ntaps * 3
    else:
        edge = padlen

    # x's 'axis' dimension must be bigger than edge.
    if x.shape[axis] <= edge:
        raise ValueError("The length of the input vector x must be greater than "
                         "padlen, which is %d." % edge)

    if padtype is not None and edge > 0:
        # Make an extension of length `edge` at each
        # end of the input array.
        if padtype == 'even':
            ext = even_ext(x, edge, axis=axis)
        elif padtype == 'odd':
            ext = odd_ext(x, edge, axis=axis)
        else:
            ext = const_ext(x, edge, axis=axis)
    else:
        ext = x
    return edge, ext


def sosfilt(sos, x, axis=-1, zi=None):
    """
    Filter data along one dimension using cascaded second-order sections.

    Filter a data sequence, `x`, using a digital IIR filter defined by
    `sos`. This is implemented by performing `lfilter` for each
    second-order section.  See `lfilter` for details.

    Parameters
    ----------
    sos : array_like
        Array of second-order filter coefficients, must have shape
        ``(n_sections, 6)``. Each row corresponds to a second-order
        section, with the first three columns providing the numerator
        coefficients and the last three providing the denominator
        coefficients.
    x : array_like
        An N-dimensional input array.
    axis : int, optional
        The axis of the input data array along which to apply the
        linear filter. The filter is applied to each subarray along
        this axis.  Default is -1.
    zi : array_like, optional
        Initial conditions for the cascaded filter delays.  It is a (at
        least 2D) vector of shape ``(n_sections, ..., 2, ...)``, where
        ``..., 2, ...`` denotes the shape of `x`, but with ``x.shape[axis]``
        replaced by 2.  If `zi` is None or is not given then initial rest
        (i.e. all zeros) is assumed.
        Note that these initial conditions are *not* the same as the initial
        conditions given by `lfiltic` or `lfilter_zi`.

    Returns
    -------
    y : ndarray
        The output of the digital filter.
    zf : ndarray, optional
        If `zi` is None, this is not returned, otherwise, `zf` holds the
        final filter delay values.

    See Also
    --------
    zpk2sos, sos2zpk, sosfilt_zi, sosfiltfilt, sosfreqz

    Notes
    -----
    The filter function is implemented as a series of second-order filters
    with direct-form II transposed structure. It is designed to minimize
    numerical precision errors for high-order filters.

    .. versionadded:: 0.16.0

    Examples
    --------
    Plot a 13th-order filter's impulse response using both `lfilter` and
    `sosfilt`, showing the instability that results from trying to do a
    13th-order filter in a single stage (the numerical error pushes some poles
    outside of the unit circle):

    >>> import matplotlib.pyplot as plt
    >>> from scipy import signal
    >>> b, a = signal.ellip(13, 0.009, 80, 0.05, output='ba')
    >>> sos = signal.ellip(13, 0.009, 80, 0.05, output='sos')
    >>> x = signal.unit_impulse(700)
    >>> y_tf = signal.lfilter(b, a, x)
    >>> y_sos = signal.sosfilt(sos, x)
    >>> plt.plot(y_tf, 'r', label='TF')
    >>> plt.plot(y_sos, 'k', label='SOS')
    >>> plt.legend(loc='best')
    >>> plt.show()

    """
    x = np.asarray(x)
    sos, n_sections = _validate_sos(sos)
    use_zi = zi is not None
    if use_zi:
        zi = np.asarray(zi)
        x_zi_shape = list(x.shape)
        x_zi_shape[axis] = 2
        x_zi_shape = tuple([n_sections] + x_zi_shape)
        if zi.shape != x_zi_shape:
            raise ValueError('Invalid zi shape. With axis=%r, an input with '
                             'shape %r, and an sos array with %d sections, zi '
                             'must have shape %r, got %r.' %
                             (axis, x.shape, n_sections, x_zi_shape, zi.shape))
        zf = zeros_like(zi)

    for section in range(n_sections):
        if use_zi:
            x, zf[section] = lfilter(sos[section, :3], sos[section, 3:],
                                     x, axis, zi=zi[section])
        else:
            x = lfilter(sos[section, :3], sos[section, 3:], x, axis)
    out = (x, zf) if use_zi else x
    return out


def sosfiltfilt(sos, x, axis=-1, padtype='odd', padlen=None):
    """
    A forward-backward digital filter using cascaded second-order sections.

    See `filtfilt` for more complete information about this method.

    Parameters
    ----------
    sos : array_like
        Array of second-order filter coefficients, must have shape
        ``(n_sections, 6)``. Each row corresponds to a second-order
        section, with the first three columns providing the numerator
        coefficients and the last three providing the denominator
        coefficients.
    x : array_like
        The array of data to be filtered.
    axis : int, optional
        The axis of `x` to which the filter is applied.
        Default is -1.
    padtype : str or None, optional
        Must be 'odd', 'even', 'constant', or None.  This determines the
        type of extension to use for the padded signal to which the filter
        is applied.  If `padtype` is None, no padding is used.  The default
        is 'odd'.
    padlen : int or None, optional
        The number of elements by which to extend `x` at both ends of
        `axis` before applying the filter.  This value must be less than
        ``x.shape[axis] - 1``.  ``padlen=0`` implies no padding.
        The default value is::

            3 * (2 * len(sos) + 1 - min((sos[:, 2] == 0).sum(),
                                        (sos[:, 5] == 0).sum()))

        The extra subtraction at the end attempts to compensate for poles
        and zeros at the origin (e.g. for odd-order filters) to yield
        equivalent estimates of `padlen` to those of `filtfilt` for
        second-order section filters built with `scipy.signal` functions.

    Returns
    -------
    y : ndarray
        The filtered output with the same shape as `x`.

    See Also
    --------
    filtfilt, sosfilt, sosfilt_zi, sosfreqz

    Notes
    -----
    .. versionadded:: 0.18.0

    Examples
    --------
    >>> from scipy.signal import sosfiltfilt, butter
    >>> import matplotlib.pyplot as plt

    Create an interesting signal to filter.

    >>> n = 201
    >>> t = np.linspace(0, 1, n)
    >>> np.random.seed(123)
    >>> x = 1 + (t < 0.5) - 0.25*t**2 + 0.05*np.random.randn(n)

    Create a lowpass Butterworth filter, and use it to filter `x`.

    >>> sos = butter(4, 0.125, output='sos')
    >>> y = sosfiltfilt(sos, x)

    For comparison, apply an 8th order filter using `sosfilt`.  The filter
    is initialized using the mean of the first four values of `x`.

    >>> from scipy.signal import sosfilt, sosfilt_zi
    >>> sos8 = butter(8, 0.125, output='sos')
    >>> zi = x[:4].mean() * sosfilt_zi(sos8)
    >>> y2, zo = sosfilt(sos8, x, zi=zi)

    Plot the results.  Note that the phase of `y` matches the input, while
    `y2` has a significant phase delay.

    >>> plt.plot(t, x, alpha=0.5, label='x(t)')
    >>> plt.plot(t, y, label='y(t)')
    >>> plt.plot(t, y2, label='y2(t)')
    >>> plt.legend(framealpha=1, shadow=True)
    >>> plt.grid(alpha=0.25)
    >>> plt.xlabel('t')
    >>> plt.show()

    """
    sos, n_sections = _validate_sos(sos)

    # `method` is "pad"...
    ntaps = 2 * n_sections + 1
    ntaps -= min((sos[:, 2] == 0).sum(), (sos[:, 5] == 0).sum())
    edge, ext = _validate_pad(padtype, padlen, x, axis,
                              ntaps=ntaps)

    # These steps follow the same form as filtfilt with modifications
    zi = sosfilt_zi(sos)  # shape (n_sections, 2) --> (n_sections, ..., 2, ...)
    zi_shape = [1] * x.ndim
    zi_shape[axis] = 2
    zi.shape = [n_sections] + zi_shape
    x_0 = axis_slice(ext, stop=1, axis=axis)
    (y, zf) = sosfilt(sos, ext, axis=axis, zi=zi * x_0)
    y_0 = axis_slice(y, start=-1, axis=axis)
    (y, zf) = sosfilt(sos, axis_reverse(y, axis=axis), axis=axis, zi=zi * y_0)
    y = axis_reverse(y, axis=axis)
    if edge > 0:
        y = axis_slice(y, start=edge, stop=-edge, axis=axis)
    return y


def decimate(x, q, n=None, ftype='iir', axis=-1, zero_phase=True):
    """
    Downsample the signal after applying an anti-aliasing filter.

    By default, an order 8 Chebyshev type I filter is used. A 30 point FIR
    filter with Hamming window is used if `ftype` is 'fir'.

    Parameters
    ----------
    x : array_like
        The signal to be downsampled, as an N-dimensional array.
    q : int
        The downsampling factor. When using IIR downsampling, it is recommended
        to call `decimate` multiple times for downsampling factors higher than
        13.
    n : int, optional
        The order of the filter (1 less than the length for 'fir'). Defaults to
        8 for 'iir' and 20 times the downsampling factor for 'fir'.
    ftype : str {'iir', 'fir'} or ``dlti`` instance, optional
        If 'iir' or 'fir', specifies the type of lowpass filter. If an instance
        of an `dlti` object, uses that object to filter before downsampling.
    axis : int, optional
        The axis along which to decimate.
    zero_phase : bool, optional
        Prevent phase shift by filtering with `filtfilt` instead of `lfilter`
        when using an IIR filter, and shifting the outputs back by the filter's
        group delay when using an FIR filter. The default value of ``True`` is
        recommended, since a phase shift is generally not desired.

        .. versionadded:: 0.18.0

    Returns
    -------
    y : ndarray
        The down-sampled signal.

    See Also
    --------
    resample : Resample up or down using the FFT method.
    resample_poly : Resample using polyphase filtering and an FIR filter.

    Notes
    -----
    The ``zero_phase`` keyword was added in 0.18.0.
    The possibility to use instances of ``dlti`` as ``ftype`` was added in
    0.18.0.
    """

    x = asarray(x)
    q = operator.index(q)

    if n is not None:
        n = operator.index(n)

    if ftype == 'fir':
        if n is None:
            half_len = 10 * q  # reasonable cutoff for our sinc-like function
            n = 2 * half_len
        b, a = firwin(n+1, 1. / q, window='hamming'), 1.
    elif ftype == 'iir':
        if n is None:
            n = 8
        system = dlti(*cheby1(n, 0.05, 0.8 / q))
        b, a = system.num, system.den
    elif isinstance(ftype, dlti):
        system = ftype._as_tf()  # Avoids copying if already in TF form
        b, a = system.num, system.den
    else:
        raise ValueError('invalid ftype')

    sl = [slice(None)] * x.ndim
    a = np.asarray(a)

    if a.size == 1:  # FIR case
        b = b / a
        if zero_phase:
            y = resample_poly(x, 1, q, axis=axis, window=b)
        else:
            # upfirdn is generally faster than lfilter by a factor equal to the
            # downsampling factor, since it only calculates the needed outputs
            n_out = x.shape[axis] // q + bool(x.shape[axis] % q)
            y = upfirdn(b, x, up=1, down=q, axis=axis)
            sl[axis] = slice(None, n_out, None)

    else:  # IIR case
        if zero_phase:
            y = filtfilt(b, a, x, axis=axis)
        else:
            y = lfilter(b, a, x, axis=axis)
        sl[axis] = slice(None, None, q)

    return y[tuple(sl)]
