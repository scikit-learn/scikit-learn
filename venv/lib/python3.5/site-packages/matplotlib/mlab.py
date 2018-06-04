"""

Numerical python functions written for compatibility with MATLAB
commands with the same names.

MATLAB compatible functions
---------------------------

:func:`cohere`
    Coherence (normalized cross spectral density)

:func:`csd`
    Cross spectral density using Welch's average periodogram

:func:`detrend`
    Remove the mean or best fit line from an array

:func:`find`
    Return the indices where some condition is true;
    numpy.nonzero is similar but more general.

:func:`griddata`
    Interpolate irregularly distributed data to a
    regular grid.

:func:`prctile`
    Find the percentiles of a sequence

:func:`prepca`
    Principal Component Analysis

:func:`psd`
    Power spectral density using Welch's average periodogram

:func:`rk4`
    A 4th order runge kutta integrator for 1D or ND systems

:func:`specgram`
    Spectrogram (spectrum over segments of time)

Miscellaneous functions
-----------------------

Functions that don't exist in MATLAB, but are useful anyway:

:func:`cohere_pairs`
    Coherence over all pairs.  This is not a MATLAB function, but we
    compute coherence a lot in my lab, and we compute it for a lot of
    pairs.  This function is optimized to do this efficiently by
    caching the direct FFTs.

:func:`rk4`
    A 4th order Runge-Kutta ODE integrator in case you ever find
    yourself stranded without scipy (and the far superior
    scipy.integrate tools)

:func:`contiguous_regions`
    Return the indices of the regions spanned by some logical mask

:func:`cross_from_below`
    Return the indices where a 1D array crosses a threshold from below

:func:`cross_from_above`
    Return the indices where a 1D array crosses a threshold from above

:func:`complex_spectrum`
    Return the complex-valued frequency spectrum of a signal

:func:`magnitude_spectrum`
    Return the magnitude of the frequency spectrum of a signal

:func:`angle_spectrum`
    Return the angle (wrapped phase) of the frequency spectrum of a signal

:func:`phase_spectrum`
    Return the phase (unwrapped angle) of the frequency spectrum of a signal

:func:`detrend_mean`
    Remove the mean from a line.

:func:`demean`
    Remove the mean from a line. This function is the same as
    :func:`detrend_mean` except for the default *axis*.

:func:`detrend_linear`
    Remove the best fit line from a line.

:func:`detrend_none`
    Return the original line.

:func:`stride_windows`
    Get all windows in an array in a memory-efficient manner

:func:`stride_repeat`
    Repeat an array in a memory-efficient manner

:func:`apply_window`
    Apply a window along a given axis


record array helper functions
-----------------------------

A collection of helper methods for numpyrecord arrays

.. _htmlonly:

    See :ref:`misc-examples-index`

:func:`rec2txt`
    Pretty print a record array

:func:`rec2csv`
    Store record array in CSV file

:func:`csv2rec`
    Import record array from CSV file with type inspection

:func:`rec_append_fields`
    Adds  field(s)/array(s) to record array

:func:`rec_drop_fields`
    Drop fields from record array

:func:`rec_join`
    Join two record arrays on sequence of fields

:func:`recs_join`
    A simple join of multiple recarrays using a single column as a key

:func:`rec_groupby`
    Summarize data by groups (similar to SQL GROUP BY)

:func:`rec_summarize`
    Helper code to filter rec array fields into new fields

For the rec viewer functions(e rec2csv), there are a bunch of Format
objects you can pass into the functions that will do things like color
negative values red, set percent formatting and scaling, etc.

Example usage::

    r = csv2rec('somefile.csv', checkrows=0)

    formatd = dict(
        weight = FormatFloat(2),
        change = FormatPercent(2),
        cost   = FormatThousands(2),
        )


    rec2excel(r, 'test.xls', formatd=formatd)
    rec2csv(r, 'test.csv', formatd=formatd)
    scroll = rec2gtk(r, formatd=formatd)

    win = gtk.Window()
    win.set_size_request(600,800)
    win.add(scroll)
    win.show_all()
    gtk.main()


"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves import map, xrange, zip

import copy
import csv
import operator
import os
import warnings

import numpy as np

import matplotlib.cbook as cbook
from matplotlib import docstring
from matplotlib.path import Path
import math


if six.PY3:
    long = int


@cbook.deprecated("2.2", alternative='numpy.logspace or numpy.geomspace')
def logspace(xmin, xmax, N):
    '''
    Return N values logarithmically spaced between xmin and xmax.

    '''
    return np.exp(np.linspace(np.log(xmin), np.log(xmax), N))


def window_hanning(x):
    '''
    Return x times the hanning window of len(x).

    See Also
    --------
    :func:`window_none`
        :func:`window_none` is another window algorithm.
    '''
    return np.hanning(len(x))*x


def window_none(x):
    '''
    No window function; simply return x.

    See Also
    --------
    :func:`window_hanning`
        :func:`window_hanning` is another window algorithm.
    '''
    return x


def apply_window(x, window, axis=0, return_window=None):
    '''
    Apply the given window to the given 1D or 2D array along the given axis.

    Parameters
    ----------
    x : 1D or 2D array or sequence
        Array or sequence containing the data.

    window : function or array.
        Either a function to generate a window or an array with length
        *x*.shape[*axis*]

    axis : integer
        The axis over which to do the repetition.
        Must be 0 or 1.  The default is 0

    return_window : bool
        If true, also return the 1D values of the window that was applied
    '''
    x = np.asarray(x)

    if x.ndim < 1 or x.ndim > 2:
        raise ValueError('only 1D or 2D arrays can be used')
    if axis+1 > x.ndim:
        raise ValueError('axis(=%s) out of bounds' % axis)

    xshape = list(x.shape)
    xshapetarg = xshape.pop(axis)

    if cbook.iterable(window):
        if len(window) != xshapetarg:
            raise ValueError('The len(window) must be the same as the shape '
                             'of x for the chosen axis')
        windowVals = window
    else:
        windowVals = window(np.ones(xshapetarg, dtype=x.dtype))

    if x.ndim == 1:
        if return_window:
            return windowVals * x, windowVals
        else:
            return windowVals * x

    xshapeother = xshape.pop()

    otheraxis = (axis+1) % 2

    windowValsRep = stride_repeat(windowVals, xshapeother, axis=otheraxis)

    if return_window:
        return windowValsRep * x, windowVals
    else:
        return windowValsRep * x


def detrend(x, key=None, axis=None):
    '''
    Return x with its trend removed.

    Parameters
    ----------
    x : array or sequence
        Array or sequence containing the data.

    key : [ 'default' | 'constant' | 'mean' | 'linear' | 'none'] or function
        Specifies the detrend algorithm to use. 'default' is 'mean', which is
        the same as :func:`detrend_mean`. 'constant' is the same. 'linear' is
        the same as :func:`detrend_linear`. 'none' is the same as
        :func:`detrend_none`. The default is 'mean'. See the corresponding
        functions for more details regarding the algorithms. Can also be a
        function that carries out the detrend operation.

    axis : integer
        The axis along which to do the detrending.

    See Also
    --------
    :func:`detrend_mean`
        :func:`detrend_mean` implements the 'mean' algorithm.

    :func:`detrend_linear`
        :func:`detrend_linear` implements the 'linear' algorithm.

    :func:`detrend_none`
        :func:`detrend_none` implements the 'none' algorithm.
    '''
    if key is None or key in ['constant', 'mean', 'default']:
        return detrend(x, key=detrend_mean, axis=axis)
    elif key == 'linear':
        return detrend(x, key=detrend_linear, axis=axis)
    elif key == 'none':
        return detrend(x, key=detrend_none, axis=axis)
    elif isinstance(key, six.string_types):
        raise ValueError("Unknown value for key %s, must be one of: "
                         "'default', 'constant', 'mean', "
                         "'linear', or a function" % key)

    if not callable(key):
        raise ValueError("Unknown value for key %s, must be one of: "
                         "'default', 'constant', 'mean', "
                         "'linear', or a function" % key)

    x = np.asarray(x)

    if axis is not None and axis+1 > x.ndim:
        raise ValueError('axis(=%s) out of bounds' % axis)

    if (axis is None and x.ndim == 0) or (not axis and x.ndim == 1):
        return key(x)

    # try to use the 'axis' argument if the function supports it,
    # otherwise use apply_along_axis to do it
    try:
        return key(x, axis=axis)
    except TypeError:
        return np.apply_along_axis(key, axis=axis, arr=x)


def demean(x, axis=0):
    '''
    Return x minus its mean along the specified axis.

    Parameters
    ----------
    x : array or sequence
        Array or sequence containing the data
        Can have any dimensionality

    axis : integer
        The axis along which to take the mean.  See numpy.mean for a
        description of this argument.

    See Also
    --------
    :func:`delinear`

    :func:`denone`
        :func:`delinear` and :func:`denone` are other detrend algorithms.

    :func:`detrend_mean`
        This function is the same as :func:`detrend_mean` except for the
        default *axis*.
    '''
    return detrend_mean(x, axis=axis)


def detrend_mean(x, axis=None):
    '''
    Return x minus the mean(x).

    Parameters
    ----------
    x : array or sequence
        Array or sequence containing the data
        Can have any dimensionality

    axis : integer
        The axis along which to take the mean.  See numpy.mean for a
        description of this argument.

    See Also
    --------
    :func:`demean`
        This function is the same as :func:`demean` except for the default
        *axis*.

    :func:`detrend_linear`

    :func:`detrend_none`
        :func:`detrend_linear` and :func:`detrend_none` are other detrend
        algorithms.

    :func:`detrend`
        :func:`detrend` is a wrapper around all the detrend algorithms.
    '''
    x = np.asarray(x)

    if axis is not None and axis+1 > x.ndim:
        raise ValueError('axis(=%s) out of bounds' % axis)

    # short-circuit 0-D array.
    if not x.ndim:
        return np.array(0., dtype=x.dtype)

    # short-circuit simple operations
    if axis == 0 or axis is None or x.ndim <= 1:
        return x - x.mean(axis)

    ind = [slice(None)] * x.ndim
    ind[axis] = np.newaxis
    return x - x.mean(axis)[ind]


def detrend_none(x, axis=None):
    '''
    Return x: no detrending.

    Parameters
    ----------
    x : any object
        An object containing the data

    axis : integer
        This parameter is ignored.
        It is included for compatibility with detrend_mean

    See Also
    --------
    :func:`denone`
        This function is the same as :func:`denone` except for the default
        *axis*, which has no effect.

    :func:`detrend_mean`

    :func:`detrend_linear`
        :func:`detrend_mean` and :func:`detrend_linear` are other detrend
        algorithms.

    :func:`detrend`
        :func:`detrend` is a wrapper around all the detrend algorithms.
    '''
    return x


def detrend_linear(y):
    '''
    Return x minus best fit line; 'linear' detrending.

    Parameters
    ----------
    y : 0-D or 1-D array or sequence
        Array or sequence containing the data

    axis : integer
        The axis along which to take the mean.  See numpy.mean for a
        description of this argument.

    See Also
    --------
    :func:`delinear`
        This function is the same as :func:`delinear` except for the default
        *axis*.

    :func:`detrend_mean`

    :func:`detrend_none`
        :func:`detrend_mean` and :func:`detrend_none` are other detrend
        algorithms.

    :func:`detrend`
        :func:`detrend` is a wrapper around all the detrend algorithms.
    '''
    # This is faster than an algorithm based on linalg.lstsq.
    y = np.asarray(y)

    if y.ndim > 1:
        raise ValueError('y cannot have ndim > 1')

    # short-circuit 0-D array.
    if not y.ndim:
        return np.array(0., dtype=y.dtype)

    x = np.arange(y.size, dtype=float)

    C = np.cov(x, y, bias=1)
    b = C[0, 1]/C[0, 0]

    a = y.mean() - b*x.mean()
    return y - (b*x + a)


def stride_windows(x, n, noverlap=None, axis=0):
    '''
    Get all windows of x with length n as a single array,
    using strides to avoid data duplication.

    .. warning::

        It is not safe to write to the output array.  Multiple
        elements may point to the same piece of memory,
        so modifying one value may change others.

    Parameters
    ----------
    x : 1D array or sequence
        Array or sequence containing the data.

    n : integer
        The number of data points in each window.

    noverlap : integer
        The overlap between adjacent windows.
        Default is 0 (no overlap)

    axis : integer
        The axis along which the windows will run.

    References
    ----------
    `stackoverflow: Rolling window for 1D arrays in Numpy?
    <http://stackoverflow.com/a/6811241>`_
    `stackoverflow: Using strides for an efficient moving average filter
    <http://stackoverflow.com/a/4947453>`_
    '''
    if noverlap is None:
        noverlap = 0

    if noverlap >= n:
        raise ValueError('noverlap must be less than n')
    if n < 1:
        raise ValueError('n cannot be less than 1')

    x = np.asarray(x)

    if x.ndim != 1:
        raise ValueError('only 1-dimensional arrays can be used')
    if n == 1 and noverlap == 0:
        if axis == 0:
            return x[np.newaxis]
        else:
            return x[np.newaxis].transpose()
    if n > x.size:
        raise ValueError('n cannot be greater than the length of x')

    # np.lib.stride_tricks.as_strided easily leads to memory corruption for
    # non integer shape and strides, i.e. noverlap or n. See #3845.
    noverlap = int(noverlap)
    n = int(n)

    step = n - noverlap
    if axis == 0:
        shape = (n, (x.shape[-1]-noverlap)//step)
        strides = (x.strides[0], step*x.strides[0])
    else:
        shape = ((x.shape[-1]-noverlap)//step, n)
        strides = (step*x.strides[0], x.strides[0])
    return np.lib.stride_tricks.as_strided(x, shape=shape, strides=strides)


def stride_repeat(x, n, axis=0):
    '''
    Repeat the values in an array in a memory-efficient manner.  Array x is
    stacked vertically n times.

    .. warning::

        It is not safe to write to the output array.  Multiple
        elements may point to the same piece of memory, so
        modifying one value may change others.

    Parameters
    ----------
    x : 1D array or sequence
        Array or sequence containing the data.

    n : integer
        The number of time to repeat the array.

    axis : integer
        The axis along which the data will run.

    References
    ----------
    `stackoverflow: Repeat NumPy array without replicating data?
    <http://stackoverflow.com/a/5568169>`_
    '''
    if axis not in [0, 1]:
        raise ValueError('axis must be 0 or 1')
    x = np.asarray(x)
    if x.ndim != 1:
        raise ValueError('only 1-dimensional arrays can be used')

    if n == 1:
        if axis == 0:
            return np.atleast_2d(x)
        else:
            return np.atleast_2d(x).T
    if n < 1:
        raise ValueError('n cannot be less than 1')

    # np.lib.stride_tricks.as_strided easily leads to memory corruption for
    # non integer shape and strides, i.e. n. See #3845.
    n = int(n)

    if axis == 0:
        shape = (n, x.size)
        strides = (0, x.strides[0])
    else:
        shape = (x.size, n)
        strides = (x.strides[0], 0)

    return np.lib.stride_tricks.as_strided(x, shape=shape, strides=strides)


def _spectral_helper(x, y=None, NFFT=None, Fs=None, detrend_func=None,
                     window=None, noverlap=None, pad_to=None,
                     sides=None, scale_by_freq=None, mode=None):
    '''
    This is a helper function that implements the commonality between the
    psd, csd, spectrogram and complex, magnitude, angle, and phase spectrums.
    It is *NOT* meant to be used outside of mlab and may change at any time.
    '''
    if y is None:
        # if y is None use x for y
        same_data = True
    else:
        # The checks for if y is x are so that we can use the same function to
        # implement the core of psd(), csd(), and spectrogram() without doing
        # extra calculations.  We return the unaveraged Pxy, freqs, and t.
        same_data = y is x

    if Fs is None:
        Fs = 2
    if noverlap is None:
        noverlap = 0
    if detrend_func is None:
        detrend_func = detrend_none
    if window is None:
        window = window_hanning

    # if NFFT is set to None use the whole signal
    if NFFT is None:
        NFFT = 256

    if mode is None or mode == 'default':
        mode = 'psd'
    elif mode not in ['psd', 'complex', 'magnitude', 'angle', 'phase']:
        raise ValueError("Unknown value for mode %s, must be one of: "
                         "'default', 'psd', 'complex', "
                         "'magnitude', 'angle', 'phase'" % mode)

    if not same_data and mode != 'psd':
        raise ValueError("x and y must be equal if mode is not 'psd'")

    # Make sure we're dealing with a numpy array. If y and x were the same
    # object to start with, keep them that way
    x = np.asarray(x)
    if not same_data:
        y = np.asarray(y)

    if sides is None or sides == 'default':
        if np.iscomplexobj(x):
            sides = 'twosided'
        else:
            sides = 'onesided'
    elif sides not in ['onesided', 'twosided']:
        raise ValueError("Unknown value for sides %s, must be one of: "
                         "'default', 'onesided', or 'twosided'" % sides)

    # zero pad x and y up to NFFT if they are shorter than NFFT
    if len(x) < NFFT:
        n = len(x)
        x = np.resize(x, (NFFT,))
        x[n:] = 0

    if not same_data and len(y) < NFFT:
        n = len(y)
        y = np.resize(y, (NFFT,))
        y[n:] = 0

    if pad_to is None:
        pad_to = NFFT

    if mode != 'psd':
        scale_by_freq = False
    elif scale_by_freq is None:
        scale_by_freq = True

    # For real x, ignore the negative frequencies unless told otherwise
    if sides == 'twosided':
        numFreqs = pad_to
        if pad_to % 2:
            freqcenter = (pad_to - 1)//2 + 1
        else:
            freqcenter = pad_to//2
        scaling_factor = 1.
    elif sides == 'onesided':
        if pad_to % 2:
            numFreqs = (pad_to + 1)//2
        else:
            numFreqs = pad_to//2 + 1
        scaling_factor = 2.

    result = stride_windows(x, NFFT, noverlap, axis=0)
    result = detrend(result, detrend_func, axis=0)
    result, windowVals = apply_window(result, window, axis=0,
                                      return_window=True)
    result = np.fft.fft(result, n=pad_to, axis=0)[:numFreqs, :]
    freqs = np.fft.fftfreq(pad_to, 1/Fs)[:numFreqs]

    if not same_data:
        # if same_data is False, mode must be 'psd'
        resultY = stride_windows(y, NFFT, noverlap)
        resultY = detrend(resultY, detrend_func, axis=0)
        resultY = apply_window(resultY, window, axis=0)
        resultY = np.fft.fft(resultY, n=pad_to, axis=0)[:numFreqs, :]
        result = np.conj(result) * resultY
    elif mode == 'psd':
        result = np.conj(result) * result
    elif mode == 'magnitude':
        result = np.abs(result) / np.abs(windowVals).sum()
    elif mode == 'angle' or mode == 'phase':
        # we unwrap the phase later to handle the onesided vs. twosided case
        result = np.angle(result)
    elif mode == 'complex':
        result /= np.abs(windowVals).sum()

    if mode == 'psd':

        # Also include scaling factors for one-sided densities and dividing by
        # the sampling frequency, if desired. Scale everything, except the DC
        # component and the NFFT/2 component:

        # if we have a even number of frequencies, don't scale NFFT/2
        if not NFFT % 2:
            slc = slice(1, -1, None)
        # if we have an odd number, just don't scale DC
        else:
            slc = slice(1, None, None)

        result[slc] *= scaling_factor

        # MATLAB divides by the sampling frequency so that density function
        # has units of dB/Hz and can be integrated by the plotted frequency
        # values. Perform the same scaling here.
        if scale_by_freq:
            result /= Fs
            # Scale the spectrum by the norm of the window to compensate for
            # windowing loss; see Bendat & Piersol Sec 11.5.2.
            result /= (np.abs(windowVals)**2).sum()
        else:
            # In this case, preserve power in the segment, not amplitude
            result /= np.abs(windowVals).sum()**2

    t = np.arange(NFFT/2, len(x) - NFFT/2 + 1, NFFT - noverlap)/Fs

    if sides == 'twosided':
        # center the frequency range at zero
        freqs = np.concatenate((freqs[freqcenter:], freqs[:freqcenter]))
        result = np.concatenate((result[freqcenter:, :],
                                 result[:freqcenter, :]), 0)
    elif not pad_to % 2:
        # get the last value correctly, it is negative otherwise
        freqs[-1] *= -1

    # we unwrap the phase here to handle the onesided vs. twosided case
    if mode == 'phase':
        result = np.unwrap(result, axis=0)

    return result, freqs, t


def _single_spectrum_helper(x, mode, Fs=None, window=None, pad_to=None,
                            sides=None):
    '''
    This is a helper function that implements the commonality between the
    complex, magnitude, angle, and phase spectrums.
    It is *NOT* meant to be used outside of mlab and may change at any time.
    '''
    if mode is None or mode == 'psd' or mode == 'default':
        raise ValueError('_single_spectrum_helper does not work with %s mode'
                         % mode)

    if pad_to is None:
        pad_to = len(x)

    spec, freqs, _ = _spectral_helper(x=x, y=None, NFFT=len(x), Fs=Fs,
                                      detrend_func=detrend_none, window=window,
                                      noverlap=0, pad_to=pad_to,
                                      sides=sides,
                                      scale_by_freq=False,
                                      mode=mode)
    if mode != 'complex':
        spec = spec.real

    if spec.ndim == 2 and spec.shape[1] == 1:
        spec = spec[:, 0]

    return spec, freqs


# Split out these keyword docs so that they can be used elsewhere
docstring.interpd.update(Spectral=cbook.dedent("""
    Fs : scalar
        The sampling frequency (samples per time unit).  It is used
        to calculate the Fourier frequencies, freqs, in cycles per time
        unit. The default value is 2.

    window : callable or ndarray
        A function or a vector of length *NFFT*. To create window
        vectors see :func:`window_hanning`, :func:`window_none`,
        :func:`numpy.blackman`, :func:`numpy.hamming`,
        :func:`numpy.bartlett`, :func:`scipy.signal`,
        :func:`scipy.signal.get_window`, etc. The default is
        :func:`window_hanning`.  If a function is passed as the
        argument, it must take a data segment as an argument and
        return the windowed version of the segment.

    sides : [ 'default' | 'onesided' | 'twosided' ]
        Specifies which sides of the spectrum to return.  Default gives the
        default behavior, which returns one-sided for real data and both
        for complex data.  'onesided' forces the return of a one-sided
        spectrum, while 'twosided' forces two-sided.
"""))


docstring.interpd.update(Single_Spectrum=cbook.dedent("""
    pad_to : integer
        The number of points to which the data segment is padded when
        performing the FFT.  While not increasing the actual resolution of
        the spectrum (the minimum distance between resolvable peaks),
        this can give more points in the plot, allowing for more
        detail. This corresponds to the *n* parameter in the call to fft().
        The default is None, which sets *pad_to* equal to the length of the
        input signal (i.e. no padding).
"""))


docstring.interpd.update(PSD=cbook.dedent("""
    pad_to : integer
        The number of points to which the data segment is padded when
        performing the FFT.  This can be different from *NFFT*, which
        specifies the number of data points used.  While not increasing
        the actual resolution of the spectrum (the minimum distance between
        resolvable peaks), this can give more points in the plot,
        allowing for more detail. This corresponds to the *n* parameter
        in the call to fft(). The default is None, which sets *pad_to*
        equal to *NFFT*

    NFFT : integer
        The number of data points used in each block for the FFT.
        A power 2 is most efficient.  The default value is 256.
        This should *NOT* be used to get zero padding, or the scaling of the
        result will be incorrect. Use *pad_to* for this instead.

    detrend : {'default', 'constant', 'mean', 'linear', 'none'} or callable
        The function applied to each segment before fft-ing,
        designed to remove the mean or linear trend.  Unlike in
        MATLAB, where the *detrend* parameter is a vector, in
        matplotlib is it a function.  The :mod:`~matplotlib.pylab`
        module defines :func:`~matplotlib.pylab.detrend_none`,
        :func:`~matplotlib.pylab.detrend_mean`, and
        :func:`~matplotlib.pylab.detrend_linear`, but you can use
        a custom function as well.  You can also use a string to choose
        one of the functions.  'default', 'constant', and 'mean' call
        :func:`~matplotlib.pylab.detrend_mean`.  'linear' calls
        :func:`~matplotlib.pylab.detrend_linear`.  'none' calls
        :func:`~matplotlib.pylab.detrend_none`.

    scale_by_freq : boolean, optional
        Specifies whether the resulting density values should be scaled
        by the scaling frequency, which gives density in units of Hz^-1.
        This allows for integration over the returned frequency values.
        The default is True for MATLAB compatibility.
"""))


@docstring.dedent_interpd
def psd(x, NFFT=None, Fs=None, detrend=None, window=None,
        noverlap=None, pad_to=None, sides=None, scale_by_freq=None):
    r"""
    Compute the power spectral density.

    Call signature::

        psd(x, NFFT=256, Fs=2, detrend=mlab.detrend_none,
            window=mlab.window_hanning, noverlap=0, pad_to=None,
            sides='default', scale_by_freq=None)

    The power spectral density :math:`P_{xx}` by Welch's average
    periodogram method.  The vector *x* is divided into *NFFT* length
    segments.  Each segment is detrended by function *detrend* and
    windowed by function *window*.  *noverlap* gives the length of
    the overlap between segments.  The :math:`|\mathrm{fft}(i)|^2`
    of each segment :math:`i` are averaged to compute :math:`P_{xx}`.

    If len(*x*) < *NFFT*, it will be zero padded to *NFFT*.

    Parameters
    ----------
    x : 1-D array or sequence
        Array or sequence containing the data

    %(Spectral)s

    %(PSD)s

    noverlap : integer
        The number of points of overlap between segments.
        The default value is 0 (no overlap).

    Returns
    -------
    Pxx : 1-D array
        The values for the power spectrum `P_{xx}` (real valued)

    freqs : 1-D array
        The frequencies corresponding to the elements in *Pxx*

    References
    ----------
    Bendat & Piersol -- Random Data: Analysis and Measurement Procedures, John
    Wiley & Sons (1986)

    See Also
    --------
    :func:`specgram`
        :func:`specgram` differs in the default overlap; in not returning the
        mean of the segment periodograms; and in returning the times of the
        segments.

    :func:`magnitude_spectrum`
        :func:`magnitude_spectrum` returns the magnitude spectrum.

    :func:`csd`
        :func:`csd` returns the spectral density between two signals.
    """
    Pxx, freqs = csd(x=x, y=None, NFFT=NFFT, Fs=Fs, detrend=detrend,
                     window=window, noverlap=noverlap, pad_to=pad_to,
                     sides=sides, scale_by_freq=scale_by_freq)
    return Pxx.real, freqs


@docstring.dedent_interpd
def csd(x, y, NFFT=None, Fs=None, detrend=None, window=None,
        noverlap=None, pad_to=None, sides=None, scale_by_freq=None):
    """
    Compute the cross-spectral density.

    Call signature::

        csd(x, y, NFFT=256, Fs=2, detrend=mlab.detrend_none,
            window=mlab.window_hanning, noverlap=0, pad_to=None,
            sides='default', scale_by_freq=None)

    The cross spectral density :math:`P_{xy}` by Welch's average
    periodogram method.  The vectors *x* and *y* are divided into
    *NFFT* length segments.  Each segment is detrended by function
    *detrend* and windowed by function *window*.  *noverlap* gives
    the length of the overlap between segments.  The product of
    the direct FFTs of *x* and *y* are averaged over each segment
    to compute :math:`P_{xy}`, with a scaling to correct for power
    loss due to windowing.

    If len(*x*) < *NFFT* or len(*y*) < *NFFT*, they will be zero
    padded to *NFFT*.

    Parameters
    ----------
    x, y : 1-D arrays or sequences
        Arrays or sequences containing the data

    %(Spectral)s

    %(PSD)s

    noverlap : integer
        The number of points of overlap between segments.
        The default value is 0 (no overlap).

    Returns
    -------
    Pxy : 1-D array
        The values for the cross spectrum `P_{xy}` before scaling (real valued)

    freqs : 1-D array
        The frequencies corresponding to the elements in *Pxy*

    References
    ----------
    Bendat & Piersol -- Random Data: Analysis and Measurement Procedures, John
    Wiley & Sons (1986)

    See Also
    --------
    :func:`psd`
        :func:`psd` is the equivalent to setting y=x.
    """
    if NFFT is None:
        NFFT = 256
    Pxy, freqs, _ = _spectral_helper(x=x, y=y, NFFT=NFFT, Fs=Fs,
                                     detrend_func=detrend, window=window,
                                     noverlap=noverlap, pad_to=pad_to,
                                     sides=sides, scale_by_freq=scale_by_freq,
                                     mode='psd')

    if Pxy.ndim == 2:
        if Pxy.shape[1] > 1:
            Pxy = Pxy.mean(axis=1)
        else:
            Pxy = Pxy[:, 0]
    return Pxy, freqs


@docstring.dedent_interpd
def complex_spectrum(x, Fs=None, window=None, pad_to=None,
                     sides=None):
    """
    Compute the complex-valued frequency spectrum of *x*.  Data is padded to a
    length of *pad_to* and the windowing function *window* is applied to the
    signal.

    Parameters
    ----------
    x : 1-D array or sequence
        Array or sequence containing the data

    %(Spectral)s

    %(Single_Spectrum)s

    Returns
    -------
    spectrum : 1-D array
        The values for the complex spectrum (complex valued)

    freqs : 1-D array
        The frequencies corresponding to the elements in *spectrum*

    See Also
    --------
    :func:`magnitude_spectrum`
        :func:`magnitude_spectrum` returns the absolute value of this function.

    :func:`angle_spectrum`
        :func:`angle_spectrum` returns the angle of this function.

    :func:`phase_spectrum`
        :func:`phase_spectrum` returns the phase (unwrapped angle) of this
        function.

    :func:`specgram`
        :func:`specgram` can return the complex spectrum of segments within the
        signal.
    """
    return _single_spectrum_helper(x=x, Fs=Fs, window=window, pad_to=pad_to,
                                   sides=sides, mode='complex')


@docstring.dedent_interpd
def magnitude_spectrum(x, Fs=None, window=None, pad_to=None,
                       sides=None):
    """
    Compute the magnitude (absolute value) of the frequency spectrum of
    *x*.  Data is padded to a length of *pad_to* and the windowing function
    *window* is applied to the signal.

    Parameters
    ----------
    x : 1-D array or sequence
        Array or sequence containing the data

    %(Spectral)s

    %(Single_Spectrum)s

    Returns
    -------
    spectrum : 1-D array
        The values for the magnitude spectrum (real valued)

    freqs : 1-D array
        The frequencies corresponding to the elements in *spectrum*

    See Also
    --------
    :func:`psd`
        :func:`psd` returns the power spectral density.

    :func:`complex_spectrum`
        This function returns the absolute value of :func:`complex_spectrum`.

    :func:`angle_spectrum`
        :func:`angle_spectrum` returns the angles of the corresponding
        frequencies.

    :func:`phase_spectrum`
        :func:`phase_spectrum` returns the phase (unwrapped angle) of the
        corresponding frequencies.

    :func:`specgram`
        :func:`specgram` can return the magnitude spectrum of segments within
        the signal.
    """
    return _single_spectrum_helper(x=x, Fs=Fs, window=window, pad_to=pad_to,
                                   sides=sides, mode='magnitude')


@docstring.dedent_interpd
def angle_spectrum(x, Fs=None, window=None, pad_to=None,
                   sides=None):
    """
    Compute the angle of the frequency spectrum (wrapped phase spectrum) of
    *x*.  Data is padded to a length of *pad_to* and the windowing function
    *window* is applied to the signal.

    Parameters
    ----------
    x : 1-D array or sequence
        Array or sequence containing the data

    %(Spectral)s

    %(Single_Spectrum)s

    Returns
    -------
    spectrum : 1-D array
        The values for the angle spectrum in radians (real valued)

    freqs : 1-D array
        The frequencies corresponding to the elements in *spectrum*

    See Also
    --------
    :func:`complex_spectrum`
        This function returns the angle value of :func:`complex_spectrum`.

    :func:`magnitude_spectrum`
        :func:`angle_spectrum` returns the magnitudes of the corresponding
        frequencies.

    :func:`phase_spectrum`
        :func:`phase_spectrum` returns the unwrapped version of this function.

    :func:`specgram`
        :func:`specgram` can return the angle spectrum of segments within the
        signal.
    """
    return _single_spectrum_helper(x=x, Fs=Fs, window=window, pad_to=pad_to,
                                   sides=sides, mode='angle')


@docstring.dedent_interpd
def phase_spectrum(x, Fs=None, window=None, pad_to=None,
                   sides=None):
    """
    Compute the phase of the frequency spectrum (unwrapped angle spectrum) of
    *x*.  Data is padded to a length of *pad_to* and the windowing function
    *window* is applied to the signal.

    Parameters
    ----------
    x : 1-D array or sequence
        Array or sequence containing the data

    %(Spectral)s

    %(Single_Spectrum)s

    Returns
    -------
    spectrum : 1-D array
        The values for the phase spectrum in radians (real valued)

    freqs : 1-D array
        The frequencies corresponding to the elements in *spectrum*

    See Also
    --------
    :func:`complex_spectrum`
        This function returns the angle value of :func:`complex_spectrum`.

    :func:`magnitude_spectrum`
        :func:`magnitude_spectrum` returns the magnitudes of the corresponding
        frequencies.

    :func:`angle_spectrum`
        :func:`angle_spectrum` returns the wrapped version of this function.

    :func:`specgram`
        :func:`specgram` can return the phase spectrum of segments within the
        signal.
    """
    return _single_spectrum_helper(x=x, Fs=Fs, window=window, pad_to=pad_to,
                                   sides=sides, mode='phase')


@docstring.dedent_interpd
def specgram(x, NFFT=None, Fs=None, detrend=None, window=None,
             noverlap=None, pad_to=None, sides=None, scale_by_freq=None,
             mode=None):
    """
    Compute a spectrogram.

    Compute and plot a spectrogram of data in x.  Data are split into
    NFFT length segments and the spectrum of each section is
    computed.  The windowing function window is applied to each
    segment, and the amount of overlap of each segment is
    specified with noverlap.

    Parameters
    ----------
    x : array_like
        1-D array or sequence.

    %(Spectral)s

    %(PSD)s

    noverlap : int, optional
        The number of points of overlap between blocks.  The default
        value is 128.
    mode : str, optional
        What sort of spectrum to use, default is 'psd'.
            'psd'
                Returns the power spectral density.

            'complex'
                Returns the complex-valued frequency spectrum.

            'magnitude'
                Returns the magnitude spectrum.

            'angle'
                Returns the phase spectrum without unwrapping.

            'phase'
                Returns the phase spectrum with unwrapping.

    Returns
    -------
    spectrum : array_like
        2-D array, columns are the periodograms of successive segments.

    freqs : array_like
        1-D array, frequencies corresponding to the rows in *spectrum*.

    t : array_like
        1-D array, the times corresponding to midpoints of segments
        (i.e the columns in *spectrum*).

    See Also
    --------
    psd : differs in the overlap and in the return values.
    complex_spectrum : similar, but with complex valued frequencies.
    magnitude_spectrum : similar single segment when mode is 'magnitude'.
    angle_spectrum : similar to single segment when mode is 'angle'.
    phase_spectrum : similar to single segment when mode is 'phase'.

    Notes
    -----
    detrend and scale_by_freq only apply when *mode* is set to 'psd'.

    """
    if noverlap is None:
        noverlap = 128  # default in _spectral_helper() is noverlap = 0
    if NFFT is None:
        NFFT = 256  # same default as in _spectral_helper()
    if len(x) <= NFFT:
        warnings.warn("Only one segment is calculated since parameter NFFT " +
                      "(=%d) >= signal length (=%d)." % (NFFT, len(x)))

    spec, freqs, t = _spectral_helper(x=x, y=None, NFFT=NFFT, Fs=Fs,
                                      detrend_func=detrend, window=window,
                                      noverlap=noverlap, pad_to=pad_to,
                                      sides=sides,
                                      scale_by_freq=scale_by_freq,
                                      mode=mode)

    if mode != 'complex':
        spec = spec.real  # Needed since helper implements generically

    return spec, freqs, t


_coh_error = """Coherence is calculated by averaging over *NFFT*
length segments.  Your signal is too short for your choice of *NFFT*.
"""


@docstring.dedent_interpd
def cohere(x, y, NFFT=256, Fs=2, detrend=detrend_none, window=window_hanning,
           noverlap=0, pad_to=None, sides='default', scale_by_freq=None):
    """
    The coherence between *x* and *y*.  Coherence is the normalized
    cross spectral density:

    .. math::

        C_{xy} = \\frac{|P_{xy}|^2}{P_{xx}P_{yy}}

    Parameters
    ----------
    x, y
        Array or sequence containing the data

    %(Spectral)s

    %(PSD)s

    noverlap : integer
        The number of points of overlap between blocks.  The default value
        is 0 (no overlap).

    Returns
    -------
    The return value is the tuple (*Cxy*, *f*), where *f* are the
    frequencies of the coherence vector. For cohere, scaling the
    individual densities by the sampling frequency has no effect,
    since the factors cancel out.

    See Also
    --------
    :func:`psd`, :func:`csd` :
        For information about the methods used to compute :math:`P_{xy}`,
        :math:`P_{xx}` and :math:`P_{yy}`.
    """

    if len(x) < 2 * NFFT:
        raise ValueError(_coh_error)
    Pxx, f = psd(x, NFFT, Fs, detrend, window, noverlap, pad_to, sides,
                 scale_by_freq)
    Pyy, f = psd(y, NFFT, Fs, detrend, window, noverlap, pad_to, sides,
                 scale_by_freq)
    Pxy, f = csd(x, y, NFFT, Fs, detrend, window, noverlap, pad_to, sides,
                 scale_by_freq)
    Cxy = np.abs(Pxy) ** 2 / (Pxx * Pyy)
    return Cxy, f


@cbook.deprecated('2.2')
def donothing_callback(*args):
    pass


@cbook.deprecated('2.2', 'scipy.signal.coherence')
def cohere_pairs(X, ij, NFFT=256, Fs=2, detrend=detrend_none,
                 window=window_hanning, noverlap=0,
                 preferSpeedOverMemory=True,
                 progressCallback=donothing_callback,
                 returnPxx=False):

    """
    Compute the coherence and phase for all pairs *ij*, in *X*.

    *X* is a *numSamples* * *numCols* array

    *ij* is a list of tuples.  Each tuple is a pair of indexes into
    the columns of X for which you want to compute coherence.  For
    example, if *X* has 64 columns, and you want to compute all
    nonredundant pairs, define *ij* as::

      ij = []
      for i in range(64):
          for j in range(i+1,64):
              ij.append( (i,j) )

    *preferSpeedOverMemory* is an optional bool. Defaults to true. If
    False, limits the caching by only making one, rather than two,
    complex cache arrays. This is useful if memory becomes critical.
    Even when *preferSpeedOverMemory* is False, :func:`cohere_pairs`
    will still give significant performance gains over calling
    :func:`cohere` for each pair, and will use subtantially less
    memory than if *preferSpeedOverMemory* is True.  In my tests with
    a 43000,64 array over all nonredundant pairs,
    *preferSpeedOverMemory* = True delivered a 33% performance boost
    on a 1.7GHZ Athlon with 512MB RAM compared with
    *preferSpeedOverMemory* = False.  But both solutions were more
    than 10x faster than naively crunching all possible pairs through
    :func:`cohere`.

    Returns
    -------
    Cxy : dictionary of (*i*, *j*) tuples -> coherence vector for
        that pair.  i.e., ``Cxy[(i,j) = cohere(X[:,i], X[:,j])``.
        Number of dictionary keys is ``len(ij)``.

    Phase : dictionary of phases of the cross spectral density at
        each frequency for each pair.  Keys are (*i*, *j*).

    freqs : vector of frequencies, equal in length to either the
         coherence or phase vectors for any (*i*, *j*) key.

    e.g., to make a coherence Bode plot::

          subplot(211)
          plot( freqs, Cxy[(12,19)])
          subplot(212)
          plot( freqs, Phase[(12,19)])

    For a large number of pairs, :func:`cohere_pairs` can be much more
    efficient than just calling :func:`cohere` for each pair, because
    it caches most of the intensive computations.  If :math:`N` is the
    number of pairs, this function is :math:`O(N)` for most of the
    heavy lifting, whereas calling cohere for each pair is
    :math:`O(N^2)`.  However, because of the caching, it is also more
    memory intensive, making 2 additional complex arrays with
    approximately the same number of elements as *X*.

    See :file:`test/cohere_pairs_test.py` in the src tree for an
    example script that shows that this :func:`cohere_pairs` and
    :func:`cohere` give the same results for a given pair.

    See Also
    --------
    :func:`psd`
        For information about the methods used to compute :math:`P_{xy}`,
        :math:`P_{xx}` and :math:`P_{yy}`.
    """
    numRows, numCols = X.shape

    # zero pad if X is too short
    if numRows < NFFT:
        tmp = X
        X = np.zeros((NFFT, numCols), X.dtype)
        X[:numRows, :] = tmp
        del tmp

    numRows, numCols = X.shape
    # get all the columns of X that we are interested in by checking
    # the ij tuples
    allColumns = set()
    for i, j in ij:
        allColumns.add(i)
        allColumns.add(j)
    Ncols = len(allColumns)

    # for real X, ignore the negative frequencies
    if np.iscomplexobj(X):
        numFreqs = NFFT
    else:
        numFreqs = NFFT//2+1

    # cache the FFT of every windowed, detrended NFFT length segment
    # of every channel.  If preferSpeedOverMemory, cache the conjugate
    # as well
    if cbook.iterable(window):
        if len(window) != NFFT:
            raise ValueError("The length of the window must be equal to NFFT")
        windowVals = window
    else:
        windowVals = window(np.ones(NFFT, X.dtype))
    ind = list(xrange(0, numRows-NFFT+1, NFFT-noverlap))
    numSlices = len(ind)
    FFTSlices = {}
    FFTConjSlices = {}
    Pxx = {}
    slices = range(numSlices)
    normVal = np.linalg.norm(windowVals)**2
    for iCol in allColumns:
        progressCallback(i/Ncols, 'Cacheing FFTs')
        Slices = np.zeros((numSlices, numFreqs), dtype=np.complex_)
        for iSlice in slices:
            thisSlice = X[ind[iSlice]:ind[iSlice]+NFFT, iCol]
            thisSlice = windowVals*detrend(thisSlice)
            Slices[iSlice, :] = np.fft.fft(thisSlice)[:numFreqs]

        FFTSlices[iCol] = Slices
        if preferSpeedOverMemory:
            FFTConjSlices[iCol] = np.conj(Slices)
        Pxx[iCol] = np.divide(np.mean(abs(Slices)**2, axis=0), normVal)
    del Slices, ind, windowVals

    # compute the coherences and phases for all pairs using the
    # cached FFTs
    Cxy = {}
    Phase = {}
    count = 0
    N = len(ij)
    for i, j in ij:
        count += 1
        if count % 10 == 0:
            progressCallback(count/N, 'Computing coherences')

        if preferSpeedOverMemory:
            Pxy = FFTSlices[i] * FFTConjSlices[j]
        else:
            Pxy = FFTSlices[i] * np.conj(FFTSlices[j])
        if numSlices > 1:
            Pxy = np.mean(Pxy, axis=0)
#       Pxy = np.divide(Pxy, normVal)
        Pxy /= normVal
#       Cxy[(i,j)] = np.divide(np.absolute(Pxy)**2, Pxx[i]*Pxx[j])
        Cxy[i, j] = abs(Pxy)**2 / (Pxx[i]*Pxx[j])
        Phase[i, j] = np.arctan2(Pxy.imag, Pxy.real)

    freqs = Fs/NFFT*np.arange(numFreqs)
    if returnPxx:
        return Cxy, Phase, freqs, Pxx
    else:
        return Cxy, Phase, freqs


@cbook.deprecated('2.2', 'scipy.stats.entropy')
def entropy(y, bins):
    r"""
    Return the entropy of the data in *y* in units of nat.

    .. math::

      -\sum p_i \ln(p_i)

    where :math:`p_i` is the probability of observing *y* in the
    :math:`i^{th}` bin of *bins*.  *bins* can be a number of bins or a
    range of bins; see :func:`numpy.histogram`.

    Compare *S* with analytic calculation for a Gaussian::

      x = mu + sigma * randn(200000)
      Sanalytic = 0.5 * ( 1.0 + log(2*pi*sigma**2.0) )
    """
    n, bins = np.histogram(y, bins)
    n = n.astype(float)

    n = np.take(n, np.nonzero(n)[0])         # get the positive

    p = np.divide(n, len(y))

    delta = bins[1] - bins[0]
    S = -1.0 * np.sum(p * np.log(p)) + np.log(delta)
    return S


@cbook.deprecated('2.2', 'scipy.stats.norm.pdf')
def normpdf(x, *args):
    "Return the normal pdf evaluated at *x*; args provides *mu*, *sigma*"
    mu, sigma = args
    return 1./(np.sqrt(2*np.pi)*sigma)*np.exp(-0.5 * (1./sigma*(x - mu))**2)


@cbook.deprecated('2.2')
def find(condition):
    "Return the indices where ravel(condition) is true"
    res, = np.nonzero(np.ravel(condition))
    return res


@cbook.deprecated('2.2')
def longest_contiguous_ones(x):
    """
    Return the indices of the longest stretch of contiguous ones in *x*,
    assuming *x* is a vector of zeros and ones.  If there are two
    equally long stretches, pick the first.
    """
    x = np.ravel(x)
    if len(x) == 0:
        return np.array([])

    ind = (x == 0).nonzero()[0]
    if len(ind) == 0:
        return np.arange(len(x))
    if len(ind) == len(x):
        return np.array([])

    y = np.zeros((len(x)+2,), x.dtype)
    y[1:-1] = x
    dif = np.diff(y)
    up = (dif == 1).nonzero()[0]
    dn = (dif == -1).nonzero()[0]
    i = (dn-up == max(dn - up)).nonzero()[0][0]
    ind = np.arange(up[i], dn[i])

    return ind


@cbook.deprecated('2.2')
def longest_ones(x):
    '''alias for longest_contiguous_ones'''
    return longest_contiguous_ones(x)


@cbook.deprecated('2.2')
class PCA(object):
    def __init__(self, a, standardize=True):
        """
        compute the SVD of a and store data for PCA.  Use project to
        project the data onto a reduced set of dimensions

        Parameters
        ----------
        a : np.ndarray
            A numobservations x numdims array
        standardize : bool
            True if input data are to be standardized. If False, only centering
            will be carried out.

        Attributes
        ----------
        a
            A centered unit sigma version of input ``a``.

        numrows, numcols
            The dimensions of ``a``.

        mu
            A numdims array of means of ``a``. This is the vector that points
            to the origin of PCA space.

        sigma
            A numdims array of standard deviation of ``a``.

        fracs
            The proportion of variance of each of the principal components.

        s
            The actual eigenvalues of the decomposition.

        Wt
            The weight vector for projecting a numdims point or array into
            PCA space.

        Y
            A projected into PCA space.

        Notes
        -----
        The factor loadings are in the ``Wt`` factor, i.e., the factor loadings
        for the first principal component are given by ``Wt[0]``. This row is
        also the first eigenvector.

        """
        n, m = a.shape
        if n < m:
            raise RuntimeError('we assume data in a is organized with '
                               'numrows>numcols')

        self.numrows, self.numcols = n, m
        self.mu = a.mean(axis=0)
        self.sigma = a.std(axis=0)
        self.standardize = standardize

        a = self.center(a)

        self.a = a

        U, s, Vh = np.linalg.svd(a, full_matrices=False)

        # Note: .H indicates the conjugate transposed / Hermitian.

        # The SVD is commonly written as a = U s V.H.
        # If U is a unitary matrix, it means that it satisfies U.H = inv(U).

        # The rows of Vh are the eigenvectors of a.H a.
        # The columns of U are the eigenvectors of a a.H.
        # For row i in Vh and column i in U, the corresponding eigenvalue is
        # s[i]**2.

        self.Wt = Vh

        # save the transposed coordinates
        Y = np.dot(Vh, a.T).T
        self.Y = Y

        # save the eigenvalues
        self.s = s**2

        # and now the contribution of the individual components
        vars = self.s / len(s)
        self.fracs = vars/vars.sum()

    def project(self, x, minfrac=0.):
        '''
        project x onto the principle axes, dropping any axes where fraction
        of variance<minfrac
        '''
        x = np.asarray(x)
        if x.shape[-1] != self.numcols:
            raise ValueError('Expected an array with dims[-1]==%d' %
                             self.numcols)
        Y = np.dot(self.Wt, self.center(x).T).T
        mask = self.fracs >= minfrac
        if x.ndim == 2:
            Yreduced = Y[:, mask]
        else:
            Yreduced = Y[mask]
        return Yreduced

    def center(self, x):
        '''
        center and optionally standardize the data using the mean and sigma
        from training set a
        '''
        if self.standardize:
            return (x - self.mu)/self.sigma
        else:
            return (x - self.mu)

    @staticmethod
    def _get_colinear():
        c0 = np.array([
            0.19294738,  0.6202667,   0.45962655,  0.07608613,  0.135818,
            0.83580842,  0.07218851,  0.48318321,  0.84472463,  0.18348462,
            0.81585306,  0.96923926,  0.12835919,  0.35075355,  0.15807861,
            0.837437,    0.10824303,  0.1723387,   0.43926494,  0.83705486])

        c1 = np.array([
            -1.17705601, -0.513883,   -0.26614584,  0.88067144,  1.00474954,
            -1.1616545,   0.0266109,   0.38227157,  1.80489433,  0.21472396,
            -1.41920399, -2.08158544, -0.10559009,  1.68999268,  0.34847107,
            -0.4685737,   1.23980423, -0.14638744, -0.35907697,  0.22442616])

        c2 = c0 + 2*c1
        c3 = -3*c0 + 4*c1
        a = np.array([c3, c0, c1, c2]).T
        return a


@cbook.deprecated('2.2', 'numpy.percentile')
def prctile(x, p=(0.0, 25.0, 50.0, 75.0, 100.0)):
    """
    Return the percentiles of *x*.  *p* can either be a sequence of
    percentile values or a scalar.  If *p* is a sequence, the ith
    element of the return sequence is the *p*(i)-th percentile of *x*.
    If *p* is a scalar, the largest value of *x* less than or equal to
    the *p* percentage point in the sequence is returned.
    """

    # This implementation derived from scipy.stats.scoreatpercentile
    def _interpolate(a, b, fraction):
        """Returns the point at the given fraction between a and b, where
        'fraction' must be between 0 and 1.
        """
        return a + (b - a) * fraction

    per = np.array(p)
    values = np.sort(x, axis=None)

    idxs = per / 100 * (values.shape[0] - 1)
    ai = idxs.astype(int)
    bi = ai + 1
    frac = idxs % 1

    # handle cases where attempting to interpolate past last index
    cond = bi >= len(values)
    if per.ndim:
        ai[cond] -= 1
        bi[cond] -= 1
        frac[cond] += 1
    else:
        if cond:
            ai -= 1
            bi -= 1
            frac += 1

    return _interpolate(values[ai], values[bi], frac)


@cbook.deprecated('2.2')
def prctile_rank(x, p):
    """
    Return the rank for each element in *x*, return the rank
    0..len(*p*).  e.g., if *p* = (25, 50, 75), the return value will be a
    len(*x*) array with values in [0,1,2,3] where 0 indicates the
    value is less than the 25th percentile, 1 indicates the value is
    >= the 25th and < 50th percentile, ... and 3 indicates the value
    is above the 75th percentile cutoff.

    *p* is either an array of percentiles in [0..100] or a scalar which
    indicates how many quantiles of data you want ranked.
    """

    if not cbook.iterable(p):
        p = np.arange(100.0/p, 100.0, 100.0/p)
    else:
        p = np.asarray(p)

    if p.max() <= 1 or p.min() < 0 or p.max() > 100:
        raise ValueError('percentiles should be in range 0..100, not 0..1')

    ptiles = prctile(x, p)
    return np.searchsorted(ptiles, x)


@cbook.deprecated('2.2')
def center_matrix(M, dim=0):
    """
    Return the matrix *M* with each row having zero mean and unit std.

    If *dim* = 1 operate on columns instead of rows.  (*dim* is
    opposite to the numpy axis kwarg.)
    """
    M = np.asarray(M, float)
    if dim:
        M = (M - M.mean(axis=0)) / M.std(axis=0)
    else:
        M = (M - M.mean(axis=1)[:, np.newaxis])
        M = M / M.std(axis=1)[:, np.newaxis]
    return M


@cbook.deprecated('2.2', 'scipy.integrate.ode')
def rk4(derivs, y0, t):
    """
    Integrate 1D or ND system of ODEs using 4-th order Runge-Kutta.
    This is a toy implementation which may be useful if you find
    yourself stranded on a system w/o scipy.  Otherwise use
    :func:`scipy.integrate`.

    Parameters
    ----------
    y0
        initial state vector

    t
        sample times

    derivs
        returns the derivative of the system and has the
        signature ``dy = derivs(yi, ti)``

    Examples
    --------

    A 2D system::

        def derivs6(x,t):
            d1 =  x[0] + 2*x[1]
            d2 =  -3*x[0] + 4*x[1]
            return (d1, d2)
        dt = 0.0005
        t = arange(0.0, 2.0, dt)
        y0 = (1,2)
        yout = rk4(derivs6, y0, t)

    A 1D system::

        alpha = 2
        def derivs(x,t):
            return -alpha*x + exp(-t)

        y0 = 1
        yout = rk4(derivs, y0, t)

    If you have access to scipy, you should probably be using the
    scipy.integrate tools rather than this function.
    """

    try:
        Ny = len(y0)
    except TypeError:
        yout = np.zeros((len(t),), float)
    else:
        yout = np.zeros((len(t), Ny), float)

    yout[0] = y0
    i = 0

    for i in np.arange(len(t)-1):

        thist = t[i]
        dt = t[i+1] - thist
        dt2 = dt/2.0
        y0 = yout[i]

        k1 = np.asarray(derivs(y0, thist))
        k2 = np.asarray(derivs(y0 + dt2*k1, thist+dt2))
        k3 = np.asarray(derivs(y0 + dt2*k2, thist+dt2))
        k4 = np.asarray(derivs(y0 + dt*k3, thist+dt))
        yout[i+1] = y0 + dt/6.0*(k1 + 2*k2 + 2*k3 + k4)
    return yout


@cbook.deprecated('2.2')
def bivariate_normal(X, Y, sigmax=1.0, sigmay=1.0,
                     mux=0.0, muy=0.0, sigmaxy=0.0):
    """
    Bivariate Gaussian distribution for equal shape *X*, *Y*.

    See `bivariate normal
    <http://mathworld.wolfram.com/BivariateNormalDistribution.html>`_
    at mathworld.
    """
    Xmu = X-mux
    Ymu = Y-muy

    rho = sigmaxy/(sigmax*sigmay)
    z = Xmu**2/sigmax**2 + Ymu**2/sigmay**2 - 2*rho*Xmu*Ymu/(sigmax*sigmay)
    denom = 2*np.pi*sigmax*sigmay*np.sqrt(1-rho**2)
    return np.exp(-z/(2*(1-rho**2))) / denom


@cbook.deprecated('2.2')
def get_xyz_where(Z, Cond):
    """
    *Z* and *Cond* are *M* x *N* matrices.  *Z* are data and *Cond* is
    a boolean matrix where some condition is satisfied.  Return value
    is (*x*, *y*, *z*) where *x* and *y* are the indices into *Z* and
    *z* are the values of *Z* at those indices.  *x*, *y*, and *z* are
    1D arrays.
    """
    X, Y = np.indices(Z.shape)
    return X[Cond], Y[Cond], Z[Cond]


@cbook.deprecated('2.2')
def get_sparse_matrix(M, N, frac=0.1):
    """
    Return a *M* x *N* sparse matrix with *frac* elements randomly
    filled.
    """
    data = np.zeros((M, N))*0.
    for i in range(int(M*N*frac)):
        x = np.random.randint(0, M-1)
        y = np.random.randint(0, N-1)
        data[x, y] = np.random.rand()
    return data


@cbook.deprecated('2.2', 'numpy.hypot')
def dist(x, y):
    """
    Return the distance between two points.
    """
    d = x-y
    return np.sqrt(np.dot(d, d))


@cbook.deprecated('2.2')
def dist_point_to_segment(p, s0, s1):
    """
    Get the distance of a point to a segment.

      *p*, *s0*, *s1* are *xy* sequences

    This algorithm from
    http://geomalgorithms.com/a02-_lines.html
    """
    p = np.asarray(p, float)
    s0 = np.asarray(s0, float)
    s1 = np.asarray(s1, float)
    v = s1 - s0
    w = p - s0

    c1 = np.dot(w, v)
    if c1 <= 0:
        return dist(p, s0)

    c2 = np.dot(v, v)
    if c2 <= c1:
        return dist(p, s1)

    b = c1 / c2
    pb = s0 + b * v
    return dist(p, pb)


@cbook.deprecated('2.2')
def segments_intersect(s1, s2):
    """
    Return *True* if *s1* and *s2* intersect.
    *s1* and *s2* are defined as::

      s1: (x1, y1), (x2, y2)
      s2: (x3, y3), (x4, y4)
    """
    (x1, y1), (x2, y2) = s1
    (x3, y3), (x4, y4) = s2

    den = ((y4-y3) * (x2-x1)) - ((x4-x3)*(y2-y1))

    n1 = ((x4-x3) * (y1-y3)) - ((y4-y3)*(x1-x3))
    n2 = ((x2-x1) * (y1-y3)) - ((y2-y1)*(x1-x3))

    if den == 0:
        # lines parallel
        return False

    u1 = n1/den
    u2 = n2/den

    return 0.0 <= u1 <= 1.0 and 0.0 <= u2 <= 1.0


@cbook.deprecated('2.2')
def fftsurr(x, detrend=detrend_none, window=window_none):
    """
    Compute an FFT phase randomized surrogate of *x*.
    """
    if cbook.iterable(window):
        x = window*detrend(x)
    else:
        x = window(detrend(x))
    z = np.fft.fft(x)
    a = 2.*np.pi*1j
    phase = a * np.random.rand(len(x))
    z = z*np.exp(phase)
    return np.fft.ifft(z).real


@cbook.deprecated('2.2')
def movavg(x, n):
    """
    Compute the len(*n*) moving average of *x*.
    """
    w = np.empty((n,), dtype=float)
    w[:] = 1.0/n
    return np.convolve(x, w, mode='valid')


# the following code was written and submitted by Fernando Perez
# from the ipython numutils package under a BSD license
# begin fperez functions

"""
A set of convenient utilities for numerical work.

Most of this module requires numpy or is meant to be used with it.

Copyright (c) 2001-2004, Fernando Perez. <Fernando.Perez@colorado.edu>
All rights reserved.

This license was generated from the BSD license template as found in:
http://www.opensource.org/licenses/bsd-license.php

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    * Neither the name of the IPython project nor the names of its
      contributors may be used to endorse or promote products derived from
      this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""


# *****************************************************************************
# Globals
# ****************************************************************************
# function definitions
exp_safe_MIN = math.log(2.2250738585072014e-308)
exp_safe_MAX = 1.7976931348623157e+308


@cbook.deprecated("2.2", 'numpy.exp')
def exp_safe(x):
    """
    Compute exponentials which safely underflow to zero.

    Slow, but convenient to use. Note that numpy provides proper
    floating point exception handling with access to the underlying
    hardware.
    """

    if type(x) is np.ndarray:
        return np.exp(np.clip(x, exp_safe_MIN, exp_safe_MAX))
    else:
        return math.exp(x)


@cbook.deprecated("2.2", alternative='numpy.array(list(map(...)))')
def amap(fn, *args):
    """
    amap(function, sequence[, sequence, ...]) -> array.

    Works like :func:`map`, but it returns an array.  This is just a
    convenient shorthand for ``numpy.array(map(...))``.
    """
    return np.array(list(map(fn, *args)))


@cbook.deprecated("2.2")
def rms_flat(a):
    """
    Return the root mean square of all the elements of *a*, flattened out.
    """
    return np.sqrt(np.mean(np.abs(a) ** 2))


@cbook.deprecated("2.2", alternative='numpy.linalg.norm(a, ord=1)')
def l1norm(a):
    """
    Return the *l1* norm of *a*, flattened out.

    Implemented as a separate function (not a call to :func:`norm` for speed).
    """
    return np.sum(np.abs(a))


@cbook.deprecated("2.2", alternative='numpy.linalg.norm(a, ord=2)')
def l2norm(a):
    """
    Return the *l2* norm of *a*, flattened out.

    Implemented as a separate function (not a call to :func:`norm` for speed).
    """
    return np.sqrt(np.sum(np.abs(a) ** 2))


@cbook.deprecated("2.2", alternative='numpy.linalg.norm(a.flat, ord=p)')
def norm_flat(a, p=2):
    """
    norm(a,p=2) -> l-p norm of a.flat

    Return the l-p norm of *a*, considered as a flat array.  This is NOT a true
    matrix norm, since arrays of arbitrary rank are always flattened.

    *p* can be a number or the string 'Infinity' to get the L-infinity norm.
    """
    # This function was being masked by a more general norm later in
    # the file.  We may want to simply delete it.
    if p == 'Infinity':
        return np.max(np.abs(a))
    else:
        return np.sum(np.abs(a) ** p) ** (1 / p)


@cbook.deprecated("2.2", 'numpy.arange')
def frange(xini, xfin=None, delta=None, **kw):
    """
    frange([start,] stop[, step, keywords]) -> array of floats

    Return a numpy ndarray containing a progression of floats. Similar to
    :func:`numpy.arange`, but defaults to a closed interval.

    ``frange(x0, x1)`` returns ``[x0, x0+1, x0+2, ..., x1]``; *start*
    defaults to 0, and the endpoint *is included*. This behavior is
    different from that of :func:`range` and
    :func:`numpy.arange`. This is deliberate, since :func:`frange`
    will probably be more useful for generating lists of points for
    function evaluation, and endpoints are often desired in this
    use. The usual behavior of :func:`range` can be obtained by
    setting the keyword *closed* = 0, in this case, :func:`frange`
    basically becomes :func:numpy.arange`.

    When *step* is given, it specifies the increment (or
    decrement). All arguments can be floating point numbers.

    ``frange(x0,x1,d)`` returns ``[x0,x0+d,x0+2d,...,xfin]`` where
    *xfin* <= *x1*.

    :func:`frange` can also be called with the keyword *npts*. This
    sets the number of points the list should contain (and overrides
    the value *step* might have been given). :func:`numpy.arange`
    doesn't offer this option.

    Examples::

      >>> frange(3)
      array([ 0.,  1.,  2.,  3.])
      >>> frange(3,closed=0)
      array([ 0.,  1.,  2.])
      >>> frange(1,6,2)
      array([1, 3, 5])   or 1,3,5,7, depending on floating point vagueries
      >>> frange(1,6.5,npts=5)
      array([ 1.   ,  2.375,  3.75 ,  5.125,  6.5  ])
    """

    # defaults
    kw.setdefault('closed', 1)
    endpoint = kw['closed'] != 0

    # funny logic to allow the *first* argument to be optional (like range())
    # This was modified with a simpler version from a similar frange() found
    # at http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/66472
    if xfin is None:
        xfin = xini + 0.0
        xini = 0.0

    if delta is None:
        delta = 1.0

    # compute # of points, spacing and return final list
    try:
        npts = kw['npts']
        delta = (xfin-xini) / (npts-endpoint)
    except KeyError:
        npts = int(np.round((xfin-xini)/delta)) + endpoint
        # round finds the nearest, so the endpoint can be up to
        # delta/2 larger than xfin.

    return np.arange(npts)*delta+xini
# end frange()


@cbook.deprecated("2.2", 'numpy.identity')
def identity(n, rank=2, dtype='l', typecode=None):
    """
    Returns the identity matrix of shape (*n*, *n*, ..., *n*) (rank *r*).

    For ranks higher than 2, this object is simply a multi-index Kronecker
    delta::

                            /  1  if i0=i1=...=iR,
        id[i0,i1,...,iR] = -|
                            \\  0  otherwise.

    Optionally a *dtype* (or typecode) may be given (it defaults to 'l').

    Since rank defaults to 2, this function behaves in the default case (when
    only *n* is given) like ``numpy.identity(n)`` -- but surprisingly, it is
    much faster.
    """
    if typecode is not None:
        dtype = typecode
    iden = np.zeros((n,)*rank, dtype)
    for i in range(n):
        idx = (i,)*rank
        iden[idx] = 1
    return iden


@cbook.deprecated("2.2")
def base_repr(number, base=2, padding=0):
    """
    Return the representation of a *number* in any given *base*.
    """
    chars = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    if number < base:
        return (padding - 1) * chars[0] + chars[int(number)]
    max_exponent = int(math.log(number)/math.log(base))
    max_power = long(base) ** max_exponent
    lead_digit = int(number/max_power)
    return (chars[lead_digit] +
            base_repr(number - max_power * lead_digit, base,
                      max(padding - 1, max_exponent)))


@cbook.deprecated("2.2")
def binary_repr(number, max_length=1025):
    """
    Return the binary representation of the input *number* as a
    string.

    This is more efficient than using :func:`base_repr` with base 2.

    Increase the value of max_length for very large numbers. Note that
    on 32-bit machines, 2**1023 is the largest integer power of 2
    which can be converted to a Python float.
    """

#   assert number < 2L << max_length
    shifts = map(operator.rshift, max_length * [number],
                 range(max_length - 1, -1, -1))
    digits = list(map(operator.mod, shifts, max_length * [2]))
    if not digits.count(1):
        return 0
    digits = digits[digits.index(1):]
    return ''.join(map(repr, digits)).replace('L', '')


@cbook.deprecated("2.2", 'numpy.log2')
def log2(x, ln2=math.log(2.0)):
    """
    Return the log(*x*) in base 2.

    This is a _slow_ function but which is guaranteed to return the correct
    integer value if the input is an integer exact power of 2.
    """
    try:
        bin_n = binary_repr(x)[1:]
    except (AssertionError, TypeError):
        return math.log(x)/ln2
    else:
        if '1' in bin_n:
            return math.log(x)/ln2
        else:
            return len(bin_n)


@cbook.deprecated("2.2")
def ispower2(n):
    """
    Returns the log base 2 of *n* if *n* is a power of 2, zero otherwise.

    Note the potential ambiguity if *n* == 1: 2**0 == 1, interpret accordingly.
    """

    bin_n = binary_repr(n)[1:]
    if '1' in bin_n:
        return 0
    else:
        return len(bin_n)


@cbook.deprecated("2.2")
def isvector(X):
    """
    Like the MATLAB function with the same name, returns *True*
    if the supplied numpy array or matrix *X* looks like a vector,
    meaning it has a one non-singleton axis (i.e., it can have
    multiple axes, but all must have length 1, except for one of
    them).

    If you just want to see if the array has 1 axis, use X.ndim == 1.
    """
    return np.prod(X.shape) == np.max(X.shape)

# end fperez numutils code


# helpers for loading, saving, manipulating and viewing numpy record arrays
@cbook.deprecated("2.2", 'numpy.isnan')
def safe_isnan(x):
    ':func:`numpy.isnan` for arbitrary types'
    if isinstance(x, six.string_types):
        return False
    try:
        b = np.isnan(x)
    except NotImplementedError:
        return False
    except TypeError:
        return False
    else:
        return b


@cbook.deprecated("2.2", 'numpy.isinf')
def safe_isinf(x):
    ':func:`numpy.isinf` for arbitrary types'
    if isinstance(x, six.string_types):
        return False
    try:
        b = np.isinf(x)
    except NotImplementedError:
        return False
    except TypeError:
        return False
    else:
        return b


@cbook.deprecated("2.2")
def rec_append_fields(rec, names, arrs, dtypes=None):
    """
    Return a new record array with field names populated with data
    from arrays in *arrs*.  If appending a single field, then *names*,
    *arrs* and *dtypes* do not have to be lists. They can just be the
    values themselves.
    """
    if (not isinstance(names, six.string_types) and cbook.iterable(names)
            and len(names) and isinstance(names[0], six.string_types)):
        if len(names) != len(arrs):
            raise ValueError("number of arrays do not match number of names")
    else:  # we have only 1 name and 1 array
        names = [names]
        arrs = [arrs]
    arrs = list(map(np.asarray, arrs))
    if dtypes is None:
        dtypes = [a.dtype for a in arrs]
    elif not cbook.iterable(dtypes):
        dtypes = [dtypes]
    if len(arrs) != len(dtypes):
        if len(dtypes) == 1:
            dtypes = dtypes * len(arrs)
        else:
            raise ValueError("dtypes must be None, a single dtype or a list")
    old_dtypes = rec.dtype.descr
    if six.PY2:
        old_dtypes = [(name.encode('utf-8'), dt) for name, dt in old_dtypes]
    newdtype = np.dtype(old_dtypes + list(zip(names, dtypes)))
    newrec = np.recarray(rec.shape, dtype=newdtype)
    for field in rec.dtype.fields:
        newrec[field] = rec[field]
    for name, arr in zip(names, arrs):
        newrec[name] = arr
    return newrec


@cbook.deprecated("2.2")
def rec_drop_fields(rec, names):
    """
    Return a new numpy record array with fields in *names* dropped.
    """

    names = set(names)

    newdtype = np.dtype([(name, rec.dtype[name]) for name in rec.dtype.names
                         if name not in names])

    newrec = np.recarray(rec.shape, dtype=newdtype)
    for field in newdtype.names:
        newrec[field] = rec[field]

    return newrec


@cbook.deprecated("2.2")
def rec_keep_fields(rec, names):
    """
    Return a new numpy record array with only fields listed in names
    """

    if isinstance(names, six.string_types):
        names = names.split(',')

    arrays = []
    for name in names:
        arrays.append(rec[name])

    return np.rec.fromarrays(arrays, names=names)


@cbook.deprecated("2.2")
def rec_groupby(r, groupby, stats):
    """
    *r* is a numpy record array

    *groupby* is a sequence of record array attribute names that
    together form the grouping key.  e.g., ('date', 'productcode')

    *stats* is a sequence of (*attr*, *func*, *outname*) tuples which
    will call ``x = func(attr)`` and assign *x* to the record array
    output with attribute *outname*.  For example::

      stats = ( ('sales', len, 'numsales'), ('sales', np.mean, 'avgsale') )

    Return record array has *dtype* names for each attribute name in
    the *groupby* argument, with the associated group values, and
    for each outname name in the *stats* argument, with the associated
    stat summary output.
    """
    # build a dictionary from groupby keys-> list of indices into r with
    # those keys
    rowd = {}
    for i, row in enumerate(r):
        key = tuple([row[attr] for attr in groupby])
        rowd.setdefault(key, []).append(i)

    rows = []
    # sort the output by groupby keys
    for key in sorted(rowd):
        row = list(key)
        # get the indices for this groupby key
        ind = rowd[key]
        thisr = r[ind]
        # call each stat function for this groupby slice
        row.extend([func(thisr[attr]) for attr, func, outname in stats])
        rows.append(row)

    # build the output record array with groupby and outname attributes
    attrs, funcs, outnames = list(zip(*stats))
    names = list(groupby)
    names.extend(outnames)
    return np.rec.fromrecords(rows, names=names)


@cbook.deprecated("2.2")
def rec_summarize(r, summaryfuncs):
    """
    *r* is a numpy record array

    *summaryfuncs* is a list of (*attr*, *func*, *outname*) tuples
    which will apply *func* to the array *r*[attr] and assign the
    output to a new attribute name *outname*.  The returned record
    array is identical to *r*, with extra arrays for each element in
    *summaryfuncs*.

    """

    names = list(r.dtype.names)
    arrays = [r[name] for name in names]

    for attr, func, outname in summaryfuncs:
        names.append(outname)
        arrays.append(np.asarray(func(r[attr])))

    return np.rec.fromarrays(arrays, names=names)


@cbook.deprecated("2.2")
def rec_join(key, r1, r2, jointype='inner', defaults=None, r1postfix='1',
             r2postfix='2'):
    """
    Join record arrays *r1* and *r2* on *key*; *key* is a tuple of
    field names -- if *key* is a string it is assumed to be a single
    attribute name. If *r1* and *r2* have equal values on all the keys
    in the *key* tuple, then their fields will be merged into a new
    record array containing the intersection of the fields of *r1* and
    *r2*.

    *r1* (also *r2*) must not have any duplicate keys.

    The *jointype* keyword can be 'inner', 'outer', 'leftouter'.  To
    do a rightouter join just reverse *r1* and *r2*.

    The *defaults* keyword is a dictionary filled with
    ``{column_name:default_value}`` pairs.

    The keywords *r1postfix* and *r2postfix* are postfixed to column names
    (other than keys) that are both in *r1* and *r2*.
    """

    if isinstance(key, six.string_types):
        key = (key, )

    for name in key:
        if name not in r1.dtype.names:
            raise ValueError('r1 does not have key field %s' % name)
        if name not in r2.dtype.names:
            raise ValueError('r2 does not have key field %s' % name)

    def makekey(row):
        return tuple([row[name] for name in key])

    r1d = {makekey(row): i for i, row in enumerate(r1)}
    r2d = {makekey(row): i for i, row in enumerate(r2)}

    r1keys = set(r1d)
    r2keys = set(r2d)

    common_keys = r1keys & r2keys

    r1ind = np.array([r1d[k] for k in common_keys])
    r2ind = np.array([r2d[k] for k in common_keys])

    common_len = len(common_keys)
    left_len = right_len = 0
    if jointype == "outer" or jointype == "leftouter":
        left_keys = r1keys.difference(r2keys)
        left_ind = np.array([r1d[k] for k in left_keys])
        left_len = len(left_ind)
    if jointype == "outer":
        right_keys = r2keys.difference(r1keys)
        right_ind = np.array([r2d[k] for k in right_keys])
        right_len = len(right_ind)

    def key_desc(name):
        '''
        if name is a string key, use the larger size of r1 or r2 before
        merging
        '''
        dt1 = r1.dtype[name]
        if dt1.type != np.string_:
            return (name, dt1.descr[0][1])

        dt2 = r2.dtype[name]
        if dt1 != dt2:
            raise ValueError("The '{}' fields in arrays 'r1' and 'r2' must "
                             "have the same dtype".format(name))
        if dt1.num > dt2.num:
            return (name, dt1.descr[0][1])
        else:
            return (name, dt2.descr[0][1])

    keydesc = [key_desc(name) for name in key]

    def mapped_r1field(name):
        """
        The column name in *newrec* that corresponds to the column in *r1*.
        """
        if name in key or name not in r2.dtype.names:
            return name
        else:
            return name + r1postfix

    def mapped_r2field(name):
        """
        The column name in *newrec* that corresponds to the column in *r2*.
        """
        if name in key or name not in r1.dtype.names:
            return name
        else:
            return name + r2postfix

    r1desc = [(mapped_r1field(desc[0]), desc[1]) for desc in r1.dtype.descr
              if desc[0] not in key]
    r2desc = [(mapped_r2field(desc[0]), desc[1]) for desc in r2.dtype.descr
              if desc[0] not in key]
    all_dtypes = keydesc + r1desc + r2desc
    if six.PY2:
        all_dtypes = [(name.encode('utf-8'), dt) for name, dt in all_dtypes]
    newdtype = np.dtype(all_dtypes)
    newrec = np.recarray((common_len + left_len + right_len,), dtype=newdtype)

    if defaults is not None:
        for thiskey in defaults:
            if thiskey not in newdtype.names:
                warnings.warn('rec_join defaults key="%s" not in new dtype '
                              'names "%s"' % (thiskey, newdtype.names))

    for name in newdtype.names:
        dt = newdtype[name]
        if dt.kind in ('f', 'i'):
            newrec[name] = 0

    if jointype != 'inner' and defaults is not None:
        # fill in the defaults enmasse
        newrec_fields = list(newrec.dtype.fields)
        for k, v in six.iteritems(defaults):
            if k in newrec_fields:
                newrec[k] = v

    for field in r1.dtype.names:
        newfield = mapped_r1field(field)
        if common_len:
            newrec[newfield][:common_len] = r1[field][r1ind]
        if (jointype == "outer" or jointype == "leftouter") and left_len:
            newrec[newfield][common_len:(common_len+left_len)] = (
                r1[field][left_ind]
            )

    for field in r2.dtype.names:
        newfield = mapped_r2field(field)
        if field not in key and common_len:
            newrec[newfield][:common_len] = r2[field][r2ind]
        if jointype == "outer" and right_len:
            newrec[newfield][-right_len:] = r2[field][right_ind]

    newrec.sort(order=key)

    return newrec


@cbook.deprecated("2.2")
def recs_join(key, name, recs, jointype='outer', missing=0., postfixes=None):
    """
    Join a sequence of record arrays on single column key.

    This function only joins a single column of the multiple record arrays

    *key*
      is the column name that acts as a key

    *name*
      is the name of the column that we want to join

    *recs*
      is a list of record arrays to join

    *jointype*
      is a string 'inner' or 'outer'

    *missing*
      is what any missing field is replaced by

    *postfixes*
      if not None, a len recs sequence of postfixes

    returns a record array with columns [rowkey, name0, name1, ... namen-1].
    or if postfixes [PF0, PF1, ..., PFN-1] are supplied,
    [rowkey, namePF0, namePF1, ... namePFN-1].

    Example::

      r = recs_join("date", "close", recs=[r0, r1], missing=0.)

    """
    results = []
    aligned_iters = cbook.align_iterators(operator.attrgetter(key),
                                          *[iter(r) for r in recs])

    def extract(r):
        if r is None:
            return missing
        else:
            return r[name]

    if jointype == "outer":
        for rowkey, row in aligned_iters:
            results.append([rowkey] + list(map(extract, row)))
    elif jointype == "inner":
        for rowkey, row in aligned_iters:
            if None not in row:  # throw out any Nones
                results.append([rowkey] + list(map(extract, row)))

    if postfixes is None:
        postfixes = ['%d' % i for i in range(len(recs))]
    names = ",".join([key] + ["%s%s" % (name, postfix)
                              for postfix in postfixes])
    return np.rec.fromrecords(results, names=names)


@cbook.deprecated("2.2")
def csv2rec(fname, comments='#', skiprows=0, checkrows=0, delimiter=',',
            converterd=None, names=None, missing='', missingd=None,
            use_mrecords=False, dayfirst=False, yearfirst=False):
    """
    Load data from comma/space/tab delimited file in *fname* into a
    numpy record array and return the record array.

    If *names* is *None*, a header row is required to automatically
    assign the recarray names.  The headers will be lower cased,
    spaces will be converted to underscores, and illegal attribute
    name characters removed.  If *names* is not *None*, it is a
    sequence of names to use for the column names.  In this case, it
    is assumed there is no header row.


    - *fname*: can be a filename or a file handle.  Support for gzipped
      files is automatic, if the filename ends in '.gz'

    - *comments*: the character used to indicate the start of a comment
      in the file, or *None* to switch off the removal of comments

    - *skiprows*: is the number of rows from the top to skip

    - *checkrows*: is the number of rows to check to validate the column
      data type.  When set to zero all rows are validated.

    - *converterd*: if not *None*, is a dictionary mapping column number or
      munged column name to a converter function.

    - *names*: if not None, is a list of header names.  In this case, no
      header will be read from the file

    - *missingd* is a dictionary mapping munged column names to field values
      which signify that the field does not contain actual data and should
      be masked, e.g., '0000-00-00' or 'unused'

    - *missing*: a string whose value signals a missing field regardless of
      the column it appears in

    - *use_mrecords*: if True, return an mrecords.fromrecords record array if
      any of the data are missing

    - *dayfirst*: default is False so that MM-DD-YY has precedence over
      DD-MM-YY.  See
      http://labix.org/python-dateutil#head-b95ce2094d189a89f80f5ae52a05b4ab7b41af47
      for further information.

    - *yearfirst*: default is False so that MM-DD-YY has precedence over
      YY-MM-DD. See
      http://labix.org/python-dateutil#head-b95ce2094d189a89f80f5ae52a05b4ab7b41af47
      for further information.

      If no rows are found, *None* is returned
    """

    if converterd is None:
        converterd = dict()

    if missingd is None:
        missingd = {}

    import dateutil.parser
    import datetime

    fh = cbook.to_filehandle(fname)

    delimiter = str(delimiter)

    class FH:
        """
        For space-delimited files, we want different behavior than
        comma or tab.  Generally, we want multiple spaces to be
        treated as a single separator, whereas with comma and tab we
        want multiple commas to return multiple (empty) fields.  The
        join/strip trick below effects this.
        """
        def __init__(self, fh):
            self.fh = fh

        def close(self):
            self.fh.close()

        def seek(self, arg):
            self.fh.seek(arg)

        def fix(self, s):
            return ' '.join(s.split())

        def __next__(self):
            return self.fix(next(self.fh))

        def __iter__(self):
            for line in self.fh:
                yield self.fix(line)

    if delimiter == ' ':
        fh = FH(fh)

    reader = csv.reader(fh, delimiter=delimiter)

    def process_skiprows(reader):
        if skiprows:
            for i, row in enumerate(reader):
                if i >= (skiprows-1):
                    break

        return fh, reader

    process_skiprows(reader)

    def ismissing(name, val):
        "Should the value val in column name be masked?"
        return val == missing or val == missingd.get(name) or val == ''

    def with_default_value(func, default):
        def newfunc(name, val):
            if ismissing(name, val):
                return default
            else:
                return func(val)
        return newfunc

    def mybool(x):
        if x == 'True':
            return True
        elif x == 'False':
            return False
        else:
            raise ValueError('invalid bool')

    dateparser = dateutil.parser.parse

    def mydateparser(x):
        # try and return a datetime object
        d = dateparser(x, dayfirst=dayfirst, yearfirst=yearfirst)
        return d

    mydateparser = with_default_value(mydateparser, datetime.datetime(1, 1, 1))

    myfloat = with_default_value(float, np.nan)
    myint = with_default_value(int, -1)
    mystr = with_default_value(str, '')
    mybool = with_default_value(mybool, None)

    def mydate(x):
        # try and return a date object
        d = dateparser(x, dayfirst=dayfirst, yearfirst=yearfirst)

        if d.hour > 0 or d.minute > 0 or d.second > 0:
            raise ValueError('not a date')
        return d.date()
    mydate = with_default_value(mydate, datetime.date(1, 1, 1))

    def get_func(name, item, func):
        # promote functions in this order
        funcs = [mybool, myint, myfloat, mydate, mydateparser, mystr]
        for func in funcs[funcs.index(func):]:
            try:
                func(name, item)
            except Exception:
                continue
            return func
        raise ValueError('Could not find a working conversion function')

    # map column names that clash with builtins -- TODO - extend this list
    itemd = {
        'return': 'return_',
        'file':   'file_',
        'print':  'print_',
        }

    def get_converters(reader, comments):

        converters = None
        i = 0
        for row in reader:
            if (len(row) and comments is not None and
                    row[0].startswith(comments)):
                continue
            if i == 0:
                converters = [mybool]*len(row)
            if checkrows and i > checkrows:
                break
            i += 1

            for j, (name, item) in enumerate(zip(names, row)):
                func = converterd.get(j)
                if func is None:
                    func = converterd.get(name)
                if func is None:
                    func = converters[j]
                    if len(item.strip()):
                        func = get_func(name, item, func)
                else:
                    # how should we handle custom converters and defaults?
                    func = with_default_value(func, None)
                converters[j] = func
        return converters

    # Get header and remove invalid characters
    needheader = names is None

    if needheader:
        for row in reader:
            if (len(row) and comments is not None and
                    row[0].startswith(comments)):
                continue
            headers = row
            break

        # remove these chars
        delete = set(r"""~!@#$%^&*()-=+~\|}[]{';: /?.>,<""")
        delete.add('"')

        names = []
        seen = dict()
        for i, item in enumerate(headers):
            item = item.strip().lower().replace(' ', '_')
            item = ''.join([c for c in item if c not in delete])
            if not len(item):
                item = 'column%d' % i

            item = itemd.get(item, item)
            cnt = seen.get(item, 0)
            if cnt > 0:
                names.append(item + '_%d' % cnt)
            else:
                names.append(item)
            seen[item] = cnt+1

    else:
        if isinstance(names, six.string_types):
            names = [n.strip() for n in names.split(',')]

    # get the converter functions by inspecting checkrows
    converters = get_converters(reader, comments)
    if converters is None:
        raise ValueError('Could not find any valid data in CSV file')

    # reset the reader and start over
    fh.seek(0)
    reader = csv.reader(fh, delimiter=delimiter)
    process_skiprows(reader)

    if needheader:
        while True:
            # skip past any comments and consume one line of column header
            row = next(reader)
            if (len(row) and comments is not None and
                    row[0].startswith(comments)):
                continue
            break

    # iterate over the remaining rows and convert the data to date
    # objects, ints, or floats as appropriate
    rows = []
    rowmasks = []
    for i, row in enumerate(reader):
        if not len(row):
            continue
        if comments is not None and row[0].startswith(comments):
            continue
        # Ensure that the row returned always has the same nr of elements
        row.extend([''] * (len(converters) - len(row)))
        rows.append([func(name, val)
                     for func, name, val in zip(converters, names, row)])
        rowmasks.append([ismissing(name, val)
                         for name, val in zip(names, row)])
    fh.close()

    if not len(rows):
        return None

    if use_mrecords and np.any(rowmasks):
        r = np.ma.mrecords.fromrecords(rows, names=names, mask=rowmasks)
    else:
        r = np.rec.fromrecords(rows, names=names)
    return r


# a series of classes for describing the format intentions of various rec views
@cbook.deprecated("2.2")
class FormatObj(object):
    def tostr(self, x):
        return self.toval(x)

    def toval(self, x):
        return str(x)

    def fromstr(self, s):
        return s

    def __hash__(self):
        """
        override the hash function of any of the formatters, so that we don't
        create duplicate excel format styles
        """
        return hash(self.__class__)


@cbook.deprecated("2.2")
class FormatString(FormatObj):
    def tostr(self, x):
        val = repr(x)
        return val[1:-1]


@cbook.deprecated("2.2")
class FormatFormatStr(FormatObj):
    def __init__(self, fmt):
        self.fmt = fmt

    def tostr(self, x):
        if x is None:
            return 'None'
        return self.fmt % self.toval(x)


@cbook.deprecated("2.2")
class FormatFloat(FormatFormatStr):
    def __init__(self, precision=4, scale=1.):
        FormatFormatStr.__init__(self, '%%1.%df' % precision)
        self.precision = precision
        self.scale = scale

    def __hash__(self):
        return hash((self.__class__, self.precision, self.scale))

    def toval(self, x):
        if x is not None:
            x = x * self.scale
        return x

    def fromstr(self, s):
        return float(s)/self.scale


@cbook.deprecated("2.2")
class FormatInt(FormatObj):

    def tostr(self, x):
        return '%d' % int(x)

    def toval(self, x):
        return int(x)

    def fromstr(self, s):
        return int(s)


@cbook.deprecated("2.2")
class FormatBool(FormatObj):
    def toval(self, x):
        return str(x)

    def fromstr(self, s):
        return bool(s)


@cbook.deprecated("2.2")
class FormatPercent(FormatFloat):
    def __init__(self, precision=4):
        FormatFloat.__init__(self, precision, scale=100.)


@cbook.deprecated("2.2")
class FormatThousands(FormatFloat):
    def __init__(self, precision=4):
        FormatFloat.__init__(self, precision, scale=1e-3)


@cbook.deprecated("2.2")
class FormatMillions(FormatFloat):
    def __init__(self, precision=4):
        FormatFloat.__init__(self, precision, scale=1e-6)


@cbook.deprecated("2.2", alternative='date.strftime')
class FormatDate(FormatObj):
    def __init__(self, fmt):
        self.fmt = fmt

    def __hash__(self):
        return hash((self.__class__, self.fmt))

    def toval(self, x):
        if x is None:
            return 'None'
        return x.strftime(self.fmt)

    def fromstr(self, x):
        import dateutil.parser
        return dateutil.parser.parse(x).date()


@cbook.deprecated("2.2", alternative='datetime.strftime')
class FormatDatetime(FormatDate):
    def __init__(self, fmt='%Y-%m-%d %H:%M:%S'):
        FormatDate.__init__(self, fmt)

    def fromstr(self, x):
        import dateutil.parser
        return dateutil.parser.parse(x)


@cbook.deprecated("2.2")
def get_formatd(r, formatd=None):
    'build a formatd guaranteed to have a key for every dtype name'
    defaultformatd = {
        np.bool_: FormatBool(),
        np.int16: FormatInt(),
        np.int32: FormatInt(),
        np.int64: FormatInt(),
        np.float32: FormatFloat(),
        np.float64: FormatFloat(),
        np.object_: FormatObj(),
        np.string_: FormatString()}

    if formatd is None:
        formatd = dict()

    for i, name in enumerate(r.dtype.names):
        dt = r.dtype[name]
        format = formatd.get(name)
        if format is None:
            format = defaultformatd.get(dt.type, FormatObj())
        formatd[name] = format
    return formatd


@cbook.deprecated("2.2")
def csvformat_factory(format):
    format = copy.deepcopy(format)
    if isinstance(format, FormatFloat):
        format.scale = 1.  # override scaling for storage
        format.fmt = '%r'
    return format


@cbook.deprecated("2.2", alternative='numpy.recarray.tofile')
def rec2txt(r, header=None, padding=3, precision=3, fields=None):
    """
    Returns a textual representation of a record array.

    Parameters
    ----------
    r: numpy recarray

    header: list
        column headers

    padding:
        space between each column

    precision: number of decimal places to use for floats.
        Set to an integer to apply to all floats.  Set to a
        list of integers to apply precision individually.
        Precision for non-floats is simply ignored.

    fields : list
        If not None, a list of field names to print.  fields
        can be a list of strings like ['field1', 'field2'] or a single
        comma separated string like 'field1,field2'

    Examples
    --------

    For ``precision=[0,2,3]``, the output is ::

      ID    Price   Return
      ABC   12.54    0.234
      XYZ    6.32   -0.076
    """

    if fields is not None:
        r = rec_keep_fields(r, fields)

    if cbook.is_numlike(precision):
        precision = [precision]*len(r.dtype)

    def get_type(item, atype=int):
        tdict = {None: int, int: float, float: str}
        try:
            atype(str(item))
        except:
            return get_type(item, tdict[atype])
        return atype

    def get_justify(colname, column, precision):
        ntype = column.dtype

        if np.issubdtype(ntype, np.character):
            fixed_width = int(ntype.str[2:])
            length = max(len(colname), fixed_width)
            return 0, length+padding, "%s"  # left justify

        if np.issubdtype(ntype, np.integer):
            length = max(len(colname),
                         np.max(list(map(len, list(map(str, column))))))
            return 1, length+padding, "%d"  # right justify

        if np.issubdtype(ntype, np.floating):
            fmt = "%." + str(precision) + "f"
            length = max(
                len(colname),
                np.max(list(map(len, list(map(lambda x: fmt % x, column)))))
            )
            return 1, length+padding, fmt   # right justify

        return (0,
                max(len(colname),
                    np.max(list(map(len, list(map(str, column))))))+padding,
                "%s")

    if header is None:
        header = r.dtype.names

    justify_pad_prec = [get_justify(header[i], r.__getitem__(colname),
                                    precision[i])
                        for i, colname in enumerate(r.dtype.names)]

    justify_pad_prec_spacer = []
    for i in range(len(justify_pad_prec)):
        just, pad, prec = justify_pad_prec[i]
        if i == 0:
            justify_pad_prec_spacer.append((just, pad, prec, 0))
        else:
            pjust, ppad, pprec = justify_pad_prec[i-1]
            if pjust == 0 and just == 1:
                justify_pad_prec_spacer.append((just, pad-padding, prec, 0))
            elif pjust == 1 and just == 0:
                justify_pad_prec_spacer.append((just, pad, prec, padding))
            else:
                justify_pad_prec_spacer.append((just, pad, prec, 0))

    def format(item, just_pad_prec_spacer):
        just, pad, prec, spacer = just_pad_prec_spacer
        if just == 0:
            return spacer*' ' + str(item).ljust(pad)
        else:
            if get_type(item) == float:
                item = (prec % float(item))
            elif get_type(item) == int:
                item = (prec % int(item))

            return item.rjust(pad)

    textl = []
    textl.append(''.join([format(colitem, justify_pad_prec_spacer[j])
                          for j, colitem in enumerate(header)]))
    for i, row in enumerate(r):
        textl.append(''.join([format(colitem, justify_pad_prec_spacer[j])
                              for j, colitem in enumerate(row)]))
        if i == 0:
            textl[0] = textl[0].rstrip()

    text = os.linesep.join(textl)
    return text


@cbook.deprecated("2.2", alternative='numpy.recarray.tofile')
def rec2csv(r, fname, delimiter=',', formatd=None, missing='',
            missingd=None, withheader=True):
    """
    Save the data from numpy recarray *r* into a
    comma-/space-/tab-delimited file.  The record array dtype names
    will be used for column headers.

    *fname*: can be a filename or a file handle.  Support for gzipped
      files is automatic, if the filename ends in '.gz'

    *withheader*: if withheader is False, do not write the attribute
      names in the first row

    for formatd type FormatFloat, we override the precision to store
    full precision floats in the CSV file

    See Also
    --------
    :func:`csv2rec`
        For information about *missing* and *missingd*, which can be used to
        fill in masked values into your CSV file.
    """

    delimiter = str(delimiter)

    if missingd is None:
        missingd = dict()

    def with_mask(func):
        def newfunc(val, mask, mval):
            if mask:
                return mval
            else:
                return func(val)
        return newfunc

    if r.ndim != 1:
        raise ValueError('rec2csv only operates on 1 dimensional recarrays')

    formatd = get_formatd(r, formatd)
    funcs = []
    for i, name in enumerate(r.dtype.names):
        funcs.append(with_mask(csvformat_factory(formatd[name]).tostr))

    fh, opened = cbook.to_filehandle(fname, 'wb', return_opened=True)
    writer = csv.writer(fh, delimiter=delimiter)
    header = r.dtype.names
    if withheader:
        writer.writerow(header)

    # Our list of specials for missing values
    mvals = []
    for name in header:
        mvals.append(missingd.get(name, missing))

    ismasked = False
    if len(r):
        row = r[0]
        ismasked = hasattr(row, '_fieldmask')

    for row in r:
        if ismasked:
            row, rowmask = row.item(), row._fieldmask.item()
        else:
            rowmask = [False] * len(row)
        writer.writerow([func(val, mask, mval) for func, val, mask, mval
                         in zip(funcs, row, rowmask, mvals)])
    if opened:
        fh.close()


@cbook.deprecated('2.2')
def griddata(x, y, z, xi, yi, interp='nn'):
    """
    Interpolates from a nonuniformly spaced grid to some other grid.

    Fits a surface of the form z = f(`x`, `y`) to the data in the
    (usually) nonuniformly spaced vectors (`x`, `y`, `z`), then
    interpolates this surface at the points specified by
    (`xi`, `yi`) to produce `zi`.

    Parameters
    ----------
    x, y, z : 1d array_like
        Coordinates of grid points to interpolate from.
    xi, yi : 1d or 2d array_like
        Coordinates of grid points to interpolate to.
    interp : string key from {'nn', 'linear'}
        Interpolation algorithm, either 'nn' for natural neighbor, or
        'linear' for linear interpolation.

    Returns
    -------
    2d float array
        Array of values interpolated at (`xi`, `yi`) points.  Array
        will be masked is any of (`xi`, `yi`) are outside the convex
        hull of (`x`, `y`).

    Notes
    -----
    If `interp` is 'nn' (the default), uses natural neighbor
    interpolation based on Delaunay triangulation.  This option is
    only available if the mpl_toolkits.natgrid module is installed.
    This can be downloaded from https://github.com/matplotlib/natgrid.
    The (`xi`, `yi`) grid must be regular and monotonically increasing
    in this case.

    If `interp` is 'linear', linear interpolation is used via
    matplotlib.tri.LinearTriInterpolator.

    Instead of using `griddata`, more flexible functionality and other
    interpolation options are available using a
    matplotlib.tri.Triangulation and a matplotlib.tri.TriInterpolator.
    """
    # Check input arguments.
    x = np.asanyarray(x, dtype=np.float64)
    y = np.asanyarray(y, dtype=np.float64)
    z = np.asanyarray(z, dtype=np.float64)
    if x.shape != y.shape or x.shape != z.shape or x.ndim != 1:
        raise ValueError("x, y and z must be equal-length 1-D arrays")

    xi = np.asanyarray(xi, dtype=np.float64)
    yi = np.asanyarray(yi, dtype=np.float64)
    if xi.ndim != yi.ndim:
        raise ValueError("xi and yi must be arrays with the same number of "
                         "dimensions (1 or 2)")
    if xi.ndim == 2 and xi.shape != yi.shape:
        raise ValueError("if xi and yi are 2D arrays, they must have the same "
                         "shape")
    if xi.ndim == 1:
        xi, yi = np.meshgrid(xi, yi)

    if interp == 'nn':
        use_nn_interpolation = True
    elif interp == 'linear':
        use_nn_interpolation = False
    else:
        raise ValueError("interp keyword must be one of 'linear' (for linear "
                         "interpolation) or 'nn' (for natural neighbor "
                         "interpolation).  Default is 'nn'.")

    # Remove masked points.
    mask = np.ma.getmask(z)
    if mask is not np.ma.nomask:
        x = x.compress(~mask)
        y = y.compress(~mask)
        z = z.compressed()

    if use_nn_interpolation:
        try:
            from mpl_toolkits.natgrid import _natgrid
        except ImportError:
            raise RuntimeError(
                "To use interp='nn' (Natural Neighbor interpolation) in "
                "griddata, natgrid must be installed. Either install it "
                "from http://github.com/matplotlib/natgrid or use "
                "interp='linear' instead.")

        if xi.ndim == 2:
            # natgrid expects 1D xi and yi arrays.
            xi = xi[0, :]
            yi = yi[:, 0]

        # Override default natgrid internal parameters.
        _natgrid.seti(b'ext', 0)
        _natgrid.setr(b'nul', np.nan)

        if np.min(np.diff(xi)) < 0 or np.min(np.diff(yi)) < 0:
            raise ValueError("Output grid defined by xi,yi must be monotone "
                             "increasing")

        # Allocate array for output (buffer will be overwritten by natgridd)
        zi = np.empty((yi.shape[0], xi.shape[0]), np.float64)

        # Natgrid requires each array to be contiguous rather than e.g. a view
        # that is a non-contiguous slice of another array.  Use numpy.require
        # to deal with this, which will copy if necessary.
        x = np.require(x, requirements=['C'])
        y = np.require(y, requirements=['C'])
        z = np.require(z, requirements=['C'])
        xi = np.require(xi, requirements=['C'])
        yi = np.require(yi, requirements=['C'])
        _natgrid.natgridd(x, y, z, xi, yi, zi)

        # Mask points on grid outside convex hull of input data.
        if np.any(np.isnan(zi)):
            zi = np.ma.masked_where(np.isnan(zi), zi)
        return zi
    else:
        # Linear interpolation performed using a matplotlib.tri.Triangulation
        # and a matplotlib.tri.LinearTriInterpolator.
        from .tri import Triangulation, LinearTriInterpolator
        triang = Triangulation(x, y)
        interpolator = LinearTriInterpolator(triang, z)
        return interpolator(xi, yi)


##################################################
# Linear interpolation algorithms
##################################################
@cbook.deprecated("2.2", alternative="numpy.interp")
def less_simple_linear_interpolation(x, y, xi, extrap=False):
    """
    This function provides simple (but somewhat less so than
    :func:`cbook.simple_linear_interpolation`) linear interpolation.
    :func:`simple_linear_interpolation` will give a list of point
    between a start and an end, while this does true linear
    interpolation at an arbitrary set of points.

    This is very inefficient linear interpolation meant to be used
    only for a small number of points in relatively non-intensive use
    cases.  For real linear interpolation, use scipy.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    xi = np.atleast_1d(xi)

    s = list(y.shape)
    s[0] = len(xi)
    yi = np.tile(np.nan, s)

    for ii, xx in enumerate(xi):
        bb = x == xx
        if np.any(bb):
            jj, = np.nonzero(bb)
            yi[ii] = y[jj[0]]
        elif xx < x[0]:
            if extrap:
                yi[ii] = y[0]
        elif xx > x[-1]:
            if extrap:
                yi[ii] = y[-1]
        else:
            jj, = np.nonzero(x < xx)
            jj = max(jj)

            yi[ii] = y[jj] + (xx-x[jj])/(x[jj+1]-x[jj]) * (y[jj+1]-y[jj])

    return yi


@cbook.deprecated("2.2")
def slopes(x, y):
    """
    :func:`slopes` calculates the slope *y*'(*x*)

    The slope is estimated using the slope obtained from that of a
    parabola through any three consecutive points.

    This method should be superior to that described in the appendix
    of A CONSISTENTLY WELL BEHAVED METHOD OF INTERPOLATION by Russel
    W. Stineman (Creative Computing July 1980) in at least one aspect:

      Circles for interpolation demand a known aspect ratio between
      *x*- and *y*-values.  For many functions, however, the abscissa
      are given in different dimensions, so an aspect ratio is
      completely arbitrary.

    The parabola method gives very similar results to the circle
    method for most regular cases but behaves much better in special
    cases.

    Norbert Nemec, Institute of Theoretical Physics, University or
    Regensburg, April 2006 Norbert.Nemec at physik.uni-regensburg.de

    (inspired by a original implementation by Halldor Bjornsson,
    Icelandic Meteorological Office, March 2006 halldor at vedur.is)
    """
    # Cast key variables as float.
    x = np.asarray(x, float)
    y = np.asarray(y, float)

    yp = np.zeros(y.shape, float)

    dx = x[1:] - x[:-1]
    dy = y[1:] - y[:-1]
    dydx = dy/dx
    yp[1:-1] = (dydx[:-1] * dx[1:] + dydx[1:] * dx[:-1])/(dx[1:] + dx[:-1])
    yp[0] = 2.0 * dy[0]/dx[0] - yp[1]
    yp[-1] = 2.0 * dy[-1]/dx[-1] - yp[-2]
    return yp


@cbook.deprecated("2.2")
def stineman_interp(xi, x, y, yp=None):
    """
    Given data vectors *x* and *y*, the slope vector *yp* and a new
    abscissa vector *xi*, the function :func:`stineman_interp` uses
    Stineman interpolation to calculate a vector *yi* corresponding to
    *xi*.

    Here's an example that generates a coarse sine curve, then
    interpolates over a finer abscissa::

      x = linspace(0,2*pi,20);  y = sin(x); yp = cos(x)
      xi = linspace(0,2*pi,40);
      yi = stineman_interp(xi,x,y,yp);
      plot(x,y,'o',xi,yi)

    The interpolation method is described in the article A
    CONSISTENTLY WELL BEHAVED METHOD OF INTERPOLATION by Russell
    W. Stineman. The article appeared in the July 1980 issue of
    Creative Computing with a note from the editor stating that while
    they were:

      not an academic journal but once in a while something serious
      and original comes in adding that this was
      "apparently a real solution" to a well known problem.

    For *yp* = *None*, the routine automatically determines the slopes
    using the :func:`slopes` routine.

    *x* is assumed to be sorted in increasing order.

    For values ``xi[j] < x[0]`` or ``xi[j] > x[-1]``, the routine
    tries an extrapolation.  The relevance of the data obtained from
    this, of course, is questionable...

    Original implementation by Halldor Bjornsson, Icelandic
    Meteorolocial Office, March 2006 halldor at vedur.is

    Completely reworked and optimized for Python by Norbert Nemec,
    Institute of Theoretical Physics, University or Regensburg, April
    2006 Norbert.Nemec at physik.uni-regensburg.de
    """

    # Cast key variables as float.
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if x.shape != y.shape:
        raise ValueError("'x' and 'y' must be of same shape")

    if yp is None:
        yp = slopes(x, y)
    else:
        yp = np.asarray(yp, float)

    xi = np.asarray(xi, float)
    yi = np.zeros(xi.shape, float)

    # calculate linear slopes
    dx = x[1:] - x[:-1]
    dy = y[1:] - y[:-1]
    s = dy/dx  # note length of s is N-1 so last element is #N-2

    # find the segment each xi is in
    # this line actually is the key to the efficiency of this implementation
    idx = np.searchsorted(x[1:-1], xi)

    # now we have generally: x[idx[j]] <= xi[j] <= x[idx[j]+1]
    # except at the boundaries, where it may be that xi[j] < x[0] or
    # xi[j] > x[-1]

    # the y-values that would come out from a linear interpolation:
    sidx = s.take(idx)
    xidx = x.take(idx)
    yidx = y.take(idx)
    xidxp1 = x.take(idx+1)
    yo = yidx + sidx * (xi - xidx)

    # the difference that comes when using the slopes given in yp
    # using the yp slope of the left point
    dy1 = (yp.take(idx) - sidx) * (xi - xidx)
    # using the yp slope of the right point
    dy2 = (yp.take(idx+1)-sidx) * (xi - xidxp1)

    dy1dy2 = dy1*dy2
    # The following is optimized for Python. The solution actually
    # does more calculations than necessary but exploiting the power
    # of numpy, this is far more efficient than coding a loop by hand
    # in Python
    yi = yo + dy1dy2 * np.choose(np.array(np.sign(dy1dy2), np.int32)+1,
                                 ((2*xi-xidx-xidxp1)/((dy1-dy2)*(xidxp1-xidx)),
                                  0.0,
                                  1/(dy1+dy2),))
    return yi


class GaussianKDE(object):
    """
    Representation of a kernel-density estimate using Gaussian kernels.

    Parameters
    ----------
    dataset : array_like
        Datapoints to estimate from. In case of univariate data this is a 1-D
        array, otherwise a 2-D array with shape (# of dims, # of data).

    bw_method : str, scalar or callable, optional
        The method used to calculate the estimator bandwidth.  This can be
        'scott', 'silverman', a scalar constant or a callable.  If a
        scalar, this will be used directly as `kde.factor`.  If a
        callable, it should take a `GaussianKDE` instance as only
        parameter and return a scalar. If None (default), 'scott' is used.

    Attributes
    ----------
    dataset : ndarray
        The dataset with which `gaussian_kde` was initialized.

    dim : int
        Number of dimensions.

    num_dp : int
        Number of datapoints.

    factor : float
        The bandwidth factor, obtained from `kde.covariance_factor`, with which
        the covariance matrix is multiplied.

    covariance : ndarray
        The covariance matrix of `dataset`, scaled by the calculated bandwidth
        (`kde.factor`).

    inv_cov : ndarray
        The inverse of `covariance`.

    Methods
    -------
    kde.evaluate(points) : ndarray
        Evaluate the estimated pdf on a provided set of points.

    kde(points) : ndarray
        Same as kde.evaluate(points)

    """

    # This implementation with minor modification was too good to pass up.
    # from scipy: https://github.com/scipy/scipy/blob/master/scipy/stats/kde.py

    def __init__(self, dataset, bw_method=None):
        self.dataset = np.atleast_2d(dataset)
        if not np.array(self.dataset).size > 1:
            raise ValueError("`dataset` input should have multiple elements.")

        self.dim, self.num_dp = np.array(self.dataset).shape
        isString = isinstance(bw_method, six.string_types)

        if bw_method is None:
            pass
        elif (isString and bw_method == 'scott'):
            self.covariance_factor = self.scotts_factor
        elif (isString and bw_method == 'silverman'):
            self.covariance_factor = self.silverman_factor
        elif (np.isscalar(bw_method) and not isString):
                self._bw_method = 'use constant'
                self.covariance_factor = lambda: bw_method
        elif callable(bw_method):
            self._bw_method = bw_method
            self.covariance_factor = lambda: self._bw_method(self)
        else:
            raise ValueError("`bw_method` should be 'scott', 'silverman', a "
                             "scalar or a callable")

        # Computes the covariance matrix for each Gaussian kernel using
        # covariance_factor().

        self.factor = self.covariance_factor()
        # Cache covariance and inverse covariance of the data
        if not hasattr(self, '_data_inv_cov'):
            self.data_covariance = np.atleast_2d(
                np.cov(
                    self.dataset,
                    rowvar=1,
                    bias=False))
            self.data_inv_cov = np.linalg.inv(self.data_covariance)

        self.covariance = self.data_covariance * self.factor ** 2
        self.inv_cov = self.data_inv_cov / self.factor ** 2
        self.norm_factor = np.sqrt(
            np.linalg.det(
                2 * np.pi * self.covariance)) * self.num_dp

    def scotts_factor(self):
        return np.power(self.num_dp, -1. / (self.dim + 4))

    def silverman_factor(self):
        return np.power(
            self.num_dp * (self.dim + 2.0) / 4.0, -1. / (self.dim + 4))

    #  Default method to calculate bandwidth, can be overwritten by subclass
    covariance_factor = scotts_factor

    def evaluate(self, points):
        """Evaluate the estimated pdf on a set of points.

        Parameters
        ----------
        points : (# of dimensions, # of points)-array
            Alternatively, a (# of dimensions,) vector can be passed in and
            treated as a single point.

        Returns
        -------
        values : (# of points,)-array
            The values at each point.

        Raises
        ------
        ValueError : if the dimensionality of the input points is different
                     than the dimensionality of the KDE.

        """
        points = np.atleast_2d(points)

        dim, num_m = np.array(points).shape
        if dim != self.dim:
            raise ValueError("points have dimension {}, dataset has dimension "
                             "{}".format(dim, self.dim))

        result = np.zeros((num_m,), dtype=float)

        if num_m >= self.num_dp:
            # there are more points than data, so loop over data
            for i in range(self.num_dp):
                diff = self.dataset[:, i, np.newaxis] - points
                tdiff = np.dot(self.inv_cov, diff)
                energy = np.sum(diff * tdiff, axis=0) / 2.0
                result = result + np.exp(-energy)
        else:
            # loop over points
            for i in range(num_m):
                diff = self.dataset - points[:, i, np.newaxis]
                tdiff = np.dot(self.inv_cov, diff)
                energy = np.sum(diff * tdiff, axis=0) / 2.0
                result[i] = np.sum(np.exp(-energy), axis=0)

        result = result / self.norm_factor

        return result

    __call__ = evaluate


##################################################
# Code related to things in and around polygons
##################################################
@cbook.deprecated("2.2")
def inside_poly(points, verts):
    """
    *points* is a sequence of *x*, *y* points.
    *verts* is a sequence of *x*, *y* vertices of a polygon.

    Return value is a sequence of indices into points for the points
    that are inside the polygon.
    """
    # Make a closed polygon path
    poly = Path(verts)

    # Check to see which points are contained within the Path
    return [idx for idx, p in enumerate(points) if poly.contains_point(p)]


@cbook.deprecated("2.2")
def poly_below(xmin, xs, ys):
    """
    Given a sequence of *xs* and *ys*, return the vertices of a
    polygon that has a horizontal base at *xmin* and an upper bound at
    the *ys*.  *xmin* is a scalar.

    Intended for use with :meth:`matplotlib.axes.Axes.fill`, e.g.,::

      xv, yv = poly_below(0, x, y)
      ax.fill(xv, yv)
    """
    if any(isinstance(var, np.ma.MaskedArray) for var in [xs, ys]):
        numpy = np.ma
    else:
        numpy = np

    xs = numpy.asarray(xs)
    ys = numpy.asarray(ys)
    Nx = len(xs)
    Ny = len(ys)
    if Nx != Ny:
        raise ValueError("'xs' and 'ys' must have the same length")
    x = xmin*numpy.ones(2*Nx)
    y = numpy.ones(2*Nx)
    x[:Nx] = xs
    y[:Nx] = ys
    y[Nx:] = ys[::-1]
    return x, y


@cbook.deprecated("2.2")
def poly_between(x, ylower, yupper):
    """
    Given a sequence of *x*, *ylower* and *yupper*, return the polygon
    that fills the regions between them.  *ylower* or *yupper* can be
    scalar or iterable.  If they are iterable, they must be equal in
    length to *x*.

    Return value is *x*, *y* arrays for use with
    :meth:`matplotlib.axes.Axes.fill`.
    """
    if any(isinstance(var, np.ma.MaskedArray) for var in [ylower, yupper, x]):
        numpy = np.ma
    else:
        numpy = np

    Nx = len(x)
    if not cbook.iterable(ylower):
        ylower = ylower*numpy.ones(Nx)

    if not cbook.iterable(yupper):
        yupper = yupper*numpy.ones(Nx)

    x = numpy.concatenate((x, x[::-1]))
    y = numpy.concatenate((yupper, ylower[::-1]))
    return x, y


@cbook.deprecated('2.2')
def is_closed_polygon(X):
    """
    Tests whether first and last object in a sequence are the same.  These are
    presumably coordinates on a polygonal curve, in which case this function
    tests if that curve is closed.
    """
    return np.all(X[0] == X[-1])


@cbook.deprecated("2.2", message='Moved to matplotlib.cbook')
def contiguous_regions(mask):
    """
    return a list of (ind0, ind1) such that mask[ind0:ind1].all() is
    True and we cover all such regions
    """
    return cbook.contiguous_regions(mask)


@cbook.deprecated("2.2")
def cross_from_below(x, threshold):
    """
    return the indices into *x* where *x* crosses some threshold from
    below, e.g., the i's where::

      x[i-1]<threshold and x[i]>=threshold

    Example code::

        import matplotlib.pyplot as plt

        t = np.arange(0.0, 2.0, 0.1)
        s = np.sin(2*np.pi*t)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(t, s, '-o')
        ax.axhline(0.5)
        ax.axhline(-0.5)

        ind = cross_from_below(s, 0.5)
        ax.vlines(t[ind], -1, 1)

        ind = cross_from_above(s, -0.5)
        ax.vlines(t[ind], -1, 1)

        plt.show()

    See Also
    --------
    :func:`cross_from_above` and :func:`contiguous_regions`

    """
    x = np.asarray(x)
    ind = np.nonzero((x[:-1] < threshold) & (x[1:] >= threshold))[0]
    if len(ind):
        return ind+1
    else:
        return ind


@cbook.deprecated("2.2")
def cross_from_above(x, threshold):
    """
    return the indices into *x* where *x* crosses some threshold from
    below, e.g., the i's where::

      x[i-1]>threshold and x[i]<=threshold

    See Also
    --------
    :func:`cross_from_below` and :func:`contiguous_regions`

    """
    x = np.asarray(x)
    ind = np.nonzero((x[:-1] >= threshold) & (x[1:] < threshold))[0]
    if len(ind):
        return ind+1
    else:
        return ind


##################################################
# Vector and path length geometry calculations
##################################################
@cbook.deprecated('2.2')
def vector_lengths(X, P=2., axis=None):
    """
    Finds the length of a set of vectors in *n* dimensions.  This is
    like the :func:`numpy.norm` function for vectors, but has the ability to
    work over a particular axis of the supplied array or matrix.

    Computes ``(sum((x_i)^P))^(1/P)`` for each ``{x_i}`` being the
    elements of *X* along the given axis.  If *axis* is *None*,
    compute over all elements of *X*.
    """
    X = np.asarray(X)
    return (np.sum(X**(P), axis=axis))**(1./P)


@cbook.deprecated('2.2')
def distances_along_curve(X):
    """
    Computes the distance between a set of successive points in *N* dimensions.

    Where *X* is an *M* x *N* array or matrix.  The distances between
    successive rows is computed.  Distance is the standard Euclidean
    distance.
    """
    X = np.diff(X, axis=0)
    return vector_lengths(X, axis=1)


@cbook.deprecated('2.2')
def path_length(X):
    """
    Computes the distance travelled along a polygonal curve in *N* dimensions.

    Where *X* is an *M* x *N* array or matrix.  Returns an array of
    length *M* consisting of the distance along the curve at each point
    (i.e., the rows of *X*).
    """
    X = distances_along_curve(X)
    return np.concatenate((np.zeros(1), np.cumsum(X)))


@cbook.deprecated('2.2')
def quad2cubic(q0x, q0y, q1x, q1y, q2x, q2y):
    """
    Converts a quadratic Bezier curve to a cubic approximation.

    The inputs are the *x* and *y* coordinates of the three control
    points of a quadratic curve, and the output is a tuple of *x* and
    *y* coordinates of the four control points of the cubic curve.
    """
    # TODO: Candidate for deprecation -- no longer used internally

    # c0x, c0y = q0x, q0y
    c1x, c1y = q0x + 2./3. * (q1x - q0x), q0y + 2./3. * (q1y - q0y)
    c2x, c2y = c1x + 1./3. * (q2x - q0x), c1y + 1./3. * (q2y - q0y)
    # c3x, c3y = q2x, q2y
    return q0x, q0y, c1x, c1y, c2x, c2y, q2x, q2y


@cbook.deprecated("2.2")
def offset_line(y, yerr):
    """
    Offsets an array *y* by +/- an error and returns a tuple
    (y - err, y + err).

    The error term can be:

    * A scalar. In this case, the returned tuple is obvious.
    * A vector of the same length as *y*. The quantities y +/- err are computed
      component-wise.
    * A tuple of length 2. In this case, yerr[0] is the error below *y* and
      yerr[1] is error above *y*. For example::

        from pylab import *
        x = linspace(0, 2*pi, num=100, endpoint=True)
        y = sin(x)
        y_minus, y_plus = mlab.offset_line(y, 0.1)
        plot(x, y)
        fill_between(x, ym, y2=yp)
        show()

    """
    if cbook.is_numlike(yerr) or (cbook.iterable(yerr) and
                                  len(yerr) == len(y)):
        ymin = y - yerr
        ymax = y + yerr
    elif len(yerr) == 2:
        ymin, ymax = y - yerr[0], y + yerr[1]
    else:
        raise ValueError("yerr must be scalar, 1xN or 2xN")
    return ymin, ymax
