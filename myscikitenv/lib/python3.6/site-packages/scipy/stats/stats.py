# Copyright 2002 Gary Strangman.  All rights reserved
# Copyright 2002-2016 The SciPy Developers
#
# The original code from Gary Strangman was heavily adapted for
# use in SciPy by Travis Oliphant.  The original code came with the
# following disclaimer:
#
# This software is provided "as-is".  There are no expressed or implied
# warranties of any kind, including, but not limited to, the warranties
# of merchantability and fitness for a given application.  In no event
# shall Gary Strangman be liable for any direct, indirect, incidental,
# special, exemplary or consequential damages (including, but not limited
# to, loss of use, data or profits, or business interruption) however
# caused and on any theory of liability, whether in contract, strict
# liability or tort (including negligence or otherwise) arising in any way
# out of the use of this software, even if advised of the possibility of
# such damage.

"""
A collection of basic statistical functions for Python.  The function
names appear below.

 Some scalar functions defined here are also available in the scipy.special
 package where they work on arbitrary sized arrays.

Disclaimers:  The function list is obviously incomplete and, worse, the
functions are not optimized.  All functions have been tested (some more
so than others), but they are far from bulletproof.  Thus, as with any
free software, no warranty or guarantee is expressed or implied. :-)  A
few extra functions that don't appear in the list below can be found by
interested treasure-hunters.  These functions don't necessarily have
both list and array versions but were deemed useful.

Central Tendency
----------------
.. autosummary::
   :toctree: generated/

    gmean
    hmean
    mode

Moments
-------
.. autosummary::
   :toctree: generated/

    moment
    variation
    skew
    kurtosis
    normaltest

Altered Versions
----------------
.. autosummary::
   :toctree: generated/

    tmean
    tvar
    tstd
    tsem
    describe

Frequency Stats
---------------
.. autosummary::
   :toctree: generated/

    itemfreq
    scoreatpercentile
    percentileofscore
    cumfreq
    relfreq

Variability
-----------
.. autosummary::
   :toctree: generated/

    obrientransform
    sem
    zmap
    zscore
    gstd
    iqr
    median_absolute_deviation

Trimming Functions
------------------
.. autosummary::
   :toctree: generated/

   trimboth
   trim1

Correlation Functions
---------------------
.. autosummary::
   :toctree: generated/

   pearsonr
   fisher_exact
   spearmanr
   pointbiserialr
   kendalltau
   weightedtau
   linregress
   theilslopes

Inferential Stats
-----------------
.. autosummary::
   :toctree: generated/

   ttest_1samp
   ttest_ind
   ttest_ind_from_stats
   ttest_rel
   chisquare
   power_divergence
   ks_2samp
   epps_singleton_2samp
   mannwhitneyu
   ranksums
   wilcoxon
   kruskal
   friedmanchisquare
   brunnermunzel
   combine_pvalues

Statistical Distances
---------------------
.. autosummary::
   :toctree: generated/

   wasserstein_distance
   energy_distance

ANOVA Functions
---------------
.. autosummary::
   :toctree: generated/

   f_oneway

Support Functions
-----------------
.. autosummary::
   :toctree: generated/

   rankdata
   rvs_ratio_uniforms

References
----------
.. [CRCProbStat2000] Zwillinger, D. and Kokoska, S. (2000). CRC Standard
   Probability and Statistics Tables and Formulae. Chapman & Hall: New
   York. 2000.

"""

from __future__ import division, print_function, absolute_import

import warnings
import sys
import math
if sys.version_info.major >= 3 and sys.version_info.minor >= 5:
    from math import gcd
else:
    from fractions import gcd
from collections import namedtuple

import numpy as np
from numpy import array, asarray, ma

from scipy._lib.six import callable, string_types
from scipy._lib._version import NumpyVersion
from scipy._lib._util import _lazywhere
import scipy.special as special
from scipy import linalg
from . import distributions
from . import mstats_basic
from ._stats_mstats_common import (_find_repeats, linregress, theilslopes,
                                   siegelslopes)
from ._stats import _kendall_dis, _toint64, _weightedrankedtau
from ._rvs_sampling import rvs_ratio_uniforms
from ._hypotests import epps_singleton_2samp


__all__ = ['find_repeats', 'gmean', 'hmean', 'mode', 'tmean', 'tvar',
           'tmin', 'tmax', 'tstd', 'tsem', 'moment', 'variation',
           'skew', 'kurtosis', 'describe', 'skewtest', 'kurtosistest',
           'normaltest', 'jarque_bera', 'itemfreq',
           'scoreatpercentile', 'percentileofscore',
           'cumfreq', 'relfreq', 'obrientransform',
           'sem', 'zmap', 'zscore', 'iqr', 'gstd', 'median_absolute_deviation',
           'sigmaclip', 'trimboth', 'trim1', 'trim_mean', 'f_oneway',
           'PearsonRConstantInputWarning', 'PearsonRNearConstantInputWarning',
           'pearsonr', 'fisher_exact', 'spearmanr', 'pointbiserialr',
           'kendalltau', 'weightedtau',
           'linregress', 'siegelslopes', 'theilslopes', 'ttest_1samp',
           'ttest_ind', 'ttest_ind_from_stats', 'ttest_rel', 'kstest',
           'chisquare', 'power_divergence', 'ks_2samp', 'mannwhitneyu',
           'tiecorrect', 'ranksums', 'kruskal', 'friedmanchisquare',
           'rankdata', 'rvs_ratio_uniforms',
           'combine_pvalues', 'wasserstein_distance', 'energy_distance',
           'brunnermunzel', 'epps_singleton_2samp']


def _chk_asarray(a, axis):
    if axis is None:
        a = np.ravel(a)
        outaxis = 0
    else:
        a = np.asarray(a)
        outaxis = axis

    if a.ndim == 0:
        a = np.atleast_1d(a)

    return a, outaxis


def _chk2_asarray(a, b, axis):
    if axis is None:
        a = np.ravel(a)
        b = np.ravel(b)
        outaxis = 0
    else:
        a = np.asarray(a)
        b = np.asarray(b)
        outaxis = axis

    if a.ndim == 0:
        a = np.atleast_1d(a)
    if b.ndim == 0:
        b = np.atleast_1d(b)

    return a, b, outaxis


def _contains_nan(a, nan_policy='propagate'):
    policies = ['propagate', 'raise', 'omit']
    if nan_policy not in policies:
        raise ValueError("nan_policy must be one of {%s}" %
                         ', '.join("'%s'" % s for s in policies))
    try:
        # Calling np.sum to avoid creating a huge array into memory
        # e.g. np.isnan(a).any()
        with np.errstate(invalid='ignore'):
            contains_nan = np.isnan(np.sum(a))
    except TypeError:
        # This can happen when attempting to sum things which are not
        # numbers (e.g. as in the function `mode`). Try an alternative method:
        try:
            contains_nan = np.nan in set(a.ravel())
        except TypeError:
            # Don't know what to do. Fall back to omitting nan values and
            # issue a warning.
            contains_nan = False
            nan_policy = 'omit'
            warnings.warn("The input array could not be properly checked for nan "
                          "values. nan values will be ignored.", RuntimeWarning)

    if contains_nan and nan_policy == 'raise':
        raise ValueError("The input contains nan values")

    return (contains_nan, nan_policy)


def gmean(a, axis=0, dtype=None):
    """
    Compute the geometric mean along the specified axis.

    Return the geometric average of the array elements.
    That is:  n-th root of (x1 * x2 * ... * xn)

    Parameters
    ----------
    a : array_like
        Input array or object that can be converted to an array.
    axis : int or None, optional
        Axis along which the geometric mean is computed. Default is 0.
        If None, compute over the whole array `a`.
    dtype : dtype, optional
        Type of the returned array and of the accumulator in which the
        elements are summed. If dtype is not specified, it defaults to the
        dtype of a, unless a has an integer dtype with a precision less than
        that of the default platform integer. In that case, the default
        platform integer is used.

    Returns
    -------
    gmean : ndarray
        see dtype parameter above

    See Also
    --------
    numpy.mean : Arithmetic average
    numpy.average : Weighted average
    hmean : Harmonic mean

    Notes
    -----
    The geometric average is computed over a single dimension of the input
    array, axis=0 by default, or all values in the array if axis=None.
    float64 intermediate and return values are used for integer inputs.

    Use masked arrays to ignore any non-finite values in the input or that
    arise in the calculations such as Not a Number and infinity because masked
    arrays automatically mask any non-finite values.

    Examples
    --------
    >>> from scipy.stats import gmean
    >>> gmean([1, 4])
    2.0
    >>> gmean([1, 2, 3, 4, 5, 6, 7])
    3.3800151591412964
    """
    if not isinstance(a, np.ndarray):
        # if not an ndarray object attempt to convert it
        log_a = np.log(np.array(a, dtype=dtype))
    elif dtype:
        # Must change the default dtype allowing array type
        if isinstance(a, np.ma.MaskedArray):
            log_a = np.log(np.ma.asarray(a, dtype=dtype))
        else:
            log_a = np.log(np.asarray(a, dtype=dtype))
    else:
        log_a = np.log(a)
    return np.exp(log_a.mean(axis=axis))


def hmean(a, axis=0, dtype=None):
    """
    Calculate the harmonic mean along the specified axis.

    That is:  n / (1/x1 + 1/x2 + ... + 1/xn)

    Parameters
    ----------
    a : array_like
        Input array, masked array or object that can be converted to an array.
    axis : int or None, optional
        Axis along which the harmonic mean is computed. Default is 0.
        If None, compute over the whole array `a`.
    dtype : dtype, optional
        Type of the returned array and of the accumulator in which the
        elements are summed. If `dtype` is not specified, it defaults to the
        dtype of `a`, unless `a` has an integer `dtype` with a precision less
        than that of the default platform integer. In that case, the default
        platform integer is used.

    Returns
    -------
    hmean : ndarray
        see `dtype` parameter above

    See Also
    --------
    numpy.mean : Arithmetic average
    numpy.average : Weighted average
    gmean : Geometric mean

    Notes
    -----
    The harmonic mean is computed over a single dimension of the input
    array, axis=0 by default, or all values in the array if axis=None.
    float64 intermediate and return values are used for integer inputs.

    Use masked arrays to ignore any non-finite values in the input or that
    arise in the calculations such as Not a Number and infinity.

    Examples
    --------
    >>> from scipy.stats import hmean
    >>> hmean([1, 4])
    1.6000000000000001
    >>> hmean([1, 2, 3, 4, 5, 6, 7])
    2.6997245179063363
    """
    if not isinstance(a, np.ndarray):
        a = np.array(a, dtype=dtype)
    if np.all(a > 0):
        # Harmonic mean only defined if greater than zero
        if isinstance(a, np.ma.MaskedArray):
            size = a.count(axis)
        else:
            if axis is None:
                a = a.ravel()
                size = a.shape[0]
            else:
                size = a.shape[axis]
        return size / np.sum(1.0 / a, axis=axis, dtype=dtype)
    else:
        raise ValueError("Harmonic mean only defined if all elements greater "
                         "than zero")


ModeResult = namedtuple('ModeResult', ('mode', 'count'))


def mode(a, axis=0, nan_policy='propagate'):
    """
    Return an array of the modal (most common) value in the passed array.

    If there is more than one such value, only the smallest is returned.
    The bin-count for the modal bins is also returned.

    Parameters
    ----------
    a : array_like
        n-dimensional array of which to find mode(s).
    axis : int or None, optional
        Axis along which to operate. Default is 0. If None, compute over
        the whole array `a`.
    nan_policy : {'propagate', 'raise', 'omit'}, optional
        Defines how to handle when input contains nan. 'propagate' returns nan,
        'raise' throws an error, 'omit' performs the calculations ignoring nan
        values. Default is 'propagate'.

    Returns
    -------
    mode : ndarray
        Array of modal values.
    count : ndarray
        Array of counts for each mode.

    Examples
    --------
    >>> a = np.array([[6, 8, 3, 0],
    ...               [3, 2, 1, 7],
    ...               [8, 1, 8, 4],
    ...               [5, 3, 0, 5],
    ...               [4, 7, 5, 9]])
    >>> from scipy import stats
    >>> stats.mode(a)
    (array([[3, 1, 0, 0]]), array([[1, 1, 1, 1]]))

    To get mode of whole array, specify ``axis=None``:

    >>> stats.mode(a, axis=None)
    (array([3]), array([3]))

    """
    a, axis = _chk_asarray(a, axis)
    if a.size == 0:
        return ModeResult(np.array([]), np.array([]))

    contains_nan, nan_policy = _contains_nan(a, nan_policy)

    if contains_nan and nan_policy == 'omit':
        a = ma.masked_invalid(a)
        return mstats_basic.mode(a, axis)

    if a.dtype == object and np.nan in set(a.ravel()):
        # Fall back to a slower method since np.unique does not work with NaN
        scores = set(np.ravel(a))  # get ALL unique values
        testshape = list(a.shape)
        testshape[axis] = 1
        oldmostfreq = np.zeros(testshape, dtype=a.dtype)
        oldcounts = np.zeros(testshape, dtype=int)

        for score in scores:
            template = (a == score)
            counts = np.expand_dims(np.sum(template, axis), axis)
            mostfrequent = np.where(counts > oldcounts, score, oldmostfreq)
            oldcounts = np.maximum(counts, oldcounts)
            oldmostfreq = mostfrequent

        return ModeResult(mostfrequent, oldcounts)

    def _mode1D(a):
        vals, cnts = np.unique(a, return_counts=True)
        return vals[cnts.argmax()], cnts.max()

    # np.apply_along_axis will convert the _mode1D tuples to a numpy array, casting types in the process
    # This recreates the results without that issue
    # View of a, rotated so the requested axis is last
    in_dims = list(range(a.ndim))
    a_view = np.transpose(a, in_dims[:axis] + in_dims[axis+1:] + [axis])

    inds = np.ndindex(a_view.shape[:-1])
    modes = np.empty(a_view.shape[:-1], dtype=a.dtype)
    counts = np.zeros(a_view.shape[:-1], dtype=np.int)
    for ind in inds:
        modes[ind], counts[ind] = _mode1D(a_view[ind])
    newshape = list(a.shape)
    newshape[axis] = 1
    return ModeResult(modes.reshape(newshape), counts.reshape(newshape))


def _mask_to_limits(a, limits, inclusive):
    """Mask an array for values outside of given limits.

    This is primarily a utility function.

    Parameters
    ----------
    a : array
    limits : (float or None, float or None)
        A tuple consisting of the (lower limit, upper limit).  Values in the
        input array less than the lower limit or greater than the upper limit
        will be masked out. None implies no limit.
    inclusive : (bool, bool)
        A tuple consisting of the (lower flag, upper flag).  These flags
        determine whether values exactly equal to lower or upper are allowed.

    Returns
    -------
    A MaskedArray.

    Raises
    ------
    A ValueError if there are no values within the given limits.
    """
    lower_limit, upper_limit = limits
    lower_include, upper_include = inclusive
    am = ma.MaskedArray(a)
    if lower_limit is not None:
        if lower_include:
            am = ma.masked_less(am, lower_limit)
        else:
            am = ma.masked_less_equal(am, lower_limit)

    if upper_limit is not None:
        if upper_include:
            am = ma.masked_greater(am, upper_limit)
        else:
            am = ma.masked_greater_equal(am, upper_limit)

    if am.count() == 0:
        raise ValueError("No array values within given limits")

    return am


def tmean(a, limits=None, inclusive=(True, True), axis=None):
    """
    Compute the trimmed mean.

    This function finds the arithmetic mean of given values, ignoring values
    outside the given `limits`.

    Parameters
    ----------
    a : array_like
        Array of values.
    limits : None or (lower limit, upper limit), optional
        Values in the input array less than the lower limit or greater than the
        upper limit will be ignored.  When limits is None (default), then all
        values are used.  Either of the limit values in the tuple can also be
        None representing a half-open interval.
    inclusive : (bool, bool), optional
        A tuple consisting of the (lower flag, upper flag).  These flags
        determine whether values exactly equal to the lower or upper limits
        are included.  The default value is (True, True).
    axis : int or None, optional
        Axis along which to compute test. Default is None.

    Returns
    -------
    tmean : float

    See also
    --------
    trim_mean : returns mean after trimming a proportion from both tails.

    Examples
    --------
    >>> from scipy import stats
    >>> x = np.arange(20)
    >>> stats.tmean(x)
    9.5
    >>> stats.tmean(x, (3,17))
    10.0

    """
    a = asarray(a)
    if limits is None:
        return np.mean(a, None)

    am = _mask_to_limits(a.ravel(), limits, inclusive)
    return am.mean(axis=axis)


def tvar(a, limits=None, inclusive=(True, True), axis=0, ddof=1):
    """
    Compute the trimmed variance.

    This function computes the sample variance of an array of values,
    while ignoring values which are outside of given `limits`.

    Parameters
    ----------
    a : array_like
        Array of values.
    limits : None or (lower limit, upper limit), optional
        Values in the input array less than the lower limit or greater than the
        upper limit will be ignored. When limits is None, then all values are
        used. Either of the limit values in the tuple can also be None
        representing a half-open interval.  The default value is None.
    inclusive : (bool, bool), optional
        A tuple consisting of the (lower flag, upper flag).  These flags
        determine whether values exactly equal to the lower or upper limits
        are included.  The default value is (True, True).
    axis : int or None, optional
        Axis along which to operate. Default is 0. If None, compute over the
        whole array `a`.
    ddof : int, optional
        Delta degrees of freedom.  Default is 1.

    Returns
    -------
    tvar : float
        Trimmed variance.

    Notes
    -----
    `tvar` computes the unbiased sample variance, i.e. it uses a correction
    factor ``n / (n - 1)``.

    Examples
    --------
    >>> from scipy import stats
    >>> x = np.arange(20)
    >>> stats.tvar(x)
    35.0
    >>> stats.tvar(x, (3,17))
    20.0

    """
    a = asarray(a)
    a = a.astype(float).ravel()
    if limits is None:
        n = len(a)
        return a.var() * n / (n - 1.)
    am = _mask_to_limits(a, limits, inclusive)
    return np.ma.var(am, ddof=ddof, axis=axis)


def tmin(a, lowerlimit=None, axis=0, inclusive=True, nan_policy='propagate'):
    """
    Compute the trimmed minimum.

    This function finds the miminum value of an array `a` along the
    specified axis, but only considering values greater than a specified
    lower limit.

    Parameters
    ----------
    a : array_like
        array of values
    lowerlimit : None or float, optional
        Values in the input array less than the given limit will be ignored.
        When lowerlimit is None, then all values are used. The default value
        is None.
    axis : int or None, optional
        Axis along which to operate. Default is 0. If None, compute over the
        whole array `a`.
    inclusive : {True, False}, optional
        This flag determines whether values exactly equal to the lower limit
        are included.  The default value is True.
    nan_policy : {'propagate', 'raise', 'omit'}, optional
        Defines how to handle when input contains nan. 'propagate' returns nan,
        'raise' throws an error, 'omit' performs the calculations ignoring nan
        values. Default is 'propagate'.

    Returns
    -------
    tmin : float, int or ndarray

    Examples
    --------
    >>> from scipy import stats
    >>> x = np.arange(20)
    >>> stats.tmin(x)
    0

    >>> stats.tmin(x, 13)
    13

    >>> stats.tmin(x, 13, inclusive=False)
    14

    """
    a, axis = _chk_asarray(a, axis)
    am = _mask_to_limits(a, (lowerlimit, None), (inclusive, False))

    contains_nan, nan_policy = _contains_nan(am, nan_policy)

    if contains_nan and nan_policy == 'omit':
        am = ma.masked_invalid(am)

    res = ma.minimum.reduce(am, axis).data
    if res.ndim == 0:
        return res[()]
    return res


def tmax(a, upperlimit=None, axis=0, inclusive=True, nan_policy='propagate'):
    """
    Compute the trimmed maximum.

    This function computes the maximum value of an array along a given axis,
    while ignoring values larger than a specified upper limit.

    Parameters
    ----------
    a : array_like
        array of values
    upperlimit : None or float, optional
        Values in the input array greater than the given limit will be ignored.
        When upperlimit is None, then all values are used. The default value
        is None.
    axis : int or None, optional
        Axis along which to operate. Default is 0. If None, compute over the
        whole array `a`.
    inclusive : {True, False}, optional
        This flag determines whether values exactly equal to the upper limit
        are included.  The default value is True.
    nan_policy : {'propagate', 'raise', 'omit'}, optional
        Defines how to handle when input contains nan. 'propagate' returns nan,
        'raise' throws an error, 'omit' performs the calculations ignoring nan
        values. Default is 'propagate'.

    Returns
    -------
    tmax : float, int or ndarray

    Examples
    --------
    >>> from scipy import stats
    >>> x = np.arange(20)
    >>> stats.tmax(x)
    19

    >>> stats.tmax(x, 13)
    13

    >>> stats.tmax(x, 13, inclusive=False)
    12

    """
    a, axis = _chk_asarray(a, axis)
    am = _mask_to_limits(a, (None, upperlimit), (False, inclusive))

    contains_nan, nan_policy = _contains_nan(am, nan_policy)

    if contains_nan and nan_policy == 'omit':
        am = ma.masked_invalid(am)

    res = ma.maximum.reduce(am, axis).data
    if res.ndim == 0:
        return res[()]
    return res


def tstd(a, limits=None, inclusive=(True, True), axis=0, ddof=1):
    """
    Compute the trimmed sample standard deviation.

    This function finds the sample standard deviation of given values,
    ignoring values outside the given `limits`.

    Parameters
    ----------
    a : array_like
        array of values
    limits : None or (lower limit, upper limit), optional
        Values in the input array less than the lower limit or greater than the
        upper limit will be ignored. When limits is None, then all values are
        used. Either of the limit values in the tuple can also be None
        representing a half-open interval.  The default value is None.
    inclusive : (bool, bool), optional
        A tuple consisting of the (lower flag, upper flag).  These flags
        determine whether values exactly equal to the lower or upper limits
        are included.  The default value is (True, True).
    axis : int or None, optional
        Axis along which to operate. Default is 0. If None, compute over the
        whole array `a`.
    ddof : int, optional
        Delta degrees of freedom.  Default is 1.

    Returns
    -------
    tstd : float

    Notes
    -----
    `tstd` computes the unbiased sample standard deviation, i.e. it uses a
    correction factor ``n / (n - 1)``.

    Examples
    --------
    >>> from scipy import stats
    >>> x = np.arange(20)
    >>> stats.tstd(x)
    5.9160797830996161
    >>> stats.tstd(x, (3,17))
    4.4721359549995796

    """
    return np.sqrt(tvar(a, limits, inclusive, axis, ddof))


def tsem(a, limits=None, inclusive=(True, True), axis=0, ddof=1):
    """
    Compute the trimmed standard error of the mean.

    This function finds the standard error of the mean for given
    values, ignoring values outside the given `limits`.

    Parameters
    ----------
    a : array_like
        array of values
    limits : None or (lower limit, upper limit), optional
        Values in the input array less than the lower limit or greater than the
        upper limit will be ignored. When limits is None, then all values are
        used. Either of the limit values in the tuple can also be None
        representing a half-open interval.  The default value is None.
    inclusive : (bool, bool), optional
        A tuple consisting of the (lower flag, upper flag).  These flags
        determine whether values exactly equal to the lower or upper limits
        are included.  The default value is (True, True).
    axis : int or None, optional
        Axis along which to operate. Default is 0. If None, compute over the
        whole array `a`.
    ddof : int, optional
        Delta degrees of freedom.  Default is 1.

    Returns
    -------
    tsem : float

    Notes
    -----
    `tsem` uses unbiased sample standard deviation, i.e. it uses a
    correction factor ``n / (n - 1)``.

    Examples
    --------
    >>> from scipy import stats
    >>> x = np.arange(20)
    >>> stats.tsem(x)
    1.3228756555322954
    >>> stats.tsem(x, (3,17))
    1.1547005383792515

    """
    a = np.asarray(a).ravel()
    if limits is None:
        return a.std(ddof=ddof) / np.sqrt(a.size)

    am = _mask_to_limits(a, limits, inclusive)
    sd = np.sqrt(np.ma.var(am, ddof=ddof, axis=axis))
    return sd / np.sqrt(am.count())


#####################################
#              MOMENTS              #
#####################################

def moment(a, moment=1, axis=0, nan_policy='propagate'):
    r"""
    Calculate the nth moment about the mean for a sample.

    A moment is a specific quantitative measure of the shape of a set of
    points. It is often used to calculate coefficients of skewness and kurtosis
    due to its close relationship with them.


    Parameters
    ----------
    a : array_like
       data
    moment : int or array_like of ints, optional
       order of central moment that is returned. Default is 1.
    axis : int or None, optional
       Axis along which the central moment is computed. Default is 0.
       If None, compute over the whole array `a`.
    nan_policy : {'propagate', 'raise', 'omit'}, optional
        Defines how to handle when input contains nan. 'propagate' returns nan,
        'raise' throws an error, 'omit' performs the calculations ignoring nan
        values. Default is 'propagate'.

    Returns
    -------
    n-th central moment : ndarray or float
       The appropriate moment along the given axis or over all values if axis
       is None. The denominator for the moment calculation is the number of
       observations, no degrees of freedom correction is done.

    See also
    --------
    kurtosis, skew, describe

    Notes
    -----
    The k-th central moment of a data sample is:

    .. math::

        m_k = \frac{1}{n} \sum_{i = 1}^n (x_i - \bar{x})^k

    Where n is the number of samples and x-bar is the mean. This function uses
    exponentiation by squares [1]_ for efficiency.

    References
    ----------
    .. [1] https://eli.thegreenplace.net/2009/03/21/efficient-integer-exponentiation-algorithms

    Examples
    --------
    >>> from scipy.stats import moment
    >>> moment([1, 2, 3, 4, 5], moment=1)
    0.0
    >>> moment([1, 2, 3, 4, 5], moment=2)
    2.0
    """
    a, axis = _chk_asarray(a, axis)

    contains_nan, nan_policy = _contains_nan(a, nan_policy)

    if contains_nan and nan_policy == 'omit':
        a = ma.masked_invalid(a)
        return mstats_basic.moment(a, moment, axis)

    if a.size == 0:
        # empty array, return nan(s) with shape matching `moment`
        if np.isscalar(moment):
            return np.nan
        else:
            return np.full(np.asarray(moment).shape, np.nan, dtype=np.float64)

    # for array_like moment input, return a value for each.
    if not np.isscalar(moment):
        mmnt = [_moment(a, i, axis) for i in moment]
        return np.array(mmnt)
    else:
        return _moment(a, moment, axis)


def _moment(a, moment, axis):
    if np.abs(moment - np.round(moment)) > 0:
        raise ValueError("All moment parameters must be integers")

    if moment == 0:
        # When moment equals 0, the result is 1, by definition.
        shape = list(a.shape)
        del shape[axis]
        if shape:
            # return an actual array of the appropriate shape
            return np.ones(shape, dtype=float)
        else:
            # the input was 1D, so return a scalar instead of a rank-0 array
            return 1.0

    elif moment == 1:
        # By definition the first moment about the mean is 0.
        shape = list(a.shape)
        del shape[axis]
        if shape:
            # return an actual array of the appropriate shape
            return np.zeros(shape, dtype=float)
        else:
            # the input was 1D, so return a scalar instead of a rank-0 array
            return np.float64(0.0)
    else:
        # Exponentiation by squares: form exponent sequence
        n_list = [moment]
        current_n = moment
        while current_n > 2:
            if current_n % 2:
                current_n = (current_n - 1) / 2
            else:
                current_n /= 2
            n_list.append(current_n)

        # Starting point for exponentiation by squares
        a_zero_mean = a - np.expand_dims(np.mean(a, axis), axis)
        if n_list[-1] == 1:
            s = a_zero_mean.copy()
        else:
            s = a_zero_mean**2

        # Perform multiplications
        for n in n_list[-2::-1]:
            s = s**2
            if n % 2:
                s *= a_zero_mean
        return np.mean(s, axis)


def variation(a, axis=0, nan_policy='propagate'):
    """
    Compute the coefficient of variation, the ratio of the biased standard
    deviation to the mean.

    Parameters
    ----------
    a : array_like
        Input array.
    axis : int or None, optional
        Axis along which to calculate the coefficient of variation. Default
        is 0. If None, compute over the whole array `a`.
    nan_policy : {'propagate', 'raise', 'omit'}, optional
        Defines how to handle when input contains nan. 'propagate' returns nan,
        'raise' throws an error, 'omit' performs the calculations ignoring nan
        values. Default is 'propagate'.

    Returns
    -------
    variation : ndarray
        The calculated variation along the requested axis.

    References
    ----------
    .. [1] Zwillinger, D. and Kokoska, S. (2000). CRC Standard
       Probability and Statistics Tables and Formulae. Chapman & Hall: New
       York. 2000.

    Examples
    --------
    >>> from scipy.stats import variation
    >>> variation([1, 2, 3, 4, 5])
    0.47140452079103173
    """
    a, axis = _chk_asarray(a, axis)

    contains_nan, nan_policy = _contains_nan(a, nan_policy)

    if contains_nan and nan_policy == 'omit':
        a = ma.masked_invalid(a)
        return mstats_basic.variation(a, axis)

    return a.std(axis) / a.mean(axis)


def skew(a, axis=0, bias=True, nan_policy='propagate'):
    """
    Compute the sample skewness of a data set.

    For normally distributed data, the skewness should be about 0. For
    unimodal continuous distributions, a skewness value > 0 means that
    there is more weight in the right tail of the distribution. The
    function `skewtest` can be used to determine if the skewness value
    is close enough to 0, statistically speaking.

    Parameters
    ----------
    a : ndarray
        data
    axis : int or None, optional
        Axis along which skewness is calculated. Default is 0.
        If None, compute over the whole array `a`.
    bias : bool, optional
        If False, then the calculations are corrected for statistical bias.
    nan_policy : {'propagate', 'raise', 'omit'}, optional
        Defines how to handle when input contains nan. 'propagate' returns nan,
        'raise' throws an error, 'omit' performs the calculations ignoring nan
        values. Default is 'propagate'.

    Returns
    -------
    skewness : ndarray
        The skewness of values along an axis, returning 0 where all values are
        equal.

    Notes
    -----
    The sample skewness is computed as the Fisher-Pearson coefficient
    of skewness, i.e.

    .. math:: g_1=\\frac{m_3}{m_2^{3/2}}

    where

    .. math:: m_i=\\frac{1}{N}\\sum_{n=1}^N(x[n]-\\bar{x})^i

    is the biased sample :math:`i\\texttt{th}` central moment, and :math:`\\bar{x}` is
    the sample mean.  If ``bias`` is False, the calculations are
    corrected for bias and the value computed is the adjusted
    Fisher-Pearson standardized moment coefficient, i.e.

    .. math:: G_1=\\frac{k_3}{k_2^{3/2}}=
                  \\frac{\\sqrt{N(N-1)}}{N-2}\\frac{m_3}{m_2^{3/2}}.

    References
    ----------

    .. [1] Zwillinger, D. and Kokoska, S. (2000). CRC Standard
       Probability and Statistics Tables and Formulae. Chapman & Hall: New
       York. 2000.
       Section 2.2.24.1

    Examples
    --------
    >>> from scipy.stats import skew
    >>> skew([1, 2, 3, 4, 5])
    0.0
    >>> skew([2, 8, 0, 4, 1, 9, 9, 0])
    0.2650554122698573
    """
    a, axis = _chk_asarray(a, axis)
    n = a.shape[axis]

    contains_nan, nan_policy = _contains_nan(a, nan_policy)

    if contains_nan and nan_policy == 'omit':
        a = ma.masked_invalid(a)
        return mstats_basic.skew(a, axis, bias)

    m2 = moment(a, 2, axis)
    m3 = moment(a, 3, axis)
    zero = (m2 == 0)
    vals = _lazywhere(~zero, (m2, m3),
                      lambda m2, m3: m3 / m2**1.5,
                      0.)
    if not bias:
        can_correct = (n > 2) & (m2 > 0)
        if can_correct.any():
            m2 = np.extract(can_correct, m2)
            m3 = np.extract(can_correct, m3)
            nval = np.sqrt((n - 1.0) * n) / (n - 2.0) * m3 / m2**1.5
            np.place(vals, can_correct, nval)

    if vals.ndim == 0:
        return vals.item()

    return vals


def kurtosis(a, axis=0, fisher=True, bias=True, nan_policy='propagate'):
    """
    Compute the kurtosis (Fisher or Pearson) of a dataset.

    Kurtosis is the fourth central moment divided by the square of the
    variance. If Fisher's definition is used, then 3.0 is subtracted from
    the result to give 0.0 for a normal distribution.

    If bias is False then the kurtosis is calculated using k statistics to
    eliminate bias coming from biased moment estimators

    Use `kurtosistest` to see if result is close enough to normal.

    Parameters
    ----------
    a : array
        data for which the kurtosis is calculated
    axis : int or None, optional
        Axis along which the kurtosis is calculated. Default is 0.
        If None, compute over the whole array `a`.
    fisher : bool, optional
        If True, Fisher's definition is used (normal ==> 0.0). If False,
        Pearson's definition is used (normal ==> 3.0).
    bias : bool, optional
        If False, then the calculations are corrected for statistical bias.
    nan_policy : {'propagate', 'raise', 'omit'}, optional
        Defines how to handle when input contains nan. 'propagate' returns nan,
        'raise' throws an error, 'omit' performs the calculations ignoring nan
        values. Default is 'propagate'.

    Returns
    -------
    kurtosis : array
        The kurtosis of values along an axis. If all values are equal,
        return -3 for Fisher's definition and 0 for Pearson's definition.

    References
    ----------
    .. [1] Zwillinger, D. and Kokoska, S. (2000). CRC Standard
       Probability and Statistics Tables and Formulae. Chapman & Hall: New
       York. 2000.

    Examples
    --------
    >>> from scipy.stats import kurtosis
    >>> kurtosis([1, 2, 3, 4, 5])
    -1.3
    """
    a, axis = _chk_asarray(a, axis)

    contains_nan, nan_policy = _contains_nan(a, nan_policy)

    if contains_nan and nan_policy == 'omit':
        a = ma.masked_invalid(a)
        return mstats_basic.kurtosis(a, axis, fisher, bias)

    n = a.shape[axis]
    m2 = moment(a, 2, axis)
    m4 = moment(a, 4, axis)
    zero = (m2 == 0)
    olderr = np.seterr(all='ignore')
    try:
        vals = np.where(zero, 0, m4 / m2**2.0)
    finally:
        np.seterr(**olderr)

    if not bias:
        can_correct = (n > 3) & (m2 > 0)
        if can_correct.any():
            m2 = np.extract(can_correct, m2)
            m4 = np.extract(can_correct, m4)
            nval = 1.0/(n-2)/(n-3) * ((n**2-1.0)*m4/m2**2.0 - 3*(n-1)**2.0)
            np.place(vals, can_correct, nval + 3.0)

    if vals.ndim == 0:
        vals = vals.item()  # array scalar

    return vals - 3 if fisher else vals


DescribeResult = namedtuple('DescribeResult',
                            ('nobs', 'minmax', 'mean', 'variance', 'skewness',
                             'kurtosis'))


def describe(a, axis=0, ddof=1, bias=True, nan_policy='propagate'):
    """
    Compute several descriptive statistics of the passed array.

    Parameters
    ----------
    a : array_like
       Input data.
    axis : int or None, optional
       Axis along which statistics are calculated. Default is 0.
       If None, compute over the whole array `a`.
    ddof : int, optional
        Delta degrees of freedom (only for variance).  Default is 1.
    bias : bool, optional
        If False, then the skewness and kurtosis calculations are corrected for
        statistical bias.
    nan_policy : {'propagate', 'raise', 'omit'}, optional
        Defines how to handle when input contains nan. 'propagate' returns nan,
        'raise' throws an error, 'omit' performs the calculations ignoring nan
        values. Default is 'propagate'.

    Returns
    -------
    nobs : int or ndarray of ints
       Number of observations (length of data along `axis`).
       When 'omit' is chosen as nan_policy, each column is counted separately.
    minmax: tuple of ndarrays or floats
       Minimum and maximum value of data array.
    mean : ndarray or float
       Arithmetic mean of data along axis.
    variance : ndarray or float
       Unbiased variance of the data along axis, denominator is number of
       observations minus one.
    skewness : ndarray or float
       Skewness, based on moment calculations with denominator equal to
       the number of observations, i.e. no degrees of freedom correction.
    kurtosis : ndarray or float
       Kurtosis (Fisher).  The kurtosis is normalized so that it is
       zero for the normal distribution.  No degrees of freedom are used.

    See Also
    --------
    skew, kurtosis

    Examples
    --------
    >>> from scipy import stats
    >>> a = np.arange(10)
    >>> stats.describe(a)
    DescribeResult(nobs=10, minmax=(0, 9), mean=4.5, variance=9.166666666666666,
                   skewness=0.0, kurtosis=-1.2242424242424244)
    >>> b = [[1, 2], [3, 4]]
    >>> stats.describe(b)
    DescribeResult(nobs=2, minmax=(array([1, 2]), array([3, 4])),
                   mean=array([2., 3.]), variance=array([2., 2.]),
                   skewness=array([0., 0.]), kurtosis=array([-2., -2.]))

    """
    a, axis = _chk_asarray(a, axis)

    contains_nan, nan_policy = _contains_nan(a, nan_policy)

    if contains_nan and nan_policy == 'omit':
        a = ma.masked_invalid(a)
        return mstats_basic.describe(a, axis, ddof, bias)

    if a.size == 0:
        raise ValueError("The input must not be empty.")
    n = a.shape[axis]
    mm = (np.min(a, axis=axis), np.max(a, axis=axis))
    m = np.mean(a, axis=axis)
    v = np.var(a, axis=axis, ddof=ddof)
    sk = skew(a, axis, bias=bias)
    kurt = kurtosis(a, axis, bias=bias)

    return DescribeResult(n, mm, m, v, sk, kurt)

#####################################
#         NORMALITY TESTS           #
#####################################


SkewtestResult = namedtuple('SkewtestResult', ('statistic', 'pvalue'))


def skewtest(a, axis=0, nan_policy='propagate'):
    """
    Test whether the skew is different from the normal distribution.

    This function tests the null hypothesis that the skewness of
    the population that the sample was drawn from is the same
    as that of a corresponding normal distribution.

    Parameters
    ----------
    a : array
        The data to be tested
    axis : int or None, optional
       Axis along which statistics are calculated. Default is 0.
       If None, compute over the whole array `a`.
    nan_policy : {'propagate', 'raise', 'omit'}, optional
        Defines how to handle when input contains nan. 'propagate' returns nan,
        'raise' throws an error, 'omit' performs the calculations ignoring nan
        values. Default is 'propagate'.

    Returns
    -------
    statistic : float
        The computed z-score for this test.
    pvalue : float
        a 2-sided p-value for the hypothesis test

    Notes
    -----
    The sample size must be at least 8.

    References
    ----------
    .. [1] R. B. D'Agostino, A. J. Belanger and R. B. D'Agostino Jr.,
            "A suggestion for using powerful and informative tests of
            normality", American Statistician 44, pp. 316-321, 1990.

    Examples
    --------
    >>> from scipy.stats import skewtest
    >>> skewtest([1, 2, 3, 4, 5, 6, 7, 8])
    SkewtestResult(statistic=1.0108048609177787, pvalue=0.3121098361421897)
    >>> skewtest([2, 8, 0, 4, 1, 9, 9, 0])
    SkewtestResult(statistic=0.44626385374196975, pvalue=0.6554066631275459)
    >>> skewtest([1, 2, 3, 4, 5, 6, 7, 8000])
    SkewtestResult(statistic=3.571773510360407, pvalue=0.0003545719905823133)
    >>> skewtest([100, 100, 100, 100, 100, 100, 100, 101])
    SkewtestResult(statistic=3.5717766638478072, pvalue=0.000354567720281634)
    """
    a, axis = _chk_asarray(a, axis)

    contains_nan, nan_policy = _contains_nan(a, nan_policy)

    if contains_nan and nan_policy == 'omit':
        a = ma.masked_invalid(a)
        return mstats_basic.skewtest(a, axis)

    if axis is None:
        a = np.ravel(a)
        axis = 0
    b2 = skew(a, axis)
    n = a.shape[axis]
    if n < 8:
        raise ValueError(
            "skewtest is not valid with less than 8 samples; %i samples"
            " were given." % int(n))
    y = b2 * math.sqrt(((n + 1) * (n + 3)) / (6.0 * (n - 2)))
    beta2 = (3.0 * (n**2 + 27*n - 70) * (n+1) * (n+3) /
             ((n-2.0) * (n+5) * (n+7) * (n+9)))
    W2 = -1 + math.sqrt(2 * (beta2 - 1))
    delta = 1 / math.sqrt(0.5 * math.log(W2))
    alpha = math.sqrt(2.0 / (W2 - 1))
    y = np.where(y == 0, 1, y)
    Z = delta * np.log(y / alpha + np.sqrt((y / alpha)**2 + 1))

    return SkewtestResult(Z, 2 * distributions.norm.sf(np.abs(Z)))


KurtosistestResult = namedtuple('KurtosistestResult', ('statistic', 'pvalue'))


def kurtosistest(a, axis=0, nan_policy='propagate'):
    """
    Test whether a dataset has normal kurtosis.

    This function tests the null hypothesis that the kurtosis
    of the population from which the sample was drawn is that
    of the normal distribution: ``kurtosis = 3(n-1)/(n+1)``.

    Parameters
    ----------
    a : array
        array of the sample data
    axis : int or None, optional
       Axis along which to compute test. Default is 0. If None,
       compute over the whole array `a`.
    nan_policy : {'propagate', 'raise', 'omit'}, optional
        Defines how to handle when input contains nan. 'propagate' returns nan,
        'raise' throws an error, 'omit' performs the calculations ignoring nan
        values. Default is 'propagate'.

    Returns
    -------
    statistic : float
        The computed z-score for this test.
    pvalue : float
        The 2-sided p-value for the hypothesis test

    Notes
    -----
    Valid only for n>20. This function uses the method described in [1]_.

    References
    ----------
    .. [1] see e.g. F. J. Anscombe, W. J. Glynn, "Distribution of the kurtosis
       statistic b2 for normal samples", Biometrika, vol. 70, pp. 227-234, 1983.

    Examples
    --------
    >>> from scipy.stats import kurtosistest
    >>> kurtosistest(list(range(20)))
    KurtosistestResult(statistic=-1.7058104152122062, pvalue=0.08804338332528348)

    >>> np.random.seed(28041990)
    >>> s = np.random.normal(0, 1, 1000)
    >>> kurtosistest(s)
    KurtosistestResult(statistic=1.2317590987707365, pvalue=0.21803908613450895)
    """
    a, axis = _chk_asarray(a, axis)

    contains_nan, nan_policy = _contains_nan(a, nan_policy)

    if contains_nan and nan_policy == 'omit':
        a = ma.masked_invalid(a)
        return mstats_basic.kurtosistest(a, axis)

    n = a.shape[axis]
    if n < 5:
        raise ValueError(
            "kurtosistest requires at least 5 observations; %i observations"
            " were given." % int(n))
    if n < 20:
        warnings.warn("kurtosistest only valid for n>=20 ... continuing "
                      "anyway, n=%i" % int(n))
    b2 = kurtosis(a, axis, fisher=False)

    E = 3.0*(n-1) / (n+1)
    varb2 = 24.0*n*(n-2)*(n-3) / ((n+1)*(n+1.)*(n+3)*(n+5))  # [1]_ Eq. 1
    x = (b2-E) / np.sqrt(varb2)  # [1]_ Eq. 4
    # [1]_ Eq. 2:
    sqrtbeta1 = 6.0*(n*n-5*n+2)/((n+7)*(n+9)) * np.sqrt((6.0*(n+3)*(n+5)) /
                                                        (n*(n-2)*(n-3)))
    # [1]_ Eq. 3:
    A = 6.0 + 8.0/sqrtbeta1 * (2.0/sqrtbeta1 + np.sqrt(1+4.0/(sqrtbeta1**2)))
    term1 = 1 - 2/(9.0*A)
    denom = 1 + x*np.sqrt(2/(A-4.0))
    term2 = np.sign(denom) * np.where(denom == 0.0, np.nan,
                                      np.power((1-2.0/A)/np.abs(denom), 1/3.0))
    if np.any(denom == 0):
        msg = "Test statistic not defined in some cases due to division by " \
              "zero. Return nan in that case..."
        warnings.warn(msg, RuntimeWarning)

    Z = (term1 - term2) / np.sqrt(2/(9.0*A))  # [1]_ Eq. 5
    if Z.ndim == 0:
        Z = Z[()]

    # zprob uses upper tail, so Z needs to be positive
    return KurtosistestResult(Z, 2 * distributions.norm.sf(np.abs(Z)))


NormaltestResult = namedtuple('NormaltestResult', ('statistic', 'pvalue'))


def normaltest(a, axis=0, nan_policy='propagate'):
    """
    Test whether a sample differs from a normal distribution.

    This function tests the null hypothesis that a sample comes
    from a normal distribution.  It is based on D'Agostino and
    Pearson's [1]_, [2]_ test that combines skew and kurtosis to
    produce an omnibus test of normality.


    Parameters
    ----------
    a : array_like
        The array containing the sample to be tested.
    axis : int or None, optional
        Axis along which to compute test. Default is 0. If None,
        compute over the whole array `a`.
    nan_policy : {'propagate', 'raise', 'omit'}, optional
        Defines how to handle when input contains nan. 'propagate' returns nan,
        'raise' throws an error, 'omit' performs the calculations ignoring nan
        values. Default is 'propagate'.

    Returns
    -------
    statistic : float or array
        ``s^2 + k^2``, where ``s`` is the z-score returned by `skewtest` and
        ``k`` is the z-score returned by `kurtosistest`.
    pvalue : float or array
       A 2-sided chi squared probability for the hypothesis test.

    References
    ----------
    .. [1] D'Agostino, R. B. (1971), "An omnibus test of normality for
           moderate and large sample size", Biometrika, 58, 341-348

    .. [2] D'Agostino, R. and Pearson, E. S. (1973), "Tests for departure from
           normality", Biometrika, 60, 613-622

    Examples
    --------
    >>> from scipy import stats
    >>> pts = 1000
    >>> np.random.seed(28041990)
    >>> a = np.random.normal(0, 1, size=pts)
    >>> b = np.random.normal(2, 1, size=pts)
    >>> x = np.concatenate((a, b))
    >>> k2, p = stats.normaltest(x)
    >>> alpha = 1e-3
    >>> print("p = {:g}".format(p))
    p = 3.27207e-11
    >>> if p < alpha:  # null hypothesis: x comes from a normal distribution
    ...     print("The null hypothesis can be rejected")
    ... else:
    ...     print("The null hypothesis cannot be rejected")
    The null hypothesis can be rejected
    """
    a, axis = _chk_asarray(a, axis)

    contains_nan, nan_policy = _contains_nan(a, nan_policy)

    if contains_nan and nan_policy == 'omit':
        a = ma.masked_invalid(a)
        return mstats_basic.normaltest(a, axis)

    s, _ = skewtest(a, axis)
    k, _ = kurtosistest(a, axis)
    k2 = s*s + k*k

    return NormaltestResult(k2, distributions.chi2.sf(k2, 2))


def jarque_bera(x):
    """
    Perform the Jarque-Bera goodness of fit test on sample data.

    The Jarque-Bera test tests whether the sample data has the skewness and
    kurtosis matching a normal distribution.

    Note that this test only works for a large enough number of data samples
    (>2000) as the test statistic asymptotically has a Chi-squared distribution
    with 2 degrees of freedom.

    Parameters
    ----------
    x : array_like
        Observations of a random variable.

    Returns
    -------
    jb_value : float
        The test statistic.
    p : float
        The p-value for the hypothesis test.

    References
    ----------
    .. [1] Jarque, C. and Bera, A. (1980) "Efficient tests for normality,
           homoscedasticity and serial independence of regression residuals",
           6 Econometric Letters 255-259.

    Examples
    --------
    >>> from scipy import stats
    >>> np.random.seed(987654321)
    >>> x = np.random.normal(0, 1, 100000)
    >>> y = np.random.rayleigh(1, 100000)
    >>> stats.jarque_bera(x)
    (4.7165707989581342, 0.09458225503041906)
    >>> stats.jarque_bera(y)
    (6713.7098548143422, 0.0)

    """
    x = np.asarray(x)
    n = x.size
    if n == 0:
        raise ValueError('At least one observation is required.')

    mu = x.mean()
    diffx = x - mu
    skewness = (1 / n * np.sum(diffx**3)) / (1 / n * np.sum(diffx**2))**(3 / 2.)
    kurtosis = (1 / n * np.sum(diffx**4)) / (1 / n * np.sum(diffx**2))**2
    jb_value = n / 6 * (skewness**2 + (kurtosis - 3)**2 / 4)
    p = 1 - distributions.chi2.cdf(jb_value, 2)

    return jb_value, p


#####################################
#        FREQUENCY FUNCTIONS        #
#####################################

@np.deprecate(message="`itemfreq` is deprecated and will be removed in a "
              "future version. Use instead `np.unique(..., return_counts=True)`")
def itemfreq(a):
    """
    Return a 2-D array of item frequencies.

    Parameters
    ----------
    a : (N,) array_like
        Input array.

    Returns
    -------
    itemfreq : (K, 2) ndarray
        A 2-D frequency table.  Column 1 contains sorted, unique values from
        `a`, column 2 contains their respective counts.

    Examples
    --------
    >>> from scipy import stats
    >>> a = np.array([1, 1, 5, 0, 1, 2, 2, 0, 1, 4])
    >>> stats.itemfreq(a)
    array([[ 0.,  2.],
           [ 1.,  4.],
           [ 2.,  2.],
           [ 4.,  1.],
           [ 5.,  1.]])
    >>> np.bincount(a)
    array([2, 4, 2, 0, 1, 1])

    >>> stats.itemfreq(a/10.)
    array([[ 0. ,  2. ],
           [ 0.1,  4. ],
           [ 0.2,  2. ],
           [ 0.4,  1. ],
           [ 0.5,  1. ]])

    """
    items, inv = np.unique(a, return_inverse=True)
    freq = np.bincount(inv)
    return np.array([items, freq]).T


def scoreatpercentile(a, per, limit=(), interpolation_method='fraction',
                      axis=None):
    """
    Calculate the score at a given percentile of the input sequence.

    For example, the score at `per=50` is the median. If the desired quantile
    lies between two data points, we interpolate between them, according to
    the value of `interpolation`. If the parameter `limit` is provided, it
    should be a tuple (lower, upper) of two values.

    Parameters
    ----------
    a : array_like
        A 1-D array of values from which to extract score.
    per : array_like
        Percentile(s) at which to extract score.  Values should be in range
        [0,100].
    limit : tuple, optional
        Tuple of two scalars, the lower and upper limits within which to
        compute the percentile. Values of `a` outside
        this (closed) interval will be ignored.
    interpolation_method : {'fraction', 'lower', 'higher'}, optional
        This optional parameter specifies the interpolation method to use,
        when the desired quantile lies between two data points `i` and `j`

          - fraction: ``i + (j - i) * fraction`` where ``fraction`` is the
            fractional part of the index surrounded by ``i`` and ``j``.
          - lower: ``i``.
          - higher: ``j``.

    axis : int, optional
        Axis along which the percentiles are computed. Default is None. If
        None, compute over the whole array `a`.

    Returns
    -------
    score : float or ndarray
        Score at percentile(s).

    See Also
    --------
    percentileofscore, numpy.percentile

    Notes
    -----
    This function will become obsolete in the future.
    For NumPy 1.9 and higher, `numpy.percentile` provides all the functionality
    that `scoreatpercentile` provides.  And it's significantly faster.
    Therefore it's recommended to use `numpy.percentile` for users that have
    numpy >= 1.9.

    Examples
    --------
    >>> from scipy import stats
    >>> a = np.arange(100)
    >>> stats.scoreatpercentile(a, 50)
    49.5

    """
    # adapted from NumPy's percentile function.  When we require numpy >= 1.8,
    # the implementation of this function can be replaced by np.percentile.
    a = np.asarray(a)
    if a.size == 0:
        # empty array, return nan(s) with shape matching `per`
        if np.isscalar(per):
            return np.nan
        else:
            return np.full(np.asarray(per).shape, np.nan, dtype=np.float64)

    if limit:
        a = a[(limit[0] <= a) & (a <= limit[1])]

    sorted_ = np.sort(a, axis=axis)
    if axis is None:
        axis = 0

    return _compute_qth_percentile(sorted_, per, interpolation_method, axis)


# handle sequence of per's without calling sort multiple times
def _compute_qth_percentile(sorted_, per, interpolation_method, axis):
    if not np.isscalar(per):
        score = [_compute_qth_percentile(sorted_, i,
                                         interpolation_method, axis)
                 for i in per]
        return np.array(score)

    if not (0 <= per <= 100):
        raise ValueError("percentile must be in the range [0, 100]")

    indexer = [slice(None)] * sorted_.ndim
    idx = per / 100. * (sorted_.shape[axis] - 1)

    if int(idx) != idx:
        # round fractional indices according to interpolation method
        if interpolation_method == 'lower':
            idx = int(np.floor(idx))
        elif interpolation_method == 'higher':
            idx = int(np.ceil(idx))
        elif interpolation_method == 'fraction':
            pass  # keep idx as fraction and interpolate
        else:
            raise ValueError("interpolation_method can only be 'fraction', "
                             "'lower' or 'higher'")

    i = int(idx)
    if i == idx:
        indexer[axis] = slice(i, i + 1)
        weights = array(1)
        sumval = 1.0
    else:
        indexer[axis] = slice(i, i + 2)
        j = i + 1
        weights = array([(j - idx), (idx - i)], float)
        wshape = [1] * sorted_.ndim
        wshape[axis] = 2
        weights.shape = wshape
        sumval = weights.sum()

    # Use np.add.reduce (== np.sum but a little faster) to coerce data type
    return np.add.reduce(sorted_[tuple(indexer)] * weights, axis=axis) / sumval


def percentileofscore(a, score, kind='rank'):
    """
    The percentile rank of a score relative to a list of scores.

    A `percentileofscore` of, for example, 80% means that 80% of the
    scores in `a` are below the given score. In the case of gaps or
    ties, the exact definition depends on the optional keyword, `kind`.

    Parameters
    ----------
    a : array_like
        Array of scores to which `score` is compared.
    score : int or float
        Score that is compared to the elements in `a`.
    kind : {'rank', 'weak', 'strict', 'mean'}, optional
        This optional parameter specifies the interpretation of the
        resulting score:

        - "rank": Average percentage ranking of score.  In case of
                  multiple matches, average the percentage rankings of
                  all matching scores.
        - "weak": This kind corresponds to the definition of a cumulative
                  distribution function.  A percentileofscore of 80%
                  means that 80% of values are less than or equal
                  to the provided score.
        - "strict": Similar to "weak", except that only values that are
                    strictly less than the given score are counted.
        - "mean": The average of the "weak" and "strict" scores, often used in
                  testing.  See

                  https://en.wikipedia.org/wiki/Percentile_rank

    Returns
    -------
    pcos : float
        Percentile-position of score (0-100) relative to `a`.

    See Also
    --------
    numpy.percentile

    Examples
    --------
    Three-quarters of the given values lie below a given score:

    >>> from scipy import stats
    >>> stats.percentileofscore([1, 2, 3, 4], 3)
    75.0

    With multiple matches, note how the scores of the two matches, 0.6
    and 0.8 respectively, are averaged:

    >>> stats.percentileofscore([1, 2, 3, 3, 4], 3)
    70.0

    Only 2/5 values are strictly less than 3:

    >>> stats.percentileofscore([1, 2, 3, 3, 4], 3, kind='strict')
    40.0

    But 4/5 values are less than or equal to 3:

    >>> stats.percentileofscore([1, 2, 3, 3, 4], 3, kind='weak')
    80.0

    The average between the weak and the strict scores is

    >>> stats.percentileofscore([1, 2, 3, 3, 4], 3, kind='mean')
    60.0

    """
    if np.isnan(score):
        return np.nan
    a = np.asarray(a)
    n = len(a)
    if n == 0:
        return 100.0

    if kind == 'rank':
        left = np.count_nonzero(a < score)
        right = np.count_nonzero(a <= score)
        pct = (right + left + (1 if right > left else 0)) * 50.0/n
        return pct
    elif kind == 'strict':
        return np.count_nonzero(a < score) / n * 100
    elif kind == 'weak':
        return np.count_nonzero(a <= score) / n * 100
    elif kind == 'mean':
        pct = (np.count_nonzero(a < score) + np.count_nonzero(a <= score)) / n * 50
        return pct
    else:
        raise ValueError("kind can only be 'rank', 'strict', 'weak' or 'mean'")


HistogramResult = namedtuple('HistogramResult',
                             ('count', 'lowerlimit', 'binsize', 'extrapoints'))


def _histogram(a, numbins=10, defaultlimits=None, weights=None, printextras=False):
    """
    Separate the range into several bins and return the number of instances
    in each bin.

    Parameters
    ----------
    a : array_like
        Array of scores which will be put into bins.
    numbins : int, optional
        The number of bins to use for the histogram. Default is 10.
    defaultlimits : tuple (lower, upper), optional
        The lower and upper values for the range of the histogram.
        If no value is given, a range slightly larger than the range of the
        values in a is used. Specifically ``(a.min() - s, a.max() + s)``,
        where ``s = (1/2)(a.max() - a.min()) / (numbins - 1)``.
    weights : array_like, optional
        The weights for each value in `a`. Default is None, which gives each
        value a weight of 1.0
    printextras : bool, optional
        If True, if there are extra points (i.e. the points that fall outside
        the bin limits) a warning is raised saying how many of those points
        there are.  Default is False.

    Returns
    -------
    count : ndarray
        Number of points (or sum of weights) in each bin.
    lowerlimit : float
        Lowest value of histogram, the lower limit of the first bin.
    binsize : float
        The size of the bins (all bins have the same size).
    extrapoints : int
        The number of points outside the range of the histogram.

    See Also
    --------
    numpy.histogram

    Notes
    -----
    This histogram is based on numpy's histogram but has a larger range by
    default if default limits is not set.

    """
    a = np.ravel(a)
    if defaultlimits is None:
        if a.size == 0:
            # handle empty arrays. Undetermined range, so use 0-1.
            defaultlimits = (0, 1)
        else:
            # no range given, so use values in `a`
            data_min = a.min()
            data_max = a.max()
            # Have bins extend past min and max values slightly
            s = (data_max - data_min) / (2. * (numbins - 1.))
            defaultlimits = (data_min - s, data_max + s)

    # use numpy's histogram method to compute bins
    hist, bin_edges = np.histogram(a, bins=numbins, range=defaultlimits,
                                   weights=weights)
    # hist are not always floats, convert to keep with old output
    hist = np.array(hist, dtype=float)
    # fixed width for bins is assumed, as numpy's histogram gives
    # fixed width bins for int values for 'bins'
    binsize = bin_edges[1] - bin_edges[0]
    # calculate number of extra points
    extrapoints = len([v for v in a
                       if defaultlimits[0] > v or v > defaultlimits[1]])
    if extrapoints > 0 and printextras:
        warnings.warn("Points outside given histogram range = %s"
                      % extrapoints)

    return HistogramResult(hist, defaultlimits[0], binsize, extrapoints)


CumfreqResult = namedtuple('CumfreqResult',
                           ('cumcount', 'lowerlimit', 'binsize',
                            'extrapoints'))


def cumfreq(a, numbins=10, defaultreallimits=None, weights=None):
    """
    Return a cumulative frequency histogram, using the histogram function.

    A cumulative histogram is a mapping that counts the cumulative number of
    observations in all of the bins up to the specified bin.

    Parameters
    ----------
    a : array_like
        Input array.
    numbins : int, optional
        The number of bins to use for the histogram. Default is 10.
    defaultreallimits : tuple (lower, upper), optional
        The lower and upper values for the range of the histogram.
        If no value is given, a range slightly larger than the range of the
        values in `a` is used. Specifically ``(a.min() - s, a.max() + s)``,
        where ``s = (1/2)(a.max() - a.min()) / (numbins - 1)``.
    weights : array_like, optional
        The weights for each value in `a`. Default is None, which gives each
        value a weight of 1.0

    Returns
    -------
    cumcount : ndarray
        Binned values of cumulative frequency.
    lowerlimit : float
        Lower real limit
    binsize : float
        Width of each bin.
    extrapoints : int
        Extra points.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from scipy import stats
    >>> x = [1, 4, 2, 1, 3, 1]
    >>> res = stats.cumfreq(x, numbins=4, defaultreallimits=(1.5, 5))
    >>> res.cumcount
    array([ 1.,  2.,  3.,  3.])
    >>> res.extrapoints
    3

    Create a normal distribution with 1000 random values

    >>> rng = np.random.RandomState(seed=12345)
    >>> samples = stats.norm.rvs(size=1000, random_state=rng)

    Calculate cumulative frequencies

    >>> res = stats.cumfreq(samples, numbins=25)

    Calculate space of values for x

    >>> x = res.lowerlimit + np.linspace(0, res.binsize*res.cumcount.size,
    ...                                  res.cumcount.size)

    Plot histogram and cumulative histogram

    >>> fig = plt.figure(figsize=(10, 4))
    >>> ax1 = fig.add_subplot(1, 2, 1)
    >>> ax2 = fig.add_subplot(1, 2, 2)
    >>> ax1.hist(samples, bins=25)
    >>> ax1.set_title('Histogram')
    >>> ax2.bar(x, res.cumcount, width=res.binsize)
    >>> ax2.set_title('Cumulative histogram')
    >>> ax2.set_xlim([x.min(), x.max()])

    >>> plt.show()

    """
    h, l, b, e = _histogram(a, numbins, defaultreallimits, weights=weights)
    cumhist = np.cumsum(h * 1, axis=0)
    return CumfreqResult(cumhist, l, b, e)


RelfreqResult = namedtuple('RelfreqResult',
                           ('frequency', 'lowerlimit', 'binsize',
                            'extrapoints'))


def relfreq(a, numbins=10, defaultreallimits=None, weights=None):
    """
    Return a relative frequency histogram, using the histogram function.

    A relative frequency  histogram is a mapping of the number of
    observations in each of the bins relative to the total of observations.

    Parameters
    ----------
    a : array_like
        Input array.
    numbins : int, optional
        The number of bins to use for the histogram. Default is 10.
    defaultreallimits : tuple (lower, upper), optional
        The lower and upper values for the range of the histogram.
        If no value is given, a range slightly larger than the range of the
        values in a is used. Specifically ``(a.min() - s, a.max() + s)``,
        where ``s = (1/2)(a.max() - a.min()) / (numbins - 1)``.
    weights : array_like, optional
        The weights for each value in `a`. Default is None, which gives each
        value a weight of 1.0

    Returns
    -------
    frequency : ndarray
        Binned values of relative frequency.
    lowerlimit : float
        Lower real limit
    binsize : float
        Width of each bin.
    extrapoints : int
        Extra points.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from scipy import stats
    >>> a = np.array([2, 4, 1, 2, 3, 2])
    >>> res = stats.relfreq(a, numbins=4)
    >>> res.frequency
    array([ 0.16666667, 0.5       , 0.16666667,  0.16666667])
    >>> np.sum(res.frequency)  # relative frequencies should add up to 1
    1.0

    Create a normal distribution with 1000 random values

    >>> rng = np.random.RandomState(seed=12345)
    >>> samples = stats.norm.rvs(size=1000, random_state=rng)

    Calculate relative frequencies

    >>> res = stats.relfreq(samples, numbins=25)

    Calculate space of values for x

    >>> x = res.lowerlimit + np.linspace(0, res.binsize*res.frequency.size,
    ...                                  res.frequency.size)

    Plot relative frequency histogram

    >>> fig = plt.figure(figsize=(5, 4))
    >>> ax = fig.add_subplot(1, 1, 1)
    >>> ax.bar(x, res.frequency, width=res.binsize)
    >>> ax.set_title('Relative frequency histogram')
    >>> ax.set_xlim([x.min(), x.max()])

    >>> plt.show()

    """
    a = np.asanyarray(a)
    h, l, b, e = _histogram(a, numbins, defaultreallimits, weights=weights)
    h = h / a.shape[0]

    return RelfreqResult(h, l, b, e)


#####################################
#        VARIABILITY FUNCTIONS      #
#####################################

def obrientransform(*args):
    """
    Compute the O'Brien transform on input data (any number of arrays).

    Used to test for homogeneity of variance prior to running one-way stats.
    Each array in ``*args`` is one level of a factor.
    If `f_oneway` is run on the transformed data and found significant,
    the variances are unequal.  From Maxwell and Delaney [1]_, p.112.

    Parameters
    ----------
    args : tuple of array_like
        Any number of arrays.

    Returns
    -------
    obrientransform : ndarray
        Transformed data for use in an ANOVA.  The first dimension
        of the result corresponds to the sequence of transformed
        arrays.  If the arrays given are all 1-D of the same length,
        the return value is a 2-D array; otherwise it is a 1-D array
        of type object, with each element being an ndarray.

    References
    ----------
    .. [1] S. E. Maxwell and H. D. Delaney, "Designing Experiments and
           Analyzing Data: A Model Comparison Perspective", Wadsworth, 1990.

    Examples
    --------
    We'll test the following data sets for differences in their variance.

    >>> x = [10, 11, 13, 9, 7, 12, 12, 9, 10]
    >>> y = [13, 21, 5, 10, 8, 14, 10, 12, 7, 15]

    Apply the O'Brien transform to the data.

    >>> from scipy.stats import obrientransform
    >>> tx, ty = obrientransform(x, y)

    Use `scipy.stats.f_oneway` to apply a one-way ANOVA test to the
    transformed data.

    >>> from scipy.stats import f_oneway
    >>> F, p = f_oneway(tx, ty)
    >>> p
    0.1314139477040335

    If we require that ``p < 0.05`` for significance, we cannot conclude
    that the variances are different.
    """
    TINY = np.sqrt(np.finfo(float).eps)

    # `arrays` will hold the transformed arguments.
    arrays = []

    for arg in args:
        a = np.asarray(arg)
        n = len(a)
        mu = np.mean(a)
        sq = (a - mu)**2
        sumsq = sq.sum()

        # The O'Brien transform.
        t = ((n - 1.5) * n * sq - 0.5 * sumsq) / ((n - 1) * (n - 2))

        # Check that the mean of the transformed data is equal to the
        # original variance.
        var = sumsq / (n - 1)
        if abs(var - np.mean(t)) > TINY:
            raise ValueError('Lack of convergence in obrientransform.')

        arrays.append(t)

    return np.array(arrays)


def sem(a, axis=0, ddof=1, nan_policy='propagate'):
    """
    Calculate the standard error of the mean (or standard error of
    measurement) of the values in the input array.

    Parameters
    ----------
    a : array_like
        An array containing the values for which the standard error is
        returned.
    axis : int or None, optional
        Axis along which to operate. Default is 0. If None, compute over
        the whole array `a`.
    ddof : int, optional
        Delta degrees-of-freedom. How many degrees of freedom to adjust
        for bias in limited samples relative to the population estimate
        of variance. Defaults to 1.
    nan_policy : {'propagate', 'raise', 'omit'}, optional
        Defines how to handle when input contains nan. 'propagate' returns nan,
        'raise' throws an error, 'omit' performs the calculations ignoring nan
        values. Default is 'propagate'.

    Returns
    -------
    s : ndarray or float
        The standard error of the mean in the sample(s), along the input axis.

    Notes
    -----
    The default value for `ddof` is different to the default (0) used by other
    ddof containing routines, such as np.std and np.nanstd.

    Examples
    --------
    Find standard error along the first axis:

    >>> from scipy import stats
    >>> a = np.arange(20).reshape(5,4)
    >>> stats.sem(a)
    array([ 2.8284,  2.8284,  2.8284,  2.8284])

    Find standard error across the whole array, using n degrees of freedom:

    >>> stats.sem(a, axis=None, ddof=0)
    1.2893796958227628

    """
    a, axis = _chk_asarray(a, axis)

    contains_nan, nan_policy = _contains_nan(a, nan_policy)

    if contains_nan and nan_policy == 'omit':
        a = ma.masked_invalid(a)
        return mstats_basic.sem(a, axis, ddof)

    n = a.shape[axis]
    s = np.std(a, axis=axis, ddof=ddof) / np.sqrt(n)
    return s


def zscore(a, axis=0, ddof=0):
    """
    Calculate the z score of each value in the sample, relative to the
    sample mean and standard deviation.

    Parameters
    ----------
    a : array_like
        An array like object containing the sample data.
    axis : int or None, optional
        Axis along which to operate. Default is 0. If None, compute over
        the whole array `a`.
    ddof : int, optional
        Degrees of freedom correction in the calculation of the
        standard deviation. Default is 0.

    Returns
    -------
    zscore : array_like
        The z-scores, standardized by mean and standard deviation of
        input array `a`.

    Notes
    -----
    This function preserves ndarray subclasses, and works also with
    matrices and masked arrays (it uses `asanyarray` instead of
    `asarray` for parameters).

    Examples
    --------
    >>> a = np.array([ 0.7972,  0.0767,  0.4383,  0.7866,  0.8091,
    ...                0.1954,  0.6307,  0.6599,  0.1065,  0.0508])
    >>> from scipy import stats
    >>> stats.zscore(a)
    array([ 1.1273, -1.247 , -0.0552,  1.0923,  1.1664, -0.8559,  0.5786,
            0.6748, -1.1488, -1.3324])

    Computing along a specified axis, using n-1 degrees of freedom
    (``ddof=1``) to calculate the standard deviation:

    >>> b = np.array([[ 0.3148,  0.0478,  0.6243,  0.4608],
    ...               [ 0.7149,  0.0775,  0.6072,  0.9656],
    ...               [ 0.6341,  0.1403,  0.9759,  0.4064],
    ...               [ 0.5918,  0.6948,  0.904 ,  0.3721],
    ...               [ 0.0921,  0.2481,  0.1188,  0.1366]])
    >>> stats.zscore(b, axis=1, ddof=1)
    array([[-0.19264823, -1.28415119,  1.07259584,  0.40420358],
           [ 0.33048416, -1.37380874,  0.04251374,  1.00081084],
           [ 0.26796377, -1.12598418,  1.23283094, -0.37481053],
           [-0.22095197,  0.24468594,  1.19042819, -1.21416216],
           [-0.82780366,  1.4457416 , -0.43867764, -0.1792603 ]])
    """
    a = np.asanyarray(a)
    mns = a.mean(axis=axis)
    sstd = a.std(axis=axis, ddof=ddof)
    if axis and mns.ndim < a.ndim:
        return ((a - np.expand_dims(mns, axis=axis)) /
                np.expand_dims(sstd, axis=axis))
    else:
        return (a - mns) / sstd


def zmap(scores, compare, axis=0, ddof=0):
    """
    Calculate the relative z-scores.

    Return an array of z-scores, i.e., scores that are standardized to
    zero mean and unit variance, where mean and variance are calculated
    from the comparison array.

    Parameters
    ----------
    scores : array_like
        The input for which z-scores are calculated.
    compare : array_like
        The input from which the mean and standard deviation of the
        normalization are taken; assumed to have the same dimension as
        `scores`.
    axis : int or None, optional
        Axis over which mean and variance of `compare` are calculated.
        Default is 0. If None, compute over the whole array `scores`.
    ddof : int, optional
        Degrees of freedom correction in the calculation of the
        standard deviation. Default is 0.

    Returns
    -------
    zscore : array_like
        Z-scores, in the same shape as `scores`.

    Notes
    -----
    This function preserves ndarray subclasses, and works also with
    matrices and masked arrays (it uses `asanyarray` instead of
    `asarray` for parameters).

    Examples
    --------
    >>> from scipy.stats import zmap
    >>> a = [0.5, 2.0, 2.5, 3]
    >>> b = [0, 1, 2, 3, 4]
    >>> zmap(a, b)
    array([-1.06066017,  0.        ,  0.35355339,  0.70710678])
    """
    scores, compare = map(np.asanyarray, [scores, compare])
    mns = compare.mean(axis=axis)
    sstd = compare.std(axis=axis, ddof=ddof)
    if axis and mns.ndim < compare.ndim:
        return ((scores - np.expand_dims(mns, axis=axis)) /
                np.expand_dims(sstd, axis=axis))
    else:
        return (scores - mns) / sstd


def gstd(a, axis=0, ddof=1):
    """Calculate the geometric standard deviation of an array

    The geometric standard deviation describes the spread of a set of numbers
    where the geometric mean is preferred. It is a multiplicative factor, and
    so a dimensionless quantity.

    It is defined as the exponent of the standard deviation of ``log(a)``.
    Mathematically the population geometric standard deviation can be
    evaluated as::

        gstd = exp(std(log(a)))

    .. versionadded:: 1.3.0

    Parameters
    ----------
    a : array_like
        An array like object containing the sample data.
    axis : int, tuple or None, optional
        Axis along which to operate. Default is 0. If None, compute over
        the whole array `a`.
    ddof : int, optional
        Degree of freedom correction in the calculation of the
        geometric standard deviation. Default is 1.

    Returns
    -------
    ndarray or float
        An array of the geometric standard deviation. If `axis` is None or `a`
        is a 1d array a float is returned.

    Notes
    -----
    As the calculation requires the use of logarithms the geometric standard
    deviation only supports strictly positive values. Any non-positive or
    infinite values will raise a `ValueError`.
    The geometric standard deviation is sometimes confused with the exponent of
    the standard deviation, ``exp(std(a))``. Instead the geometric standard
    deviation is ``exp(std(log(a)))``.
    The default value for `ddof` is different to the default value (0) used
    by other ddof containing functions, such as ``np.std`` and ``np.nanstd``.

    Examples
    --------
    Find the geometric standard deviation of a log-normally distributed sample.
    Note that the standard deviation of the distribution is one, on a
    log scale this evaluates to approximately ``exp(1)``.

    >>> from scipy.stats import gstd
    >>> np.random.seed(123)
    >>> sample = np.random.lognormal(mean=0, sigma=1, size=1000)
    >>> gstd(sample)
    2.7217860664589946

    Compute the geometric standard deviation of a multidimensional array and
    of a given axis.

    >>> a = np.arange(1, 25).reshape(2, 3, 4)
    >>> gstd(a, axis=None)
    2.2944076136018947
    >>> gstd(a, axis=2)
    array([[1.82424757, 1.22436866, 1.13183117],
           [1.09348306, 1.07244798, 1.05914985]])
    >>> gstd(a, axis=(1,2))
    array([2.12939215, 1.22120169])

    The geometric standard deviation further handles masked arrays.

    >>> a = np.arange(1, 25).reshape(2, 3, 4)
    >>> ma = np.ma.masked_where(a > 16, a)
    >>> ma
    masked_array(
      data=[[[1, 2, 3, 4],
             [5, 6, 7, 8],
             [9, 10, 11, 12]],
            [[13, 14, 15, 16],
             [--, --, --, --],
             [--, --, --, --]]],
      mask=[[[False, False, False, False],
             [False, False, False, False],
             [False, False, False, False]],
            [[False, False, False, False],
             [ True,  True,  True,  True],
             [ True,  True,  True,  True]]],
      fill_value=999999)
    >>> gstd(ma, axis=2)
    masked_array(
      data=[[1.8242475707663655, 1.2243686572447428, 1.1318311657788478],
            [1.0934830582350938, --, --]],
      mask=[[False, False, False],
            [False,  True,  True]],
      fill_value=999999)
    """
    a = np.asanyarray(a)
    log = ma.log if isinstance(a, ma.MaskedArray) else np.log

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            return np.exp(np.std(log(a), axis=axis, ddof=ddof))
    except RuntimeWarning as w:
        if np.isinf(a).any():
            raise ValueError(
                'Infinite value encountered. The geometric standard deviation '
                'is defined for strictly positive values only.')
        elif np.less_equal(a, 0).any():
            raise ValueError(
                'Non positive value encountered. The geometric standard '
                'deviation is defined for strictly positive values only.')
        elif 'Degrees of freedom <= 0 for slice' == str(w):
            raise ValueError(w)
        else:
            #  Remaining warnings don't need to be exceptions.
            warnings.warn(w)
    except TypeError:
        raise ValueError(
            'Invalid array input. The inputs could not be '
            'safely coerced to any supported types')


# Private dictionary initialized only once at module level
# See https://en.wikipedia.org/wiki/Robust_measures_of_scale
_scale_conversions = {'raw': 1.0,
                      'normal': special.erfinv(0.5) * 2.0 * math.sqrt(2.0)}


def iqr(x, axis=None, rng=(25, 75), scale='raw', nan_policy='propagate',
        interpolation='linear', keepdims=False):
    """
    Compute the interquartile range of the data along the specified axis.

    The interquartile range (IQR) is the difference between the 75th and
    25th percentile of the data. It is a measure of the dispersion
    similar to standard deviation or variance, but is much more robust
    against outliers [2]_.

    The ``rng`` parameter allows this function to compute other
    percentile ranges than the actual IQR. For example, setting
    ``rng=(0, 100)`` is equivalent to `numpy.ptp`.

    The IQR of an empty array is `np.nan`.

    .. versionadded:: 0.18.0

    Parameters
    ----------
    x : array_like
        Input array or object that can be converted to an array.
    axis : int or sequence of int, optional
        Axis along which the range is computed. The default is to
        compute the IQR for the entire array.
    rng : Two-element sequence containing floats in range of [0,100] optional
        Percentiles over which to compute the range. Each must be
        between 0 and 100, inclusive. The default is the true IQR:
        `(25, 75)`. The order of the elements is not important.
    scale : scalar or str, optional
        The numerical value of scale will be divided out of the final
        result. The following string values are recognized:

          'raw' : No scaling, just return the raw IQR.
          'normal' : Scale by :math:`2 \\sqrt{2} erf^{-1}(\\frac{1}{2}) \\approx 1.349`.

        The default is 'raw'. Array-like scale is also allowed, as long
        as it broadcasts correctly to the output such that
        ``out / scale`` is a valid operation. The output dimensions
        depend on the input array, `x`, the `axis` argument, and the
        `keepdims` flag.
    nan_policy : {'propagate', 'raise', 'omit'}, optional
        Defines how to handle when input contains nan. 'propagate'
        returns nan, 'raise' throws an error, 'omit' performs the
        calculations ignoring nan values. Default is 'propagate'.
    interpolation : {'linear', 'lower', 'higher', 'midpoint', 'nearest'}, optional
        Specifies the interpolation method to use when the percentile
        boundaries lie between two data points `i` and `j`:

          * 'linear' : `i + (j - i) * fraction`, where `fraction` is the
              fractional part of the index surrounded by `i` and `j`.
          * 'lower' : `i`.
          * 'higher' : `j`.
          * 'nearest' : `i` or `j` whichever is nearest.
          * 'midpoint' : `(i + j) / 2`.

        Default is 'linear'.
    keepdims : bool, optional
        If this is set to `True`, the reduced axes are left in the
        result as dimensions with size one. With this option, the result
        will broadcast correctly against the original array `x`.

    Returns
    -------
    iqr : scalar or ndarray
        If ``axis=None``, a scalar is returned. If the input contains
        integers or floats of smaller precision than ``np.float64``, then the
        output data-type is ``np.float64``. Otherwise, the output data-type is
        the same as that of the input.

    See Also
    --------
    numpy.std, numpy.var

    Examples
    --------
    >>> from scipy.stats import iqr
    >>> x = np.array([[10, 7, 4], [3, 2, 1]])
    >>> x
    array([[10,  7,  4],
           [ 3,  2,  1]])
    >>> iqr(x)
    4.0
    >>> iqr(x, axis=0)
    array([ 3.5,  2.5,  1.5])
    >>> iqr(x, axis=1)
    array([ 3.,  1.])
    >>> iqr(x, axis=1, keepdims=True)
    array([[ 3.],
           [ 1.]])

    Notes
    -----
    This function is heavily dependent on the version of `numpy` that is
    installed. Versions greater than 1.11.0b3 are highly recommended, as they
    include a number of enhancements and fixes to `numpy.percentile` and
    `numpy.nanpercentile` that affect the operation of this function. The
    following modifications apply:

    Below 1.10.0 : `nan_policy` is poorly defined.
        The default behavior of `numpy.percentile` is used for 'propagate'. This
        is a hybrid of 'omit' and 'propagate' that mostly yields a skewed
        version of 'omit' since NaNs are sorted to the end of the data. A
        warning is raised if there are NaNs in the data.
    Below 1.9.0: `numpy.nanpercentile` does not exist.
        This means that `numpy.percentile` is used regardless of `nan_policy`
        and a warning is issued. See previous item for a description of the
        behavior.
    Below 1.9.0: `keepdims` and `interpolation` are not supported.
        The keywords get ignored with a warning if supplied with non-default
        values. However, multiple axes are still supported.

    References
    ----------
    .. [1] "Interquartile range" https://en.wikipedia.org/wiki/Interquartile_range
    .. [2] "Robust measures of scale" https://en.wikipedia.org/wiki/Robust_measures_of_scale
    .. [3] "Quantile" https://en.wikipedia.org/wiki/Quantile
    """
    x = asarray(x)

    # This check prevents percentile from raising an error later. Also, it is
    # consistent with `np.var` and `np.std`.
    if not x.size:
        return np.nan

    # An error may be raised here, so fail-fast, before doing lengthy
    # computations, even though `scale` is not used until later
    if isinstance(scale, string_types):
        scale_key = scale.lower()
        if scale_key not in _scale_conversions:
            raise ValueError("{0} not a valid scale for `iqr`".format(scale))
        scale = _scale_conversions[scale_key]

    # Select the percentile function to use based on nans and policy
    contains_nan, nan_policy = _contains_nan(x, nan_policy)

    if contains_nan and nan_policy == 'omit':
        percentile_func = _iqr_nanpercentile
    else:
        percentile_func = _iqr_percentile

    if len(rng) != 2:
        raise TypeError("quantile range must be two element sequence")

    if np.isnan(rng).any():
        raise ValueError("range must not contain NaNs")

    rng = sorted(rng)
    pct = percentile_func(x, rng, axis=axis, interpolation=interpolation,
                          keepdims=keepdims, contains_nan=contains_nan)
    out = np.subtract(pct[1], pct[0])

    if scale != 1.0:
        out /= scale

    return out


def median_absolute_deviation(x, axis=0, center=np.median, scale=1.4826,
                              nan_policy='propagate'):
    """
    Compute the median absolute deviation of the data along the given axis.

    The median absolute deviation (MAD, [1]_) computes the median over the
    absolute deviations from the median. It is a measure of dispersion
    similar to the standard deviation, but is more robust to outliers [2]_.

    The MAD of an empty array is ``np.nan``.

    .. versionadded:: 1.3.0

    Parameters
    ----------
    x : array_like
        Input array or object that can be converted to an array.
    axis : int or None, optional
        Axis along which the range is computed. Default is 0. If None, compute
        the MAD over the entire array.
    center : callable, optional
        A function that will return the central value. The default is to use
        np.median. Any user defined function used will need to have the function
        signature ``func(arr, axis)``.
    scale : int, optional
        The scaling factor applied to the MAD. The default scale (1.4826)
        ensures consistency with the standard deviation for normally distributed
        data.
    nan_policy : {'propagate', 'raise', 'omit'}, optional
        Defines how to handle when input contains nan. 'propagate'
        returns nan, 'raise' throws an error, 'omit' performs the
        calculations ignoring nan values. Default is 'propagate'.

    Returns
    -------
    mad : scalar or ndarray
        If ``axis=None``, a scalar is returned. If the input contains
        integers or floats of smaller precision than ``np.float64``, then the
        output data-type is ``np.float64``. Otherwise, the output data-type is
        the same as that of the input.

    See Also
    --------
    numpy.std, numpy.var, numpy.median, scipy.stats.iqr, scipy.stats.tmean,
    scipy.stats.tstd, scipy.stats.tvar

    Notes
    -----
    The `center` argument only affects the calculation of the central value
    around which the MAD is calculated. That is, passing in ``center=np.mean``
    will calculate the MAD around the mean - it will not calculate the *mean*
    absolute deviation.

    References
    ----------
    .. [1] "Median absolute deviation" https://en.wikipedia.org/wiki/Median_absolute_deviation
    .. [2] "Robust measures of scale" https://en.wikipedia.org/wiki/Robust_measures_of_scale

    Examples
    --------
    When comparing the behavior of `median_absolute_deviation` with ``np.std``,
    the latter is affected when we change a single value of an array to have an
    outlier value while the MAD hardly changes:

    >>> from scipy import stats
    >>> x = stats.norm.rvs(size=100, scale=1, random_state=123456)
    >>> x.std()
    0.9973906394005013
    >>> stats.median_absolute_deviation(x)
    1.2280762773108278
    >>> x[0] = 345.6
    >>> x.std()
    34.42304872314415
    >>> stats.median_absolute_deviation(x)
    1.2340335571164334

    Axis handling example:

    >>> x = np.array([[10, 7, 4], [3, 2, 1]])
    >>> x
    array([[10,  7,  4],
           [ 3,  2,  1]])
    >>> stats.median_absolute_deviation(x)
    array([5.1891, 3.7065, 2.2239])
    >>> stats.median_absolute_deviation(x, axis=None)
    2.9652

    """
    x = asarray(x)

    # Consistent with `np.var` and `np.std`.
    if not x.size:
        return np.nan

    contains_nan, nan_policy = _contains_nan(x, nan_policy)

    if contains_nan and nan_policy == 'propagate':
        return np.nan

    if contains_nan and nan_policy == 'omit':
        # Way faster than carrying the masks around
        arr = ma.masked_invalid(x).compressed()
    else:
        arr = x

    if axis is None:
        med = center(arr)
        mad = np.median(np.abs(arr - med))
    else:
        med = np.apply_over_axes(center, arr, axis)
        mad = np.median(np.abs(arr - med), axis=axis)

    return scale * mad


def _iqr_percentile(x, q, axis=None, interpolation='linear', keepdims=False, contains_nan=False):
    """
    Private wrapper that works around older versions of `numpy`.

    While this function is pretty much necessary for the moment, it
    should be removed as soon as the minimum supported numpy version
    allows.
    """
    if contains_nan and NumpyVersion(np.__version__) < '1.10.0a':
        # I see no way to avoid the version check to ensure that the corrected
        # NaN behavior has been implemented except to call `percentile` on a
        # small array.
        msg = "Keyword nan_policy='propagate' not correctly supported for " \
              "numpy versions < 1.10.x. The default behavior of " \
              "`numpy.percentile` will be used."
        warnings.warn(msg, RuntimeWarning)

    try:
        # For older versions of numpy, there are two things that can cause a
        # problem here: missing keywords and non-scalar axis. The former can be
        # partially handled with a warning, the latter can be handled fully by
        # hacking in an implementation similar to numpy's function for
        # providing multi-axis functionality
        # (`numpy.lib.function_base._ureduce` for the curious).
        result = np.percentile(x, q, axis=axis, keepdims=keepdims,
                               interpolation=interpolation)
    except TypeError:
        if interpolation != 'linear' or keepdims:
            # At time or writing, this means np.__version__ < 1.9.0
            warnings.warn("Keywords interpolation and keepdims not supported "
                          "for your version of numpy", RuntimeWarning)
        try:
            # Special processing if axis is an iterable
            original_size = len(axis)
        except TypeError:
            # Axis is a scalar at this point
            pass
        else:
            axis = np.unique(np.asarray(axis) % x.ndim)
            if original_size > axis.size:
                # mimic numpy if axes are duplicated
                raise ValueError("duplicate value in axis")
            if axis.size == x.ndim:
                # axis includes all axes: revert to None
                axis = None
            elif axis.size == 1:
                # no rolling necessary
                axis = axis[0]
            else:
                # roll multiple axes to the end and flatten that part out
                for ax in axis[::-1]:
                    x = np.rollaxis(x, ax, x.ndim)
                x = x.reshape(x.shape[:-axis.size] +
                              (np.prod(x.shape[-axis.size:]),))
                axis = -1
        result = np.percentile(x, q, axis=axis)

    return result


def _iqr_nanpercentile(x, q, axis=None, interpolation='linear', keepdims=False,
                       contains_nan=False):
    """
    Private wrapper that works around the following:

      1. A bug in `np.nanpercentile` that was around until numpy version
         1.11.0.
      2. A bug in `np.percentile` NaN handling that was fixed in numpy
         version 1.10.0.
      3. The non-existence of `np.nanpercentile` before numpy version
         1.9.0.

    While this function is pretty much necessary for the moment, it
    should be removed as soon as the minimum supported numpy version
    allows.
    """
    if hasattr(np, 'nanpercentile'):
        # At time or writing, this means np.__version__ < 1.9.0
        result = np.nanpercentile(x, q, axis=axis,
                                  interpolation=interpolation,
                                  keepdims=keepdims)
        # If non-scalar result and nanpercentile does not do proper axis roll.
        # I see no way of avoiding the version test since dimensions may just
        # happen to match in the data.
        if result.ndim > 1 and NumpyVersion(np.__version__) < '1.11.0a':
            axis = np.asarray(axis)
            if axis.size == 1:
                # If only one axis specified, reduction happens along that dimension
                if axis.ndim == 0:
                    axis = axis[None]
                result = np.rollaxis(result, axis[0])
            else:
                # If multiple axes, reduced dimeision is last
                result = np.rollaxis(result, -1)
    else:
        msg = "Keyword nan_policy='omit' not correctly supported for numpy " \
              "versions < 1.9.x. The default behavior of  numpy.percentile " \
              "will be used."
        warnings.warn(msg, RuntimeWarning)
        result = _iqr_percentile(x, q, axis=axis)

    return result


#####################################
#         TRIMMING FUNCTIONS        #
#####################################

SigmaclipResult = namedtuple('SigmaclipResult', ('clipped', 'lower', 'upper'))


def sigmaclip(a, low=4., high=4.):
    """
    Iterative sigma-clipping of array elements.

    Starting from the full sample, all elements outside the critical range are
    removed, i.e. all elements of the input array `c` that satisfy either of
    the following conditions ::

        c < mean(c) - std(c)*low
        c > mean(c) + std(c)*high

    The iteration continues with the updated sample until no
    elements are outside the (updated) range.

    Parameters
    ----------
    a : array_like
        Data array, will be raveled if not 1-D.
    low : float, optional
        Lower bound factor of sigma clipping. Default is 4.
    high : float, optional
        Upper bound factor of sigma clipping. Default is 4.

    Returns
    -------
    clipped : ndarray
        Input array with clipped elements removed.
    lower : float
        Lower threshold value use for clipping.
    upper : float
        Upper threshold value use for clipping.

    Examples
    --------
    >>> from scipy.stats import sigmaclip
    >>> a = np.concatenate((np.linspace(9.5, 10.5, 31),
    ...                     np.linspace(0, 20, 5)))
    >>> fact = 1.5
    >>> c, low, upp = sigmaclip(a, fact, fact)
    >>> c
    array([  9.96666667,  10.        ,  10.03333333,  10.        ])
    >>> c.var(), c.std()
    (0.00055555555555555165, 0.023570226039551501)
    >>> low, c.mean() - fact*c.std(), c.min()
    (9.9646446609406727, 9.9646446609406727, 9.9666666666666668)
    >>> upp, c.mean() + fact*c.std(), c.max()
    (10.035355339059327, 10.035355339059327, 10.033333333333333)

    >>> a = np.concatenate((np.linspace(9.5, 10.5, 11),
    ...                     np.linspace(-100, -50, 3)))
    >>> c, low, upp = sigmaclip(a, 1.8, 1.8)
    >>> (c == np.linspace(9.5, 10.5, 11)).all()
    True

    """
    c = np.asarray(a).ravel()
    delta = 1
    while delta:
        c_std = c.std()
        c_mean = c.mean()
        size = c.size
        critlower = c_mean - c_std * low
        critupper = c_mean + c_std * high
        c = c[(c >= critlower) & (c <= critupper)]
        delta = size - c.size

    return SigmaclipResult(c, critlower, critupper)


def trimboth(a, proportiontocut, axis=0):
    """
    Slices off a proportion of items from both ends of an array.

    Slices off the passed proportion of items from both ends of the passed
    array (i.e., with `proportiontocut` = 0.1, slices leftmost 10% **and**
    rightmost 10% of scores). The trimmed values are the lowest and
    highest ones.
    Slices off less if proportion results in a non-integer slice index (i.e.,
    conservatively slices off`proportiontocut`).

    Parameters
    ----------
    a : array_like
        Data to trim.
    proportiontocut : float
        Proportion (in range 0-1) of total data set to trim of each end.
    axis : int or None, optional
        Axis along which to trim data. Default is 0. If None, compute over
        the whole array `a`.

    Returns
    -------
    out : ndarray
        Trimmed version of array `a`. The order of the trimmed content
        is undefined.

    See Also
    --------
    trim_mean

    Examples
    --------
    >>> from scipy import stats
    >>> a = np.arange(20)
    >>> b = stats.trimboth(a, 0.1)
    >>> b.shape
    (16,)

    """
    a = np.asarray(a)

    if a.size == 0:
        return a

    if axis is None:
        a = a.ravel()
        axis = 0

    nobs = a.shape[axis]
    lowercut = int(proportiontocut * nobs)
    uppercut = nobs - lowercut
    if (lowercut >= uppercut):
        raise ValueError("Proportion too big.")

    atmp = np.partition(a, (lowercut, uppercut - 1), axis)

    sl = [slice(None)] * atmp.ndim
    sl[axis] = slice(lowercut, uppercut)
    return atmp[tuple(sl)]


def trim1(a, proportiontocut, tail='right', axis=0):
    """
    Slices off a proportion from ONE end of the passed array distribution.

    If `proportiontocut` = 0.1, slices off 'leftmost' or 'rightmost'
    10% of scores. The lowest or highest values are trimmed (depending on
    the tail).
    Slices off less if proportion results in a non-integer slice index
    (i.e., conservatively slices off `proportiontocut` ).

    Parameters
    ----------
    a : array_like
        Input array
    proportiontocut : float
        Fraction to cut off of 'left' or 'right' of distribution
    tail : {'left', 'right'}, optional
        Defaults to 'right'.
    axis : int or None, optional
        Axis along which to trim data. Default is 0. If None, compute over
        the whole array `a`.

    Returns
    -------
    trim1 : ndarray
        Trimmed version of array `a`. The order of the trimmed content is
        undefined.

    """
    a = np.asarray(a)
    if axis is None:
        a = a.ravel()
        axis = 0

    nobs = a.shape[axis]

    # avoid possible corner case
    if proportiontocut >= 1:
        return []

    if tail.lower() == 'right':
        lowercut = 0
        uppercut = nobs - int(proportiontocut * nobs)

    elif tail.lower() == 'left':
        lowercut = int(proportiontocut * nobs)
        uppercut = nobs

    atmp = np.partition(a, (lowercut, uppercut - 1), axis)

    return atmp[lowercut:uppercut]


def trim_mean(a, proportiontocut, axis=0):
    """
    Return mean of array after trimming distribution from both tails.

    If `proportiontocut` = 0.1, slices off 'leftmost' and 'rightmost' 10% of
    scores. The input is sorted before slicing. Slices off less if proportion
    results in a non-integer slice index (i.e., conservatively slices off
    `proportiontocut` ).

    Parameters
    ----------
    a : array_like
        Input array
    proportiontocut : float
        Fraction to cut off of both tails of the distribution
    axis : int or None, optional
        Axis along which the trimmed means are computed. Default is 0.
        If None, compute over the whole array `a`.

    Returns
    -------
    trim_mean : ndarray
        Mean of trimmed array.

    See Also
    --------
    trimboth
    tmean : compute the trimmed mean ignoring values outside given `limits`.

    Examples
    --------
    >>> from scipy import stats
    >>> x = np.arange(20)
    >>> stats.trim_mean(x, 0.1)
    9.5
    >>> x2 = x.reshape(5, 4)
    >>> x2
    array([[ 0,  1,  2,  3],
           [ 4,  5,  6,  7],
           [ 8,  9, 10, 11],
           [12, 13, 14, 15],
           [16, 17, 18, 19]])
    >>> stats.trim_mean(x2, 0.25)
    array([  8.,   9.,  10.,  11.])
    >>> stats.trim_mean(x2, 0.25, axis=1)
    array([  1.5,   5.5,   9.5,  13.5,  17.5])

    """
    a = np.asarray(a)

    if a.size == 0:
        return np.nan

    if axis is None:
        a = a.ravel()
        axis = 0

    nobs = a.shape[axis]
    lowercut = int(proportiontocut * nobs)
    uppercut = nobs - lowercut
    if (lowercut > uppercut):
        raise ValueError("Proportion too big.")

    atmp = np.partition(a, (lowercut, uppercut - 1), axis)

    sl = [slice(None)] * atmp.ndim
    sl[axis] = slice(lowercut, uppercut)
    return np.mean(atmp[tuple(sl)], axis=axis)


F_onewayResult = namedtuple('F_onewayResult', ('statistic', 'pvalue'))


def f_oneway(*args):
    """
    Performs a 1-way ANOVA.

    The one-way ANOVA tests the null hypothesis that two or more groups have
    the same population mean.  The test is applied to samples from two or
    more groups, possibly with differing sizes.

    Parameters
    ----------
    sample1, sample2, ... : array_like
        The sample measurements for each group.

    Returns
    -------
    statistic : float
        The computed F-value of the test.
    pvalue : float
        The associated p-value from the F-distribution.

    Notes
    -----
    The ANOVA test has important assumptions that must be satisfied in order
    for the associated p-value to be valid.

    1. The samples are independent.
    2. Each sample is from a normally distributed population.
    3. The population standard deviations of the groups are all equal.  This
       property is known as homoscedasticity.

    If these assumptions are not true for a given set of data, it may still be
    possible to use the Kruskal-Wallis H-test (`scipy.stats.kruskal`) although
    with some loss of power.

    The algorithm is from Heiman[2], pp.394-7.


    References
    ----------
    .. [1] R. Lowry, "Concepts and Applications of Inferential Statistics",
           Chapter 14, 2014, http://vassarstats.net/textbook/

    .. [2] G.W. Heiman, "Understanding research methods and statistics: An
           integrated introduction for psychology", Houghton, Mifflin and
           Company, 2001.

    .. [3] G.H. McDonald, "Handbook of Biological Statistics", One-way ANOVA.
           http://www.biostathandbook.com/onewayanova.html

    Examples
    --------
    >>> import scipy.stats as stats

    [3]_ Here are some data on a shell measurement (the length of the anterior
    adductor muscle scar, standardized by dividing by length) in the mussel
    Mytilus trossulus from five locations: Tillamook, Oregon; Newport, Oregon;
    Petersburg, Alaska; Magadan, Russia; and Tvarminne, Finland, taken from a
    much larger data set used in McDonald et al. (1991).

    >>> tillamook = [0.0571, 0.0813, 0.0831, 0.0976, 0.0817, 0.0859, 0.0735,
    ...              0.0659, 0.0923, 0.0836]
    >>> newport = [0.0873, 0.0662, 0.0672, 0.0819, 0.0749, 0.0649, 0.0835,
    ...            0.0725]
    >>> petersburg = [0.0974, 0.1352, 0.0817, 0.1016, 0.0968, 0.1064, 0.105]
    >>> magadan = [0.1033, 0.0915, 0.0781, 0.0685, 0.0677, 0.0697, 0.0764,
    ...            0.0689]
    >>> tvarminne = [0.0703, 0.1026, 0.0956, 0.0973, 0.1039, 0.1045]
    >>> stats.f_oneway(tillamook, newport, petersburg, magadan, tvarminne)
    (7.1210194716424473, 0.00028122423145345439)

    """
    args = [np.asarray(arg, dtype=float) for arg in args]
    # ANOVA on N groups, each in its own array
    num_groups = len(args)
    alldata = np.concatenate(args)
    bign = len(alldata)

    # Determine the mean of the data, and subtract that from all inputs to a
    # variance (via sum_of_sq / sq_of_sum) calculation.  Variance is invariance
    # to a shift in location, and centering all data around zero vastly
    # improves numerical stability.
    offset = alldata.mean()
    alldata -= offset

    sstot = _sum_of_squares(alldata) - (_square_of_sums(alldata) / bign)
    ssbn = 0
    for a in args:
        ssbn += _square_of_sums(a - offset) / len(a)

    # Naming: variables ending in bn/b are for "between treatments", wn/w are
    # for "within treatments"
    ssbn -= _square_of_sums(alldata) / bign
    sswn = sstot - ssbn
    dfbn = num_groups - 1
    dfwn = bign - num_groups
    msb = ssbn / dfbn
    msw = sswn / dfwn
    f = msb / msw

    prob = special.fdtrc(dfbn, dfwn, f)   # equivalent to stats.f.sf

    return F_onewayResult(f, prob)


class PearsonRConstantInputWarning(RuntimeWarning):
    """
    Warning generated by `pearsonr` when an input is constant.
    """

    def __init__(self, msg=None):
        if msg is None:
            msg = ("An input array is constant; the correlation coefficent "
                   "is not defined.")
        self.args = (msg,)


class PearsonRNearConstantInputWarning(RuntimeWarning):
    """
    Warning generated by `pearsonr` when an input is nearly constant.
    """

    def __init__(self, msg=None):
        if msg is None:
            msg = ("An input array is nearly constant; the computed "
                   "correlation coefficent may be inaccurate.")
        self.args = (msg,)


def pearsonr(x, y):
    r"""
    Pearson correlation coefficient and p-value for testing non-correlation.

    The Pearson correlation coefficient [1]_ measures the linear relationship
    between two datasets.  The calculation of the p-value relies on the
    assumption that each dataset is normally distributed.  (See Kowalski [3]_
    for a discussion of the effects of non-normality of the input on the
    distribution of the correlation coefficient.)  Like other correlation
    coefficients, this one varies between -1 and +1 with 0 implying no
    correlation. Correlations of -1 or +1 imply an exact linear relationship.
    Positive correlations imply that as x increases, so does y. Negative
    correlations imply that as x increases, y decreases.

    The p-value roughly indicates the probability of an uncorrelated system
    producing datasets that have a Pearson correlation at least as extreme
    as the one computed from these datasets.

    Parameters
    ----------
    x : (N,) array_like
        Input
    y : (N,) array_like
        Input

    Returns
    -------
    r : float
        Pearson's correlation coefficient
    p-value : float
        two-tailed p-value

    Warns
    -----
    PearsonRConstantInputWarning
        Raised if an input is a constant array.  The correlation coefficient
        is not defined in this case, so ``np.nan`` is returned.

    PearsonRNearConstantInputWarning
        Raised if an input is "nearly" constant.  The array ``x`` is considered
        nearly constant if ``norm(x - mean(x)) < 1e-13 * abs(mean(x))``.
        Numerical errors in the calculation ``x - mean(x)`` in this case might
        result in an inaccurate calculation of r.

    See Also
    --------
    spearmanr : Spearman rank-order correlation coefficient.
    kendalltau : Kendall's tau, a correlation measure for ordinal data.

    Notes
    -----

    The correlation coefficient is calculated as follows:

    .. math::

        r = \frac{\sum (x - m_x) (y - m_y)}
                 {\sqrt{\sum (x - m_x)^2 \sum (y - m_y)^2}}

    where :math:`m_x` is the mean of the vector :math:`x` and :math:`m_y` is
    the mean of the vector :math:`y`.

    Under the assumption that x and y are drawn from independent normal
    distributions (so the population correlation coefficient is 0), the
    probability density function of the sample correlation coefficient r
    is ([1]_, [2]_)::

               (1 - r**2)**(n/2 - 2)
        f(r) = ---------------------
                  B(1/2, n/2 - 1)

    where n is the number of samples, and B is the beta function.  This
    is sometimes referred to as the exact distribution of r.  This is
    the distribution that is used in `pearsonr` to compute the p-value.
    The distribution is a beta distribution on the interval [-1, 1],
    with equal shape parameters a = b = n/2 - 1.  In terms of SciPy's
    implementation of the beta distribution, the distribution of r is::

        dist = scipy.stats.beta(n/2 - 1, n/2 - 1, loc=-1, scale=2)

    The p-value returned by `pearsonr` is a two-sided p-value.  For a
    given sample with correlation coefficient r, the p-value is
    the probability that abs(r') of a random sample x' and y' drawn from
    the population with zero correlation would be greater than or equal
    to abs(r).  In terms of the object `dist` shown above, the p-value
    for a given r and length n can be computed as::

        p = 2*dist.cdf(-abs(r))

    When n is 2, the above continuous distribution is not well-defined.
    One can interpret the limit of the beta distribution as the shape
    parameters a and b approach a = b = 0 as a discrete distribution with
    equal probability masses at r = 1 and r = -1.  More directly, one
    can observe that, given the data x = [x1, x2] and y = [y1, y2], and
    assuming x1 != x2 and y1 != y2, the only possible values for r are 1
    and -1.  Because abs(r') for any sample x' and y' with length 2 will
    be 1, the two-sided p-value for a sample of length 2 is always 1.

    References
    ----------
    .. [1] "Pearson correlation coefficient", Wikipedia,
           https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
    .. [2] Student, "Probable error of a correlation coefficient",
           Biometrika, Volume 6, Issue 2-3, 1 September 1908, pp. 302-310.
    .. [3] C. J. Kowalski, "On the Effects of Non-Normality on the Distribution
           of the Sample Product-Moment Correlation Coefficient"
           Journal of the Royal Statistical Society. Series C (Applied
           Statistics), Vol. 21, No. 1 (1972), pp. 1-12.

    Examples
    --------
    >>> from scipy import stats
    >>> a = np.array([0, 0, 0, 1, 1, 1, 1])
    >>> b = np.arange(7)
    >>> stats.pearsonr(a, b)
    (0.8660254037844386, 0.011724811003954649)

    >>> stats.pearsonr([1, 2, 3, 4, 5], [10, 9, 2.5, 6, 4])
    (-0.7426106572325057, 0.1505558088534455)

    """
    n = len(x)
    if n != len(y):
        raise ValueError('x and y must have the same length.')

    if n < 2:
        raise ValueError('x and y must have length at least 2.')

    x = np.asarray(x)
    y = np.asarray(y)

    # If an input is constant, the correlation coefficient is not defined.
    if (x == x[0]).all() or (y == y[0]).all():
        warnings.warn(PearsonRConstantInputWarning())
        return np.nan, np.nan

    # dtype is the data type for the calculations.  This expression ensures
    # that the data type is at least 64 bit floating point.  It might have
    # more precision if the input is, for example, np.longdouble.
    dtype = type(1.0 + x[0] + y[0])

    if n == 2:
        return dtype(np.sign(x[1] - x[0])*np.sign(y[1] - y[0])), 1.0

    xmean = x.mean(dtype=dtype)
    ymean = y.mean(dtype=dtype)

    # By using `astype(dtype)`, we ensure that the intermediate calculations
    # use at least 64 bit floating point.
    xm = x.astype(dtype) - xmean
    ym = y.astype(dtype) - ymean

    # Unlike np.linalg.norm or the expression sqrt((xm*xm).sum()),
    # scipy.linalg.norm(xm) does not overflow if xm is, for example,
    # [-5e210, 5e210, 3e200, -3e200]
    normxm = linalg.norm(xm)
    normym = linalg.norm(ym)

    threshold = 1e-13
    if normxm < threshold*abs(xmean) or normym < threshold*abs(ymean):
        # If all the values in x (likewise y) are very close to the mean,
        # the loss of precision that occurs in the subtraction xm = x - xmean
        # might result in large errors in r.
        warnings.warn(PearsonRNearConstantInputWarning())

    r = np.dot(xm/normxm, ym/normym)

    # Presumably, if abs(r) > 1, then it is only some small artifact of
    # floating point arithmetic.
    r = max(min(r, 1.0), -1.0)

    # As explained in the docstring, the p-value can be computed as
    #     p = 2*dist.cdf(-abs(r))
    # where dist is the beta distribution on [-1, 1] with shape parameters
    # a = b = n/2 - 1.  `special.btdtr` is the CDF for the beta distribution
    # on [0, 1].  To use it, we make the transformation  x = (r + 1)/2; the
    # shape parameters do not change.  Then -abs(r) used in `cdf(-abs(r))`
    # becomes x = (-abs(r) + 1)/2 = 0.5*(1 - abs(r)).  (r is cast to float64
    # to avoid a TypeError raised by btdtr when r is higher precision.)
    ab = n/2 - 1
    prob = 2*special.btdtr(ab, ab, 0.5*(1 - abs(np.float64(r))))

    return r, prob


def fisher_exact(table, alternative='two-sided'):
    """Performs a Fisher exact test on a 2x2 contingency table.

    Parameters
    ----------
    table : array_like of ints
        A 2x2 contingency table.  Elements should be non-negative integers.
    alternative : {'two-sided', 'less', 'greater'}, optional
        Which alternative hypothesis to the null hypothesis the test uses.
        Default is 'two-sided'.

    Returns
    -------
    oddsratio : float
        This is prior odds ratio and not a posterior estimate.
    p_value : float
        P-value, the probability of obtaining a distribution at least as
        extreme as the one that was actually observed, assuming that the
        null hypothesis is true.

    See Also
    --------
    chi2_contingency : Chi-square test of independence of variables in a
        contingency table.

    Notes
    -----
    The calculated odds ratio is different from the one R uses. This scipy
    implementation returns the (more common) "unconditional Maximum
    Likelihood Estimate", while R uses the "conditional Maximum Likelihood
    Estimate".

    For tables with large numbers, the (inexact) chi-square test implemented
    in the function `chi2_contingency` can also be used.

    Examples
    --------
    Say we spend a few days counting whales and sharks in the Atlantic and
    Indian oceans. In the Atlantic ocean we find 8 whales and 1 shark, in the
    Indian ocean 2 whales and 5 sharks. Then our contingency table is::

                Atlantic  Indian
        whales     8        2
        sharks     1        5

    We use this table to find the p-value:

    >>> import scipy.stats as stats
    >>> oddsratio, pvalue = stats.fisher_exact([[8, 2], [1, 5]])
    >>> pvalue
    0.0349...

    The probability that we would observe this or an even more imbalanced ratio
    by chance is about 3.5%.  A commonly used significance level is 5%--if we
    adopt that, we can therefore conclude that our observed imbalance is
    statistically significant; whales prefer the Atlantic while sharks prefer
    the Indian ocean.

    """
    hypergeom = distributions.hypergeom
    c = np.asarray(table, dtype=np.int64)  # int32 is not enough for the algorithm
    if not c.shape == (2, 2):
        raise ValueError("The input `table` must be of shape (2, 2).")

    if np.any(c < 0):
        raise ValueError("All values in `table` must be nonnegative.")

    if 0 in c.sum(axis=0) or 0 in c.sum(axis=1):
        # If both values in a row or column are zero, the p-value is 1 and
        # the odds ratio is NaN.
        return np.nan, 1.0

    if c[1, 0] > 0 and c[0, 1] > 0:
        oddsratio = c[0, 0] * c[1, 1] / (c[1, 0] * c[0, 1])
    else:
        oddsratio = np.inf

    n1 = c[0, 0] + c[0, 1]
    n2 = c[1, 0] + c[1, 1]
    n = c[0, 0] + c[1, 0]

    def binary_search(n, n1, n2, side):
        """Binary search for where to begin lower/upper halves in two-sided
        test.
        """
        if side == "upper":
            minval = mode
            maxval = n
        else:
            minval = 0
            maxval = mode
        guess = -1
        while maxval - minval > 1:
            if maxval == minval + 1 and guess == minval:
                guess = maxval
            else:
                guess = (maxval + minval) // 2
            pguess = hypergeom.pmf(guess, n1 + n2, n1, n)
            if side == "upper":
                ng = guess - 1
            else:
                ng = guess + 1
            if pguess <= pexact < hypergeom.pmf(ng, n1 + n2, n1, n):
                break
            elif pguess < pexact:
                maxval = guess
            else:
                minval = guess
        if guess == -1:
            guess = minval
        if side == "upper":
            while guess > 0 and hypergeom.pmf(guess, n1 + n2, n1, n) < pexact * epsilon:
                guess -= 1
            while hypergeom.pmf(guess, n1 + n2, n1, n) > pexact / epsilon:
                guess += 1
        else:
            while hypergeom.pmf(guess, n1 + n2, n1, n) < pexact * epsilon:
                guess += 1
            while guess > 0 and hypergeom.pmf(guess, n1 + n2, n1, n) > pexact / epsilon:
                guess -= 1
        return guess

    if alternative == 'less':
        pvalue = hypergeom.cdf(c[0, 0], n1 + n2, n1, n)
    elif alternative == 'greater':
        # Same formula as the 'less' case, but with the second column.
        pvalue = hypergeom.cdf(c[0, 1], n1 + n2, n1, c[0, 1] + c[1, 1])
    elif alternative == 'two-sided':
        mode = int((n + 1) * (n1 + 1) / (n1 + n2 + 2))
        pexact = hypergeom.pmf(c[0, 0], n1 + n2, n1, n)
        pmode = hypergeom.pmf(mode, n1 + n2, n1, n)

        epsilon = 1 - 1e-4
        if np.abs(pexact - pmode) / np.maximum(pexact, pmode) <= 1 - epsilon:
            return oddsratio, 1.

        elif c[0, 0] < mode:
            plower = hypergeom.cdf(c[0, 0], n1 + n2, n1, n)
            if hypergeom.pmf(n, n1 + n2, n1, n) > pexact / epsilon:
                return oddsratio, plower

            guess = binary_search(n, n1, n2, "upper")
            pvalue = plower + hypergeom.sf(guess - 1, n1 + n2, n1, n)
        else:
            pupper = hypergeom.sf(c[0, 0] - 1, n1 + n2, n1, n)
            if hypergeom.pmf(0, n1 + n2, n1, n) > pexact / epsilon:
                return oddsratio, pupper

            guess = binary_search(n, n1, n2, "lower")
            pvalue = pupper + hypergeom.cdf(guess, n1 + n2, n1, n)
    else:
        msg = "`alternative` should be one of {'two-sided', 'less', 'greater'}"
        raise ValueError(msg)

    pvalue = min(pvalue, 1.0)

    return oddsratio, pvalue


SpearmanrResult = namedtuple('SpearmanrResult', ('correlation', 'pvalue'))


def spearmanr(a, b=None, axis=0, nan_policy='propagate'):
    """
    Calculate a Spearman rank-order correlation coefficient and the p-value
    to test for non-correlation.

    The Spearman correlation is a nonparametric measure of the monotonicity
    of the relationship between two datasets. Unlike the Pearson correlation,
    the Spearman correlation does not assume that both datasets are normally
    distributed. Like other correlation coefficients, this one varies
    between -1 and +1 with 0 implying no correlation. Correlations of -1 or
    +1 imply an exact monotonic relationship. Positive correlations imply that
    as x increases, so does y. Negative correlations imply that as x
    increases, y decreases.

    The p-value roughly indicates the probability of an uncorrelated system
    producing datasets that have a Spearman correlation at least as extreme
    as the one computed from these datasets. The p-values are not entirely
    reliable but are probably reasonable for datasets larger than 500 or so.

    Parameters
    ----------
    a, b : 1D or 2D array_like, b is optional
        One or two 1-D or 2-D arrays containing multiple variables and
        observations. When these are 1-D, each represents a vector of
        observations of a single variable. For the behavior in the 2-D case,
        see under ``axis``, below.
        Both arrays need to have the same length in the ``axis`` dimension.
    axis : int or None, optional
        If axis=0 (default), then each column represents a variable, with
        observations in the rows. If axis=1, the relationship is transposed:
        each row represents a variable, while the columns contain observations.
        If axis=None, then both arrays will be raveled.
    nan_policy : {'propagate', 'raise', 'omit'}, optional
        Defines how to handle when input contains nan. 'propagate' returns nan,
        'raise' throws an error, 'omit' performs the calculations ignoring nan
        values. Default is 'propagate'.

    Returns
    -------
    correlation : float or ndarray (2-D square)
        Spearman correlation matrix or correlation coefficient (if only 2
        variables are given as parameters. Correlation matrix is square with
        length equal to total number of variables (columns or rows) in ``a``
        and ``b`` combined.
    pvalue : float
        The two-sided p-value for a hypothesis test whose null hypothesis is
        that two sets of data are uncorrelated, has same dimension as rho.

    References
    ----------

    .. [1] Zwillinger, D. and Kokoska, S. (2000). CRC Standard
       Probability and Statistics Tables and Formulae. Chapman & Hall: New
       York. 2000.
       Section  14.7

    Examples
    --------
    >>> from scipy import stats
    >>> stats.spearmanr([1,2,3,4,5], [5,6,7,8,7])
    (0.82078268166812329, 0.088587005313543798)
    >>> np.random.seed(1234321)
    >>> x2n = np.random.randn(100, 2)
    >>> y2n = np.random.randn(100, 2)
    >>> stats.spearmanr(x2n)
    (0.059969996999699973, 0.55338590803773591)
    >>> stats.spearmanr(x2n[:,0], x2n[:,1])
    (0.059969996999699973, 0.55338590803773591)
    >>> rho, pval = stats.spearmanr(x2n, y2n)
    >>> rho
    array([[ 1.        ,  0.05997   ,  0.18569457,  0.06258626],
           [ 0.05997   ,  1.        ,  0.110003  ,  0.02534653],
           [ 0.18569457,  0.110003  ,  1.        ,  0.03488749],
           [ 0.06258626,  0.02534653,  0.03488749,  1.        ]])
    >>> pval
    array([[ 0.        ,  0.55338591,  0.06435364,  0.53617935],
           [ 0.55338591,  0.        ,  0.27592895,  0.80234077],
           [ 0.06435364,  0.27592895,  0.        ,  0.73039992],
           [ 0.53617935,  0.80234077,  0.73039992,  0.        ]])
    >>> rho, pval = stats.spearmanr(x2n.T, y2n.T, axis=1)
    >>> rho
    array([[ 1.        ,  0.05997   ,  0.18569457,  0.06258626],
           [ 0.05997   ,  1.        ,  0.110003  ,  0.02534653],
           [ 0.18569457,  0.110003  ,  1.        ,  0.03488749],
           [ 0.06258626,  0.02534653,  0.03488749,  1.        ]])
    >>> stats.spearmanr(x2n, y2n, axis=None)
    (0.10816770419260482, 0.1273562188027364)
    >>> stats.spearmanr(x2n.ravel(), y2n.ravel())
    (0.10816770419260482, 0.1273562188027364)

    >>> xint = np.random.randint(10, size=(100, 2))
    >>> stats.spearmanr(xint)
    (0.052760927029710199, 0.60213045837062351)

    """
    a, axisout = _chk_asarray(a, axis)
    if a.ndim > 2:
        raise ValueError("spearmanr only handles 1-D or 2-D arrays")

    if b is None:
        if a.ndim < 2:
            raise ValueError("`spearmanr` needs at least 2 variables to compare")
    else:
        # Concatenate a and b, so that we now only have to handle the case
        # of a 2-D `a`.
        b, _ = _chk_asarray(b, axis)
        if axisout == 0:
            a = np.column_stack((a, b))
        else:
            a = np.row_stack((a, b))

    n_vars = a.shape[1 - axisout]
    n_obs = a.shape[axisout]
    if n_obs <= 1:
        # Handle empty arrays or single observations.
        return SpearmanrResult(np.nan, np.nan)

    a_contains_nan, nan_policy = _contains_nan(a, nan_policy)
    variable_has_nan = np.zeros(n_vars, dtype=bool)
    if a_contains_nan:
        if nan_policy == 'omit':
            return mstats_basic.spearmanr(a, axis=axis, nan_policy=nan_policy)
        elif nan_policy == 'propagate':
            if a.ndim == 1 or n_vars <= 2:
                return SpearmanrResult(np.nan, np.nan)
            else:
                # Keep track of variables with NaNs, set the outputs to NaN
                # only for those variables
                variable_has_nan = np.isnan(a).sum(axis=axisout)

    a_ranked = np.apply_along_axis(rankdata, axisout, a)
    rs = np.corrcoef(a_ranked, rowvar=axisout)
    dof = n_obs - 2  # degrees of freedom

    # rs can have elements equal to 1, so avoid zero division warnings
    olderr = np.seterr(divide='ignore')
    try:
        # clip the small negative values possibly caused by rounding
        # errors before taking the square root
        t = rs * np.sqrt((dof/((rs+1.0)*(1.0-rs))).clip(0))
    finally:
        np.seterr(**olderr)

    prob = 2 * distributions.t.sf(np.abs(t), dof)

    # For backwards compatibility, return scalars when comparing 2 columns
    if rs.shape == (2, 2):
        return SpearmanrResult(rs[1, 0], prob[1, 0])
    else:
        rs[variable_has_nan, :] = np.nan
        rs[:, variable_has_nan] = np.nan
        return SpearmanrResult(rs, prob)


PointbiserialrResult = namedtuple('PointbiserialrResult',
                                  ('correlation', 'pvalue'))


def pointbiserialr(x, y):
    r"""
    Calculate a point biserial correlation coefficient and its p-value.

    The point biserial correlation is used to measure the relationship
    between a binary variable, x, and a continuous variable, y. Like other
    correlation coefficients, this one varies between -1 and +1 with 0
    implying no correlation. Correlations of -1 or +1 imply a determinative
    relationship.

    This function uses a shortcut formula but produces the same result as
    `pearsonr`.

    Parameters
    ----------
    x : array_like of bools
        Input array.
    y : array_like
        Input array.

    Returns
    -------
    correlation : float
        R value
    pvalue : float
        2-tailed p-value

    Notes
    -----
    `pointbiserialr` uses a t-test with ``n-1`` degrees of freedom.
    It is equivalent to `pearsonr.`

    The value of the point-biserial correlation can be calculated from:

    .. math::

        r_{pb} = \frac{\overline{Y_{1}} -
                 \overline{Y_{0}}}{s_{y}}\sqrt{\frac{N_{1} N_{2}}{N (N - 1))}}

    Where :math:`Y_{0}` and :math:`Y_{1}` are means of the metric
    observations coded 0 and 1 respectively; :math:`N_{0}` and :math:`N_{1}`
    are number of observations coded 0 and 1 respectively; :math:`N` is the
    total number of observations and :math:`s_{y}` is the standard
    deviation of all the metric observations.

    A value of :math:`r_{pb}` that is significantly different from zero is
    completely equivalent to a significant difference in means between the two
    groups. Thus, an independent groups t Test with :math:`N-2` degrees of
    freedom may be used to test whether :math:`r_{pb}` is nonzero. The
    relation between the t-statistic for comparing two independent groups and
    :math:`r_{pb}` is given by:

    .. math::

        t = \sqrt{N - 2}\frac{r_{pb}}{\sqrt{1 - r^{2}_{pb}}}

    References
    ----------
    .. [1] J. Lev, "The Point Biserial Coefficient of Correlation", Ann. Math.
           Statist., Vol. 20, no.1, pp. 125-126, 1949.

    .. [2] R.F. Tate, "Correlation Between a Discrete and a Continuous
           Variable. Point-Biserial Correlation.", Ann. Math. Statist., Vol. 25,
           np. 3, pp. 603-607, 1954.

    .. [3] D. Kornbrot "Point Biserial Correlation", In Wiley StatsRef:
           Statistics Reference Online (eds N. Balakrishnan, et al.), 2014.
           https://doi.org/10.1002/9781118445112.stat06227

    Examples
    --------
    >>> from scipy import stats
    >>> a = np.array([0, 0, 0, 1, 1, 1, 1])
    >>> b = np.arange(7)
    >>> stats.pointbiserialr(a, b)
    (0.8660254037844386, 0.011724811003954652)
    >>> stats.pearsonr(a, b)
    (0.86602540378443871, 0.011724811003954626)
    >>> np.corrcoef(a, b)
    array([[ 1.       ,  0.8660254],
           [ 0.8660254,  1.       ]])

    """
    rpb, prob = pearsonr(x, y)
    return PointbiserialrResult(rpb, prob)


KendalltauResult = namedtuple('KendalltauResult', ('correlation', 'pvalue'))


def kendalltau(x, y, initial_lexsort=None, nan_policy='propagate', method='auto'):
    """
    Calculate Kendall's tau, a correlation measure for ordinal data.

    Kendall's tau is a measure of the correspondence between two rankings.
    Values close to 1 indicate strong agreement, values close to -1 indicate
    strong disagreement.  This is the 1945 "tau-b" version of Kendall's
    tau [2]_, which can account for ties and which reduces to the 1938 "tau-a"
    version [1]_ in absence of ties.

    Parameters
    ----------
    x, y : array_like
        Arrays of rankings, of the same shape. If arrays are not 1-D, they will
        be flattened to 1-D.
    initial_lexsort : bool, optional
        Unused (deprecated).
    nan_policy : {'propagate', 'raise', 'omit'}, optional
        Defines how to handle when input contains nan. 'propagate' returns nan,
        'raise' throws an error, 'omit' performs the calculations ignoring nan
        values. Default is 'propagate'. Note that if the input contains nan
        'omit' delegates to mstats_basic.kendalltau(), which has a different
        implementation.
    method: {'auto', 'asymptotic', 'exact'}, optional
        Defines which method is used to calculate the p-value [5]_.
        'asymptotic' uses a normal approximation valid for large samples.
        'exact' computes the exact p-value, but can only be used if no ties
        are present. 'auto' is the default and selects the appropriate
        method based on a trade-off between speed and accuracy.

    Returns
    -------
    correlation : float
       The tau statistic.
    pvalue : float
       The two-sided p-value for a hypothesis test whose null hypothesis is
       an absence of association, tau = 0.

    See also
    --------
    spearmanr : Calculates a Spearman rank-order correlation coefficient.
    theilslopes : Computes the Theil-Sen estimator for a set of points (x, y).
    weightedtau : Computes a weighted version of Kendall's tau.

    Notes
    -----
    The definition of Kendall's tau that is used is [2]_::

      tau = (P - Q) / sqrt((P + Q + T) * (P + Q + U))

    where P is the number of concordant pairs, Q the number of discordant
    pairs, T the number of ties only in `x`, and U the number of ties only in
    `y`.  If a tie occurs for the same pair in both `x` and `y`, it is not
    added to either T or U.

    References
    ----------
    .. [1] Maurice G. Kendall, "A New Measure of Rank Correlation", Biometrika
           Vol. 30, No. 1/2, pp. 81-93, 1938.
    .. [2] Maurice G. Kendall, "The treatment of ties in ranking problems",
           Biometrika Vol. 33, No. 3, pp. 239-251. 1945.
    .. [3] Gottfried E. Noether, "Elements of Nonparametric Statistics", John
           Wiley & Sons, 1967.
    .. [4] Peter M. Fenwick, "A new data structure for cumulative frequency
           tables", Software: Practice and Experience, Vol. 24, No. 3,
           pp. 327-336, 1994.
    .. [5] Maurice G. Kendall, "Rank Correlation Methods" (4th Edition),
           Charles Griffin & Co., 1970.

    Examples
    --------
    >>> from scipy import stats
    >>> x1 = [12, 2, 1, 12, 2]
    >>> x2 = [1, 4, 7, 1, 0]
    >>> tau, p_value = stats.kendalltau(x1, x2)
    >>> tau
    -0.47140452079103173
    >>> p_value
    0.2827454599327748

    """
    x = np.asarray(x).ravel()
    y = np.asarray(y).ravel()

    if x.size != y.size:
        raise ValueError("All inputs to `kendalltau` must be of the same size, "
                         "found x-size %s and y-size %s" % (x.size, y.size))
    elif not x.size or not y.size:
        return KendalltauResult(np.nan, np.nan)  # Return NaN if arrays are empty

    # check both x and y
    cnx, npx = _contains_nan(x, nan_policy)
    cny, npy = _contains_nan(y, nan_policy)
    contains_nan = cnx or cny
    if npx == 'omit' or npy == 'omit':
        nan_policy = 'omit'

    if contains_nan and nan_policy == 'propagate':
        return KendalltauResult(np.nan, np.nan)

    elif contains_nan and nan_policy == 'omit':
        x = ma.masked_invalid(x)
        y = ma.masked_invalid(y)
        return mstats_basic.kendalltau(x, y, method=method)

    if initial_lexsort is not None:  # deprecate to drop!
        warnings.warn('"initial_lexsort" is gone!')

    def count_rank_tie(ranks):
        cnt = np.bincount(ranks).astype('int64', copy=False)
        cnt = cnt[cnt > 1]
        return ((cnt * (cnt - 1) // 2).sum(),
            (cnt * (cnt - 1.) * (cnt - 2)).sum(),
            (cnt * (cnt - 1.) * (2*cnt + 5)).sum())

    size = x.size
    perm = np.argsort(y)  # sort on y and convert y to dense ranks
    x, y = x[perm], y[perm]
    y = np.r_[True, y[1:] != y[:-1]].cumsum(dtype=np.intp)

    # stable sort on x and convert x to dense ranks
    perm = np.argsort(x, kind='mergesort')
    x, y = x[perm], y[perm]
    x = np.r_[True, x[1:] != x[:-1]].cumsum(dtype=np.intp)

    dis = _kendall_dis(x, y)  # discordant pairs

    obs = np.r_[True, (x[1:] != x[:-1]) | (y[1:] != y[:-1]), True]
    cnt = np.diff(np.nonzero(obs)[0]).astype('int64', copy=False)

    ntie = (cnt * (cnt - 1) // 2).sum()  # joint ties
    xtie, x0, x1 = count_rank_tie(x)     # ties in x, stats
    ytie, y0, y1 = count_rank_tie(y)     # ties in y, stats

    tot = (size * (size - 1)) // 2

    if xtie == tot or ytie == tot:
        return KendalltauResult(np.nan, np.nan)

    # Note that tot = con + dis + (xtie - ntie) + (ytie - ntie) + ntie
    #               = con + dis + xtie + ytie - ntie
    con_minus_dis = tot - xtie - ytie + ntie - 2 * dis
    tau = con_minus_dis / np.sqrt(tot - xtie) / np.sqrt(tot - ytie)
    # Limit range to fix computational errors
    tau = min(1., max(-1., tau))

    if method == 'exact' and (xtie != 0 or ytie != 0):
        raise ValueError("Ties found, exact method cannot be used.")

    if method == 'auto':
        if (xtie == 0 and ytie == 0) and (size <= 33 or min(dis, tot-dis) <= 1):
            method = 'exact'
        else:
            method = 'asymptotic'

    if xtie == 0 and ytie == 0 and method == 'exact':
        # Exact p-value, see Maurice G. Kendall, "Rank Correlation Methods" (4th Edition), Charles Griffin & Co., 1970.
        c = min(dis, tot-dis)
        if size <= 0:
            raise ValueError
        elif c < 0 or 2*c > size*(size-1):
            raise ValueError
        elif size == 1:
            pvalue = 1.0
        elif size == 2:
            pvalue = 1.0
        elif c == 0:
            pvalue = 2.0/math.factorial(size) if size < 171 else 0.0
        elif c == 1:
            pvalue = 2.0/math.factorial(size-1) if (size-1) < 171 else 0.0
        else:
            new = [0.0]*(c+1)
            new[0] = 1.0
            new[1] = 1.0
            for j in range(3,size+1):
                old = new[:]
                for k in range(1,min(j,c+1)):
                    new[k] += new[k-1]
                for k in range(j,c+1):
                    new[k] += new[k-1] - old[k-j]
            pvalue = 2.0*sum(new)/math.factorial(size) if size < 171 else 0.0

    elif method == 'asymptotic':
        # con_minus_dis is approx normally distributed with this variance [3]_
        var = (size * (size - 1) * (2.*size + 5) - x1 - y1) / 18. + (
            2. * xtie * ytie) / (size * (size - 1)) + x0 * y0 / (9. *
            size * (size - 1) * (size - 2))
        pvalue = special.erfc(np.abs(con_minus_dis) / np.sqrt(var) / np.sqrt(2))
    else:
        raise ValueError("Unknown method "+str(method)+" specified, please use auto, exact or asymptotic.")

    return KendalltauResult(tau, pvalue)


WeightedTauResult = namedtuple('WeightedTauResult', ('correlation', 'pvalue'))


def weightedtau(x, y, rank=True, weigher=None, additive=True):
    r"""
    Compute a weighted version of Kendall's :math:`\tau`.

    The weighted :math:`\tau` is a weighted version of Kendall's
    :math:`\tau` in which exchanges of high weight are more influential than
    exchanges of low weight. The default parameters compute the additive
    hyperbolic version of the index, :math:`\tau_\mathrm h`, which has
    been shown to provide the best balance between important and
    unimportant elements [1]_.

    The weighting is defined by means of a rank array, which assigns a
    nonnegative rank to each element, and a weigher function, which
    assigns a weight based from the rank to each element. The weight of an
    exchange is then the sum or the product of the weights of the ranks of
    the exchanged elements. The default parameters compute
    :math:`\tau_\mathrm h`: an exchange between elements with rank
    :math:`r` and :math:`s` (starting from zero) has weight
    :math:`1/(r+1) + 1/(s+1)`.

    Specifying a rank array is meaningful only if you have in mind an
    external criterion of importance. If, as it usually happens, you do
    not have in mind a specific rank, the weighted :math:`\tau` is
    defined by averaging the values obtained using the decreasing
    lexicographical rank by (`x`, `y`) and by (`y`, `x`). This is the
    behavior with default parameters.

    Note that if you are computing the weighted :math:`\tau` on arrays of
    ranks, rather than of scores (i.e., a larger value implies a lower
    rank) you must negate the ranks, so that elements of higher rank are
    associated with a larger value.

    Parameters
    ----------
    x, y : array_like
        Arrays of scores, of the same shape. If arrays are not 1-D, they will
        be flattened to 1-D.
    rank: array_like of ints or bool, optional
        A nonnegative rank assigned to each element. If it is None, the
        decreasing lexicographical rank by (`x`, `y`) will be used: elements of
        higher rank will be those with larger `x`-values, using `y`-values to
        break ties (in particular, swapping `x` and `y` will give a different
        result). If it is False, the element indices will be used
        directly as ranks. The default is True, in which case this
        function returns the average of the values obtained using the
        decreasing lexicographical rank by (`x`, `y`) and by (`y`, `x`).
    weigher : callable, optional
        The weigher function. Must map nonnegative integers (zero
        representing the most important element) to a nonnegative weight.
        The default, None, provides hyperbolic weighing, that is,
        rank :math:`r` is mapped to weight :math:`1/(r+1)`.
    additive : bool, optional
        If True, the weight of an exchange is computed by adding the
        weights of the ranks of the exchanged elements; otherwise, the weights
        are multiplied. The default is True.

    Returns
    -------
    correlation : float
       The weighted :math:`\tau` correlation index.
    pvalue : float
       Presently ``np.nan``, as the null statistics is unknown (even in the
       additive hyperbolic case).

    See also
    --------
    kendalltau : Calculates Kendall's tau.
    spearmanr : Calculates a Spearman rank-order correlation coefficient.
    theilslopes : Computes the Theil-Sen estimator for a set of points (x, y).

    Notes
    -----
    This function uses an :math:`O(n \log n)`, mergesort-based algorithm
    [1]_ that is a weighted extension of Knight's algorithm for Kendall's
    :math:`\tau` [2]_. It can compute Shieh's weighted :math:`\tau` [3]_
    between rankings without ties (i.e., permutations) by setting
    `additive` and `rank` to False, as the definition given in [1]_ is a
    generalization of Shieh's.

    NaNs are considered the smallest possible score.

    .. versionadded:: 0.19.0

    References
    ----------
    .. [1] Sebastiano Vigna, "A weighted correlation index for rankings with
           ties", Proceedings of the 24th international conference on World
           Wide Web, pp. 1166-1176, ACM, 2015.
    .. [2] W.R. Knight, "A Computer Method for Calculating Kendall's Tau with
           Ungrouped Data", Journal of the American Statistical Association,
           Vol. 61, No. 314, Part 1, pp. 436-439, 1966.
    .. [3] Grace S. Shieh. "A weighted Kendall's tau statistic", Statistics &
           Probability Letters, Vol. 39, No. 1, pp. 17-24, 1998.

    Examples
    --------
    >>> from scipy import stats
    >>> x = [12, 2, 1, 12, 2]
    >>> y = [1, 4, 7, 1, 0]
    >>> tau, p_value = stats.weightedtau(x, y)
    >>> tau
    -0.56694968153682723
    >>> p_value
    nan
    >>> tau, p_value = stats.weightedtau(x, y, additive=False)
    >>> tau
    -0.62205716951801038

    NaNs are considered the smallest possible score:

    >>> x = [12, 2, 1, 12, 2]
    >>> y = [1, 4, 7, 1, np.nan]
    >>> tau, _ = stats.weightedtau(x, y)
    >>> tau
    -0.56694968153682723

    This is exactly Kendall's tau:

    >>> x = [12, 2, 1, 12, 2]
    >>> y = [1, 4, 7, 1, 0]
    >>> tau, _ = stats.weightedtau(x, y, weigher=lambda x: 1)
    >>> tau
    -0.47140452079103173

    >>> x = [12, 2, 1, 12, 2]
    >>> y = [1, 4, 7, 1, 0]
    >>> stats.weightedtau(x, y, rank=None)
    WeightedTauResult(correlation=-0.4157652301037516, pvalue=nan)
    >>> stats.weightedtau(y, x, rank=None)
    WeightedTauResult(correlation=-0.7181341329699028, pvalue=nan)

    """
    x = np.asarray(x).ravel()
    y = np.asarray(y).ravel()

    if x.size != y.size:
        raise ValueError("All inputs to `weightedtau` must be of the same size, "
                         "found x-size %s and y-size %s" % (x.size, y.size))
    if not x.size:
        return WeightedTauResult(np.nan, np.nan)  # Return NaN if arrays are empty

    # If there are NaNs we apply _toint64()
    if np.isnan(np.sum(x)):
        x = _toint64(x)
    if np.isnan(np.sum(x)):
        y = _toint64(y)

    # Reduce to ranks unsupported types
    if x.dtype != y.dtype:
        if x.dtype != np.int64:
            x = _toint64(x)
        if y.dtype != np.int64:
            y = _toint64(y)
    else:
        if x.dtype not in (np.int32, np.int64, np.float32, np.float64):
            x = _toint64(x)
            y = _toint64(y)

    if rank is True:
        return WeightedTauResult((
            _weightedrankedtau(x, y, None, weigher, additive) +
            _weightedrankedtau(y, x, None, weigher, additive)
            ) / 2, np.nan)

    if rank is False:
        rank = np.arange(x.size, dtype=np.intp)
    elif rank is not None:
        rank = np.asarray(rank).ravel()
        if rank.size != x.size:
            raise ValueError("All inputs to `weightedtau` must be of the same size, "
                         "found x-size %s and rank-size %s" % (x.size, rank.size))

    return WeightedTauResult(_weightedrankedtau(x, y, rank, weigher, additive), np.nan)


#####################################
#       INFERENTIAL STATISTICS      #
#####################################

Ttest_1sampResult = namedtuple('Ttest_1sampResult', ('statistic', 'pvalue'))


def ttest_1samp(a, popmean, axis=0, nan_policy='propagate'):
    """
    Calculate the T-test for the mean of ONE group of scores.

    This is a two-sided test for the null hypothesis that the expected value
    (mean) of a sample of independent observations `a` is equal to the given
    population mean, `popmean`.

    Parameters
    ----------
    a : array_like
        sample observation
    popmean : float or array_like
        expected value in null hypothesis. If array_like, then it must have the
        same shape as `a` excluding the axis dimension
    axis : int or None, optional
        Axis along which to compute test. If None, compute over the whole
        array `a`.
    nan_policy : {'propagate', 'raise', 'omit'}, optional
        Defines how to handle when input contains nan. 'propagate' returns nan,
        'raise' throws an error, 'omit' performs the calculations ignoring nan
        values. Default is 'propagate'.

    Returns
    -------
    statistic : float or array
        t-statistic
    pvalue : float or array
        two-tailed p-value

    Examples
    --------
    >>> from scipy import stats

    >>> np.random.seed(7654567)  # fix seed to get the same result
    >>> rvs = stats.norm.rvs(loc=5, scale=10, size=(50,2))

    Test if mean of random sample is equal to true mean, and different mean.
    We reject the null hypothesis in the second case and don't reject it in
    the first case.

    >>> stats.ttest_1samp(rvs,5.0)
    (array([-0.68014479, -0.04323899]), array([ 0.49961383,  0.96568674]))
    >>> stats.ttest_1samp(rvs,0.0)
    (array([ 2.77025808,  4.11038784]), array([ 0.00789095,  0.00014999]))

    Examples using axis and non-scalar dimension for population mean.

    >>> stats.ttest_1samp(rvs,[5.0,0.0])
    (array([-0.68014479,  4.11038784]), array([  4.99613833e-01,   1.49986458e-04]))
    >>> stats.ttest_1samp(rvs.T,[5.0,0.0],axis=1)
    (array([-0.68014479,  4.11038784]), array([  4.99613833e-01,   1.49986458e-04]))
    >>> stats.ttest_1samp(rvs,[[5.0],[0.0]])
    (array([[-0.68014479, -0.04323899],
           [ 2.77025808,  4.11038784]]), array([[  4.99613833e-01,   9.65686743e-01],
           [  7.89094663e-03,   1.49986458e-04]]))

    """
    a, axis = _chk_asarray(a, axis)

    contains_nan, nan_policy = _contains_nan(a, nan_policy)

    if contains_nan and nan_policy == 'omit':
        a = ma.masked_invalid(a)
        return mstats_basic.ttest_1samp(a, popmean, axis)

    n = a.shape[axis]
    df = n - 1

    d = np.mean(a, axis) - popmean
    v = np.var(a, axis, ddof=1)
    denom = np.sqrt(v / n)

    with np.errstate(divide='ignore', invalid='ignore'):
        t = np.divide(d, denom)
    t, prob = _ttest_finish(df, t)

    return Ttest_1sampResult(t, prob)


def _ttest_finish(df, t):
    """Common code between all 3 t-test functions."""
    prob = distributions.t.sf(np.abs(t), df) * 2  # use np.abs to get upper tail
    if t.ndim == 0:
        t = t[()]

    return t, prob


def _ttest_ind_from_stats(mean1, mean2, denom, df):

    d = mean1 - mean2
    with np.errstate(divide='ignore', invalid='ignore'):
        t = np.divide(d, denom)
    t, prob = _ttest_finish(df, t)

    return (t, prob)


def _unequal_var_ttest_denom(v1, n1, v2, n2):
    vn1 = v1 / n1
    vn2 = v2 / n2
    with np.errstate(divide='ignore', invalid='ignore'):
        df = (vn1 + vn2)**2 / (vn1**2 / (n1 - 1) + vn2**2 / (n2 - 1))

    # If df is undefined, variances are zero (assumes n1 > 0 & n2 > 0).
    # Hence it doesn't matter what df is as long as it's not NaN.
    df = np.where(np.isnan(df), 1, df)
    denom = np.sqrt(vn1 + vn2)
    return df, denom


def _equal_var_ttest_denom(v1, n1, v2, n2):
    df = n1 + n2 - 2.0
    svar = ((n1 - 1) * v1 + (n2 - 1) * v2) / df
    denom = np.sqrt(svar * (1.0 / n1 + 1.0 / n2))
    return df, denom


Ttest_indResult = namedtuple('Ttest_indResult', ('statistic', 'pvalue'))


def ttest_ind_from_stats(mean1, std1, nobs1, mean2, std2, nobs2,
                         equal_var=True):
    """
    T-test for means of two independent samples from descriptive statistics.

    This is a two-sided test for the null hypothesis that two independent
    samples have identical average (expected) values.

    Parameters
    ----------
    mean1 : array_like
        The mean(s) of sample 1.
    std1 : array_like
        The standard deviation(s) of sample 1.
    nobs1 : array_like
        The number(s) of observations of sample 1.
    mean2 : array_like
        The mean(s) of sample 2
    std2 : array_like
        The standard deviations(s) of sample 2.
    nobs2 : array_like
        The number(s) of observations of sample 2.
    equal_var : bool, optional
        If True (default), perform a standard independent 2 sample test
        that assumes equal population variances [1]_.
        If False, perform Welch's t-test, which does not assume equal
        population variance [2]_.

    Returns
    -------
    statistic : float or array
        The calculated t-statistics
    pvalue : float or array
        The two-tailed p-value.

    See Also
    --------
    scipy.stats.ttest_ind

    Notes
    -----

    .. versionadded:: 0.16.0

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/T-test#Independent_two-sample_t-test

    .. [2] https://en.wikipedia.org/wiki/Welch%27s_t-test

    Examples
    --------
    Suppose we have the summary data for two samples, as follows::

                         Sample   Sample
                   Size   Mean   Variance
        Sample 1    13    15.0     87.5
        Sample 2    11    12.0     39.0

    Apply the t-test to this data (with the assumption that the population
    variances are equal):

    >>> from scipy.stats import ttest_ind_from_stats
    >>> ttest_ind_from_stats(mean1=15.0, std1=np.sqrt(87.5), nobs1=13,
    ...                      mean2=12.0, std2=np.sqrt(39.0), nobs2=11)
    Ttest_indResult(statistic=0.9051358093310269, pvalue=0.3751996797581487)

    For comparison, here is the data from which those summary statistics
    were taken.  With this data, we can compute the same result using
    `scipy.stats.ttest_ind`:

    >>> a = np.array([1, 3, 4, 6, 11, 13, 15, 19, 22, 24, 25, 26, 26])
    >>> b = np.array([2, 4, 6, 9, 11, 13, 14, 15, 18, 19, 21])
    >>> from scipy.stats import ttest_ind
    >>> ttest_ind(a, b)
    Ttest_indResult(statistic=0.905135809331027, pvalue=0.3751996797581486)

    """
    if equal_var:
        df, denom = _equal_var_ttest_denom(std1**2, nobs1, std2**2, nobs2)
    else:
        df, denom = _unequal_var_ttest_denom(std1**2, nobs1,
                                             std2**2, nobs2)

    res = _ttest_ind_from_stats(mean1, mean2, denom, df)
    return Ttest_indResult(*res)


def ttest_ind(a, b, axis=0, equal_var=True, nan_policy='propagate'):
    """
    Calculate the T-test for the means of *two independent* samples of scores.

    This is a two-sided test for the null hypothesis that 2 independent samples
    have identical average (expected) values. This test assumes that the
    populations have identical variances by default.

    Parameters
    ----------
    a, b : array_like
        The arrays must have the same shape, except in the dimension
        corresponding to `axis` (the first, by default).
    axis : int or None, optional
        Axis along which to compute test. If None, compute over the whole
        arrays, `a`, and `b`.
    equal_var : bool, optional
        If True (default), perform a standard independent 2 sample test
        that assumes equal population variances [1]_.
        If False, perform Welch's t-test, which does not assume equal
        population variance [2]_.

        .. versionadded:: 0.11.0
    nan_policy : {'propagate', 'raise', 'omit'}, optional
        Defines how to handle when input contains nan. 'propagate' returns nan,
        'raise' throws an error, 'omit' performs the calculations ignoring nan
        values. Default is 'propagate'.


    Returns
    -------
    statistic : float or array
        The calculated t-statistic.
    pvalue : float or array
        The two-tailed p-value.

    Notes
    -----
    We can use this test, if we observe two independent samples from
    the same or different population, e.g. exam scores of boys and
    girls or of two ethnic groups. The test measures whether the
    average (expected) value differs significantly across samples. If
    we observe a large p-value, for example larger than 0.05 or 0.1,
    then we cannot reject the null hypothesis of identical average scores.
    If the p-value is smaller than the threshold, e.g. 1%, 5% or 10%,
    then we reject the null hypothesis of equal averages.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/T-test#Independent_two-sample_t-test

    .. [2] https://en.wikipedia.org/wiki/Welch%27s_t-test

    Examples
    --------
    >>> from scipy import stats
    >>> np.random.seed(12345678)

    Test with sample with identical means:

    >>> rvs1 = stats.norm.rvs(loc=5,scale=10,size=500)
    >>> rvs2 = stats.norm.rvs(loc=5,scale=10,size=500)
    >>> stats.ttest_ind(rvs1,rvs2)
    (0.26833823296239279, 0.78849443369564776)
    >>> stats.ttest_ind(rvs1,rvs2, equal_var = False)
    (0.26833823296239279, 0.78849452749500748)

    `ttest_ind` underestimates p for unequal variances:

    >>> rvs3 = stats.norm.rvs(loc=5, scale=20, size=500)
    >>> stats.ttest_ind(rvs1, rvs3)
    (-0.46580283298287162, 0.64145827413436174)
    >>> stats.ttest_ind(rvs1, rvs3, equal_var = False)
    (-0.46580283298287162, 0.64149646246569292)

    When n1 != n2, the equal variance t-statistic is no longer equal to the
    unequal variance t-statistic:

    >>> rvs4 = stats.norm.rvs(loc=5, scale=20, size=100)
    >>> stats.ttest_ind(rvs1, rvs4)
    (-0.99882539442782481, 0.3182832709103896)
    >>> stats.ttest_ind(rvs1, rvs4, equal_var = False)
    (-0.69712570584654099, 0.48716927725402048)

    T-test with different means, variance, and n:

    >>> rvs5 = stats.norm.rvs(loc=8, scale=20, size=100)
    >>> stats.ttest_ind(rvs1, rvs5)
    (-1.4679669854490653, 0.14263895620529152)
    >>> stats.ttest_ind(rvs1, rvs5, equal_var = False)
    (-0.94365973617132992, 0.34744170334794122)

    """
    a, b, axis = _chk2_asarray(a, b, axis)

    # check both a and b
    cna, npa = _contains_nan(a, nan_policy)
    cnb, npb = _contains_nan(b, nan_policy)
    contains_nan = cna or cnb
    if npa == 'omit' or npb == 'omit':
        nan_policy = 'omit'

    if contains_nan and nan_policy == 'omit':
        a = ma.masked_invalid(a)
        b = ma.masked_invalid(b)
        return mstats_basic.ttest_ind(a, b, axis, equal_var)

    if a.size == 0 or b.size == 0:
        return Ttest_indResult(np.nan, np.nan)

    v1 = np.var(a, axis, ddof=1)
    v2 = np.var(b, axis, ddof=1)
    n1 = a.shape[axis]
    n2 = b.shape[axis]

    if equal_var:
        df, denom = _equal_var_ttest_denom(v1, n1, v2, n2)
    else:
        df, denom = _unequal_var_ttest_denom(v1, n1, v2, n2)

    res = _ttest_ind_from_stats(np.mean(a, axis), np.mean(b, axis), denom, df)

    return Ttest_indResult(*res)


Ttest_relResult = namedtuple('Ttest_relResult', ('statistic', 'pvalue'))


def ttest_rel(a, b, axis=0, nan_policy='propagate'):
    """
    Calculate the T-test on TWO RELATED samples of scores, a and b.

    This is a two-sided test for the null hypothesis that 2 related or
    repeated samples have identical average (expected) values.

    Parameters
    ----------
    a, b : array_like
        The arrays must have the same shape.
    axis : int or None, optional
        Axis along which to compute test. If None, compute over the whole
        arrays, `a`, and `b`.
    nan_policy : {'propagate', 'raise', 'omit'}, optional
        Defines how to handle when input contains nan. 'propagate' returns nan,
        'raise' throws an error, 'omit' performs the calculations ignoring nan
        values. Default is 'propagate'.

    Returns
    -------
    statistic : float or array
        t-statistic
    pvalue : float or array
        two-tailed p-value

    Notes
    -----
    Examples for the use are scores of the same set of student in
    different exams, or repeated sampling from the same units. The
    test measures whether the average score differs significantly
    across samples (e.g. exams). If we observe a large p-value, for
    example greater than 0.05 or 0.1 then we cannot reject the null
    hypothesis of identical average scores. If the p-value is smaller
    than the threshold, e.g. 1%, 5% or 10%, then we reject the null
    hypothesis of equal averages. Small p-values are associated with
    large t-statistics.

    References
    ----------
    https://en.wikipedia.org/wiki/T-test#Dependent_t-test_for_paired_samples

    Examples
    --------
    >>> from scipy import stats
    >>> np.random.seed(12345678) # fix random seed to get same numbers

    >>> rvs1 = stats.norm.rvs(loc=5,scale=10,size=500)
    >>> rvs2 = (stats.norm.rvs(loc=5,scale=10,size=500) +
    ...         stats.norm.rvs(scale=0.2,size=500))
    >>> stats.ttest_rel(rvs1,rvs2)
    (0.24101764965300962, 0.80964043445811562)
    >>> rvs3 = (stats.norm.rvs(loc=8,scale=10,size=500) +
    ...         stats.norm.rvs(scale=0.2,size=500))
    >>> stats.ttest_rel(rvs1,rvs3)
    (-3.9995108708727933, 7.3082402191726459e-005)

    """
    a, b, axis = _chk2_asarray(a, b, axis)

    cna, npa = _contains_nan(a, nan_policy)
    cnb, npb = _contains_nan(b, nan_policy)
    contains_nan = cna or cnb
    if npa == 'omit' or npb == 'omit':
        nan_policy = 'omit'

    if contains_nan and nan_policy == 'omit':
        a = ma.masked_invalid(a)
        b = ma.masked_invalid(b)
        m = ma.mask_or(ma.getmask(a), ma.getmask(b))
        aa = ma.array(a, mask=m, copy=True)
        bb = ma.array(b, mask=m, copy=True)
        return mstats_basic.ttest_rel(aa, bb, axis)

    if a.shape[axis] != b.shape[axis]:
        raise ValueError('unequal length arrays')

    if a.size == 0 or b.size == 0:
        return np.nan, np.nan

    n = a.shape[axis]
    df = n - 1

    d = (a - b).astype(np.float64)
    v = np.var(d, axis, ddof=1)
    dm = np.mean(d, axis)
    denom = np.sqrt(v / n)

    with np.errstate(divide='ignore', invalid='ignore'):
        t = np.divide(dm, denom)
    t, prob = _ttest_finish(df, t)

    return Ttest_relResult(t, prob)


KstestResult = namedtuple('KstestResult', ('statistic', 'pvalue'))


def kstest(rvs, cdf, args=(), N=20, alternative='two-sided', mode='approx'):
    """
    Perform the Kolmogorov-Smirnov test for goodness of fit.

    This performs a test of the distribution F(x) of an observed
    random variable against a given distribution G(x). Under the null
    hypothesis the two distributions are identical, F(x)=G(x). The
    alternative hypothesis can be either 'two-sided' (default), 'less'
    or 'greater'. The KS test is only valid for continuous distributions.

    Parameters
    ----------
    rvs : str, array or callable
        If a string, it should be the name of a distribution in `scipy.stats`.
        If an array, it should be a 1-D array of observations of random
        variables.
        If a callable, it should be a function to generate random variables;
        it is required to have a keyword argument `size`.
    cdf : str or callable
        If a string, it should be the name of a distribution in `scipy.stats`.
        If `rvs` is a string then `cdf` can be False or the same as `rvs`.
        If a callable, that callable is used to calculate the cdf.
    args : tuple, sequence, optional
        Distribution parameters, used if `rvs` or `cdf` are strings.
    N : int, optional
        Sample size if `rvs` is string or callable.  Default is 20.
    alternative : {'two-sided', 'less','greater'}, optional
        Defines the alternative hypothesis (see explanation above).
        Default is 'two-sided'.
    mode : 'approx' (default) or 'asymp', optional
        Defines the distribution used for calculating the p-value.

          - 'approx' : use approximation to exact distribution of test statistic
          - 'asymp' : use asymptotic distribution of test statistic

    Returns
    -------
    statistic : float
        KS test statistic, either D, D+ or D-.
    pvalue :  float
        One-tailed or two-tailed p-value.

    Notes
    -----
    In the one-sided test, the alternative is that the empirical
    cumulative distribution function of the random variable is "less"
    or "greater" than the cumulative distribution function G(x) of the
    hypothesis, ``F(x)<=G(x)``, resp. ``F(x)>=G(x)``.

    Examples
    --------
    >>> from scipy import stats

    >>> x = np.linspace(-15, 15, 9)
    >>> stats.kstest(x, 'norm')
    (0.44435602715924361, 0.038850142705171065)

    >>> np.random.seed(987654321) # set random seed to get the same result
    >>> stats.kstest('norm', False, N=100)
    (0.058352892479417884, 0.88531190944151261)

    The above lines are equivalent to:

    >>> np.random.seed(987654321)
    >>> stats.kstest(stats.norm.rvs(size=100), 'norm')
    (0.058352892479417884, 0.88531190944151261)

    *Test against one-sided alternative hypothesis*

    Shift distribution to larger values, so that ``cdf_dgp(x) < norm.cdf(x)``:

    >>> np.random.seed(987654321)
    >>> x = stats.norm.rvs(loc=0.2, size=100)
    >>> stats.kstest(x,'norm', alternative = 'less')
    (0.12464329735846891, 0.040989164077641749)

    Reject equal distribution against alternative hypothesis: less

    >>> stats.kstest(x,'norm', alternative = 'greater')
    (0.0072115233216311081, 0.98531158590396395)

    Don't reject equal distribution against alternative hypothesis: greater

    >>> stats.kstest(x,'norm', mode='asymp')
    (0.12464329735846891, 0.08944488871182088)

    *Testing t distributed random variables against normal distribution*

    With 100 degrees of freedom the t distribution looks close to the normal
    distribution, and the K-S test does not reject the hypothesis that the
    sample came from the normal distribution:

    >>> np.random.seed(987654321)
    >>> stats.kstest(stats.t.rvs(100,size=100),'norm')
    (0.072018929165471257, 0.67630062862479168)

    With 3 degrees of freedom the t distribution looks sufficiently different
    from the normal distribution, that we can reject the hypothesis that the
    sample came from the normal distribution at the 10% level:

    >>> np.random.seed(987654321)
    >>> stats.kstest(stats.t.rvs(3,size=100),'norm')
    (0.131016895759829, 0.058826222555312224)

    """
    if isinstance(rvs, string_types):
        if (not cdf) or (cdf == rvs):
            cdf = getattr(distributions, rvs).cdf
            rvs = getattr(distributions, rvs).rvs
        else:
            raise AttributeError("if rvs is string, cdf has to be the "
                                 "same distribution")

    if isinstance(cdf, string_types):
        cdf = getattr(distributions, cdf).cdf
    if callable(rvs):
        kwds = {'size': N}
        vals = np.sort(rvs(*args, **kwds))
    else:
        vals = np.sort(rvs)
        N = len(vals)
    cdfvals = cdf(vals, *args)

    # to not break compatibility with existing code
    if alternative == 'two_sided':
        alternative = 'two-sided'

    if alternative in ['two-sided', 'greater']:
        Dplus = (np.arange(1.0, N + 1)/N - cdfvals).max()
        if alternative == 'greater':
            return KstestResult(Dplus, distributions.ksone.sf(Dplus, N))

    if alternative in ['two-sided', 'less']:
        Dmin = (cdfvals - np.arange(0.0, N)/N).max()
        if alternative == 'less':
            return KstestResult(Dmin, distributions.ksone.sf(Dmin, N))

    if alternative == 'two-sided':
        D = np.max([Dplus, Dmin])
        if mode == 'asymp':
            return KstestResult(D, distributions.kstwobign.sf(D * np.sqrt(N)))
        if mode == 'approx':
            pval_two = distributions.kstwobign.sf(D * np.sqrt(N))
            if N > 2666 or pval_two > 0.80 - N*0.3/1000:
                return KstestResult(D, pval_two)
            else:
                return KstestResult(D, 2 * distributions.ksone.sf(D, N))


# Map from names to lambda_ values used in power_divergence().
_power_div_lambda_names = {
    "pearson": 1,
    "log-likelihood": 0,
    "freeman-tukey": -0.5,
    "mod-log-likelihood": -1,
    "neyman": -2,
    "cressie-read": 2/3,
}


def _count(a, axis=None):
    """
    Count the number of non-masked elements of an array.

    This function behaves like np.ma.count(), but is much faster
    for ndarrays.
    """
    if hasattr(a, 'count'):
        num = a.count(axis=axis)
        if isinstance(num, np.ndarray) and num.ndim == 0:
            # In some cases, the `count` method returns a scalar array (e.g.
            # np.array(3)), but we want a plain integer.
            num = int(num)
    else:
        if axis is None:
            num = a.size
        else:
            num = a.shape[axis]
    return num


Power_divergenceResult = namedtuple('Power_divergenceResult',
                                    ('statistic', 'pvalue'))

def power_divergence(f_obs, f_exp=None, ddof=0, axis=0, lambda_=None):
    """
    Cressie-Read power divergence statistic and goodness of fit test.

    This function tests the null hypothesis that the categorical data
    has the given frequencies, using the Cressie-Read power divergence
    statistic.

    Parameters
    ----------
    f_obs : array_like
        Observed frequencies in each category.
    f_exp : array_like, optional
        Expected frequencies in each category.  By default the categories are
        assumed to be equally likely.
    ddof : int, optional
        "Delta degrees of freedom": adjustment to the degrees of freedom
        for the p-value.  The p-value is computed using a chi-squared
        distribution with ``k - 1 - ddof`` degrees of freedom, where `k`
        is the number of observed frequencies.  The default value of `ddof`
        is 0.
    axis : int or None, optional
        The axis of the broadcast result of `f_obs` and `f_exp` along which to
        apply the test.  If axis is None, all values in `f_obs` are treated
        as a single data set.  Default is 0.
    lambda_ : float or str, optional
        `lambda_` gives the power in the Cressie-Read power divergence
        statistic.  The default is 1.  For convenience, `lambda_` may be
        assigned one of the following strings, in which case the
        corresponding numerical value is used::

            String              Value   Description
            "pearson"             1     Pearson's chi-squared statistic.
                                        In this case, the function is
                                        equivalent to `stats.chisquare`.
            "log-likelihood"      0     Log-likelihood ratio. Also known as
                                        the G-test [3]_.
            "freeman-tukey"      -1/2   Freeman-Tukey statistic.
            "mod-log-likelihood" -1     Modified log-likelihood ratio.
            "neyman"             -2     Neyman's statistic.
            "cressie-read"        2/3   The power recommended in [5]_.

    Returns
    -------
    statistic : float or ndarray
        The Cressie-Read power divergence test statistic.  The value is
        a float if `axis` is None or if` `f_obs` and `f_exp` are 1-D.
    pvalue : float or ndarray
        The p-value of the test.  The value is a float if `ddof` and the
        return value `stat` are scalars.

    See Also
    --------
    chisquare

    Notes
    -----
    This test is invalid when the observed or expected frequencies in each
    category are too small.  A typical rule is that all of the observed
    and expected frequencies should be at least 5.

    When `lambda_` is less than zero, the formula for the statistic involves
    dividing by `f_obs`, so a warning or error may be generated if any value
    in `f_obs` is 0.

    Similarly, a warning or error may be generated if any value in `f_exp` is
    zero when `lambda_` >= 0.

    The default degrees of freedom, k-1, are for the case when no parameters
    of the distribution are estimated. If p parameters are estimated by
    efficient maximum likelihood then the correct degrees of freedom are
    k-1-p. If the parameters are estimated in a different way, then the
    dof can be between k-1-p and k-1. However, it is also possible that
    the asymptotic distribution is not a chisquare, in which case this
    test is not appropriate.

    This function handles masked arrays.  If an element of `f_obs` or `f_exp`
    is masked, then data at that position is ignored, and does not count
    towards the size of the data set.

    .. versionadded:: 0.13.0

    References
    ----------
    .. [1] Lowry, Richard.  "Concepts and Applications of Inferential
           Statistics". Chapter 8.
           https://web.archive.org/web/20171015035606/http://faculty.vassar.edu/lowry/ch8pt1.html
    .. [2] "Chi-squared test", https://en.wikipedia.org/wiki/Chi-squared_test
    .. [3] "G-test", https://en.wikipedia.org/wiki/G-test
    .. [4] Sokal, R. R. and Rohlf, F. J. "Biometry: the principles and
           practice of statistics in biological research", New York: Freeman
           (1981)
    .. [5] Cressie, N. and Read, T. R. C., "Multinomial Goodness-of-Fit
           Tests", J. Royal Stat. Soc. Series B, Vol. 46, No. 3 (1984),
           pp. 440-464.

    Examples
    --------

    (See `chisquare` for more examples.)

    When just `f_obs` is given, it is assumed that the expected frequencies
    are uniform and given by the mean of the observed frequencies.  Here we
    perform a G-test (i.e. use the log-likelihood ratio statistic):

    >>> from scipy.stats import power_divergence
    >>> power_divergence([16, 18, 16, 14, 12, 12], lambda_='log-likelihood')
    (2.006573162632538, 0.84823476779463769)

    The expected frequencies can be given with the `f_exp` argument:

    >>> power_divergence([16, 18, 16, 14, 12, 12],
    ...                  f_exp=[16, 16, 16, 16, 16, 8],
    ...                  lambda_='log-likelihood')
    (3.3281031458963746, 0.6495419288047497)

    When `f_obs` is 2-D, by default the test is applied to each column.

    >>> obs = np.array([[16, 18, 16, 14, 12, 12], [32, 24, 16, 28, 20, 24]]).T
    >>> obs.shape
    (6, 2)
    >>> power_divergence(obs, lambda_="log-likelihood")
    (array([ 2.00657316,  6.77634498]), array([ 0.84823477,  0.23781225]))

    By setting ``axis=None``, the test is applied to all data in the array,
    which is equivalent to applying the test to the flattened array.

    >>> power_divergence(obs, axis=None)
    (23.31034482758621, 0.015975692534127565)
    >>> power_divergence(obs.ravel())
    (23.31034482758621, 0.015975692534127565)

    `ddof` is the change to make to the default degrees of freedom.

    >>> power_divergence([16, 18, 16, 14, 12, 12], ddof=1)
    (2.0, 0.73575888234288467)

    The calculation of the p-values is done by broadcasting the
    test statistic with `ddof`.

    >>> power_divergence([16, 18, 16, 14, 12, 12], ddof=[0,1,2])
    (2.0, array([ 0.84914504,  0.73575888,  0.5724067 ]))

    `f_obs` and `f_exp` are also broadcast.  In the following, `f_obs` has
    shape (6,) and `f_exp` has shape (2, 6), so the result of broadcasting
    `f_obs` and `f_exp` has shape (2, 6).  To compute the desired chi-squared
    statistics, we must use ``axis=1``:

    >>> power_divergence([16, 18, 16, 14, 12, 12],
    ...                  f_exp=[[16, 16, 16, 16, 16, 8],
    ...                         [8, 20, 20, 16, 12, 12]],
    ...                  axis=1)
    (array([ 3.5 ,  9.25]), array([ 0.62338763,  0.09949846]))

    """
    # Convert the input argument `lambda_` to a numerical value.
    if isinstance(lambda_, string_types):
        if lambda_ not in _power_div_lambda_names:
            names = repr(list(_power_div_lambda_names.keys()))[1:-1]
            raise ValueError("invalid string for lambda_: {0!r}.  Valid strings "
                             "are {1}".format(lambda_, names))
        lambda_ = _power_div_lambda_names[lambda_]
    elif lambda_ is None:
        lambda_ = 1

    f_obs = np.asanyarray(f_obs)

    if f_exp is not None:
        f_exp = np.atleast_1d(np.asanyarray(f_exp))
    else:
        # Compute the equivalent of
        #   f_exp = f_obs.mean(axis=axis, keepdims=True)
        # Older versions of numpy do not have the 'keepdims' argument, so
        # we have to do a little work to achieve the same result.
        # Ignore 'invalid' errors so the edge case of a data set with length 0
        # is handled without spurious warnings.
        with np.errstate(invalid='ignore'):
            f_exp = np.atleast_1d(f_obs.mean(axis=axis))
        if axis is not None:
            reduced_shape = list(f_obs.shape)
            reduced_shape[axis] = 1
            f_exp.shape = reduced_shape

    # `terms` is the array of terms that are summed along `axis` to create
    # the test statistic.  We use some specialized code for a few special
    # cases of lambda_.
    if lambda_ == 1:
        # Pearson's chi-squared statistic
        terms = (f_obs - f_exp)**2 / f_exp
    elif lambda_ == 0:
        # Log-likelihood ratio (i.e. G-test)
        terms = 2.0 * special.xlogy(f_obs, f_obs / f_exp)
    elif lambda_ == -1:
        # Modified log-likelihood ratio
        terms = 2.0 * special.xlogy(f_exp, f_exp / f_obs)
    else:
        # General Cressie-Read power divergence.
        terms = f_obs * ((f_obs / f_exp)**lambda_ - 1)
        terms /= 0.5 * lambda_ * (lambda_ + 1)

    stat = terms.sum(axis=axis)

    num_obs = _count(terms, axis=axis)
    ddof = asarray(ddof)
    p = distributions.chi2.sf(stat, num_obs - 1 - ddof)

    return Power_divergenceResult(stat, p)


def chisquare(f_obs, f_exp=None, ddof=0, axis=0):
    """
    Calculate a one-way chi square test.

    The chi square test tests the null hypothesis that the categorical data
    has the given frequencies.

    Parameters
    ----------
    f_obs : array_like
        Observed frequencies in each category.
    f_exp : array_like, optional
        Expected frequencies in each category.  By default the categories are
        assumed to be equally likely.
    ddof : int, optional
        "Delta degrees of freedom": adjustment to the degrees of freedom
        for the p-value.  The p-value is computed using a chi-squared
        distribution with ``k - 1 - ddof`` degrees of freedom, where `k`
        is the number of observed frequencies.  The default value of `ddof`
        is 0.
    axis : int or None, optional
        The axis of the broadcast result of `f_obs` and `f_exp` along which to
        apply the test.  If axis is None, all values in `f_obs` are treated
        as a single data set.  Default is 0.

    Returns
    -------
    chisq : float or ndarray
        The chi-squared test statistic.  The value is a float if `axis` is
        None or `f_obs` and `f_exp` are 1-D.
    p : float or ndarray
        The p-value of the test.  The value is a float if `ddof` and the
        return value `chisq` are scalars.

    See Also
    --------
    scipy.stats.power_divergence

    Notes
    -----
    This test is invalid when the observed or expected frequencies in each
    category are too small.  A typical rule is that all of the observed
    and expected frequencies should be at least 5.

    The default degrees of freedom, k-1, are for the case when no parameters
    of the distribution are estimated. If p parameters are estimated by
    efficient maximum likelihood then the correct degrees of freedom are
    k-1-p. If the parameters are estimated in a different way, then the
    dof can be between k-1-p and k-1. However, it is also possible that
    the asymptotic distribution is not a chisquare, in which case this
    test is not appropriate.

    References
    ----------
    .. [1] Lowry, Richard.  "Concepts and Applications of Inferential
           Statistics". Chapter 8.
           https://web.archive.org/web/20171022032306/http://vassarstats.net:80/textbook/ch8pt1.html
    .. [2] "Chi-squared test", https://en.wikipedia.org/wiki/Chi-squared_test

    Examples
    --------
    When just `f_obs` is given, it is assumed that the expected frequencies
    are uniform and given by the mean of the observed frequencies.

    >>> from scipy.stats import chisquare
    >>> chisquare([16, 18, 16, 14, 12, 12])
    (2.0, 0.84914503608460956)

    With `f_exp` the expected frequencies can be given.

    >>> chisquare([16, 18, 16, 14, 12, 12], f_exp=[16, 16, 16, 16, 16, 8])
    (3.5, 0.62338762774958223)

    When `f_obs` is 2-D, by default the test is applied to each column.

    >>> obs = np.array([[16, 18, 16, 14, 12, 12], [32, 24, 16, 28, 20, 24]]).T
    >>> obs.shape
    (6, 2)
    >>> chisquare(obs)
    (array([ 2.        ,  6.66666667]), array([ 0.84914504,  0.24663415]))

    By setting ``axis=None``, the test is applied to all data in the array,
    which is equivalent to applying the test to the flattened array.

    >>> chisquare(obs, axis=None)
    (23.31034482758621, 0.015975692534127565)
    >>> chisquare(obs.ravel())
    (23.31034482758621, 0.015975692534127565)

    `ddof` is the change to make to the default degrees of freedom.

    >>> chisquare([16, 18, 16, 14, 12, 12], ddof=1)
    (2.0, 0.73575888234288467)

    The calculation of the p-values is done by broadcasting the
    chi-squared statistic with `ddof`.

    >>> chisquare([16, 18, 16, 14, 12, 12], ddof=[0,1,2])
    (2.0, array([ 0.84914504,  0.73575888,  0.5724067 ]))

    `f_obs` and `f_exp` are also broadcast.  In the following, `f_obs` has
    shape (6,) and `f_exp` has shape (2, 6), so the result of broadcasting
    `f_obs` and `f_exp` has shape (2, 6).  To compute the desired chi-squared
    statistics, we use ``axis=1``:

    >>> chisquare([16, 18, 16, 14, 12, 12],
    ...           f_exp=[[16, 16, 16, 16, 16, 8], [8, 20, 20, 16, 12, 12]],
    ...           axis=1)
    (array([ 3.5 ,  9.25]), array([ 0.62338763,  0.09949846]))

    """
    return power_divergence(f_obs, f_exp=f_exp, ddof=ddof, axis=axis,
                            lambda_="pearson")


Ks_2sampResult = namedtuple('Ks_2sampResult', ('statistic', 'pvalue'))


def _compute_prob_inside_method(m, n, g, h):
    """Count the proportion of paths that stay strictly inside two diagonal lines.

    Parameters
    ----------
    m : integer
        m > 0
    n : integer
        n > 0
    g : integer
        g is greatest common divisor of m and n
    h : integer
        0 <= h <= lcm(m,n)

    Returns
    -------
    p : float
        The proportion of paths that stay inside the two lines.


    Count the integer lattice paths from (0, 0) to (m, n) which satisfy
    |x/m - y/n| < h / lcm(m, n).
    The paths make steps of size +1 in either positive x or positive y directions.

    We generally follow Hodges' treatment of Drion/Gnedenko/Korolyuk.
    Hodges, J.L. Jr.,
    "The Significance Probability of the Smirnov Two-Sample Test,"
    Arkiv fiur Matematik, 3, No. 43 (1958), 469-86.
    """
    # Probability is symmetrical in m, n.  Computation below uses m >= n.
    if m < n:
        m, n = n, m
    mg = m // g
    ng = n // g

    # Count the integer lattice paths from (0, 0) to (m, n) which satisfy
    # |nx/g - my/g| < h.
    # Compute matrix A such that:
    #  A(x, 0) = A(0, y) = 1
    #  A(x, y) = A(x, y-1) + A(x-1, y), for x,y>=1, except that
    #  A(x, y) = 0 if |x/m - y/n|>= h
    # Probability is A(m, n)/binom(m+n, n)
    # Optimizations exist for m==n, m==n*p.
    # Only need to preserve a single column of A, and only a sliding window of it.
    # minj keeps track of the slide.
    minj, maxj = 0, min(int(np.ceil(h / mg)), n + 1)
    curlen = maxj - minj
    # Make a vector long enough to hold maximum window needed.
    lenA = min(2 * maxj + 2, n + 1)
    # This is an integer calculation, but the entries are essentially
    # binomial coefficients, hence grow quickly.
    # Scaling after each column is computed avoids dividing by a
    # large binomial coefficent at the end. Instead it is incorporated
    # one factor at a time during the computation.
    dtype = np.float64
    A = np.zeros(lenA, dtype=dtype)
    # Initialize the first column
    A[minj:maxj] = 1
    for i in range(1, m + 1):
        # Generate the next column.
        # First calculate the sliding window
        lastminj, lastmaxj, lastlen = minj, maxj, curlen
        minj = max(int(np.floor((ng * i - h) / mg)) + 1, 0)
        minj = min(minj, n)
        maxj = min(int(np.ceil((ng * i + h) / mg)), n + 1)
        if maxj <= minj:
            return 0
        # Now fill in the values
        A[0:maxj - minj] = np.cumsum(A[minj - lastminj:maxj - lastminj])
        curlen = maxj - minj
        if lastlen > curlen:
            # Set some carried-over elements to 0
            A[maxj - minj:maxj - minj + (lastlen - curlen)] = 0
        # Peel off one term from each of top and bottom of the binomial coefficient.
        scaling_factor = i * 1.0 / (n + i)
        A *= scaling_factor
    return A[maxj - minj - 1]


def _compute_prob_outside_square(n, h):
    """Compute the proportion of paths that pass outside the two diagonal lines.

    Parameters
    ----------
    n : integer
        n > 0
    h : integer
        0 <= h <= n

    Returns
    -------
    p : float
        The proportion of paths that pass outside the lines x-y = +/-h.

    """
    # Compute Pr(D_{n,n} >= h/n)
    # Prob = 2 * ( binom(2n, n-h) - binom(2n, n-2a) + binom(2n, n-3a) - ... )  / binom(2n, n)
    # This formulation exhibits subtractive cancellation.
    # Instead divide each term by binom(2n, n), then factor common terms
    # and use a Horner-like algorithm
    # P = 2 * A0 * (1 - A1*(1 - A2*(1 - A3*(1 - A4*(...)))))

    P = 0.0
    k = int(np.floor(n / h))
    while k >= 0:
        p1 = 1.0
        # Each of the Ai terms has numerator and denominator with h simple terms.
        for j in range(h):
            p1 = (n - k * h - j) * p1 / (n + k * h + j + 1)
        P = p1 * (1.0 - P)
        k -= 1
    return 2 * P


def _count_paths_outside_method(m, n, g, h):
    """Count the number of paths that pass outside the specified diagonal.

    Parameters
    ----------
    m : integer
        m > 0
    n : integer
        n > 0
    g : integer
        g is greatest common divisor of m and n
    h : integer
        0 <= h <= lcm(m,n)

    Returns
    -------
    p : float
        The number of paths that go low.
        The calculation may overflow - check for a finite answer.

    Exceptions
    ----------
    FloatingPointError: Raised if the intermediate computation goes outside
    the range of a float.

    Notes
    -----
    Count the integer lattice paths from (0, 0) to (m, n), which at some
    point (x, y) along the path, satisfy:
      m*y <= n*x - h*g
    The paths make steps of size +1 in either positive x or positive y directions.

    We generally follow Hodges' treatment of Drion/Gnedenko/Korolyuk.
    Hodges, J.L. Jr.,
    "The Significance Probability of the Smirnov Two-Sample Test,"
    Arkiv fiur Matematik, 3, No. 43 (1958), 469-86.
    """
    # Compute #paths which stay lower than x/m-y/n = h/lcm(m,n)
    # B(x, y) = #{paths from (0,0) to (x,y) without previously crossing the boundary}
    #         = binom(x, y) - #{paths which already reached the boundary}
    # Multiply by the number of path extensions going from (x, y) to (m, n)
    # Sum.

    # Probability is symmetrical in m, n.  Computation below assumes m >= n.
    if m < n:
        m, n = n, m
    mg = m // g
    ng = n // g

    #  0 <= x_j <= m is the smallest integer for which n*x_j - m*j < g*h
    xj = [int(np.ceil((h + mg * j)/ng)) for j in range(n+1)]
    xj = [_ for _ in xj if _ <= m]
    lxj = len(xj)
    # B is an array just holding a few values of B(x,y), the ones needed.
    # B[j] == B(x_j, j)
    if lxj == 0:
        return np.round(special.binom(m + n, n))
    B = np.zeros(lxj)
    B[0] = 1
    # Compute the B(x, y) terms
    # The binomial coefficient is an integer, but special.binom() may return a float.
    # Round it to the nearest integer.
    for j in range(1, lxj):
        Bj = np.round(special.binom(xj[j] + j, j))
        if not np.isfinite(Bj):
            raise FloatingPointError()
        for i in range(j):
            bin = np.round(special.binom(xj[j] - xj[i] + j - i, j-i))
            dec = bin * B[i]
            Bj -= dec
        B[j] = Bj
        if not np.isfinite(Bj):
            raise FloatingPointError()
    # Compute the number of path extensions...
    num_paths = 0
    for j in range(lxj):
        bin = np.round(special.binom((m-xj[j]) + (n - j), n-j))
        term = B[j] * bin
        if not np.isfinite(term):
            raise FloatingPointError()
        num_paths += term
    return np.round(num_paths)


def ks_2samp(data1, data2, alternative='two-sided', mode='auto'):
    """
    Compute the Kolmogorov-Smirnov statistic on 2 samples.

    This is a two-sided test for the null hypothesis that 2 independent samples
    are drawn from the same continuous distribution.  The
    alternative hypothesis can be either 'two-sided' (default), 'less'
    or 'greater'.

    Parameters
    ----------
    data1, data2 : sequence of 1-D ndarrays
        Two arrays of sample observations assumed to be drawn from a continuous
        distribution, sample sizes can be different.
    alternative : {'two-sided', 'less', 'greater'}, optional
        Defines the alternative hypothesis (see explanation above).
        Default is 'two-sided'.
    mode : {'auto', 'exact', 'asymp'}, optional
        Defines the method used for calculating the p-value.
        Default is 'auto'.

        - 'exact' : use approximation to exact distribution of test statistic
        - 'asymp' : use asymptotic distribution of test statistic
        - 'auto' : use 'exact' for small size arrays, 'asymp' for large.

    Returns
    -------
    statistic : float
        KS statistic
    pvalue : float
        two-tailed p-value

    Notes
    -----
    This tests whether 2 samples are drawn from the same distribution. Note
    that, like in the case of the one-sample K-S test, the distribution is
    assumed to be continuous.

    In the one-sided test, the alternative is that the empirical
    cumulative distribution function F(x) of the data1 variable is "less"
    or "greater" than the empirical cumulative distribution function G(x)
    of the data2 variable, ``F(x)<=G(x)``, resp. ``F(x)>=G(x)``.

    If the K-S statistic is small or the p-value is high, then we cannot
    reject the hypothesis that the distributions of the two samples
    are the same.

    If the mode is 'auto', the computation is exact if the sample sizes are
    less than 10000.  For larger sizes, the computation uses the
    Kolmogorov-Smirnov distributions to compute an approximate value.

    We generally follow Hodges' treatment of Drion/Gnedenko/Korolyuk [1]_.

    References
    ----------
    .. [1] Hodges, J.L. Jr.,  "The Significance Probability of the Smirnov
           Two-Sample Test," Arkiv fiur Matematik, 3, No. 43 (1958), 469-86.


    Examples
    --------
    >>> from scipy import stats
    >>> np.random.seed(12345678)  #fix random seed to get the same result
    >>> n1 = 200  # size of first sample
    >>> n2 = 300  # size of second sample

    For a different distribution, we can reject the null hypothesis since the
    pvalue is below 1%:

    >>> rvs1 = stats.norm.rvs(size=n1, loc=0., scale=1)
    >>> rvs2 = stats.norm.rvs(size=n2, loc=0.5, scale=1.5)
    >>> stats.ks_2samp(rvs1, rvs2)
    (0.20833333333333334, 5.129279597781977e-05)

    For a slightly different distribution, we cannot reject the null hypothesis
    at a 10% or lower alpha since the p-value at 0.144 is higher than 10%

    >>> rvs3 = stats.norm.rvs(size=n2, loc=0.01, scale=1.0)
    >>> stats.ks_2samp(rvs1, rvs3)
    (0.10333333333333333, 0.14691437867433876)

    For an identical distribution, we cannot reject the null hypothesis since
    the p-value is high, 41%:

    >>> rvs4 = stats.norm.rvs(size=n2, loc=0.0, scale=1.0)
    >>> stats.ks_2samp(rvs1, rvs4)
    (0.07999999999999996, 0.41126949729859719)

    """
    LARGE_N = 10000  # 'auto' will attempt to be exact if n1,n2 <= LARGE_N
    data1 = np.sort(data1)
    data2 = np.sort(data2)
    n1 = data1.shape[0]
    n2 = data2.shape[0]
    if min(n1, n2) == 0:
        raise ValueError('Data passed to ks_2samp must not be empty')

    data_all = np.concatenate([data1, data2])
    # using searchsorted solves equal data problem
    cdf1 = np.searchsorted(data1, data_all, side='right') / n1
    cdf2 = np.searchsorted(data2, data_all, side='right') / n2
    cddiffs = cdf1 - cdf2
    minS = -np.min(cddiffs)
    maxS = np.max(cddiffs)
    alt2Dvalue = {'less': minS, 'greater': maxS, 'two-sided': max(minS, maxS)}
    d = alt2Dvalue[alternative]
    g = gcd(n1, n2)
    n1g = n1 // g
    n2g = n2 // g
    prob = -np.inf
    original_mode = mode
    if mode == 'auto':
        if max(n1, n2) <= LARGE_N:
            mode = 'exact'
        else:
            mode = 'asymp'
    elif mode == 'exact':
        # If lcm(n1, n2) is too big, switch from exact to asymp
        if n1g >= np.iinfo(np.int).max / n2g:
            mode = 'asymp'
            warnings.warn(
                "Exact ks_2samp calculation not possible with samples sizes "
                "%d and %d. Switching to 'asymp' " % (n1, n2), RuntimeWarning)

    saw_fp_error = False
    if mode == 'exact':
        lcm = (n1 // g) * n2
        h = int(np.round(d * lcm))
        d = h * 1.0 / lcm
        if h == 0:
            prob = 1.0
        else:
            try:
                if alternative == 'two-sided':
                    if n1 == n2:
                        prob = _compute_prob_outside_square(n1, h)
                    else:
                        prob = 1 - _compute_prob_inside_method(n1, n2, g, h)
                else:
                    if n1 == n2:
                        # prob = binom(2n, n-h) / binom(2n, n)
                        # Evaluating in that form incurs roundoff errors
                        # from special.binom. Instead calculate directly
                        prob = 1.0
                        for j in range(h):
                            prob = (n1 - j) * prob / (n1 + j + 1)
                    else:
                        num_paths = _count_paths_outside_method(n1, n2, g, h)
                        bin = special.binom(n1 + n2, n1)
                        if not np.isfinite(bin) or not np.isfinite(num_paths) or num_paths > bin:
                            raise FloatingPointError()
                        prob = num_paths / bin

            except FloatingPointError:
                # Switch mode
                mode = 'asymp'
                saw_fp_error = True
                # Can't raise warning here, inside the try
            finally:
                if saw_fp_error:
                    if original_mode == 'exact':
                        warnings.warn(
                            "ks_2samp: Exact calculation overflowed. "
                            "Switching to mode=%s" % mode, RuntimeWarning)
                else:
                    if prob > 1 or prob < 0:
                        mode = 'asymp'
                        if original_mode == 'exact':
                            warnings.warn(
                                "ks_2samp: Exact calculation incurred large"
                                " rounding error. Switching to mode=%s" % mode,
                                RuntimeWarning)

    if mode == 'asymp':
        # The product n1*n2 is large.  Use Smirnov's asymptoptic formula.
        if alternative == 'two-sided':
            en = np.sqrt(n1 * n2 / (n1 + n2))
            # Switch to using kstwo.sf() when it becomes available.
            # prob = distributions.kstwo.sf(d, int(np.round(en)))
            prob = distributions.kstwobign.sf(en * d)
        else:
            m, n = max(n1, n2), min(n1, n2)
            z = np.sqrt(m*n/(m+n)) * d
            # Use Hodges' suggested approximation Eqn 5.3
            expt = -2 * z**2 - 2 * z * (m + 2*n)/np.sqrt(m*n*(m+n))/3.0
            prob = np.exp(expt)

    prob = (0 if prob < 0 else (1 if prob > 1 else prob))
    return Ks_2sampResult(d, prob)


def tiecorrect(rankvals):
    """
    Tie correction factor for ties in the Mann-Whitney U and
    Kruskal-Wallis H tests.

    Parameters
    ----------
    rankvals : array_like
        A 1-D sequence of ranks.  Typically this will be the array
        returned by `~scipy.stats.rankdata`.

    Returns
    -------
    factor : float
        Correction factor for U or H.

    See Also
    --------
    rankdata : Assign ranks to the data
    mannwhitneyu : Mann-Whitney rank test
    kruskal : Kruskal-Wallis H test

    References
    ----------
    .. [1] Siegel, S. (1956) Nonparametric Statistics for the Behavioral
           Sciences.  New York: McGraw-Hill.

    Examples
    --------
    >>> from scipy.stats import tiecorrect, rankdata
    >>> tiecorrect([1, 2.5, 2.5, 4])
    0.9
    >>> ranks = rankdata([1, 3, 2, 4, 5, 7, 2, 8, 4])
    >>> ranks
    array([ 1. ,  4. ,  2.5,  5.5,  7. ,  8. ,  2.5,  9. ,  5.5])
    >>> tiecorrect(ranks)
    0.9833333333333333

    """
    arr = np.sort(rankvals)
    idx = np.nonzero(np.r_[True, arr[1:] != arr[:-1], True])[0]
    cnt = np.diff(idx).astype(np.float64)

    size = np.float64(arr.size)
    return 1.0 if size < 2 else 1.0 - (cnt**3 - cnt).sum() / (size**3 - size)


MannwhitneyuResult = namedtuple('MannwhitneyuResult', ('statistic', 'pvalue'))

def mannwhitneyu(x, y, use_continuity=True, alternative=None):
    """
    Compute the Mann-Whitney rank test on samples x and y.

    Parameters
    ----------
    x, y : array_like
        Array of samples, should be one-dimensional.
    use_continuity : bool, optional
            Whether a continuity correction (1/2.) should be taken into
            account. Default is True.
    alternative : None (deprecated), 'less', 'two-sided', or 'greater'
            Whether to get the p-value for the one-sided hypothesis ('less'
            or 'greater') or for the two-sided hypothesis ('two-sided').
            Defaults to None, which results in a p-value half the size of
            the 'two-sided' p-value and a different U statistic. The
            default behavior is not the same as using 'less' or 'greater':
            it only exists for backward compatibility and is deprecated.

    Returns
    -------
    statistic : float
        The Mann-Whitney U statistic, equal to min(U for x, U for y) if
        `alternative` is equal to None (deprecated; exists for backward
        compatibility), and U for y otherwise.
    pvalue : float
        p-value assuming an asymptotic normal distribution. One-sided or
        two-sided, depending on the choice of `alternative`.

    Notes
    -----
    Use only when the number of observation in each sample is > 20 and
    you have 2 independent samples of ranks. Mann-Whitney U is
    significant if the u-obtained is LESS THAN or equal to the critical
    value of U.

    This test corrects for ties and by default uses a continuity correction.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Mann-Whitney_U_test

    .. [2] H.B. Mann and D.R. Whitney, "On a Test of Whether one of Two Random
           Variables is Stochastically Larger than the Other," The Annals of
           Mathematical Statistics, vol. 18, no. 1, pp. 50-60, 1947.

    """
    if alternative is None:
        warnings.warn("Calling `mannwhitneyu` without specifying "
                      "`alternative` is deprecated.", DeprecationWarning)

    x = np.asarray(x)
    y = np.asarray(y)
    n1 = len(x)
    n2 = len(y)
    ranked = rankdata(np.concatenate((x, y)))
    rankx = ranked[0:n1]  # get the x-ranks
    u1 = n1*n2 + (n1*(n1+1))/2.0 - np.sum(rankx, axis=0)  # calc U for x
    u2 = n1*n2 - u1  # remainder is U for y
    T = tiecorrect(ranked)
    if T == 0:
        raise ValueError('All numbers are identical in mannwhitneyu')
    sd = np.sqrt(T * n1 * n2 * (n1+n2+1) / 12.0)

    meanrank = n1*n2/2.0 + 0.5 * use_continuity
    if alternative is None or alternative == 'two-sided':
        bigu = max(u1, u2)
    elif alternative == 'less':
        bigu = u1
    elif alternative == 'greater':
        bigu = u2
    else:
        raise ValueError("alternative should be None, 'less', 'greater' "
                         "or 'two-sided'")

    z = (bigu - meanrank) / sd
    if alternative is None:
        # This behavior, equal to half the size of the two-sided
        # p-value, is deprecated.
        p = distributions.norm.sf(abs(z))
    elif alternative == 'two-sided':
        p = 2 * distributions.norm.sf(abs(z))
    else:
        p = distributions.norm.sf(z)

    u = u2
    # This behavior is deprecated.
    if alternative is None:
        u = min(u1, u2)
    return MannwhitneyuResult(u, p)


RanksumsResult = namedtuple('RanksumsResult', ('statistic', 'pvalue'))


def ranksums(x, y):
    """
    Compute the Wilcoxon rank-sum statistic for two samples.

    The Wilcoxon rank-sum test tests the null hypothesis that two sets
    of measurements are drawn from the same distribution.  The alternative
    hypothesis is that values in one sample are more likely to be
    larger than the values in the other sample.

    This test should be used to compare two samples from continuous
    distributions.  It does not handle ties between measurements
    in x and y.  For tie-handling and an optional continuity correction
    see `scipy.stats.mannwhitneyu`.

    Parameters
    ----------
    x,y : array_like
        The data from the two samples

    Returns
    -------
    statistic : float
        The test statistic under the large-sample approximation that the
        rank sum statistic is normally distributed
    pvalue : float
        The two-sided p-value of the test

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Wilcoxon_rank-sum_test

    """
    x, y = map(np.asarray, (x, y))
    n1 = len(x)
    n2 = len(y)
    alldata = np.concatenate((x, y))
    ranked = rankdata(alldata)
    x = ranked[:n1]
    s = np.sum(x, axis=0)
    expected = n1 * (n1+n2+1) / 2.0
    z = (s - expected) / np.sqrt(n1*n2*(n1+n2+1)/12.0)
    prob = 2 * distributions.norm.sf(abs(z))

    return RanksumsResult(z, prob)


KruskalResult = namedtuple('KruskalResult', ('statistic', 'pvalue'))


def kruskal(*args, **kwargs):
    """
    Compute the Kruskal-Wallis H-test for independent samples

    The Kruskal-Wallis H-test tests the null hypothesis that the population
    median of all of the groups are equal.  It is a non-parametric version of
    ANOVA.  The test works on 2 or more independent samples, which may have
    different sizes.  Note that rejecting the null hypothesis does not
    indicate which of the groups differs.  Post-hoc comparisons between
    groups are required to determine which groups are different.

    Parameters
    ----------
    sample1, sample2, ... : array_like
       Two or more arrays with the sample measurements can be given as
       arguments.
    nan_policy : {'propagate', 'raise', 'omit'}, optional
        Defines how to handle when input contains nan. 'propagate' returns nan,
        'raise' throws an error, 'omit' performs the calculations ignoring nan
        values. Default is 'propagate'.

    Returns
    -------
    statistic : float
       The Kruskal-Wallis H statistic, corrected for ties
    pvalue : float
       The p-value for the test using the assumption that H has a chi
       square distribution

    See Also
    --------
    f_oneway : 1-way ANOVA
    mannwhitneyu : Mann-Whitney rank test on two samples.
    friedmanchisquare : Friedman test for repeated measurements

    Notes
    -----
    Due to the assumption that H has a chi square distribution, the number
    of samples in each group must not be too small.  A typical rule is
    that each sample must have at least 5 measurements.

    References
    ----------
    .. [1] W. H. Kruskal & W. W. Wallis, "Use of Ranks in
       One-Criterion Variance Analysis", Journal of the American Statistical
       Association, Vol. 47, Issue 260, pp. 583-621, 1952.
    .. [2] https://en.wikipedia.org/wiki/Kruskal-Wallis_one-way_analysis_of_variance

    Examples
    --------
    >>> from scipy import stats
    >>> x = [1, 3, 5, 7, 9]
    >>> y = [2, 4, 6, 8, 10]
    >>> stats.kruskal(x, y)
    KruskalResult(statistic=0.2727272727272734, pvalue=0.6015081344405895)

    >>> x = [1, 1, 1]
    >>> y = [2, 2, 2]
    >>> z = [2, 2]
    >>> stats.kruskal(x, y, z)
    KruskalResult(statistic=7.0, pvalue=0.0301973834223185)

    """
    args = list(map(np.asarray, args))
    num_groups = len(args)
    if num_groups < 2:
        raise ValueError("Need at least two groups in stats.kruskal()")

    for arg in args:
        if arg.size == 0:
            return KruskalResult(np.nan, np.nan)
    n = np.asarray(list(map(len, args)))

    if 'nan_policy' in kwargs.keys():
        if kwargs['nan_policy'] not in ('propagate', 'raise', 'omit'):
            raise ValueError("nan_policy must be 'propagate', "
                             "'raise' or'omit'")
        else:
            nan_policy = kwargs['nan_policy']
    else:
        nan_policy = 'propagate'

    contains_nan = False
    for arg in args:
        cn = _contains_nan(arg, nan_policy)
        if cn[0]:
            contains_nan = True
            break

    if contains_nan and nan_policy == 'omit':
        for a in args:
            a = ma.masked_invalid(a)
        return mstats_basic.kruskal(*args)

    if contains_nan and nan_policy == 'propagate':
        return KruskalResult(np.nan, np.nan)

    alldata = np.concatenate(args)
    ranked = rankdata(alldata)
    ties = tiecorrect(ranked)
    if ties == 0:
        raise ValueError('All numbers are identical in kruskal')

    # Compute sum^2/n for each group and sum
    j = np.insert(np.cumsum(n), 0, 0)
    ssbn = 0
    for i in range(num_groups):
        ssbn += _square_of_sums(ranked[j[i]:j[i+1]]) / n[i]

    totaln = np.sum(n)
    h = 12.0 / (totaln * (totaln + 1)) * ssbn - 3 * (totaln + 1)
    df = num_groups - 1
    h /= ties

    return KruskalResult(h, distributions.chi2.sf(h, df))


FriedmanchisquareResult = namedtuple('FriedmanchisquareResult',
                                     ('statistic', 'pvalue'))


def friedmanchisquare(*args):
    """
    Compute the Friedman test for repeated measurements

    The Friedman test tests the null hypothesis that repeated measurements of
    the same individuals have the same distribution.  It is often used
    to test for consistency among measurements obtained in different ways.
    For example, if two measurement techniques are used on the same set of
    individuals, the Friedman test can be used to determine if the two
    measurement techniques are consistent.

    Parameters
    ----------
    measurements1, measurements2, measurements3... : array_like
        Arrays of measurements.  All of the arrays must have the same number
        of elements.  At least 3 sets of measurements must be given.

    Returns
    -------
    statistic : float
        the test statistic, correcting for ties
    pvalue : float
        the associated p-value assuming that the test statistic has a chi
        squared distribution

    Notes
    -----
    Due to the assumption that the test statistic has a chi squared
    distribution, the p-value is only reliable for n > 10 and more than
    6 repeated measurements.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Friedman_test

    """
    k = len(args)
    if k < 3:
        raise ValueError('Less than 3 levels.  Friedman test not appropriate.')

    n = len(args[0])
    for i in range(1, k):
        if len(args[i]) != n:
            raise ValueError('Unequal N in friedmanchisquare.  Aborting.')

    # Rank data
    data = np.vstack(args).T
    data = data.astype(float)
    for i in range(len(data)):
        data[i] = rankdata(data[i])

    # Handle ties
    ties = 0
    for i in range(len(data)):
        replist, repnum = find_repeats(array(data[i]))
        for t in repnum:
            ties += t * (t*t - 1)
    c = 1 - ties / (k*(k*k - 1)*n)

    ssbn = np.sum(data.sum(axis=0)**2)
    chisq = (12.0 / (k*n*(k+1)) * ssbn - 3*n*(k+1)) / c

    return FriedmanchisquareResult(chisq, distributions.chi2.sf(chisq, k - 1))


BrunnerMunzelResult = namedtuple('BrunnerMunzelResult',
                                 ('statistic', 'pvalue'))


def brunnermunzel(x, y, alternative="two-sided", distribution="t",
                  nan_policy='propagate'):
    """
    Computes the Brunner-Munzel test on samples x and y

    The Brunner-Munzel test is a nonparametric test of the null hypothesis that
    when values are taken one by one from each group, the probabilities of
    getting large values in both groups are equal.
    Unlike the Wilcoxon-Mann-Whitney's U test, this does not require the
    assumption of equivariance of two groups. Note that this does not assume
    the distributions are same. This test works on two independent samples,
    which may have different sizes.

    Parameters
    ----------
    x, y : array_like
        Array of samples, should be one-dimensional.
    alternative :  'less', 'two-sided', or 'greater', optional
        Whether to get the p-value for the one-sided hypothesis ('less'
        or 'greater') or for the two-sided hypothesis ('two-sided').
        Defaults value is 'two-sided' .
    distribution: 't' or 'normal', optional
        Whether to get the p-value by t-distribution or by standard normal
        distribution.
        Defaults value is 't' .
    nan_policy : {'propagate', 'raise', 'omit'}, optional
        Defines how to handle when input contains nan. 'propagate' returns nan,
        'raise' throws an error, 'omit' performs the calculations ignoring nan
        values. Default is 'propagate'.

    Returns
    -------
    statistic : float
        The Brunner-Munzer W statistic.
    pvalue : float
        p-value assuming an t distribution. One-sided or
        two-sided, depending on the choice of `alternative` and `distribution`.

    See Also
    --------
    mannwhitneyu : Mann-Whitney rank test on two samples.

    Notes
    -------
    Brunner and Munzel recommended to estimate the p-value by t-distribution
    when the size of data is 50 or less. If the size is lower than 10, it would
    be better to use permuted Brunner Munzel test (see [2]_).

    References
    ----------
    .. [1] Brunner, E. and Munzel, U. "The nonparametric Benhrens-Fisher
           problem: Asymptotic theory and a small-sample approximation".
           Biometrical Journal. Vol. 42(2000): 17-25.
    .. [2] Neubert, K. and Brunner, E. "A studentized permutation test for the
           non-parametric Behrens-Fisher problem". Computational Statistics and
           Data Analysis. Vol. 51(2007): 5192-5204.

    Examples
    --------
    >>> from scipy import stats
    >>> x1 = [1,2,1,1,1,1,1,1,1,1,2,4,1,1]
    >>> x2 = [3,3,4,3,1,2,3,1,1,5,4]
    >>> w, p_value = stats.brunnermunzel(x1, x2)
    >>> w
    3.1374674823029505
    >>> p_value
    0.0057862086661515377

    """
    x = np.asarray(x)
    y = np.asarray(y)

    # check both x and y
    cnx, npx = _contains_nan(x, nan_policy)
    cny, npy = _contains_nan(y, nan_policy)
    contains_nan = cnx or cny
    if npx == "omit" or npy == "omit":
        nan_policy = "omit"

    if contains_nan and nan_policy == "propagate":
        return BrunnerMunzelResult(np.nan, np.nan)
    elif contains_nan and nan_policy == "omit":
        x = ma.masked_invalid(x)
        y = ma.masked_invalid(y)
        return mstats_basic.brunnermunzel(x, y, alternative, distribution)

    nx = len(x)
    ny = len(y)
    if nx == 0 or ny == 0:
        return BrunnerMunzelResult(np.nan, np.nan)
    rankc = rankdata(np.concatenate((x, y)))
    rankcx = rankc[0:nx]
    rankcy = rankc[nx:nx+ny]
    rankcx_mean = np.mean(rankcx)
    rankcy_mean = np.mean(rankcy)
    rankx = rankdata(x)
    ranky = rankdata(y)
    rankx_mean = np.mean(rankx)
    ranky_mean = np.mean(ranky)

    Sx = np.sum(np.power(rankcx - rankx - rankcx_mean + rankx_mean, 2.0))
    Sx /= nx - 1
    Sy = np.sum(np.power(rankcy - ranky - rankcy_mean + ranky_mean, 2.0))
    Sy /= ny - 1

    wbfn = nx * ny * (rankcy_mean - rankcx_mean)
    wbfn /= (nx + ny) * np.sqrt(nx * Sx + ny * Sy)

    if distribution == "t":
        df_numer = np.power(nx * Sx + ny * Sy, 2.0)
        df_denom = np.power(nx * Sx, 2.0) / (nx - 1)
        df_denom += np.power(ny * Sy, 2.0) / (ny - 1)
        df = df_numer / df_denom
        p = distributions.t.cdf(wbfn, df)
    elif distribution == "normal":
        p = distributions.norm.cdf(wbfn)
    else:
        raise ValueError(
            "distribution should be 't' or 'normal'")

    if alternative == "greater":
        p = p
    elif alternative == "less":
        p = 1 - p
    elif alternative == "two-sided":
        p = 2 * np.min([p, 1-p])
    else:
        raise ValueError(
            "alternative should be 'less', 'greater' or 'two-sided'")

    return BrunnerMunzelResult(wbfn, p)


def combine_pvalues(pvalues, method='fisher', weights=None):
    """
    Methods for combining the p-values of independent tests bearing upon the
    same hypothesis.

    Parameters
    ----------
    pvalues : array_like, 1-D
        Array of p-values assumed to come from independent tests.
    method : {'fisher', 'pearson', 'tippett', 'stouffer', 'mudholkar_george'},
    optional.
        Name of method to use to combine p-values. The following methods are
        available:

        - "fisher": Fisher's method (Fisher's combined probability test),
          the default, the sum of the logarithm of the p-values.
        - "pearson": Pearson's method (similar to Fisher's but uses sum of the
          complement of the p-values inside the logarithms).
        - "tippett": Tippett's method (minimum of p-values).
        - "stouffer": Stouffer's Z-score method.
        - "mudholkar_george": the difference of Fisher's and Pearson's methods
           divided by 2.
    weights : array_like, 1-D, optional
        Optional array of weights used only for Stouffer's Z-score method.

    Returns
    -------
    statistic: float
        The statistic calculated by the specified method.
    pval: float
        The combined p-value.

    Notes
    -----
    Fisher's method (also known as Fisher's combined probability test) [1]_ uses
    a chi-squared statistic to compute a combined p-value. The closely related
    Stouffer's Z-score method [2]_ uses Z-scores rather than p-values. The
    advantage of Stouffer's method is that it is straightforward to introduce
    weights, which can make Stouffer's method more powerful than Fisher's
    method when the p-values are from studies of different size [6]_ [7]_.
    The Pearson's method uses :math:`log(1-p_i)` inside the sum whereas Fisher's
    method uses :math:`log(p_i)` [4]_. For Fisher's and Pearson's method, the
    sum of the logarithms is multiplied by -2 in the implementation. This
    quantity has a chisquare distribution that determines the p-value. The
    `mudholkar_george` method is the difference of the Fisher's and Pearson's
    test statistics, each of which include the -2 factor [4]_. However, the
    `mudholkar_george` method does not include these -2 factors. The test
    statistic of `mudholkar_george` is the sum of logisitic random variables and
    equation 3.6 in [3]_ is used to approximate the p-value based on Student's
    t-distribution.

    Fisher's method may be extended to combine p-values from dependent tests
    [5]_. Extensions such as Brown's method and Kost's method are not currently
    implemented.

    .. versionadded:: 0.15.0

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Fisher%27s_method
    .. [2] https://en.wikipedia.org/wiki/Fisher%27s_method#Relation_to_Stouffer.27s_Z-score_method
    .. [3] George, E. O., and G. S. Mudholkar. "On the convolution of logistic
           random variables." Metrika 30.1 (1983): 1-13.
    .. [4] Heard, N. and Rubin-Delanchey, P. "Choosing between methods of
           combining p-values."  Biometrika 105.1 (2018): 239-246.
    .. [5] Whitlock, M. C. "Combining probability from independent tests: the
           weighted Z-method is superior to Fisher's approach." Journal of
           Evolutionary Biology 18, no. 5 (2005): 1368-1373.
    .. [6] Zaykin, Dmitri V. "Optimally weighted Z-test is a powerful method
           for combining probabilities in meta-analysis." Journal of
           Evolutionary Biology 24, no. 8 (2011): 1836-1841.
    .. [7] https://en.wikipedia.org/wiki/Extensions_of_Fisher%27s_method

    """
    pvalues = np.asarray(pvalues)
    if pvalues.ndim != 1:
        raise ValueError("pvalues is not 1-D")

    if method == 'fisher':
        statistic = -2 * np.sum(np.log(pvalues))
        pval = distributions.chi2.sf(statistic, 2 * len(pvalues))
    elif method == 'pearson':
        statistic = -2 * np.sum(np.log1p(-pvalues))
        pval = distributions.chi2.sf(statistic, 2 * len(pvalues))
    elif method == 'mudholkar_george':
        statistic = -np.sum(np.log(pvalues)) + np.sum(np.log1p(-pvalues))
        nu = 5 * len(pvalues) + 4
        approx_factor = np.sqrt(nu / (nu - 2))
        pval = distributions.t.sf(statistic * approx_factor, nu)
    elif method == 'tippett':
        statistic = np.min(pvalues)
        pval = distributions.beta.sf(statistic, 1, len(pvalues))
    elif method == 'stouffer':
        if weights is None:
            weights = np.ones_like(pvalues)
        elif len(weights) != len(pvalues):
            raise ValueError("pvalues and weights must be of the same size.")

        weights = np.asarray(weights)
        if weights.ndim != 1:
            raise ValueError("weights is not 1-D")

        Zi = distributions.norm.isf(pvalues)
        statistic = np.dot(weights, Zi) / np.linalg.norm(weights)
        pval = distributions.norm.sf(statistic)

    else:
        raise ValueError(
            "Invalid method '%s'. Options are 'fisher', 'pearson', \
            'mudholkar_george', 'tippett', 'or 'stouffer'", method)

    return (statistic, pval)


#####################################
#       STATISTICAL DISTANCES       #
#####################################

def wasserstein_distance(u_values, v_values, u_weights=None, v_weights=None):
    r"""
    Compute the first Wasserstein distance between two 1D distributions.

    This distance is also known as the earth mover's distance, since it can be
    seen as the minimum amount of "work" required to transform :math:`u` into
    :math:`v`, where "work" is measured as the amount of distribution weight
    that must be moved, multiplied by the distance it has to be moved.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    u_values, v_values : array_like
        Values observed in the (empirical) distribution.
    u_weights, v_weights : array_like, optional
        Weight for each value. If unspecified, each value is assigned the same
        weight.
        `u_weights` (resp. `v_weights`) must have the same length as
        `u_values` (resp. `v_values`). If the weight sum differs from 1, it
        must still be positive and finite so that the weights can be normalized
        to sum to 1.

    Returns
    -------
    distance : float
        The computed distance between the distributions.

    Notes
    -----
    The first Wasserstein distance between the distributions :math:`u` and
    :math:`v` is:

    .. math::

        l_1 (u, v) = \inf_{\pi \in \Gamma (u, v)} \int_{\mathbb{R} \times
        \mathbb{R}} |x-y| \mathrm{d} \pi (x, y)

    where :math:`\Gamma (u, v)` is the set of (probability) distributions on
    :math:`\mathbb{R} \times \mathbb{R}` whose marginals are :math:`u` and
    :math:`v` on the first and second factors respectively.

    If :math:`U` and :math:`V` are the respective CDFs of :math:`u` and
    :math:`v`, this distance also equals to:

    .. math::

        l_1(u, v) = \int_{-\infty}^{+\infty} |U-V|

    See [2]_ for a proof of the equivalence of both definitions.

    The input distributions can be empirical, therefore coming from samples
    whose values are effectively inputs of the function, or they can be seen as
    generalized functions, in which case they are weighted sums of Dirac delta
    functions located at the specified values.

    References
    ----------
    .. [1] "Wasserstein metric", https://en.wikipedia.org/wiki/Wasserstein_metric
    .. [2] Ramdas, Garcia, Cuturi "On Wasserstein Two Sample Testing and Related
           Families of Nonparametric Tests" (2015). :arXiv:`1509.02237`.

    Examples
    --------
    >>> from scipy.stats import wasserstein_distance
    >>> wasserstein_distance([0, 1, 3], [5, 6, 8])
    5.0
    >>> wasserstein_distance([0, 1], [0, 1], [3, 1], [2, 2])
    0.25
    >>> wasserstein_distance([3.4, 3.9, 7.5, 7.8], [4.5, 1.4],
    ...                      [1.4, 0.9, 3.1, 7.2], [3.2, 3.5])
    4.0781331438047861
    """
    return _cdf_distance(1, u_values, v_values, u_weights, v_weights)


def energy_distance(u_values, v_values, u_weights=None, v_weights=None):
    r"""
    Compute the energy distance between two 1D distributions.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    u_values, v_values : array_like
        Values observed in the (empirical) distribution.
    u_weights, v_weights : array_like, optional
        Weight for each value. If unspecified, each value is assigned the same
        weight.
        `u_weights` (resp. `v_weights`) must have the same length as
        `u_values` (resp. `v_values`). If the weight sum differs from 1, it
        must still be positive and finite so that the weights can be normalized
        to sum to 1.

    Returns
    -------
    distance : float
        The computed distance between the distributions.

    Notes
    -----
    The energy distance between two distributions :math:`u` and :math:`v`, whose
    respective CDFs are :math:`U` and :math:`V`, equals to:

    .. math::

        D(u, v) = \left( 2\mathbb E|X - Y| - \mathbb E|X - X'| -
        \mathbb E|Y - Y'| \right)^{1/2}

    where :math:`X` and :math:`X'` (resp. :math:`Y` and :math:`Y'`) are
    independent random variables whose probability distribution is :math:`u`
    (resp. :math:`v`).

    As shown in [2]_, for one-dimensional real-valued variables, the energy
    distance is linked to the non-distribution-free version of the Cramer-von
    Mises distance:

    .. math::

        D(u, v) = \sqrt{2} l_2(u, v) = \left( 2 \int_{-\infty}^{+\infty} (U-V)^2
        \right)^{1/2}

    Note that the common Cramer-von Mises criterion uses the distribution-free
    version of the distance. See [2]_ (section 2), for more details about both
    versions of the distance.

    The input distributions can be empirical, therefore coming from samples
    whose values are effectively inputs of the function, or they can be seen as
    generalized functions, in which case they are weighted sums of Dirac delta
    functions located at the specified values.

    References
    ----------
    .. [1] "Energy distance", https://en.wikipedia.org/wiki/Energy_distance
    .. [2] Szekely "E-statistics: The energy of statistical samples." Bowling
           Green State University, Department of Mathematics and Statistics,
           Technical Report 02-16 (2002).
    .. [3] Rizzo, Szekely "Energy distance." Wiley Interdisciplinary Reviews:
           Computational Statistics, 8(1):27-38 (2015).
    .. [4] Bellemare, Danihelka, Dabney, Mohamed, Lakshminarayanan, Hoyer,
           Munos "The Cramer Distance as a Solution to Biased Wasserstein
           Gradients" (2017). :arXiv:`1705.10743`.

    Examples
    --------
    >>> from scipy.stats import energy_distance
    >>> energy_distance([0], [2])
    2.0000000000000004
    >>> energy_distance([0, 8], [0, 8], [3, 1], [2, 2])
    1.0000000000000002
    >>> energy_distance([0.7, 7.4, 2.4, 6.8], [1.4, 8. ],
    ...                 [2.1, 4.2, 7.4, 8. ], [7.6, 8.8])
    0.88003340976158217
    """
    return np.sqrt(2) * _cdf_distance(2, u_values, v_values,
                                      u_weights, v_weights)


def _cdf_distance(p, u_values, v_values, u_weights=None, v_weights=None):
    r"""
    Compute, between two one-dimensional distributions :math:`u` and
    :math:`v`, whose respective CDFs are :math:`U` and :math:`V`, the
    statistical distance that is defined as:

    .. math::

        l_p(u, v) = \left( \int_{-\infty}^{+\infty} |U-V|^p \right)^{1/p}

    p is a positive parameter; p = 1 gives the Wasserstein distance, p = 2
    gives the energy distance.

    Parameters
    ----------
    u_values, v_values : array_like
        Values observed in the (empirical) distribution.
    u_weights, v_weights : array_like, optional
        Weight for each value. If unspecified, each value is assigned the same
        weight.
        `u_weights` (resp. `v_weights`) must have the same length as
        `u_values` (resp. `v_values`). If the weight sum differs from 1, it
        must still be positive and finite so that the weights can be normalized
        to sum to 1.

    Returns
    -------
    distance : float
        The computed distance between the distributions.

    Notes
    -----
    The input distributions can be empirical, therefore coming from samples
    whose values are effectively inputs of the function, or they can be seen as
    generalized functions, in which case they are weighted sums of Dirac delta
    functions located at the specified values.

    References
    ----------
    .. [1] Bellemare, Danihelka, Dabney, Mohamed, Lakshminarayanan, Hoyer,
           Munos "The Cramer Distance as a Solution to Biased Wasserstein
           Gradients" (2017). :arXiv:`1705.10743`.
    """
    u_values, u_weights = _validate_distribution(u_values, u_weights)
    v_values, v_weights = _validate_distribution(v_values, v_weights)

    u_sorter = np.argsort(u_values)
    v_sorter = np.argsort(v_values)

    all_values = np.concatenate((u_values, v_values))
    all_values.sort(kind='mergesort')

    # Compute the differences between pairs of successive values of u and v.
    deltas = np.diff(all_values)

    # Get the respective positions of the values of u and v among the values of
    # both distributions.
    u_cdf_indices = u_values[u_sorter].searchsorted(all_values[:-1], 'right')
    v_cdf_indices = v_values[v_sorter].searchsorted(all_values[:-1], 'right')

    # Calculate the CDFs of u and v using their weights, if specified.
    if u_weights is None:
        u_cdf = u_cdf_indices / u_values.size
    else:
        u_sorted_cumweights = np.concatenate(([0],
                                              np.cumsum(u_weights[u_sorter])))
        u_cdf = u_sorted_cumweights[u_cdf_indices] / u_sorted_cumweights[-1]

    if v_weights is None:
        v_cdf = v_cdf_indices / v_values.size
    else:
        v_sorted_cumweights = np.concatenate(([0],
                                              np.cumsum(v_weights[v_sorter])))
        v_cdf = v_sorted_cumweights[v_cdf_indices] / v_sorted_cumweights[-1]

    # Compute the value of the integral based on the CDFs.
    # If p = 1 or p = 2, we avoid using np.power, which introduces an overhead
    # of about 15%.
    if p == 1:
        return np.sum(np.multiply(np.abs(u_cdf - v_cdf), deltas))
    if p == 2:
        return np.sqrt(np.sum(np.multiply(np.square(u_cdf - v_cdf), deltas)))
    return np.power(np.sum(np.multiply(np.power(np.abs(u_cdf - v_cdf), p),
                                       deltas)), 1/p)


def _validate_distribution(values, weights):
    """
    Validate the values and weights from a distribution input of `cdf_distance`
    and return them as ndarray objects.

    Parameters
    ----------
    values : array_like
        Values observed in the (empirical) distribution.
    weights : array_like
        Weight for each value.

    Returns
    -------
    values : ndarray
        Values as ndarray.
    weights : ndarray
        Weights as ndarray.
    """
    # Validate the value array.
    values = np.asarray(values, dtype=float)
    if len(values) == 0:
        raise ValueError("Distribution can't be empty.")

    # Validate the weight array, if specified.
    if weights is not None:
        weights = np.asarray(weights, dtype=float)
        if len(weights) != len(values):
            raise ValueError('Value and weight array-likes for the same '
                             'empirical distribution must be of the same size.')
        if np.any(weights < 0):
            raise ValueError('All weights must be non-negative.')
        if not 0 < np.sum(weights) < np.inf:
            raise ValueError('Weight array-like sum must be positive and '
                             'finite. Set as None for an equal distribution of '
                             'weight.')

        return values, weights

    return values, None


#####################################
#         SUPPORT FUNCTIONS         #
#####################################

RepeatedResults = namedtuple('RepeatedResults', ('values', 'counts'))


def find_repeats(arr):
    """
    Find repeats and repeat counts.

    Parameters
    ----------
    arr : array_like
        Input array. This is cast to float64.

    Returns
    -------
    values : ndarray
        The unique values from the (flattened) input that are repeated.

    counts : ndarray
        Number of times the corresponding 'value' is repeated.

    Notes
    -----
    In numpy >= 1.9 `numpy.unique` provides similar functionality. The main
    difference is that `find_repeats` only returns repeated values.

    Examples
    --------
    >>> from scipy import stats
    >>> stats.find_repeats([2, 1, 2, 3, 2, 2, 5])
    RepeatedResults(values=array([2.]), counts=array([4]))

    >>> stats.find_repeats([[10, 20, 1, 2], [5, 5, 4, 4]])
    RepeatedResults(values=array([4.,  5.]), counts=array([2, 2]))

    """
    # Note: always copies.
    return RepeatedResults(*_find_repeats(np.array(arr, dtype=np.float64)))


def _sum_of_squares(a, axis=0):
    """
    Square each element of the input array, and return the sum(s) of that.

    Parameters
    ----------
    a : array_like
        Input array.
    axis : int or None, optional
        Axis along which to calculate. Default is 0. If None, compute over
        the whole array `a`.

    Returns
    -------
    sum_of_squares : ndarray
        The sum along the given axis for (a**2).

    See also
    --------
    _square_of_sums : The square(s) of the sum(s) (the opposite of
    `_sum_of_squares`).
    """
    a, axis = _chk_asarray(a, axis)
    return np.sum(a*a, axis)


def _square_of_sums(a, axis=0):
    """
    Sum elements of the input array, and return the square(s) of that sum.

    Parameters
    ----------
    a : array_like
        Input array.
    axis : int or None, optional
        Axis along which to calculate. Default is 0. If None, compute over
        the whole array `a`.

    Returns
    -------
    square_of_sums : float or ndarray
        The square of the sum over `axis`.

    See also
    --------
    _sum_of_squares : The sum of squares (the opposite of `square_of_sums`).
    """
    a, axis = _chk_asarray(a, axis)
    s = np.sum(a, axis)
    if not np.isscalar(s):
        return s.astype(float) * s
    else:
        return float(s) * s


def rankdata(a, method='average'):
    """
    Assign ranks to data, dealing with ties appropriately.

    Ranks begin at 1.  The `method` argument controls how ranks are assigned
    to equal values.  See [1]_ for further discussion of ranking methods.

    Parameters
    ----------
    a : array_like
        The array of values to be ranked.  The array is first flattened.
    method : str, optional
        The method used to assign ranks to tied elements.
        The options are 'average', 'min', 'max', 'dense' and 'ordinal'.

        'average':
            The average of the ranks that would have been assigned to
            all the tied values is assigned to each value.
        'min':
            The minimum of the ranks that would have been assigned to all
            the tied values is assigned to each value.  (This is also
            referred to as "competition" ranking.)
        'max':
            The maximum of the ranks that would have been assigned to all
            the tied values is assigned to each value.
        'dense':
            Like 'min', but the rank of the next highest element is assigned
            the rank immediately after those assigned to the tied elements.
        'ordinal':
            All values are given a distinct rank, corresponding to the order
            that the values occur in `a`.

        The default is 'average'.

    Returns
    -------
    ranks : ndarray
         An array of length equal to the size of `a`, containing rank
         scores.

    References
    ----------
    .. [1] "Ranking", https://en.wikipedia.org/wiki/Ranking

    Examples
    --------
    >>> from scipy.stats import rankdata
    >>> rankdata([0, 2, 3, 2])
    array([ 1. ,  2.5,  4. ,  2.5])
    >>> rankdata([0, 2, 3, 2], method='min')
    array([ 1,  2,  4,  2])
    >>> rankdata([0, 2, 3, 2], method='max')
    array([ 1,  3,  4,  3])
    >>> rankdata([0, 2, 3, 2], method='dense')
    array([ 1,  2,  3,  2])
    >>> rankdata([0, 2, 3, 2], method='ordinal')
    array([ 1,  2,  4,  3])
    """
    if method not in ('average', 'min', 'max', 'dense', 'ordinal'):
        raise ValueError('unknown method "{0}"'.format(method))

    arr = np.ravel(np.asarray(a))
    algo = 'mergesort' if method == 'ordinal' else 'quicksort'
    sorter = np.argsort(arr, kind=algo)

    inv = np.empty(sorter.size, dtype=np.intp)
    inv[sorter] = np.arange(sorter.size, dtype=np.intp)

    if method == 'ordinal':
        return inv + 1

    arr = arr[sorter]
    obs = np.r_[True, arr[1:] != arr[:-1]]
    dense = obs.cumsum()[inv]

    if method == 'dense':
        return dense

    # cumulative counts of each unique value
    count = np.r_[np.nonzero(obs)[0], len(obs)]

    if method == 'max':
        return count[dense]

    if method == 'min':
        return count[dense - 1] + 1

    # average method
    return .5 * (count[dense] + count[dense - 1] + 1)
