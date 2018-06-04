"""
Statistical functions and tests, following scipy.stats.

Some differences

- We don't handle missing values at all

"""
# This is lightly adapted from scipy.stats 0.19
# http://github.com/scipy/scipy/blob/v0.19.0/scipy/stats/stats.py
# The original copyright notice follows:

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
import math

import numpy as np
import dask.array as da
from dask.array.random import doc_wraps
from dask.array.ufunc import wrap_elemwise
from dask import delayed
try:
    import scipy.stats
except ImportError:
    raise ImportError("`dask.array.stats` requires `scipy` to be installed.")
from scipy.stats import distributions
from scipy import special
from scipy.stats.stats import (Ttest_indResult, Ttest_1sampResult,
                               Ttest_relResult,
                               Power_divergenceResult,
                               NormaltestResult,
                               SkewtestResult,
                               KurtosistestResult,
                               F_onewayResult)


__all__ = ['ttest_ind', 'ttest_1samp', 'ttest_rel', 'chisquare',
           'power_divergence', 'skew', 'skewtest', 'kurtosis',
           'kurtosistest', 'normaltest', 'f_oneway', 'moment']

# -----------------
# Statistical Tests
# -----------------


@doc_wraps(scipy.stats.ttest_ind)
def ttest_ind(a, b, axis=0, equal_var=True):
    v1 = da.var(a, axis, ddof=1)  # XXX: np -> da
    v2 = da.var(b, axis, ddof=1)  # XXX: np -> da
    n1 = a.shape[axis]
    n2 = b.shape[axis]

    if equal_var:
        df, denom = _equal_var_ttest_denom(v1, n1, v2, n2)
    else:
        df, denom = _unequal_var_ttest_denom(v1, n1, v2, n2)

    res = _ttest_ind_from_stats(da.mean(a, axis), da.mean(b, axis), denom, df)

    return delayed(Ttest_indResult, nout=2)(*res)


@doc_wraps(scipy.stats.ttest_1samp)
def ttest_1samp(a, popmean, axis=0, nan_policy='propagate'):
    if nan_policy != 'propagate':
        raise NotImplementedError("`nan_policy` other than 'propagate' "
                                  "have not been implemented.")
    n = a.shape[axis]
    df = n - 1

    d = da.mean(a, axis) - popmean
    v = da.var(a, axis, ddof=1)
    denom = da.sqrt(v / float(n))

    with np.errstate(divide='ignore', invalid='ignore'):
        t = da.divide(d, denom)
    t, prob = _ttest_finish(df, t)
    return delayed(Ttest_1sampResult, nout=2)(t, prob)


@doc_wraps(scipy.stats.ttest_rel)
def ttest_rel(a, b, axis=0, nan_policy='propagate'):
    if nan_policy != 'propagate':
        raise NotImplementedError("`nan_policy` other than 'propagate' "
                                  "have not been implemented.")

    n = a.shape[axis]
    df = float(n - 1)

    d = (a - b).astype(np.float64)
    v = da.var(d, axis, ddof=1)
    dm = da.mean(d, axis)
    denom = da.sqrt(v / float(n))

    with np.errstate(divide='ignore', invalid='ignore'):
        t = da.divide(dm, denom)
    t, prob = _ttest_finish(df, t)

    return delayed(Ttest_relResult, nout=2)(t, prob)


@doc_wraps(scipy.stats.chisquare)
def chisquare(f_obs, f_exp=None, ddof=0, axis=0):
    return power_divergence(f_obs, f_exp=f_exp, ddof=ddof, axis=axis,
                            lambda_="pearson")


@doc_wraps(scipy.stats.power_divergence)
def power_divergence(f_obs, f_exp=None, ddof=0, axis=0, lambda_=None):

    if isinstance(lambda_, str):
        # TODO: public api
        if lambda_ not in scipy.stats.stats._power_div_lambda_names:
            names = repr(list(scipy.stats.stats._power_div_lambda_names.keys()))[1:-1]
            raise ValueError("invalid string for lambda_: {0!r}.  Valid strings "
                             "are {1}".format(lambda_, names))
        lambda_ = scipy.stats.stats._power_div_lambda_names[lambda_]
    elif lambda_ is None:
        lambda_ = 1

    if f_exp is not None:
        # f_exp = np.atleast_1d(np.asanyarray(f_exp))
        pass
    else:
        f_exp = f_obs.mean(axis=axis, keepdims=True)

    # `terms` is the array of terms that are summed along `axis` to create
    # the test statistic.  We use some specialized code for a few special
    # cases of lambda_.
    if lambda_ == 1:
        # Pearson's chi-squared statistic
        terms = (f_obs - f_exp)**2 / f_exp
    elif lambda_ == 0:
        # Log-likelihood ratio (i.e. G-test)
        terms = 2.0 * _xlogy(f_obs, f_obs / f_exp)
    elif lambda_ == -1:
        # Modified log-likelihood ratio
        terms = 2.0 * _xlogy(f_exp, f_exp / f_obs)
    else:
        # General Cressie-Read power divergence.
        terms = f_obs * ((f_obs / f_exp)**lambda_ - 1)
        terms /= 0.5 * lambda_ * (lambda_ + 1)

    stat = terms.sum(axis=axis)

    num_obs = _count(terms, axis=axis)
    # ddof = asarray(ddof)
    p = delayed(distributions.chi2.sf)(stat, num_obs - 1 - ddof)

    return delayed(Power_divergenceResult, nout=2)(stat, p)


@doc_wraps(scipy.stats.skew)
def skew(a, axis=0, bias=True, nan_policy='propagate'):
    if nan_policy != 'propagate':
        raise NotImplementedError("`nan_policy` other than 'propagate' "
                                  "have not been implemented.")

    n = a.shape[axis]  # noqa; for bias
    m2 = moment(a, 2, axis)
    m3 = moment(a, 3, axis)
    zero = (m2 == 0)
    vals = da.where(~zero, m3 / m2**1.5, 0.)
    # vals = da.where(~zero, (m2, m3),
    #                 lambda m2, m3: m3 / m2**1.5,
    #                 0.)
    if not bias:
        # Need a version of np.place
        raise NotImplementedError("bias=False is not implemented.")

    if vals.ndim == 0:
        return vals
        # TODO: scalar
        # return vals.item()

    return vals


@doc_wraps(scipy.stats.skewtest)
def skewtest(a, axis=0, nan_policy='propagate'):
    if nan_policy != 'propagate':
        raise NotImplementedError("`nan_policy` other than 'propagate' "
                                  "have not been implemented.")

    b2 = skew(a, axis)
    n = float(a.shape[axis])
    if n < 8:
        raise ValueError(
            "skewtest is not valid with less than 8 samples; %i samples"
            " were given." % int(n))
    y = b2 * math.sqrt(((n + 1) * (n + 3)) / (6.0 * (n - 2)))
    beta2 = (3.0 * (n**2 + 27 * n - 70) * (n + 1) * (n + 3) /
             ((n - 2.0) * (n + 5) * (n + 7) * (n + 9)))
    W2 = -1 + math.sqrt(2 * (beta2 - 1))
    delta = 1 / math.sqrt(0.5 * math.log(W2))
    alpha = math.sqrt(2.0 / (W2 - 1))
    y = np.where(y == 0, 1, y)
    Z = delta * np.log(y / alpha + np.sqrt((y / alpha)**2 + 1))

    return delayed(SkewtestResult, nout=2)(Z, 2 * distributions.norm.sf(np.abs(Z)))


@doc_wraps(scipy.stats.kurtosis)
def kurtosis(a, axis=0, fisher=True, bias=True, nan_policy='propagate'):
    if nan_policy != 'propagate':
        raise NotImplementedError("`nan_policy` other than 'propagate' "
                                  "have not been implemented.")
    n = a.shape[axis]  # noqa; for bias
    m2 = moment(a, 2, axis)
    m4 = moment(a, 4, axis)
    zero = (m2 == 0)
    olderr = np.seterr(all='ignore')
    try:
        vals = da.where(zero, 0, m4 / m2**2.0)
    finally:
        np.seterr(**olderr)

    if not bias:
        # need a version of np.place
        raise NotImplementedError("bias=False is not implemented.")

    if vals.ndim == 0:
        return vals  # TODO: scalar
        # vals = vals.item()  # array scalar

    if fisher:
        return vals - 3
    else:
        return vals
        # TODO: scalar; vals = vals.item()  # array scalar


@doc_wraps(scipy.stats.kurtosistest)
def kurtosistest(a, axis=0, nan_policy='propagate'):
    if nan_policy != 'propagate':
        raise NotImplementedError("`nan_policy` other than 'propagate' "
                                  "have not been implemented.")

    n = float(a.shape[axis])
    b2 = kurtosis(a, axis, fisher=False)

    E = 3.0 * (n - 1) / (n + 1)
    varb2 = 24.0 * n * (n - 2) * (n - 3) / ((n + 1) * (n + 1.) * (n + 3) * (n + 5))  # [1]_ Eq. 1
    x = (b2 - E) / np.sqrt(varb2)  # [1]_ Eq. 4
    # [1]_ Eq. 2:
    sqrtbeta1 = 6.0 * (n * n - 5 * n + 2) / ((n + 7) * (n + 9)) * np.sqrt((6.0 * (n + 3) * (n + 5)) /
                                                                          (n * (n - 2) * (n - 3)))
    # [1]_ Eq. 3:
    A = 6.0 + 8.0 / sqrtbeta1 * (2.0 / sqrtbeta1 + np.sqrt(1 + 4.0 / (sqrtbeta1**2)))
    term1 = 1 - 2 / (9.0 * A)
    denom = 1 + x * np.sqrt(2 / (A - 4.0))
    denom = np.where(denom < 0, 99, denom)
    term2 = np.where(denom < 0, term1, np.power((1 - 2.0 / A) / denom, 1 / 3.0))
    Z = (term1 - term2) / np.sqrt(2 / (9.0 * A))  # [1]_ Eq. 5
    Z = np.where(denom == 99, 0, Z)
    if Z.ndim == 0:
        Z = Z[()]

    # zprob uses upper tail, so Z needs to be positive
    return delayed(KurtosistestResult, nout=2)(Z, 2 * distributions.norm.sf(np.abs(Z)))


@doc_wraps(scipy.stats.normaltest)
def normaltest(a, axis=0, nan_policy='propagate'):
    if nan_policy != 'propagate':
        raise NotImplementedError("`nan_policy` other than 'propagate' "
                                  "have not been implemented.")

    s, _ = skewtest(a, axis)
    k, _ = kurtosistest(a, axis)
    k2 = s * s + k * k
    return delayed(NormaltestResult, nout=2)(k2, delayed(distributions.chi2.sf)(k2, 2))


@doc_wraps(scipy.stats.f_oneway)
def f_oneway(*args):
    # args = [np.asarray(arg, dtype=float) for arg in args]
    # ANOVA on N groups, each in its own array
    num_groups = len(args)
    alldata = da.concatenate(args)
    bign = len(alldata)

    # Determine the mean of the data, and subtract that from all inputs to a
    # variance (via sum_of_sq / sq_of_sum) calculation.  Variance is invariance
    # to a shift in location, and centering all data around zero vastly
    # improves numerical stability.
    offset = alldata.mean()
    alldata -= offset

    sstot = _sum_of_squares(alldata) - (_square_of_sums(alldata) / float(bign))
    ssbn = 0
    for a in args:
        ssbn += _square_of_sums(a - offset) / float(len(a))

    # Naming: variables ending in bn/b are for "between treatments", wn/w are
    # for "within treatments"
    ssbn -= (_square_of_sums(alldata) / float(bign))
    sswn = sstot - ssbn
    dfbn = num_groups - 1
    dfwn = bign - num_groups
    msb = ssbn / float(dfbn)
    msw = sswn / float(dfwn)
    f = msb / msw

    prob = _fdtrc(dfbn, dfwn, f)   # equivalent to stats.f.sf

    return delayed(F_onewayResult, nout=2)(f, prob)


@doc_wraps(scipy.stats.moment)
def moment(a, moment=1, axis=0, nan_policy='propagate'):
    if nan_policy != 'propagate':
        raise NotImplementedError("`nan_policy` other than 'propagate' "
                                  "have not been implemented.")
    return da.moment(a, moment, axis=axis)

# -------
# Helpers
# -------
# Don't really want to do all of scipy.special (or do we?)


_xlogy = wrap_elemwise(special.xlogy)
_fdtrc = wrap_elemwise(special.fdtrc)


def _equal_var_ttest_denom(v1, n1, v2, n2):
    df = n1 + n2 - 2.0
    svar = ((n1 - 1) * v1 + (n2 - 1) * v2) / df
    denom = da.sqrt(svar * (1.0 / n1 + 1.0 / n2))  # XXX: np -> da
    return df, denom


def _unequal_var_ttest_denom(v1, n1, v2, n2):
    vn1 = v1 / n1
    vn2 = v2 / n2
    with np.errstate(divide='ignore', invalid='ignore'):
        df = (vn1 + vn2)**2 / (vn1**2 / (n1 - 1) + vn2**2 / (n2 - 1))

    # If df is undefined, variances are zero (assumes n1 > 0 & n2 > 0).
    # Hence it doesn't matter what df is as long as it's not NaN.
    df = da.where(da.isnan(df), 1, df)  # XXX: np -> da
    denom = da.sqrt(vn1 + vn2)
    return df, denom


def _ttest_ind_from_stats(mean1, mean2, denom, df):

    d = mean1 - mean2
    with np.errstate(divide='ignore', invalid='ignore'):
        t = da.divide(d, denom)
    t, prob = _ttest_finish(df, t)

    return (t, prob)


def _ttest_finish(df, t):
    """Common code between all 3 t-test functions."""
    # XXX: np.abs -> da.absolute
    # XXX: delayed(distributions.t.sf)
    prob = delayed(distributions.t.sf)(da.absolute(t),
                                       df) * 2  # use np.abs to get upper tail
    if t.ndim == 0:
        t = t[()]

    return t, prob


def _count(x, axis=None):
    if axis is None:
        return x.size
    else:
        return x.shape[axis]


def _sum_of_squares(a, axis=0):
    """
    Squares each element of the input array, and returns the sum(s) of that.
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
    return da.sum(a * a, axis)


def _square_of_sums(a, axis=0):
    """
    Sums elements of the input array, and returns the square(s) of that sum.
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
    s = da.sum(a, axis)
    return s * s
