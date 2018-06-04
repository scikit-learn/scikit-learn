from collections import namedtuple

import numpy as np

from . import distributions


__all__ = ['_find_repeats', 'linregress', 'theilslopes']

LinregressResult = namedtuple('LinregressResult', ('slope', 'intercept',
                                                   'rvalue', 'pvalue',
                                                   'stderr'))

def linregress(x, y=None):
    """
    Calculate a linear least-squares regression for two sets of measurements.

    Parameters
    ----------
    x, y : array_like
        Two sets of measurements.  Both arrays should have the same length.
        If only x is given (and y=None), then it must be a two-dimensional
        array where one dimension has length 2.  The two sets of measurements
        are then found by splitting the array along the length-2 dimension.

    Returns
    -------
    slope : float
        slope of the regression line
    intercept : float
        intercept of the regression line
    rvalue : float
        correlation coefficient
    pvalue : float
        two-sided p-value for a hypothesis test whose null hypothesis is
        that the slope is zero, using Wald Test with t-distribution of
        the test statistic.
    stderr : float
        Standard error of the estimated gradient.

    See also
    --------
    :func:`scipy.optimize.curve_fit` : Use non-linear
     least squares to fit a function to data.
    :func:`scipy.optimize.leastsq` : Minimize the sum of
     squares of a set of equations.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from scipy import stats
    >>> np.random.seed(12345678)
    >>> x = np.random.random(10)
    >>> y = np.random.random(10)
    >>> slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

    To get coefficient of determination (r_squared)

    >>> print("r-squared:", r_value**2)
    r-squared: 0.08040226853902833

    Plot the data along with the fitted line

    >>> plt.plot(x, y, 'o', label='original data')
    >>> plt.plot(x, intercept + slope*x, 'r', label='fitted line')
    >>> plt.legend()
    >>> plt.show()

    """
    TINY = 1.0e-20
    if y is None:  # x is a (2, N) or (N, 2) shaped array_like
        x = np.asarray(x)
        if x.shape[0] == 2:
            x, y = x
        elif x.shape[1] == 2:
            x, y = x.T
        else:
            msg = ("If only `x` is given as input, it has to be of shape "
                   "(2, N) or (N, 2), provided shape was %s" % str(x.shape))
            raise ValueError(msg)
    else:
        x = np.asarray(x)
        y = np.asarray(y)

    if x.size == 0 or y.size == 0:
        raise ValueError("Inputs must not be empty.")

    n = len(x)
    xmean = np.mean(x, None)
    ymean = np.mean(y, None)

    # average sum of squares:
    ssxm, ssxym, ssyxm, ssym = np.cov(x, y, bias=1).flat
    r_num = ssxym
    r_den = np.sqrt(ssxm * ssym)
    if r_den == 0.0:
        r = 0.0
    else:
        r = r_num / r_den
        # test for numerical error propagation
        if r > 1.0:
            r = 1.0
        elif r < -1.0:
            r = -1.0

    df = n - 2
    slope = r_num / ssxm
    intercept = ymean - slope*xmean
    if n == 2:
        # handle case when only two points are passed in
        if y[0] == y[1]:
            prob = 1.0
        else:
            prob = 0.0
        sterrest = 0.0
    else:
        t = r * np.sqrt(df / ((1.0 - r + TINY)*(1.0 + r + TINY)))
        prob = 2 * distributions.t.sf(np.abs(t), df)
        sterrest = np.sqrt((1 - r**2) * ssym / ssxm / df)

    return LinregressResult(slope, intercept, r, prob, sterrest)


def theilslopes(y, x=None, alpha=0.95):
    r"""
    Computes the Theil-Sen estimator for a set of points (x, y).

    `theilslopes` implements a method for robust linear regression.  It
    computes the slope as the median of all slopes between paired values.

    Parameters
    ----------
    y : array_like
        Dependent variable.
    x : array_like or None, optional
        Independent variable. If None, use ``arange(len(y))`` instead.
    alpha : float, optional
        Confidence degree between 0 and 1. Default is 95% confidence.
        Note that `alpha` is symmetric around 0.5, i.e. both 0.1 and 0.9 are
        interpreted as "find the 90% confidence interval".

    Returns
    -------
    medslope : float
        Theil slope.
    medintercept : float
        Intercept of the Theil line, as ``median(y) - medslope*median(x)``.
    lo_slope : float
        Lower bound of the confidence interval on `medslope`.
    up_slope : float
        Upper bound of the confidence interval on `medslope`.

    Notes
    -----
    The implementation of `theilslopes` follows [1]_. The intercept is
    not defined in [1]_, and here it is defined as ``median(y) -
    medslope*median(x)``, which is given in [3]_. Other definitions of
    the intercept exist in the literature. A confidence interval for
    the intercept is not given as this question is not addressed in
    [1]_.

    References
    ----------
    .. [1] P.K. Sen, "Estimates of the regression coefficient based on Kendall's tau",
           J. Am. Stat. Assoc., Vol. 63, pp. 1379-1389, 1968.
    .. [2] H. Theil, "A rank-invariant method of linear and polynomial
           regression analysis I, II and III",  Nederl. Akad. Wetensch., Proc.
           53:, pp. 386-392, pp. 521-525, pp. 1397-1412, 1950.
    .. [3] W.L. Conover, "Practical nonparametric statistics", 2nd ed.,
           John Wiley and Sons, New York, pp. 493.

    Examples
    --------
    >>> from scipy import stats
    >>> import matplotlib.pyplot as plt

    >>> x = np.linspace(-5, 5, num=150)
    >>> y = x + np.random.normal(size=x.size)
    >>> y[11:15] += 10  # add outliers
    >>> y[-5:] -= 7

    Compute the slope, intercept and 90% confidence interval.  For comparison,
    also compute the least-squares fit with `linregress`:

    >>> res = stats.theilslopes(y, x, 0.90)
    >>> lsq_res = stats.linregress(x, y)

    Plot the results. The Theil-Sen regression line is shown in red, with the
    dashed red lines illustrating the confidence interval of the slope (note
    that the dashed red lines are not the confidence interval of the regression
    as the confidence interval of the intercept is not included). The green
    line shows the least-squares fit for comparison.

    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111)
    >>> ax.plot(x, y, 'b.')
    >>> ax.plot(x, res[1] + res[0] * x, 'r-')
    >>> ax.plot(x, res[1] + res[2] * x, 'r--')
    >>> ax.plot(x, res[1] + res[3] * x, 'r--')
    >>> ax.plot(x, lsq_res[1] + lsq_res[0] * x, 'g-')
    >>> plt.show()

    """
    # We copy both x and y so we can use _find_repeats.
    y = np.array(y).flatten()
    if x is None:
        x = np.arange(len(y), dtype=float)
    else:
        x = np.array(x, dtype=float).flatten()
        if len(x) != len(y):
            raise ValueError("Incompatible lengths ! (%s<>%s)" % (len(y), len(x)))

    # Compute sorted slopes only when deltax > 0
    deltax = x[:, np.newaxis] - x
    deltay = y[:, np.newaxis] - y
    slopes = deltay[deltax > 0] / deltax[deltax > 0]
    slopes.sort()
    medslope = np.median(slopes)
    medinter = np.median(y) - medslope * np.median(x)
    # Now compute confidence intervals
    if alpha > 0.5:
        alpha = 1. - alpha

    z = distributions.norm.ppf(alpha / 2.)
    # This implements (2.6) from Sen (1968)
    _, nxreps = _find_repeats(x)
    _, nyreps = _find_repeats(y)
    nt = len(slopes)       # N in Sen (1968)
    ny = len(y)            # n in Sen (1968)
    # Equation 2.6 in Sen (1968):
    sigsq = 1/18. * (ny * (ny-1) * (2*ny+5) -
                     sum(k * (k-1) * (2*k + 5) for k in nxreps) -
                     sum(k * (k-1) * (2*k + 5) for k in nyreps))
    # Find the confidence interval indices in `slopes`
    sigma = np.sqrt(sigsq)
    Ru = min(int(np.round((nt - z*sigma)/2.)), len(slopes)-1)
    Rl = max(int(np.round((nt + z*sigma)/2.)) - 1, 0)
    delta = slopes[[Rl, Ru]]
    return medslope, medinter, delta[0], delta[1]


def _find_repeats(arr):
    # This function assumes it may clobber its input.
    if len(arr) == 0:
        return np.array(0, np.float64), np.array(0, np.intp)

    # XXX This cast was previously needed for the Fortran implementation,
    # should we ditch it?
    arr = np.asarray(arr, np.float64).ravel()
    arr.sort()

    # Taken from NumPy 1.9's np.unique.
    change = np.concatenate(([True], arr[1:] != arr[:-1]))
    unique = arr[change]
    change_idx = np.concatenate(np.nonzero(change) + ([arr.size],))
    freq = np.diff(change_idx)
    atleast2 = freq > 1
    return unique[atleast2], freq[atleast2]
