import math
from scipy import stats
from scipy._lib._array_api import (xp_capabilities, array_namespace, xp_promote,
                                   xp_result_type)
from scipy.stats._stats_py import _SimpleNormal, SignificanceResult, _get_pvalue
from scipy.stats._axis_nan_policy import _axis_nan_policy_factory


__all__ = ['chatterjeexi', 'spearmanrho']


def _xi_statistic(x, y, y_continuous, xp):
    # Compute xi correlation statistic

    # `axis=-1` is guaranteed by _axis_nan_policy decorator
    n = x.shape[-1]

    # "Rearrange the data as (X(1), Y(1)), . . . ,(X(n), Y(n))
    # such that X(1) ≤ ··· ≤ X(n)"
    j = xp.argsort(x, axis=-1)
    j, y = xp.broadcast_arrays(j, y)
    y = xp.take_along_axis(y, j, axis=-1)

    # "Let ri be the rank of Y(i), that is, the number of j such that Y(j) ≤ Y(i)"
    r = stats.rankdata(y, method='max', axis=-1)
    # " additionally define li to be the number of j such that Y(j) ≥ Y(i)"
    # Could probably compute this from r, but that can be an enhancement
    l = stats.rankdata(-y, method='max', axis=-1)
    r, l = xp.astype(r, x.dtype), xp.astype(l, x.dtype)

    num = xp.sum(xp.abs(xp.diff(r, axis=-1)), axis=-1)
    if y_continuous:  # [1] Eq. 1.1
        statistic = 1 - 3 * num / (n ** 2 - 1)
    else:  # [1] Eq. 1.2
        den = 2 * xp.sum((n - l) * l, axis=-1)
        statistic = 1 - n * num / den

    return statistic, r, l


def _xi_std(r, l, y_continuous, xp):
    # Compute asymptotic standard deviation of xi under null hypothesis of independence

    # `axis=-1` is guaranteed by _axis_nan_policy decorator
    n = r.shape[-1]

    # "Suppose that X and Y are independent and Y is continuous. Then
    # √n·ξn(X, Y) → N(0, 2/5) in distribution as n → ∞"
    if y_continuous:  # [1] Theorem 2.1
        return xp.asarray(math.sqrt(2 / 5) / math.sqrt(n), dtype=r.dtype)

    # "Suppose that X and Y are independent. Then √n·ξn(X, Y)
    # converges to N(0, τ²) in distribution as n → ∞
    # [1] Eq. 2.2 and surrounding math
    i = xp.arange(1, n + 1, dtype=r.dtype)
    u = xp.sort(r, axis=-1)
    v = xp.cumulative_sum(u, axis=-1)
    an = 1 / n**4 * xp.sum((2*n - 2*i + 1) * u**2, axis=-1)
    bn = 1 / n**5 * xp.sum((v + (n - i)*u)**2, axis=-1)
    cn = 1 / n**3 * xp.sum((2*n - 2*i + 1) * u, axis=-1)
    dn = 1 / n**3 * xp.sum((l * (n - l)), axis=-1)
    tau2 = (an - 2*bn + cn**2) / dn**2

    return xp.sqrt(tau2) / math.sqrt(n)


def _chatterjeexi_iv(y_continuous, method):
    # Input validation for `chatterjeexi`
    # x, y, `axis` input validation taken care of by decorator

    if y_continuous not in {True, False}:
        raise ValueError('`y_continuous` must be boolean.')

    if not isinstance(method, stats.PermutationMethod):
        method = method.lower()
        message = "`method` must be 'asymptotic' or a `PermutationMethod` instance."
        if method != 'asymptotic':
            raise ValueError(message)

    return y_continuous, method


def _unpack(res, _):
    return res.statistic, res.pvalue


@xp_capabilities(skip_backends=[('dask.array', 'no take_along_axis'),
                                ('cupy', 'no rankdata (xp.repeats limitation)')])
@_axis_nan_policy_factory(SignificanceResult, paired=True, n_samples=2,
                          result_to_tuple=_unpack, n_outputs=2, too_small=1)
def chatterjeexi(x, y, *, axis=0, y_continuous=False, method='asymptotic'):
    r"""Compute the xi correlation and perform a test of independence

    The xi correlation coefficient is a measure of association between two
    variables; the value tends to be close to zero when the variables are
    independent and close to 1 when there is a strong association. Unlike
    other correlation coefficients, the xi correlation is effective even
    when the association is not monotonic.

    Parameters
    ----------
    x, y : array-like
        The samples: corresponding observations of the independent and
        dependent variable. The (N-d) arrays must be broadcastable.
    axis : int, default: 0
        Axis along which to perform the test.
    method : 'asymptotic' or `PermutationMethod` instance, optional
        Selects the method used to calculate the *p*-value.
        Default is 'asymptotic'. The following options are available.

        * ``'asymptotic'``: compares the standardized test statistic
          against the normal distribution.
        * `PermutationMethod` instance. In this case, the p-value
          is computed using `permutation_test` with the provided
          configuration options and other appropriate settings.

    y_continuous : bool, default: False
        Whether `y` is assumed to be drawn from a continuous distribution.
        If `y` is drawn from a continuous distribution, results are valid
        whether this is assumed or not, but enabling this assumption will
        result in faster computation and typically produce similar results.

    Returns
    -------
    res : SignificanceResult
        An object containing attributes:

        statistic : float
            The xi correlation statistic.
        pvalue : float
            The associated *p*-value: the probability of a statistic at least as
            high as the observed value under the null hypothesis of independence.

    See Also
    --------
    scipy.stats.pearsonr, scipy.stats.spearmanr, scipy.stats.kendalltau

    Notes
    -----
    There is currently no special handling of ties in `x`; they are broken arbitrarily
    by the implementation. [1]_ recommends: "if there are ties among the Xi's, then
    choose an increasing rearrangement as above by breaking ties uniformly at random."
    This is easily accomplished by adding a small amount of random noise to `x`; see
    examples.

    [1]_ notes that the statistic is not symmetric in `x` and `y` *by design*:
    "...we may want to understand if :math:`Y` is a function :math:`X`, and not just
    if one of the variables is a function of the other." See [1]_ Remark 1.

    References
    ----------
    .. [1] Chatterjee, Sourav. "A new coefficient of correlation." Journal of
           the American Statistical Association 116.536 (2021): 2009-2022.
           :doi:`10.1080/01621459.2020.1758115`.

    Examples
    --------
    Generate perfectly correlated data, and observe that the xi correlation is
    nearly 1.0.

    >>> import numpy as np
    >>> from scipy import stats
    >>> rng = np.random.default_rng(348932549825235)
    >>> x = rng.uniform(0, 10, size=100)
    >>> y = np.sin(x)
    >>> res = stats.chatterjeexi(x, y)
    >>> res.statistic
    np.float64(0.9012901290129013)

    The probability of observing such a high value of the statistic under the
    null hypothesis of independence is very low.

    >>> res.pvalue
    np.float64(2.2206974648177804e-46)

    As noise is introduced, the correlation coefficient decreases.

    >>> noise = rng.normal(scale=[[0.1], [0.5], [1]], size=(3, 100))
    >>> res = stats.chatterjeexi(x, y + noise, axis=-1)
    >>> res.statistic
    array([0.79507951, 0.41824182, 0.16651665])

    Because the distribution of `y` is continuous, it is valid to pass
    ``y_continuous=True``. The statistic is identical, and the p-value
    (not shown) is only slightly different.

    >>> stats.chatterjeexi(x, y + noise, y_continuous=True, axis=-1).statistic
    array([0.79507951, 0.41824182, 0.16651665])

    Consider a case in which there are ties in `x`.

    >>> x = rng.integers(10, size=1000)
    >>> y = rng.integers(10, size=1000)

    [1]_ recommends breaking the ties uniformly at random.

    >>> d = rng.uniform(1e-5, size=x.size)
    >>> res = stats.chatterjeexi(x + d, y)
    >>> res.statistic
    -0.029919991638798438

    Since this gives a randomized estimate of the statistic, [1]_ also suggests
    considering the average over all possibilities of breaking ties. This is
    computationally infeasible when there are many ties, but a randomized estimate of
    *this* quantity can be obtained by considering many random possibilities of breaking
    ties.

    >>> d = rng.uniform(1e-5, size=(9999, x.size))
    >>> res = stats.chatterjeexi(x + d, y, axis=1)
    >>> np.mean(res.statistic)
    0.001186895213756626

    """
    xp = array_namespace(x, y)

    # x, y, `axis` input validation taken care of by decorator
    # In fact, `axis` is guaranteed to be -1
    y_continuous, method = _chatterjeexi_iv(y_continuous, method)
    x, y = xp_promote(x, y, force_floating=True, xp=xp)

    # A highly negative statistic is possible, e.g.
    # x = np.arange(100.), y = (x % 2 == 0)
    # Unclear whether we should expose `alternative`, though.
    alternative = 'greater'

    if method == 'asymptotic':
        xi, r, l = _xi_statistic(x, y, y_continuous, xp=xp)
        std = _xi_std(r, l, y_continuous, xp=xp)
        norm = _SimpleNormal()
        pvalue = _get_pvalue(xi / std, norm, alternative=alternative, xp=xp)
    elif isinstance(method, stats.PermutationMethod):
        res = stats.permutation_test(
            # Could be faster if we just permuted the ranks; for now, keep it simple.
            data=(y,),
            statistic=lambda y, axis: _xi_statistic(x, y, y_continuous, xp=xp)[0],
            alternative=alternative, permutation_type='pairings', **method._asdict(),
            axis=-1)  # `axis=-1` is guaranteed by _axis_nan_policy decorator

        xi, pvalue = res.statistic, res.pvalue

    xi = xi[()] if xi.ndim == 0 else xi
    pvalue = pvalue[()] if pvalue.ndim == 0 else pvalue
    return SignificanceResult(xi, pvalue)


@xp_capabilities(cpu_only=True, exceptions=['jax.numpy'],
    skip_backends=[('dask.array', 'not supported by rankdata (take_along_axis)')]
)
@_axis_nan_policy_factory(SignificanceResult, paired=True, n_samples=2,
                          result_to_tuple=_unpack, n_outputs=2, too_small=1)
def spearmanrho(x, y, /, *, alternative='two-sided', method=None, axis=0):
    r"""Calculate a Spearman rho correlation coefficient with associated p-value.

    The Spearman rank-order correlation coefficient is a nonparametric measure
    of the monotonicity of the relationship between two datasets.
    Like other correlation coefficients, it varies between -1 and +1 with 0
    implying no correlation. Coefficients of -1 or +1 are associated with an exact
    monotonic relationship.  Positive correlations indicate that as `x` increases,
    so does `y`; negative correlations indicate that as `x` increases, `y` decreases.
    The p-value is the probability of an uncorrelated system producing datasets
    with a Spearman correlation at least as extreme as the one computed from the
    observed dataset.

    Parameters
    ----------
    x, y : array-like
        The samples: corresponding observations of the independent and
        dependent variable. The (N-d) arrays must be broadcastable.
    alternative : {'two-sided', 'less', 'greater'}, optional
        Defines the alternative hypothesis. Default is 'two-sided'.
        The following options are available:

        * 'two-sided': the correlation is nonzero
        * 'less': the correlation is negative (less than zero)
        * 'greater':  the correlation is positive (greater than zero)

    method : ResamplingMethod, optional
        Defines the method used to compute the p-value. If `method` is an
        instance of `PermutationMethod`/`MonteCarloMethod`, the p-value is
        computed using
        `scipy.stats.permutation_test`/`scipy.stats.monte_carlo_test` with the
        provided configuration options and other appropriate settings.
        Otherwise, the p-value is computed using an asymptotic approximation of
        the null distribution.
    axis : int or None, optional
        If axis=0 (default), then each column represents a variable, with
        observations in the rows. If axis=1, the relationship is transposed:
        each row represents a variable, while the columns contain observations.
        If axis=None, then both arrays will be raveled.
        Like other `scipy.stats` functions, `axis` is interpreted after the
        arrays are broadcasted.

    Returns
    -------
    res : SignificanceResult
        An object containing attributes:

        statistic : floating point array or NumPy scalar
            Spearman correlation coefficient
        pvalue : floating point array NumPy scalar
            The p-value - the probabilitiy of realizing such an extreme statistic
            value under the null hypothesis that two samples have no ordinal
            correlation. See `alternative` above for alternative hypotheses.

    Warns
    -----
    `~scipy.stats.ConstantInputWarning`
        Raised if an input is a constant array.  The correlation coefficient
        is not defined in this case, so ``np.nan`` is returned.

    Notes
    -----
    `spearmanrho` was created to make improvements to SciPy's implementation of
    the Spearman correlation test without making backward-incompatible changes
    to `spearmanr`. Advantages of `spearmanrho` over `spearmanr` include:

    - `spearmanrho` follows standard array broadcasting rules.
    - `spearmanrho` is compatible with some non-NumPy arrays.
    - `spearmanrho` can compute exact p-values, even in the presence of ties,
      when an appropriate instance of `PermutationMethod` is provided via the
      `method` argument.

    References
    ----------
    .. [1] Zwillinger, D. and Kokoska, S. (2000). CRC Standard
       Probability and Statistics Tables and Formulae. Chapman & Hall: New
       York. 2000.
       Section  14.7
    .. [2] Kendall, M. G. and Stuart, A. (1973).
       The Advanced Theory of Statistics, Volume 2: Inference and Relationship.
       Griffin. 1973.
       Section 31.18

    Examples
    --------
    Univariate samples, approximate p-value.

    >>> import numpy as np
    >>> from scipy import stats
    >>> x = [1, 2, 3, 4, 5]
    >>> y = [5, 6, 7, 8, 7]
    >>> res = stats.spearmanrho(x, y)
    >>> res.statistic
    np.float64(0.8207826816681233)
    >>> res.pvalue
    np.float64(0.08858700531354405)

    Univariate samples, exact p-value.

    >>> res = stats.spearmanrho(x, y, method=stats.PermutationMethod())
    >>> res.statistic
    np.float64(0.8207826816681233)
    >>> res.pvalue
    np.float64(0.13333333333333333)

    Batch of univariate samples, one vectorized call.

    >>> rng = np.random.default_rng(98145152315484)
    >>> x2 = rng.standard_normal((2, 100))
    >>> y2 = rng.standard_normal((2, 100))
    >>> res = stats.spearmanrho(x2, y2, axis=-1)
    >>> res.statistic
    array([ 0.16585659, -0.12151215])
    >>> res.pvalue
    array([0.0991155 , 0.22846869])

    Bivariate samples using standard broadcasting rules.

    >>> res = stats.spearmanrho(x2[np.newaxis, :], x2[:, np.newaxis], axis=-1)
    >>> res.statistic
    array([[ 1.        , -0.14670267],
           [-0.14670267,  1.        ]])
    >>> res.pvalue
    array([[0.        , 0.14526128],
           [0.14526128, 0.        ]])

    """
    xp = array_namespace(x, y)
    dtype = xp_result_type(x, y, force_floating=True, xp=xp)
    rx = stats.rankdata(x, axis=axis)
    ry = stats.rankdata(y, axis=axis)
    rx = xp.astype(rx, dtype, copy=False)
    ry = xp.astype(ry, dtype, copy=False)
    res = stats.pearsonr(rx, ry, method=method, alternative=alternative, axis=axis)
    return SignificanceResult(res.statistic, res.pvalue)
