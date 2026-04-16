import numpy as np

from scipy import stats
from ._stats_py import _get_pvalue, _rankdata, _SimpleNormal
from . import _morestats
from ._axis_nan_policy import _broadcast_arrays
from ._hypotests import _get_wilcoxon_distr
from scipy._lib._util import _get_nan
from scipy._lib._array_api import array_namespace, xp_promote, xp_size
import scipy._lib.array_api_extra as xpx


class WilcoxonDistribution:

    def __init__(self, n):
        n = np.asarray(n).astype(int, copy=False)
        self.n = n
        self._dists = {ni: _get_wilcoxon_distr(ni) for ni in np.unique(n)}

    def _cdf1(self, k, n):
        pmfs = self._dists[n]
        return pmfs[:k + 1].sum()

    def _cdf(self, k, n):
        return np.vectorize(self._cdf1, otypes=[float])(k, n)

    def _sf1(self, k, n):
        pmfs = self._dists[n]
        return pmfs[k:].sum()

    def _sf(self, k, n):
        return np.vectorize(self._sf1, otypes=[float])(k, n)

    def mean(self):
        return self.n * (self.n + 1) / 4

    def _prep(self, k):
        k = np.asarray(k).astype(int, copy=False)
        mn = self.mean()
        out = np.empty(k.shape, dtype=np.float64)
        return k, mn, out

    def cdf(self, k):
        k, mn, out = self._prep(k)
        return xpx.apply_where(
            k <= mn, (k, self.n),
            self._cdf,
            lambda k, n: 1 - self._sf(k+1, n))[()]

    def sf(self, k):
        k, mn, out = self._prep(k)
        return xpx.apply_where(
            k <= mn, (k, self.n),
            self._sf,
            lambda k, n: 1 - self._cdf(k-1, n))[()]


def _wilcoxon_iv(x, y, zero_method, correction, alternative, method, axis):
    xp = array_namespace(x, y)
    x, y = xp_promote(x, y, force_floating=True, xp=xp)

    axis = np.asarray(axis)[()]  # OK to use NumPy for input validation
    message = "`axis` must be an integer."
    if not np.issubdtype(axis.dtype, np.integer) or axis.ndim != 0:
        raise ValueError(message)
    axis = int(axis)

    message = '`axis` must be compatible with the shape(s) of `x` (and `y`)'
    AxisError = getattr(np, 'AxisError', None) or np.exceptions.AxisError
    try:
        if y is None:
            d = x
        else:
            x, y = _broadcast_arrays((x, y), axis=axis, xp=xp)
            d = x - y
        d = xp.moveaxis(d, axis, -1)
    except AxisError as e:
        raise AxisError(message) from e

    message = "`x` and `y` must have the same length along `axis`."
    if y is not None and x.shape[axis] != y.shape[axis]:
        raise ValueError(message)

    message = "`x` (and `y`, if provided) must be an array of real numbers."
    if not xp.isdtype(d.dtype, "real floating"):
        raise ValueError(message)

    zero_method = str(zero_method).lower()
    zero_methods = {"wilcox", "pratt", "zsplit"}
    message = f"`zero_method` must be one of {zero_methods}."
    if zero_method not in zero_methods:
        raise ValueError(message)

    corrections = {True, False}
    message = f"`correction` must be one of {corrections}."
    if correction not in corrections:
        raise ValueError(message)

    alternative = str(alternative).lower()
    alternatives = {"two-sided", "less", "greater"}
    message = f"`alternative` must be one of {alternatives}."
    if alternative not in alternatives:
        raise ValueError(message)

    if not isinstance(method, stats.PermutationMethod):
        methods = {"auto", "asymptotic", "exact"}
        message = (f"`method` must be one of {methods} or "
                   "an instance of `stats.PermutationMethod`.")
        if method not in methods:
            raise ValueError(message)
    output_z = True if method == 'asymptotic' else False

    # For small samples, we decide later whether to perform an exact test or a
    # permutation test. The reason is that the presence of ties is not
    # known at the input validation stage.
    n_zero = xp.count_nonzero(d == 0, axis=None)
    if method == "auto" and d.shape[-1] > 50:
        method = "asymptotic"

    return d, zero_method, correction, alternative, method, axis, output_z, n_zero, xp


def _wilcoxon_statistic(d, method, zero_method='wilcox', *, xp):
    dtype = d.dtype
    i_zeros = (d == 0)

    if zero_method == 'wilcox':
        # Wilcoxon's method for treating zeros was to remove them from
        # the calculation. We do this by replacing 0s with NaNs, which
        # are ignored anyway.
        # Copy required for array-api-strict. See data-apis/array-api-extra#506.
        d = xpx.at(d)[i_zeros].set(xp.nan, copy=True)

    i_nan = xp.isnan(d)
    n_nan = xp.count_nonzero(i_nan, axis=-1)
    count = xp.astype(d.shape[-1] - n_nan, dtype)

    r, t = _rankdata(xp.abs(d), 'average', return_ties=True, xp=xp)
    r, t = xp.astype(r, dtype, copy=False), xp.astype(t, dtype, copy=False)

    r_plus = xp.sum(xp.astype(d > 0, dtype) * r, axis=-1)
    r_minus = xp.sum(xp.astype(d < 0, dtype) * r, axis=-1)

    has_ties = xp.any(t == 0)

    if zero_method == "zsplit":
        # The "zero-split" method for treating zeros is to add half their contribution
        # to r_plus and half to r_minus.
        # See gh-2263 for the origin of this method.
        r_zero_2 = xp.sum(xp.astype(i_zeros, dtype) * r, axis=-1) / 2
        r_plus = xpx.at(r_plus)[...].add(r_zero_2)
        r_minus = xpx.at(r_minus)[...].add(r_zero_2)

    mn = count * (count + 1.) * 0.25
    se = count * (count + 1.) * (2. * count + 1.)

    if zero_method == "pratt":
        # Pratt's method for treating zeros was just to modify the z-statistic.

        # normal approximation needs to be adjusted, see Cureton (1967)
        n_zero = xp.astype(xp.count_nonzero(i_zeros, axis=-1), dtype)
        mn = xpx.at(mn)[...].subtract(n_zero * (n_zero + 1.) * 0.25)
        se = xpx.at(se)[...].subtract(n_zero * (n_zero + 1.) * (2. * n_zero + 1.))

        # zeros are not to be included in tie-correction.
        # any tie counts corresponding with zeros are in the 0th column
        # t[xp.any(i_zeros, axis=-1), 0] = 0
        t_i_zeros = xp.zeros_like(i_zeros)
        t_i_zeros = xpx.at(t_i_zeros)[..., 0].set(xp.any(i_zeros, axis=-1))
        t = xpx.at(t)[t_i_zeros].set(0.)

    tie_correct = xp.sum(t**3 - t, axis=-1)
    se = xp.sqrt((se - tie_correct/2) / 24)

    # se = 0 means that no non-zero values are left in d. we only need z
    # if method is asymptotic. however, if method="auto", the switch to
    # asymptotic might only happen after the statistic is calculated, so z
    # needs to be computed. in all other cases, avoid division by zero warning
    # (z is not needed anyways)
    if method in ["asymptotic", "auto"]:
        z = (r_plus - mn) / se
    else:
        z = xp.nan

    return r_plus, r_minus, se, z, count, has_ties


def _correction_sign(z, alternative, xp):
    if alternative == 'greater':
        return 1
    elif alternative == 'less':
        return -1
    else:
        return xp.sign(z)


def _wilcoxon_nd(x, y=None, zero_method='wilcox', correction=True,
                 alternative='two-sided', method='auto', axis=0):

    temp = _wilcoxon_iv(x, y, zero_method, correction, alternative, method, axis)
    d, zero_method, correction, alternative, method, axis, output_z, n_zero, xp = temp

    if xp_size(d) == 0:
        NaN = _get_nan(d, xp=xp)
        res = _morestats.WilcoxonResult(statistic=NaN, pvalue=NaN)
        if method == 'asymptotic':
            res.zstatistic = NaN
        return res

    r_plus, r_minus, se, z, count, has_ties = _wilcoxon_statistic(
        d, method, zero_method, xp=xp
    )

    # we only know if there are ties after computing the statistic and not
    # at the input validation stage. if the original method was auto and
    # the decision was to use an exact test, we override this to
    # a permutation test now (since method='exact' is not exact in the
    # presence of ties)
    if method == "auto":
        if not (has_ties or n_zero > 0):
            method = "exact"
        elif d.shape[-1] <= 13:
            # the possible outcomes to be simulated by the permutation test
            # are 2**n, where n is the sample size.
            # if n <= 13, the p-value is deterministic since 2**13 is less
            # than 9999, the default number of n_resamples
            method = stats.PermutationMethod()
        else:
            # if there are ties and the sample size is too large to
            # run a deterministic permutation test, fall back to asymptotic
            method = "asymptotic"

    if method == 'asymptotic':
        if correction:
            sign = _correction_sign(z, alternative, xp=xp)
            z = xpx.at(z)[...].subtract(sign * 0.5 / se)
        p = _get_pvalue(z, _SimpleNormal(), alternative, xp=xp)
    elif method == 'exact':
        dist = WilcoxonDistribution(count)
        # The null distribution in `dist` is exact only if there are no ties
        # or zeros. If there are ties or zeros, the statistic can be non-
        # integral, but the null distribution is only defined for integral
        # values of the statistic. Therefore, we're conservative: round
        # non-integral statistic up before computing CDF and down before
        # computing SF. This preserves symmetry w.r.t. alternatives and
        # order of the input arguments. See gh-19872.
        r_plus_np = np.asarray(r_plus)
        if alternative == 'less':
            p = dist.cdf(np.ceil(r_plus_np))
        elif alternative == 'greater':
            p = dist.sf(np.floor(r_plus_np))
        else:
            p = 2 * np.minimum(dist.sf(np.floor(r_plus_np)),
                               dist.cdf(np.ceil(r_plus_np)))
            p = np.clip(p, 0, 1)
        p = xp.asarray(p, dtype=d.dtype)
    else:  # `PermutationMethod` instance (already validated)
        p = stats.permutation_test(
            (d,), lambda d: _wilcoxon_statistic(d, method, zero_method, xp=xp)[0],
            permutation_type='samples', **method._asdict(),
            alternative=alternative, axis=-1).pvalue

    # for backward compatibility...
    statistic = xp.minimum(r_plus, r_minus) if alternative=='two-sided' else r_plus
    z = -xp.abs(z) if (alternative == 'two-sided' and method == 'asymptotic') else z

    statistic = statistic[()] if statistic.ndim == 0 else statistic
    p = p[()] if p.ndim == 0 else p
    res = _morestats.WilcoxonResult(statistic=statistic, pvalue=p)
    if output_z:
        res.zstatistic = z[()] if z.ndim == 0 else z
    return res
