#
# Author:  Travis Oliphant  2002-2011 with contributions from
#          SciPy Developers 2004-2011
#
from __future__ import division, print_function, absolute_import

import warnings

import numpy as np

from scipy.misc.doccer import (extend_notes_in_docstring,
                               replace_notes_in_docstring)
from scipy import optimize
from scipy import integrate
import scipy.special as sc
from scipy._lib._numpy_compat import broadcast_to

from . import _stats
from ._tukeylambda_stats import (tukeylambda_variance as _tlvar,
                                 tukeylambda_kurtosis as _tlkurt)
from ._distn_infrastructure import (get_distribution_names, _kurtosis,
                                    _lazyselect, _lazywhere, _ncx2_cdf,
                                    _ncx2_log_pdf, _ncx2_pdf,
                                    rv_continuous, _skew, valarray)
from ._constants import _XMIN, _EULER, _ZETA3, _XMAX, _LOGXMAX


# In numpy 1.12 and above, np.power refuses to raise integers to negative
# powers, and `np.float_power` is a new replacement.
try:
    float_power = np.float_power
except AttributeError:
    float_power = np.power


## Kolmogorov-Smirnov one-sided and two-sided test statistics
class ksone_gen(rv_continuous):
    """General Kolmogorov-Smirnov one-sided test.

    %(default)s

    """
    def _cdf(self, x, n):
        return 1.0 - sc.smirnov(n, x)

    def _ppf(self, q, n):
        return sc.smirnovi(n, 1.0 - q)


ksone = ksone_gen(a=0.0, name='ksone')


class kstwobign_gen(rv_continuous):
    """Kolmogorov-Smirnov two-sided test for large N.

    %(default)s

    """
    def _cdf(self, x):
        return 1.0 - sc.kolmogorov(x)

    def _sf(self, x):
        return sc.kolmogorov(x)

    def _ppf(self, q):
        return sc.kolmogi(1.0 - q)


kstwobign = kstwobign_gen(a=0.0, name='kstwobign')


## Normal distribution

# loc = mu, scale = std
# Keep these implementations out of the class definition so they can be reused
# by other distributions.
_norm_pdf_C = np.sqrt(2*np.pi)
_norm_pdf_logC = np.log(_norm_pdf_C)


def _norm_pdf(x):
    return np.exp(-x**2/2.0) / _norm_pdf_C


def _norm_logpdf(x):
    return -x**2 / 2.0 - _norm_pdf_logC


def _norm_cdf(x):
    return sc.ndtr(x)


def _norm_logcdf(x):
    return sc.log_ndtr(x)


def _norm_ppf(q):
    return sc.ndtri(q)


def _norm_sf(x):
    return _norm_cdf(-x)


def _norm_logsf(x):
    return _norm_logcdf(-x)


def _norm_isf(q):
    return -_norm_ppf(q)


class norm_gen(rv_continuous):
    r"""A normal continuous random variable.

    The location (loc) keyword specifies the mean.
    The scale (scale) keyword specifies the standard deviation.

    %(before_notes)s

    Notes
    -----
    The probability density function for `norm` is:

    .. math::

        f(x) = \frac{\exp(-x^2/2)}{\sqrt{2\pi}}

    The survival function, ``norm.sf``, is also referred to as the
    Q-function in some contexts (see, e.g.,
    `Wikipedia's <https://en.wikipedia.org/wiki/Q-function>`_ definition).

    %(after_notes)s

    %(example)s

    """
    def _rvs(self):
        return self._random_state.standard_normal(self._size)

    def _pdf(self, x):
        # norm.pdf(x) = exp(-x**2/2)/sqrt(2*pi)
        return _norm_pdf(x)

    def _logpdf(self, x):
        return _norm_logpdf(x)

    def _cdf(self, x):
        return _norm_cdf(x)

    def _logcdf(self, x):
        return _norm_logcdf(x)

    def _sf(self, x):
        return _norm_sf(x)

    def _logsf(self, x):
        return _norm_logsf(x)

    def _ppf(self, q):
        return _norm_ppf(q)

    def _isf(self, q):
        return _norm_isf(q)

    def _stats(self):
        return 0.0, 1.0, 0.0, 0.0

    def _entropy(self):
        return 0.5*(np.log(2*np.pi)+1)

    @replace_notes_in_docstring(rv_continuous, notes="""\
        This function uses explicit formulas for the maximum likelihood
        estimation of the normal distribution parameters, so the
        `optimizer` argument is ignored.\n\n""")
    def fit(self, data, **kwds):
        floc = kwds.get('floc', None)
        fscale = kwds.get('fscale', None)

        if floc is not None and fscale is not None:
            # This check is for consistency with `rv_continuous.fit`.
            # Without this check, this function would just return the
            # parameters that were given.
            raise ValueError("All parameters fixed. There is nothing to "
                             "optimize.")

        data = np.asarray(data)

        if floc is None:
            loc = data.mean()
        else:
            loc = floc

        if fscale is None:
            scale = np.sqrt(((data - loc)**2).mean())
        else:
            scale = fscale

        return loc, scale


norm = norm_gen(name='norm')


class alpha_gen(rv_continuous):
    r"""An alpha continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `alpha` is:

    .. math::

        f(x, a) = \frac{1}{x^2 \Phi(a) \sqrt{2\pi}} *
                  \exp(-\frac{1}{2} (a-1/x)^2)

    where ``Phi(alpha)`` is the normal CDF, ``x > 0``, and ``a > 0``.

    `alpha` takes ``a`` as a shape parameter.

    %(after_notes)s

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _pdf(self, x, a):
        # alpha.pdf(x, a) = 1/(x**2*Phi(a)*sqrt(2*pi)) * exp(-1/2 * (a-1/x)**2)
        return 1.0/(x**2)/_norm_cdf(a)*_norm_pdf(a-1.0/x)

    def _logpdf(self, x, a):
        return -2*np.log(x) + _norm_logpdf(a-1.0/x) - np.log(_norm_cdf(a))

    def _cdf(self, x, a):
        return _norm_cdf(a-1.0/x) / _norm_cdf(a)

    def _ppf(self, q, a):
        return 1.0/np.asarray(a-sc.ndtri(q*_norm_cdf(a)))

    def _stats(self, a):
        return [np.inf]*2 + [np.nan]*2


alpha = alpha_gen(a=0.0, name='alpha')


class anglit_gen(rv_continuous):
    r"""An anglit continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `anglit` is:

    .. math::

        f(x) = \sin(2x + \pi/2) = \cos(2x)

    for :math:`-\pi/4 \le x \le \pi/4`.

    %(after_notes)s

    %(example)s

    """
    def _pdf(self, x):
        # anglit.pdf(x) = sin(2*x + \pi/2) = cos(2*x)
        return np.cos(2*x)

    def _cdf(self, x):
        return np.sin(x+np.pi/4)**2.0

    def _ppf(self, q):
        return np.arcsin(np.sqrt(q))-np.pi/4

    def _stats(self):
        return 0.0, np.pi*np.pi/16-0.5, 0.0, -2*(np.pi**4 - 96)/(np.pi*np.pi-8)**2

    def _entropy(self):
        return 1-np.log(2)


anglit = anglit_gen(a=-np.pi/4, b=np.pi/4, name='anglit')


class arcsine_gen(rv_continuous):
    r"""An arcsine continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `arcsine` is:

    .. math::

        f(x) = \frac{1}{\pi \sqrt{x (1-x)}}

    for :math:`0 \le x \le 1`.

    %(after_notes)s

    %(example)s

    """
    def _pdf(self, x):
        # arcsine.pdf(x) = 1/(pi*sqrt(x*(1-x)))
        return 1.0/np.pi/np.sqrt(x*(1-x))

    def _cdf(self, x):
        return 2.0/np.pi*np.arcsin(np.sqrt(x))

    def _ppf(self, q):
        return np.sin(np.pi/2.0*q)**2.0

    def _stats(self):
        mu = 0.5
        mu2 = 1.0/8
        g1 = 0
        g2 = -3.0/2.0
        return mu, mu2, g1, g2

    def _entropy(self):
        return -0.24156447527049044468


arcsine = arcsine_gen(a=0.0, b=1.0, name='arcsine')


class FitDataError(ValueError):
    # This exception is raised by, for example, beta_gen.fit when both floc
    # and fscale  are fixed and there are values in the data not in the open
    # interval (floc, floc+fscale).
    def __init__(self, distr, lower, upper):
        self.args = (
            "Invalid values in `data`.  Maximum likelihood "
            "estimation with {distr!r} requires that {lower!r} < x "
            "< {upper!r} for each x in `data`.".format(
                distr=distr, lower=lower, upper=upper),
        )


class FitSolverError(RuntimeError):
    # This exception is raised by, for example, beta_gen.fit when
    # optimize.fsolve returns with ier != 1.
    def __init__(self, mesg):
        emsg = "Solver for the MLE equations failed to converge: "
        emsg += mesg.replace('\n', '')
        self.args = (emsg,)


def _beta_mle_a(a, b, n, s1):
    # The zeros of this function give the MLE for `a`, with
    # `b`, `n` and `s1` given.  `s1` is the sum of the logs of
    # the data. `n` is the number of data points.
    psiab = sc.psi(a + b)
    func = s1 - n * (-psiab + sc.psi(a))
    return func


def _beta_mle_ab(theta, n, s1, s2):
    # Zeros of this function are critical points of
    # the maximum likelihood function.  Solving this system
    # for theta (which contains a and b) gives the MLE for a and b
    # given `n`, `s1` and `s2`.  `s1` is the sum of the logs of the data,
    # and `s2` is the sum of the logs of 1 - data.  `n` is the number
    # of data points.
    a, b = theta
    psiab = sc.psi(a + b)
    func = [s1 - n * (-psiab + sc.psi(a)),
            s2 - n * (-psiab + sc.psi(b))]
    return func


class beta_gen(rv_continuous):
    r"""A beta continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `beta` is:

    .. math::

        f(x, a, b) = \frac{\gamma(a+b) x^{a-1} (1-x)^{b-1}}
                          {\gamma(a) \gamma(b)}

    for :math:`0 < x < 1`, :math:`a > 0`, :math:`b > 0`, where
    :math:`\gamma(z)` is the gamma function (`scipy.special.gamma`).

    `beta` takes :math:`a` and :math:`b` as shape parameters.

    %(after_notes)s

    %(example)s

    """
    def _rvs(self, a, b):
        return self._random_state.beta(a, b, self._size)

    def _pdf(self, x, a, b):
        #                     gamma(a+b) * x**(a-1) * (1-x)**(b-1)
        # beta.pdf(x, a, b) = ------------------------------------
        #                              gamma(a)*gamma(b)
        return np.exp(self._logpdf(x, a, b))

    def _logpdf(self, x, a, b):
        lPx = sc.xlog1py(b - 1.0, -x) + sc.xlogy(a - 1.0, x)
        lPx -= sc.betaln(a, b)
        return lPx

    def _cdf(self, x, a, b):
        return sc.btdtr(a, b, x)

    def _ppf(self, q, a, b):
        return sc.btdtri(a, b, q)

    def _stats(self, a, b):
        mn = a*1.0 / (a + b)
        var = (a*b*1.0)/(a+b+1.0)/(a+b)**2.0
        g1 = 2.0*(b-a)*np.sqrt((1.0+a+b)/(a*b)) / (2+a+b)
        g2 = 6.0*(a**3 + a**2*(1-2*b) + b**2*(1+b) - 2*a*b*(2+b))
        g2 /= a*b*(a+b+2)*(a+b+3)
        return mn, var, g1, g2

    def _fitstart(self, data):
        g1 = _skew(data)
        g2 = _kurtosis(data)

        def func(x):
            a, b = x
            sk = 2*(b-a)*np.sqrt(a + b + 1) / (a + b + 2) / np.sqrt(a*b)
            ku = a**3 - a**2*(2*b-1) + b**2*(b+1) - 2*a*b*(b+2)
            ku /= a*b*(a+b+2)*(a+b+3)
            ku *= 6
            return [sk-g1, ku-g2]
        a, b = optimize.fsolve(func, (1.0, 1.0))
        return super(beta_gen, self)._fitstart(data, args=(a, b))

    @extend_notes_in_docstring(rv_continuous, notes="""\
        In the special case where both `floc` and `fscale` are given, a
        `ValueError` is raised if any value `x` in `data` does not satisfy
        `floc < x < floc + fscale`.\n\n""")
    def fit(self, data, *args, **kwds):
        # Override rv_continuous.fit, so we can more efficiently handle the
        # case where floc and fscale are given.

        f0 = (kwds.get('f0', None) or kwds.get('fa', None) or
              kwds.get('fix_a', None))
        f1 = (kwds.get('f1', None) or kwds.get('fb', None) or
              kwds.get('fix_b', None))
        floc = kwds.get('floc', None)
        fscale = kwds.get('fscale', None)

        if floc is None or fscale is None:
            # do general fit
            return super(beta_gen, self).fit(data, *args, **kwds)

        if f0 is not None and f1 is not None:
            # This check is for consistency with `rv_continuous.fit`.
            raise ValueError("All parameters fixed. There is nothing to "
                             "optimize.")

        # Special case: loc and scale are constrained, so we are fitting
        # just the shape parameters.  This can be done much more efficiently
        # than the method used in `rv_continuous.fit`.  (See the subsection
        # "Two unknown parameters" in the section "Maximum likelihood" of
        # the Wikipedia article on the Beta distribution for the formulas.)

        # Normalize the data to the interval [0, 1].
        data = (np.ravel(data) - floc) / fscale
        if np.any(data <= 0) or np.any(data >= 1):
            raise FitDataError("beta", lower=floc, upper=floc + fscale)
        xbar = data.mean()

        if f0 is not None or f1 is not None:
            # One of the shape parameters is fixed.

            if f0 is not None:
                # The shape parameter a is fixed, so swap the parameters
                # and flip the data.  We always solve for `a`.  The result
                # will be swapped back before returning.
                b = f0
                data = 1 - data
                xbar = 1 - xbar
            else:
                b = f1

            # Initial guess for a.  Use the formula for the mean of the beta
            # distribution, E[x] = a / (a + b), to generate a reasonable
            # starting point based on the mean of the data and the given
            # value of b.
            a = b * xbar / (1 - xbar)

            # Compute the MLE for `a` by solving _beta_mle_a.
            theta, info, ier, mesg = optimize.fsolve(
                _beta_mle_a, a,
                args=(b, len(data), np.log(data).sum()),
                full_output=True
            )
            if ier != 1:
                raise FitSolverError(mesg=mesg)
            a = theta[0]

            if f0 is not None:
                # The shape parameter a was fixed, so swap back the
                # parameters.
                a, b = b, a

        else:
            # Neither of the shape parameters is fixed.

            # s1 and s2 are used in the extra arguments passed to _beta_mle_ab
            # by optimize.fsolve.
            s1 = np.log(data).sum()
            s2 = sc.log1p(-data).sum()

            # Use the "method of moments" to estimate the initial
            # guess for a and b.
            fac = xbar * (1 - xbar) / data.var(ddof=0) - 1
            a = xbar * fac
            b = (1 - xbar) * fac

            # Compute the MLE for a and b by solving _beta_mle_ab.
            theta, info, ier, mesg = optimize.fsolve(
                _beta_mle_ab, [a, b],
                args=(len(data), s1, s2),
                full_output=True
            )
            if ier != 1:
                raise FitSolverError(mesg=mesg)
            a, b = theta

        return a, b, floc, fscale


beta = beta_gen(a=0.0, b=1.0, name='beta')


class betaprime_gen(rv_continuous):
    r"""A beta prime continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `betaprime` is:

    .. math::

        f(x, a, b) = \frac{x^{a-1} (1+x)^{-a-b}}{\beta(a, b)}

    for ``x > 0``, ``a > 0``, ``b > 0``, where ``beta(a, b)`` is the beta
    function (see `scipy.special.beta`).

    `betaprime` takes ``a`` and ``b`` as shape parameters.

    %(after_notes)s

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _rvs(self, a, b):
        sz, rndm = self._size, self._random_state
        u1 = gamma.rvs(a, size=sz, random_state=rndm)
        u2 = gamma.rvs(b, size=sz, random_state=rndm)
        return u1 / u2

    def _pdf(self, x, a, b):
        # betaprime.pdf(x, a, b) = x**(a-1) * (1+x)**(-a-b) / beta(a, b)
        return np.exp(self._logpdf(x, a, b))

    def _logpdf(self, x, a, b):
        return sc.xlogy(a - 1.0, x) - sc.xlog1py(a + b, x) - sc.betaln(a, b)

    def _cdf(self, x, a, b):
        return sc.betainc(a, b, x/(1.+x))

    def _munp(self, n, a, b):
        if n == 1.0:
            return np.where(b > 1,
                            a/(b-1.0),
                            np.inf)
        elif n == 2.0:
            return np.where(b > 2,
                            a*(a+1.0)/((b-2.0)*(b-1.0)),
                            np.inf)
        elif n == 3.0:
            return np.where(b > 3,
                            a*(a+1.0)*(a+2.0)/((b-3.0)*(b-2.0)*(b-1.0)),
                            np.inf)
        elif n == 4.0:
            return np.where(b > 4,
                            (a*(a + 1.0)*(a + 2.0)*(a + 3.0) /
                             ((b - 4.0)*(b - 3.0)*(b - 2.0)*(b - 1.0))),
                            np.inf)
        else:
            raise NotImplementedError


betaprime = betaprime_gen(a=0.0, name='betaprime')


class bradford_gen(rv_continuous):
    r"""A Bradford continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `bradford` is:

    .. math::

        f(x, c) = \frac{c}{k (1+cx)}

    for :math:`0 < x < 1`, :math:`c > 0` and :math:`k = \log(1+c)`.

    `bradford` takes :math:`c` as a shape parameter.

    %(after_notes)s

    %(example)s

    """
    def _pdf(self, x, c):
        # bradford.pdf(x, c) = c / (k * (1+c*x))
        return c / (c*x + 1.0) / sc.log1p(c)

    def _cdf(self, x, c):
        return sc.log1p(c*x) / sc.log1p(c)

    def _ppf(self, q, c):
        return sc.expm1(q * sc.log1p(c)) / c

    def _stats(self, c, moments='mv'):
        k = np.log(1.0+c)
        mu = (c-k)/(c*k)
        mu2 = ((c+2.0)*k-2.0*c)/(2*c*k*k)
        g1 = None
        g2 = None
        if 's' in moments:
            g1 = np.sqrt(2)*(12*c*c-9*c*k*(c+2)+2*k*k*(c*(c+3)+3))
            g1 /= np.sqrt(c*(c*(k-2)+2*k))*(3*c*(k-2)+6*k)
        if 'k' in moments:
            g2 = (c**3*(k-3)*(k*(3*k-16)+24)+12*k*c*c*(k-4)*(k-3) +
                  6*c*k*k*(3*k-14) + 12*k**3)
            g2 /= 3*c*(c*(k-2)+2*k)**2
        return mu, mu2, g1, g2

    def _entropy(self, c):
        k = np.log(1+c)
        return k/2.0 - np.log(c/k)


bradford = bradford_gen(a=0.0, b=1.0, name='bradford')


class burr_gen(rv_continuous):
    r"""A Burr (Type III) continuous random variable.

    %(before_notes)s

    See Also
    --------
    fisk : a special case of either `burr` or ``burr12`` with ``d = 1``
    burr12 : Burr Type XII distribution

    Notes
    -----
    The probability density function for `burr` is:

    .. math::

        f(x, c, d) = c d x^{-c-1} (1+x^{-c})^{-d-1}

    for :math:`x > 0`.

    `burr` takes :math:`c` and :math:`d` as shape parameters.

    This is the PDF corresponding to the third CDF given in Burr's list;
    specifically, it is equation (11) in Burr's paper [1]_.

    %(after_notes)s

    References
    ----------
    .. [1] Burr, I. W. "Cumulative frequency functions", Annals of
       Mathematical Statistics, 13(2), pp 215-232 (1942).

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _pdf(self, x, c, d):
        # burr.pdf(x, c, d) = c * d * x**(-c-1) * (1+x**(-c))**(-d-1)
        return c * d * (x**(-c - 1.0)) * ((1 + x**(-c))**(-d - 1.0))

    def _cdf(self, x, c, d):
        return (1 + x**(-c))**(-d)

    def _ppf(self, q, c, d):
        return (q**(-1.0/d) - 1)**(-1.0/c)

    def _munp(self, n, c, d):
        nc = 1. * n / c
        return d * sc.beta(1.0 - nc, d + nc)


burr = burr_gen(a=0.0, name='burr')


class burr12_gen(rv_continuous):
    r"""A Burr (Type XII) continuous random variable.

    %(before_notes)s

    See Also
    --------
    fisk : a special case of either `burr` or ``burr12`` with ``d = 1``
    burr : Burr Type III distribution

    Notes
    -----
    The probability density function for `burr` is:

    .. math::

        f(x, c, d) = c d x^{c-1} (1+x^c)^{-d-1}

    for :math:`x > 0`.

    `burr12` takes :math:`c` and :math:`d` as shape parameters.

    This is the PDF corresponding to the twelfth CDF given in Burr's list;
    specifically, it is equation (20) in Burr's paper [1]_.

    %(after_notes)s

    The Burr type 12 distribution is also sometimes referred to as
    the Singh-Maddala distribution from NIST [2]_.

    References
    ----------
    .. [1] Burr, I. W. "Cumulative frequency functions", Annals of
       Mathematical Statistics, 13(2), pp 215-232 (1942).

    .. [2] http://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/b12pdf.htm

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _pdf(self, x, c, d):
        # burr12.pdf(x, c, d) = c * d * x**(c-1) * (1+x**(c))**(-d-1)
        return np.exp(self._logpdf(x, c, d))

    def _logpdf(self, x, c, d):
        return np.log(c) + np.log(d) + sc.xlogy(c - 1, x) + sc.xlog1py(-d-1, x**c)

    def _cdf(self, x, c, d):
        return -sc.expm1(self._logsf(x, c, d))

    def _logcdf(self, x, c, d):
        return sc.log1p(-(1 + x**c)**(-d))

    def _sf(self, x, c, d):
        return np.exp(self._logsf(x, c, d))

    def _logsf(self, x, c, d):
        return sc.xlog1py(-d, x**c)

    def _ppf(self, q, c, d):
        # The following is an implementation of
        #   ((1 - q)**(-1.0/d) - 1)**(1.0/c)
        # that does a better job handling small values of q.
        return sc.expm1(-1/d * sc.log1p(-q))**(1/c)

    def _munp(self, n, c, d):
        nc = 1. * n / c
        return d * sc.beta(1.0 + nc, d - nc)


burr12 = burr12_gen(a=0.0, name='burr12')


class fisk_gen(burr_gen):
    r"""A Fisk continuous random variable.

    The Fisk distribution is also known as the log-logistic distribution, and
    equals the Burr distribution with ``d == 1``.

    `fisk` takes :math:`c` as a shape parameter.

    %(before_notes)s

    Notes
    -----
    The probability density function for `fisk` is:

    .. math::

        f(x, c) = c x^{-c-1} (1 + x^{-c})^{-2}

    for :math:`x > 0`.

    `fisk` takes :math:`c` as a shape parameters.

    %(after_notes)s

    See Also
    --------
    burr

    %(example)s

    """
    def _pdf(self, x, c):
        # fisk.pdf(x, c) = c * x**(-c-1) * (1 + x**(-c))**(-2)
        return burr_gen._pdf(self, x, c, 1.0)

    def _cdf(self, x, c):
        return burr_gen._cdf(self, x, c, 1.0)

    def _ppf(self, x, c):
        return burr_gen._ppf(self, x, c, 1.0)

    def _munp(self, n, c):
        return burr_gen._munp(self, n, c, 1.0)

    def _entropy(self, c):
        return 2 - np.log(c)


fisk = fisk_gen(a=0.0, name='fisk')


# median = loc
class cauchy_gen(rv_continuous):
    r"""A Cauchy continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `cauchy` is:

    .. math::

        f(x) = \frac{1}{\pi (1 + x^2)}

    %(after_notes)s

    %(example)s

    """
    def _pdf(self, x):
        # cauchy.pdf(x) = 1 / (pi * (1 + x**2))
        return 1.0/np.pi/(1.0+x*x)

    def _cdf(self, x):
        return 0.5 + 1.0/np.pi*np.arctan(x)

    def _ppf(self, q):
        return np.tan(np.pi*q-np.pi/2.0)

    def _sf(self, x):
        return 0.5 - 1.0/np.pi*np.arctan(x)

    def _isf(self, q):
        return np.tan(np.pi/2.0-np.pi*q)

    def _stats(self):
        return np.nan, np.nan, np.nan, np.nan

    def _entropy(self):
        return np.log(4*np.pi)

    def _fitstart(self, data, args=None):
        # Initialize ML guesses using quartiles instead of moments.
        p25, p50, p75 = np.percentile(data, [25, 50, 75])
        return p50, (p75 - p25)/2


cauchy = cauchy_gen(name='cauchy')


class chi_gen(rv_continuous):
    r"""A chi continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `chi` is:

    .. math::

        f(x, df) = \frac{x^{df-1} \exp(-x^2/2)}{2^{df/2-1} \gamma(df/2)}

    for :math:`x > 0`.

    Special cases of `chi` are:

        - ``chi(1, loc, scale)`` is equivalent to `halfnorm`
        - ``chi(2, 0, scale)`` is equivalent to `rayleigh`
        - ``chi(3, 0, scale)`` is equivalent to `maxwell`

    `chi` takes ``df`` as a shape parameter.

    %(after_notes)s

    %(example)s

    """

    def _rvs(self, df):
        sz, rndm = self._size, self._random_state
        return np.sqrt(chi2.rvs(df, size=sz, random_state=rndm))

    def _pdf(self, x, df):
        #                   x**(df-1) * exp(-x**2/2)
        # chi.pdf(x, df) =  -------------------------
        #                   2**(df/2-1) * gamma(df/2)
        return np.exp(self._logpdf(x, df))

    def _logpdf(self, x, df):
        l = np.log(2) - .5*np.log(2)*df - sc.gammaln(.5*df)
        return l + sc.xlogy(df - 1., x) - .5*x**2

    def _cdf(self, x, df):
        return sc.gammainc(.5*df, .5*x**2)

    def _ppf(self, q, df):
        return np.sqrt(2*sc.gammaincinv(.5*df, q))

    def _stats(self, df):
        mu = np.sqrt(2)*sc.gamma(df/2.0+0.5)/sc.gamma(df/2.0)
        mu2 = df - mu*mu
        g1 = (2*mu**3.0 + mu*(1-2*df))/np.asarray(np.power(mu2, 1.5))
        g2 = 2*df*(1.0-df)-6*mu**4 + 4*mu**2 * (2*df-1)
        g2 /= np.asarray(mu2**2.0)
        return mu, mu2, g1, g2


chi = chi_gen(a=0.0, name='chi')


## Chi-squared (gamma-distributed with loc=0 and scale=2 and shape=df/2)
class chi2_gen(rv_continuous):
    r"""A chi-squared continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `chi2` is:

    .. math::

        f(x, df) = \frac{1}{(2 \gamma(df/2)} (x/2)^{df/2-1} \exp(-x/2)

    `chi2` takes ``df`` as a shape parameter.

    %(after_notes)s

    %(example)s

    """
    def _rvs(self, df):
        return self._random_state.chisquare(df, self._size)

    def _pdf(self, x, df):
        # chi2.pdf(x, df) = 1 / (2*gamma(df/2)) * (x/2)**(df/2-1) * exp(-x/2)
        return np.exp(self._logpdf(x, df))

    def _logpdf(self, x, df):
        return sc.xlogy(df/2.-1, x) - x/2. - sc.gammaln(df/2.) - (np.log(2)*df)/2.

    def _cdf(self, x, df):
        return sc.chdtr(df, x)

    def _sf(self, x, df):
        return sc.chdtrc(df, x)

    def _isf(self, p, df):
        return sc.chdtri(df, p)

    def _ppf(self, p, df):
        return self._isf(1.0-p, df)

    def _stats(self, df):
        mu = df
        mu2 = 2*df
        g1 = 2*np.sqrt(2.0/df)
        g2 = 12.0/df
        return mu, mu2, g1, g2


chi2 = chi2_gen(a=0.0, name='chi2')


class cosine_gen(rv_continuous):
    r"""A cosine continuous random variable.

    %(before_notes)s

    Notes
    -----
    The cosine distribution is an approximation to the normal distribution.
    The probability density function for `cosine` is:

    .. math::

        f(x) = \frac{1}{2\pi} (1+\cos(x))

    for :math:`-\pi \le x \le \pi`.

    %(after_notes)s

    %(example)s

    """
    def _pdf(self, x):
        # cosine.pdf(x) = 1/(2*pi) * (1+cos(x))
        return 1.0/2/np.pi*(1+np.cos(x))

    def _cdf(self, x):
        return 1.0/2/np.pi*(np.pi + x + np.sin(x))

    def _stats(self):
        return 0.0, np.pi*np.pi/3.0-2.0, 0.0, -6.0*(np.pi**4-90)/(5.0*(np.pi*np.pi-6)**2)

    def _entropy(self):
        return np.log(4*np.pi)-1.0


cosine = cosine_gen(a=-np.pi, b=np.pi, name='cosine')


class dgamma_gen(rv_continuous):
    r"""A double gamma continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `dgamma` is:

    .. math::

        f(x, a) = \frac{1}{2\gamma(a)} |x|^{a-1} \exp(-|x|)

    for :math:`a > 0`.

    `dgamma` takes :math:`a` as a shape parameter.

    %(after_notes)s

    %(example)s

    """
    def _rvs(self, a):
        sz, rndm = self._size, self._random_state
        u = rndm.random_sample(size=sz)
        gm = gamma.rvs(a, size=sz, random_state=rndm)
        return gm * np.where(u >= 0.5, 1, -1)

    def _pdf(self, x, a):
        # dgamma.pdf(x, a) = 1 / (2*gamma(a)) * abs(x)**(a-1) * exp(-abs(x))
        ax = abs(x)
        return 1.0/(2*sc.gamma(a))*ax**(a-1.0) * np.exp(-ax)

    def _logpdf(self, x, a):
        ax = abs(x)
        return sc.xlogy(a - 1.0, ax) - ax - np.log(2) - sc.gammaln(a)

    def _cdf(self, x, a):
        fac = 0.5*sc.gammainc(a, abs(x))
        return np.where(x > 0, 0.5 + fac, 0.5 - fac)

    def _sf(self, x, a):
        fac = 0.5*sc.gammainc(a, abs(x))
        return np.where(x > 0, 0.5-fac, 0.5+fac)

    def _ppf(self, q, a):
        fac = sc.gammainccinv(a, 1-abs(2*q-1))
        return np.where(q > 0.5, fac, -fac)

    def _stats(self, a):
        mu2 = a*(a+1.0)
        return 0.0, mu2, 0.0, (a+2.0)*(a+3.0)/mu2-3.0


dgamma = dgamma_gen(name='dgamma')


class dweibull_gen(rv_continuous):
    r"""A double Weibull continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `dweibull` is:

    .. math::

        f(x, c) = c / 2 |x|^{c-1} \exp(-|x|^c)

    `dweibull` takes :math:`d` as a shape parameter.

    %(after_notes)s

    %(example)s

    """
    def _rvs(self, c):
        sz, rndm = self._size, self._random_state
        u = rndm.random_sample(size=sz)
        w = weibull_min.rvs(c, size=sz, random_state=rndm)
        return w * (np.where(u >= 0.5, 1, -1))

    def _pdf(self, x, c):
        # dweibull.pdf(x, c) = c / 2 * abs(x)**(c-1) * exp(-abs(x)**c)
        ax = abs(x)
        Px = c / 2.0 * ax**(c-1.0) * np.exp(-ax**c)
        return Px

    def _logpdf(self, x, c):
        ax = abs(x)
        return np.log(c) - np.log(2.0) + sc.xlogy(c - 1.0, ax) - ax**c

    def _cdf(self, x, c):
        Cx1 = 0.5 * np.exp(-abs(x)**c)
        return np.where(x > 0, 1 - Cx1, Cx1)

    def _ppf(self, q, c):
        fac = 2. * np.where(q <= 0.5, q, 1. - q)
        fac = np.power(-np.log(fac), 1.0 / c)
        return np.where(q > 0.5, fac, -fac)

    def _munp(self, n, c):
        return (1 - (n % 2)) * sc.gamma(1.0 + 1.0 * n / c)

    # since we know that all odd moments are zeros, return them at once.
    # returning Nones from _stats makes the public stats call _munp
    # so overall we're saving one or two gamma function evaluations here.
    def _stats(self, c):
        return 0, None, 0, None


dweibull = dweibull_gen(name='dweibull')


## Exponential (gamma distributed with a=1.0, loc=loc and scale=scale)
class expon_gen(rv_continuous):
    r"""An exponential continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `expon` is:

    .. math::

        f(x) = \exp(-x)

    for :math:`x \ge 0`.

    %(after_notes)s

    A common parameterization for `expon` is in terms of the rate parameter
    ``lambda``, such that ``pdf = lambda * exp(-lambda * x)``. This
    parameterization corresponds to using ``scale = 1 / lambda``.

    %(example)s

    """
    def _rvs(self):
        return self._random_state.standard_exponential(self._size)

    def _pdf(self, x):
        # expon.pdf(x) = exp(-x)
        return np.exp(-x)

    def _logpdf(self, x):
        return -x

    def _cdf(self, x):
        return -sc.expm1(-x)

    def _ppf(self, q):
        return -sc.log1p(-q)

    def _sf(self, x):
        return np.exp(-x)

    def _logsf(self, x):
        return -x

    def _isf(self, q):
        return -np.log(q)

    def _stats(self):
        return 1.0, 1.0, 2.0, 6.0

    def _entropy(self):
        return 1.0

    @replace_notes_in_docstring(rv_continuous, notes="""\
        This function uses explicit formulas for the maximum likelihood
        estimation of the exponential distribution parameters, so the
        `optimizer`, `loc` and `scale` keyword arguments are ignored.\n\n""")
    def fit(self, data, *args, **kwds):
        if len(args) > 0:
            raise TypeError("Too many arguments.")

        floc = kwds.pop('floc', None)
        fscale = kwds.pop('fscale', None)

        # Ignore the optimizer-related keyword arguments, if given.
        kwds.pop('loc', None)
        kwds.pop('scale', None)
        kwds.pop('optimizer', None)
        if kwds:
            raise TypeError("Unknown arguments: %s." % kwds)

        if floc is not None and fscale is not None:
            # This check is for consistency with `rv_continuous.fit`.
            raise ValueError("All parameters fixed. There is nothing to "
                             "optimize.")

        data = np.asarray(data)
        data_min = data.min()
        if floc is None:
            # ML estimate of the location is the minimum of the data.
            loc = data_min
        else:
            loc = floc
            if data_min < loc:
                # There are values that are less than the specified loc.
                raise FitDataError("expon", lower=floc, upper=np.inf)

        if fscale is None:
            # ML estimate of the scale is the shifted mean.
            scale = data.mean() - loc
        else:
            scale = fscale

        # We expect the return values to be floating point, so ensure it
        # by explicitly converting to float.
        return float(loc), float(scale)


expon = expon_gen(a=0.0, name='expon')


## Exponentially Modified Normal (exponential distribution
##  convolved with a Normal).
## This is called an exponentially modified gaussian on wikipedia
class exponnorm_gen(rv_continuous):
    r"""An exponentially modified Normal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `exponnorm` is:

    .. math::

        f(x, K) = \frac{1}{2K} \exp\left(\frac{1}{2 K^2}\right) \exp(-x / K)
                  \text{erfc}\left(-\frac{x - 1/K}{\sqrt{2}}\right)

    where the shape parameter :math:`K > 0`.

    It can be thought of as the sum of a normally distributed random
    value with mean ``loc`` and sigma ``scale`` and an exponentially
    distributed random number with a pdf proportional to ``exp(-lambda * x)``
    where ``lambda = (K * scale)**(-1)``.

    %(after_notes)s

    An alternative parameterization of this distribution (for example, in
    `Wikipedia <http://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution>`_)
    involves three parameters, :math:`\mu`, :math:`\lambda` and
    :math:`\sigma`.
    In the present parameterization this corresponds to having ``loc`` and
    ``scale`` equal to :math:`\mu` and :math:`\sigma`, respectively, and
    shape parameter :math:`K = 1/(\sigma\lambda)`.

    .. versionadded:: 0.16.0

    %(example)s

    """
    def _rvs(self, K):
        expval = self._random_state.standard_exponential(self._size) * K
        gval = self._random_state.standard_normal(self._size)
        return expval + gval

    def _pdf(self, x, K):
        # exponnorm.pdf(x, K) =
        #     1/(2*K) exp(1/(2 * K**2)) exp(-x / K) * erfc-(x - 1/K) / sqrt(2))
        invK = 1.0 / K
        exparg = 0.5 * invK**2 - invK * x
        # Avoid overflows; setting np.exp(exparg) to the max float works
        #  all right here
        expval = _lazywhere(exparg < _LOGXMAX, (exparg,), np.exp, _XMAX)
        return 0.5 * invK * expval * sc.erfc(-(x - invK) / np.sqrt(2))

    def _logpdf(self, x, K):
        invK = 1.0 / K
        exparg = 0.5 * invK**2 - invK * x
        return exparg + np.log(0.5 * invK * sc.erfc(-(x - invK) / np.sqrt(2)))

    def _cdf(self, x, K):
        invK = 1.0 / K
        expval = invK * (0.5 * invK - x)
        return _norm_cdf(x) - np.exp(expval) * _norm_cdf(x - invK)

    def _sf(self, x, K):
        invK = 1.0 / K
        expval = invK * (0.5 * invK - x)
        return _norm_cdf(-x) + np.exp(expval) * _norm_cdf(x - invK)

    def _stats(self, K):
        K2 = K * K
        opK2 = 1.0 + K2
        skw = 2 * K**3 * opK2**(-1.5)
        krt = 6.0 * K2 * K2 * opK2**(-2)
        return K, opK2, skw, krt


exponnorm = exponnorm_gen(name='exponnorm')


class exponweib_gen(rv_continuous):
    r"""An exponentiated Weibull continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `exponweib` is:

    .. math::

        f(x, a, c) = a c (1-\exp(-x^c))^{a-1} \exp(-x^c) x^{c-1}

    for :math:`x > 0`, :math:`a > 0`, :math:`c > 0`.

    `exponweib` takes :math:`a` and :math:`c` as shape parameters.

    %(after_notes)s

    %(example)s

    """
    def _pdf(self, x, a, c):
        # exponweib.pdf(x, a, c) =
        #     a * c * (1-exp(-x**c))**(a-1) * exp(-x**c)*x**(c-1)
        return np.exp(self._logpdf(x, a, c))

    def _logpdf(self, x, a, c):
        negxc = -x**c
        exm1c = -sc.expm1(negxc)
        logp = (np.log(a) + np.log(c) + sc.xlogy(a - 1.0, exm1c) +
                negxc + sc.xlogy(c - 1.0, x))
        return logp

    def _cdf(self, x, a, c):
        exm1c = -sc.expm1(-x**c)
        return exm1c**a

    def _ppf(self, q, a, c):
        return (-sc.log1p(-q**(1.0/a)))**np.asarray(1.0/c)


exponweib = exponweib_gen(a=0.0, name='exponweib')


class exponpow_gen(rv_continuous):
    r"""An exponential power continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `exponpow` is:

    .. math::

        f(x, b) = b x^{b-1} \exp(1 + x^b - \exp(x^b))

    for :math:`x \ge 0`, :math:`b > 0``.  Note that this is a different
    distribution from the exponential power distribution that is also known
    under the names "generalized normal" or "generalized Gaussian".

    `exponpow` takes :math:`b` as a shape parameter.

    %(after_notes)s

    References
    ----------
    http://www.math.wm.edu/~leemis/chart/UDR/PDFs/Exponentialpower.pdf

    %(example)s

    """
    def _pdf(self, x, b):
        # exponpow.pdf(x, b) = b * x**(b-1) * exp(1 + x**b - exp(x**b))
        return np.exp(self._logpdf(x, b))

    def _logpdf(self, x, b):
        xb = x**b
        f = 1 + np.log(b) + sc.xlogy(b - 1.0, x) + xb - np.exp(xb)
        return f

    def _cdf(self, x, b):
        return -sc.expm1(-sc.expm1(x**b))

    def _sf(self, x, b):
        return np.exp(-sc.expm1(x**b))

    def _isf(self, x, b):
        return (sc.log1p(-np.log(x)))**(1./b)

    def _ppf(self, q, b):
        return pow(sc.log1p(-sc.log1p(-q)), 1.0/b)


exponpow = exponpow_gen(a=0.0, name='exponpow')


class fatiguelife_gen(rv_continuous):
    r"""A fatigue-life (Birnbaum-Saunders) continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `fatiguelife` is:

    .. math::

        f(x, c) = \frac{x+1}{ 2c\sqrt{2\pi x^3} \exp(-\frac{(x-1)^2}{2x c^2}}

    for :math:`x > 0`.

    `fatiguelife` takes :math:`c` as a shape parameter.

    %(after_notes)s

    References
    ----------
    .. [1] "Birnbaum-Saunders distribution",
           http://en.wikipedia.org/wiki/Birnbaum-Saunders_distribution

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _rvs(self, c):
        z = self._random_state.standard_normal(self._size)
        x = 0.5*c*z
        x2 = x*x
        t = 1.0 + 2*x2 + 2*x*np.sqrt(1 + x2)
        return t

    def _pdf(self, x, c):
        # fatiguelife.pdf(x, c) =
        #     (x+1) / (2*c*sqrt(2*pi*x**3)) * exp(-(x-1)**2/(2*x*c**2))
        return np.exp(self._logpdf(x, c))

    def _logpdf(self, x, c):
        return (np.log(x+1) - (x-1)**2 / (2.0*x*c**2) - np.log(2*c) -
                0.5*(np.log(2*np.pi) + 3*np.log(x)))

    def _cdf(self, x, c):
        return _norm_cdf(1.0 / c * (np.sqrt(x) - 1.0/np.sqrt(x)))

    def _ppf(self, q, c):
        tmp = c*sc.ndtri(q)
        return 0.25 * (tmp + np.sqrt(tmp**2 + 4))**2

    def _stats(self, c):
        # NB: the formula for kurtosis in wikipedia seems to have an error:
        # it's 40, not 41. At least it disagrees with the one from Wolfram
        # Alpha.  And the latter one, below, passes the tests, while the wiki
        # one doesn't So far I didn't have the guts to actually check the
        # coefficients from the expressions for the raw moments.
        c2 = c*c
        mu = c2 / 2.0 + 1.0
        den = 5.0 * c2 + 4.0
        mu2 = c2*den / 4.0
        g1 = 4 * c * (11*c2 + 6.0) / np.power(den, 1.5)
        g2 = 6 * c2 * (93*c2 + 40.0) / den**2.0
        return mu, mu2, g1, g2


fatiguelife = fatiguelife_gen(a=0.0, name='fatiguelife')


class foldcauchy_gen(rv_continuous):
    r"""A folded Cauchy continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `foldcauchy` is:

    .. math::

        f(x, c) = \frac{1}{\pi (1+(x-c)^2)} + \frac{1}{\pi (1+(x+c)^2)}

    for :math:`x \ge 0``.

    `foldcauchy` takes :math:`c` as a shape parameter.

    %(example)s

    """
    def _rvs(self, c):
        return abs(cauchy.rvs(loc=c, size=self._size,
                              random_state=self._random_state))

    def _pdf(self, x, c):
        # foldcauchy.pdf(x, c) = 1/(pi*(1+(x-c)**2)) + 1/(pi*(1+(x+c)**2))
        return 1.0/np.pi*(1.0/(1+(x-c)**2) + 1.0/(1+(x+c)**2))

    def _cdf(self, x, c):
        return 1.0/np.pi*(np.arctan(x-c) + np.arctan(x+c))

    def _stats(self, c):
        return np.inf, np.inf, np.nan, np.nan


foldcauchy = foldcauchy_gen(a=0.0, name='foldcauchy')


class f_gen(rv_continuous):
    r"""An F continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `f` is:

    .. math::

        f(x, df_1, df_2) = \frac{df_2^{df_2/2} df_1^{df_1/2} x^{df_1 / 2-1}}
                                {(df_2+df_1 x)^{(df_1+df_2)/2}
                                 B(df_1/2, df_2/2)}

    for :math:`x > 0`.

    `f` takes ``dfn`` and ``dfd`` as shape parameters.

    %(after_notes)s

    %(example)s

    """
    def _rvs(self, dfn, dfd):
        return self._random_state.f(dfn, dfd, self._size)

    def _pdf(self, x, dfn, dfd):
        #                      df2**(df2/2) * df1**(df1/2) * x**(df1/2-1)
        # F.pdf(x, df1, df2) = --------------------------------------------
        #                      (df2+df1*x)**((df1+df2)/2) * B(df1/2, df2/2)
        return np.exp(self._logpdf(x, dfn, dfd))

    def _logpdf(self, x, dfn, dfd):
        n = 1.0 * dfn
        m = 1.0 * dfd
        lPx = m/2 * np.log(m) + n/2 * np.log(n) + (n/2 - 1) * np.log(x)
        lPx -= ((n+m)/2) * np.log(m + n*x) + sc.betaln(n/2, m/2)
        return lPx

    def _cdf(self, x, dfn, dfd):
        return sc.fdtr(dfn, dfd, x)

    def _sf(self, x, dfn, dfd):
        return sc.fdtrc(dfn, dfd, x)

    def _ppf(self, q, dfn, dfd):
        return sc.fdtri(dfn, dfd, q)

    def _stats(self, dfn, dfd):
        v1, v2 = 1. * dfn, 1. * dfd
        v2_2, v2_4, v2_6, v2_8 = v2 - 2., v2 - 4., v2 - 6., v2 - 8.

        mu = _lazywhere(
            v2 > 2, (v2, v2_2),
            lambda v2, v2_2: v2 / v2_2,
            np.inf)

        mu2 = _lazywhere(
            v2 > 4, (v1, v2, v2_2, v2_4),
            lambda v1, v2, v2_2, v2_4:
            2 * v2 * v2 * (v1 + v2_2) / (v1 * v2_2**2 * v2_4),
            np.inf)

        g1 = _lazywhere(
            v2 > 6, (v1, v2_2, v2_4, v2_6),
            lambda v1, v2_2, v2_4, v2_6:
            (2 * v1 + v2_2) / v2_6 * np.sqrt(v2_4 / (v1 * (v1 + v2_2))),
            np.nan)
        g1 *= np.sqrt(8.)

        g2 = _lazywhere(
            v2 > 8, (g1, v2_6, v2_8),
            lambda g1, v2_6, v2_8: (8 + g1 * g1 * v2_6) / v2_8,
            np.nan)
        g2 *= 3. / 2.

        return mu, mu2, g1, g2


f = f_gen(a=0.0, name='f')


## Folded Normal
##   abs(Z) where (Z is normal with mu=L and std=S so that c=abs(L)/S)
##
##  note: regress docs have scale parameter correct, but first parameter
##    he gives is a shape parameter A = c * scale

##  Half-normal is folded normal with shape-parameter c=0.

class foldnorm_gen(rv_continuous):
    r"""A folded normal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `foldnorm` is:

    .. math::

        f(x, c) = \sqrt{2/\pi} cosh(c x) \exp(-\frac{x^2+c^2}{2})

    for :math:`c \ge 0`.

    `foldnorm` takes :math:`c` as a shape parameter.

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, c):
        return c >= 0

    def _rvs(self, c):
        return abs(self._random_state.standard_normal(self._size) + c)

    def _pdf(self, x, c):
        # foldnormal.pdf(x, c) = sqrt(2/pi) * cosh(c*x) * exp(-(x**2+c**2)/2)
        return _norm_pdf(x + c) + _norm_pdf(x-c)

    def _cdf(self, x, c):
        return _norm_cdf(x-c) + _norm_cdf(x+c) - 1.0

    def _stats(self, c):
        # Regina C. Elandt, Technometrics 3, 551 (1961)
        # http://www.jstor.org/stable/1266561
        #
        c2 = c*c
        expfac = np.exp(-0.5*c2) / np.sqrt(2.*np.pi)

        mu = 2.*expfac + c * sc.erf(c/np.sqrt(2))
        mu2 = c2 + 1 - mu*mu

        g1 = 2. * (mu*mu*mu - c2*mu - expfac)
        g1 /= np.power(mu2, 1.5)

        g2 = c2 * (c2 + 6.) + 3 + 8.*expfac*mu
        g2 += (2. * (c2 - 3.) - 3. * mu**2) * mu**2
        g2 = g2 / mu2**2.0 - 3.

        return mu, mu2, g1, g2


foldnorm = foldnorm_gen(a=0.0, name='foldnorm')


class weibull_min_gen(rv_continuous):
    r"""Weibull minimum continuous random variable.

    %(before_notes)s

    See Also
    --------
    weibull_max

    Notes
    -----
    The probability density function for `weibull_min` is:

    .. math::

        f(x, c) = c x^{c-1} \exp(-x^c)

    for :math:`x > 0`, :math:`c > 0`.

    `weibull_min` takes ``c`` as a shape parameter.

    %(after_notes)s

    %(example)s

    """

    def _pdf(self, x, c):
        # frechet_r.pdf(x, c) = c * x**(c-1) * exp(-x**c)
        return c*pow(x, c-1)*np.exp(-pow(x, c))

    def _logpdf(self, x, c):
        return np.log(c) + sc.xlogy(c - 1, x) - pow(x, c)

    def _cdf(self, x, c):
        return -sc.expm1(-pow(x, c))

    def _sf(self, x, c):
        return np.exp(-pow(x, c))

    def _logsf(self, x, c):
        return -pow(x, c)

    def _ppf(self, q, c):
        return pow(-sc.log1p(-q), 1.0/c)

    def _munp(self, n, c):
        return sc.gamma(1.0+n*1.0/c)

    def _entropy(self, c):
        return -_EULER / c - np.log(c) + _EULER + 1


weibull_min = weibull_min_gen(a=0.0, name='weibull_min')


class weibull_max_gen(rv_continuous):
    r"""Weibull maximum continuous random variable.

    %(before_notes)s

    See Also
    --------
    weibull_min

    Notes
    -----
    The probability density function for `weibull_max` is:

    .. math::

        f(x, c) = c (-x)^{c-1} \exp(-(-x)^c)

    for :math:`x < 0`, :math:`c > 0`.

    `weibull_max` takes ``c`` as a shape parameter.

    %(after_notes)s

    %(example)s

    """
    def _pdf(self, x, c):
        # frechet_l.pdf(x, c) = c * (-x)**(c-1) * exp(-(-x)**c)
        return c*pow(-x, c-1)*np.exp(-pow(-x, c))

    def _logpdf(self, x, c):
        return np.log(c) + sc.xlogy(c-1, -x) - pow(-x, c)

    def _cdf(self, x, c):
        return np.exp(-pow(-x, c))

    def _logcdf(self, x, c):
        return -pow(-x, c)

    def _sf(self, x, c):
        return -sc.expm1(-pow(-x, c))

    def _ppf(self, q, c):
        return -pow(-np.log(q), 1.0/c)

    def _munp(self, n, c):
        val = sc.gamma(1.0+n*1.0/c)
        if int(n) % 2:
            sgn = -1
        else:
            sgn = 1
        return sgn * val

    def _entropy(self, c):
        return -_EULER / c - np.log(c) + _EULER + 1


weibull_max = weibull_max_gen(b=0.0, name='weibull_max')

# Public methods to be deprecated in frechet_r and frechet_l:
# ['__call__', 'cdf', 'entropy', 'expect', 'fit', 'fit_loc_scale', 'freeze',
#  'interval', 'isf', 'logcdf', 'logpdf', 'logsf', 'mean', 'median', 'moment',
#  'nnlf', 'pdf', 'ppf', 'rvs', 'sf', 'stats', 'std', 'var']

_frechet_r_deprec_msg = """\
The distribution `frechet_r` is a synonym for `weibull_min`; this historical
usage is deprecated because of possible confusion with the (quite different)
Frechet distribution.  To preserve the existing behavior of the program, use
`scipy.stats.weibull_min`.  For the Frechet distribution (i.e. the Type II
extreme value distribution), use `scipy.stats.invweibull`."""

class frechet_r_gen(weibull_min_gen):

    @np.deprecate(old_name='frechet_r', message=_frechet_r_deprec_msg)
    def __call__(self, *args, **kwargs):
        return weibull_min_gen.__call__(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_r', message=_frechet_r_deprec_msg)
    def cdf(self, *args, **kwargs):
        return weibull_min_gen.cdf(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_r', message=_frechet_r_deprec_msg)
    def entropy(self, *args, **kwargs):
        return weibull_min_gen.entropy(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_r', message=_frechet_r_deprec_msg)
    def expect(self, *args, **kwargs):
        return weibull_min_gen.expect(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_r', message=_frechet_r_deprec_msg)
    def fit(self, *args, **kwargs):
        return weibull_min_gen.fit(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_r', message=_frechet_r_deprec_msg)
    def fit_loc_scale(self, *args, **kwargs):
        return weibull_min_gen.fit_loc_scale(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_r', message=_frechet_r_deprec_msg)
    def freeze(self, *args, **kwargs):
        return weibull_min_gen.freeze(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_r', message=_frechet_r_deprec_msg)
    def interval(self, *args, **kwargs):
        return weibull_min_gen.interval(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_r', message=_frechet_r_deprec_msg)
    def isf(self, *args, **kwargs):
        return weibull_min_gen.isf(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_r', message=_frechet_r_deprec_msg)
    def logcdf(self, *args, **kwargs):
        return weibull_min_gen.logcdf(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_r', message=_frechet_r_deprec_msg)
    def logpdf(self, *args, **kwargs):
        return weibull_min_gen.logpdf(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_r', message=_frechet_r_deprec_msg)
    def logsf(self, *args, **kwargs):
        return weibull_min_gen.logsf(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_r', message=_frechet_r_deprec_msg)
    def mean(self, *args, **kwargs):
        return weibull_min_gen.mean(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_r', message=_frechet_r_deprec_msg)
    def median(self, *args, **kwargs):
        return weibull_min_gen.median(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_r', message=_frechet_r_deprec_msg)
    def moment(self, *args, **kwargs):
        return weibull_min_gen.moment(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_r', message=_frechet_r_deprec_msg)
    def nnlf(self, *args, **kwargs):
        return weibull_min_gen.nnlf(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_r', message=_frechet_r_deprec_msg)
    def pdf(self, *args, **kwargs):
        return weibull_min_gen.pdf(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_r', message=_frechet_r_deprec_msg)
    def ppf(self, *args, **kwargs):
        return weibull_min_gen.ppf(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_r', message=_frechet_r_deprec_msg)
    def rvs(self, *args, **kwargs):
        return weibull_min_gen.rvs(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_r', message=_frechet_r_deprec_msg)
    def sf(self, *args, **kwargs):
        return weibull_min_gen.sf(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_r', message=_frechet_r_deprec_msg)
    def stats(self, *args, **kwargs):
        return weibull_min_gen.stats(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_r', message=_frechet_r_deprec_msg)
    def std(self, *args, **kwargs):
        return weibull_min_gen.std(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_r', message=_frechet_r_deprec_msg)
    def var(self, *args, **kwargs):
        return weibull_min_gen.var(self, *args, **kwargs)


frechet_r = frechet_r_gen(a=0.0, name='frechet_r')


_frechet_l_deprec_msg = """\
The distribution `frechet_l` is a synonym for `weibull_max`; this historical
usage is deprecated because of possible confusion with the (quite different)
Frechet distribution.  To preserve the existing behavior of the program, use
`scipy.stats.weibull_max`.  For the Frechet distribution (i.e. the Type II
extreme value distribution), use `scipy.stats.invweibull`."""

class frechet_l_gen(weibull_max_gen):

    @np.deprecate(old_name='frechet_l', message=_frechet_l_deprec_msg)
    def __call__(self, *args, **kwargs):
        return weibull_max_gen.__call__(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_l', message=_frechet_l_deprec_msg)
    def cdf(self, *args, **kwargs):
        return weibull_max_gen.cdf(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_l', message=_frechet_l_deprec_msg)
    def entropy(self, *args, **kwargs):
        return weibull_max_gen.entropy(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_l', message=_frechet_l_deprec_msg)
    def expect(self, *args, **kwargs):
        return weibull_max_gen.expect(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_l', message=_frechet_l_deprec_msg)
    def fit(self, *args, **kwargs):
        return weibull_max_gen.fit(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_l', message=_frechet_l_deprec_msg)
    def fit_loc_scale(self, *args, **kwargs):
        return weibull_max_gen.fit_loc_scale(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_l', message=_frechet_l_deprec_msg)
    def freeze(self, *args, **kwargs):
        return weibull_max_gen.freeze(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_l', message=_frechet_l_deprec_msg)
    def interval(self, *args, **kwargs):
        return weibull_max_gen.interval(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_l', message=_frechet_l_deprec_msg)
    def isf(self, *args, **kwargs):
        return weibull_max_gen.isf(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_l', message=_frechet_l_deprec_msg)
    def logcdf(self, *args, **kwargs):
        return weibull_max_gen.logcdf(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_l', message=_frechet_l_deprec_msg)
    def logpdf(self, *args, **kwargs):
        return weibull_max_gen.logpdf(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_l', message=_frechet_l_deprec_msg)
    def logsf(self, *args, **kwargs):
        return weibull_max_gen.logsf(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_l', message=_frechet_l_deprec_msg)
    def mean(self, *args, **kwargs):
        return weibull_max_gen.mean(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_l', message=_frechet_l_deprec_msg)
    def median(self, *args, **kwargs):
        return weibull_max_gen.median(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_l', message=_frechet_l_deprec_msg)
    def moment(self, *args, **kwargs):
        return weibull_max_gen.moment(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_l', message=_frechet_l_deprec_msg)
    def nnlf(self, *args, **kwargs):
        return weibull_max_gen.nnlf(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_l', message=_frechet_l_deprec_msg)
    def pdf(self, *args, **kwargs):
        return weibull_max_gen.pdf(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_l', message=_frechet_l_deprec_msg)
    def ppf(self, *args, **kwargs):
        return weibull_max_gen.ppf(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_l', message=_frechet_l_deprec_msg)
    def rvs(self, *args, **kwargs):
        return weibull_max_gen.rvs(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_l', message=_frechet_l_deprec_msg)
    def sf(self, *args, **kwargs):
        return weibull_max_gen.sf(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_l', message=_frechet_l_deprec_msg)
    def stats(self, *args, **kwargs):
        return weibull_max_gen.stats(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_l', message=_frechet_l_deprec_msg)
    def std(self, *args, **kwargs):
        return weibull_max_gen.std(self, *args, **kwargs)

    @np.deprecate(old_name='frechet_l', message=_frechet_l_deprec_msg)
    def var(self, *args, **kwargs):
        return weibull_max_gen.var(self, *args, **kwargs)


frechet_l = frechet_l_gen(b=0.0, name='frechet_l')


class genlogistic_gen(rv_continuous):
    r"""A generalized logistic continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `genlogistic` is:

    .. math::

        f(x, c) = c \frac{\exp(-x)}
                         {(1 + \exp(-x))^{c+1}}

    for :math:`x > 0`, :math:`c > 0`.

    `genlogistic` takes :math:`c` as a shape parameter.

    %(after_notes)s

    %(example)s

    """
    def _pdf(self, x, c):
        # genlogistic.pdf(x, c) = c * exp(-x) / (1 + exp(-x))**(c+1)
        return np.exp(self._logpdf(x, c))

    def _logpdf(self, x, c):
        return np.log(c) - x - (c+1.0)*sc.log1p(np.exp(-x))

    def _cdf(self, x, c):
        Cx = (1+np.exp(-x))**(-c)
        return Cx

    def _ppf(self, q, c):
        vals = -np.log(pow(q, -1.0/c)-1)
        return vals

    def _stats(self, c):
        mu = _EULER + sc.psi(c)
        mu2 = np.pi*np.pi/6.0 + sc.zeta(2, c)
        g1 = -2*sc.zeta(3, c) + 2*_ZETA3
        g1 /= np.power(mu2, 1.5)
        g2 = np.pi**4/15.0 + 6*sc.zeta(4, c)
        g2 /= mu2**2.0
        return mu, mu2, g1, g2


genlogistic = genlogistic_gen(name='genlogistic')


class genpareto_gen(rv_continuous):
    r"""A generalized Pareto continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `genpareto` is:

    .. math::

        f(x, c) = (1 + c x)^{-1 - 1/c}

    defined for :math:`x \ge 0` if :math:`c \ge 0`, and for
    :math:`0 \le x \le -1/c` if :math:`c < 0`.

    `genpareto` takes :math:`c` as a shape parameter.

    For ``c == 0``, `genpareto` reduces to the exponential
    distribution, `expon`:

    .. math::

        f(x, c=0) = \exp(-x)

    For ``c == -1``, `genpareto` is uniform on ``[0, 1]``:

    .. math::

        f(x, c=-1) = x

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, c):
        c = np.asarray(c)
        self.b = _lazywhere(c < 0, (c,),
                            lambda c: -1. / c,
                            np.inf)
        return True

    def _pdf(self, x, c):
        # genpareto.pdf(x, c) = (1 + c * x)**(-1 - 1/c)
        return np.exp(self._logpdf(x, c))

    def _logpdf(self, x, c):
        return _lazywhere((x == x) & (c != 0), (x, c),
                          lambda x, c: -sc.xlog1py(c + 1., c*x) / c,
                          -x)

    def _cdf(self, x, c):
        return -sc.inv_boxcox1p(-x, -c)

    def _sf(self, x, c):
        return sc.inv_boxcox(-x, -c)

    def _logsf(self, x, c):
        return _lazywhere((x == x) & (c != 0), (x, c),
                          lambda x, c: -sc.log1p(c*x) / c,
                          -x)

    def _ppf(self, q, c):
        return -sc.boxcox1p(-q, -c)

    def _isf(self, q, c):
        return -sc.boxcox(q, -c)

    def _munp(self, n, c):
        def __munp(n, c):
            val = 0.0
            k = np.arange(0, n + 1)
            for ki, cnk in zip(k, sc.comb(n, k)):
                val = val + cnk * (-1) ** ki / (1.0 - c * ki)
            return np.where(c * n < 1, val * (-1.0 / c) ** n, np.inf)
        return _lazywhere(c != 0, (c,),
                          lambda c: __munp(n, c),
                          sc.gamma(n + 1))

    def _entropy(self, c):
        return 1. + c


genpareto = genpareto_gen(a=0.0, name='genpareto')


class genexpon_gen(rv_continuous):
    r"""A generalized exponential continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `genexpon` is:

    .. math::

        f(x, a, b, c) = (a + b (1 - \exp(-c x)))
                        \exp(-a x - b x + \frac{b}{c}  (1-\exp(-c x)))

    for :math:`x \ge 0`, :math:`a, b, c > 0`.

    `genexpon` takes :math:`a`, :math:`b` and :math:`c` as shape parameters.

    %(after_notes)s

    References
    ----------
    H.K. Ryu, "An Extension of Marshall and Olkin's Bivariate Exponential
    Distribution", Journal of the American Statistical Association, 1993.

    N. Balakrishnan, "The Exponential Distribution: Theory, Methods and
    Applications", Asit P. Basu.

    %(example)s

    """
    def _pdf(self, x, a, b, c):
        # genexpon.pdf(x, a, b, c) = (a + b * (1 - exp(-c*x))) * \
        #                            exp(-a*x - b*x + b/c * (1-exp(-c*x)))
        return (a + b*(-sc.expm1(-c*x)))*np.exp((-a-b)*x +
                                                b*(-sc.expm1(-c*x))/c)

    def _cdf(self, x, a, b, c):
        return -sc.expm1((-a-b)*x + b*(-sc.expm1(-c*x))/c)

    def _logpdf(self, x, a, b, c):
        return np.log(a+b*(-sc.expm1(-c*x))) + (-a-b)*x+b*(-sc.expm1(-c*x))/c


genexpon = genexpon_gen(a=0.0, name='genexpon')


class genextreme_gen(rv_continuous):
    r"""A generalized extreme value continuous random variable.

    %(before_notes)s

    See Also
    --------
    gumbel_r

    Notes
    -----
    For :math:`c=0`, `genextreme` is equal to `gumbel_r`.
    The probability density function for `genextreme` is:

    .. math::

        f(x, c) = \begin{cases}
                    \exp(-\exp(-x)) \exp(-x)              &\text{for } c = 0\\
                    \exp(-(1-c x)^{1/c}) (1-c x)^{1/c-1}  &\text{for }
                                                            x \le 1/c, c > 0
                  \end{cases}


    Note that several sources and software packages use the opposite
    convention for the sign of the shape parameter :math:`c`.

    `genextreme` takes :math:`c` as a shape parameter.

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, c):
        self.b = np.where(c > 0, 1.0 / np.maximum(c, _XMIN), np.inf)
        self.a = np.where(c < 0, 1.0 / np.minimum(c, -_XMIN), -np.inf)
        return np.where(abs(c) == np.inf, 0, 1)

    def _loglogcdf(self, x, c):
        return _lazywhere((x == x) & (c != 0), (x, c),
                          lambda x, c: sc.log1p(-c*x)/c, -x)

    def _pdf(self, x, c):
        # genextreme.pdf(x, c) =
        #     exp(-exp(-x))*exp(-x),                    for c==0
        #     exp(-(1-c*x)**(1/c))*(1-c*x)**(1/c-1),    for x \le 1/c, c > 0
        return np.exp(self._logpdf(x, c))

    def _logpdf(self, x, c):
        cx = _lazywhere((x == x) & (c != 0), (x, c), lambda x, c: c*x, 0.0)
        logex2 = sc.log1p(-cx)
        logpex2 = self._loglogcdf(x, c)
        pex2 = np.exp(logpex2)
        # Handle special cases
        np.putmask(logpex2, (c == 0) & (x == -np.inf), 0.0)
        logpdf = np.where((cx == 1) | (cx == -np.inf),
                          -np.inf,
                          -pex2+logpex2-logex2)
        np.putmask(logpdf, (c == 1) & (x == 1), 0.0)
        return logpdf

    def _logcdf(self, x, c):
        return -np.exp(self._loglogcdf(x, c))

    def _cdf(self, x, c):
        return np.exp(self._logcdf(x, c))

    def _sf(self, x, c):
        return -sc.expm1(self._logcdf(x, c))

    def _ppf(self, q, c):
        x = -np.log(-np.log(q))
        return _lazywhere((x == x) & (c != 0), (x, c),
                          lambda x, c: -sc.expm1(-c * x) / c, x)

    def _isf(self, q, c):
        x = -np.log(-sc.log1p(-q))
        return _lazywhere((x == x) & (c != 0), (x, c),
                          lambda x, c: -sc.expm1(-c * x) / c, x)

    def _stats(self, c):
        g = lambda n: sc.gamma(n*c + 1)
        g1 = g(1)
        g2 = g(2)
        g3 = g(3)
        g4 = g(4)
        g2mg12 = np.where(abs(c) < 1e-7, (c*np.pi)**2.0/6.0, g2-g1**2.0)
        gam2k = np.where(abs(c) < 1e-7, np.pi**2.0/6.0,
                         sc.expm1(sc.gammaln(2.0*c+1.0)-2*sc.gammaln(c + 1.0))/c**2.0)
        eps = 1e-14
        gamk = np.where(abs(c) < eps, -_EULER, sc.expm1(sc.gammaln(c + 1))/c)

        m = np.where(c < -1.0, np.nan, -gamk)
        v = np.where(c < -0.5, np.nan, g1**2.0*gam2k)

        # skewness
        sk1 = _lazywhere(c >= -1./3,
                         (c, g1, g2, g3, g2mg12),
                         lambda c, g1, g2, g3, g2gm12:
                             np.sign(c)*(-g3 + (g2 + 2*g2mg12)*g1)/g2mg12**1.5,
                         fillvalue=np.nan)
        sk = np.where(abs(c) <= eps**0.29, 12*np.sqrt(6)*_ZETA3/np.pi**3, sk1)

        # kurtosis
        ku1 = _lazywhere(c >= -1./4,
                         (g1, g2, g3, g4, g2mg12),
                         lambda g1, g2, g3, g4, g2mg12:
                             (g4 + (-4*g3 + 3*(g2 + g2mg12)*g1)*g1)/g2mg12**2,
                         fillvalue=np.nan)
        ku = np.where(abs(c) <= (eps)**0.23, 12.0/5.0, ku1-3.0)
        return m, v, sk, ku

    def _fitstart(self, data):
        # This is better than the default shape of (1,).
        g = _skew(data)
        if g < 0:
            a = 0.5
        else:
            a = -0.5
        return super(genextreme_gen, self)._fitstart(data, args=(a,))

    def _munp(self, n, c):
        k = np.arange(0, n+1)
        vals = 1.0/c**n * np.sum(
            sc.comb(n, k) * (-1)**k * sc.gamma(c*k + 1),
            axis=0)
        return np.where(c*n > -1, vals, np.inf)

    def _entropy(self, c):
        return _EULER*(1 - c) + 1


genextreme = genextreme_gen(name='genextreme')


def _digammainv(y):
    # Inverse of the digamma function (real positive arguments only).
    # This function is used in the `fit` method of `gamma_gen`.
    # The function uses either optimize.fsolve or optimize.newton
    # to solve `sc.digamma(x) - y = 0`.  There is probably room for
    # improvement, but currently it works over a wide range of y:
    #    >>> y = 64*np.random.randn(1000000)
    #    >>> y.min(), y.max()
    #    (-311.43592651416662, 351.77388222276869)
    #    x = [_digammainv(t) for t in y]
    #    np.abs(sc.digamma(x) - y).max()
    #    1.1368683772161603e-13
    #
    _em = 0.5772156649015328606065120
    func = lambda x: sc.digamma(x) - y
    if y > -0.125:
        x0 = np.exp(y) + 0.5
        if y < 10:
            # Some experimentation shows that newton reliably converges
            # must faster than fsolve in this y range.  For larger y,
            # newton sometimes fails to converge.
            value = optimize.newton(func, x0, tol=1e-10)
            return value
    elif y > -3:
        x0 = np.exp(y/2.332) + 0.08661
    else:
        x0 = 1.0 / (-y - _em)

    value, info, ier, mesg = optimize.fsolve(func, x0, xtol=1e-11,
                                             full_output=True)
    if ier != 1:
        raise RuntimeError("_digammainv: fsolve failed, y = %r" % y)

    return value[0]


## Gamma (Use MATLAB and MATHEMATICA (b=theta=scale, a=alpha=shape) definition)

## gamma(a, loc, scale)  with a an integer is the Erlang distribution
## gamma(1, loc, scale)  is the Exponential distribution
## gamma(df/2, 0, 2) is the chi2 distribution with df degrees of freedom.

class gamma_gen(rv_continuous):
    r"""A gamma continuous random variable.

    %(before_notes)s

    See Also
    --------
    erlang, expon

    Notes
    -----
    The probability density function for `gamma` is:

    .. math::

        f(x, a) = \frac{x^{a-1} \exp(-x)}{\Gamma(a)}

    for :math:`x \ge 0`, :math:`a > 0`. Here :math:`\Gamma(a)` refers to the
    gamma function.

    `gamma` has a shape parameter `a` which needs to be set explicitly.

    When :math:`a` is an integer, `gamma` reduces to the Erlang
    distribution, and when :math:`a=1` to the exponential distribution.

    %(after_notes)s

    %(example)s

    """
    def _rvs(self, a):
        return self._random_state.standard_gamma(a, self._size)

    def _pdf(self, x, a):
        # gamma.pdf(x, a) = x**(a-1) * exp(-x) / gamma(a)
        return np.exp(self._logpdf(x, a))

    def _logpdf(self, x, a):
        return sc.xlogy(a-1.0, x) - x - sc.gammaln(a)

    def _cdf(self, x, a):
        return sc.gammainc(a, x)

    def _sf(self, x, a):
        return sc.gammaincc(a, x)

    def _ppf(self, q, a):
        return sc.gammaincinv(a, q)

    def _stats(self, a):
        return a, a, 2.0/np.sqrt(a), 6.0/a

    def _entropy(self, a):
        return sc.psi(a)*(1-a) + a + sc.gammaln(a)

    def _fitstart(self, data):
        # The skewness of the gamma distribution is `4 / np.sqrt(a)`.
        # We invert that to estimate the shape `a` using the skewness
        # of the data.  The formula is regularized with 1e-8 in the
        # denominator to allow for degenerate data where the skewness
        # is close to 0.
        a = 4 / (1e-8 + _skew(data)**2)
        return super(gamma_gen, self)._fitstart(data, args=(a,))

    @extend_notes_in_docstring(rv_continuous, notes="""\
        When the location is fixed by using the argument `floc`, this
        function uses explicit formulas or solves a simpler numerical
        problem than the full ML optimization problem.  So in that case,
        the `optimizer`, `loc` and `scale` arguments are ignored.\n\n""")
    def fit(self, data, *args, **kwds):
        f0 = (kwds.get('f0', None) or kwds.get('fa', None) or
              kwds.get('fix_a', None))
        floc = kwds.get('floc', None)
        fscale = kwds.get('fscale', None)

        if floc is None:
            # loc is not fixed.  Use the default fit method.
            return super(gamma_gen, self).fit(data, *args, **kwds)

        # Special case: loc is fixed.

        if f0 is not None and fscale is not None:
            # This check is for consistency with `rv_continuous.fit`.
            # Without this check, this function would just return the
            # parameters that were given.
            raise ValueError("All parameters fixed. There is nothing to "
                             "optimize.")

        # Fixed location is handled by shifting the data.
        data = np.asarray(data)
        if np.any(data <= floc):
            raise FitDataError("gamma", lower=floc, upper=np.inf)
        if floc != 0:
            # Don't do the subtraction in-place, because `data` might be a
            # view of the input array.
            data = data - floc
        xbar = data.mean()

        # Three cases to handle:
        # * shape and scale both free
        # * shape fixed, scale free
        # * shape free, scale fixed

        if fscale is None:
            # scale is free
            if f0 is not None:
                # shape is fixed
                a = f0
            else:
                # shape and scale are both free.
                # The MLE for the shape parameter `a` is the solution to:
                # np.log(a) - sc.digamma(a) - np.log(xbar) +
                #                             np.log(data.mean) = 0
                s = np.log(xbar) - np.log(data).mean()
                func = lambda a: np.log(a) - sc.digamma(a) - s
                aest = (3-s + np.sqrt((s-3)**2 + 24*s)) / (12*s)
                xa = aest*(1-0.4)
                xb = aest*(1+0.4)
                a = optimize.brentq(func, xa, xb, disp=0)

            # The MLE for the scale parameter is just the data mean
            # divided by the shape parameter.
            scale = xbar / a
        else:
            # scale is fixed, shape is free
            # The MLE for the shape parameter `a` is the solution to:
            # sc.digamma(a) - np.log(data).mean() + np.log(fscale) = 0
            c = np.log(data).mean() - np.log(fscale)
            a = _digammainv(c)
            scale = fscale

        return a, floc, scale


gamma = gamma_gen(a=0.0, name='gamma')


class erlang_gen(gamma_gen):
    """An Erlang continuous random variable.

    %(before_notes)s

    See Also
    --------
    gamma

    Notes
    -----
    The Erlang distribution is a special case of the Gamma distribution, with
    the shape parameter `a` an integer.  Note that this restriction is not
    enforced by `erlang`. It will, however, generate a warning the first time
    a non-integer value is used for the shape parameter.

    Refer to `gamma` for examples.

    """

    def _argcheck(self, a):
        allint = np.all(np.floor(a) == a)
        allpos = np.all(a > 0)
        if not allint:
            # An Erlang distribution shouldn't really have a non-integer
            # shape parameter, so warn the user.
            warnings.warn(
                'The shape parameter of the erlang distribution '
                'has been given a non-integer value %r.' % (a,),
                RuntimeWarning)
        return allpos

    def _fitstart(self, data):
        # Override gamma_gen_fitstart so that an integer initial value is
        # used.  (Also regularize the division, to avoid issues when
        # _skew(data) is 0 or close to 0.)
        a = int(4.0 / (1e-8 + _skew(data)**2))
        return super(gamma_gen, self)._fitstart(data, args=(a,))

    # Trivial override of the fit method, so we can monkey-patch its
    # docstring.
    def fit(self, data, *args, **kwds):
        return super(erlang_gen, self).fit(data, *args, **kwds)

    if fit.__doc__ is not None:
        fit.__doc__ = (rv_continuous.fit.__doc__ +
            """
            Notes
            -----
            The Erlang distribution is generally defined to have integer values
            for the shape parameter.  This is not enforced by the `erlang` class.
            When fitting the distribution, it will generally return a non-integer
            value for the shape parameter.  By using the keyword argument
            `f0=<integer>`, the fit method can be constrained to fit the data to
            a specific integer shape parameter.
            """)


erlang = erlang_gen(a=0.0, name='erlang')


class gengamma_gen(rv_continuous):
    r"""A generalized gamma continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `gengamma` is:

    .. math::

        f(x, a, c) = \frac{|c| x^{c a-1} \exp(-x^c)}{\gamma(a)}

    for :math:`x \ge 0`, :math:`a > 0`, and :math:`c \ne 0`.

    `gengamma` takes :math:`a` and :math:`c` as shape parameters.

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, a, c):
        return (a > 0) & (c != 0)

    def _pdf(self, x, a, c):
        # gengamma.pdf(x, a, c) = abs(c) * x**(c*a-1) * exp(-x**c) / gamma(a)
        return np.exp(self._logpdf(x, a, c))

    def _logpdf(self, x, a, c):
        return np.log(abs(c)) + sc.xlogy(c*a - 1, x) - x**c - sc.gammaln(a)

    def _cdf(self, x, a, c):
        xc = x**c
        val1 = sc.gammainc(a, xc)
        val2 = sc.gammaincc(a, xc)
        return np.where(c > 0, val1, val2)

    def _sf(self, x, a, c):
        xc = x**c
        val1 = sc.gammainc(a, xc)
        val2 = sc.gammaincc(a, xc)
        return np.where(c > 0, val2, val1)

    def _ppf(self, q, a, c):
        val1 = sc.gammaincinv(a, q)
        val2 = sc.gammainccinv(a, q)
        return np.where(c > 0, val1, val2)**(1.0/c)

    def _isf(self, q, a, c):
        val1 = sc.gammaincinv(a, q)
        val2 = sc.gammainccinv(a, q)
        return np.where(c > 0, val2, val1)**(1.0/c)

    def _munp(self, n, a, c):
        # Pochhammer symbol: sc.pocha,n) = gamma(a+n)/gamma(a)
        return sc.poch(a, n*1.0/c)

    def _entropy(self, a, c):
        val = sc.psi(a)
        return a*(1-val) + 1.0/c*val + sc.gammaln(a) - np.log(abs(c))


gengamma = gengamma_gen(a=0.0, name='gengamma')


class genhalflogistic_gen(rv_continuous):
    r"""A generalized half-logistic continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `genhalflogistic` is:

    .. math::

        f(x, c) = \frac{2 (1 - c x)^{1/(c-1)}}{[1 + (1 - c x)^{1/c}]^2}

    for :math:`0 \le x \le 1/c`, and :math:`c > 0`.

    `genhalflogistic` takes :math:`c` as a shape parameter.

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, c):
        self.b = 1.0 / c
        return c > 0

    def _pdf(self, x, c):
        # genhalflogistic.pdf(x, c) =
        #    2 * (1-c*x)**(1/c-1) / (1+(1-c*x)**(1/c))**2
        limit = 1.0/c
        tmp = np.asarray(1-c*x)
        tmp0 = tmp**(limit-1)
        tmp2 = tmp0*tmp
        return 2*tmp0 / (1+tmp2)**2

    def _cdf(self, x, c):
        limit = 1.0/c
        tmp = np.asarray(1-c*x)
        tmp2 = tmp**(limit)
        return (1.0-tmp2) / (1+tmp2)

    def _ppf(self, q, c):
        return 1.0/c*(1-((1.0-q)/(1.0+q))**c)

    def _entropy(self, c):
        return 2 - (2*c+1)*np.log(2)


genhalflogistic = genhalflogistic_gen(a=0.0, name='genhalflogistic')


class gompertz_gen(rv_continuous):
    r"""A Gompertz (or truncated Gumbel) continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `gompertz` is:

    .. math::

        f(x, c) = c \exp(x) \exp(-c (e^x-1))

    for :math:`x \ge 0`, :math:`c > 0`.

    `gompertz` takes :math:`c` as a shape parameter.

    %(after_notes)s

    %(example)s

    """
    def _pdf(self, x, c):
        # gompertz.pdf(x, c) = c * exp(x) * exp(-c*(exp(x)-1))
        return np.exp(self._logpdf(x, c))

    def _logpdf(self, x, c):
        return np.log(c) + x - c * sc.expm1(x)

    def _cdf(self, x, c):
        return -sc.expm1(-c * sc.expm1(x))

    def _ppf(self, q, c):
        return sc.log1p(-1.0 / c * sc.log1p(-q))

    def _entropy(self, c):
        return 1.0 - np.log(c) - np.exp(c)*sc.expn(1, c)


gompertz = gompertz_gen(a=0.0, name='gompertz')


class gumbel_r_gen(rv_continuous):
    r"""A right-skewed Gumbel continuous random variable.

    %(before_notes)s

    See Also
    --------
    gumbel_l, gompertz, genextreme

    Notes
    -----
    The probability density function for `gumbel_r` is:

    .. math::

        f(x) = \exp(-(x + e^{-x}))

    The Gumbel distribution is sometimes referred to as a type I Fisher-Tippett
    distribution.  It is also related to the extreme value distribution,
    log-Weibull and Gompertz distributions.

    %(after_notes)s

    %(example)s

    """
    def _pdf(self, x):
        # gumbel_r.pdf(x) = exp(-(x + exp(-x)))
        return np.exp(self._logpdf(x))

    def _logpdf(self, x):
        return -x - np.exp(-x)

    def _cdf(self, x):
        return np.exp(-np.exp(-x))

    def _logcdf(self, x):
        return -np.exp(-x)

    def _ppf(self, q):
        return -np.log(-np.log(q))

    def _stats(self):
        return _EULER, np.pi*np.pi/6.0, 12*np.sqrt(6)/np.pi**3 * _ZETA3, 12.0/5

    def _entropy(self):
        # http://en.wikipedia.org/wiki/Gumbel_distribution
        return _EULER + 1.


gumbel_r = gumbel_r_gen(name='gumbel_r')


class gumbel_l_gen(rv_continuous):
    r"""A left-skewed Gumbel continuous random variable.

    %(before_notes)s

    See Also
    --------
    gumbel_r, gompertz, genextreme

    Notes
    -----
    The probability density function for `gumbel_l` is:

    .. math::

        f(x) = \exp(x - e^x)

    The Gumbel distribution is sometimes referred to as a type I Fisher-Tippett
    distribution.  It is also related to the extreme value distribution,
    log-Weibull and Gompertz distributions.

    %(after_notes)s

    %(example)s

    """
    def _pdf(self, x):
        # gumbel_l.pdf(x) = exp(x - exp(x))
        return np.exp(self._logpdf(x))

    def _logpdf(self, x):
        return x - np.exp(x)

    def _cdf(self, x):
        return -sc.expm1(-np.exp(x))

    def _ppf(self, q):
        return np.log(-sc.log1p(-q))

    def _logsf(self, x):
        return -np.exp(x)

    def _sf(self, x):
        return np.exp(-np.exp(x))

    def _isf(self, x):
        return np.log(-np.log(x))

    def _stats(self):
        return -_EULER, np.pi*np.pi/6.0, \
               -12*np.sqrt(6)/np.pi**3 * _ZETA3, 12.0/5

    def _entropy(self):
        return _EULER + 1.


gumbel_l = gumbel_l_gen(name='gumbel_l')


class halfcauchy_gen(rv_continuous):
    r"""A Half-Cauchy continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `halfcauchy` is:

    .. math::

        f(x) = \frac{2}{\pi (1 + x^2)}

    for :math:`x \ge 0`.

    %(after_notes)s

    %(example)s

    """
    def _pdf(self, x):
        # halfcauchy.pdf(x) = 2 / (pi * (1 + x**2))
        return 2.0/np.pi/(1.0+x*x)

    def _logpdf(self, x):
        return np.log(2.0/np.pi) - sc.log1p(x*x)

    def _cdf(self, x):
        return 2.0/np.pi*np.arctan(x)

    def _ppf(self, q):
        return np.tan(np.pi/2*q)

    def _stats(self):
        return np.inf, np.inf, np.nan, np.nan

    def _entropy(self):
        return np.log(2*np.pi)


halfcauchy = halfcauchy_gen(a=0.0, name='halfcauchy')


class halflogistic_gen(rv_continuous):
    r"""A half-logistic continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `halflogistic` is:

    .. math::

        f(x) = \frac{ 2 e^{-x} }{ (1+e^{-x})^2 } = \frac{1}{2} sech(x/2)^2

    for :math:`x \ge 0`.

    %(after_notes)s

    %(example)s

    """
    def _pdf(self, x):
        # halflogistic.pdf(x) = 2 * exp(-x) / (1+exp(-x))**2
        #                     = 1/2 * sech(x/2)**2
        return np.exp(self._logpdf(x))

    def _logpdf(self, x):
        return np.log(2) - x - 2. * sc.log1p(np.exp(-x))

    def _cdf(self, x):
        return np.tanh(x/2.0)

    def _ppf(self, q):
        return 2*np.arctanh(q)

    def _munp(self, n):
        if n == 1:
            return 2*np.log(2)
        if n == 2:
            return np.pi*np.pi/3.0
        if n == 3:
            return 9*_ZETA3
        if n == 4:
            return 7*np.pi**4 / 15.0
        return 2*(1-pow(2.0, 1-n))*sc.gamma(n+1)*sc.zeta(n, 1)

    def _entropy(self):
        return 2-np.log(2)


halflogistic = halflogistic_gen(a=0.0, name='halflogistic')


class halfnorm_gen(rv_continuous):
    r"""A half-normal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `halfnorm` is:

    .. math::

        f(x) = \sqrt{2/\pi} e^{-\frac{x^2}{2}}

    for :math:`x > 0`.

    `halfnorm` is a special case of :math`\chi` with ``df == 1``.

    %(after_notes)s

    %(example)s

    """
    def _rvs(self):
        return abs(self._random_state.standard_normal(size=self._size))

    def _pdf(self, x):
        # halfnorm.pdf(x) = sqrt(2/pi) * exp(-x**2/2)
        return np.sqrt(2.0/np.pi)*np.exp(-x*x/2.0)

    def _logpdf(self, x):
        return 0.5 * np.log(2.0/np.pi) - x*x/2.0

    def _cdf(self, x):
        return _norm_cdf(x)*2-1.0

    def _ppf(self, q):
        return sc.ndtri((1+q)/2.0)

    def _stats(self):
        return (np.sqrt(2.0/np.pi),
                1-2.0/np.pi,
                np.sqrt(2)*(4-np.pi)/(np.pi-2)**1.5,
                8*(np.pi-3)/(np.pi-2)**2)

    def _entropy(self):
        return 0.5*np.log(np.pi/2.0)+0.5


halfnorm = halfnorm_gen(a=0.0, name='halfnorm')


class hypsecant_gen(rv_continuous):
    r"""A hyperbolic secant continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `hypsecant` is:

    .. math::

        f(x) = \frac{1}{\pi} sech(x)

    %(after_notes)s

    %(example)s

    """
    def _pdf(self, x):
        # hypsecant.pdf(x) = 1/pi * sech(x)
        return 1.0/(np.pi*np.cosh(x))

    def _cdf(self, x):
        return 2.0/np.pi*np.arctan(np.exp(x))

    def _ppf(self, q):
        return np.log(np.tan(np.pi*q/2.0))

    def _stats(self):
        return 0, np.pi*np.pi/4, 0, 2

    def _entropy(self):
        return np.log(2*np.pi)


hypsecant = hypsecant_gen(name='hypsecant')


class gausshyper_gen(rv_continuous):
    r"""A Gauss hypergeometric continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `gausshyper` is:

    .. math::

        f(x, a, b, c, z) = C x^{a-1} (1-x)^{b-1} (1+zx)^{-c}

    for :math:`0 \le x \le 1`, :math:`a > 0`, :math:`b > 0`, and
    :math:`C = \frac{1}{B(a, b) F[2, 1](c, a; a+b; -z)}`

    `gausshyper` takes :math:`a`, :math:`b`, :math:`c` and :math:`z` as shape
    parameters.

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, a, b, c, z):
        return (a > 0) & (b > 0) & (c == c) & (z == z)

    def _pdf(self, x, a, b, c, z):
        # gausshyper.pdf(x, a, b, c, z) =
        #   C * x**(a-1) * (1-x)**(b-1) * (1+z*x)**(-c)
        Cinv = sc.gamma(a)*sc.gamma(b)/sc.gamma(a+b)*sc.hyp2f1(c, a, a+b, -z)
        return 1.0/Cinv * x**(a-1.0) * (1.0-x)**(b-1.0) / (1.0+z*x)**c

    def _munp(self, n, a, b, c, z):
        fac = sc.beta(n+a, b) / sc.beta(a, b)
        num = sc.hyp2f1(c, a+n, a+b+n, -z)
        den = sc.hyp2f1(c, a, a+b, -z)
        return fac*num / den


gausshyper = gausshyper_gen(a=0.0, b=1.0, name='gausshyper')


class invgamma_gen(rv_continuous):
    r"""An inverted gamma continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `invgamma` is:

    .. math::

        f(x, a) = \frac{x^{-a-1}}{\gamma(a)} \exp(-\frac{1}{x})

    for :math:`x > 0`, :math:`a > 0`.

    `invgamma` takes :math:`a` as a shape parameter.

    `invgamma` is a special case of `gengamma` with ``c == -1``.

    %(after_notes)s

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _pdf(self, x, a):
        # invgamma.pdf(x, a) = x**(-a-1) / gamma(a) * exp(-1/x)
        return np.exp(self._logpdf(x, a))

    def _logpdf(self, x, a):
        return -(a+1) * np.log(x) - sc.gammaln(a) - 1.0/x

    def _cdf(self, x, a):
        return sc.gammaincc(a, 1.0 / x)

    def _ppf(self, q, a):
        return 1.0 / sc.gammainccinv(a, q)

    def _sf(self, x, a):
        return sc.gammainc(a, 1.0 / x)

    def _isf(self, q, a):
        return 1.0 / sc.gammaincinv(a, q)

    def _stats(self, a, moments='mvsk'):
        m1 = _lazywhere(a > 1, (a,), lambda x: 1. / (x - 1.), np.inf)
        m2 = _lazywhere(a > 2, (a,), lambda x: 1. / (x - 1.)**2 / (x - 2.),
                        np.inf)

        g1, g2 = None, None
        if 's' in moments:
            g1 = _lazywhere(
                a > 3, (a,),
                lambda x: 4. * np.sqrt(x - 2.) / (x - 3.), np.nan)
        if 'k' in moments:
            g2 = _lazywhere(
                a > 4, (a,),
                lambda x: 6. * (5. * x - 11.) / (x - 3.) / (x - 4.), np.nan)
        return m1, m2, g1, g2

    def _entropy(self, a):
        return a - (a+1.0) * sc.psi(a) + sc.gammaln(a)


invgamma = invgamma_gen(a=0.0, name='invgamma')


# scale is gamma from DATAPLOT and B from Regress
class invgauss_gen(rv_continuous):
    r"""An inverse Gaussian continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `invgauss` is:

    .. math::

        f(x, \mu) = \frac{1}{\sqrt{2 \pi x^3}}
                    \exp(-\frac{(x-\mu)^2}{2 x \mu^2})

    for :math:`x > 0`.

    `invgauss` takes :math:`\mu` as a shape parameter.

    %(after_notes)s

    When :math:`\mu` is too small, evaluating the cumulative distribution
    function will be inaccurate due to ``cdf(mu -> 0) = inf * 0``.
    NaNs are returned for :math:`\mu \le 0.0028`.

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _rvs(self, mu):
        return self._random_state.wald(mu, 1.0, size=self._size)

    def _pdf(self, x, mu):
        # invgauss.pdf(x, mu) =
        #                  1 / sqrt(2*pi*x**3) * exp(-(x-mu)**2/(2*x*mu**2))
        return 1.0/np.sqrt(2*np.pi*x**3.0)*np.exp(-1.0/(2*x)*((x-mu)/mu)**2)

    def _logpdf(self, x, mu):
        return -0.5*np.log(2*np.pi) - 1.5*np.log(x) - ((x-mu)/mu)**2/(2*x)

    def _cdf(self, x, mu):
        fac = np.sqrt(1.0/x)
        # Numerical accuracy for small `mu` is bad.  See #869.
        C1 = _norm_cdf(fac*(x-mu)/mu)
        C1 += np.exp(1.0/mu) * _norm_cdf(-fac*(x+mu)/mu) * np.exp(1.0/mu)
        return C1

    def _stats(self, mu):
        return mu, mu**3.0, 3*np.sqrt(mu), 15*mu


invgauss = invgauss_gen(a=0.0, name='invgauss')


class norminvgauss_gen(rv_continuous):
    r"""A Normal Inverse Gaussian continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `norminvgauss` is:

    .. math::

        f(x; a, b) = (a \exp(\sqrt{a^2 - b^2} + b x)) /
                     (\pi \sqrt{1 + x^2} \, K_1(a * \sqrt{1 + x^2}))

    where `x` is a real number, the parameter `a` is the tail heaviness
    and `b` is the asymmetry parameter satisfying `a > 0` and `abs(b) <= a`.
    `K_1` is the modified Bessel function of second kind (`scipy.special.k1`).

    %(after_notes)s

    A normal inverse Gaussian random variable with parameters `a` and `b` can
    be expressed  as `X = b * V + sqrt(V) * X` where `X` is `norm(0,1)`
    and `V` is `invgauss(mu=1/sqrt(a**2 - b**2))`. This representation is used
    to generate random variates.

    References
    ----------
    O. Barndorff-Nielsen, "Hyperbolic Distributions and Distributions on
    Hyperbolae", Scandinavian Journal of Statistics, Vol. 5(3),
    pp. 151-157, 1978.

    O. Barndorff-Nielsen, "Normal Inverse Gaussian Distributions and Stochastic
    Volatility Modelling", Scandinavian Journal of Statistics, Vol. 24,
    pp. 1–13, 1997.

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _argcheck(self, a, b):
        return (a > 0) & (np.absolute(b) < a)

    def _pdf(self, x, a, b):
        gamma = np.sqrt(a**2 - b**2)
        fac1 = a / np.pi * np.exp(gamma)
        sq = np.hypot(1, x)  # reduce overflows
        return fac1 * sc.k1e(a * sq) * np.exp(b*x - a*sq) / sq

    def _rvs(self, a, b):
        # note: X = b * V + sqrt(V) * X is norminvgaus(a,b) if X is standard
        # normal and V is invgauss(mu=1/sqrt(a**2 - b**2))
        gamma = np.sqrt(a**2 - b**2)
        sz, rndm = self._size, self._random_state
        ig = invgauss.rvs(mu=1/gamma, size=sz, random_state=rndm)
        return b * ig + np.sqrt(ig) * norm.rvs(size=sz, random_state=rndm)

    def _stats(self, a, b):
        gamma = np.sqrt(a**2 - b**2)
        mean = b / gamma
        variance = a**2 / gamma**3
        skewness = 3.0 * b / (a * np.sqrt(gamma))
        kurtosis = 3.0 * (1 + 4 * b**2 / a**2) / gamma
        return mean, variance, skewness, kurtosis


norminvgauss = norminvgauss_gen(name="norminvgauss")


class invweibull_gen(rv_continuous):
    r"""An inverted Weibull continuous random variable.

    This distribution is also known as the Fréchet distribution or the
    type II extreme value distribution.

    %(before_notes)s

    Notes
    -----
    The probability density function for `invweibull` is:

    .. math::

        f(x, c) = c x^{-c-1} \exp(-x^{-c})

    for :math:`x > 0``, :math:`c > 0``.

    `invweibull` takes :math:`c`` as a shape parameter.

    %(after_notes)s

    References
    ----------
    F.R.S. de Gusmao, E.M.M Ortega and G.M. Cordeiro, "The generalized inverse
    Weibull distribution", Stat. Papers, vol. 52, pp. 591-619, 2011.

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _pdf(self, x, c):
        # invweibull.pdf(x, c) = c * x**(-c-1) * exp(-x**(-c))
        xc1 = np.power(x, -c - 1.0)
        xc2 = np.power(x, -c)
        xc2 = np.exp(-xc2)
        return c * xc1 * xc2

    def _cdf(self, x, c):
        xc1 = np.power(x, -c)
        return np.exp(-xc1)

    def _ppf(self, q, c):
        return np.power(-np.log(q), -1.0/c)

    def _munp(self, n, c):
        return sc.gamma(1 - n / c)

    def _entropy(self, c):
        return 1+_EULER + _EULER / c - np.log(c)


invweibull = invweibull_gen(a=0, name='invweibull')


class johnsonsb_gen(rv_continuous):
    r"""A Johnson SB continuous random variable.

    %(before_notes)s

    See Also
    --------
    johnsonsu

    Notes
    -----
    The probability density function for `johnsonsb` is:

    .. math::

        f(x, a, b) = \frac{b}{x(1-x)}  \phi(a + b \log \frac{x}{1-x} )

    for :math:`0 < x < 1` and :math:`a, b > 0`, and :math:`\phi` is the normal
    pdf.

    `johnsonsb` takes :math:`a` and :math:`b` as shape parameters.

    %(after_notes)s

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _argcheck(self, a, b):
        return (b > 0) & (a == a)

    def _pdf(self, x, a, b):
        # johnsonsb.pdf(x, a, b) = b / (x*(1-x)) * phi(a + b * log(x/(1-x)))
        trm = _norm_pdf(a + b*np.log(x/(1.0-x)))
        return b*1.0/(x*(1-x))*trm

    def _cdf(self, x, a, b):
        return _norm_cdf(a + b*np.log(x/(1.0-x)))

    def _ppf(self, q, a, b):
        return 1.0 / (1 + np.exp(-1.0 / b * (_norm_ppf(q) - a)))


johnsonsb = johnsonsb_gen(a=0.0, b=1.0, name='johnsonsb')


class johnsonsu_gen(rv_continuous):
    r"""A Johnson SU continuous random variable.

    %(before_notes)s

    See Also
    --------
    johnsonsb

    Notes
    -----
    The probability density function for `johnsonsu` is:

    .. math::

        f(x, a, b) = \frac{b}{\sqrt{x^2 + 1}}
                     \phi(a + b \log(x + \sqrt{x^2 + 1}))

    for all :math:`x, a, b > 0`, and :math:`\phi` is the normal pdf.

    `johnsonsu` takes :math:`a` and :math:`b` as shape parameters.

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, a, b):
        return (b > 0) & (a == a)

    def _pdf(self, x, a, b):
        # johnsonsu.pdf(x, a, b) = b / sqrt(x**2 + 1) *
        #                          phi(a + b * log(x + sqrt(x**2 + 1)))
        x2 = x*x
        trm = _norm_pdf(a + b * np.log(x + np.sqrt(x2+1)))
        return b*1.0/np.sqrt(x2+1.0)*trm

    def _cdf(self, x, a, b):
        return _norm_cdf(a + b * np.log(x + np.sqrt(x*x + 1)))

    def _ppf(self, q, a, b):
        return np.sinh((_norm_ppf(q) - a) / b)


johnsonsu = johnsonsu_gen(name='johnsonsu')


class laplace_gen(rv_continuous):
    r"""A Laplace continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `laplace` is:

    .. math::

        f(x) = \frac{1}{2} \exp(-|x|)

    %(after_notes)s

    %(example)s

    """
    def _rvs(self):
        return self._random_state.laplace(0, 1, size=self._size)

    def _pdf(self, x):
        # laplace.pdf(x) = 1/2 * exp(-abs(x))
        return 0.5*np.exp(-abs(x))

    def _cdf(self, x):
        return np.where(x > 0, 1.0-0.5*np.exp(-x), 0.5*np.exp(x))

    def _ppf(self, q):
        return np.where(q > 0.5, -np.log(2*(1-q)), np.log(2*q))

    def _stats(self):
        return 0, 2, 0, 3

    def _entropy(self):
        return np.log(2)+1


laplace = laplace_gen(name='laplace')


class levy_gen(rv_continuous):
    r"""A Levy continuous random variable.

    %(before_notes)s

    See Also
    --------
    levy_stable, levy_l

    Notes
    -----
    The probability density function for `levy` is:

    .. math::

        f(x) = \frac{1}{x \sqrt{2\pi x}) \exp(-\frac{1}{2x}}

    for :math:`x > 0`.

    This is the same as the Levy-stable distribution with :math:`a=1/2` and
    :math:`b=1`.

    %(after_notes)s

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _pdf(self, x):
        # levy.pdf(x) = 1 / (x * sqrt(2*pi*x)) * exp(-1/(2*x))
        return 1 / np.sqrt(2*np.pi*x) / x * np.exp(-1/(2*x))

    def _cdf(self, x):
        # Equivalent to 2*norm.sf(np.sqrt(1/x))
        return sc.erfc(np.sqrt(0.5 / x))

    def _ppf(self, q):
        # Equivalent to 1.0/(norm.isf(q/2)**2) or 0.5/(erfcinv(q)**2)
        val = -sc.ndtri(q/2)
        return 1.0 / (val * val)

    def _stats(self):
        return np.inf, np.inf, np.nan, np.nan


levy = levy_gen(a=0.0, name="levy")


class levy_l_gen(rv_continuous):
    r"""A left-skewed Levy continuous random variable.

    %(before_notes)s

    See Also
    --------
    levy, levy_stable

    Notes
    -----
    The probability density function for `levy_l` is:

    .. math::

        f(x) = \frac{1}{|x| \sqrt{2\pi |x|}} \exp(-\frac{1}{2 |x|})

    for :math:`x < 0`.

    This is the same as the Levy-stable distribution with :math:`a=1/2` and
    :math:`b=-1`.

    %(after_notes)s

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _pdf(self, x):
        # levy_l.pdf(x) = 1 / (abs(x) * sqrt(2*pi*abs(x))) * exp(-1/(2*abs(x)))
        ax = abs(x)
        return 1/np.sqrt(2*np.pi*ax)/ax*np.exp(-1/(2*ax))

    def _cdf(self, x):
        ax = abs(x)
        return 2 * _norm_cdf(1 / np.sqrt(ax)) - 1

    def _ppf(self, q):
        val = _norm_ppf((q + 1.0) / 2)
        return -1.0 / (val * val)

    def _stats(self):
        return np.inf, np.inf, np.nan, np.nan


levy_l = levy_l_gen(b=0.0, name="levy_l")


class levy_stable_gen(rv_continuous):
    r"""A Levy-stable continuous random variable.

    %(before_notes)s

    See Also
    --------
    levy, levy_l

    Notes
    -----
    Levy-stable distribution (only random variates available -- ignore other
    docs)

    """

    def _rvs(self, alpha, beta):

        def alpha1func(alpha, beta, TH, aTH, bTH, cosTH, tanTH, W):
            return (2/np.pi*(np.pi/2 + bTH)*tanTH -
                    beta*np.log((np.pi/2*W*cosTH)/(np.pi/2 + bTH)))

        def beta0func(alpha, beta, TH, aTH, bTH, cosTH, tanTH, W):
            return (W/(cosTH/np.tan(aTH) + np.sin(TH)) *
                    ((np.cos(aTH) + np.sin(aTH)*tanTH)/W)**(1.0/alpha))

        def otherwise(alpha, beta, TH, aTH, bTH, cosTH, tanTH, W):
            # alpha is not 1 and beta is not 0
            val0 = beta*np.tan(np.pi*alpha/2)
            th0 = np.arctan(val0)/alpha
            val3 = W/(cosTH/np.tan(alpha*(th0 + TH)) + np.sin(TH))
            res3 = val3*((np.cos(aTH) + np.sin(aTH)*tanTH -
                          val0*(np.sin(aTH) - np.cos(aTH)*tanTH))/W)**(1.0/alpha)
            return res3

        def alphanot1func(alpha, beta, TH, aTH, bTH, cosTH, tanTH, W):
            res = _lazywhere(beta == 0,
                             (alpha, beta, TH, aTH, bTH, cosTH, tanTH, W),
                             beta0func, f2=otherwise)
            return res

        sz = self._size
        alpha = broadcast_to(alpha, sz)
        beta = broadcast_to(beta, sz)
        TH = uniform.rvs(loc=-np.pi/2.0, scale=np.pi, size=sz,
                         random_state=self._random_state)
        W = expon.rvs(size=sz, random_state=self._random_state)
        aTH = alpha*TH
        bTH = beta*TH
        cosTH = np.cos(TH)
        tanTH = np.tan(TH)
        res = _lazywhere(alpha == 1,
                         (alpha, beta, TH, aTH, bTH, cosTH, tanTH, W),
                         alpha1func, f2=alphanot1func)
        return res

    def _argcheck(self, alpha, beta):
        return (alpha > 0) & (alpha <= 2) & (beta <= 1) & (beta >= -1)

    def _pdf(self, x, alpha, beta):
        raise NotImplementedError


levy_stable = levy_stable_gen(name='levy_stable')


class logistic_gen(rv_continuous):
    r"""A logistic (or Sech-squared) continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `logistic` is:

    .. math::

        f(x) = \frac{\exp(-x)}
                    {(1+exp(-x))^2}

    `logistic` is a special case of `genlogistic` with ``c == 1``.

    %(after_notes)s

    %(example)s

    """
    def _rvs(self):
        return self._random_state.logistic(size=self._size)

    def _pdf(self, x):
        # logistic.pdf(x) = exp(-x) / (1+exp(-x))**2
        return np.exp(self._logpdf(x))

    def _logpdf(self, x):
        return -x - 2. * sc.log1p(np.exp(-x))

    def _cdf(self, x):
        return sc.expit(x)

    def _ppf(self, q):
        return sc.logit(q)

    def _sf(self, x):
        return sc.expit(-x)

    def _isf(self, q):
        return -sc.logit(q)

    def _stats(self):
        return 0, np.pi*np.pi/3.0, 0, 6.0/5.0

    def _entropy(self):
        # http://en.wikipedia.org/wiki/Logistic_distribution
        return 2.0


logistic = logistic_gen(name='logistic')


class loggamma_gen(rv_continuous):
    r"""A log gamma continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `loggamma` is:

    .. math::

        f(x, c) = \frac{\exp(c x - \exp(x))}
                       {\gamma(c)}

    for all :math:`x, c > 0`.

    `loggamma` takes :math:`c` as a shape parameter.

    %(after_notes)s

    %(example)s

    """
    def _rvs(self, c):
        return np.log(self._random_state.gamma(c, size=self._size))

    def _pdf(self, x, c):
        # loggamma.pdf(x, c) = exp(c*x-exp(x)) / gamma(c)
        return np.exp(c*x-np.exp(x)-sc.gammaln(c))

    def _cdf(self, x, c):
        return sc.gammainc(c, np.exp(x))

    def _ppf(self, q, c):
        return np.log(sc.gammaincinv(c, q))

    def _stats(self, c):
        # See, for example, "A Statistical Study of Log-Gamma Distribution", by
        # Ping Shing Chan (thesis, McMaster University, 1993).
        mean = sc.digamma(c)
        var = sc.polygamma(1, c)
        skewness = sc.polygamma(2, c) / np.power(var, 1.5)
        excess_kurtosis = sc.polygamma(3, c) / (var*var)
        return mean, var, skewness, excess_kurtosis


loggamma = loggamma_gen(name='loggamma')


class loglaplace_gen(rv_continuous):
    r"""A log-Laplace continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `loglaplace` is:

    .. math::

        f(x, c) = \begin{cases}\frac{c}{2} x^{ c-1}  &\text{for } 0 < x < 1\\
                               \frac{c}{2} x^{-c-1}  &\text{for } x \ge 1
                  \end{cases}

    for ``c > 0``.

    `loglaplace` takes ``c`` as a shape parameter.

    %(after_notes)s

    References
    ----------
    T.J. Kozubowski and K. Podgorski, "A log-Laplace growth rate model",
    The Mathematical Scientist, vol. 28, pp. 49-60, 2003.

    %(example)s

    """
    def _pdf(self, x, c):
        # loglaplace.pdf(x, c) = c / 2 * x**(c-1),   for 0 < x < 1
        #                      = c / 2 * x**(-c-1),  for x >= 1
        cd2 = c/2.0
        c = np.where(x < 1, c, -c)
        return cd2*x**(c-1)

    def _cdf(self, x, c):
        return np.where(x < 1, 0.5*x**c, 1-0.5*x**(-c))

    def _ppf(self, q, c):
        return np.where(q < 0.5, (2.0*q)**(1.0/c), (2*(1.0-q))**(-1.0/c))

    def _munp(self, n, c):
        return c**2 / (c**2 - n**2)

    def _entropy(self, c):
        return np.log(2.0/c) + 1.0


loglaplace = loglaplace_gen(a=0.0, name='loglaplace')


def _lognorm_logpdf(x, s):
    return _lazywhere(x != 0, (x, s),
                      lambda x, s: -np.log(x)**2 / (2*s**2) - np.log(s*x*np.sqrt(2*np.pi)),
                      -np.inf)


class lognorm_gen(rv_continuous):
    r"""A lognormal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `lognorm` is:

    .. math::

        f(x, s) = \frac{1}{s x \sqrt{2\pi}}
                  \exp(-\frac{1}{2} (\frac{\log(x)}{s})^2)

    for ``x > 0``, ``s > 0``.

    `lognorm` takes ``s`` as a shape parameter.

    %(after_notes)s

    A common parametrization for a lognormal random variable ``Y`` is in
    terms of the mean, ``mu``, and standard deviation, ``sigma``, of the
    unique normally distributed random variable ``X`` such that exp(X) = Y.
    This parametrization corresponds to setting ``s = sigma`` and ``scale =
    exp(mu)``.

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _rvs(self, s):
        return np.exp(s * self._random_state.standard_normal(self._size))

    def _pdf(self, x, s):
        # lognorm.pdf(x, s) = 1 / (s*x*sqrt(2*pi)) * exp(-1/2*(log(x)/s)**2)
        return np.exp(self._logpdf(x, s))

    def _logpdf(self, x, s):
        return _lognorm_logpdf(x, s)

    def _cdf(self, x, s):
        return _norm_cdf(np.log(x) / s)

    def _logcdf(self, x, s):
        return _norm_logcdf(np.log(x) / s)

    def _ppf(self, q, s):
        return np.exp(s * _norm_ppf(q))

    def _sf(self, x, s):
        return _norm_sf(np.log(x) / s)

    def _logsf(self, x, s):
        return _norm_logsf(np.log(x) / s)

    def _stats(self, s):
        p = np.exp(s*s)
        mu = np.sqrt(p)
        mu2 = p*(p-1)
        g1 = np.sqrt((p-1))*(2+p)
        g2 = np.polyval([1, 2, 3, 0, -6.0], p)
        return mu, mu2, g1, g2

    def _entropy(self, s):
        return 0.5 * (1 + np.log(2*np.pi) + 2 * np.log(s))

    @extend_notes_in_docstring(rv_continuous, notes="""\
        When the location parameter is fixed by using the `floc` argument,
        this function uses explicit formulas for the maximum likelihood
        estimation of the log-normal shape and scale parameters, so the
        `optimizer`, `loc` and `scale` keyword arguments are ignored.\n\n""")
    def fit(self, data, *args, **kwds):
        floc = kwds.get('floc', None)
        if floc is None:
            # loc is not fixed.  Use the default fit method.
            return super(lognorm_gen, self).fit(data, *args, **kwds)

        f0 = (kwds.get('f0', None) or kwds.get('fs', None) or
              kwds.get('fix_s', None))
        fscale = kwds.get('fscale', None)

        if len(args) > 1:
            raise TypeError("Too many input arguments.")
        for name in ['f0', 'fs', 'fix_s', 'floc', 'fscale', 'loc', 'scale',
                     'optimizer']:
            kwds.pop(name, None)
        if kwds:
            raise TypeError("Unknown arguments: %s." % kwds)

        # Special case: loc is fixed.  Use the maximum likelihood formulas
        # instead of the numerical solver.

        if f0 is not None and fscale is not None:
            # This check is for consistency with `rv_continuous.fit`.
            raise ValueError("All parameters fixed. There is nothing to "
                             "optimize.")

        data = np.asarray(data)
        floc = float(floc)
        if floc != 0:
            # Shifting the data by floc. Don't do the subtraction in-place,
            # because `data` might be a view of the input array.
            data = data - floc
        if np.any(data <= 0):
            raise FitDataError("lognorm", lower=floc, upper=np.inf)
        lndata = np.log(data)

        # Three cases to handle:
        # * shape and scale both free
        # * shape fixed, scale free
        # * shape free, scale fixed

        if fscale is None:
            # scale is free.
            scale = np.exp(lndata.mean())
            if f0 is None:
                # shape is free.
                shape = lndata.std()
            else:
                # shape is fixed.
                shape = float(f0)
        else:
            # scale is fixed, shape is free
            scale = float(fscale)
            shape = np.sqrt(((lndata - np.log(scale))**2).mean())

        return shape, floc, scale


lognorm = lognorm_gen(a=0.0, name='lognorm')


class gilbrat_gen(rv_continuous):
    r"""A Gilbrat continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `gilbrat` is:

    .. math::

        f(x) = \frac{1}{x \sqrt{2\pi}} \exp(-\frac{1}{2} (\log(x))^2)

    `gilbrat` is a special case of `lognorm` with ``s = 1``.

    %(after_notes)s

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _rvs(self):
        return np.exp(self._random_state.standard_normal(self._size))

    def _pdf(self, x):
        # gilbrat.pdf(x) = 1/(x*sqrt(2*pi)) * exp(-1/2*(log(x))**2)
        return np.exp(self._logpdf(x))

    def _logpdf(self, x):
        return _lognorm_logpdf(x, 1.0)

    def _cdf(self, x):
        return _norm_cdf(np.log(x))

    def _ppf(self, q):
        return np.exp(_norm_ppf(q))

    def _stats(self):
        p = np.e
        mu = np.sqrt(p)
        mu2 = p * (p - 1)
        g1 = np.sqrt((p - 1)) * (2 + p)
        g2 = np.polyval([1, 2, 3, 0, -6.0], p)
        return mu, mu2, g1, g2

    def _entropy(self):
        return 0.5 * np.log(2 * np.pi) + 0.5


gilbrat = gilbrat_gen(a=0.0, name='gilbrat')


class maxwell_gen(rv_continuous):
    r"""A Maxwell continuous random variable.

    %(before_notes)s

    Notes
    -----
    A special case of a `chi` distribution,  with ``df = 3``, ``loc = 0.0``,
    and given ``scale = a``, where ``a`` is the parameter used in the
    Mathworld description [1]_.

    The probability density function for `maxwell` is:

    .. math::

        f(x) = \sqrt{2/\pi}x^2 \exp(-x^2/2)

    for ``x > 0``.

    %(after_notes)s

    References
    ----------
    .. [1] http://mathworld.wolfram.com/MaxwellDistribution.html

    %(example)s
    """
    def _rvs(self):
        return chi.rvs(3.0, size=self._size, random_state=self._random_state)

    def _pdf(self, x):
        # maxwell.pdf(x) = sqrt(2/pi)x**2 * exp(-x**2/2)
        return np.sqrt(2.0/np.pi)*x*x*np.exp(-x*x/2.0)

    def _cdf(self, x):
        return sc.gammainc(1.5, x*x/2.0)

    def _ppf(self, q):
        return np.sqrt(2*sc.gammaincinv(1.5, q))

    def _stats(self):
        val = 3*np.pi-8
        return (2*np.sqrt(2.0/np.pi),
                3-8/np.pi,
                np.sqrt(2)*(32-10*np.pi)/val**1.5,
                (-12*np.pi*np.pi + 160*np.pi - 384) / val**2.0)

    def _entropy(self):
        return _EULER + 0.5*np.log(2*np.pi)-0.5


maxwell = maxwell_gen(a=0.0, name='maxwell')


class mielke_gen(rv_continuous):
    r"""A Mielke's Beta-Kappa continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `mielke` is:

    .. math::

        f(x, k, s) = \frac{k x^{k-1}}{(1+x^s)^{1+k/s}}

    for ``x > 0``.

    `mielke` takes ``k`` and ``s`` as shape parameters.

    %(after_notes)s

    %(example)s

    """
    def _pdf(self, x, k, s):
        # mielke.pdf(x, k, s) = k * x**(k-1) / (1+x**s)**(1+k/s)
        return k*x**(k-1.0) / (1.0+x**s)**(1.0+k*1.0/s)

    def _cdf(self, x, k, s):
        return x**k / (1.0+x**s)**(k*1.0/s)

    def _ppf(self, q, k, s):
        qsk = pow(q, s*1.0/k)
        return pow(qsk/(1.0-qsk), 1.0/s)


mielke = mielke_gen(a=0.0, name='mielke')


class kappa4_gen(rv_continuous):
    r"""Kappa 4 parameter distribution.

    %(before_notes)s

    Notes
    -----
    The probability density function for kappa4 is:

    .. math::

        f(x, h, k) = (1 - k x)^{1/k - 1} (1 - h (1 - k x)^{1/k})^{1/h-1}

    if :math:`h` and :math:`k` are not equal to 0.

    If :math:`h` or :math:`k` are zero then the pdf can be simplified:

    h = 0 and k != 0::

        kappa4.pdf(x, h, k) = (1.0 - k*x)**(1.0/k - 1.0)*
                              exp(-(1.0 - k*x)**(1.0/k))

    h != 0 and k = 0::

        kappa4.pdf(x, h, k) = exp(-x)*(1.0 - h*exp(-x))**(1.0/h - 1.0)

    h = 0 and k = 0::

        kappa4.pdf(x, h, k) = exp(-x)*exp(-exp(-x))

    kappa4 takes :math:`h` and :math:`k` as shape parameters.

    The kappa4 distribution returns other distributions when certain
    :math:`h` and :math:`k` values are used.

    +------+-------------+----------------+------------------+
    | h    | k=0.0       | k=1.0          | -inf<=k<=inf     |
    +======+=============+================+==================+
    | -1.0 | Logistic    |                | Generalized      |
    |      |             |                | Logistic(1)      |
    |      |             |                |                  |
    |      | logistic(x) |                |                  |
    +------+-------------+----------------+------------------+
    |  0.0 | Gumbel      | Reverse        | Generalized      |
    |      |             | Exponential(2) | Extreme Value    |
    |      |             |                |                  |
    |      | gumbel_r(x) |                | genextreme(x, k) |
    +------+-------------+----------------+------------------+
    |  1.0 | Exponential | Uniform        | Generalized      |
    |      |             |                | Pareto           |
    |      |             |                |                  |
    |      | expon(x)    | uniform(x)     | genpareto(x, -k) |
    +------+-------------+----------------+------------------+

    (1) There are at least five generalized logistic distributions.
        Four are described here:
        https://en.wikipedia.org/wiki/Generalized_logistic_distribution
        The "fifth" one is the one kappa4 should match which currently
        isn't implemented in scipy:
        https://en.wikipedia.org/wiki/Talk:Generalized_logistic_distribution
        http://www.mathwave.com/help/easyfit/html/analyses/distributions/gen_logistic.html
    (2) This distribution is currently not in scipy.

    References
    ----------
    J.C. Finney, "Optimization of a Skewed Logistic Distribution With Respect
    to the Kolmogorov-Smirnov Test", A Dissertation Submitted to the Graduate
    Faculty of the Louisiana State University and Agricultural and Mechanical
    College, (August, 2004),
    http://digitalcommons.lsu.edu/cgi/viewcontent.cgi?article=4671&context=gradschool_dissertations

    J.R.M. Hosking, "The four-parameter kappa distribution". IBM J. Res.
    Develop. 38 (3), 25 1-258 (1994).

    B. Kumphon, A. Kaew-Man, P. Seenoi, "A Rainfall Distribution for the Lampao
    Site in the Chi River Basin, Thailand", Journal of Water Resource and
    Protection, vol. 4, 866-869, (2012).
    http://file.scirp.org/pdf/JWARP20121000009_14676002.pdf

    C. Winchester, "On Estimation of the Four-Parameter Kappa Distribution", A
    Thesis Submitted to Dalhousie University, Halifax, Nova Scotia, (March
    2000).
    http://www.nlc-bnc.ca/obj/s4/f2/dsk2/ftp01/MQ57336.pdf

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, h, k):
        condlist = [np.logical_and(h > 0, k > 0),
                    np.logical_and(h > 0, k == 0),
                    np.logical_and(h > 0, k < 0),
                    np.logical_and(h <= 0, k > 0),
                    np.logical_and(h <= 0, k == 0),
                    np.logical_and(h <= 0, k < 0)]

        def f0(h, k):
            return (1.0 - float_power(h, -k))/k

        def f1(h, k):
            return np.log(h)

        def f3(h, k):
            a = np.empty(np.shape(h))
            a[:] = -np.inf
            return a

        def f5(h, k):
            return 1.0/k

        self.a = _lazyselect(condlist,
                             [f0, f1, f0, f3, f3, f5],
                             [h, k],
                             default=np.nan)

        def f0(h, k):
            return 1.0/k

        def f1(h, k):
            a = np.empty(np.shape(h))
            a[:] = np.inf
            return a

        self.b = _lazyselect(condlist,
                             [f0, f1, f1, f0, f1, f1],
                             [h, k],
                             default=np.nan)
        return h == h

    def _pdf(self, x, h, k):
        # kappa4.pdf(x, h, k) = (1.0 - k*x)**(1.0/k - 1.0)*
        #                       (1.0 - h*(1.0 - k*x)**(1.0/k))**(1.0/h-1)
        return np.exp(self._logpdf(x, h, k))

    def _logpdf(self, x, h, k):
        condlist = [np.logical_and(h != 0, k != 0),
                    np.logical_and(h == 0, k != 0),
                    np.logical_and(h != 0, k == 0),
                    np.logical_and(h == 0, k == 0)]

        def f0(x, h, k):
            '''pdf = (1.0 - k*x)**(1.0/k - 1.0)*(
                      1.0 - h*(1.0 - k*x)**(1.0/k))**(1.0/h-1.0)
               logpdf = ...
            '''
            return (sc.xlog1py(1.0/k - 1.0, -k*x) +
                    sc.xlog1py(1.0/h - 1.0, -h*(1.0 - k*x)**(1.0/k)))

        def f1(x, h, k):
            '''pdf = (1.0 - k*x)**(1.0/k - 1.0)*np.exp(-(
                      1.0 - k*x)**(1.0/k))
               logpdf = ...
            '''
            return sc.xlog1py(1.0/k - 1.0, -k*x) - (1.0 - k*x)**(1.0/k)

        def f2(x, h, k):
            '''pdf = np.exp(-x)*(1.0 - h*np.exp(-x))**(1.0/h - 1.0)
               logpdf = ...
            '''
            return -x + sc.xlog1py(1.0/h - 1.0, -h*np.exp(-x))

        def f3(x, h, k):
            '''pdf = np.exp(-x-np.exp(-x))
               logpdf = ...
            '''
            return -x - np.exp(-x)

        return _lazyselect(condlist,
                           [f0, f1, f2, f3],
                           [x, h, k],
                           default=np.nan)

    def _cdf(self, x, h, k):
        return np.exp(self._logcdf(x, h, k))

    def _logcdf(self, x, h, k):
        condlist = [np.logical_and(h != 0, k != 0),
                    np.logical_and(h == 0, k != 0),
                    np.logical_and(h != 0, k == 0),
                    np.logical_and(h == 0, k == 0)]

        def f0(x, h, k):
            '''cdf = (1.0 - h*(1.0 - k*x)**(1.0/k))**(1.0/h)
               logcdf = ...
            '''
            return (1.0/h)*sc.log1p(-h*(1.0 - k*x)**(1.0/k))

        def f1(x, h, k):
            '''cdf = np.exp(-(1.0 - k*x)**(1.0/k))
               logcdf = ...
            '''
            return -(1.0 - k*x)**(1.0/k)

        def f2(x, h, k):
            '''cdf = (1.0 - h*np.exp(-x))**(1.0/h)
               logcdf = ...
            '''
            return (1.0/h)*sc.log1p(-h*np.exp(-x))

        def f3(x, h, k):
            '''cdf = np.exp(-np.exp(-x))
               logcdf = ...
            '''
            return -np.exp(-x)

        return _lazyselect(condlist,
                           [f0, f1, f2, f3],
                           [x, h, k],
                           default=np.nan)

    def _ppf(self, q, h, k):
        condlist = [np.logical_and(h != 0, k != 0),
                    np.logical_and(h == 0, k != 0),
                    np.logical_and(h != 0, k == 0),
                    np.logical_and(h == 0, k == 0)]

        def f0(q, h, k):
            return 1.0/k*(1.0 - ((1.0 - (q**h))/h)**k)

        def f1(q, h, k):
            return 1.0/k*(1.0 - (-np.log(q))**k)

        def f2(q, h, k):
            '''ppf = -np.log((1.0 - (q**h))/h)
            '''
            return -sc.log1p(-(q**h)) + np.log(h)

        def f3(q, h, k):
            return -np.log(-np.log(q))

        return _lazyselect(condlist,
                           [f0, f1, f2, f3],
                           [q, h, k],
                           default=np.nan)

    def _stats(self, h, k):
        if h >= 0 and k >= 0:
            maxr = 5
        elif h < 0 and k >= 0:
            maxr = int(-1.0/h*k)
        elif k < 0:
            maxr = int(-1.0/k)
        else:
            maxr = 5

        outputs = [None if r < maxr else np.nan for r in range(1, 5)]
        return outputs[:]


kappa4 = kappa4_gen(name='kappa4')


class kappa3_gen(rv_continuous):
    r"""Kappa 3 parameter distribution.

    %(before_notes)s

    Notes
    -----
    The probability density function for `kappa` is:

    .. math::

        f(x, a) = \begin{cases}
                    a [a + x^a]^{-(a + 1)/a},     &\text{for } x > 0\\
                    0.0,                          &\text{for } x \le 0
                  \end{cases}

    `kappa3` takes :math:`a` as a shape parameter and :math:`a > 0`.

    References
    ----------
    P.W. Mielke and E.S. Johnson, "Three-Parameter Kappa Distribution Maximum
    Likelihood and Likelihood Ratio Tests", Methods in Weather Research,
    701-707, (September, 1973),
    http://docs.lib.noaa.gov/rescue/mwr/101/mwr-101-09-0701.pdf

    B. Kumphon, "Maximum Entropy and Maximum Likelihood Estimation for the
    Three-Parameter Kappa Distribution", Open Journal of Statistics, vol 2,
    415-419 (2012)
    http://file.scirp.org/pdf/OJS20120400011_95789012.pdf

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, a):
        return a > 0

    def _pdf(self, x, a):
        # kappa3.pdf(x, a) =
        #     a*[a + x**a]**(-(a + 1)/a),     for x > 0
        #     0.0,                            for x <= 0
        return a*(a + x**a)**(-1.0/a-1)

    def _cdf(self, x, a):
        return x*(a + x**a)**(-1.0/a)

    def _ppf(self, q, a):
        return (a/(q**-a - 1.0))**(1.0/a)

    def _stats(self, a):
        outputs = [None if i < a else np.nan for i in range(1, 5)]
        return outputs[:]


kappa3 = kappa3_gen(a=0.0, name='kappa3')

class moyal_gen(rv_continuous):
    r"""A Moyal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `moyal` is:

    .. math::

        f(x) = \exp(-(x + \exp(-x))/2) / \sqrt{2\pi}

    %(after_notes)s

    This distribution has utility in high-energy physics and radiation
    detection. It describes the energy loss of a charged relativistic
    particle due to ionization of the medium [1]_. It also provides an
    approximation for the Landau distribution. For an in depth description
    see [2]_. For additional description, see [3]_.

    References
    ----------
    .. [1] J.E. Moyal, "XXX. Theory of ionization fluctuations",
           The London, Edinburgh, and Dublin Philosophical Magazine
           and Journal of Science, vol 46, 263-280, (1955).
           https://doi.org/10.1080/14786440308521076 (gated)
    .. [2] G. Cordeiro et al., "The beta Moyal: a useful skew distribution",
           International Journal of Research and Reviews in Applied Sciences,
           vol 10, 171-192, (2012).
           http://www.arpapress.com/Volumes/Vol10Issue2/IJRRAS_10_2_02.pdf
    .. [3] C. Walck, "Handbook on Statistical Distributions for
           Experimentalists; International Report SUF-PFY/96-01", Chapter 26,
           University of Stockholm: Stockholm, Sweden, (2007).
           www.stat.rice.edu/~dobelman/textfiles/DistributionsHandbook.pdf

    .. versionadded:: 1.1.0

    %(example)s

    """
    def _rvs(self):
        sz, rndm = self._size, self._random_state
        u1 = gamma.rvs(a = 0.5, scale = 2, size=sz, random_state=rndm)
        return -np.log(u1)

    def _pdf(self, x):
        return np.exp(-0.5 * (x + np.exp(-x))) / np.sqrt(2*np.pi)

    def _cdf(self, x):
        return sc.erfc(np.exp(-0.5 * x) / np.sqrt(2))

    def _sf(self, x):
        return sc.erf(np.exp(-0.5 * x) / np.sqrt(2))

    def _ppf(self, x):
        return -np.log(2 * sc.erfcinv(x)**2)

    def _stats(self):
        mu = np.log(2) + np.euler_gamma
        mu2 = np.pi**2 / 2
        g1 = 28 * np.sqrt(2) * sc.zeta(3) / np.pi**3
        g2 = 4.
        return mu, mu2, g1, g2

    def _munp(self, n):
        if n == 1.0:
            return np.log(2) + np.euler_gamma
        elif n == 2.0:
            return np.pi**2 / 2 + (np.log(2) + np.euler_gamma)**2
        elif n == 3.0:
            tmp1 = 1.5 * np.pi**2 * (np.log(2)+np.euler_gamma)
            tmp2 = (np.log(2)+np.euler_gamma)**3
            tmp3 = 14 * sc.zeta(3)
            return tmp1 + tmp2 + tmp3
        elif n == 4.0:
            tmp1 = 4 * 14 * sc.zeta(3) * (np.log(2) + np.euler_gamma)
            tmp2 = 3 * np.pi**2 * (np.log(2) + np.euler_gamma)**2
            tmp3 = (np.log(2) + np.euler_gamma)**4
            tmp4 = 7 * np.pi**4 / 4
            return tmp1 + tmp2 + tmp3 + tmp4
        else:
            # return generic for higher moments
            # return rv_continuous._mom1_sc(self, n, b)
            return self._mom1_sc(n)


moyal = moyal_gen(name="moyal")


class nakagami_gen(rv_continuous):
    r"""A Nakagami continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `nakagami` is:

    .. math::

        f(x, nu) = \frac{2 \nu^\nu}{\Gamma(\nu)} x^{2\nu-1} \exp(-\nu x^2)

    for ``x > 0``, ``nu > 0``.

    `nakagami` takes ``nu`` as a shape parameter.

    %(after_notes)s

    %(example)s

    """
    def _pdf(self, x, nu):
        # nakagami.pdf(x, nu) = 2 * nu**nu / gamma(nu) *
        #                       x**(2*nu-1) * exp(-nu*x**2)
        return 2*nu**nu/sc.gamma(nu)*(x**(2*nu-1.0))*np.exp(-nu*x*x)

    def _cdf(self, x, nu):
        return sc.gammainc(nu, nu*x*x)

    def _ppf(self, q, nu):
        return np.sqrt(1.0/nu*sc.gammaincinv(nu, q))

    def _stats(self, nu):
        mu = sc.gamma(nu+0.5)/sc.gamma(nu)/np.sqrt(nu)
        mu2 = 1.0-mu*mu
        g1 = mu * (1 - 4*nu*mu2) / 2.0 / nu / np.power(mu2, 1.5)
        g2 = -6*mu**4*nu + (8*nu-2)*mu**2-2*nu + 1
        g2 /= nu*mu2**2.0
        return mu, mu2, g1, g2


nakagami = nakagami_gen(a=0.0, name="nakagami")


class ncx2_gen(rv_continuous):
    r"""A non-central chi-squared continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `ncx2` is:

    .. math::

        f(x, df, nc) = \exp(-\frac{nc+x}{2}) \frac{1}{2} (x/nc)^{(df-2)/4}
                       I[(df-2)/2](\sqrt{nc x})

    for :math:`x > 0`.

    `ncx2` takes ``df`` and ``nc`` as shape parameters.

    %(after_notes)s

    %(example)s

    """
    def _rvs(self, df, nc):
        return self._random_state.noncentral_chisquare(df, nc, self._size)

    def _logpdf(self, x, df, nc):
        return _ncx2_log_pdf(x, df, nc)

    def _pdf(self, x, df, nc):
        # ncx2.pdf(x, df, nc) = exp(-(nc+x)/2) * 1/2 * (x/nc)**((df-2)/4)
        #                       * I[(df-2)/2](sqrt(nc*x))
        return _ncx2_pdf(x, df, nc)

    def _cdf(self, x, df, nc):
        return _ncx2_cdf(x, df, nc)

    def _ppf(self, q, df, nc):
        return sc.chndtrix(q, df, nc)

    def _stats(self, df, nc):
        val = df + 2.0*nc
        return (df + nc,
                2*val,
                np.sqrt(8)*(val+nc)/val**1.5,
                12.0*(val+2*nc)/val**2.0)


ncx2 = ncx2_gen(a=0.0, name='ncx2')


class ncf_gen(rv_continuous):
    r"""A non-central F distribution continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `ncf` is:

    .. math::

        f(x, n_1, n_2, \lambda) =
                          \exp(\frac{\lambda}{2} + \lambda n_1 \frac{x}{2(n_1 x+n_2)})
                          n_1^{n_1/2} n_2^{n_2/2} x^{n_1/2 - 1} \\
                          (n_2+n_1 x)^{-(n_1+n_2)/2}
                          \gamma(n_1/2) \gamma(1+n_2/2) \\
                         \frac{L^{\frac{v_1}{2}-1}_{v_2/2}
                               (-\lambda v_1 \frac{x}{2(v_1 x+v_2)})}
                              {(B(v_1/2, v_2/2)  \gamma(\frac{v_1+v_2}{2})}

    for :math:`n_1 > 1`, :math:`n_2, \lambda > 0`.  Here :math:`n_1` is the
    degrees of freedom in the numerator, :math:`n_2` the degrees of freedom in
    the denominator, :math:`\lambda` the non-centrality parameter,
    :math:`\gamma` is the logarithm of the Gamma function, :math:`L_n^k` is a
    generalized Laguerre polynomial and :math:`B` is the beta function.

    `ncf` takes ``df1``, ``df2`` and ``nc`` as shape parameters.

    %(after_notes)s

    %(example)s

    """
    def _rvs(self, dfn, dfd, nc):
        return self._random_state.noncentral_f(dfn, dfd, nc, self._size)

    def _pdf_skip(self, x, dfn, dfd, nc):
        # ncf.pdf(x, df1, df2, nc) = exp(nc/2 + nc*df1*x/(2*(df1*x+df2))) *
        #             df1**(df1/2) * df2**(df2/2) * x**(df1/2-1) *
        #             (df2+df1*x)**(-(df1+df2)/2) *
        #             gamma(df1/2)*gamma(1+df2/2) *
        #             L^{v1/2-1}^{v2/2}(-nc*v1*x/(2*(v1*x+v2))) /
        #             (B(v1/2, v2/2) * gamma((v1+v2)/2))
        n1, n2 = dfn, dfd
        term = -nc/2+nc*n1*x/(2*(n2+n1*x)) + sc.gammaln(n1/2.)+sc.gammaln(1+n2/2.)
        term -= sc.gammaln((n1+n2)/2.0)
        Px = np.exp(term)
        Px *= n1**(n1/2) * n2**(n2/2) * x**(n1/2-1)
        Px *= (n2+n1*x)**(-(n1+n2)/2)
        Px *= sc.assoc_laguerre(-nc*n1*x/(2.0*(n2+n1*x)), n2/2, n1/2-1)
        Px /= sc.beta(n1/2, n2/2)
        # This function does not have a return.  Drop it for now, the generic
        # function seems to work OK.

    def _cdf(self, x, dfn, dfd, nc):
        return sc.ncfdtr(dfn, dfd, nc, x)

    def _ppf(self, q, dfn, dfd, nc):
        return sc.ncfdtri(dfn, dfd, nc, q)

    def _munp(self, n, dfn, dfd, nc):
        val = (dfn * 1.0/dfd)**n
        term = sc.gammaln(n+0.5*dfn) + sc.gammaln(0.5*dfd-n) - sc.gammaln(dfd*0.5)
        val *= np.exp(-nc / 2.0+term)
        val *= sc.hyp1f1(n+0.5*dfn, 0.5*dfn, 0.5*nc)
        return val

    def _stats(self, dfn, dfd, nc):
        mu = np.where(dfd <= 2, np.inf, dfd / (dfd-2.0)*(1+nc*1.0/dfn))
        mu2 = np.where(dfd <= 4, np.inf, 2*(dfd*1.0/dfn)**2.0 *
                       ((dfn+nc/2.0)**2.0 + (dfn+nc)*(dfd-2.0)) /
                       ((dfd-2.0)**2.0 * (dfd-4.0)))
        return mu, mu2, None, None


ncf = ncf_gen(a=0.0, name='ncf')


class t_gen(rv_continuous):
    r"""A Student's T continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `t` is:

    .. math::

        f(x, df) = \frac{\gamma((df+1)/2)}
                        {\sqrt{\pi*df} \gamma(df/2) (1+x^2/df)^{(df+1)/2}}

    for ``df > 0``.

    `t` takes ``df`` as a shape parameter.

    %(after_notes)s

    %(example)s

    """
    def _rvs(self, df):
        return self._random_state.standard_t(df, size=self._size)

    def _pdf(self, x, df):
        #                                gamma((df+1)/2)
        # t.pdf(x, df) = ---------------------------------------------------
        #                sqrt(pi*df) * gamma(df/2) * (1+x**2/df)**((df+1)/2)
        r = np.asarray(df*1.0)
        Px = np.exp(sc.gammaln((r+1)/2)-sc.gammaln(r/2))
        Px /= np.sqrt(r*np.pi)*(1+(x**2)/r)**((r+1)/2)
        return Px

    def _logpdf(self, x, df):
        r = df*1.0
        lPx = sc.gammaln((r+1)/2)-sc.gammaln(r/2)
        lPx -= 0.5*np.log(r*np.pi) + (r+1)/2*np.log(1+(x**2)/r)
        return lPx

    def _cdf(self, x, df):
        return sc.stdtr(df, x)

    def _sf(self, x, df):
        return sc.stdtr(df, -x)

    def _ppf(self, q, df):
        return sc.stdtrit(df, q)

    def _isf(self, q, df):
        return -sc.stdtrit(df, q)

    def _stats(self, df):
        mu2 = _lazywhere(df > 2, (df,),
                         lambda df: df / (df-2.0),
                         np.inf)
        g1 = np.where(df > 3, 0.0, np.nan)
        g2 = _lazywhere(df > 4, (df,),
                        lambda df: 6.0 / (df-4.0),
                        np.nan)
        return 0, mu2, g1, g2


t = t_gen(name='t')


class nct_gen(rv_continuous):
    r"""A non-central Student's T continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `nct` is:

    .. math::

        f(x, df, nc) = \frac{df^{df/2} \gamma(df+1)}{2^{df}
                       \exp(nc^2 / 2) (df+x^2)^{df/2} \gamma(df/2)}

    for ``df > 0``.

    `nct` takes ``df`` and ``nc`` as shape parameters.

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, df, nc):
        return (df > 0) & (nc == nc)

    def _rvs(self, df, nc):
        sz, rndm = self._size, self._random_state
        n = norm.rvs(loc=nc, size=sz, random_state=rndm)
        c2 = chi2.rvs(df, size=sz, random_state=rndm)
        return n * np.sqrt(df) / np.sqrt(c2)

    def _pdf(self, x, df, nc):
        # nct.pdf(x, df, nc) =
        #                               df**(df/2) * gamma(df+1)
        #                ----------------------------------------------------
        #                2**df*exp(nc**2/2) * (df+x**2)**(df/2) * gamma(df/2)
        n = df*1.0
        nc = nc*1.0
        x2 = x*x
        ncx2 = nc*nc*x2
        fac1 = n + x2
        trm1 = n/2.*np.log(n) + sc.gammaln(n+1)
        trm1 -= n*np.log(2)+nc*nc/2.+(n/2.)*np.log(fac1)+sc.gammaln(n/2.)
        Px = np.exp(trm1)
        valF = ncx2 / (2*fac1)
        trm1 = np.sqrt(2)*nc*x*sc.hyp1f1(n/2+1, 1.5, valF)
        trm1 /= np.asarray(fac1*sc.gamma((n+1)/2))
        trm2 = sc.hyp1f1((n+1)/2, 0.5, valF)
        trm2 /= np.asarray(np.sqrt(fac1)*sc.gamma(n/2+1))
        Px *= trm1+trm2
        return Px

    def _cdf(self, x, df, nc):
        return sc.nctdtr(df, nc, x)

    def _ppf(self, q, df, nc):
        return sc.nctdtrit(df, nc, q)

    def _stats(self, df, nc, moments='mv'):
        #
        # See D. Hogben, R.S. Pinkham, and M.B. Wilk,
        # 'The moments of the non-central t-distribution'
        # Biometrika 48, p. 465 (2961).
        # e.g. http://www.jstor.org/stable/2332772 (gated)
        #
        mu, mu2, g1, g2 = None, None, None, None

        gfac = sc.gamma(df/2.-0.5) / sc.gamma(df/2.)
        c11 = np.sqrt(df/2.) * gfac
        c20 = df / (df-2.)
        c22 = c20 - c11*c11
        mu = np.where(df > 1, nc*c11, np.inf)
        mu2 = np.where(df > 2, c22*nc*nc + c20, np.inf)
        if 's' in moments:
            c33t = df * (7.-2.*df) / (df-2.) / (df-3.) + 2.*c11*c11
            c31t = 3.*df / (df-2.) / (df-3.)
            mu3 = (c33t*nc*nc + c31t) * c11*nc
            g1 = np.where(df > 3, mu3 / np.power(mu2, 1.5), np.nan)
        # kurtosis
        if 'k' in moments:
            c44 = df*df / (df-2.) / (df-4.)
            c44 -= c11*c11 * 2.*df*(5.-df) / (df-2.) / (df-3.)
            c44 -= 3.*c11**4
            c42 = df / (df-4.) - c11*c11 * (df-1.) / (df-3.)
            c42 *= 6.*df / (df-2.)
            c40 = 3.*df*df / (df-2.) / (df-4.)

            mu4 = c44 * nc**4 + c42*nc**2 + c40
            g2 = np.where(df > 4, mu4/mu2**2 - 3., np.nan)
        return mu, mu2, g1, g2


nct = nct_gen(name="nct")


class pareto_gen(rv_continuous):
    r"""A Pareto continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `pareto` is:

    .. math::

        f(x, b) = \frac{b}{x^{b+1}}

    for :math:`x \ge 1`, :math:`b > 0`.

    `pareto` takes :math:`b` as a shape parameter.

    %(after_notes)s

    %(example)s

    """
    def _pdf(self, x, b):
        # pareto.pdf(x, b) = b / x**(b+1)
        return b * x**(-b-1)

    def _cdf(self, x, b):
        return 1 - x**(-b)

    def _ppf(self, q, b):
        return pow(1-q, -1.0/b)

    def _sf(self, x, b):
        return x**(-b)

    def _stats(self, b, moments='mv'):
        mu, mu2, g1, g2 = None, None, None, None
        if 'm' in moments:
            mask = b > 1
            bt = np.extract(mask, b)
            mu = valarray(np.shape(b), value=np.inf)
            np.place(mu, mask, bt / (bt-1.0))
        if 'v' in moments:
            mask = b > 2
            bt = np.extract(mask, b)
            mu2 = valarray(np.shape(b), value=np.inf)
            np.place(mu2, mask, bt / (bt-2.0) / (bt-1.0)**2)
        if 's' in moments:
            mask = b > 3
            bt = np.extract(mask, b)
            g1 = valarray(np.shape(b), value=np.nan)
            vals = 2 * (bt + 1.0) * np.sqrt(bt - 2.0) / ((bt - 3.0) * np.sqrt(bt))
            np.place(g1, mask, vals)
        if 'k' in moments:
            mask = b > 4
            bt = np.extract(mask, b)
            g2 = valarray(np.shape(b), value=np.nan)
            vals = (6.0*np.polyval([1.0, 1.0, -6, -2], bt) /
                    np.polyval([1.0, -7.0, 12.0, 0.0], bt))
            np.place(g2, mask, vals)
        return mu, mu2, g1, g2

    def _entropy(self, c):
        return 1 + 1.0/c - np.log(c)


pareto = pareto_gen(a=1.0, name="pareto")


class lomax_gen(rv_continuous):
    r"""A Lomax (Pareto of the second kind) continuous random variable.

    %(before_notes)s

    Notes
    -----
    The Lomax distribution is a special case of the Pareto distribution, with
    (loc=-1.0).

    The probability density function for `lomax` is:

    .. math::

        f(x, c) = \frac{c}{(1+x)^{c+1}}

    for :math:`x \ge 0`, ``c > 0``.

    `lomax` takes :math:`c` as a shape parameter.

    %(after_notes)s

    %(example)s

    """
    def _pdf(self, x, c):
        # lomax.pdf(x, c) = c / (1+x)**(c+1)
        return c*1.0/(1.0+x)**(c+1.0)

    def _logpdf(self, x, c):
        return np.log(c) - (c+1)*sc.log1p(x)

    def _cdf(self, x, c):
        return -sc.expm1(-c*sc.log1p(x))

    def _sf(self, x, c):
        return np.exp(-c*sc.log1p(x))

    def _logsf(self, x, c):
        return -c*sc.log1p(x)

    def _ppf(self, q, c):
        return sc.expm1(-sc.log1p(-q)/c)

    def _stats(self, c):
        mu, mu2, g1, g2 = pareto.stats(c, loc=-1.0, moments='mvsk')
        return mu, mu2, g1, g2

    def _entropy(self, c):
        return 1+1.0/c-np.log(c)


lomax = lomax_gen(a=0.0, name="lomax")


class pearson3_gen(rv_continuous):
    r"""A pearson type III continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `pearson3` is:

    .. math::

        f(x, skew) = \frac{|\beta|}{\gamma(\alpha)}
                     (\beta (x - \zeta))^{alpha - 1} \exp(-\beta (x - \zeta))

    where:

    .. math::

            \beta = \frac{2}{skew  stddev}
            \alpha = (stddev \beta)^2
            \zeta = loc - \frac{\alpha}{\beta}

    `pearson3` takes ``skew`` as a shape parameter.

    %(after_notes)s

    %(example)s

    References
    ----------
    R.W. Vogel and D.E. McMartin, "Probability Plot Goodness-of-Fit and
    Skewness Estimation Procedures for the Pearson Type 3 Distribution", Water
    Resources Research, Vol.27, 3149-3158 (1991).

    L.R. Salvosa, "Tables of Pearson's Type III Function", Ann. Math. Statist.,
    Vol.1, 191-198 (1930).

    "Using Modern Computing Tools to Fit the Pearson Type III Distribution to
    Aviation Loads Data", Office of Aviation Research (2003).

    """
    def _preprocess(self, x, skew):
        # The real 'loc' and 'scale' are handled in the calling pdf(...). The
        # local variables 'loc' and 'scale' within pearson3._pdf are set to
        # the defaults just to keep them as part of the equations for
        # documentation.
        loc = 0.0
        scale = 1.0

        # If skew is small, return _norm_pdf. The divide between pearson3
        # and norm was found by brute force and is approximately a skew of
        # 0.000016.  No one, I hope, would actually use a skew value even
        # close to this small.
        norm2pearson_transition = 0.000016

        ans, x, skew = np.broadcast_arrays([1.0], x, skew)
        ans = ans.copy()

        # mask is True where skew is small enough to use the normal approx.
        mask = np.absolute(skew) < norm2pearson_transition
        invmask = ~mask

        beta = 2.0 / (skew[invmask] * scale)
        alpha = (scale * beta)**2
        zeta = loc - alpha / beta

        transx = beta * (x[invmask] - zeta)
        return ans, x, transx, mask, invmask, beta, alpha, zeta

    def _argcheck(self, skew):
        # The _argcheck function in rv_continuous only allows positive
        # arguments.  The skew argument for pearson3 can be zero (which I want
        # to handle inside pearson3._pdf) or negative.  So just return True
        # for all skew args.
        return np.ones(np.shape(skew), dtype=bool)

    def _stats(self, skew):
        _, _, _, _, _, beta, alpha, zeta = (
            self._preprocess([1], skew))
        m = zeta + alpha / beta
        v = alpha / (beta**2)
        s = 2.0 / (alpha**0.5) * np.sign(beta)
        k = 6.0 / alpha
        return m, v, s, k

    def _pdf(self, x, skew):
        # pearson3.pdf(x, skew) = abs(beta) / gamma(alpha) *
        #     (beta * (x - zeta))**(alpha - 1) * exp(-beta*(x - zeta))
        # Do the calculation in _logpdf since helps to limit
        # overflow/underflow problems
        ans = np.exp(self._logpdf(x, skew))
        if ans.ndim == 0:
            if np.isnan(ans):
                return 0.0
            return ans
        ans[np.isnan(ans)] = 0.0
        return ans

    def _logpdf(self, x, skew):
        #   PEARSON3 logpdf                           GAMMA logpdf
        #   np.log(abs(beta))
        # + (alpha - 1)*np.log(beta*(x - zeta))          + (a - 1)*np.log(x)
        # - beta*(x - zeta)                           - x
        # - sc.gammalnalpha)                              - sc.gammalna)
        ans, x, transx, mask, invmask, beta, alpha, _ = (
            self._preprocess(x, skew))

        ans[mask] = np.log(_norm_pdf(x[mask]))
        ans[invmask] = np.log(abs(beta)) + gamma._logpdf(transx, alpha)
        return ans

    def _cdf(self, x, skew):
        ans, x, transx, mask, invmask, _, alpha, _ = (
            self._preprocess(x, skew))

        ans[mask] = _norm_cdf(x[mask])
        ans[invmask] = gamma._cdf(transx, alpha)
        return ans

    def _rvs(self, skew):
        skew = broadcast_to(skew, self._size)
        ans, _, _, mask, invmask, beta, alpha, zeta = (
            self._preprocess([0], skew))

        nsmall = mask.sum()
        nbig = mask.size - nsmall
        ans[mask] = self._random_state.standard_normal(nsmall)
        ans[invmask] = (self._random_state.standard_gamma(alpha, nbig)/beta +
                        zeta)

        if self._size == ():
            ans = ans[0]
        return ans

    def _ppf(self, q, skew):
        ans, q, _, mask, invmask, beta, alpha, zeta = (
            self._preprocess(q, skew))
        ans[mask] = _norm_ppf(q[mask])
        ans[invmask] = sc.gammaincinv(alpha, q[invmask])/beta + zeta
        return ans


pearson3 = pearson3_gen(name="pearson3")


class powerlaw_gen(rv_continuous):
    r"""A power-function continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `powerlaw` is:

    .. math::

        f(x, a) = a x^{a-1}

    for :math:`0 \le x \le 1`, :math:`a > 0`.

    `powerlaw` takes :math:`a` as a shape parameter.

    %(after_notes)s

    `powerlaw` is a special case of `beta` with ``b == 1``.

    %(example)s

    """
    def _pdf(self, x, a):
        # powerlaw.pdf(x, a) = a * x**(a-1)
        return a*x**(a-1.0)

    def _logpdf(self, x, a):
        return np.log(a) + sc.xlogy(a - 1, x)

    def _cdf(self, x, a):
        return x**(a*1.0)

    def _logcdf(self, x, a):
        return a*np.log(x)

    def _ppf(self, q, a):
        return pow(q, 1.0/a)

    def _stats(self, a):
        return (a / (a + 1.0),
                a / (a + 2.0) / (a + 1.0) ** 2,
                -2.0 * ((a - 1.0) / (a + 3.0)) * np.sqrt((a + 2.0) / a),
                6 * np.polyval([1, -1, -6, 2], a) / (a * (a + 3.0) * (a + 4)))

    def _entropy(self, a):
        return 1 - 1.0/a - np.log(a)


powerlaw = powerlaw_gen(a=0.0, b=1.0, name="powerlaw")


class powerlognorm_gen(rv_continuous):
    r"""A power log-normal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `powerlognorm` is:

    .. math::

        f(x, c, s) = \frac{c}{x s} \phi(\log(x)/s)
                     (\Phi(-\log(x)/s))^{c-1}

    where :math:`\phi` is the normal pdf, and :math:`\Phi` is the normal cdf,
    and :math:`x > 0`, :math:`s, c > 0`.

    `powerlognorm` takes :math:`c` and :math:`s` as shape parameters.

    %(after_notes)s

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _pdf(self, x, c, s):
        # powerlognorm.pdf(x, c, s) = c / (x*s) * phi(log(x)/s) *
        #                                         (Phi(-log(x)/s))**(c-1),
        return (c/(x*s) * _norm_pdf(np.log(x)/s) *
                pow(_norm_cdf(-np.log(x)/s), c*1.0-1.0))

    def _cdf(self, x, c, s):
        return 1.0 - pow(_norm_cdf(-np.log(x)/s), c*1.0)

    def _ppf(self, q, c, s):
        return np.exp(-s * _norm_ppf(pow(1.0 - q, 1.0 / c)))


powerlognorm = powerlognorm_gen(a=0.0, name="powerlognorm")


class powernorm_gen(rv_continuous):
    r"""A power normal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `powernorm` is:

    .. math::

        f(x, c) = c \phi(x) (\Phi(-x))^{c-1}

    where :math:`\phi` is the normal pdf, and :math:`\Phi` is the normal cdf,
    and :math:`x > 0`, :math:`c > 0`.

    `powernorm` takes :math:`c` as a shape parameter.

    %(after_notes)s

    %(example)s

    """
    def _pdf(self, x, c):
        # powernorm.pdf(x, c) = c * phi(x) * (Phi(-x))**(c-1)
        return c*_norm_pdf(x) * (_norm_cdf(-x)**(c-1.0))

    def _logpdf(self, x, c):
        return np.log(c) + _norm_logpdf(x) + (c-1)*_norm_logcdf(-x)

    def _cdf(self, x, c):
        return 1.0-_norm_cdf(-x)**(c*1.0)

    def _ppf(self, q, c):
        return -_norm_ppf(pow(1.0 - q, 1.0 / c))


powernorm = powernorm_gen(name='powernorm')


class rdist_gen(rv_continuous):
    r"""An R-distributed continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `rdist` is:

    .. math::

        f(x, c) = \frac{(1-x^2)^{c/2-1}}{B(1/2, c/2)}

    for :math:`-1 \le x \le 1`, :math:`c > 0`.

    `rdist` takes :math:`c` as a shape parameter.

    This distribution includes the following distribution kernels as
    special cases::

        c = 2:  uniform
        c = 4:  Epanechnikov (parabolic)
        c = 6:  quartic (biweight)
        c = 8:  triweight

    %(after_notes)s

    %(example)s

    """
    def _pdf(self, x, c):
        # rdist.pdf(x, c) = (1-x**2)**(c/2-1) / B(1/2, c/2)
        return np.power((1.0 - x**2), c / 2.0 - 1) / sc.beta(0.5, c / 2.0)

    def _cdf(self, x, c):
        term1 = x / sc.beta(0.5, c / 2.0)
        res = 0.5 + term1 * sc.hyp2f1(0.5, 1 - c / 2.0, 1.5, x**2)
        # There's an issue with hyp2f1, it returns nans near x = +-1, c > 100.
        # Use the generic implementation in that case.  See gh-1285 for
        # background.
        if np.any(np.isnan(res)):
            return rv_continuous._cdf(self, x, c)
        return res

    def _munp(self, n, c):
        numerator = (1 - (n % 2)) * sc.beta((n + 1.0) / 2, c / 2.0)
        return numerator / sc.beta(1. / 2, c / 2.)


rdist = rdist_gen(a=-1.0, b=1.0, name="rdist")


class rayleigh_gen(rv_continuous):
    r"""A Rayleigh continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `rayleigh` is:

    .. math::

        f(r) = r \exp(-r^2/2)

    for :math:`x \ge 0`.

    `rayleigh` is a special case of `chi` with ``df == 2``.

    %(after_notes)s

    %(example)s

    """
    _support_mask = rv_continuous._open_support_mask

    def _rvs(self):
        return chi.rvs(2, size=self._size, random_state=self._random_state)

    def _pdf(self, r):
        # rayleigh.pdf(r) = r * exp(-r**2/2)
        return np.exp(self._logpdf(r))

    def _logpdf(self, r):
        return np.log(r) - 0.5 * r * r

    def _cdf(self, r):
        return -sc.expm1(-0.5 * r**2)

    def _ppf(self, q):
        return np.sqrt(-2 * sc.log1p(-q))

    def _sf(self, r):
        return np.exp(self._logsf(r))

    def _logsf(self, r):
        return -0.5 * r * r

    def _isf(self, q):
        return np.sqrt(-2 * np.log(q))

    def _stats(self):
        val = 4 - np.pi
        return (np.sqrt(np.pi/2),
                val/2,
                2*(np.pi-3)*np.sqrt(np.pi)/val**1.5,
                6*np.pi/val-16/val**2)

    def _entropy(self):
        return _EULER/2.0 + 1 - 0.5*np.log(2)


rayleigh = rayleigh_gen(a=0.0, name="rayleigh")


class reciprocal_gen(rv_continuous):
    r"""A reciprocal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `reciprocal` is:

    .. math::

        f(x, a, b) = \frac{1}{x \log(b/a)}

    for :math:`a \le x \le b`, :math:`a, b > 0`.

    `reciprocal` takes :math:`a` and :math:`b` as shape parameters.

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, a, b):
        self.a = a
        self.b = b
        self.d = np.log(b*1.0 / a)
        return (a > 0) & (b > 0) & (b > a)

    def _pdf(self, x, a, b):
        # reciprocal.pdf(x, a, b) = 1 / (x*log(b/a))
        return 1.0 / (x * self.d)

    def _logpdf(self, x, a, b):
        return -np.log(x) - np.log(self.d)

    def _cdf(self, x, a, b):
        return (np.log(x)-np.log(a)) / self.d

    def _ppf(self, q, a, b):
        return a*pow(b*1.0/a, q)

    def _munp(self, n, a, b):
        return 1.0/self.d / n * (pow(b*1.0, n) - pow(a*1.0, n))

    def _entropy(self, a, b):
        return 0.5*np.log(a*b)+np.log(np.log(b/a))


reciprocal = reciprocal_gen(name="reciprocal")


class rice_gen(rv_continuous):
    r"""A Rice continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `rice` is:

    .. math::

        f(x, b) = x \exp(- \frac{x^2 + b^2}{2}) I[0](x b)

    for :math:`x > 0`, :math:`b > 0`.

    `rice` takes :math:`b` as a shape parameter.

    %(after_notes)s

    The Rice distribution describes the length, :math:`r`, of a 2-D vector with
    components :math:`(U+u, V+v)`, where :math:`U, V` are constant, :math:`u,
    v` are independent Gaussian random variables with standard deviation
    :math:`s`.  Let :math:`R = \sqrt{U^2 + V^2}`. Then the pdf of :math:`r` is
    ``rice.pdf(x, R/s, scale=s)``.

    %(example)s

    """
    def _argcheck(self, b):
        return b >= 0

    def _rvs(self, b):
        # http://en.wikipedia.org/wiki/Rice_distribution
        t = b/np.sqrt(2) + self._random_state.standard_normal(size=(2,) +
                                                              self._size)
        return np.sqrt((t*t).sum(axis=0))

    def _cdf(self, x, b):
        return sc.chndtr(np.square(x), 2, np.square(b))

    def _ppf(self, q, b):
        return np.sqrt(sc.chndtrix(q, 2, np.square(b)))

    def _pdf(self, x, b):
        # rice.pdf(x, b) = x * exp(-(x**2+b**2)/2) * I[0](x*b)
        #
        # We use (x**2 + b**2)/2 = ((x-b)**2)/2 + xb.
        # The factor of np.exp(-xb) is then included in the i0e function
        # in place of the modified Bessel function, i0, improving
        # numerical stability for large values of xb.
        return x * np.exp(-(x-b)*(x-b)/2.0) * sc.i0e(x*b)

    def _munp(self, n, b):
        nd2 = n/2.0
        n1 = 1 + nd2
        b2 = b*b/2.0
        return (2.0**(nd2) * np.exp(-b2) * sc.gamma(n1) *
                sc.hyp1f1(n1, 1, b2))


rice = rice_gen(a=0.0, name="rice")


# FIXME: PPF does not work.
class recipinvgauss_gen(rv_continuous):
    r"""A reciprocal inverse Gaussian continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `recipinvgauss` is:

    .. math::

        f(x, \mu) = \frac{1}{\sqrt{2\pi x}} \frac{\exp(-(1-\mu x)^2}{2x\mu^2)}

    for :math:`x \ge 0`.

    `recipinvgauss` takes :math:`\mu` as a shape parameter.

    %(after_notes)s

    %(example)s

    """

    def _pdf(self, x, mu):
        # recipinvgauss.pdf(x, mu) =
        #                     1/sqrt(2*pi*x) * exp(-(1-mu*x)**2/(2*x*mu**2))
        return 1.0/np.sqrt(2*np.pi*x)*np.exp(-(1-mu*x)**2.0 / (2*x*mu**2.0))

    def _logpdf(self, x, mu):
        return -(1-mu*x)**2.0 / (2*x*mu**2.0) - 0.5*np.log(2*np.pi*x)

    def _cdf(self, x, mu):
        trm1 = 1.0/mu - x
        trm2 = 1.0/mu + x
        isqx = 1.0/np.sqrt(x)
        return 1.0-_norm_cdf(isqx*trm1)-np.exp(2.0/mu)*_norm_cdf(-isqx*trm2)

    def _rvs(self, mu):
        return 1.0/self._random_state.wald(mu, 1.0, size=self._size)


recipinvgauss = recipinvgauss_gen(a=0.0, name='recipinvgauss')


class semicircular_gen(rv_continuous):
    r"""A semicircular continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `semicircular` is:

    .. math::

        f(x) = \frac{2}{\pi} \sqrt{1-x^2}

    for :math:`-1 \le x \le 1`.

    %(after_notes)s

    %(example)s

    """
    def _pdf(self, x):
        # semicircular.pdf(x) = 2/pi * sqrt(1-x**2)
        return 2.0/np.pi*np.sqrt(1-x*x)

    def _cdf(self, x):
        return 0.5+1.0/np.pi*(x*np.sqrt(1-x*x) + np.arcsin(x))

    def _stats(self):
        return 0, 0.25, 0, -1.0

    def _entropy(self):
        return 0.64472988584940017414


semicircular = semicircular_gen(a=-1.0, b=1.0, name="semicircular")


class skew_norm_gen(rv_continuous):
    r"""A skew-normal random variable.

    %(before_notes)s

    Notes
    -----
    The pdf is::

        skewnorm.pdf(x, a) = 2 * norm.pdf(x) * norm.cdf(a*x)

    `skewnorm` takes :math:`a` as a skewness parameter
    When ``a = 0`` the distribution is identical to a normal distribution.
    rvs implements the method of [1]_.

    %(after_notes)s

    %(example)s

    References
    ----------
    .. [1] A. Azzalini and A. Capitanio (1999). Statistical applications of the
        multivariate skew-normal distribution. J. Roy. Statist. Soc., B 61, 579-602.
        http://azzalini.stat.unipd.it/SN/faq-r.html

    """
    def _argcheck(self, a):
        return np.isfinite(a)

    def _pdf(self, x, a):
        return 2.*_norm_pdf(x)*_norm_cdf(a*x)

    def _cdf_single(self, x, *args):
        if x <= 0:
            cdf = integrate.quad(self._pdf, self.a, x, args=args)[0]
        else:
            t1 = integrate.quad(self._pdf, self.a, 0, args=args)[0]
            t2 = integrate.quad(self._pdf, 0, x, args=args)[0]
            cdf = t1 + t2
        if cdf > 1:
            # Presumably numerical noise, e.g. 1.0000000000000002
            cdf = 1.0
        return cdf

    def _sf(self, x, a):
        return self._cdf(-x, -a)

    def _rvs(self, a):
        u0 = self._random_state.normal(size=self._size)
        v = self._random_state.normal(size=self._size)
        d = a/np.sqrt(1 + a**2)
        u1 = d*u0 + v*np.sqrt(1 - d**2)
        return np.where(u0 >= 0, u1, -u1)

    def _stats(self, a, moments='mvsk'):
        output = [None, None, None, None]
        const = np.sqrt(2/np.pi) * a/np.sqrt(1 + a**2)

        if 'm' in moments:
            output[0] = const
        if 'v' in moments:
            output[1] = 1 - const**2
        if 's' in moments:
            output[2] = ((4 - np.pi)/2) * (const/np.sqrt(1 - const**2))**3
        if 'k' in moments:
            output[3] = (2*(np.pi - 3)) * (const**4/(1 - const**2)**2)

        return output


skewnorm = skew_norm_gen(name='skewnorm')


class trapz_gen(rv_continuous):
    r"""A trapezoidal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The trapezoidal distribution can be represented with an up-sloping line
    from ``loc`` to ``(loc + c*scale)``, then constant to ``(loc + d*scale)``
    and then downsloping from ``(loc + d*scale)`` to ``(loc+scale)``.

    `trapz` takes :math:`c` and :math:`d` as shape parameters.

    %(after_notes)s

    The standard form is in the range [0, 1] with c the mode.
    The location parameter shifts the start to `loc`.
    The scale parameter changes the width from 1 to `scale`.

    %(example)s

    """
    def _argcheck(self, c, d):
        return (c >= 0) & (c <= 1) & (d >= 0) & (d <= 1) & (d >= c)

    def _pdf(self, x, c, d):
        u = 2 / (d-c+1)

        return _lazyselect([x < c,
                            (c <= x) & (x <= d),
                            x > d],
                           [lambda x, c, d, u: u * x / c,
                            lambda x, c, d, u: u,
                            lambda x, c, d, u: u * (1-x) / (1-d)],
                            (x, c, d, u))

    def _cdf(self, x, c, d):
        return _lazyselect([x < c,
                            (c <= x) & (x <= d),
                            x > d],
                           [lambda x, c, d: x**2 / c / (d-c+1),
                            lambda x, c, d: (c + 2 * (x-c)) / (d-c+1),
                            lambda x, c, d: 1-((1-x) ** 2
                                               / (d-c+1) / (1-d))],
                            (x, c, d))

    def _ppf(self, q, c, d):
        qc, qd = self._cdf(c, c, d), self._cdf(d, c, d)
        condlist = [q < qc, q <= qd, q > qd]
        choicelist = [np.sqrt(q * c * (1 + d - c)),
                      0.5 * q * (1 + d - c) + 0.5 * c,
                      1 - np.sqrt((1 - q) * (d - c + 1) * (1 - d))]
        return np.select(condlist, choicelist)


trapz = trapz_gen(a=0.0, b=1.0, name="trapz")


class triang_gen(rv_continuous):
    r"""A triangular continuous random variable.

    %(before_notes)s

    Notes
    -----
    The triangular distribution can be represented with an up-sloping line from
    ``loc`` to ``(loc + c*scale)`` and then downsloping for ``(loc + c*scale)``
    to ``(loc+scale)``.

    `triang` takes :math:`c` as a shape parameter.

    %(after_notes)s

    The standard form is in the range [0, 1] with c the mode.
    The location parameter shifts the start to `loc`.
    The scale parameter changes the width from 1 to `scale`.

    %(example)s

    """
    def _rvs(self, c):
        return self._random_state.triangular(0, c, 1, self._size)

    def _argcheck(self, c):
        return (c >= 0) & (c <= 1)

    def _pdf(self, x, c):
        # 0: edge case where c=0
        # 1: generalised case for x < c, don't use x <= c, as it doesn't cope
        #    with c = 0.
        # 2: generalised case for x >= c, but doesn't cope with c = 1
        # 3: edge case where c=1
        r = _lazyselect([c == 0,
                         x < c,
                         (x >= c) & (c != 1),
                         c == 1],
                        [lambda x, c: 2 - 2 * x,
                         lambda x, c: 2 * x / c,
                         lambda x, c: 2 * (1 - x) / (1 - c),
                         lambda x, c: 2 * x],
                        (x, c))
        return r

    def _cdf(self, x, c):
        r = _lazyselect([c == 0,
                         x < c,
                         (x >= c) & (c != 1),
                         c == 1],
                        [lambda x, c: 2*x - x*x,
                         lambda x, c: x * x / c,
                         lambda x, c: (x*x - 2*x + c) / (c-1),
                         lambda x, c: x * x],
                        (x, c))
        return r

    def _ppf(self, q, c):
        return np.where(q < c, np.sqrt(c * q), 1-np.sqrt((1-c) * (1-q)))

    def _stats(self, c):
        return ((c+1.0)/3.0,
                (1.0-c+c*c)/18,
                np.sqrt(2)*(2*c-1)*(c+1)*(c-2) / (5*np.power((1.0-c+c*c), 1.5)),
                -3.0/5.0)

    def _entropy(self, c):
        return 0.5-np.log(2)


triang = triang_gen(a=0.0, b=1.0, name="triang")


class truncexpon_gen(rv_continuous):
    r"""A truncated exponential continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `truncexpon` is:

    .. math::

        f(x, b) = \frac{\exp(-x)}{1 - \exp(-b)}

    for :math:`0 < x < b`.

    `truncexpon` takes :math:`b` as a shape parameter.

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, b):
        self.b = b
        return b > 0

    def _pdf(self, x, b):
        # truncexpon.pdf(x, b) = exp(-x) / (1-exp(-b))
        return np.exp(-x)/(-sc.expm1(-b))

    def _logpdf(self, x, b):
        return -x - np.log(-sc.expm1(-b))

    def _cdf(self, x, b):
        return sc.expm1(-x)/sc.expm1(-b)

    def _ppf(self, q, b):
        return -sc.log1p(q*sc.expm1(-b))

    def _munp(self, n, b):
        # wrong answer with formula, same as in continuous.pdf
        # return sc.gamman+1)-sc.gammainc1+n, b)
        if n == 1:
            return (1-(b+1)*np.exp(-b))/(-sc.expm1(-b))
        elif n == 2:
            return 2*(1-0.5*(b*b+2*b+2)*np.exp(-b))/(-sc.expm1(-b))
        else:
            # return generic for higher moments
            # return rv_continuous._mom1_sc(self, n, b)
            return self._mom1_sc(n, b)

    def _entropy(self, b):
        eB = np.exp(b)
        return np.log(eB-1)+(1+eB*(b-1.0))/(1.0-eB)


truncexpon = truncexpon_gen(a=0.0, name='truncexpon')


class truncnorm_gen(rv_continuous):
    r"""A truncated normal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The standard form of this distribution is a standard normal truncated to
    the range [a, b] --- notice that a and b are defined over the domain of the
    standard normal.  To convert clip values for a specific mean and standard
    deviation, use::

        a, b = (myclip_a - my_mean) / my_std, (myclip_b - my_mean) / my_std

    `truncnorm` takes :math:`a` and :math:`b` as shape parameters.

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, a, b):
        self.a = a
        self.b = b
        self._nb = _norm_cdf(b)
        self._na = _norm_cdf(a)
        self._sb = _norm_sf(b)
        self._sa = _norm_sf(a)
        self._delta = np.where(self.a > 0,
                               -(self._sb - self._sa),
                               self._nb - self._na)
        self._logdelta = np.log(self._delta)
        return a != b

    def _pdf(self, x, a, b):
        return _norm_pdf(x) / self._delta

    def _logpdf(self, x, a, b):
        return _norm_logpdf(x) - self._logdelta

    def _cdf(self, x, a, b):
        return (_norm_cdf(x) - self._na) / self._delta

    def _ppf(self, q, a, b):
        # XXX Use _lazywhere...
        ppf = np.where(self.a > 0,
                       _norm_isf(q*self._sb + self._sa*(1.0-q)),
                       _norm_ppf(q*self._nb + self._na*(1.0-q)))
        return ppf

    def _stats(self, a, b):
        nA, nB = self._na, self._nb
        d = nB - nA
        pA, pB = _norm_pdf(a), _norm_pdf(b)
        mu = (pA - pB) / d   # correction sign
        mu2 = 1 + (a*pA - b*pB) / d - mu*mu
        return mu, mu2, None, None


truncnorm = truncnorm_gen(name='truncnorm')


# FIXME: RVS does not work.
class tukeylambda_gen(rv_continuous):
    r"""A Tukey-Lamdba continuous random variable.

    %(before_notes)s

    Notes
    -----
    A flexible distribution, able to represent and interpolate between the
    following distributions:

    - Cauchy                (lam=-1)
    - logistic              (lam=0.0)
    - approx Normal         (lam=0.14)
    - u-shape               (lam = 0.5)
    - uniform from -1 to 1  (lam = 1)

    `tukeylambda` takes ``lam`` as a shape parameter.

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, lam):
        return np.ones(np.shape(lam), dtype=bool)

    def _pdf(self, x, lam):
        Fx = np.asarray(sc.tklmbda(x, lam))
        Px = Fx**(lam-1.0) + (np.asarray(1-Fx))**(lam-1.0)
        Px = 1.0/np.asarray(Px)
        return np.where((lam <= 0) | (abs(x) < 1.0/np.asarray(lam)), Px, 0.0)

    def _cdf(self, x, lam):
        return sc.tklmbda(x, lam)

    def _ppf(self, q, lam):
        return sc.boxcox(q, lam) - sc.boxcox1p(-q, lam)

    def _stats(self, lam):
        return 0, _tlvar(lam), 0, _tlkurt(lam)

    def _entropy(self, lam):
        def integ(p):
            return np.log(pow(p, lam-1)+pow(1-p, lam-1))
        return integrate.quad(integ, 0, 1)[0]


tukeylambda = tukeylambda_gen(name='tukeylambda')


class FitUniformFixedScaleDataError(FitDataError):
    def __init__(self, ptp, fscale):
        self.args = (
            "Invalid values in `data`.  Maximum likelihood estimation with "
            "the uniform distribution and fixed scale requires that "
            "data.ptp() <= fscale, but data.ptp() = %r and fscale = %r." %
            (ptp, fscale),
        )


class uniform_gen(rv_continuous):
    r"""A uniform continuous random variable.

    This distribution is constant between `loc` and ``loc + scale``.

    %(before_notes)s

    %(example)s

    """
    def _rvs(self):
        return self._random_state.uniform(0.0, 1.0, self._size)

    def _pdf(self, x):
        return 1.0*(x == x)

    def _cdf(self, x):
        return x

    def _ppf(self, q):
        return q

    def _stats(self):
        return 0.5, 1.0/12, 0, -1.2

    def _entropy(self):
        return 0.0

    def fit(self, data, *args, **kwds):
        """
        Maximum likelihood estimate for the location and scale parameters.

        `uniform.fit` uses only the following parameters.  Because exact
        formulas are used, the parameters related to optimization that are
        available in the `fit` method of other distributions are ignored
        here.  The only positional argument accepted is `data`.

        Parameters
        ----------
        data : array_like
            Data to use in calculating the maximum likelihood estimate.
        floc : float, optional
            Hold the location parameter fixed to the specified value.
        fscale : float, optional
            Hold the scale parameter fixed to the specified value.

        Returns
        -------
        loc, scale : float
            Maximum likelihood estimates for the location and scale.

        Notes
        -----
        An error is raised if `floc` is given and any values in `data` are
        less than `floc`, or if `fscale` is given and `fscale` is less
        than ``data.max() - data.min()``.  An error is also raised if both
        `floc` and `fscale` are given.

        Examples
        --------
        >>> from scipy.stats import uniform

        We'll fit the uniform distribution to `x`:

        >>> x = np.array([2, 2.5, 3.1, 9.5, 13.0])

        For a uniform distribution MLE, the location is the minimum of the
        data, and the scale is the maximum minus the minimum.

        >>> loc, scale = uniform.fit(x)
        >>> loc
        2.0
        >>> scale
        11.0

        If we know the data comes from a uniform distribution where the support
        starts at 0, we can use `floc=0`:

        >>> loc, scale = uniform.fit(x, floc=0)
        >>> loc
        0.0
        >>> scale
        13.0

        Alternatively, if we know the length of the support is 12, we can use
        `fscale=12`:

        >>> loc, scale = uniform.fit(x, fscale=12)
        >>> loc
        1.5
        >>> scale
        12.0

        In that last example, the support interval is [1.5, 13.5].  This
        solution is not unique.  For example, the distribution with ``loc=2``
        and ``scale=12`` has the same likelihood as the one above.  When
        `fscale` is given and it is larger than ``data.max() - data.min()``,
        the parameters returned by the `fit` method center the support over
        the interval ``[data.min(), data.max()]``.

        """
        if len(args) > 0:
            raise TypeError("Too many arguments.")

        floc = kwds.pop('floc', None)
        fscale = kwds.pop('fscale', None)

        # Ignore the optimizer-related keyword arguments, if given.
        kwds.pop('loc', None)
        kwds.pop('scale', None)
        kwds.pop('optimizer', None)
        if kwds:
            raise TypeError("Unknown arguments: %s." % kwds)

        if floc is not None and fscale is not None:
            # This check is for consistency with `rv_continuous.fit`.
            raise ValueError("All parameters fixed. There is nothing to "
                             "optimize.")

        data = np.asarray(data)

        # MLE for the uniform distribution
        # --------------------------------
        # The PDF is
        #
        #     f(x, loc, scale) = {1/scale  for loc <= x <= loc + scale
        #                        {0        otherwise}
        #
        # The likelihood function is
        #     L(x, loc, scale) = (1/scale)**n
        # where n is len(x), assuming loc <= x <= loc + scale for all x.
        # The log-likelihood is
        #     l(x, loc, scale) = -n*log(scale)
        # The log-likelihood is maximized by making scale as small as possible,
        # while keeping loc <= x <= loc + scale.   So if neither loc nor scale
        # are fixed, the log-likelihood is maximized by choosing
        #     loc = x.min()
        #     scale = x.ptp()
        # If loc is fixed, it must be less than or equal to x.min(), and then
        # the scale is
        #     scale = x.max() - loc
        # If scale is fixed, it must not be less than x.ptp().  If scale is
        # greater than x.ptp(), the solution is not unique.  Note that the
        # likelihood does not depend on loc, except for the requirement that
        # loc <= x <= loc + scale.  All choices of loc for which
        #     x.max() - scale <= loc <= x.min()
        # have the same log-likelihood.  In this case, we choose loc such that
        # the support is centered over the interval [data.min(), data.max()]:
        #     loc = x.min() = 0.5*(scale - x.ptp())

        if fscale is None:
            # scale is not fixed.
            if floc is None:
                # loc is not fixed, scale is not fixed.
                loc = data.min()
                scale = data.ptp()
            else:
                # loc is fixed, scale is not fixed.
                loc = floc
                scale = data.max() - loc
                if data.min() < loc:
                    raise FitDataError("uniform", lower=loc, upper=loc + scale)
        else:
            # loc is not fixed, scale is fixed.
            ptp = data.ptp()
            if ptp > fscale:
                raise FitUniformFixedScaleDataError(ptp=ptp, fscale=fscale)
            # If ptp < fscale, the ML estimate is not unique; see the comments
            # above.  We choose the distribution for which the support is
            # centered over the interval [data.min(), data.max()].
            loc = data.min() - 0.5*(fscale - ptp)
            scale = fscale

        # We expect the return values to be floating point, so ensure it
        # by explicitly converting to float.
        return float(loc), float(scale)


uniform = uniform_gen(a=0.0, b=1.0, name='uniform')


class vonmises_gen(rv_continuous):
    r"""A Von Mises continuous random variable.

    %(before_notes)s

    Notes
    -----
    If `x` is not in range or `loc` is not in range it assumes they are angles
    and converts them to [-\pi, \pi] equivalents.

    The probability density function for `vonmises` is:

    .. math::

        f(x, \kappa) = \frac{ \exp(\kappa \cos(x)) }{ 2 \pi I[0](\kappa) }

    for :math:`-\pi \le x \le \pi`, :math:`\kappa > 0`.

    `vonmises` takes :math:`\kappa` as a shape parameter.

    %(after_notes)s

    See Also
    --------
    vonmises_line : The same distribution, defined on a [-\pi, \pi] segment
                    of the real line.

    %(example)s

    """
    def _rvs(self, kappa):
        return self._random_state.vonmises(0.0, kappa, size=self._size)

    def _pdf(self, x, kappa):
        # vonmises.pdf(x, \kappa) = exp(\kappa * cos(x)) / (2*pi*I[0](\kappa))
        return np.exp(kappa * np.cos(x)) / (2*np.pi*sc.i0(kappa))

    def _cdf(self, x, kappa):
        return _stats.von_mises_cdf(kappa, x)

    def _stats_skip(self, kappa):
        return 0, None, 0, None

    def _entropy(self, kappa):
        return (-kappa * sc.i1(kappa) / sc.i0(kappa) +
                np.log(2 * np.pi * sc.i0(kappa)))


vonmises = vonmises_gen(name='vonmises')
vonmises_line = vonmises_gen(a=-np.pi, b=np.pi, name='vonmises_line')


class wald_gen(invgauss_gen):
    r"""A Wald continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `wald` is:

    .. math::

        f(x) = \frac{1}{\sqrt{2\pi x^3}} \exp(- \frac{ (x-1)^2 }{ 2x })

    for :math:`x > 0`.

    `wald` is a special case of `invgauss` with ``mu == 1``.

    %(after_notes)s

    %(example)s
    """
    _support_mask = rv_continuous._open_support_mask

    def _rvs(self):
        return self._random_state.wald(1.0, 1.0, size=self._size)

    def _pdf(self, x):
        # wald.pdf(x) = 1/sqrt(2*pi*x**3) * exp(-(x-1)**2/(2*x))
        return invgauss._pdf(x, 1.0)

    def _logpdf(self, x):
        return invgauss._logpdf(x, 1.0)

    def _cdf(self, x):
        return invgauss._cdf(x, 1.0)

    def _stats(self):
        return 1.0, 1.0, 3.0, 15.0


wald = wald_gen(a=0.0, name="wald")


class wrapcauchy_gen(rv_continuous):
    r"""A wrapped Cauchy continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `wrapcauchy` is:

    .. math::

        f(x, c) = \frac{1-c^2}{2\pi (1+c^2 - 2c \cos(x))}

    for :math:`0 \le x \le 2\pi`, :math:`0 < c < 1`.

    `wrapcauchy` takes :math:`c` as a shape parameter.

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, c):
        return (c > 0) & (c < 1)

    def _pdf(self, x, c):
        # wrapcauchy.pdf(x, c) = (1-c**2) / (2*pi*(1+c**2-2*c*cos(x)))
        return (1.0-c*c)/(2*np.pi*(1+c*c-2*c*np.cos(x)))

    def _cdf(self, x, c):
        output = np.zeros(x.shape, dtype=x.dtype)
        val = (1.0+c)/(1.0-c)
        c1 = x < np.pi
        c2 = 1-c1
        xp = np.extract(c1, x)
        xn = np.extract(c2, x)
        if np.any(xn):
            valn = np.extract(c2, np.ones_like(x)*val)
            xn = 2*np.pi - xn
            yn = np.tan(xn/2.0)
            on = 1.0-1.0/np.pi*np.arctan(valn*yn)
            np.place(output, c2, on)
        if np.any(xp):
            valp = np.extract(c1, np.ones_like(x)*val)
            yp = np.tan(xp/2.0)
            op = 1.0/np.pi*np.arctan(valp*yp)
            np.place(output, c1, op)
        return output

    def _ppf(self, q, c):
        val = (1.0-c)/(1.0+c)
        rcq = 2*np.arctan(val*np.tan(np.pi*q))
        rcmq = 2*np.pi-2*np.arctan(val*np.tan(np.pi*(1-q)))
        return np.where(q < 1.0/2, rcq, rcmq)

    def _entropy(self, c):
        return np.log(2*np.pi*(1-c*c))


wrapcauchy = wrapcauchy_gen(a=0.0, b=2*np.pi, name='wrapcauchy')


class gennorm_gen(rv_continuous):
    r"""A generalized normal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `gennorm` is [1]_::

                                     beta
        gennorm.pdf(x, beta) =  ---------------  exp(-|x|**beta)
                                2 gamma(1/beta)

    `gennorm` takes :math:`\beta` as a shape parameter.
    For :math:`\beta = 1`, it is identical to a Laplace distribution.
    For ``\beta = 2``, it is identical to a normal distribution
    (with :math:`scale=1/\sqrt{2}`).

    See Also
    --------
    laplace : Laplace distribution
    norm : normal distribution

    References
    ----------

    .. [1] "Generalized normal distribution, Version 1",
           https://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1

    %(example)s

    """

    def _pdf(self, x, beta):
        return np.exp(self._logpdf(x, beta))

    def _logpdf(self, x, beta):
        return np.log(0.5*beta) - sc.gammaln(1.0/beta) - abs(x)**beta

    def _cdf(self, x, beta):
        c = 0.5 * np.sign(x)
        # evaluating (.5 + c) first prevents numerical cancellation
        return (0.5 + c) - c * sc.gammaincc(1.0/beta, abs(x)**beta)

    def _ppf(self, x, beta):
        c = np.sign(x - 0.5)
        # evaluating (1. + c) first prevents numerical cancellation
        return c * sc.gammainccinv(1.0/beta, (1.0 + c) - 2.0*c*x)**(1.0/beta)

    def _sf(self, x, beta):
        return self._cdf(-x, beta)

    def _isf(self, x, beta):
        return -self._ppf(x, beta)

    def _stats(self, beta):
        c1, c3, c5 = sc.gammaln([1.0/beta, 3.0/beta, 5.0/beta])
        return 0., np.exp(c3 - c1), 0., np.exp(c5 + c1 - 2.0*c3) - 3.

    def _entropy(self, beta):
        return 1. / beta - np.log(.5 * beta) + sc.gammaln(1. / beta)


gennorm = gennorm_gen(name='gennorm')


class halfgennorm_gen(rv_continuous):
    r"""The upper half of a generalized normal continuous random variable.

    %(before_notes)s

    Notes
    -----
    The probability density function for `halfgennorm` is:

    .. math::

        f(x, \beta) = \frac{\beta}{\gamma(1/\beta)} \exp(-|x|^\beta)

    `gennorm` takes :math:`\beta` as a shape parameter.
    For :math:`\beta = 1`, it is identical to an exponential distribution.
    For :math:`\beta = 2`, it is identical to a half normal distribution
    (with :math:`scale=1/\sqrt{2}`).

    See Also
    --------
    gennorm : generalized normal distribution
    expon : exponential distribution
    halfnorm : half normal distribution

    References
    ----------

    .. [1] "Generalized normal distribution, Version 1",
           https://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1

    %(example)s

    """

    def _pdf(self, x, beta):
        #                                 beta
        # halfgennorm.pdf(x, beta) =  -------------  exp(-|x|**beta)
        #                             gamma(1/beta)
        return np.exp(self._logpdf(x, beta))

    def _logpdf(self, x, beta):
        return np.log(beta) - sc.gammaln(1.0/beta) - x**beta

    def _cdf(self, x, beta):
        return sc.gammainc(1.0/beta, x**beta)

    def _ppf(self, x, beta):
        return sc.gammaincinv(1.0/beta, x)**(1.0/beta)

    def _sf(self, x, beta):
        return sc.gammaincc(1.0/beta, x**beta)

    def _isf(self, x, beta):
        return sc.gammainccinv(1.0/beta, x)**(1.0/beta)

    def _entropy(self, beta):
        return 1.0/beta - np.log(beta) + sc.gammaln(1.0/beta)


halfgennorm = halfgennorm_gen(a=0, name='halfgennorm')


class crystalball_gen(rv_continuous):
    r"""
    Crystalball distribution

    %(before_notes)s

    Notes
    -----
    The probability density function for `crystalball` is:

    .. math::

        f(x, \beta, m) =  \begin{cases}
                            N \exp(-x^2 / 2),  &\text{for } x > -\beta\\
                            N A (B - x)^{-m}  &\text{for } x \le -\beta
                          \end{cases}

    where :math:`A = (m / |beta|)**n * exp(-beta**2 / 2)`,
    :math:`B = m/|beta| - |beta|` and :math:`N` is a normalisation constant.

    `crystalball` takes :math:`\beta` and :math:`m` as shape parameters.
    :math:`\beta` defines the point where the pdf changes from a power-law to a
    gaussian distribution :math:`m` is power of the power-law tail.

    References
    ----------
    .. [1] "Crystal Ball Function",
           https://en.wikipedia.org/wiki/Crystal_Ball_function

    %(after_notes)s

    .. versionadded:: 0.19.0

    %(example)s
    """
    def _pdf(self, x, beta, m):
        """
        Return PDF of the crystalball function.

                                            --
                                           | exp(-x**2 / 2),  for x > -beta
        crystalball.pdf(x, beta, m) =  N * |
                                           | A * (B - x)**(-m), for x <= -beta
                                            --
        """
        N = 1.0 / (m/beta / (m-1) * np.exp(-beta**2 / 2.0) + _norm_pdf_C * _norm_cdf(beta))
        rhs = lambda x, beta, m: np.exp(-x**2 / 2)
        lhs = lambda x, beta, m: (m/beta)**m * np.exp(-beta**2 / 2.0) * (m/beta - beta - x)**(-m)
        return N * _lazywhere(np.atleast_1d(x > -beta), (x, beta, m), f=rhs, f2=lhs)

    def _cdf(self, x, beta, m):
        """
        Return CDF of the crystalball function
        """
        N = 1.0 / (m/beta / (m-1) * np.exp(-beta**2 / 2.0) + _norm_pdf_C * _norm_cdf(beta))
        rhs = lambda x, beta, m: (m/beta) * np.exp(-beta**2 / 2.0) / (m-1) + _norm_pdf_C * (_norm_cdf(x) - _norm_cdf(-beta))
        lhs = lambda x, beta, m: (m/beta)**m * np.exp(-beta**2 / 2.0) * (m/beta - beta - x)**(-m+1) / (m-1)
        return N * _lazywhere(np.atleast_1d(x > -beta), (x, beta, m), f=rhs, f2=lhs)

    def _munp(self, n, beta, m):
        """
        Returns the n-th non-central moment of the crystalball function.
        """
        N = 1.0 / (m/beta / (m-1) * np.exp(-beta**2 / 2.0) + _norm_pdf_C * _norm_cdf(beta))

        def n_th_moment(n, beta, m):
            """
            Returns n-th moment. Defined only if n+1 < m
            Function cannot broadcast due to the loop over n
            """
            A = (m/beta)**m * np.exp(-beta**2 / 2.0)
            B = m/beta - beta
            rhs = 2**((n-1)/2.0) * sc.gamma((n+1)/2) * (1.0 + (-1)**n * sc.gammainc((n+1)/2, beta**2 / 2))
            lhs = np.zeros(rhs.shape)
            for k in range(n + 1):
                lhs += sc.binom(n, k) * B**(n-k) * (-1)**k / (m - k - 1) * (m/beta)**(-m + k + 1)
            return A * lhs + rhs

        return N * _lazywhere(np.atleast_1d(n + 1 < m),
                              (n, beta, m),
                              np.vectorize(n_th_moment, otypes=[np.float]),
                              np.inf)

    def _argcheck(self, beta, m):
        """
        In HEP crystal-ball is also defined for m = 1 (see plot on wikipedia)
        But the function doesn't have a finite integral in this corner case,
        and isn't a PDF anymore (but can still be used on a finite range).
        Here we restrict the function to m > 1.
        In addition we restrict beta to be positive
        """
        return (m > 1) & (beta > 0)


crystalball = crystalball_gen(name='crystalball', longname="A Crystalball Function")


def _argus_phi(chi):
    """
    Utility function for the argus distribution
    used in the CDF and norm of the Argus Funktion
    """
    return _norm_cdf(chi) - chi * _norm_pdf(chi) - 0.5


class argus_gen(rv_continuous):
    r"""
    Argus distribution

    %(before_notes)s

    Notes
    -----
    The probability density function for `argus` is:

    .. math::

        f(x, \chi) = \frac{\chi^3}{\sqrt{2\pi} \Psi(\chi)} x \sqrt{1-x^2}
                     \exp(- 0.5 \chi^2 (1 - x^2))

        where:

    .. math::

        \Psi(\chi) = \Phi(\chi) - \chi \phi(\chi) - 1/2

    with :math:`\Phi` and :math:`\phi` being the CDF and PDF of a standard
    normal distribution, respectively.

    `argus` takes :math:`\chi` as shape a parameter.

    References
    ----------

    .. [1] "ARGUS distribution",
           https://en.wikipedia.org/wiki/ARGUS_distribution

    %(after_notes)s

    .. versionadded:: 0.19.0

    %(example)s
    """
    def _pdf(self, x, chi):
        """
        Return PDF of the argus function

        argus.pdf(x, chi) = chi**3 / (sqrt(2*pi) * Psi(chi)) * x *
                            sqrt(1-x**2) * exp(- 0.5 * chi**2 * (1 - x**2))
        """
        y = 1.0 - x**2
        return chi**3 / (_norm_pdf_C * _argus_phi(chi)) * x * np.sqrt(y) * np.exp(-chi**2 * y / 2)

    def _cdf(self, x, chi):
        """
        Return CDF of the argus function
        """
        return 1.0 - self._sf(x, chi)

    def _sf(self, x, chi):
        """
        Return survival function of the argus function
        """
        return _argus_phi(chi * np.sqrt(1 - x**2)) / _argus_phi(chi)


argus = argus_gen(name='argus', longname="An Argus Function", a=0.0, b=1.0)


class rv_histogram(rv_continuous):
    """
    Generates a distribution given by a histogram.
    This is useful to generate a template distribution from a binned
    datasample.

    As a subclass of the `rv_continuous` class, `rv_histogram` inherits from it
    a collection of generic methods (see `rv_continuous` for the full list),
    and implements them based on the properties of the provided binned
    datasample.

    Parameters
    ----------
    histogram : tuple of array_like
      Tuple containing two array_like objects
      The first containing the content of n bins
      The second containing the (n+1) bin boundaries
      In particular the return value np.histogram is accepted

    Notes
    -----
    There are no additional shape parameters except for the loc and scale.
    The pdf is defined as a stepwise function from the provided histogram
    The cdf is a linear interpolation of the pdf.

    .. versionadded:: 0.19.0

    Examples
    --------

    Create a scipy.stats distribution from a numpy histogram

    >>> import scipy.stats
    >>> import numpy as np
    >>> data = scipy.stats.norm.rvs(size=100000, loc=0, scale=1.5, random_state=123)
    >>> hist = np.histogram(data, bins=100)
    >>> hist_dist = scipy.stats.rv_histogram(hist)

    Behaves like an ordinary scipy rv_continuous distribution

    >>> hist_dist.pdf(1.0)
    0.20538577847618705
    >>> hist_dist.cdf(2.0)
    0.90818568543056499

    PDF is zero above (below) the highest (lowest) bin of the histogram,
    defined by the max (min) of the original dataset

    >>> hist_dist.pdf(np.max(data))
    0.0
    >>> hist_dist.cdf(np.max(data))
    1.0
    >>> hist_dist.pdf(np.min(data))
    7.7591907244498314e-05
    >>> hist_dist.cdf(np.min(data))
    0.0

    PDF and CDF follow the histogram

    >>> import matplotlib.pyplot as plt
    >>> X = np.linspace(-5.0, 5.0, 100)
    >>> plt.title("PDF from Template")
    >>> plt.hist(data, density=True, bins=100)
    >>> plt.plot(X, hist_dist.pdf(X), label='PDF')
    >>> plt.plot(X, hist_dist.cdf(X), label='CDF')
    >>> plt.show()

    """
    _support_mask = rv_continuous._support_mask

    def __init__(self, histogram, *args, **kwargs):
        """
        Create a new distribution using the given histogram

        Parameters
        ----------
        histogram : tuple of array_like
          Tuple containing two array_like objects
          The first containing the content of n bins
          The second containing the (n+1) bin boundaries
          In particular the return value np.histogram is accepted
        """
        self._histogram = histogram
        if len(histogram) != 2:
            raise ValueError("Expected length 2 for parameter histogram")
        self._hpdf = np.asarray(histogram[0])
        self._hbins = np.asarray(histogram[1])
        if len(self._hpdf) + 1 != len(self._hbins):
            raise ValueError("Number of elements in histogram content "
                             "and histogram boundaries do not match, "
                             "expected n and n+1.")
        self._hbin_widths = self._hbins[1:] - self._hbins[:-1]
        self._hpdf = self._hpdf / float(np.sum(self._hpdf * self._hbin_widths))
        self._hcdf = np.cumsum(self._hpdf * self._hbin_widths)
        self._hpdf = np.hstack([0.0, self._hpdf, 0.0])
        self._hcdf = np.hstack([0.0, self._hcdf])
        # Set support
        kwargs['a'] = self._hbins[0]
        kwargs['b'] = self._hbins[-1]
        super(rv_histogram, self).__init__(*args, **kwargs)

    def _pdf(self, x):
        """
        PDF of the histogram
        """
        return self._hpdf[np.searchsorted(self._hbins, x, side='right')]

    def _cdf(self, x):
        """
        CDF calculated from the histogram
        """
        return np.interp(x, self._hbins, self._hcdf)

    def _ppf(self, x):
        """
        Percentile function calculated from the histogram
        """
        return np.interp(x, self._hcdf, self._hbins)

    def _munp(self, n):
        """Compute the n-th non-central moment."""
        integrals = (self._hbins[1:]**(n+1) - self._hbins[:-1]**(n+1)) / (n+1)
        return np.sum(self._hpdf[1:-1] * integrals)

    def _entropy(self):
        """Compute entropy of distribution"""
        res = _lazywhere(self._hpdf[1:-1] > 0.0,
                         (self._hpdf[1:-1],),
                         np.log,
                         0.0)
        return -np.sum(self._hpdf[1:-1] * res * self._hbin_widths)

    def _updated_ctor_param(self):
        """
        Set the histogram as additional constructor argument
        """
        dct = super(rv_histogram, self)._updated_ctor_param()
        dct['histogram'] = self._histogram
        return dct


# Collect names of classes and objects in this module.
pairs = list(globals().items())
_distn_names, _distn_gen_names = get_distribution_names(pairs, rv_continuous)

__all__ = _distn_names + _distn_gen_names + ['rv_histogram']
