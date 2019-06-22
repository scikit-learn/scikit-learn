#
# Author:  Travis Oliphant  2002-2011 with contributions from
#          SciPy Developers 2004-2011
#
from __future__ import division, print_function, absolute_import

from scipy import special
from scipy.special import entr, logsumexp, betaln, gammaln as gamln
from scipy._lib._numpy_compat import broadcast_to
from scipy._lib._util import _lazywhere

from numpy import floor, ceil, log, exp, sqrt, log1p, expm1, tanh, cosh, sinh

import numpy as np

from ._distn_infrastructure import (
        rv_discrete, _ncx2_pdf, _ncx2_cdf, get_distribution_names)


class binom_gen(rv_discrete):
    r"""A binomial discrete random variable.

    %(before_notes)s

    Notes
    -----
    The probability mass function for `binom` is:

    .. math::

       f(k) = \binom{n}{k} p^k (1-p)^{n-k}

    for ``k`` in ``{0, 1,..., n}``.

    `binom` takes ``n`` and ``p`` as shape parameters.

    %(after_notes)s

    %(example)s

    """
    def _rvs(self, n, p):
        return self._random_state.binomial(n, p, self._size)

    def _argcheck(self, n, p):
        return (n >= 0) & (p >= 0) & (p <= 1)

    def _get_support(self, n, p):
        return self.a, n

    def _logpmf(self, x, n, p):
        k = floor(x)
        combiln = (gamln(n+1) - (gamln(k+1) + gamln(n-k+1)))
        return combiln + special.xlogy(k, p) + special.xlog1py(n-k, -p)

    def _pmf(self, x, n, p):
        # binom.pmf(k) = choose(n, k) * p**k * (1-p)**(n-k)
        return exp(self._logpmf(x, n, p))

    def _cdf(self, x, n, p):
        k = floor(x)
        vals = special.bdtr(k, n, p)
        return vals

    def _sf(self, x, n, p):
        k = floor(x)
        return special.bdtrc(k, n, p)

    def _ppf(self, q, n, p):
        vals = ceil(special.bdtrik(q, n, p))
        vals1 = np.maximum(vals - 1, 0)
        temp = special.bdtr(vals1, n, p)
        return np.where(temp >= q, vals1, vals)

    def _stats(self, n, p, moments='mv'):
        q = 1.0 - p
        mu = n * p
        var = n * p * q
        g1, g2 = None, None
        if 's' in moments:
            g1 = (q - p) / sqrt(var)
        if 'k' in moments:
            g2 = (1.0 - 6*p*q) / var
        return mu, var, g1, g2

    def _entropy(self, n, p):
        k = np.r_[0:n + 1]
        vals = self._pmf(k, n, p)
        return np.sum(entr(vals), axis=0)


binom = binom_gen(name='binom')


class bernoulli_gen(binom_gen):
    r"""A Bernoulli discrete random variable.

    %(before_notes)s

    Notes
    -----
    The probability mass function for `bernoulli` is:

    .. math::

       f(k) = \begin{cases}1-p  &\text{if } k = 0\\
                           p    &\text{if } k = 1\end{cases}

    for :math:`k` in :math:`\{0, 1\}`.

    `bernoulli` takes :math:`p` as shape parameter.

    %(after_notes)s

    %(example)s

    """
    def _rvs(self, p):
        return binom_gen._rvs(self, 1, p)

    def _argcheck(self, p):
        return (p >= 0) & (p <= 1)

    def _get_support(self, p):
        # Overrides binom_gen._get_support!x
        return self.a, self.b

    def _logpmf(self, x, p):
        return binom._logpmf(x, 1, p)

    def _pmf(self, x, p):
        # bernoulli.pmf(k) = 1-p  if k = 0
        #                  = p    if k = 1
        return binom._pmf(x, 1, p)

    def _cdf(self, x, p):
        return binom._cdf(x, 1, p)

    def _sf(self, x, p):
        return binom._sf(x, 1, p)

    def _ppf(self, q, p):
        return binom._ppf(q, 1, p)

    def _stats(self, p):
        return binom._stats(1, p)

    def _entropy(self, p):
        return entr(p) + entr(1-p)


bernoulli = bernoulli_gen(b=1, name='bernoulli')


class nbinom_gen(rv_discrete):
    r"""A negative binomial discrete random variable.

    %(before_notes)s

    Notes
    -----
    Negative binomial distribution describes a sequence of i.i.d. Bernoulli
    trials, repeated until a predefined, non-random number of successes occurs.

    The probability mass function of the number of failures for `nbinom` is:

    .. math::

       f(k) = \binom{k+n-1}{n-1} p^n (1-p)^k

    for :math:`k \ge 0`.

    `nbinom` takes :math:`n` and :math:`p` as shape parameters where n is the
    number of successes, whereas p is the probability of a single success.

    %(after_notes)s

    %(example)s

    """
    def _rvs(self, n, p):
        return self._random_state.negative_binomial(n, p, self._size)

    def _argcheck(self, n, p):
        return (n > 0) & (p >= 0) & (p <= 1)

    def _pmf(self, x, n, p):
        # nbinom.pmf(k) = choose(k+n-1, n-1) * p**n * (1-p)**k
        return exp(self._logpmf(x, n, p))

    def _logpmf(self, x, n, p):
        coeff = gamln(n+x) - gamln(x+1) - gamln(n)
        return coeff + n*log(p) + special.xlog1py(x, -p)

    def _cdf(self, x, n, p):
        k = floor(x)
        return special.betainc(n, k+1, p)

    def _sf_skip(self, x, n, p):
        # skip because special.nbdtrc doesn't work for 0<n<1
        k = floor(x)
        return special.nbdtrc(k, n, p)

    def _ppf(self, q, n, p):
        vals = ceil(special.nbdtrik(q, n, p))
        vals1 = (vals-1).clip(0.0, np.inf)
        temp = self._cdf(vals1, n, p)
        return np.where(temp >= q, vals1, vals)

    def _stats(self, n, p):
        Q = 1.0 / p
        P = Q - 1.0
        mu = n*P
        var = n*P*Q
        g1 = (Q+P)/sqrt(n*P*Q)
        g2 = (1.0 + 6*P*Q) / (n*P*Q)
        return mu, var, g1, g2


nbinom = nbinom_gen(name='nbinom')


class geom_gen(rv_discrete):
    r"""A geometric discrete random variable.

    %(before_notes)s

    Notes
    -----
    The probability mass function for `geom` is:

    .. math::

        f(k) = (1-p)^{k-1} p

    for :math:`k \ge 1`.

    `geom` takes :math:`p` as shape parameter.

    %(after_notes)s

    See Also
    --------
    planck

    %(example)s

    """
    def _rvs(self, p):
        return self._random_state.geometric(p, size=self._size)

    def _argcheck(self, p):
        return (p <= 1) & (p >= 0)

    def _pmf(self, k, p):
        return np.power(1-p, k-1) * p

    def _logpmf(self, k, p):
        return special.xlog1py(k - 1, -p) + log(p)

    def _cdf(self, x, p):
        k = floor(x)
        return -expm1(log1p(-p)*k)

    def _sf(self, x, p):
        return np.exp(self._logsf(x, p))

    def _logsf(self, x, p):
        k = floor(x)
        return k*log1p(-p)

    def _ppf(self, q, p):
        vals = ceil(log1p(-q) / log1p(-p))
        temp = self._cdf(vals-1, p)
        return np.where((temp >= q) & (vals > 0), vals-1, vals)

    def _stats(self, p):
        mu = 1.0/p
        qr = 1.0-p
        var = qr / p / p
        g1 = (2.0-p) / sqrt(qr)
        g2 = np.polyval([1, -6, 6], p)/(1.0-p)
        return mu, var, g1, g2


geom = geom_gen(a=1, name='geom', longname="A geometric")


class hypergeom_gen(rv_discrete):
    r"""A hypergeometric discrete random variable.

    The hypergeometric distribution models drawing objects from a bin.
    `M` is the total number of objects, `n` is total number of Type I objects.
    The random variate represents the number of Type I objects in `N` drawn
    without replacement from the total population.

    %(before_notes)s

    Notes
    -----
    The symbols used to denote the shape parameters (`M`, `n`, and `N`) are not
    universally accepted.  See the Examples for a clarification of the
    definitions used here.

    The probability mass function is defined as,

    .. math:: p(k, M, n, N) = \frac{\binom{n}{k} \binom{M - n}{N - k}}
                                   {\binom{M}{N}}

    for :math:`k \in [\max(0, N - M + n), \min(n, N)]`, where the binomial
    coefficients are defined as,

    .. math:: \binom{n}{k} \equiv \frac{n!}{k! (n - k)!}.

    %(after_notes)s

    Examples
    --------
    >>> from scipy.stats import hypergeom
    >>> import matplotlib.pyplot as plt

    Suppose we have a collection of 20 animals, of which 7 are dogs.  Then if
    we want to know the probability of finding a given number of dogs if we
    choose at random 12 of the 20 animals, we can initialize a frozen
    distribution and plot the probability mass function:

    >>> [M, n, N] = [20, 7, 12]
    >>> rv = hypergeom(M, n, N)
    >>> x = np.arange(0, n+1)
    >>> pmf_dogs = rv.pmf(x)

    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111)
    >>> ax.plot(x, pmf_dogs, 'bo')
    >>> ax.vlines(x, 0, pmf_dogs, lw=2)
    >>> ax.set_xlabel('# of dogs in our group of chosen animals')
    >>> ax.set_ylabel('hypergeom PMF')
    >>> plt.show()

    Instead of using a frozen distribution we can also use `hypergeom`
    methods directly.  To for example obtain the cumulative distribution
    function, use:

    >>> prb = hypergeom.cdf(x, M, n, N)

    And to generate random numbers:

    >>> R = hypergeom.rvs(M, n, N, size=10)

    """
    def _rvs(self, M, n, N):
        return self._random_state.hypergeometric(n, M-n, N, size=self._size)

    def _get_support(self, M, n, N):
        return np.maximum(N-(M-n), 0), np.minimum(n, N)

    def _argcheck(self, M, n, N):
        cond = (M > 0) & (n >= 0) & (N >= 0)
        cond &= (n <= M) & (N <= M)
        return cond

    def _logpmf(self, k, M, n, N):
        tot, good = M, n
        bad = tot - good
        result = (betaln(good+1, 1) + betaln(bad+1, 1) + betaln(tot-N+1, N+1) -
                  betaln(k+1, good-k+1) - betaln(N-k+1, bad-N+k+1) -
                  betaln(tot+1, 1))
        return result

    def _pmf(self, k, M, n, N):
        # same as the following but numerically more precise
        # return comb(good, k) * comb(bad, N-k) / comb(tot, N)
        return exp(self._logpmf(k, M, n, N))

    def _stats(self, M, n, N):
        # tot, good, sample_size = M, n, N
        # "wikipedia".replace('N', 'M').replace('n', 'N').replace('K', 'n')
        M, n, N = 1.*M, 1.*n, 1.*N
        m = M - n
        p = n/M
        mu = N*p

        var = m*n*N*(M - N)*1.0/(M*M*(M-1))
        g1 = (m - n)*(M-2*N) / (M-2.0) * sqrt((M-1.0) / (m*n*N*(M-N)))

        g2 = M*(M+1) - 6.*N*(M-N) - 6.*n*m
        g2 *= (M-1)*M*M
        g2 += 6.*n*N*(M-N)*m*(5.*M-6)
        g2 /= n * N * (M-N) * m * (M-2.) * (M-3.)
        return mu, var, g1, g2

    def _entropy(self, M, n, N):
        k = np.r_[N - (M - n):min(n, N) + 1]
        vals = self.pmf(k, M, n, N)
        return np.sum(entr(vals), axis=0)

    def _sf(self, k, M, n, N):
        # This for loop is needed because `k` can be an array. If that's the
        # case, the sf() method makes M, n and N arrays of the same shape. We
        # therefore unpack all inputs args, so we can do the manual
        # integration.
        res = []
        for quant, tot, good, draw in zip(k, M, n, N):
            # Manual integration over probability mass function. More accurate
            # than integrate.quad.
            k2 = np.arange(quant + 1, draw + 1)
            res.append(np.sum(self._pmf(k2, tot, good, draw)))
        return np.asarray(res)

    def _logsf(self, k, M, n, N):
        res = []
        for quant, tot, good, draw in zip(k, M, n, N):
            if (quant + 0.5) * (tot + 0.5) < (good - 0.5) * (draw - 0.5):
                # Less terms to sum if we calculate log(1-cdf)
                res.append(log1p(-exp(self.logcdf(quant, tot, good, draw))))
            else:
                # Integration over probability mass function using logsumexp
                k2 = np.arange(quant + 1, draw + 1)
                res.append(logsumexp(self._logpmf(k2, tot, good, draw)))
        return np.asarray(res)

    def _logcdf(self, k, M, n, N):
        res = []
        for quant, tot, good, draw in zip(k, M, n, N):
            if (quant + 0.5) * (tot + 0.5) > (good - 0.5) * (draw - 0.5):
                # Less terms to sum if we calculate log(1-sf)
                res.append(log1p(-exp(self.logsf(quant, tot, good, draw))))
            else:
                # Integration over probability mass function using logsumexp
                k2 = np.arange(0, quant + 1)
                res.append(logsumexp(self._logpmf(k2, tot, good, draw)))
        return np.asarray(res)


hypergeom = hypergeom_gen(name='hypergeom')


# FIXME: Fails _cdfvec
class logser_gen(rv_discrete):
    r"""A Logarithmic (Log-Series, Series) discrete random variable.

    %(before_notes)s

    Notes
    -----
    The probability mass function for `logser` is:

    .. math::

        f(k) = - \frac{p^k}{k \log(1-p)}

    for :math:`k \ge 1`.

    `logser` takes :math:`p` as shape parameter.

    %(after_notes)s

    %(example)s

    """
    def _rvs(self, p):
        # looks wrong for p>0.5, too few k=1
        # trying to use generic is worse, no k=1 at all
        return self._random_state.logseries(p, size=self._size)

    def _argcheck(self, p):
        return (p > 0) & (p < 1)

    def _pmf(self, k, p):
        # logser.pmf(k) = - p**k / (k*log(1-p))
        return -np.power(p, k) * 1.0 / k / special.log1p(-p)

    def _stats(self, p):
        r = special.log1p(-p)
        mu = p / (p - 1.0) / r
        mu2p = -p / r / (p - 1.0)**2
        var = mu2p - mu*mu
        mu3p = -p / r * (1.0+p) / (1.0 - p)**3
        mu3 = mu3p - 3*mu*mu2p + 2*mu**3
        g1 = mu3 / np.power(var, 1.5)

        mu4p = -p / r * (
            1.0 / (p-1)**2 - 6*p / (p - 1)**3 + 6*p*p / (p-1)**4)
        mu4 = mu4p - 4*mu3p*mu + 6*mu2p*mu*mu - 3*mu**4
        g2 = mu4 / var**2 - 3.0
        return mu, var, g1, g2


logser = logser_gen(a=1, name='logser', longname='A logarithmic')


class poisson_gen(rv_discrete):
    r"""A Poisson discrete random variable.

    %(before_notes)s

    Notes
    -----
    The probability mass function for `poisson` is:

    .. math::

        f(k) = \exp(-\mu) \frac{\mu^k}{k!}

    for :math:`k \ge 0`.

    `poisson` takes :math:`\mu` as shape parameter.

    %(after_notes)s

    %(example)s

    """

    # Override rv_discrete._argcheck to allow mu=0.
    def _argcheck(self, mu):
        return mu >= 0

    def _rvs(self, mu):
        return self._random_state.poisson(mu, self._size)

    def _logpmf(self, k, mu):
        Pk = special.xlogy(k, mu) - gamln(k + 1) - mu
        return Pk

    def _pmf(self, k, mu):
        # poisson.pmf(k) = exp(-mu) * mu**k / k!
        return exp(self._logpmf(k, mu))

    def _cdf(self, x, mu):
        k = floor(x)
        return special.pdtr(k, mu)

    def _sf(self, x, mu):
        k = floor(x)
        return special.pdtrc(k, mu)

    def _ppf(self, q, mu):
        vals = ceil(special.pdtrik(q, mu))
        vals1 = np.maximum(vals - 1, 0)
        temp = special.pdtr(vals1, mu)
        return np.where(temp >= q, vals1, vals)

    def _stats(self, mu):
        var = mu
        tmp = np.asarray(mu)
        mu_nonzero = tmp > 0
        g1 = _lazywhere(mu_nonzero, (tmp,), lambda x: sqrt(1.0/x), np.inf)
        g2 = _lazywhere(mu_nonzero, (tmp,), lambda x: 1.0/x, np.inf)
        return mu, var, g1, g2


poisson = poisson_gen(name="poisson", longname='A Poisson')


class planck_gen(rv_discrete):
    r"""A Planck discrete exponential random variable.

    %(before_notes)s

    Notes
    -----
    The probability mass function for `planck` is:

    .. math::

        f(k) = (1-\exp(-\lambda)) \exp(-\lambda k)

    for :math:`k \ge 0` and :math:`\lambda > 0`.

    `planck` takes :math:`\lambda` as shape parameter. The Planck distribution
    can be written as a geometric distribution (`geom`) with
    :math:`p = 1 - \exp(-\lambda)` shifted by `loc = -1`.

    %(after_notes)s

    See Also
    --------
    geom

    %(example)s

    """
    def _argcheck(self, lambda_):
        return lambda_ > 0

    def _pmf(self, k, lambda_):
        return -expm1(-lambda_)*exp(-lambda_*k)

    def _cdf(self, x, lambda_):
        k = floor(x)
        return -expm1(-lambda_*(k+1))

    def _sf(self, x, lambda_):
        return exp(self._logsf(x, lambda_))

    def _logsf(self, x, lambda_):
        k = floor(x)
        return -lambda_*(k+1)

    def _ppf(self, q, lambda_):
        vals = ceil(-1.0/lambda_ * log1p(-q)-1)
        vals1 = (vals-1).clip(*(self._get_support(lambda_)))
        temp = self._cdf(vals1, lambda_)
        return np.where(temp >= q, vals1, vals)

    def _rvs(self, lambda_):
        # use relation to geometric distribution for sampling
        p = -expm1(-lambda_)
        return self._random_state.geometric(p, size=self._size) - 1.0

    def _stats(self, lambda_):
        mu = 1/expm1(lambda_)
        var = exp(-lambda_)/(expm1(-lambda_))**2
        g1 = 2*cosh(lambda_/2.0)
        g2 = 4+2*cosh(lambda_)
        return mu, var, g1, g2

    def _entropy(self, lambda_):
        C = -expm1(-lambda_)
        return lambda_*exp(-lambda_)/C - log(C)


planck = planck_gen(a=0, name='planck', longname='A discrete exponential ')


class boltzmann_gen(rv_discrete):
    r"""A Boltzmann (Truncated Discrete Exponential) random variable.

    %(before_notes)s

    Notes
    -----
    The probability mass function for `boltzmann` is:

    .. math::

        f(k) = (1-\exp(-\lambda)) \exp(-\lambda k) / (1-\exp(-\lambda N))

    for :math:`k = 0,..., N-1`.

    `boltzmann` takes :math:`\lambda > 0` and :math:`N > 0` as shape parameters.

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, lambda_, N):
        return (lambda_ > 0) & (N > 0)

    def _get_support(self, lambda_, N):
        return self.a, N - 1

    def _pmf(self, k, lambda_, N):
        # boltzmann.pmf(k) =
        #               (1-exp(-lambda_)*exp(-lambda_*k)/(1-exp(-lambda_*N))
        fact = (1-exp(-lambda_))/(1-exp(-lambda_*N))
        return fact*exp(-lambda_*k)

    def _cdf(self, x, lambda_, N):
        k = floor(x)
        return (1-exp(-lambda_*(k+1)))/(1-exp(-lambda_*N))

    def _ppf(self, q, lambda_, N):
        qnew = q*(1-exp(-lambda_*N))
        vals = ceil(-1.0/lambda_ * log(1-qnew)-1)
        vals1 = (vals-1).clip(0.0, np.inf)
        temp = self._cdf(vals1, lambda_, N)
        return np.where(temp >= q, vals1, vals)

    def _stats(self, lambda_, N):
        z = exp(-lambda_)
        zN = exp(-lambda_*N)
        mu = z/(1.0-z)-N*zN/(1-zN)
        var = z/(1.0-z)**2 - N*N*zN/(1-zN)**2
        trm = (1-zN)/(1-z)
        trm2 = (z*trm**2 - N*N*zN)
        g1 = z*(1+z)*trm**3 - N**3*zN*(1+zN)
        g1 = g1 / trm2**(1.5)
        g2 = z*(1+4*z+z*z)*trm**4 - N**4 * zN*(1+4*zN+zN*zN)
        g2 = g2 / trm2 / trm2
        return mu, var, g1, g2


boltzmann = boltzmann_gen(name='boltzmann', a=0,
                          longname='A truncated discrete exponential ')


class randint_gen(rv_discrete):
    r"""A uniform discrete random variable.

    %(before_notes)s

    Notes
    -----
    The probability mass function for `randint` is:

    .. math::

        f(k) = \frac{1}{high - low}

    for ``k = low, ..., high - 1``.

    `randint` takes ``low`` and ``high`` as shape parameters.

    %(after_notes)s

    %(example)s

    """
    def _argcheck(self, low, high):
        return (high > low)

    def _get_support(self, low, high):
        return low, high-1

    def _pmf(self, k, low, high):
        # randint.pmf(k) = 1./(high - low)
        p = np.ones_like(k) / (high - low)
        return np.where((k >= low) & (k < high), p, 0.)

    def _cdf(self, x, low, high):
        k = floor(x)
        return (k - low + 1.) / (high - low)

    def _ppf(self, q, low, high):
        vals = ceil(q * (high - low) + low) - 1
        vals1 = (vals - 1).clip(low, high)
        temp = self._cdf(vals1, low, high)
        return np.where(temp >= q, vals1, vals)

    def _stats(self, low, high):
        m2, m1 = np.asarray(high), np.asarray(low)
        mu = (m2 + m1 - 1.0) / 2
        d = m2 - m1
        var = (d*d - 1) / 12.0
        g1 = 0.0
        g2 = -6.0/5.0 * (d*d + 1.0) / (d*d - 1.0)
        return mu, var, g1, g2

    def _rvs(self, low, high):
        """An array of *size* random integers >= ``low`` and < ``high``."""
        if self._size is not None:
            # NumPy's RandomState.randint() doesn't broadcast its arguments.
            # Use `broadcast_to()` to extend the shapes of low and high
            # up to self._size.  Then we can use the numpy.vectorize'd
            # randint without needing to pass it a `size` argument.
            low = broadcast_to(low, self._size)
            high = broadcast_to(high, self._size)
        randint = np.vectorize(self._random_state.randint, otypes=[np.int_])
        return randint(low, high)

    def _entropy(self, low, high):
        return log(high - low)


randint = randint_gen(name='randint', longname='A discrete uniform '
                      '(random integer)')


# FIXME: problems sampling.
class zipf_gen(rv_discrete):
    r"""A Zipf discrete random variable.

    %(before_notes)s

    Notes
    -----
    The probability mass function for `zipf` is:

    .. math::

        f(k, a) = \frac{1}{\zeta(a) k^a}

    for :math:`k \ge 1`.

    `zipf` takes :math:`a` as shape parameter. :math:`\zeta` is the
    Riemann zeta function (`scipy.special.zeta`)

    %(after_notes)s

    %(example)s

    """
    def _rvs(self, a):
        return self._random_state.zipf(a, size=self._size)

    def _argcheck(self, a):
        return a > 1

    def _pmf(self, k, a):
        # zipf.pmf(k, a) = 1/(zeta(a) * k**a)
        Pk = 1.0 / special.zeta(a, 1) / k**a
        return Pk

    def _munp(self, n, a):
        return _lazywhere(
            a > n + 1, (a, n),
            lambda a, n: special.zeta(a - n, 1) / special.zeta(a, 1),
            np.inf)


zipf = zipf_gen(a=1, name='zipf', longname='A Zipf')


class dlaplace_gen(rv_discrete):
    r"""A  Laplacian discrete random variable.

    %(before_notes)s

    Notes
    -----
    The probability mass function for `dlaplace` is:

    .. math::

        f(k) = \tanh(a/2) \exp(-a |k|)

    for integers :math:`k` and :math:`a > 0`.

    `dlaplace` takes :math:`a` as shape parameter.

    %(after_notes)s

    %(example)s

    """
    def _pmf(self, k, a):
        # dlaplace.pmf(k) = tanh(a/2) * exp(-a*abs(k))
        return tanh(a/2.0) * exp(-a * abs(k))

    def _cdf(self, x, a):
        k = floor(x)
        f = lambda k, a: 1.0 - exp(-a * k) / (exp(a) + 1)
        f2 = lambda k, a: exp(a * (k+1)) / (exp(a) + 1)
        return _lazywhere(k >= 0, (k, a), f=f, f2=f2)

    def _ppf(self, q, a):
        const = 1 + exp(a)
        vals = ceil(np.where(q < 1.0 / (1 + exp(-a)),
                             log(q*const) / a - 1,
                             -log((1-q) * const) / a))
        vals1 = vals - 1
        return np.where(self._cdf(vals1, a) >= q, vals1, vals)

    def _stats(self, a):
        ea = exp(a)
        mu2 = 2.*ea/(ea-1.)**2
        mu4 = 2.*ea*(ea**2+10.*ea+1.) / (ea-1.)**4
        return 0., mu2, 0., mu4/mu2**2 - 3.

    def _entropy(self, a):
        return a / sinh(a) - log(tanh(a/2.0))


dlaplace = dlaplace_gen(a=-np.inf,
                        name='dlaplace', longname='A discrete Laplacian')


class skellam_gen(rv_discrete):
    r"""A  Skellam discrete random variable.

    %(before_notes)s

    Notes
    -----
    Probability distribution of the difference of two correlated or
    uncorrelated Poisson random variables.

    Let :math:`k_1` and :math:`k_2` be two Poisson-distributed r.v. with
    expected values :math:`\lambda_1` and :math:`\lambda_2`. Then,
    :math:`k_1 - k_2` follows a Skellam distribution with parameters
    :math:`\mu_1 = \lambda_1 - \rho \sqrt{\lambda_1 \lambda_2}` and
    :math:`\mu_2 = \lambda_2 - \rho \sqrt{\lambda_1 \lambda_2}`, where
    :math:`\rho` is the correlation coefficient between :math:`k_1` and
    :math:`k_2`. If the two Poisson-distributed r.v. are independent then
    :math:`\rho = 0`.

    Parameters :math:`\mu_1` and :math:`\mu_2` must be strictly positive.

    For details see: https://en.wikipedia.org/wiki/Skellam_distribution

    `skellam` takes :math:`\mu_1` and :math:`\mu_2` as shape parameters.

    %(after_notes)s

    %(example)s

    """
    def _rvs(self, mu1, mu2):
        n = self._size
        return (self._random_state.poisson(mu1, n) -
                self._random_state.poisson(mu2, n))

    def _pmf(self, x, mu1, mu2):
        px = np.where(x < 0,
                      _ncx2_pdf(2*mu2, 2*(1-x), 2*mu1)*2,
                      _ncx2_pdf(2*mu1, 2*(1+x), 2*mu2)*2)
        # ncx2.pdf() returns nan's for extremely low probabilities
        return px

    def _cdf(self, x, mu1, mu2):
        x = floor(x)
        px = np.where(x < 0,
                      _ncx2_cdf(2*mu2, -2*x, 2*mu1),
                      1 - _ncx2_cdf(2*mu1, 2*(x+1), 2*mu2))
        return px

    def _stats(self, mu1, mu2):
        mean = mu1 - mu2
        var = mu1 + mu2
        g1 = mean / sqrt((var)**3)
        g2 = 1 / var
        return mean, var, g1, g2


skellam = skellam_gen(a=-np.inf, name="skellam", longname='A Skellam')


class yulesimon_gen(rv_discrete):
    r"""A Yule-Simon discrete random variable.

    %(before_notes)s

    Notes
    -----

    The probability mass function for the `yulesimon` is:

    .. math::

        f(k) =  \alpha B(k, \alpha+1)

    for :math:`k=1,2,3,...`, where :math:`\alpha>0`.
    Here :math:`B` refers to the `scipy.special.beta` function.

    The sampling of random variates is based on pg 553, Section 6.3 of [1]_.
    Our notation maps to the referenced logic via :math:`\alpha=a-1`.

    For details see the wikipedia entry [2]_.

    References
    ----------
    .. [1] Devroye, Luc. "Non-uniform Random Variate Generation",
         (1986) Springer, New York.

    .. [2] https://en.wikipedia.org/wiki/Yule-Simon_distribution

    %(after_notes)s

    %(example)s

    """
    def _rvs(self, alpha):
        E1 = self._random_state.standard_exponential(self._size)
        E2 = self._random_state.standard_exponential(self._size)
        ans = ceil(-E1 / log1p(-exp(-E2 / alpha)))
        return ans

    def _pmf(self, x, alpha):
        return alpha * special.beta(x, alpha + 1)

    def _argcheck(self, alpha):
        return (alpha > 0)

    def _logpmf(self, x, alpha):
        return log(alpha) + special.betaln(x, alpha + 1)

    def _cdf(self, x, alpha):
        return 1 - x * special.beta(x, alpha + 1)

    def _sf(self, x, alpha):
        return x * special.beta(x, alpha + 1)

    def _logsf(self, x, alpha):
        return log(x) + special.betaln(x, alpha + 1)

    def _stats(self, alpha):
        mu = np.where(alpha <= 1, np.inf, alpha / (alpha - 1))
        mu2 = np.where(alpha > 2,
                alpha**2 / ((alpha - 2.0) * (alpha - 1)**2),
                np.inf)
        mu2 = np.where(alpha <= 1, np.nan, mu2)
        g1 = np.where(alpha > 3,
                sqrt(alpha - 2) * (alpha + 1)**2 / (alpha * (alpha - 3)),
                np.inf)
        g1 = np.where(alpha <= 2, np.nan, g1)
        g2 = np.where(alpha > 4,
                (alpha + 3) + (alpha**3 - 49 * alpha - 22) / (alpha *
                        (alpha - 4) * (alpha - 3)), np.inf)
        g2 = np.where(alpha <= 2, np.nan, g2)
        return mu, mu2, g1, g2


yulesimon = yulesimon_gen(name='yulesimon', a=1)


# Collect names of classes and objects in this module.
pairs = list(globals().items())
_distn_names, _distn_gen_names = get_distribution_names(pairs, rv_discrete)

__all__ = _distn_names + _distn_gen_names
