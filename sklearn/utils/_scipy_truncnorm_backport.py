#
# Author:  Travis Oliphant  2002-2011 with contributions from
#          SciPy Developers 2004-2011
#
import numbers

import numpy as np
from .fixes import sp_version
if sp_version < (0, 14):
    from scipy.stats.distributions import rv_frozen, rv_continuous
    from scipy.stats.distributions import (_norm_cdf, _norm_sf, _norm_pdf,
                                           _norm_logpdf, _norm_isf, _norm_ppf)
else:
    from scipy.stats._distn_infrastructure import rv_frozen, rv_continuous
    from scipy.stats._continuous_distns import (_norm_cdf, _norm_sf, _norm_pdf,
                                                _norm_logpdf, _norm_isf,
                                                _norm_ppf)

__all__ = ['truncnorm']


def check_random_state(seed):
    """Turn seed into a np.random.RandomState instance
    Parameters
    ----------
    seed : None | int | instance of RandomState
        If seed is None, return the RandomState singleton used by np.random.
        If seed is an int, return a new RandomState instance seeded with seed.
        If seed is already a RandomState instance, return it.
        Otherwise raise ValueError.
    """
    if seed is None or seed is np.random:
        return np.random.mtrand._rand
    if isinstance(seed, (numbers.Integral, np.integer)):
        return np.random.RandomState(seed)
    if isinstance(seed, np.random.RandomState):
        return seed
    raise ValueError('%r cannot be used to seed a numpy.random.RandomState'
                     ' instance' % seed)


class rv_frozen_patch(rv_frozen):
    def rvs(self, size=None, random_state=None, **kwds):
        kwds = self.kwds.copy()
        kwds.update({'size': size, 'random_state': random_state})
        return self.dist.rvs(*self.args, **self.kwds)


class truncnorm_gen(rv_continuous):
    r"""A truncated normal continuous random variable.
    Notes
    -----
    The standard form of this distribution is a standard normal truncated to
    the range [a, b] --- notice that a and b are defined over the domain of the
    standard normal.  To convert clip values for a specific mean and standard
    deviation, use::
        a, b = (myclip_a - my_mean) / my_std, (myclip_b - my_mean) / my_std
    `truncnorm` takes :math:`a` and :math:`b` as shape parameters.
    """
    def __init__(self, seed=None, name=None, **kwds):
        super(truncnorm_gen, self).__init__(**kwds)
        self._random_state = check_random_state(seed)

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

    def _rvs(self, *args):
        # This method must handle self._size being a tuple, and it must
        # properly broadcast *args and self._size.  self._size might be
        # an empty tuple, which means a scalar random variate is to be
        # generated.

        # Use basic inverse cdf algorithm for RV generation as default.
        U = self._random_state.random_sample(self._size)
        Y = self._ppf(U, *args)
        return Y

    def rvs(self, *args, **kwds):
        """
        Random variates of given type.
        Parameters
        ----------
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution (see docstring of the
            instance object for more information).
        loc : array_like, optional
            Location parameter (default=0).
        scale : array_like, optional
            Scale parameter (default=1).
        size : int or tuple of ints, optional
            Defining number of random variates (default is 1).
        random_state : None or int or ``np.random.RandomState`` instance,
            optional.
            If int or RandomState, use it for drawing the random variates.
            If None, rely on ``self.random_state``.
            Default is None.
        Returns
        -------
        rvs : ndarray or scalar
            Random variates of given `size`.
        """
        discrete = kwds.pop('discrete', None)
        rndm = kwds.pop('random_state', None)
        args, loc, scale, size = self._parse_args_rvs(*args, **kwds)
        cond = np.logical_and(self._argcheck(*args), (scale >= 0))
        if not np.all(cond):
            raise ValueError("Domain error in arguments.")

        if np.all(scale == 0):
            return loc*np.ones(size, 'd')

        # extra gymnastics needed for a custom random_state
        if rndm is not None:
            from . import check_random_state
            random_state_saved = self._random_state
            self._random_state = check_random_state(rndm)

        # `size` should just be an argument to _rvs(), but for, um,
        # historical reasons, it is made an attribute that is read
        # by _rvs().
        self._size = size
        vals = self._rvs(*args)

        vals = vals * scale + loc

        # do not forget to restore the _random_state
        if rndm is not None:
            self._random_state = random_state_saved

        # Cast to int if discrete
        if discrete:
            if size == ():
                vals = int(vals)
            else:
                vals = vals.astype(int)

        return vals

    def freeze(self, *args, **kwds):
        """Freeze the distribution for the given arguments.
        Parameters
        ----------
        arg1, arg2, arg3,... : array_like
            The shape parameter(s) for the distribution.  Should include all
            the non-optional arguments, may include ``loc`` and ``scale``.
        Returns
        -------
        rv_frozen : rv_frozen instance
            The frozen distribution.
        """
        return rv_frozen_patch(self, *args, **kwds)


truncnorm = truncnorm_gen(name='truncnorm')
