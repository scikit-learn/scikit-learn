""" Collection of Model instances for use with the odrpack fitting package.
"""
from __future__ import division, print_function, absolute_import

import numpy as np
from scipy.odr.odrpack import Model

__all__ = ['Model', 'exponential', 'multilinear', 'unilinear', 'quadratic',
           'polynomial']


def _lin_fcn(B, x):
    a, b = B[0], B[1:]
    b.shape = (b.shape[0], 1)

    return a + (x*b).sum(axis=0)


def _lin_fjb(B, x):
    a = np.ones(x.shape[-1], float)
    res = np.concatenate((a, x.ravel()))
    res.shape = (B.shape[-1], x.shape[-1])
    return res


def _lin_fjd(B, x):
    b = B[1:]
    b = np.repeat(b, (x.shape[-1],)*b.shape[-1],axis=0)
    b.shape = x.shape
    return b


def _lin_est(data):
    # Eh. The answer is analytical, so just return all ones.
    # Don't return zeros since that will interfere with
    # ODRPACK's auto-scaling procedures.

    if len(data.x.shape) == 2:
        m = data.x.shape[0]
    else:
        m = 1

    return np.ones((m + 1,), float)


def _poly_fcn(B, x, powers):
    a, b = B[0], B[1:]
    b.shape = (b.shape[0], 1)

    return a + np.sum(b * np.power(x, powers), axis=0)


def _poly_fjacb(B, x, powers):
    res = np.concatenate((np.ones(x.shape[-1], float), np.power(x,
        powers).flat))
    res.shape = (B.shape[-1], x.shape[-1])
    return res


def _poly_fjacd(B, x, powers):
    b = B[1:]
    b.shape = (b.shape[0], 1)

    b = b * powers

    return np.sum(b * np.power(x, powers-1),axis=0)


def _exp_fcn(B, x):
    return B[0] + np.exp(B[1] * x)


def _exp_fjd(B, x):
    return B[1] * np.exp(B[1] * x)


def _exp_fjb(B, x):
    res = np.concatenate((np.ones(x.shape[-1], float), x * np.exp(B[1] * x)))
    res.shape = (2, x.shape[-1])
    return res


def _exp_est(data):
    # Eh.
    return np.array([1., 1.])


multilinear = Model(_lin_fcn, fjacb=_lin_fjb,
               fjacd=_lin_fjd, estimate=_lin_est,
               meta={'name': 'Arbitrary-dimensional Linear',
                     'equ':'y = B_0 + Sum[i=1..m, B_i * x_i]',
                     'TeXequ':r'$y=\beta_0 + \sum_{i=1}^m \beta_i x_i$'})


def polynomial(order):
    """
    Factory function for a general polynomial model.

    Parameters
    ----------
    order : int or sequence
        If an integer, it becomes the order of the polynomial to fit. If
        a sequence of numbers, then these are the explicit powers in the
        polynomial.
        A constant term (power 0) is always included, so don't include 0.
        Thus, polynomial(n) is equivalent to polynomial(range(1, n+1)).

    Returns
    -------
    polynomial : Model instance
        Model instance.

    """

    powers = np.asarray(order)
    if powers.shape == ():
        # Scalar.
        powers = np.arange(1, powers + 1)

    powers.shape = (len(powers), 1)
    len_beta = len(powers) + 1

    def _poly_est(data, len_beta=len_beta):
        # Eh. Ignore data and return all ones.
        return np.ones((len_beta,), float)

    return Model(_poly_fcn, fjacd=_poly_fjacd, fjacb=_poly_fjacb,
                 estimate=_poly_est, extra_args=(powers,),
                 meta={'name': 'Sorta-general Polynomial',
                 'equ': 'y = B_0 + Sum[i=1..%s, B_i * (x**i)]' % (len_beta-1),
                 'TeXequ': r'$y=\beta_0 + \sum_{i=1}^{%s} \beta_i x^i$' %
                        (len_beta-1)})


exponential = Model(_exp_fcn, fjacd=_exp_fjd, fjacb=_exp_fjb,
                    estimate=_exp_est, meta={'name':'Exponential',
                    'equ': 'y= B_0 + exp(B_1 * x)',
                    'TeXequ': r'$y=\beta_0 + e^{\beta_1 x}$'})


def _unilin(B, x):
    return x*B[0] + B[1]


def _unilin_fjd(B, x):
    return np.ones(x.shape, float) * B[0]


def _unilin_fjb(B, x):
    _ret = np.concatenate((x, np.ones(x.shape, float)))
    _ret.shape = (2,) + x.shape

    return _ret


def _unilin_est(data):
    return (1., 1.)


def _quadratic(B, x):
    return x*(x*B[0] + B[1]) + B[2]


def _quad_fjd(B, x):
    return 2*x*B[0] + B[1]


def _quad_fjb(B, x):
    _ret = np.concatenate((x*x, x, np.ones(x.shape, float)))
    _ret.shape = (3,) + x.shape

    return _ret


def _quad_est(data):
    return (1.,1.,1.)


unilinear = Model(_unilin, fjacd=_unilin_fjd, fjacb=_unilin_fjb,
                  estimate=_unilin_est, meta={'name': 'Univariate Linear',
                  'equ': 'y = B_0 * x + B_1',
                  'TeXequ': '$y = \\beta_0 x + \\beta_1$'})

quadratic = Model(_quadratic, fjacd=_quad_fjd, fjacb=_quad_fjb,
                  estimate=_quad_est, meta={'name': 'Quadratic',
                  'equ': 'y = B_0*x**2 + B_1*x + B_2',
                  'TeXequ': '$y = \\beta_0 x^2 + \\beta_1 x + \\beta_2'})
