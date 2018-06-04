from __future__ import division, print_function, absolute_import

from numpy import zeros, asarray, eye, poly1d, hstack, r_
from scipy import linalg

__all__ = ["pade"]

def pade(an, m):
    """
    Return Pade approximation to a polynomial as the ratio of two polynomials.

    Parameters
    ----------
    an : (N,) array_like
        Taylor series coefficients.
    m : int
        The order of the returned approximating polynomials.

    Returns
    -------
    p, q : Polynomial class
        The Pade approximation of the polynomial defined by `an` is
        ``p(x)/q(x)``.

    Examples
    --------
    >>> from scipy.interpolate import pade
    >>> e_exp = [1.0, 1.0, 1.0/2.0, 1.0/6.0, 1.0/24.0, 1.0/120.0]
    >>> p, q = pade(e_exp, 2)

    >>> e_exp.reverse()
    >>> e_poly = np.poly1d(e_exp)

    Compare ``e_poly(x)`` and the Pade approximation ``p(x)/q(x)``

    >>> e_poly(1)
    2.7166666666666668

    >>> p(1)/q(1)
    2.7179487179487181

    """
    an = asarray(an)
    N = len(an) - 1
    n = N - m
    if n < 0:
        raise ValueError("Order of q <m> must be smaller than len(an)-1.")
    Akj = eye(N+1, n+1)
    Bkj = zeros((N+1, m), 'd')
    for row in range(1, m+1):
        Bkj[row,:row] = -(an[:row])[::-1]
    for row in range(m+1, N+1):
        Bkj[row,:] = -(an[row-m:row])[::-1]
    C = hstack((Akj, Bkj))
    pq = linalg.solve(C, an)
    p = pq[:n+1]
    q = r_[1.0, pq[n+1:]]
    return poly1d(p[::-1]), poly1d(q[::-1])

