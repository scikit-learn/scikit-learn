"""
:author: Gary Ruben, 2009
:license: modified BSD
"""

__all__ = ["frt2", "ifrt2"]

import numpy as np
from numpy import roll, newaxis


def frt2(a):
    """Compute the 2-dimensional finite Radon transform (FRT) for the input array.

    Parameters
    ----------
    a : ndarray of int, shape (M, M)
        Input array.

    Returns
    -------
    FRT : ndarray of int, shape (M+1, M)
        Finite Radon Transform array of coefficients.

    See Also
    --------
    ifrt2 : The two-dimensional inverse FRT.

    Notes
    -----
    The FRT has a unique inverse if and only if M is prime. [FRT]
    The idea for this algorithm is due to Vlad Negnevitski.

    Examples
    --------

    Generate a test image:
    Use a prime number for the array dimensions

    >>> SIZE = 59
    >>> img = np.tri(SIZE, dtype=np.int32)

    Apply the Finite Radon Transform:

    >>> f = frt2(img)

    References
    ----------
    .. [FRT] A. Kingston and I. Svalbe, "Projective transforms on periodic
             discrete image arrays," in P. Hawkes (Ed), Advances in Imaging
             and Electron Physics, 139 (2006)

    """
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError("Input must be a square, 2-D array")

    ai = a.copy()
    n = ai.shape[0]
    f = np.empty((n + 1, n), np.uint32)
    f[0] = ai.sum(axis=0)
    for m in range(1, n):
        # Roll the pth row of ai left by p places
        for row in range(1, n):
            ai[row] = roll(ai[row], -row)
        f[m] = ai.sum(axis=0)
    f[n] = ai.sum(axis=1)
    return f


def ifrt2(a):
    """Compute the 2-dimensional inverse finite Radon transform (iFRT) for the input array.

    Parameters
    ----------
    a : ndarray of int, shape (M+1, M)
        Input array.

    Returns
    -------
    iFRT : ndarray of int, shape (M, M)
        Inverse Finite Radon Transform coefficients.

    See Also
    --------
    frt2 : The two-dimensional FRT

    Notes
    -----
    The FRT has a unique inverse if and only if M is prime.
    See [1]_ for an overview.
    The idea for this algorithm is due to Vlad Negnevitski.

    Examples
    --------

    >>> SIZE = 59
    >>> img = np.tri(SIZE, dtype=np.int32)

    Apply the Finite Radon Transform:

    >>> f = frt2(img)

    Apply the Inverse Finite Radon Transform to recover the input

    >>> fi = ifrt2(f)

    Check that it's identical to the original

    >>> assert len(np.nonzero(img-fi)[0]) == 0

    References
    ----------
    .. [1] A. Kingston and I. Svalbe, "Projective transforms on periodic
             discrete image arrays," in P. Hawkes (Ed), Advances in Imaging
             and Electron Physics, 139 (2006)

    """
    if a.ndim != 2 or a.shape[0] != a.shape[1] + 1:
        raise ValueError("Input must be an (n+1) row x n column, 2-D array")

    ai = a.copy()[:-1]
    n = ai.shape[1]
    f = np.empty((n, n), np.uint32)
    f[0] = ai.sum(axis=0)
    for m in range(1, n):
        # Rolls the pth row of ai right by p places.
        for row in range(1, ai.shape[0]):
            ai[row] = roll(ai[row], row)
        f[m] = ai.sum(axis=0)
    f += a[-1][newaxis].T
    f = (f - ai[0].sum()) / n
    return f
