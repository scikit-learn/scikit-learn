"""Delegation layer for linear algebra functions."""

from . import _agnostic
from ._lib import _compat, _helpers
from ._lib._typing import Array, ArrayNamespace

__all__ = ["kron"]


def kron(
    a: Array | complex,
    b: Array | complex,
    /,
    *,
    xp: ArrayNamespace | None = None,
) -> Array:
    """
    Kronecker product of two arrays.

    Computes the Kronecker product, a composite array made of blocks of the
    second array scaled by the first.

    Equivalent to ``numpy.kron`` for NumPy arrays.

    Parameters
    ----------
    a, b : Array | int | float | complex
        Input arrays or scalars. At least one must be an array.
    xp : array_namespace, optional
        The standard-compatible namespace for `a` and `b`. Default: infer.

    Returns
    -------
    array
        The Kronecker product of `a` and `b`.

    Notes
    -----
    The function assumes that the number of dimensions of `a` and `b`
    are the same, if necessary prepending the smallest with ones.
    If ``a.shape = (r0,r1,..,rN)`` and ``b.shape = (s0,s1,...,sN)``,
    the Kronecker product has shape ``(r0*s0, r1*s1, ..., rN*SN)``.
    The elements are products of elements from `a` and `b`, organized
    explicitly by::

        kron(a,b)[k0,k1,...,kN] = a[i0,i1,...,iN] * b[j0,j1,...,jN]

    where::

        kt = it * st + jt,  t = 0,...,N

    In the common 2-D case (N=1), the block structure can be visualized::

        [[ a[0,0]*b,   a[0,1]*b,  ... , a[0,-1]*b  ],
         [  ...                              ...   ],
         [ a[-1,0]*b,  a[-1,1]*b, ... , a[-1,-1]*b ]]

    Examples
    --------
    >>> import array_api_strict as xp
    >>> import array_api_extra as xpx
    >>> xpx.kron(xp.asarray([1, 10, 100]), xp.asarray([5, 6, 7]), xp=xp)
    Array([  5,   6,   7,  50,  60,  70, 500,
           600, 700], dtype=array_api_strict.int64)

    >>> xpx.kron(xp.asarray([5, 6, 7]), xp.asarray([1, 10, 100]), xp=xp)
    Array([  5,  50, 500,   6,  60, 600,   7,
            70, 700], dtype=array_api_strict.int64)

    >>> xpx.kron(xp.eye(2), xp.ones((2, 2)), xp=xp)
    Array([[1., 1., 0., 0.],
           [1., 1., 0., 0.],
           [0., 0., 1., 1.],
           [0., 0., 1., 1.]], dtype=array_api_strict.float64)

    >>> a = xp.reshape(xp.arange(100), (2, 5, 2, 5))
    >>> b = xp.reshape(xp.arange(24), (2, 3, 4))
    >>> c = xpx.kron(a, b, xp=xp)
    >>> c.shape
    (2, 10, 6, 20)
    >>> I = (1, 3, 0, 2)
    >>> J = (0, 2, 1)
    >>> J1 = (0,) + J             # extend to ndim=4
    >>> S1 = (1,) + b.shape
    >>> K = tuple(xp.asarray(I) * xp.asarray(S1) + xp.asarray(J1))
    >>> c[K] == a[I]*b[J]
    Array(True, dtype=array_api_strict.bool)
    """
    if xp is None:
        xp = _compat.array_namespace(a, b)

    a, b = _helpers.asarrays(a, b, xp=xp)

    if (
        _compat.is_cupy_namespace(xp)
        or _compat.is_jax_namespace(xp)
        or _compat.is_numpy_namespace(xp)
        or _compat.is_torch_namespace(xp)
    ):
        return xp.kron(a, b)

    return _agnostic._linalg.kron(a, b, xp=xp)
