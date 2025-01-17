import numpy as np
from ._cython_nnls import _nnls


__all__ = ['nnls']


def nnls(A, b, maxiter=None, *, atol=None):
    """
    Solve ``argmin_x || Ax - b ||_2`` for ``x>=0``.

    This problem, often called as NonNegative Least Squares, is a convex
    optimization problem with convex constraints. It typically arises when
    the ``x`` models quantities for which only nonnegative values are
    attainable; weight of ingredients, component costs and so on.

    Parameters
    ----------
    A : (m, n) ndarray
        Coefficient array
    b : (m,) ndarray, float
        Right-hand side vector.
    maxiter: int, optional
        Maximum number of iterations, optional. Default value is ``3 * n``.
    atol: float
        Tolerance value used in the algorithm to assess closeness to zero in
        the projected residual ``(A.T @ (A x - b)`` entries. Increasing this
        value relaxes the solution constraints. A typical relaxation value can
        be selected as ``max(m, n) * np.linalg.norm(a, 1) * np.spacing(1.)``.
        This value is not set as default since the norm operation becomes
        expensive for large problems hence can be used only when necessary.

    Returns
    -------
    x : ndarray
        Solution vector.
    rnorm : float
        The 2-norm of the residual, ``|| Ax-b ||_2``.

    See Also
    --------
    lsq_linear : Linear least squares with bounds on the variables

    Notes
    -----
    The code is based on [2]_ which is an improved version of the classical
    algorithm of [1]_. It utilizes an active set method and solves the KKT
    (Karush-Kuhn-Tucker) conditions for the non-negative least squares problem.

    References
    ----------
    .. [1] : Lawson C., Hanson R.J., "Solving Least Squares Problems", SIAM,
       1995, :doi:`10.1137/1.9781611971217`
    .. [2] : Bro, Rasmus and de Jong, Sijmen, "A Fast Non-Negativity-
       Constrained Least Squares Algorithm", Journal Of Chemometrics, 1997,
       :doi:`10.1002/(SICI)1099-128X(199709/10)11:5<393::AID-CEM483>3.0.CO;2-L`

     Examples
    --------
    >>> import numpy as np
    >>> from scipy.optimize import nnls
    ...
    >>> A = np.array([[1, 0], [1, 0], [0, 1]])
    >>> b = np.array([2, 1, 1])
    >>> nnls(A, b)
    (array([1.5, 1. ]), 0.7071067811865475)

    >>> b = np.array([-1, -1, -1])
    >>> nnls(A, b)
    (array([0., 0.]), 1.7320508075688772)

    """

    A = np.asarray_chkfinite(A, dtype=np.float64, order='C')
    b = np.asarray_chkfinite(b, dtype=np.float64)

    if len(A.shape) != 2:
        raise ValueError("Expected a two-dimensional array (matrix)" +
                         f", but the shape of A is {A.shape}")
    if len(b.shape) != 1:
        raise ValueError("Expected a one-dimensional array (vector)" +
                         f", but the shape of b is {b.shape}")

    m, n = A.shape

    if m != b.shape[0]:
        raise ValueError(
                "Incompatible dimensions. The first dimension of " +
                f"A is {m}, while the shape of b is {(b.shape[0], )}")

    if not maxiter:
        maxiter = 3*n
    x, rnorm, info = _nnls(A, b, maxiter)
    if info == -1:
        raise RuntimeError("Maximum number of iterations reached.")

    return x, rnorm
