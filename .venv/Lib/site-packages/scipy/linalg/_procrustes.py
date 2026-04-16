"""
Solve the orthogonal Procrustes problem.

"""

from scipy._lib._util import _apply_over_batch
from ._decomp_svd import svd
from scipy._lib._array_api import array_namespace, xp_capabilities, _asarray, is_numpy


__all__ = ['orthogonal_procrustes']


@xp_capabilities(
    jax_jit=False,
    skip_backends=[("dask.array", "full_matrices=True is not supported by dask")],
)
@_apply_over_batch(('A', 2), ('B', 2))
def orthogonal_procrustes(A, B, check_finite=True):
    """
    Compute the matrix solution of the orthogonal (or unitary) Procrustes problem.

    Given matrices `A` and `B` of the same shape, find an orthogonal (or unitary in
    the case of complex input) matrix `R` that most closely maps `A` to `B` using the
    algorithm given in [1]_.

    Parameters
    ----------
    A : (M, N) array_like
        Matrix to be mapped.
    B : (M, N) array_like
        Target matrix.
    check_finite : bool, optional
        Whether to check that the input matrices contain only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.

    Returns
    -------
    R : (N, N) ndarray
        The matrix solution of the orthogonal Procrustes problem.
        Minimizes the Frobenius norm of ``(A @ R) - B``, subject to
        ``R.conj().T @ R = I``.
    scale : float
        Sum of the singular values of ``A.conj().T @ B``.

    Raises
    ------
    ValueError
        If the input array shapes don't match or if check_finite is True and
        the arrays contain Inf or NaN.

    Notes
    -----
    Note that unlike higher level Procrustes analyses of spatial data, this
    function only uses orthogonal transformations like rotations and
    reflections, and it does not use scaling or translation.

    References
    ----------
    .. [1] Peter H. Schonemann, "A generalized solution of the orthogonal
           Procrustes problem", Psychometrica -- Vol. 31, No. 1, March, 1966.
           :doi:`10.1007/BF02289451`

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.linalg import orthogonal_procrustes
    >>> A = np.array([[ 2,  0,  1], [-2,  0,  0]])

    Flip the order of columns and check for the anti-diagonal mapping

    >>> R, sca = orthogonal_procrustes(A, np.fliplr(A))
    >>> R
    array([[-5.34384992e-17,  0.00000000e+00,  1.00000000e+00],
           [ 0.00000000e+00,  1.00000000e+00,  0.00000000e+00],
           [ 1.00000000e+00,  0.00000000e+00, -7.85941422e-17]])
    >>> sca
    9.0

    As an example of the unitary Procrustes problem, generate a
    random complex matrix ``A``, a random unitary matrix ``Q``,
    and their product ``B``.

    >>> shape = (4, 4)
    >>> rng = np.random.default_rng(589234981235)
    >>> A = rng.random(shape) + rng.random(shape)*1j
    >>> Q = rng.random(shape) + rng.random(shape)*1j
    >>> Q, _ = np.linalg.qr(Q)
    >>> B = A @ Q

    `orthogonal_procrustes` recovers the unitary matrix ``Q``
    from ``A`` and ``B``.

    >>> R, _ = orthogonal_procrustes(A, B)
    >>> np.allclose(R, Q)
    True

    """
    xp = array_namespace(A, B)

    A = _asarray(A, xp=xp, check_finite=check_finite, subok=True)
    B = _asarray(B, xp=xp, check_finite=check_finite, subok=True)

    if A.ndim != 2:
        raise ValueError(f'expected ndim to be 2, but observed {A.ndim}')
    if A.shape != B.shape:
        raise ValueError(f'the shapes of A and B differ ({A.shape} vs {B.shape})')

    # Be clever with transposes, with the intention to save memory.
    # The conjugate has no effect for real inputs, but gives the correct solution
    # for complex inputs.
    if is_numpy(xp):
        u, w, vt = svd((B.T @ xp.conj(A)).T)
    else:
        u, w, vt = xp.linalg.svd((B.T @ xp.conj(A)).T)
    R = u @ vt
    scale = xp.sum(w)
    return R, scale
