"""
  Restared Krylov method for evaluating f(A)b

  The original code was written in MATLAB by Stefan Guttel.
  It was later adapted to C++ and Python by Nicolas Guidotti.
  Both authors agrees to relicense the code under the BSD license.

  Copyright (C) 2025 Nicolas Guidotti and Stefan Guttel.
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  1. Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

  3. Neither the name of the copyright holder nor the names of its contributors
  may be used to endorse or promote products derived from this software without
  specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""

import numpy as np
from scipy.linalg import norm
from ._isolve.iterative import _get_atol_rtol

__all__ = ['funm_multiply_krylov']

def _funm_multiply_krylov_arnoldi(A, b, bnorm, V, H, m):
    """
    The Arnoldi iteration for constructing the basis V and the projection H = V * A V
    for the Krylov subspace Km(A, b) of order m.

    Parameters
    ----------
    A : transposable linear operator
        The operator whose matrix function is of interest.
    b : ndarray
        The vector b to multiply the f(A) with.
    V : ndarray
        The n x (m + 1) matrix whose columns determines the basis for
        Krylov subspace Km(A, b).
    H : ndarray
        A (m + 1) x m upper Hessenberg matrix representing the projection of A
        onto Km(A, b).
    m : int
        The order of the Krylov subspace.

    Returns
    -------
    breakdown : bool
        Indicate if the Arnoldi broke down or not

    iter : int
        Returns the last valid iteration.

    """

    dotprod = np.vdot if np.iscomplexobj(b) else np.dot
    norm_tol = np.finfo(b.dtype.char).eps ** 2
    V[:, 0] = b / bnorm

    for k in range(0, m):
        V[:, k + 1] = A.dot(V[:, k])

        # Uses the modified Gram-Schmift process to orthogonalize V[:, k + 1]
        # against the previous basis vectors
        for i in range(0, k + 1):
            H[i, k] = dotprod(V[:, i], V[:, k + 1])
            V[:, k + 1] = V[:, k + 1] - H[i, k] * V[:, i]

        H[k + 1, k] = norm(V[:, k + 1])
        if H[k + 1, k] < norm_tol:
            return True, k

        V[:, k + 1] = V[:, k + 1] / H[k + 1, k]

    return False, m

def _funm_multiply_krylov_lanczos(A, b, bnorm, V, H, m):
    """
    The Lanczos iteration for constructing the basis V and the projection H = V * A V
    for the Krylov subspace Km(A, b) of order m. A must be Hermitian.

    Parameters
    ----------
    A : transposable linear operator
        The operator whose matrix function is of interest.
    b : ndarray
        The vector b to multiply the f(A) with.
    V : ndarray
        The n x (m + 1) matrix whose columns determines the basis for
        Krylov subspace Km(A, b).
    H : ndarray
        A (m + 1) x m upper Hessenberg matrix representing the projection of A
        onto Km(A, b).
    m : int
        The order of the Krylov subspace.

    Returns
    -------
    breakdown : bool
        Indicate if the Arnoldi broke down or not

    iter : int
        Returns the last valid iteration.

    """
    dotprod = np.vdot if np.iscomplexobj(b) else np.dot
    norm_tol = np.finfo(b.dtype.char).eps ** 2
    V[:, 0] = b / bnorm

    for k in range(0, m):
        if k > 0:
            V[:, k + 1] = A.dot(V[:, k]) - H[k, k - 1] * V[:, k - 1]
        else:
            V[:, k + 1] = A.dot(V[:, k])

        H[k, k] = dotprod(V[:, k + 1], V[:, k])
        V[:, k + 1] = V[:, k + 1] - H[k, k] * V[:, k]

        H[k + 1, k] = norm(V[:, k + 1])

        if H[k + 1, k] < norm_tol:
            return True, k

        V[:, k + 1] = V[:, k + 1] / H[k + 1, k]
        if k < m - 1:
            H[k, k + 1] = H[k + 1, k]

    return False, m


def funm_multiply_krylov(f, A, b, *, assume_a = "general", t = 1.0, atol = 0.0,
                         rtol = 1e-6, restart_every_m = None, max_restarts = 20):
    """
    A restarted Krylov method for evaluating ``y = f(tA) b`` from [1]_ [2]_.

    Parameters
    ----------
    f : callable
        Callable object that computes the matrix function ``F = f(X)``.

    A : {sparse array, ndarray, LinearOperator}
        A real or complex N-by-N matrix.
        Alternatively, `A` can be a linear operator which can
        produce ``Ax`` using, e.g., ``scipy.sparse.linalg.LinearOperator``.

    b : ndarray
        A vector to multiply the ``f(tA)`` with.

    assume_a : string, optional
        Indicate the structure of ``A``. The algorithm will use this information
        to select the appropriated code path. The available options are
        'hermitian'/'her' and 'general'/'gen'. If ommited, then it is assumed
        that ``A`` has a 'general' structure.

    t : float, optional
        The value to scale the matrix ``A`` with. The default is ``t = 1.0``

    atol, rtol : float, optional
        Parameters for the convergence test. For convergence,
        ``norm(||y_k - y_k-1||) <= max(rtol*norm(b), atol)`` should be satisfied.
        The default is ``atol=0.`` and ``rtol=1e-6``.

    restart_every_m : integer
        If the iteration number reaches this value a restart is triggered.
        Larger values increase iteration cost but may be necessary for convergence.
        If omitted, ``min(20, n)`` is used.

    max_restarts : int, optional
        Maximum number of restart cycles. The algorithm will stop
        after max_restarts cycles even if the specified tolerance has not been
        achieved. The default is ``max_restarts=20``

    Returns
    -------
    y : ndarray
        The result of ``f(tA) b``.

    Notes
    -----
    The convergence of the Krylov method heavily depends on the spectrum
    of ``A`` and the function ``f``. With restarting, there are only formal
    proofs for functions of order 1 (e.g., ``exp``, ``sin``, ``cos``) and
    Stieltjes functions [2]_ [3]_, while the general case remains an open problem.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse import csr_array
    >>> from scipy.sparse.linalg import funm_multiply_krylov
    >>> from scipy.linalg import expm, solve
    >>> A = csr_array([[3, 2, 0], [1, -1, 0], [0, 5, 1]], dtype=float)
    >>> b = np.array([2, 4, -1], dtype=float)
    >>> t = 0.1

    Compute ``y = exp(tA) b``.

    >>> y = funm_multiply_krylov(expm, A, b, t = t)
    >>> y
    array([3.6164913 , 3.88421511 , 0.96073457])

    >>> ref = expm(t * A.todense()) @ b
    >>> err = y - ref
    >>> err
    array([4.44089210e-16 , 0.00000000e+00 , 2.22044605e-16])

    Compute :math:`y = (A^3 - A) b`.

    >>> poly = lambda X : X @ X @ X - X
    >>> y = funm_multiply_krylov(poly, A, b)
    >>> y
    array([132. , 24. , 70.])

    >>> ref = poly(A.todense()) @ b
    >>> err = y - ref
    >>> err
    array([ 0.00000000e+00 , 7.10542736e-15 , -2.84217094e-14])

    Compute :math:`y = f(tA) b`, where  :math:`f(X) = X^{-1}(e^{X} - I)`. This is
    known as the "phi function" from the exponential integrator literature.

    >>> phim_1 = lambda X : solve(X, expm(X) - np.eye(X.shape[0]))
    >>> y = funm_multiply_krylov(phim_1, A, b, t = t)
    >>> y
    array([ 2.76984306 , 3.92769192 , -0.03111392])

    >>> ref = phim_1(t * A.todense()) @ b
    >>> err = y - ref
    >>> err
    array([ 0.00000000e+00 , 8.88178420e-16 , -4.60742555e-15])

    References
    ----------
    .. [1] M. Afanasjew, M. Eiermann, O. G. Ernst, and S. Güttel,
          "Implementation of a restarted Krylov subspace method for the
          evaluation of matrix functions," Linear Algebra and its Applications,
          vol. 429, no. 10, pp. 2293-2314, Nov. 2008, :doi:`10.1016/j.laa.2008.06.029`.

    .. [2] M. Eiermann and O. G. Ernst, "A Restarted Krylov Subspace Method
           for the Evaluation of Matrix Functions," SIAM J. Numer. Anal., vol. 44,
           no. 6, pp. 2481-2504, Jan. 2006, :doi:`10.1137/050633846`.

    .. [3] A. Frommer, S. Güttel, and M. Schweitzer, "Convergence of Restarted
           Krylov Subspace Methods for Stieltjes Functions of Matrices," SIAM J.
           Matrix Anal. Appl., vol. 35, no. 4, pp. 1602-1624,
           Jan. 2014, :doi:`10.1137/140973463`.

    """

    if assume_a not in {'hermitian', 'general', 'her', 'gen'}:
        raise ValueError(f'scipy.sparse.linalg.funm_multiply_krylov: {assume_a} '
                         'is not a recognized matrix structure')
    is_hermitian = (assume_a == 'her') or (assume_a == 'hermitian')

    if len(b.shape) != 1:
        raise ValueError("scipy.sparse.linalg.funm_multiply_krylov: "
                         "argument 'b' must be a 1D array.")
    n = b.shape[0]

    if restart_every_m is None:
        restart_every_m = min(20, n)

    restart_every_m = int(restart_every_m)
    max_restarts = int(max_restarts)

    if restart_every_m <= 0:
        raise ValueError("scipy.sparse.linalg.funm_multiply_krylov: "
                         "argument 'restart_every_m' must be positive.")

    if max_restarts <= 0:
            raise ValueError("scipy.sparse.linalg.funm_multiply_krylov: "
                             "argument 'max_restarts' must be positive.")

    m = restart_every_m
    max_restarts = min(max_restarts, int(n / m) + 1)
    mmax = m * max_restarts

    bnorm = norm(b)
    atol, _ = _get_atol_rtol("funm_multiply_krylov", bnorm, atol, rtol)

    if bnorm == 0:
        y = np.array(b)
        return y

    # Preallocate the maximum memory space.
    # Using the column major order here since we work with
    # each individual column separately.
    internal_type = np.common_type(A, b)
    V = np.zeros((n, m + 1), dtype = internal_type, order = 'F')
    H = np.zeros((mmax + 1, mmax), dtype = internal_type, order = 'F')

    restart = 1

    if is_hermitian:
        breakdown, j = _funm_multiply_krylov_lanczos(A, b, bnorm, V,
                                                        H[:m + 1, :m], m)
    else:
        breakdown, j = _funm_multiply_krylov_arnoldi(A, b, bnorm, V,
                                                        H[:m + 1, :m], m)

    fH = f(t * H[:j, :j])
    y = bnorm * V[:, :j].dot(fH[:, 0])

    if breakdown:
        return y

    update_norm = norm(bnorm * fH[:, 0])

    while restart < max_restarts and update_norm > atol:
        begin = restart * m
        end = (restart + 1) * m

        if is_hermitian:
            breakdown, j = _funm_multiply_krylov_lanczos(A, V[:, m], 1, V,
                                                         H[begin:end + 1, begin:end], m)
        else:
            breakdown, j = _funm_multiply_krylov_arnoldi(A, V[:, m], 1, V,
                                                         H[begin:end + 1, begin:end], m)

        if breakdown:
            end = begin + j
            fH = f(t * H[:end, :end])
            y[:end] = y[:end] + bnorm * V[:, :m].dot(fH[begin:end, 0])
            return y

        fH = f(t * H[:end, :end])
        y = y + bnorm * V[:, :m].dot(fH[begin:end, 0])
        update_norm = norm(bnorm * fH[begin:end, 0])
        restart += 1

    return y
