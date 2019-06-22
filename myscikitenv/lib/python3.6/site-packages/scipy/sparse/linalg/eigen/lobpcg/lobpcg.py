"""
Locally Optimal Block Preconditioned Conjugate Gradient Method (LOBPCG).

References
----------
.. [1] A. V. Knyazev (2001),
       Toward the Optimal Preconditioned Eigensolver: Locally Optimal
       Block Preconditioned Conjugate Gradient Method.
       SIAM Journal on Scientific Computing 23, no. 2,
       pp. 517-541. http://dx.doi.org/10.1137/S1064827500366124

.. [2] A. V. Knyazev, I. Lashuk, M. E. Argentati, and E. Ovchinnikov (2007),
       Block Locally Optimal Preconditioned Eigenvalue Xolvers (BLOPEX)
       in hypre and PETSc.  https://arxiv.org/abs/0705.2626

.. [3] A. V. Knyazev's C and MATLAB implementations:
       https://bitbucket.org/joseroman/blopex
"""

from __future__ import division, print_function, absolute_import
import numpy as np
from scipy.linalg import (inv, eigh, cho_factor, cho_solve, cholesky,
                          LinAlgError)
from scipy.sparse.linalg import aslinearoperator
from scipy.sparse.sputils import bmat

__all__ = ['lobpcg']


def _save(ar, fileName):
    # Used only when verbosity level > 10.
    np.savetxt(fileName, ar)


def _report_nonhermitian(M, a, b, name):
    """
    Report if `M` is not a hermitian matrix given the tolerances `a`, `b`.
    """
    from scipy.linalg import norm

    md = M - M.T.conj()

    nmd = norm(md, 1)
    tol = np.spacing(max(10**a, (10**b)*norm(M, 1)))
    if nmd > tol:
        print('matrix %s is not sufficiently Hermitian for a=%d, b=%d:'
              % (name, a, b))
        print('condition: %.e < %e' % (nmd, tol))


def _as2d(ar):
    """
    If the input array is 2D return it, if it is 1D, append a dimension,
    making it a column vector.
    """
    if ar.ndim == 2:
        return ar
    else:  # Assume 1!
        aux = np.array(ar, copy=False)
        aux.shape = (ar.shape[0], 1)
        return aux


def _makeOperator(operatorInput, expectedShape):
    """Takes a dense numpy array or a sparse matrix or
    a function and makes an operator performing matrix * blockvector
    products."""
    if operatorInput is None:
        return None
    else:
        operator = aslinearoperator(operatorInput)

    if operator.shape != expectedShape:
        raise ValueError('operator has invalid shape')

    return operator


def _applyConstraints(blockVectorV, factYBY, blockVectorBY, blockVectorY):
    """Changes blockVectorV in place."""
    gramYBV = np.dot(blockVectorBY.T.conj(), blockVectorV)
    tmp = cho_solve(factYBY, gramYBV)
    blockVectorV -= np.dot(blockVectorY, tmp)


def _b_orthonormalize(B, blockVectorV, blockVectorBV=None, retInvR=False):
    if blockVectorBV is None:
        if B is not None:
            blockVectorBV = B(blockVectorV)
        else:
            blockVectorBV = blockVectorV  # Shared data!!!
    gramVBV = np.dot(blockVectorV.T.conj(), blockVectorBV)
    gramVBV = cholesky(gramVBV)
    gramVBV = inv(gramVBV, overwrite_a=True)
    # gramVBV is now R^{-1}.
    blockVectorV = np.dot(blockVectorV, gramVBV)
    if B is not None:
        blockVectorBV = np.dot(blockVectorBV, gramVBV)
    else:
        blockVectorBV = None

    if retInvR:
        return blockVectorV, blockVectorBV, gramVBV
    else:
        return blockVectorV, blockVectorBV


def _get_indx(_lambda, num, largest):
    """Get `num` indices into `_lambda` depending on `largest` option."""
    ii = np.argsort(_lambda)
    if largest:
        ii = ii[:-num-1:-1]
    else:
        ii = ii[:num]

    return ii


def lobpcg(A, X,
           B=None, M=None, Y=None,
           tol=None, maxiter=20,
           largest=True, verbosityLevel=0,
           retLambdaHistory=False, retResidualNormsHistory=False):
    """Locally Optimal Block Preconditioned Conjugate Gradient Method (LOBPCG)

    LOBPCG is a preconditioned eigensolver for large symmetric positive
    definite (SPD) generalized eigenproblems.

    Parameters
    ----------
    A : {sparse matrix, dense matrix, LinearOperator}
        The symmetric linear operator of the problem, usually a
        sparse matrix.  Often called the "stiffness matrix".
    X : array_like
        Initial approximation to the k eigenvectors. If A has
        shape=(n,n) then X should have shape shape=(n,k).
    B : {dense matrix, sparse matrix, LinearOperator}, optional
        the right hand side operator in a generalized eigenproblem.
        by default, B = Identity
        often called the "mass matrix"
    M : {dense matrix, sparse matrix, LinearOperator}, optional
        preconditioner to A; by default M = Identity
        M should approximate the inverse of A
    Y : array_like, optional
        n-by-sizeY matrix of constraints, sizeY < n
        The iterations will be performed in the B-orthogonal complement
        of the column-space of Y. Y must be full rank.
    tol : scalar, optional
        Solver tolerance (stopping criterion)
        by default: tol=n*sqrt(eps)
    maxiter : integer, optional
        maximum number of iterations
        by default: maxiter=min(n,20)
    largest : bool, optional
        when True, solve for the largest eigenvalues, otherwise the smallest
    verbosityLevel : integer, optional
        controls solver output.  default: verbosityLevel = 0.
    retLambdaHistory : boolean, optional
        whether to return eigenvalue history
    retResidualNormsHistory : boolean, optional
        whether to return history of residual norms

    Returns
    -------
    w : array
        Array of k eigenvalues
    v : array
        An array of k eigenvectors.  V has the same shape as X.
    lambdas : list of arrays, optional
        The eigenvalue history, if `retLambdaHistory` is True.
    rnorms : list of arrays, optional
        The history of residual norms, if `retResidualNormsHistory` is True.

    Examples
    --------

    Solve A x = lambda B x with constraints and preconditioning.

    >>> from scipy.sparse import spdiags, issparse
    >>> from scipy.sparse.linalg import lobpcg, LinearOperator
    >>> n = 100
    >>> vals = [np.arange(n, dtype=np.float64) + 1]
    >>> A = spdiags(vals, 0, n, n)
    >>> A.toarray()
    array([[  1.,   0.,   0., ...,   0.,   0.,   0.],
           [  0.,   2.,   0., ...,   0.,   0.,   0.],
           [  0.,   0.,   3., ...,   0.,   0.,   0.],
           ...,
           [  0.,   0.,   0., ...,  98.,   0.,   0.],
           [  0.,   0.,   0., ...,   0.,  99.,   0.],
           [  0.,   0.,   0., ...,   0.,   0., 100.]])

    Constraints.

    >>> Y = np.eye(n, 3)

    Initial guess for eigenvectors, should have linearly independent
    columns. Column dimension = number of requested eigenvalues.

    >>> X = np.random.rand(n, 3)

    Preconditioner -- inverse of A (as an abstract linear operator).

    >>> invA = spdiags([1./vals[0]], 0, n, n)
    >>> def precond( x ):
    ...     return invA  * x
    >>> M = LinearOperator(matvec=precond, shape=(n, n), dtype=float)

    Here, ``invA`` could of course have been used directly as a preconditioner.
    Let us then solve the problem:

    >>> eigs, vecs = lobpcg(A, X, Y=Y, M=M, largest=False)
    >>> eigs
    array([4., 5., 6.])

    Note that the vectors passed in Y are the eigenvectors of the 3 smallest
    eigenvalues. The results returned are orthogonal to those.

    Notes
    -----
    If both retLambdaHistory and retResidualNormsHistory are True,
    the return tuple has the following format
    (lambda, V, lambda history, residual norms history).

    In the following ``n`` denotes the matrix size and ``m`` the number
    of required eigenvalues (smallest or largest).

    The LOBPCG code internally solves eigenproblems of the size 3``m`` on every
    iteration by calling the "standard" dense eigensolver, so if ``m`` is not
    small enough compared to ``n``, it does not make sense to call the LOBPCG
    code, but rather one should use the "standard" eigensolver,
    e.g. numpy or scipy function in this case.
    If one calls the LOBPCG algorithm for 5``m``>``n``,
    it will most likely break internally, so the code tries to call
    the standard function instead.

    It is not that n should be large for the LOBPCG to work, but rather the
    ratio ``n``/``m`` should be large. It you call LOBPCG with ``m``=1
    and ``n``=10, it works though ``n`` is small. The method is intended
    for extremely large ``n``/``m``, see e.g., reference [28] in
    https://arxiv.org/abs/0705.2626

    The convergence speed depends basically on two factors:

    1. How well relatively separated the seeking eigenvalues are from the rest
       of the eigenvalues. One can try to vary ``m`` to make this better.

    2. How well conditioned the problem is. This can be changed by using proper
       preconditioning. For example, a rod vibration test problem (under tests
       directory) is ill-conditioned for large ``n``, so convergence will be
       slow, unless efficient preconditioning is used. For this specific
       problem, a good simple preconditioner function would be a linear solve
       for A, which is easy to code since A is tridiagonal.

    *Acknowledgements*

    lobpcg.py code was written by Robert Cimrman.
    Many thanks belong to Andrew Knyazev, the author of the algorithm,
    for lots of advice and support.

    References
    ----------
    .. [1] A. V. Knyazev (2001),
           Toward the Optimal Preconditioned Eigensolver: Locally Optimal
           Block Preconditioned Conjugate Gradient Method.
           SIAM Journal on Scientific Computing 23, no. 2,
           pp. 517-541. http://dx.doi.org/10.1137/S1064827500366124

    .. [2] A. V. Knyazev, I. Lashuk, M. E. Argentati, and E. Ovchinnikov
           (2007), Block Locally Optimal Preconditioned Eigenvalue Xolvers
           (BLOPEX) in hypre and PETSc. https://arxiv.org/abs/0705.2626

    .. [3] A. V. Knyazev's C and MATLAB implementations:
           https://bitbucket.org/joseroman/blopex
    """
    blockVectorX = X
    blockVectorY = Y
    residualTolerance = tol
    maxIterations = maxiter

    if blockVectorY is not None:
        sizeY = blockVectorY.shape[1]
    else:
        sizeY = 0

    # Block size.
    if len(blockVectorX.shape) != 2:
        raise ValueError('expected rank-2 array for argument X')

    n, sizeX = blockVectorX.shape

    if verbosityLevel:
        aux = "Solving "
        if B is None:
            aux += "standard"
        else:
            aux += "generalized"
        aux += " eigenvalue problem with"
        if M is None:
            aux += "out"
        aux += " preconditioning\n\n"
        aux += "matrix size %d\n" % n
        aux += "block size %d\n\n" % sizeX
        if blockVectorY is None:
            aux += "No constraints\n\n"
        else:
            if sizeY > 1:
                aux += "%d constraints\n\n" % sizeY
            else:
                aux += "%d constraint\n\n" % sizeY
        print(aux)

    A = _makeOperator(A, (n, n))
    B = _makeOperator(B, (n, n))
    M = _makeOperator(M, (n, n))

    if (n - sizeY) < (5 * sizeX):
        # warn('The problem size is small compared to the block size.' \
        #        ' Using dense eigensolver instead of LOBPCG.')

        sizeX = min(sizeX, n)

        if blockVectorY is not None:
            raise NotImplementedError('The dense eigensolver '
                                      'does not support constraints.')

        # Define the closed range of indices of eigenvalues to return.
        if largest:
            eigvals = (n - sizeX, n-1)
        else:
            eigvals = (0, sizeX-1)

        A_dense = A(np.eye(n, dtype=A.dtype))
        B_dense = None if B is None else B(np.eye(n, dtype=B.dtype))

        vals, vecs = eigh(A_dense, B_dense, eigvals=eigvals,
                          check_finite=False)
        if largest:
            # Reverse order to be compatible with eigs() in 'LM' mode.
            vals = vals[::-1]
            vecs = vecs[:, ::-1]

        return vals, vecs

    if (residualTolerance is None) or (residualTolerance <= 0.0):
        residualTolerance = np.sqrt(1e-15) * n

    # Apply constraints to X.
    if blockVectorY is not None:

        if B is not None:
            blockVectorBY = B(blockVectorY)
        else:
            blockVectorBY = blockVectorY

        # gramYBY is a dense array.
        gramYBY = np.dot(blockVectorY.T.conj(), blockVectorBY)
        try:
            # gramYBY is a Cholesky factor from now on...
            gramYBY = cho_factor(gramYBY)
        except LinAlgError:
            raise ValueError('cannot handle linearly dependent constraints')

        _applyConstraints(blockVectorX, gramYBY, blockVectorBY, blockVectorY)

    ##
    # B-orthonormalize X.
    blockVectorX, blockVectorBX = _b_orthonormalize(B, blockVectorX)

    ##
    # Compute the initial Ritz vectors: solve the eigenproblem.
    blockVectorAX = A(blockVectorX)
    gramXAX = np.dot(blockVectorX.T.conj(), blockVectorAX)

    _lambda, eigBlockVector = eigh(gramXAX, check_finite=False)
    ii = _get_indx(_lambda, sizeX, largest)
    _lambda = _lambda[ii]

    eigBlockVector = np.asarray(eigBlockVector[:, ii])
    blockVectorX = np.dot(blockVectorX, eigBlockVector)
    blockVectorAX = np.dot(blockVectorAX, eigBlockVector)
    if B is not None:
        blockVectorBX = np.dot(blockVectorBX, eigBlockVector)

    ##
    # Active index set.
    activeMask = np.ones((sizeX,), dtype=bool)

    lambdaHistory = [_lambda]
    residualNormsHistory = []

    previousBlockSize = sizeX
    ident = np.eye(sizeX, dtype=A.dtype)
    ident0 = np.eye(sizeX, dtype=A.dtype)

    ##
    # Main iteration loop.

    blockVectorP = None  # set during iteration
    blockVectorAP = None
    blockVectorBP = None

    iterationNumber = -1
    while iterationNumber < maxIterations:
        iterationNumber += 1
        if verbosityLevel > 0:
            print('iteration %d' % iterationNumber)

        if B is not None:
            aux = blockVectorBX * _lambda[np.newaxis, :]

        else:
            aux = blockVectorX * _lambda[np.newaxis, :]

        blockVectorR = blockVectorAX - aux

        aux = np.sum(blockVectorR.conjugate() * blockVectorR, 0)
        residualNorms = np.sqrt(aux)

        residualNormsHistory.append(residualNorms)

        ii = np.where(residualNorms > residualTolerance, True, False)
        activeMask = activeMask & ii
        if verbosityLevel > 2:
            print(activeMask)

        currentBlockSize = activeMask.sum()
        if currentBlockSize != previousBlockSize:
            previousBlockSize = currentBlockSize
            ident = np.eye(currentBlockSize, dtype=A.dtype)

        if currentBlockSize == 0:
            break

        if verbosityLevel > 0:
            print('current block size:', currentBlockSize)
            print('eigenvalue:', _lambda)
            print('residual norms:', residualNorms)
        if verbosityLevel > 10:
            print(eigBlockVector)

        activeBlockVectorR = _as2d(blockVectorR[:, activeMask])

        if iterationNumber > 0:
            activeBlockVectorP = _as2d(blockVectorP[:, activeMask])
            activeBlockVectorAP = _as2d(blockVectorAP[:, activeMask])
            if B is not None:
                activeBlockVectorBP = _as2d(blockVectorBP[:, activeMask])

        if M is not None:
            # Apply preconditioner T to the active residuals.
            activeBlockVectorR = M(activeBlockVectorR)

        ##
        # Apply constraints to the preconditioned residuals.
        if blockVectorY is not None:
            _applyConstraints(activeBlockVectorR,
                              gramYBY, blockVectorBY, blockVectorY)

        ##
        # B-orthonormalize the preconditioned residuals.

        aux = _b_orthonormalize(B, activeBlockVectorR)
        activeBlockVectorR, activeBlockVectorBR = aux

        activeBlockVectorAR = A(activeBlockVectorR)

        if iterationNumber > 0:
            if B is not None:
                aux = _b_orthonormalize(B, activeBlockVectorP,
                                        activeBlockVectorBP, retInvR=True)
                activeBlockVectorP, activeBlockVectorBP, invR = aux
                activeBlockVectorAP = np.dot(activeBlockVectorAP, invR)

            else:
                aux = _b_orthonormalize(B, activeBlockVectorP, retInvR=True)
                activeBlockVectorP, _, invR = aux
                activeBlockVectorAP = np.dot(activeBlockVectorAP, invR)

        ##
        # Perform the Rayleigh Ritz Procedure:
        # Compute symmetric Gram matrices:

        if B is not None:
            xaw = np.dot(blockVectorX.T.conj(), activeBlockVectorAR)
            waw = np.dot(activeBlockVectorR.T.conj(), activeBlockVectorAR)
            xbw = np.dot(blockVectorX.T.conj(), activeBlockVectorBR)

            if iterationNumber > 0:
                xap = np.dot(blockVectorX.T.conj(), activeBlockVectorAP)
                wap = np.dot(activeBlockVectorR.T.conj(), activeBlockVectorAP)
                pap = np.dot(activeBlockVectorP.T.conj(), activeBlockVectorAP)
                xbp = np.dot(blockVectorX.T.conj(), activeBlockVectorBP)
                wbp = np.dot(activeBlockVectorR.T.conj(), activeBlockVectorBP)

                gramA = bmat([[np.diag(_lambda), xaw, xap],
                              [xaw.T.conj(), waw, wap],
                              [xap.T.conj(), wap.T.conj(), pap]])

                gramB = bmat([[ident0, xbw, xbp],
                              [xbw.T.conj(), ident, wbp],
                              [xbp.T.conj(), wbp.T.conj(), ident]])
            else:
                gramA = bmat([[np.diag(_lambda), xaw],
                              [xaw.T.conj(), waw]])
                gramB = bmat([[ident0, xbw],
                              [xbw.T.conj(), ident]])

        else:
            xaw = np.dot(blockVectorX.T.conj(), activeBlockVectorAR)
            waw = np.dot(activeBlockVectorR.T.conj(), activeBlockVectorAR)
            xbw = np.dot(blockVectorX.T.conj(), activeBlockVectorR)

            if iterationNumber > 0:
                xap = np.dot(blockVectorX.T.conj(), activeBlockVectorAP)
                wap = np.dot(activeBlockVectorR.T.conj(), activeBlockVectorAP)
                pap = np.dot(activeBlockVectorP.T.conj(), activeBlockVectorAP)
                xbp = np.dot(blockVectorX.T.conj(), activeBlockVectorP)
                wbp = np.dot(activeBlockVectorR.T.conj(), activeBlockVectorP)

                gramA = bmat([[np.diag(_lambda), xaw, xap],
                              [xaw.T.conj(), waw, wap],
                              [xap.T.conj(), wap.T.conj(), pap]])

                gramB = bmat([[ident0, xbw, xbp],
                              [xbw.T.conj(), ident, wbp],
                              [xbp.T.conj(), wbp.T.conj(), ident]])
            else:
                gramA = bmat([[np.diag(_lambda), xaw],
                              [xaw.T.conj(), waw]])
                gramB = bmat([[ident0, xbw],
                              [xbw.T.conj(), ident]])

        if verbosityLevel > 0:
            _report_nonhermitian(gramA, 3, -1, 'gramA')
            _report_nonhermitian(gramB, 3, -1, 'gramB')

        if verbosityLevel > 10:
            _save(gramA, 'gramA')
            _save(gramB, 'gramB')

        # Solve the generalized eigenvalue problem.
        _lambda, eigBlockVector = eigh(gramA, gramB, check_finite=False)
        ii = _get_indx(_lambda, sizeX, largest)

        if verbosityLevel > 10:
            print(ii)
            print(_lambda)

        _lambda = _lambda[ii]
        eigBlockVector = eigBlockVector[:, ii]

        lambdaHistory.append(_lambda)

        if verbosityLevel > 10:
            print('lambda:', _lambda)
#         # Normalize eigenvectors!
#         aux = np.sum( eigBlockVector.conjugate() * eigBlockVector, 0 )
#         eigVecNorms = np.sqrt( aux )
#         eigBlockVector = eigBlockVector / eigVecNorms[np.newaxis, :]
#         eigBlockVector, aux = _b_orthonormalize( B, eigBlockVector )

        if verbosityLevel > 10:
            print(eigBlockVector)

        # Compute Ritz vectors.
        if B is not None:
            if iterationNumber > 0:
                eigBlockVectorX = eigBlockVector[:sizeX]
                eigBlockVectorR = eigBlockVector[sizeX:sizeX+currentBlockSize]
                eigBlockVectorP = eigBlockVector[sizeX+currentBlockSize:]

                pp = np.dot(activeBlockVectorR, eigBlockVectorR)
                pp += np.dot(activeBlockVectorP, eigBlockVectorP)

                app = np.dot(activeBlockVectorAR, eigBlockVectorR)
                app += np.dot(activeBlockVectorAP, eigBlockVectorP)

                bpp = np.dot(activeBlockVectorBR, eigBlockVectorR)
                bpp += np.dot(activeBlockVectorBP, eigBlockVectorP)
            else:
                eigBlockVectorX = eigBlockVector[:sizeX]
                eigBlockVectorR = eigBlockVector[sizeX:]

                pp = np.dot(activeBlockVectorR, eigBlockVectorR)
                app = np.dot(activeBlockVectorAR, eigBlockVectorR)
                bpp = np.dot(activeBlockVectorBR, eigBlockVectorR)

            if verbosityLevel > 10:
                print(pp)
                print(app)
                print(bpp)

            blockVectorX = np.dot(blockVectorX, eigBlockVectorX) + pp
            blockVectorAX = np.dot(blockVectorAX, eigBlockVectorX) + app
            blockVectorBX = np.dot(blockVectorBX, eigBlockVectorX) + bpp

            blockVectorP, blockVectorAP, blockVectorBP = pp, app, bpp

        else:
            if iterationNumber > 0:
                eigBlockVectorX = eigBlockVector[:sizeX]
                eigBlockVectorR = eigBlockVector[sizeX:sizeX+currentBlockSize]
                eigBlockVectorP = eigBlockVector[sizeX+currentBlockSize:]

                pp = np.dot(activeBlockVectorR, eigBlockVectorR)
                pp += np.dot(activeBlockVectorP, eigBlockVectorP)

                app = np.dot(activeBlockVectorAR, eigBlockVectorR)
                app += np.dot(activeBlockVectorAP, eigBlockVectorP)
            else:
                eigBlockVectorX = eigBlockVector[:sizeX]
                eigBlockVectorR = eigBlockVector[sizeX:]

                pp = np.dot(activeBlockVectorR, eigBlockVectorR)
                app = np.dot(activeBlockVectorAR, eigBlockVectorR)

            if verbosityLevel > 10:
                print(pp)
                print(app)

            blockVectorX = np.dot(blockVectorX, eigBlockVectorX) + pp
            blockVectorAX = np.dot(blockVectorAX, eigBlockVectorX) + app

            blockVectorP, blockVectorAP = pp, app

    if B is not None:
        aux = blockVectorBX * _lambda[np.newaxis, :]

    else:
        aux = blockVectorX * _lambda[np.newaxis, :]

    blockVectorR = blockVectorAX - aux

    aux = np.sum(blockVectorR.conjugate() * blockVectorR, 0)
    residualNorms = np.sqrt(aux)

    if verbosityLevel > 0:
        print('final eigenvalue:', _lambda)
        print('final residual norms:', residualNorms)

    if retLambdaHistory:
        if retResidualNormsHistory:
            return _lambda, blockVectorX, lambdaHistory, residualNormsHistory
        else:
            return _lambda, blockVectorX, lambdaHistory
    else:
        if retResidualNormsHistory:
            return _lambda, blockVectorX, residualNormsHistory
        else:
            return _lambda, blockVectorX
