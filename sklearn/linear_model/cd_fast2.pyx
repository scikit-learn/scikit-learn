cimport numpy as np
import numpy as np
from types cimport floating, complexing
from utils cimport fused_nrm2_squared
import sys
if False:
    sys.path.append("/home/elvis/CODE/FORKED/coordescendant")
    from coordescendant_slow import coordescendant_slow as coordescendant, L11_PENALTY, L21_PENALTY
else:
    from coordescendant import coordescendant, L11_PENALTY, L21_PENALTY


def enet_coordinate_descent(np.ndarray[complexing, ndim=1] w,
                            floating alpha, floating beta,
                            np.ndarray[complexing, ndim=2, mode='fortran'] X,
                            np.ndarray[complexing, ndim=1, mode='c'] y,
                            int max_iter, floating tol, object rng,
                            bint random=0, bint positive=0):
    """Cython version of the coordinate descent algorithm
    for Elastic-Net regression

    We minimize

        (1/2) * norm(y - X w, 2)^2 + alpha norm(w, 1) + (beta/2) norm(w, 2)^2

    """
    cdef np.ndarray[complexing, ndim=2, mode="c"] W
    W = np.ascontiguousarray(w[:, None])
    cdef np.ndarray[complexing, ndim=2, mode="fortran"] Y
    Y = np.asfortranarray(y[:, None])
    W, gap, tol, n_iter = coordescendant(
        W, alpha, beta, X, Y, precomputed=False, max_iter=max_iter, tol=tol,
        penalty_model=L11_PENALTY, rng=rng, random=random, positive=positive)
    return W[:, 0], gap, tol, n_iter


def enet_coordinate_descent_gram(complexing[:] w, floating alpha, floating beta,
                                 np.ndarray[complexing, ndim=2, mode='c'] Q,
                                 np.ndarray[complexing, ndim=1, mode='c'] q,
                                 np.ndarray[complexing, ndim=1] y,
                                 int max_iter, floating tol, object rng,
                                 bint random=0, bint positive=0):
    """Cython version of the coordinate descent algorithm
    for Elastic-Net regression

    We minimize

        (1/2) * w^T Q w - q^T w + alpha norm(w, 1) + (beta/2) * norm(w, 2)^2

    which amounts to the Elastic-Net problem when:

        Q = X^T X (Gram matrix)
        q = X^T y
    """
    cdef np.ndarray[complexing, ndim=2, mode="c"] W
    W = np.ascontiguousarray(w[:, None])
    cdef np.ndarray[complexing, ndim=2, mode="fortran"] Gram
    Gram = np.asfortranarray(Q)
    cdef np.ndarray[complexing, ndim=2, mode="fortran"] Cov
    Cov = np.asfortranarray(q[:, None])
    cdef floating Y_norm2 = <floating>fused_nrm2_squared(len(y), &y[0], 1)
    W, gap, tol, n_iter = coordescendant(
        W, alpha, beta, Gram, Cov, Y_norm2=Y_norm2, precomputed=True,
        max_iter=max_iter, tol=tol, penalty_model=L11_PENALTY, rng=rng,
        random=random, positive=positive)
    return W[:, 0], gap, tol, n_iter


def enet_coordinate_descent_multi_task(complexing[::1, :] W, floating l1_reg,
                                       floating l2_reg,
                                       np.ndarray[complexing, ndim=2, mode='fortran'] X,
                                       np.ndarray[complexing, ndim=2] Y,
                                       int max_iter, floating tol, object rng,
                                       bint random=0):
    """Cython version of the coordinate descent algorithm
        for Elastic-Net mult-task regression

        We minimize

        0.5 * ||y - X w||^2_2 + l1_reg ||w||_21 + 0.5 * l2_reg ||w||^2_2

    """
    cdef np.ndarray[complexing, ndim=2, mode="c"] W_c
    W_c = np.ascontiguousarray(W.T)
    cdef np.ndarray[complexing, ndim=2, mode="fortran"] Y_f
    Y_f = np.asfortranarray(Y)
    W_c, gap, tol, n_iter = coordescendant(
        W_c, l1_reg, l2_reg, X, Y_f, precomputed=False, max_iter=max_iter,
        tol=tol, penalty_model=L21_PENALTY, rng=rng, random=random)
    return W_c.T, gap, tol, n_iter
