import warnings
cimport numpy as np
import numpy as np
from cython cimport floating
from blas_api cimport fused_dot, fused_nrm2, fused_asum
from coordescendant cimport L11_PENALTY, L21_PENALTY


cdef floating fmax(floating x, floating y) nogil:
    """Max of two real numbers"""
    if x < y:
        return y
    else:
        return x


cdef floating arr_max(int n, floating *a, int inca) nogil:
    "Maximum value of an array"
    cdef int i
    cdef floating m = a[0]
    for i in range(1, n):
        i *= inca
        m = fmax(m, a[i])
    return m


cdef floating abs_max(int n, floating *a, int inca) nogil:
    """b = np.max(np.abs(a))"""
    cdef int i
    cdef floating m = fabs(a[0])
    for i in range(1, n):
        i *= inca
        m = fmax(m, fabs(a[i]))
    return m


cdef floating diff_abs_max(int n, floating *a, int inca, floating *b,
                           int incb) nogil:
    """c = np.max(np.abs(a - b))"""
    cdef int i
    cdef double m = fabs(a[0] - b[0])
    for i in range(1, n):
        m = fmax(m, fabs(a[i * inca] - b[i * incb]))
    return m


cdef void relu(int n, floating *x, int incx) nogil:
    cdef int i
    for i in range(n):
        i *= incx
        if x[i] < 0.:
            x[i] = 0.


cdef floating fused_nrm2_squared(int N, floating *X, int incX) nogil:
    """Computes squared L2 norm of X"""
    return fused_dot(N, X, incX, X, incX)


cdef floating fsign(floating x) nogil:
    if x == 0:
        return 0
    elif x > 0:
        return 1.
    else:
        return -1.


cdef floating compute_dual_gap(int n_samples, int n_features, int n_targets,
                               floating *W_ptr,  # C-order 2d array
                               floating reg, floating l2_reg,
                               floating *X_or_Gram_ptr,  # F-order 2d
                               floating *Y_or_Cov_ptr,  # F-order 2d
                               floating *R_or_RCov_ptr,  # F-order 2d
                               floating *Grad_ptr,  # F-order 2d
                               floating Y_norm2, bint precomputed,
                               int penalty_model,
                               bint positive) nogil:
    """Compute duality gap for penalized regression problems.
    Parameters
    ----------
    W_ptr: pointer to 2d C-order array of model coefficients

    X_or_Gram_ptr: pointer to 2d F-order array of design of Gram matrix
    and points to `X^HX` in `precomputed` mode of `X` in normal mode.

    Y_or_Cov_ptr: pointer to 2d F-order array of target of covariance matrix
    and points to `X^HY` in `precomputed` mode or `Y` in normal mode.

    R_or_RCov_ptr: pointer to 2d F-order array of residuals `X^H(Y - XW)` in
    `precomputed` mode or `Y - XW` in normal mode`

    Grad_ptr: pointer to pre-allocated 2d F-order array. Contents are overwritten

    Y_norm2: Frobenius norm squared of target matrix Y. Not used if if not
    in `precomputed`, else exact value is required.
    """
    with gil:
        if penalty_model not in [L21_PENALTY, L11_PENALTY]:
            raise NotImplementedError(
                "Dual gap computation not implemented for "
                "penalty_model=%s." % penalty_model)
        if precomputed and np.isnan(Y_norm2):
            raise ValueError(
                "`precomputed=True` was specified, expected a value for "
                "`Y_norm2`")

    # matrix sizes
    cdef int W_size = n_features * n_targets
    cdef int R_or_RCov_size = n_samples * n_targets

    # counters
    cdef int i, j, k, f_ind, c_ind

    # more auxiliary variables
    cdef floating R_norm2, ry_sum, W_norm2, C, pen, gap
    cdef floating Grad_axis1norm, dual_norm_Grad

    # compute R_norm2, ry_sum, and Grad
    if precomputed:
        # R_norm2 = ||R||^2_F = ||Y||^2_F - Re(<W, R_or_RCov + Cov>_F)
        # ry_sum = <R, Y>_F = np.sum(R * Y) = ||Y||^2_F - Re(<W, Cov>_F)
        # N.B. we temporarily store R_or_RCov + Y_or_Cov in Grad
        R_norm2 = Y_norm2
        ry_sum = Y_norm2
        for j in range(n_features):
            for k in range(n_targets):
                f_ind = j + n_features * k
                Grad_ptr[f_ind] = R_or_RCov_ptr[f_ind] + Y_or_Cov_ptr[f_ind]
            R_norm2 -= fused_dot(n_targets, W_ptr + j * n_targets, 1,
                                 Grad_ptr + j, n_features)
            ry_sum -= fused_dot(n_targets, W_ptr + j * n_targets, 1,
                                Y_or_Cov_ptr + j, n_features)

        # Grad = RCov - l2_reg * W
        # XXX Is there a faster BLAS way for doing this ?
        for j in range(n_features):
            for k in range(n_targets):
                f_ind = j + n_features * k
                c_ind = k + n_targets * j
                Grad_ptr[f_ind]  = R_or_RCov_ptr[f_ind] - l2_reg * W_ptr[c_ind]
    else:
        # R_norm2 = ||R||^2_F
        R_norm2 = fused_nrm2_squared(R_or_RCov_size, R_or_RCov_ptr, 1)

        # ry_sum = <R, Y>_F
        ry_sum = fused_dot(R_or_RCov_size, R_or_RCov_ptr, 1, Y_or_Cov_ptr, 1)

        # Grad = np.dot(X.conjugate().T, R) - l2_reg * W
        for j in range(n_features):
            for k in range(n_targets):
                f_ind = j + n_features * k
                c_ind = k + n_targets * j
                Grad_ptr[f_ind]  = fused_dot(
                    n_samples, X_or_Gram_ptr + j * n_samples, 1,
                    R_or_RCov_ptr + k * n_samples, 1) - l2_reg * W_ptr[c_ind]

    # dual_norm_Grad = np.max(np.sqrt(np.sum(np.abs(Grad) ** 2, axis=1)))
    dual_norm_Grad = 0.
    for j in range(n_features):
        if penalty_model == L21_PENALTY:
            # Grad_axis1norm = np.sqrt(np.sum(Grad[j] ** 2))
            Grad_axis1norm = fused_nrm2(n_targets, Grad_ptr + j, n_features)
        elif penalty_model == L11_PENALTY:
            # Grad_axis1norm = np.abs(Grad[j]).max()
            if positive:
                Grad_axis1norm = arr_max(n_targets, Grad_ptr + j, n_features)
            else:
                Grad_axis1norm = abs_max(n_targets, Grad_ptr + j, n_features)
        if Grad_axis1norm > dual_norm_Grad:
            dual_norm_Grad = Grad_axis1norm

    # prevent numerical mysteries which may occur if reg = 0
    if -1e-12 < dual_norm_Grad < 0:
        dual_norm_Grad = 0

    # W_norm2 = ||W||^2_F
    W_norm2 = fused_nrm2_squared(W_size, W_ptr, 1)

    # R-scale Grad to ensure dual-feasibility
    if (dual_norm_Grad > reg):
        C = reg / dual_norm_Grad
    else:
        C = 1.

    # data term of primal loss + dual loss
    gap = .5 * (1 + C ** 2) * R_norm2 - C * ry_sum

    # evaluate penalty
    pen = 0.
    for j in range(n_features):
        if penalty_model == L21_PENALTY:
            # pen += np.sqrt(np.sum(W[j] ** 2))
            pen += fused_nrm2(n_targets, W_ptr + j * n_targets, 1)
        else:
            # pen += np.abs(W[j]).sum()
            pen += fused_asum(n_targets, W_ptr + j * n_targets, 1)
    pen *= reg
    pen += .5 * l2_reg * (1. + C ** 2) * W_norm2

    # complete the computation of the dual gap
    gap += pen
    return gap
