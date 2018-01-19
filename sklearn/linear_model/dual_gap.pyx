# encoding: utf-8
# cython: wraparound=False
# cython: boundscheck=False
# cython: cdivision=True
#
# Author: Elvis Dohmatob <gmdopp@gmail.com>
# License: BSD
#
# Synopsis: Helper function for computing dual gap for a some
#           well-known models

import warnings
cimport numpy as np
import numpy as np
from types cimport complexing, floating
from blas_api cimport (fused_nrm2, fused_dotu, fused_dotc)
from utils cimport l1_norm
from cd_fast2 cimport L11_PENALTY, L21_PENALTY

cdef extern from "complex.h" nogil:
    double creal(double complex)
    float crealf(float complex)


cdef floating _compute_dual_gap(int n_samples,
                                int n_features,
                                int n_targets,
                                complexing *W_ptr,  # C-order 2d array
                                floating reg,
                                floating l2_reg,
                                complexing *X_or_Gram_conj_ptr,  # F-order 2d array
                                complexing *Y_or_Cov_ptr,  # F-order 2d array
                                complexing *R_or_RCov_ptr,  # F-order 2d array
                                complexing *Grad_ptr,  # F-order 2d array
                                floating Y_norm2,
                                bint precomputed,
                                int penalty_model) nogil except *:
    """Compute duality gap for penalized regression problems.
    Parameters
    ----------
    W_ptr: pointer to 2d C-order array of model coefficients

    X_or_Gram_conj_ptr: pointer to 2d F-order array of design of Gram matrix
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
                "Dual gap computation not implemented for penalty_model=%s." % (
                    penalty_model))

        if precomputed and np.isnan(Y_norm2):
            raise ValueError(
                "`precomputed=True` was specified, expected a value for `Y_norm2`")

    # specialization of fuzed types / functions
    if complexing is float:
        real_part = crealf
    elif complexing is double:
        real_part = creal
    elif complexing is complex:
        real_part = creal
    else:
        real_part = crealf

    # matrix sizes
    cdef int W_size = n_features * n_targets
    cdef int R_or_RCov_size = n_samples * n_targets

    # counters
    cdef int i, j, k, f_index, c_index

    # more auxiliary variables
    cdef floating R_norm2, W_norm2, C, pen, gap, Grad_axis1norm, dual_norm_Grad
    cdef complexing alpha, beta, ry_sum

    # compute R_norm2, ry_sum, and Grad
    if precomputed:
        # R_norm2 = ||R||^2_F = ||Y||^2_F - Re(<W, R_or_RCov + Cov>_F)
        # ry_sum = <R, Y>_F = np.sum(R * Y) = ||Y||^2_F - Re(<W, Cov>_F)
        # N.B. we temporarily store R_or_RCov + Y_or_Cov in Grad
        alpha = Y_norm2
        beta = Y_norm2
        for j in range(n_features):
            for k in range(n_targets):
                f_index = j + n_features * k
                Grad_ptr[f_index] = R_or_RCov_ptr[f_index] + Y_or_Cov_ptr[f_index]
            alpha -= fused_dotc(n_targets,
                                W_ptr + j * n_targets,
                                1,
                                Grad_ptr + j,
                                n_features)
            beta -= fused_dotc(n_targets,
                               W_ptr + j * n_targets,
                               1,
                               Y_or_Cov_ptr + j,
                               n_features)
        R_norm2 = real_part(alpha)
        ry_sum = real_part(beta)

        # Grad = RCov - l2_reg * W
        # XXX Is there a faster BLAS way for doing this ?
        for j in range(n_features):
            for k in range(n_targets):
                f_index = j + n_features * k
                c_index = k + n_targets * j
                Grad_ptr[f_index]  = R_or_RCov_ptr[f_index] - l2_reg * W_ptr[c_index]
    else:
        # R_norm2 = ||R||^2_F
        R_norm2 = real_part(fused_dotc(R_or_RCov_size,
                                       R_or_RCov_ptr,
                                       1,
                                       R_or_RCov_ptr,
                                       1))

        # ry_sum = <R_conj, Y>_F
        ry_sum = real_part(fused_dotc(R_or_RCov_size,
                                      R_or_RCov_ptr,
                                      1,
                                      Y_or_Cov_ptr,
                                      1))

        # Grad = np.dot(X.conjugate().T, R) - l2_reg * W
        for j in range(n_features):
            for k in range(n_targets):
                f_index = j + n_features * k
                c_index = k + n_targets * j
                Grad_ptr[f_index]  = fused_dotu(n_samples,
                                                X_or_Gram_conj_ptr + j * n_samples,
                                                1,
                                                R_or_RCov_ptr + k * n_samples,
                                                1) - l2_reg * W_ptr[c_index]

    # dual_norm_Grad = np.max(np.sqrt(np.sum(np.abs(Grad) ** 2, axis=1)))
    dual_norm_Grad = 0.
    for j in range(n_features):
        if penalty_model == L21_PENALTY:
            # Grad_axis1norm = np.sqrt(np.sum(Grad[j] ** 2))
            Grad_axis1norm = fused_nrm2(n_targets,
                                        Grad_ptr + j,
                                        n_features)
        else:
            # Grad_axis1norm = np.abs(Grad[j]).sum()
            Grad_axis1norm = l1_norm(n_targets,
                                     Grad_ptr + j,
                                     n_features)
        if Grad_axis1norm > dual_norm_Grad:
            dual_norm_Grad = Grad_axis1norm

    # W_norm2 = ||W||^2_F
    W_norm2 = real_part(fused_dotc(W_size,
                                   W_ptr,
                                   1,
                                   W_ptr,
                                   1))

    # R-scale Grad to ensure dual-feasibility
    if (dual_norm_Grad > reg):
        C = reg / dual_norm_Grad
    else:
        C = 1.

    # primal loss + a piece of the dual loss
    gap = .5 * (1 + C ** 2) * R_norm2

    # evaluate penalty
    pen = 0.
    for j in range(n_features):
        if penalty_model == L21_PENALTY:
            # pen += np.sqrt(np.sum(W[j] ** 2))
            pen += fused_nrm2(n_targets,
                              W_ptr + j * n_targets,
                              1)
        else:
            # pen += np.abs(W[j]).sum()
            pen += l1_norm(n_targets,
                           W_ptr + j * n_targets,
                           1)
    pen *= reg
    pen += .5 * l2_reg * (1. + C ** 2) * W_norm2

    # complete the computation of the dual gap
    if complexing is float or complexing is double:
        gap += pen - C * ry_sum
    else:
        gap += pen - C * real_part(ry_sum)

    return gap


def _py_compute_dual_gap(np.ndarray[complexing, ndim=2, mode="c"] W,
                         floating reg,
                         floating l2_reg,
                         np.ndarray[complexing, ndim=2, mode="fortran"] X_or_Gram_conj,
                         np.ndarray[complexing, ndim=2, mode="fortran"] Y_or_Cov,
                         np.ndarray[complexing, ndim=2, mode="fortran"] R,
                         np.ndarray[complexing, ndim=2, mode="fortran"] Grad,
                         floating Y_norm2=np.nan,
                         int precomputed=True,
                         int penalty_model=L21_PENALTY):
    """Python wrapper for the purpose of testing, etc."""
    cdef int n_samples = len(X_or_Gram_conj)
    cdef int n_features = len(W)
    cdef int n_targets = len(W[0])
    return _compute_dual_gap(n_samples,
                             n_features,
                             n_targets,
                             &W[0, 0],
                             reg,
                             l2_reg,
                             &X_or_Gram_conj[0, 0],
                             &Y_or_Cov[0, 0],
                             &R[0, 0],
                             &Grad[0, 0],
                             Y_norm2,
                             precomputed,
                             penalty_model)

