import warnings
import numpy as np
from sklearn.linear_model.coordescendant import L21_PENALTY, L11_PENALTY
from .prox_slow import nrm2


def compute_dual_gap_slow(W, reg, l2_reg, X_or_Gram_conj, Y_or_Cov, R_or_RCov,
                          Grad, Y_norm2=None, precomputed=True,
                          penalty_model=L21_PENALTY, positive=False):
    """Compute duality gap for penalized regression problems"""
    if penalty_model not in [L21_PENALTY, L11_PENALTY]:
        raise NotImplementedError(
            "Dual gap computation not implemented for penalty_model=%s." % (
                penalty_model))

    if precomputed and np.isnan(Y_norm2):
        raise ValueError("`precomputed=True` was specified, expected a value for `Y_norm2`")

    n_features, n_targets = W.shape
    n_samples = len(X_or_Gram_conj)
    if len(R_or_RCov) != n_samples or len(R_or_RCov[0]) != n_targets or \
       R_or_RCov.shape != Y_or_Cov.shape:
        raise RuntimeError("Inconsistent matrix shapes")

    if precomputed:
        # R_norm2 = ||R||^2_F = ||Y||^2_F - Re(<W, R_or_RCov + Cov>_F)
        R_norm2 = Y_norm2 - np.dot(W.ravel().conjugate(),
                                   np.ravel(R_or_RCov + Y_or_Cov))

        # ry_sum = <R, Y>_F = ||Y||^2_F - Re(<W, Cov>_F)
        ry_sum = Y_norm2 - np.dot(W.ravel().conjugate(), Y_or_Cov.ravel())

        # Grad = R_or_RCov - l2_reg * W
        Grad = R_or_RCov - l2_reg * W
    else:
        # R_norm2 = ||R||^2_F = ||Y||^2_F
        # ry_sum = <R, Y>_F
        R_flat = R_or_RCov.ravel() 
        R_norm2 = np.dot(R_flat.conjugate(), R_flat)
        ry_sum = np.real(np.dot(R_flat.conjugate(), Y_or_Cov.ravel()))
        # Grad = np.dot(X.T, R) - l2_reg * W
        Grad = np.dot(X_or_Gram_conj.T, R_or_RCov) - l2_reg * W
    R_norm2 = np.real(R_norm2)
    ry_sum = np.real(ry_sum)

    # dual_norm_Grad = np.max(np.sqrt(np.sum(Grad ** 2, axis=1)))
    dual_norm_Grad = 0.
    for j in range(n_features):
        # np.sqrt(np.sum(Grad ** 2, axis=1))
        if penalty_model == L21_PENALTY:
            Grad_axis1norm = nrm2(Grad[j])
        elif penalty_model == L11_PENALTY:
            if positive:
                Grad_axis1norm = np.max(Grad[j])
            else:
                Grad_axis1norm = np.max(np.abs(Grad[j]))
        if Grad_axis1norm > dual_norm_Grad:
            dual_norm_Grad = Grad_axis1norm

    if -1e-12 < dual_norm_Grad < 0:
        dual_norm_Grad = 0

    # check dual feasibility
    if (dual_norm_Grad > reg):
        const = reg / dual_norm_Grad
    else:
        const = 1.
    gap = .5 * (1 + const ** 2) * R_norm2 - const * ry_sum

    # W_norm2 = ||W||^2_F
    W_norm2 = nrm2(W.ravel(), squared=True)

    # l21_norm = np.sqrt(np.sum(W ** 2, axis=1)).sum()
    pen = 0.
    for j in range(n_features):
        # np.sqrt(np.sum(W ** 2, axis=1))
        if penalty_model == L21_PENALTY:
            pen += nrm2(W[j])
        elif penalty_model == L11_PENALTY:
            pen += np.abs(W[j]).sum()
        else:
            raise NotImplementedError
    pen *= reg
    pen += .5 * l2_reg * (1. + const ** 2) * W_norm2
    gap += pen

    return gap
