import numpy as np
from sklearn.utils import check_random_state
from sklearn.linear_model.coordescendant import (L2INF_CONSTRAINT, L21_PENALTY,
                                                 L1INF_CONSTRAINT,
                                                 L11_PENALTY, NOP)
from .prox_slow import nrm2
from .dual_gap_slow import compute_dual_gap_slow
from .prox_slow import proj_l1_slow, proj_l2_slow, prox_l2_slow, prox_l1_slow


def coordescendant_slow(W, reg, l2_reg, X_or_Gram, Y_or_Cov, Y_norm2=np.nan,
                        precomputed=True, penalty_model=NOP, max_iter=100,
                        tol=0., random=False, rng=None, positive=False,
                        emulate_sklearn_dl=False, verbose=0):
    # some sanity checks
    if positive:
        if penalty_model != L11_PENALTY:
            raise ValueError(
                "positive=True only makes sense when L1 penalty is imposed.")
        for stuff in W, X_or_Gram, Y_or_Cov:
            if stuff.dtype not in [np.float32, np.float64]:
                raise TypeError(
                    "positive=True for complex data makes no sense")

    # misc
    rng = check_random_state(rng)
    n_features, n_targets = W.shape

    # shall we be computing dual gaps ?
    dual_gap_available = False
    if penalty_model in [L11_PENALTY, L21_PENALTY]:
        dual_gap_available = True

    # setup the proximal operator / penalty-enforcer
    if hasattr(penalty_model, "__call__"):
        prox = penalty_model
    else:
        if penalty_model == NOP:
            prox = None
        elif penalty_model == L11_PENALTY:
            prox = prox_l1_slow
        elif penalty_model == L21_PENALTY:
            prox = prox_l2_slow
        elif penalty_model == L2INF_CONSTRAINT:
            prox = proj_l2_slow
        elif penalty_model == L1INF_CONSTRAINT:
            prox = proj_l1_slow
        else:
            raise NotImplementedError("penalty_model=%s" % penalty_model)

        # sklearn only supports contraints (not penalties)
        if penalty_model > 0:
            emulate_sklearn_dl = False

    # initialize residuals
    R = np.asfortranarray(Y_or_Cov - np.dot(X_or_Gram, W))

    # misc
    Grad = np.zeros_like(W, order="F")
    X_or_Gram_conj = X_or_Gram.conjugate()

    # pre-compute column norms of design matrix
    if not precomputed:
        X_col_squared_norms = [nrm2(X_or_Gram[:, j], squared=True)
                               for j in range(n_features)]

    # the main loop
    d_W_abs_max = 0.
    W_abs_max = 0.
    n_iter = 0
    d_W_abs_tol = tol
    gap = tol + 1.
    if not precomputed:
        tol *= nrm2(Y_or_Cov.ravel(), squared=True)
    for n_iter in range(max_iter):
        W_abs_max = 0.
        d_W_abs_max = 0.
        for j in range(n_features):
            # pick a coordinate to update
            if random:
                j = rng.randint(n_features)

            # coordinate-wise Lipschitz constant
            if precomputed:
                ajj = X_or_Gram[j, j].real
            else:
                ajj = X_col_squared_norms[j]
            ajj += l2_reg

            # skip if dead coordinate
            if ajj == 0. and not emulate_sklearn_dl:
                continue

            # rank-1 update
            R += np.outer(X_or_Gram[:, j], W[j])

            # unconstrained update
            Wj = W[j].copy()
            if emulate_sklearn_dl or not precomputed:
                if emulate_sklearn_dl:
                    alpha = ajj
                else:
                    alpha = 1.
                W[j] = np.dot(X_or_Gram[:, j].conjugate(), R)
                if alpha != 1.:
                    W[j] *= alpha
            else:
                W[j] = R[j].copy()

            # proximal update
            if positive:
                W[j] = np.maximum(W[j], 0.)
            if prox is None:
                W[j] /= ajj
            else:
                prox(W[j], reg, ajj)

            # rank-1 update
            R -= np.outer(X_or_Gram[:, j], W[j])

            # update the maximum absolute coefficient update
            Wj_abs_max = np.max(np.abs(W[j]))
            W_abs_max = max(W_abs_max, Wj_abs_max)
            d_Wj_abs_max = np.max(np.abs(W[j] - Wj))
            d_W_abs_max = max(d_W_abs_max, d_Wj_abs_max)

        if (W_abs_max == 0. or d_W_abs_max / W_abs_max < d_W_abs_tol or
                n_iter == max_iter - 1):
            if dual_gap_available:
                gap = compute_dual_gap_slow(
                    W, reg, l2_reg, X_or_Gram_conj, Y_or_Cov, R, Grad,
                    Y_norm2=Y_norm2, precomputed=precomputed,
                    penalty_model=penalty_model, positive=positive)
            else:
                gap = tol

            # XXX TODO: check for bad things like negative gap, etc.
            pass

            # exit if we reached desired tolerance
            if gap < tol:
                break

    return W, gap, tol, n_iter + 1
