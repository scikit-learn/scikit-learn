# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from libc.math cimport fabs, sqrt
from libc.stdlib cimport free, malloc
from libc.string cimport memset
# from libc.time cimport time, time_t
import numpy as np

from cython cimport floating
import warnings
from scipy import sparse
from sklearn.exceptions import ConvergenceWarning

# Note: The use of BLAS can cause random results within floating point
# arithmetic, e.g. the summation order is not deterministic.
from sklearn.utils._cython_blas cimport (
    _axpy, _dot, _asum, _gemv, _nrm2, _copy, _scal
)
from sklearn.utils._cython_blas cimport ColMajor, Trans, NoTrans
from sklearn.utils._typedefs cimport float64_t, int32_t, uint8_t, uint32_t
from sklearn.utils._random cimport our_rand_r

from sklearn.linear_model._linear_loss import Multinomial_LDL_Decomposition


cdef extern from "<float.h>":
    const float FLT_EPSILON
    const double DBL_EPSILON

# The following two functions are shamelessly copied from the tree code.

cdef enum:
    # Max value for our rand_r replacement (near the bottom).
    # We don't use RAND_MAX because it's different across platforms and
    # particularly tiny on Windows/MSVC.
    # It corresponds to the maximum representable value for
    # 32-bit signed integers (i.e. 2^31 - 1).
    RAND_R_MAX = 2147483647


cdef inline uint32_t rand_int(uint32_t end, uint32_t* random_state) noexcept nogil:
    """Generate a random integer in [0; end)."""
    return our_rand_r(random_state) % end


cdef inline floating fmax(floating x, floating y) noexcept nogil:
    if x > y:
        return x
    return y


cdef inline floating fsign(floating f) noexcept nogil:
    if f == 0:
        return 0
    elif f > 0:
        return 1.0
    else:
        return -1.0


cdef inline floating abs_max(int n, const floating* a) noexcept nogil:
    """np.max(np.abs(a))"""
    cdef int i
    cdef floating m = fabs(a[0])
    cdef floating d
    for i in range(1, n):
        d = fabs(a[i])
        if d > m:
            m = d
    return m


cdef inline floating max(int n, floating* a) noexcept nogil:
    """np.max(a)"""
    cdef int i
    cdef floating m = a[0]
    cdef floating d
    for i in range(1, n):
        d = a[i]
        if d > m:
            m = d
    return m


cdef inline floating diff_abs_max(int n, const floating* a, floating* b) noexcept nogil:
    """np.max(np.abs(a - b))"""
    cdef int i
    cdef floating m = fabs(a[0] - b[0])
    cdef floating d
    for i in range(1, n):
        d = fabs(a[i] - b[i])
        if d > m:
            m = d
    return m


message_conv = (
    "Objective did not converge. You might want to increase "
    "the number of iterations, check the scale of the "
    "features or consider increasing regularisation."
)


message_ridge = (
    "Linear regression models with a zero l1 penalization "
    "strength are more efficiently fitted using one of the "
    "solvers implemented in "
    "sklearn.linear_model.Ridge/RidgeCV instead."
)


cdef inline floating dual_gap_formulation_A(
    floating alpha,  # L1 penalty
    floating beta,  # L1 penalty
    floating w_l1_norm,
    floating w_l2_norm2,
    floating R_norm2,  # R @ R
    floating Ry,  # R @ y
    floating dual_norm_XtA,
    bint gap_smaller_eps,
) noexcept nogil:
    """Compute dual gap according to formulation A."""
    cdef floating gap, primal, dual
    cdef floating scale  # Scaling factor to achieve dual feasible point.

    if floating is float:
        eps = FLT_EPSILON
    else:
        eps = DBL_EPSILON

    primal = 0.5 * (R_norm2 + beta * w_l2_norm2) + alpha * w_l1_norm

    if (dual_norm_XtA > alpha):
        scale = alpha / dual_norm_XtA
    else:
        scale = 1.0
    dual = -0.5 * (scale ** 2) * (R_norm2 + beta * w_l2_norm2) + scale * Ry
    gap = primal - dual
    if gap_smaller_eps and abs(gap) <= 2 * eps * primal:
        gap = 0.0
    return gap


cdef (floating, floating) gap_enet(
    int n_samples,
    int n_features,
    const floating[::1] w,
    floating alpha,  # L1 penalty
    floating beta,  # L2 penalty
    const floating[::1, :] X,
    const floating[::1] y,
    const floating[::1] R,  # current residuals = y - X @ w
    floating[::1] XtA,  # XtA = X.T @ R - beta * w is calculated inplace
    bint positive,
    bint gap_smaller_eps,
) noexcept nogil:
    """Compute dual gap for use in enet_coordinate_descent.

    alpha > 0:            formulation A of the duality gap
    alpha = 0 & beta > 0: formulation B of the duality gap
    alpha = beta = 0:     OLS first order condition (=gradient)
    """
    cdef floating gap, primal, dual
    cdef floating dual_norm_XtA
    cdef floating R_norm2
    cdef floating Ry
    cdef floating w_l1_norm
    cdef floating w_l2_norm2 = 0.0

    if floating is float:
        eps = FLT_EPSILON
    else:
        eps = DBL_EPSILON

    # w_l2_norm2 = w @ w
    if beta > 0:
        w_l2_norm2 = _dot(n_features, &w[0], 1, &w[0], 1)
    # R_norm2 = R @ R
    R_norm2 = _dot(n_samples, &R[0], 1, &R[0], 1)
    # Ry = R @ y
    if not (alpha == 0 and beta == 0):
        Ry = _dot(n_samples, &R[0], 1, &y[0], 1)

    if alpha == 0:
        # XtA = X.T @ R
        _gemv(
            ColMajor, Trans, n_samples, n_features, 1.0, &X[0, 0],
            n_samples, &R[0], 1, 0, &XtA[0], 1,
        )
        # ||X'R||_2^2
        dual_norm_XtA = _dot(n_features, &XtA[0], 1, &XtA[0], 1)
        if beta == 0:
            # This is OLS, no dual gap available. Resort to first order condition
            #     X'R = 0
            #     gap = ||X'R||_2^2
            # Compare with stopping criterion of LSQR.
            gap = dual_norm_XtA
            return gap, dual_norm_XtA
        # This is Ridge regression, we use formulation B for the dual gap.
        primal = 0.5 * (R_norm2 + beta * w_l2_norm2)
        dual = -0.5 * R_norm2 + Ry - 1 / (2 * beta) * dual_norm_XtA
        gap = primal - dual
        if gap_smaller_eps and abs(gap) <= 2 * eps * primal:
            gap = 0.0
        return gap, dual_norm_XtA

    # XtA = X.T @ R - beta * w
    _copy(n_features, &w[0], 1, &XtA[0], 1)
    _gemv(ColMajor, Trans, n_samples, n_features, 1.0, &X[0, 0],
          n_samples, &R[0], 1,
          -beta, &XtA[0], 1)

    # dual_norm_XtA
    if positive:
        dual_norm_XtA = max(n_features, &XtA[0])
    else:
        dual_norm_XtA = abs_max(n_features, &XtA[0])

    # w_l1_norm = np.sum(np.abs(w))
    w_l1_norm = _asum(n_features, &w[0], 1)

    gap = dual_gap_formulation_A(
        alpha=alpha,
        beta=beta,
        w_l1_norm=w_l1_norm,
        w_l2_norm2=w_l2_norm2,
        R_norm2=R_norm2,
        Ry=Ry,
        dual_norm_XtA=dual_norm_XtA,
        gap_smaller_eps=gap_smaller_eps,
    )
    return gap, dual_norm_XtA


def enet_coordinate_descent(
    floating[::1] w,
    floating alpha,
    floating beta,
    const floating[::1, :] X,
    const floating[::1] y,
    unsigned int max_iter,
    floating tol,
    object rng,
    bint random=0,
    bint positive=0,
    bint do_screening=1,
    bint early_stopping=1,
):
    """
    Cython version of the coordinate descent algorithm for Elastic-Net regression.

    The algorithm mostly follows [Friedman 2010].
    We minimize the primal

        P(w) = 1/2 ||y - X w||_2^2 + alpha ||w||_1 + beta/2 ||w||_2^2

    The dual for beta = 0, see e.g. [Fercoq 2015] with v = alpha * theta, is

        D(v) = -1/2 ||v||_2^2 + y' v    (formulation A)

    with dual feasible condition ||X^T v||_inf <= alpha.
    For beta > 0, one uses extended versions of X and y by adding n_features rows

        X -> (           X)    y -> (y)
             (sqrt(beta) I)         (0)

    Note that the residual R = y - X w is an important ingredient for the estimation of
    a dual feasible point v.
    At optimum of primal w* and dual v*, one has

        v* = y - X w*

    The duality gap is

        G(w, v) = P(w) - D(v) <= P(w) - P(w*)

    Strong duality holds: G(w*, v*) = 0.
    For testing convergence, one uses G(w, v) with current w and uses

        v = R                            if ||X^T R||_inf <= alpha
        v = R * alpha / ||X^T R||_inf    else

    The final stopping criterion is based on the duality gap

        tol ||y||_2^2 <= G(w, v)

    The tolerance here is multiplied by ||y||_2^2 to have an inequality that scales the
    same on both sides and because one has G(0, 0) = 1/2 ||y||_2^2.

    Note:
    The above dual D(v) and duality gap G require alpha > 0 because of the dual
    feasible condition.
    There is, however, an alternative dual formulation, see [Dünner 2016] 5.2.3 and
    https://github.com/scikit-learn/scikit-learn/issues/22836:

        D(v) = -1/2 ||v||_2^2 + y' v
               -1/(2 beta) sum_j (|X_j' v| - alpha)_+^2    (formulation B)

    The dual feasible set is v element real numbers. It requires beta > 0, but
    alpha = 0 is allowed. Strong duality holds and at optimum, v* = y - X w*.

    Returns
    -------
    w : ndarray of shape (n_features,)
        ElasticNet coefficients.
    gap : float
        Achieved dual gap.
    tol : float
        Equals input `tol` times `np.dot(y, y)`. The tolerance used for the dual gap.
    n_iter : int
        Number of coordinate descent iterations.

    References
    ----------
    .. [Friedman 2010]
       Jerome H. Friedman, Trevor Hastie, Rob Tibshirani. (2010)
       Regularization Paths for Generalized Linear Models via Coordinate Descent
       https://www.jstatsoft.org/article/view/v033i01

    .. [Fercoq 2015]
       Olivier Fercoq, Alexandre Gramfort, Joseph Salmon. (2015)
       Mind the duality gap: safer rules for the Lasso
       https://arxiv.org/abs/1505.03410

    .. [Dünner 2016]
       Celestine Dünner, Simon Forte, Martin Takác, Martin Jaggi. (2016).
       Primal-Dual Rates and Certificates. In ICML 2016.
       https://arxiv.org/abs/1602.05205
    """

    if floating is float:
        dtype = np.float32
    else:
        dtype = np.float64

    # get the data information into easy vars
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]

    # compute squared norms of the columns of X
    # same as norm2_cols_X = np.square(X).sum(axis=0)
    cdef floating[::1] norm2_cols_X = np.einsum(
        "ij,ij->j", X, X, dtype=dtype, order="C"
    )

    # initial value of the residuals
    cdef floating[::1] R = np.empty(n_samples, dtype=dtype)
    cdef floating[::1] XtA = np.empty(n_features, dtype=dtype)

    cdef floating d_j
    cdef floating Xj_theta
    cdef floating tmp
    cdef floating w_j
    cdef floating d_w_max
    cdef floating w_max
    cdef floating d_w_j
    cdef floating gap = tol + 1.0
    cdef floating d_w_tol = tol
    cdef floating dual_norm_XtA
    cdef unsigned int n_active = n_features
    cdef uint32_t[::1] active_set
    # TODO: use binset instead of array of bools
    cdef uint8_t[::1] excluded_set
    cdef unsigned int j
    cdef unsigned int n_iter = 0
    cdef unsigned int f_iter
    cdef uint32_t rand_r_state_seed = rng.randint(0, RAND_R_MAX)
    cdef uint32_t* rand_r_state = &rand_r_state_seed

    if alpha == 0:
        # No screeing without L1-penalty.
        do_screening = False

    if do_screening:
        active_set = np.empty(n_features, dtype=np.uint32)  # map [:n_active] -> j
        excluded_set = np.empty(n_features, dtype=np.uint8)

    with nogil:
        # R = y - np.dot(X, w)
        _copy(n_samples, &y[0], 1, &R[0], 1)
        _gemv(ColMajor, NoTrans, n_samples, n_features, -1.0, &X[0, 0],
              n_samples, &w[0], 1, 1.0, &R[0], 1)

        # tol *= np.dot(y, y)
        tol *= _dot(n_samples, &y[0], 1, &y[0], 1)

        # Check convergence before entering the main loop.
        gap, dual_norm_XtA = gap_enet(
            n_samples, n_features, w, alpha, beta, X, y, R, XtA, positive, False
        )
        if early_stopping and gap <= tol:
            with gil:
                return np.asarray(w), gap, tol, 0

        # Gap Safe Screening Rules, see https://arxiv.org/abs/1802.07481, Eq. 11
        if do_screening:
            n_active = 0
            for j in range(n_features):
                if norm2_cols_X[j] == 0:
                    w[j] = 0
                    excluded_set[j] = 1
                    continue
                Xj_theta = XtA[j] / fmax(alpha, dual_norm_XtA)  # X[:,j] @ dual_theta
                d_j = (1 - fabs(Xj_theta)) / sqrt(norm2_cols_X[j] + beta)
                if d_j <= sqrt(2 * gap) / alpha:
                    # include feature j
                    active_set[n_active] = j
                    excluded_set[j] = 0
                    n_active += 1
                else:
                    if w[j] != 0:
                        # R += w[j] * X[:,j]
                        _axpy(n_samples, w[j], &X[0, j], 1, &R[0], 1)
                        w[j] = 0
                    excluded_set[j] = 1

        for n_iter in range(max_iter):
            w_max = 0.0
            d_w_max = 0.0
            for f_iter in range(n_active):  # Loop over coordinates
                if random:
                    j = rand_int(n_active, rand_r_state)
                else:
                    j = f_iter

                if do_screening:
                    j = active_set[j]

                if norm2_cols_X[j] == 0.0:
                    continue

                w_j = w[j]  # Store previous value

                # tmp = X[:,j] @ (R + w_j * X[:,j])
                tmp = _dot(n_samples, &X[0, j], 1, &R[0], 1) + w_j * norm2_cols_X[j]

                if positive and tmp < 0:
                    w[j] = 0.0
                else:
                    w[j] = (fsign(tmp) * fmax(fabs(tmp) - alpha, 0)
                            / (norm2_cols_X[j] + beta))

                if w[j] != w_j:
                    # R -= (w[j] - w_j) * X[:,j] # Update residual
                    _axpy(n_samples, w_j - w[j], &X[0, j], 1, &R[0], 1)

                # update the maximum absolute coefficient update
                d_w_j = fabs(w[j] - w_j)
                d_w_max = fmax(d_w_max, d_w_j)

                w_max = fmax(w_max, fabs(w[j]))

            if (
                w_max == 0.0
                or d_w_max / w_max <= d_w_tol
                or n_iter == max_iter - 1
                or n_active <= 1  # We have an analytical exact solution.
            ):
                # the biggest coordinate update of this iteration was smaller
                # than the tolerance: check the duality gap as ultimate
                # stopping criterion
                gap, dual_norm_XtA = gap_enet(
                    n_samples, n_features, w, alpha, beta, X, y, R, XtA, positive, True
                )
                if gap <= tol:
                    # return if we reached desired tolerance
                    break

                # Gap Safe Screening Rules, see https://arxiv.org/abs/1802.07481, Eq. 11
                if do_screening:
                    n_active = 0
                    for j in range(n_features):
                        if excluded_set[j]:
                            continue
                        Xj_theta = XtA[j] / fmax(alpha, dual_norm_XtA)  # X @ dual_theta
                        d_j = (1 - fabs(Xj_theta)) / sqrt(norm2_cols_X[j] + beta)
                        if d_j <= sqrt(2 * gap) / alpha:
                            # include feature j
                            active_set[n_active] = j
                            excluded_set[j] = 0
                            n_active += 1
                        else:
                            if w[j] != 0:
                                # R += w[j] * X[:,j]
                                _axpy(n_samples, w[j], &X[0, j], 1, &R[0], 1)
                                w[j] = 0
                            excluded_set[j] = 1

        else:
            # for/else, runs if for doesn't end with a `break`
            with gil:
                message = (
                    message_conv +
                    f" Duality gap: {gap:.6e}, tolerance: {tol:.3e}"
                )
                if alpha < np.finfo(np.float64).eps:
                    message += "\n" + message_ridge
                warnings.warn(message, ConvergenceWarning)

    return np.asarray(w), gap, tol, n_iter + 1


cdef inline void R_plus_wj_Xj(
    unsigned int n_samples,
    floating[::1] R,  # out
    const floating[::1] X_data,
    const int[::1] X_indices,
    const int[::1] X_indptr,
    const floating[::1] X_mean,
    bint center,
    const floating[::1] sample_weight,
    bint no_sample_weights,
    floating w_j,
    unsigned int j,
) noexcept nogil:
    """R += w_j * X[:,j]"""
    cdef unsigned int startptr = X_indptr[j]
    cdef unsigned int endptr = X_indptr[j + 1]
    cdef floating sw
    cdef floating X_mean_j = X_mean[j]
    if no_sample_weights:
        for i in range(startptr, endptr):
            R[X_indices[i]] += X_data[i] * w_j
        if center:
            for i in range(n_samples):
                R[i] -= X_mean_j * w_j
    else:
        for i in range(startptr, endptr):
            sw = sample_weight[X_indices[i]]
            R[X_indices[i]] += sw * X_data[i] * w_j
        if center:
            for i in range(n_samples):
                R[i] -= sample_weight[i] * X_mean_j * w_j


cdef (floating, floating) gap_enet_sparse(
    int n_samples,
    int n_features,
    const floating[::1] w,
    floating alpha,  # L1 penalty
    floating beta,  # L2 penalty
    const floating[::1] X_data,
    const int[::1] X_indices,
    const int[::1] X_indptr,
    const floating[::1] y,
    const floating[::1] sample_weight,
    bint no_sample_weights,
    const floating[::1] X_mean,
    bint center,
    const floating[::1] R,  # current residuals = y - X @ w
    floating R_sum,
    floating[::1] XtA,  # XtA = X.T @ R - beta * w is calculated inplace
    bint positive,
    bint gap_smaller_eps,
) noexcept nogil:
    """Compute dual gap for use in sparse_enet_coordinate_descent.

    alpha > 0:            formulation A of the duality gap
    alpha = 0 & beta > 0: formulation B of the duality gap
    alpha = beta = 0:     OLS first order condition (=gradient)
    """
    cdef floating gap, primal, dual
    cdef floating dual_norm_XtA
    cdef floating R_norm2
    cdef floating Ry
    cdef floating w_l1_norm
    cdef floating w_l2_norm2 = 0.0
    cdef unsigned int i, j

    if floating is float:
        eps = FLT_EPSILON
    else:
        eps = DBL_EPSILON

    # w_l2_norm2 = w @ w
    if beta > 0:
        w_l2_norm2 = _dot(n_features, &w[0], 1, &w[0], 1)
    # R_norm2 = R @ R
    if no_sample_weights:
        R_norm2 = _dot(n_samples, &R[0], 1, &R[0], 1)
    else:
        R_norm2 = 0.0
        for i in range(n_samples):
            # R is already multiplied by sample_weight
            if sample_weight[i] != 0:
                R_norm2 += (R[i] ** 2) / sample_weight[i]
    # Ry = R @ y
    if not (alpha == 0 and beta == 0):
        # Note that with sample_weight, R equals R*sw and y is just y, such that
        # Ry = (sw * R) @ y, as it should be.
        Ry = _dot(n_samples, &R[0], 1, &y[0], 1)

    if alpha == 0:
        # XtA = X.T @ R
        for j in range(n_features):
            XtA[j] = 0.0
            for i in range(X_indptr[j], X_indptr[j + 1]):
                XtA[j] += X_data[i] * R[X_indices[i]]

            if center:
                XtA[j] -= X_mean[j] * R_sum
        # ||X'R||_2^2
        dual_norm_XtA = _dot(n_features, &XtA[0], 1, &XtA[0], 1)
        if beta == 0:
            # This is OLS, no dual gap available. Resort to first order condition
            #     X'R = 0
            #     gap = ||X'R||_2^2
            # Compare with stopping criterion of LSQR.
            gap = dual_norm_XtA
            return gap, dual_norm_XtA
        # This is Ridge regression, we use formulation B for the dual gap.
        primal = 0.5 * (R_norm2 + beta * w_l2_norm2)
        dual = -0.5 * R_norm2 + Ry - 1 / (2 * beta) * dual_norm_XtA
        gap = primal - dual
        if gap_smaller_eps and abs(gap) <= 2 * eps * primal:
            gap = 0.0
        return gap, dual_norm_XtA

    # XtA = X.T @ R - beta * w
    # sparse X.T @ dense R
    for j in range(n_features):
        XtA[j] = 0.0
        for i in range(X_indptr[j], X_indptr[j + 1]):
            XtA[j] += X_data[i] * R[X_indices[i]]

        if center:
            XtA[j] -= X_mean[j] * R_sum
        XtA[j] -= beta * w[j]

    # dual_norm_XtA
    if positive:
        dual_norm_XtA = max(n_features, &XtA[0])
    else:
        dual_norm_XtA = abs_max(n_features, &XtA[0])

    # w_l1_norm = np.sum(np.abs(w))
    w_l1_norm = _asum(n_features, &w[0], 1)

    gap = dual_gap_formulation_A(
        alpha=alpha,
        beta=beta,
        w_l1_norm=w_l1_norm,
        w_l2_norm2=w_l2_norm2,
        R_norm2=R_norm2,
        Ry=Ry,
        dual_norm_XtA=dual_norm_XtA,
        gap_smaller_eps=gap_smaller_eps,
    )
    return gap, dual_norm_XtA


def sparse_enet_coordinate_descent(
    floating[::1] w,
    floating alpha,
    floating beta,
    const floating[::1] X_data,
    const int[::1] X_indices,
    const int[::1] X_indptr,
    const floating[::1] y,
    const floating[::1] sample_weight,
    const floating[::1] X_mean,
    unsigned int max_iter,
    floating tol,
    object rng,
    bint random=0,
    bint positive=0,
    bint do_screening=1,
    bint early_stopping=1,
):
    """Cython version of the coordinate descent algorithm for Elastic-Net

    We minimize:

        1/2 * norm(y - Z w, 2)^2 + alpha * norm(w, 1) + (beta/2) * norm(w, 2)^2

    where Z = X - X_mean.
    With sample weights sw, this becomes

        1/2 * sum(sw * (y - Z w)^2, axis=0) + alpha * norm(w, 1)
        + (beta/2) * norm(w, 2)^2

    and X_mean is the weighted average of X (per column).

    The rest is the same as enet_coordinate_descent, but for sparse X.

    Returns
    -------
    w : ndarray of shape (n_features,)
        ElasticNet coefficients.
    gap : float
        Achieved dual gap.
    tol : float
        Equals input `tol` times `np.dot(y, y)`. The tolerance used for the dual gap.
    n_iter : int
        Number of coordinate descent iterations.
    """
    # Notes for sample_weight:
    # For dense X, one centers X and y and then rescales them by sqrt(sample_weight).
    # Here, for sparse X, we get the sample_weight averaged center X_mean. We take care
    # that every calculation results as if we had rescaled y and X (and therefore also
    # X_mean) by sqrt(sample_weight) without actually calculating the square root.
    # We work with:
    #     yw = sample_weight * y
    #     R = sample_weight * residual
    #     norm2_cols_X = np.sum(sample_weight * (X - X_mean)**2, axis=0)

    if floating is float:
        dtype = np.float32
    else:
        dtype = np.float64

    # get the data information into easy vars
    cdef unsigned int n_samples = y.shape[0]
    cdef unsigned int n_features = w.shape[0]

    # compute squared norms of the columns of X
    cdef floating[::1] norm2_cols_X = np.zeros(n_features, dtype=dtype)

    # initial value of the residuals
    # R = y - Zw, weighted version R = sample_weight * (y - Zw)
    cdef floating[::1] R
    cdef floating[::1] XtA = np.empty(n_features, dtype=dtype)
    cdef const floating[::1] yw

    cdef floating d_j
    cdef floating Xj_theta
    cdef floating tmp
    cdef floating w_j
    cdef floating d_w_max
    cdef floating w_max
    cdef floating d_w_j
    cdef floating gap = tol + 1.0
    cdef floating d_w_tol = tol
    cdef floating dual_norm_XtA
    cdef floating X_mean_j
    cdef floating R_sum = 0.0
    cdef floating normalize_sum
    cdef unsigned int n_active = n_features
    cdef uint32_t[::1] active_set
    # TODO: use binset insteaf of array of bools
    cdef uint8_t[::1] excluded_set
    cdef unsigned int i
    cdef unsigned int j
    cdef unsigned int n_iter = 0
    cdef unsigned int f_iter
    cdef unsigned int startptr = X_indptr[0]
    cdef unsigned int endptr
    cdef uint32_t rand_r_state_seed = rng.randint(0, RAND_R_MAX)
    cdef uint32_t* rand_r_state = &rand_r_state_seed
    cdef bint center = False
    cdef bint no_sample_weights = sample_weight is None

    if alpha == 0:
        # No screeing without L1-penalty.
        do_screening = False

    if do_screening:
        active_set = np.empty(n_features, dtype=np.uint32)  # map [:n_active] -> j
        excluded_set = np.empty(n_features, dtype=np.uint8)

    if no_sample_weights:
        yw = y
        R = y.copy()
    else:
        yw = np.multiply(sample_weight, y)
        R = yw.copy()

    with nogil:
        # center = (X_mean != 0).any()
        for j in range(n_features):
            if X_mean[j]:
                center = True
                break

        # R = y - np.dot(X, w)
        for j in range(n_features):
            X_mean_j = X_mean[j]
            endptr = X_indptr[j + 1]
            normalize_sum = 0.0
            w_j = w[j]

            if no_sample_weights:
                for i in range(startptr, endptr):
                    normalize_sum += (X_data[i] - X_mean_j) ** 2
                    R[X_indices[i]] -= X_data[i] * w_j
                norm2_cols_X[j] = normalize_sum + \
                    (n_samples - endptr + startptr) * X_mean_j ** 2
                if center:
                    for i in range(n_samples):
                        R[i] += X_mean_j * w_j
                        R_sum += R[i]
            else:
                # R = sw * (y - np.dot(X, w))
                for i in range(startptr, endptr):
                    tmp = sample_weight[X_indices[i]]
                    # second term will be subtracted by loop over range(n_samples)
                    normalize_sum += (tmp * (X_data[i] - X_mean_j) ** 2
                                      - tmp * X_mean_j ** 2)
                    R[X_indices[i]] -= tmp * X_data[i] * w_j
                if center:
                    for i in range(n_samples):
                        normalize_sum += sample_weight[i] * X_mean_j ** 2
                        R[i] += sample_weight[i] * X_mean_j * w_j
                        R_sum += R[i]
                norm2_cols_X[j] = normalize_sum
            startptr = endptr

        # Note: No need to update R_sum from here on because the update terms cancel
        # each other: w_j * np.sum(X[:,j] - X_mean[j]) = 0. R_sum is only ever
        # needed and calculated if X_mean is provided.

        # tol *= np.dot(y, y)
        # with sample weights: tol *= y @ (sw * y)
        tol *= _dot(n_samples, &y[0], 1, &yw[0], 1)

        # Check convergence before entering the main loop.
        gap, dual_norm_XtA = gap_enet_sparse(
            n_samples,
            n_features,
            w,
            alpha,
            beta,
            X_data,
            X_indices,
            X_indptr,
            y,
            sample_weight,
            no_sample_weights,
            X_mean,
            center,
            R,
            R_sum,
            XtA,
            positive,
            False,
        )
        if early_stopping and gap <= tol:
            with gil:
                return np.asarray(w), gap, tol, 0

        # Gap Safe Screening Rules, see https://arxiv.org/abs/1802.07481, Eq. 11
        if do_screening:
            n_active = 0
            for j in range(n_features):
                if norm2_cols_X[j] == 0:
                    w[j] = 0
                    excluded_set[j] = 1
                    continue
                Xj_theta = XtA[j] / fmax(alpha, dual_norm_XtA)  # X[:,j] @ dual_theta
                d_j = (1 - fabs(Xj_theta)) / sqrt(norm2_cols_X[j] + beta)
                if d_j <= sqrt(2 * gap) / alpha:
                    # include feature j
                    active_set[n_active] = j
                    excluded_set[j] = 0
                    n_active += 1
                else:
                    if w[j] != 0:
                        # R += w[j] * X[:,j]
                        R_plus_wj_Xj(
                            n_samples,
                            R,
                            X_data,
                            X_indices,
                            X_indptr,
                            X_mean,
                            center,
                            sample_weight,
                            no_sample_weights,
                            w[j],
                            j,
                        )
                        w[j] = 0
                    excluded_set[j] = 1

        for n_iter in range(max_iter):
            w_max = 0.0
            d_w_max = 0.0
            for f_iter in range(n_active):  # Loop over coordinates
                if random:
                    j = rand_int(n_active, rand_r_state)
                else:
                    j = f_iter

                if do_screening:
                    j = active_set[j]

                if norm2_cols_X[j] == 0.0:
                    continue

                startptr = X_indptr[j]
                endptr = X_indptr[j + 1]
                w_j = w[j]  # Store previous value
                X_mean_j = X_mean[j]

                # tmp = X[:,j] @ (R + w_j * X[:,j])
                tmp = 0.0
                for i in range(startptr, endptr):
                    tmp += R[X_indices[i]] * X_data[i]
                tmp += w_j * norm2_cols_X[j]

                if center:
                    tmp -= R_sum * X_mean_j

                if positive and tmp < 0.0:
                    w[j] = 0.0
                else:
                    w[j] = fsign(tmp) * fmax(fabs(tmp) - alpha, 0) \
                            / (norm2_cols_X[j] + beta)

                if w[j] != w_j:
                    # R -=  (w[j] - w_j) * X[:,j] # Update residual
                    R_plus_wj_Xj(
                        n_samples,
                        R,
                        X_data,
                        X_indices,
                        X_indptr,
                        X_mean,
                        center,
                        sample_weight,
                        no_sample_weights,
                        w_j - w[j],
                        j,
                    )

                # update the maximum absolute coefficient update
                d_w_j = fabs(w[j] - w_j)
                d_w_max = fmax(d_w_max, d_w_j)

                w_max = fmax(w_max, fabs(w[j]))

            if (
                w_max == 0.0
                or d_w_max / w_max <= d_w_tol
                or n_iter == max_iter - 1
                or n_active <= 1  # We have an analytical exact solution.
            ):
                # the biggest coordinate update of this iteration was smaller than
                # the tolerance: check the duality gap as ultimate stopping
                # criterion
                gap, dual_norm_XtA = gap_enet_sparse(
                    n_samples,
                    n_features,
                    w,
                    alpha,
                    beta,
                    X_data,
                    X_indices,
                    X_indptr,
                    y,
                    sample_weight,
                    no_sample_weights,
                    X_mean,
                    center,
                    R,
                    R_sum,
                    XtA,
                    positive,
                    True,
                )

                if gap <= tol:
                    # return if we reached desired tolerance
                    break

                # Gap Safe Screening Rules, see https://arxiv.org/abs/1802.07481, Eq. 11
                if do_screening:
                    n_active = 0
                    for j in range(n_features):
                        if excluded_set[j]:
                            continue
                        Xj_theta = XtA[j] / fmax(alpha, dual_norm_XtA)  # X @ dual_theta
                        d_j = (1 - fabs(Xj_theta)) / sqrt(norm2_cols_X[j] + beta)
                        if d_j <= sqrt(2 * gap) / alpha:
                            # include feature j
                            active_set[n_active] = j
                            excluded_set[j] = 0
                            n_active += 1
                        else:
                            if w[j] != 0:
                                # R += w[j] * X[:,j]
                                R_plus_wj_Xj(
                                    n_samples,
                                    R,
                                    X_data,
                                    X_indices,
                                    X_indptr,
                                    X_mean,
                                    center,
                                    sample_weight,
                                    no_sample_weights,
                                    w[j],
                                    j,
                                )
                                w[j] = 0
                            excluded_set[j] = 1

        else:
            # for/else, runs if for doesn't end with a `break`
            with gil:
                message = (
                    message_conv +
                    f" Duality gap: {gap:.6e}, tolerance: {tol:.3e}"
                )
                if alpha < np.finfo(np.float64).eps:
                    message += "\n" + message_ridge
                warnings.warn(message, ConvergenceWarning)

    return np.asarray(w), gap, tol, n_iter + 1


cdef (floating, floating) gap_enet_gram(
    int n_features,
    const floating[::1] w,
    floating alpha,  # L1 penalty
    floating beta,  # L2 penalty
    const floating[::1] Qw,
    const floating[::1] q,
    const floating y_norm2,
    floating[::1] XtA,  # XtA = X.T @ R - beta * w is calculated inplace
    bint positive,
    bint gap_smaller_eps,
) noexcept nogil:
    """Compute dual gap for use in enet_coordinate_descent.

    alpha > 0:            formulation A of the duality gap
    alpha = 0 & beta > 0: formulation B of the duality gap
    alpha = beta = 0:     OLS first order condition (=gradient)
    """
    cdef floating gap, primal, dual
    cdef floating dual_norm_XtA
    cdef floating R_norm2
    cdef floating Ry
    cdef floating w_l1_norm
    cdef floating w_l2_norm2 = 0.0
    cdef floating q_dot_w
    cdef floating wQw
    cdef unsigned int j

    if floating is float:
        eps = FLT_EPSILON
    else:
        eps = DBL_EPSILON

    # w_l2_norm2 = w @ w
    if beta > 0:
        w_l2_norm2 = _dot(n_features, &w[0], 1, &w[0], 1)
    # q_dot_w = w @ q
    q_dot_w = _dot(n_features, &w[0], 1, &q[0], 1)
    # wQw = w @ Q @ w
    wQw = _dot(n_features, &w[0], 1, &Qw[0], 1)
    # R_norm2 = R @ R, residual R = y - Xw
    R_norm2 = y_norm2 + wQw - 2.0 * q_dot_w
    # Ry = R @ y
    if not (alpha == 0 and beta == 0):
        # Note that R'y = (y - Xw)' y = ||y||_2^2 - w'X'y = y_norm2 - q_dot_w
        Ry = y_norm2 - q_dot_w

    if alpha == 0:
        # XtA = X'R
        for j in range(n_features):
            XtA[j] = q[j] - Qw[j]
        # ||X'R||_2^2
        dual_norm_XtA = _dot(n_features, &XtA[0], 1, &XtA[0], 1)
        if beta == 0:
            # This is OLS, no dual gap available. Resort to first order condition
            #     X'R = 0
            #     gap = ||X'R||_2^2
            # Compare with stopping criterion of LSQR.
            gap = dual_norm_XtA
            return gap, dual_norm_XtA
        # This is Ridge regression, we use formulation B for the dual gap.
        primal = 0.5 * (R_norm2 + beta * w_l2_norm2)
        dual = -0.5 * R_norm2 + Ry - 1 / (2 * beta) * dual_norm_XtA
        gap = primal - dual
        if gap_smaller_eps and abs(gap) <= 2 * eps * primal:
            gap = 0.0
        return gap, dual_norm_XtA

    # XtA = X.T @ R - beta * w = X.T @ y - X.T @ X @ w - beta * w
    for j in range(n_features):
        XtA[j] = q[j] - Qw[j] - beta * w[j]

    # dual_norm_XtA
    if positive:
        dual_norm_XtA = max(n_features, &XtA[0])
    else:
        dual_norm_XtA = abs_max(n_features, &XtA[0])

    # w_l1_norm = np.sum(np.abs(w))
    w_l1_norm = _asum(n_features, &w[0], 1)

    gap = dual_gap_formulation_A(
        alpha=alpha,
        beta=beta,
        w_l1_norm=w_l1_norm,
        w_l2_norm2=w_l2_norm2,
        R_norm2=R_norm2,
        Ry=Ry,
        dual_norm_XtA=dual_norm_XtA,
        gap_smaller_eps=gap_smaller_eps,
    )
    return gap, dual_norm_XtA


def enet_coordinate_descent_gram(
    floating[::1] w,
    floating alpha,
    floating beta,
    const floating[:, ::1] Q,
    const floating[::1] q,
    const floating[:] y,
    unsigned int max_iter,
    floating tol,
    object rng,
    bint random=0,
    bint positive=0,
    bint do_screening=1,
    bint early_stopping=1,
):
    """Cython version of the coordinate descent algorithm
        for Elastic-Net regression

        We minimize

        (1/2) * w^T Q w - q^T w + alpha norm(w, 1) + (beta/2) * norm(w, 2)^2
        +1/2 * y^T y

        which amount to the Elastic-Net problem when:
        Q = X^T X (Gram matrix)
        q = X^T y

    Returns
    -------
    w : ndarray of shape (n_features,)
        ElasticNet coefficients.
    gap : float
        Achieved dual gap.
    tol : float
        Equals input `tol` times `np.dot(y, y)`. The tolerance used for the dual gap.
    n_iter : int
        Number of coordinate descent iterations.
    """

    if floating is float:
        dtype = np.float32
    else:
        dtype = np.float64

    # get the data information into easy vars
    cdef unsigned int n_features = Q.shape[0]

    # initial value "Q w" which will be kept of up to date in the iterations
    cdef floating[::1] Qw = np.dot(Q, w)
    cdef floating[::1] XtA = np.zeros(n_features, dtype=dtype)
    cdef floating y_norm2 = np.dot(y, y)

    cdef floating d_j
    cdef floating radius
    cdef floating Xj_theta
    cdef floating tmp
    cdef floating w_j
    cdef floating d_w_max
    cdef floating w_max
    cdef floating d_w_j
    cdef floating gap = tol + 1.0
    cdef floating d_w_tol = tol
    cdef floating dual_norm_XtA
    cdef unsigned int n_active = n_features
    cdef uint32_t[::1] active_set
    # TODO: use binset insteaf of array of bools
    cdef uint8_t[::1] excluded_set
    cdef unsigned int j
    cdef unsigned int n_iter = 0
    cdef unsigned int f_iter
    cdef uint32_t rand_r_state_seed = rng.randint(0, RAND_R_MAX)
    cdef uint32_t* rand_r_state = &rand_r_state_seed

    if alpha == 0:
        # No screeing without L1-penalty.
        do_screening = False

    if do_screening:
        active_set = np.empty(n_features, dtype=np.uint32)  # map [:n_active] -> j
        excluded_set = np.empty(n_features, dtype=np.uint8)

    with nogil:
        tol *= y_norm2

        # Check convergence before entering the main loop.
        gap, dual_norm_XtA = gap_enet_gram(
            n_features, w, alpha, beta, Qw, q, y_norm2, XtA, positive, False
        )
        if early_stopping and 0 <= gap <= tol:
            # Only if gap >=0 as singular Q may cause dubious values of gap.
            with gil:
                return np.asarray(w), gap, tol, 0

        # Gap Safe Screening Rules, see https://arxiv.org/abs/1802.07481, Eq. 11
        if do_screening:
            # Due to floating point issues, gap might be negative.
            radius = sqrt(2 * fabs(gap)) / alpha
            n_active = 0
            for j in range(n_features):
                if Q[j, j] == 0:
                    w[j] = 0
                    excluded_set[j] = 1
                    continue
                Xj_theta = XtA[j] / fmax(alpha, dual_norm_XtA)  # X[:,j] @ dual_theta
                d_j = (1 - fabs(Xj_theta)) / sqrt(Q[j, j] + beta)
                if d_j <= radius:
                    # include feature j
                    active_set[n_active] = j
                    excluded_set[j] = 0
                    n_active += 1
                else:
                    if w[j] != 0:
                        # Qw -= w[j] * Q[j]  # Update Qw = Q @ w
                        _axpy(n_features, -w[j], &Q[j, 0], 1, &Qw[0], 1)
                        w[j] = 0
                    excluded_set[j] = 1

        for n_iter in range(max_iter):
            w_max = 0.0
            d_w_max = 0.0
            for f_iter in range(n_active):  # Loop over coordinates
                if random:
                    j = rand_int(n_active, rand_r_state)
                else:
                    j = f_iter

                if do_screening:
                    j = active_set[j]

                if Q[j, j] == 0.0:
                    continue

                w_j = w[j]  # Store previous value

                # if Q = X.T @ X then tmp = X[:,j] @ (y - X @ w + X[:, j] * w_j)
                tmp = q[j] - Qw[j] + w_j * Q[j, j]

                if positive and tmp < 0:
                    w[j] = 0.0
                else:
                    w[j] = fsign(tmp) * fmax(fabs(tmp) - alpha, 0) \
                        / (Q[j, j] + beta)

                if w[j] != w_j:
                    # Qw += (w[j] - w_j) * Q[j]  # Update Qw = Q @ w
                    _axpy(n_features, w[j] - w_j, &Q[j, 0], 1, &Qw[0], 1)

                # update the maximum absolute coefficient update
                d_w_j = fabs(w[j] - w_j)
                if d_w_j > d_w_max:
                    d_w_max = d_w_j

                if fabs(w[j]) > w_max:
                    w_max = fabs(w[j])

            if (
                w_max == 0.0
                or d_w_max / w_max <= d_w_tol
                or n_iter == max_iter - 1
                or n_active <= 1  # We have an analytical exact solution.
            ):
                # the biggest coordinate update of this iteration was smaller than
                # the tolerance: check the duality gap as ultimate stopping
                # criterion
                gap, dual_norm_XtA = gap_enet_gram(
                    n_features, w, alpha, beta, Qw, q, y_norm2, XtA, positive, True
                )

                if gap <= tol:
                    # return if we reached desired tolerance
                    break

                # Gap Safe Screening Rules, see https://arxiv.org/abs/1802.07481, Eq. 11
                if do_screening:
                    # Due to floating point issues, gap might be negative.
                    radius = sqrt(2 * fabs(gap)) / alpha
                    n_active = 0
                    for j in range(n_features):
                        if excluded_set[j]:
                            continue
                        Xj_theta = XtA[j] / fmax(alpha, dual_norm_XtA)  # X @ dual_theta
                        d_j = (1 - fabs(Xj_theta)) / sqrt(Q[j, j] + beta)
                        if d_j <= radius:
                            # include feature j
                            active_set[n_active] = j
                            excluded_set[j] = 0
                            n_active += 1
                        else:
                            if w[j] != 0:
                                # Qw -= w[j] * Q[j]  # Update Qw = Q @ w
                                _axpy(n_features, -w[j], &Q[j, 0], 1, &Qw[0], 1)
                                w[j] = 0
                            excluded_set[j] = 1

        else:
            # for/else, runs if for doesn't end with a `break`
            with gil:
                message = (
                    message_conv +
                    f" Duality gap: {gap:.6e}, tolerance: {tol:.3e}"
                )
                if alpha < np.finfo(np.float64).eps:
                    message += "\n" + message_ridge
                warnings.warn(message, ConvergenceWarning)

    return np.asarray(w), gap, tol, n_iter + 1


cdef (floating, floating) gap_enet_multi_task(
    int n_samples,
    int n_features,
    int n_tasks,
    const floating[::1, :] W,  # in
    floating alpha,
    floating beta,
    const floating[::1, :] X,  # in
    const floating[::1, :] Y,  # in
    const floating[::1, :] R,  # in
    floating[:, ::1] XtA,  # out
    floating[::1] XtA_row_norms,  # out
    bint gap_smaller_eps,
) noexcept nogil:
    """Compute dual gap for use in enet_coordinate_descent_multi_task.

    Parameters
    ----------
    W : memoryview of shape (n_tasks, n_features)
    X : memoryview of shape (n_samples, n_features)
    Y : memoryview of shape (n_samples, n_tasks)
    R : memoryview of shape (n_samples, n_tasks)
        Current residuals = Y - X @ W.T
    XtA : memoryview of shape (n_features, n_tasks)
        Inplace calculated as XtA = X.T @ R - beta * W.T
    XtA_row_norms : memoryview of shape n_features
        Inplace calculated as np.sqrt(np.sum(XtA ** 2, axis=1))
    """
    cdef floating gap, primal, dual
    cdef floating dual_norm_XtA
    cdef floating R_norm2
    cdef floating Ry
    cdef floating w_l21_norm
    cdef floating w_l2_norm2 = 0.0
    cdef unsigned int t, j

    if floating is float:
        eps = FLT_EPSILON
    else:
        eps = DBL_EPSILON

    # w_l2_norm2 = linalg.norm(W, ord="fro") ** 2
    if beta > 0:
        w_l2_norm2 = _dot(n_features * n_tasks, &W[0, 0], 1, &W[0, 0], 1)
    # R_norm2 = linalg.norm(R, ord="fro") ** 2
    R_norm2 = _dot(n_samples * n_tasks, &R[0, 0], 1, &R[0, 0], 1)
    # Ry = np.sum(R * Y)
    if not (alpha == 0 and beta == 0):
        Ry = _dot(n_samples * n_tasks, &R[0, 0], 1, &Y[0, 0], 1)

    if alpha == 0:
        # XtA = X.T @ R
        for j in range(n_features):
            for t in range(n_tasks):
                XtA[j, t] = _dot(n_samples, &X[0, j], 1, &R[0, t], 1)
        # ||X'R||_2^2
        dual_norm_XtA = _dot(n_features * n_tasks, &XtA[0, 0], 1, &XtA[0, 0], 1)
        if beta == 0:
            # This is OLS, no dual gap available. Resort to first order condition
            #     X'R = 0
            #     gap = ||X'R||_2^2
            # Compare with stopping criterion of LSQR.
            gap = dual_norm_XtA
            return gap, dual_norm_XtA
        # This is Ridge regression, we use formulation B for the dual gap.
        primal = 0.5 * (R_norm2 + beta * w_l2_norm2)
        dual = -0.5 * R_norm2 + Ry - 1 / (2 * beta) * dual_norm_XtA
        gap = primal - dual
        if gap_smaller_eps and abs(gap) <= 2 * eps * primal:
            gap = 0.0
        return gap, dual_norm_XtA

    # XtA = X.T @ R - beta * W.T
    for j in range(n_features):
        for t in range(n_tasks):
            XtA[j, t] = _dot(n_samples, &X[0, j], 1, &R[0, t], 1) - beta * W[t, j]

    # dual_norm_XtA = np.max(np.sqrt(np.sum(XtA ** 2, axis=1)))
    dual_norm_XtA = 0.0
    for j in range(n_features):
        # np.sqrt(np.sum(XtA ** 2, axis=1))
        XtA_row_norms[j] = _nrm2(n_tasks, &XtA[j, 0], 1)
        if XtA_row_norms[j] > dual_norm_XtA:
            dual_norm_XtA = XtA_row_norms[j]

    # w_l21_norm = np.sqrt(np.sum(W ** 2, axis=0)).sum()
    w_l21_norm = 0.0
    for ii in range(n_features):
        w_l21_norm += _nrm2(n_tasks, &W[0, ii], 1)

    gap = dual_gap_formulation_A(
        alpha=alpha,
        beta=beta,
        w_l1_norm=w_l21_norm,
        w_l2_norm2=w_l2_norm2,
        R_norm2=R_norm2,
        Ry=Ry,
        dual_norm_XtA=dual_norm_XtA,
        gap_smaller_eps=gap_smaller_eps,
    )
    return gap, dual_norm_XtA


def enet_coordinate_descent_multi_task(
    floating[::1, :] W,
    floating alpha,
    floating beta,
    const floating[::1, :] X,
    const floating[::1, :] Y,
    unsigned int max_iter,
    floating tol,
    object rng,
    bint random=0,
    bint do_screening=1,
    bint early_stopping=1,
):
    """Cython version of the coordinate descent algorithm
        for Elastic-Net multi-task regression

        We minimize

        0.5 * norm(Y - X W.T, 2)^2 + alpha * ||W.T||_21 + 0.5 * beta * norm(W.T, 2)^2

    The algorithm follows
    Noah Simon, Jerome Friedman, Trevor Hastie. 2013.
    A Blockwise Descent Algorithm for Group-penalized Multiresponse and Multinomial
    Regression
    https://doi.org/10.48550/arXiv.1311.6529

    Returns
    -------
    W : ndarray of shape (n_tasks, n_features)
        ElasticNet coefficients.
    gap : float
        Achieved dual gap.
    tol : float
        Equals input `tol` times `np.dot(y, y)`. The tolerance used for the dual gap.
    n_iter : int
        Number of coordinate descent iterations.
    """

    if floating is float:
        dtype = np.float32
    else:
        dtype = np.float64

    # get the data information into easy vars
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]
    cdef unsigned int n_tasks = Y.shape[1]

    # compute squared norms of the columns of X
    # same as norm2_cols_X = np.square(X).sum(axis=0)
    cdef floating[::1] norm2_cols_X = np.einsum(
        "ij,ij->j", X, X, dtype=dtype, order="C"
    )

    # initial value of the residuals
    cdef floating[::1, :] R = np.empty((n_samples, n_tasks), dtype=dtype, order='F')
    cdef floating[:, ::1] XtA = np.empty((n_features, n_tasks), dtype=dtype)
    cdef floating[::1] XtA_row_norms = np.empty(n_features, dtype=dtype)

    cdef floating d_j
    cdef floating Xj_theta
    cdef floating[::1] tmp = np.empty(n_tasks, dtype=dtype)
    cdef floating[::1] w_j = np.empty(n_tasks, dtype=dtype)
    cdef floating d_w_max
    cdef floating w_max
    cdef floating d_w_j
    cdef floating nn
    cdef floating W_j_abs_max
    cdef floating gap = tol + 1.0
    cdef floating d_w_tol = tol
    cdef floating dual_norm_XtA
    cdef unsigned int n_active = n_features
    cdef uint32_t[::1] active_set
    # TODO: use binset instead of array of bools
    cdef uint8_t[::1] excluded_set
    cdef unsigned int j
    cdef unsigned int t
    cdef unsigned int n_iter = 0
    cdef unsigned int f_iter
    cdef uint32_t rand_r_state_seed = rng.randint(0, RAND_R_MAX)
    cdef uint32_t* rand_r_state = &rand_r_state_seed

    if alpha == 0:
        # No screeing without L1-penalty.
        do_screening = False

    if do_screening:
        active_set = np.empty(n_features, dtype=np.uint32)  # map [:n_active] -> j
        excluded_set = np.empty(n_features, dtype=np.uint8)

    with nogil:
        # R = Y - X @ W.T
        _copy(n_samples * n_tasks, &Y[0, 0], 1, &R[0, 0], 1)
        for j in range(n_features):
            for t in range(n_tasks):
                if W[t, j] != 0:
                    _axpy(n_samples, -W[t, j], &X[0, j], 1, &R[0, t], 1)

        # tol = tol * linalg.norm(Y, ord='fro') ** 2
        tol = tol * _nrm2(n_samples * n_tasks, &Y[0, 0], 1) ** 2

        # Check convergence before entering the main loop.
        gap, dual_norm_XtA = gap_enet_multi_task(
            n_samples, n_features, n_tasks, W, alpha, beta, X, Y, R, XtA, XtA_row_norms, False
        )
        if early_stopping and gap <= tol:
            with gil:
                return np.asarray(W), gap, tol, 0

        # Gap Safe Screening Rules for multi-task Lasso, see
        # https://arxiv.org/abs/1703.07285 Eq 2.2. (also arxiv:1506.03736)
        if do_screening:
            n_active = 0
            for j in range(n_features):
                if norm2_cols_X[j] == 0:
                    for t in range(n_tasks):
                        W[t, j] = 0
                    excluded_set[j] = 1
                    continue
                # Xj_theta = ||X[:,j] @ dual_theta||_2
                Xj_theta = XtA_row_norms[j] / fmax(alpha, dual_norm_XtA)
                d_j = (1 - Xj_theta) / sqrt(norm2_cols_X[j] + beta)
                if d_j <= sqrt(2 * gap) / alpha:
                    # include feature j
                    active_set[n_active] = j
                    excluded_set[j] = 0
                    n_active += 1
                else:
                    # R += W[:, 1] * X[:, 1][:, None]
                    for t in range(n_tasks):
                        if W[t, j] != 0:
                            _axpy(n_samples, W[t, j], &X[0, j], 1, &R[0, t], 1)
                            W[t, j] = 0
                    excluded_set[j] = 1

        for n_iter in range(max_iter):
            w_max = 0.0
            d_w_max = 0.0
            for f_iter in range(n_active):  # Loop over coordinates
                if random:
                    j = rand_int(n_active, rand_r_state)
                else:
                    j = f_iter

                if do_screening:
                    j = active_set[j]

                if norm2_cols_X[j] == 0.0:
                    continue

                # w_j = W[:, j] # Store previous value
                _copy(n_tasks, &W[0, j], 1, &w_j[0], 1)

                # tmp = X[:, j] @ (R + w_j * X[:,j][:, None])
                # first part: X[:, j] @ R
                #   Using BLAS Level 2:
                #   _gemv(RowMajor, Trans, n_samples, n_tasks, 1.0, &R[0, 0],
                #         n_tasks, &X[0, j], 1, 0.0, &tmp[0], 1)
                # second part: (X[:, j] @ X[:,j]) * w_j = norm2_cols * w_j
                #   Using BLAS Level 1:
                #   _axpy(n_tasks, norm2_cols[j], &w_j[0], 1, &tmp[0], 1)
                # Using BLAS Level 1 (faster for small vectors like here):
                for t in range(n_tasks):
                    tmp[t] = _dot(n_samples, &X[0, j], 1, &R[0, t], 1)
                    # As we have the loop already, we use it to replace the second BLAS
                    # Level 1, i.e., _axpy, too.
                    tmp[t] += w_j[t] * norm2_cols_X[j]

                # nn = sqrt(np.sum(tmp ** 2))
                nn = _nrm2(n_tasks, &tmp[0], 1)

                # W[:, j] = tmp * fmax(1. - alpha / nn, 0) / (norm2_cols_X[j] + beta)
                _copy(n_tasks, &tmp[0], 1, &W[0, j], 1)
                _scal(n_tasks, fmax(1. - alpha / nn, 0) / (norm2_cols_X[j] + beta),
                      &W[0, j], 1)

                # Update residual
                # Using numpy:
                #   R -= (W[:, j] - w_j) * X[:, j][:, None]
                # Using BLAS Level 1 and 2:
                #   _axpy(n_tasks, -1.0, &W[0, j], 1, &w_j[0], 1)
                #   _ger(RowMajor, n_samples, n_tasks, 1.0,
                #        &X[0, j], 1, &w_j, 1,
                #        &R[0, 0], n_tasks)
                # Using BLAS Level 1 (faster for small vectors like here):
                for t in range(n_tasks):
                    if W[t, j] != w_j[t]:
                        _axpy(n_samples, w_j[t] - W[t, j], &X[0, j], 1, &R[0, t], 1)

                # update the maximum absolute coefficient update
                d_w_j = diff_abs_max(n_tasks, &W[0, j], &w_j[0])

                if d_w_j > d_w_max:
                    d_w_max = d_w_j

                W_j_abs_max = abs_max(n_tasks, &W[0, j])
                if W_j_abs_max > w_max:
                    w_max = W_j_abs_max

            if (
                w_max == 0.0
                or d_w_max / w_max <= d_w_tol
                or n_iter == max_iter - 1
                or n_active <= 1  # We have an analytical exact solution.
            ):
                # the biggest coordinate update of this iteration was smaller than
                # the tolerance: check the duality gap as ultimate stopping
                # criterion
                gap, dual_norm_XtA = gap_enet_multi_task(
                    n_samples, n_features, n_tasks, W, alpha, beta, X, Y, R, XtA, XtA_row_norms, True
                )
                if gap <= tol:
                    # return if we reached desired tolerance
                    break

                # Gap Safe Screening Rules for multi-task Lasso, see
                # https://arxiv.org/abs/1703.07285 Eq 2.2. (also arxiv:1506.03736)
                if do_screening:
                    n_active = 0
                    for j in range(n_features):
                        if excluded_set[j]:
                            continue
                        # Xj_theta = ||X[:,j] @ dual_theta||_2
                        Xj_theta = XtA_row_norms[j] / fmax(alpha, dual_norm_XtA)
                        d_j = (1 - Xj_theta) / sqrt(norm2_cols_X[j] + beta)
                        if d_j <= sqrt(2 * gap) / alpha:
                            # include feature j
                            active_set[n_active] = j
                            excluded_set[j] = 0
                            n_active += 1
                        else:
                            # R += W[:, 1] * X[:, 1][:, None]
                            for t in range(n_tasks):
                                if W[t, j] != 0:
                                    _axpy(n_samples, W[t, j], &X[0, j], 1, &R[0, t], 1)
                                    W[t, j] = 0
                            excluded_set[j] = 1

        else:
            # for/else, runs if for doesn't end with a `break`
            with gil:
                message = (
                    message_conv +
                    f" Duality gap: {gap:.6e}, tolerance: {tol:.3e}"
                )
                if alpha < np.finfo(np.float64).eps:
                    message += "\n" + message_ridge
                warnings.warn(message, ConvergenceWarning)

    return np.asarray(W), gap, tol, n_iter + 1


cdef void inverse_L_sqrt_D_matmul(
    const floating[::1, :] p,  # in
    const floating[::1, :] q_inv,  # in
    const floating[::1, :] sqrt_d,  # in
    floating[::1, :] x,  # out
    floating[::1] buffer,  # temporary memory buffer
) noexcept nogil:
    """Compute 1 / sqrt(D) L^(-1) x from the multinomial LDL' decomposition.

    Same as Multinomial_LDL_Decomposition.inverse_L_sqrt_D_matmul

    L^(-1) is again lower triangular and given by:
        L^(-1)_ij = f[:, i] = p[:, i] / q[:, i-1]

    For n_classes = 4, L^(-1) looks like:: text

        L^(-1) = (1   0  0  0)
                 (f1  1  0  0)
                 (f2  f2 1  0)
                 (f3  f3 f3 1)

    Note that (L sqrt(D))^(-1) = 1/sqrt(D) L^(-1).

    Parameters
    ----------
    p : memoryview of shape (n_samples, n_classes)
        Probabilities.

    q_inv : memoryview of shape (n_samples, n_classes)
        Inverse of `q = 1 - np.cumsum(self.p, axis=1)`.

    sqrt_d : memoryview of shape (n_samples, n_classes)
        Square root of diagonal D, with D_ii = p_i * q_i / q_{i-1}.

    x : memoryview of shape (n_samples, n_classes)
        This array is overwritten with the result.

    buffer : memoryview of shape (n_samples)
        Used as temporary memory buffer.
    """
    cdef unsigned int n_samples = p.shape[0]
    cdef unsigned int n_classes = p.shape[1]
    cdef unsigned int i, k, t
    cdef bint mask
    cdef floating fj
    cdef floating[::1] x_sum = buffer

    # x_sum = np.sum(x[:, :-1], axis=1)  # precomputation
    _copy(n_samples, &x[0, 0], 1, &x_sum[0], 1)
    for k in range(1, n_classes - 1):
        for t in range(n_samples):
            x_sum[t] += x[t, k]

    for i in range(n_classes - 1, 0, -1):  # row i, here i > 0.
        # fj = p[:, i] * q_inv[:, i - 1]
        # Using precomputation
        # x[:, i] += fj * x_sum
        # if i > 1:
        #     x_sum -= x[:, i - 1]
        for t in range(n_samples):
            fj = p[t, i] * q_inv[t, i - 1]
            x[t, i] += fj * x_sum[t]
        if i > 1:
            for t in range(n_samples):
                x_sum[t] -= x[t, i - 1]

    # x[:, :-1] /= sqrt_d[:, :-1]
    for k in range(n_classes - 1):
        for t in range(n_samples):
            mask = (sqrt_d[t, k] == 0)
            x[t, k] *= (not mask) / (sqrt_d[t, k] + mask)
    # Important Note:
    # Strictly speaking, the inverse of D does not exist.
    # We use 0 as the inverse of 0 and just set:
    # x[:, -1] = 0
    memset(&x[0, n_classes - 1], 0, n_samples * sizeof(floating))


cdef (floating, floating) gap_enet_multinomial(
    const floating[::1, :] w,
    floating alpha,  # L1 penalty
    floating beta,  # L2 penalty
    const floating[::1, :] X,
    bint X_is_sparse,
    floating[::1] X_data,
    int32_t[::1] X_indices,
    int32_t[::1] X_indptr,
    const floating[::1, :] y,
    const floating[::1, :] R,  # current residual b - A @ w
    const floating[::1, :] LD_R,  # LD_R = LDL.L_sqrt_D_matmul(R.copy()) * sqrt_sw[:, None]
    const floating[::1, :] H00_pinv_H0,
    bint fit_intercept,
    floating[::1, :] XtA,  # XtA = X.T @ R - beta * w is calculated inplace
    bint gap_smaller_eps,
) noexcept nogil:
    """Compute dual gap for use in enet_coordinate_descent.

    Note that X must always be replaced by A = sqrt(D) L' X
    """
    cdef floating gap, primal, dual
    cdef floating scale  # Scaling factor to achieve dual feasible point.
    cdef floating dual_norm_XtA
    cdef floating R_norm2
    cdef floating Ry
    cdef floating w_l1_norm = 0.0
    cdef floating w_l2_norm2 = 0.0
    cdef unsigned int n_classes = w.shape[0]
    cdef unsigned int n_features = w.shape[1]
    cdef unsigned int n_samples = R.shape[0]
    cdef floating* LD_R_sum0
    cdef unsigned int i, j, k
    cdef int32_t i_ind

    if floating is float:
        eps = FLT_EPSILON
    else:
        eps = DBL_EPSILON

    if beta > 0:
        # l2_norm2 = w @ w
        w_l2_norm2 = _dot(n_classes * n_features, &w[0, 0], 1, &w[0, 0], 1)
    if alpha > 0:
        # w_l1_norm = np.sum(np.abs(w))
        w_l1_norm = _asum(n_classes * n_features, &w[0, 0], 1)
    # R_norm2 = np.linalg.norm(R, ord="fro") ** 2
    R_norm2 = _dot(n_classes * n_samples, &R[0, 0], 1, &R[0, 0], 1)
    if not (alpha == 0 and beta == 0):
        # Ry = R.ravel(order="F") @ y.ravel(order="F")
        Ry = _dot(n_classes * n_samples, &R[0, 0], 1, &y[0, 0], 1)

    primal = 0.5 * R_norm2 + alpha * w_l1_norm + 0.5 * beta * w_l2_norm2

    # XtA = X'R -> A'R = X' L sqrt(D) R
    # XtA[:, :] = (X.T @ LD_R).T
    if not X_is_sparse:
        for k in range(n_classes):
            # XtA[:, k] = X.T @ LD_R[:, k]
            _gemv(
                ColMajor, Trans, n_samples, n_features, 1.0, &X[0, 0],
                n_samples, &LD_R[0, k], 1, 0.0, &XtA[k, 0], n_classes,
            )
    else:
        # XtA[:, :] = 0
        memset(&XtA[0, 0], 0, n_classes * n_features * sizeof(floating))
        # XtA[:, :] = X.T @ LD_R
        for j in range(n_features):
            for i_ind in range(X_indptr[j], X_indptr[j + 1]):
                for k in range(n_classes):
                    XtA[k, j] += X_data[i_ind] * LD_R[X_indices[i_ind], k]

    if fit_intercept:
        # temporary memory
        LD_R_sum0 = <floating *> malloc(sizeof(floating) * (n_classes - 1))
        # XtA -= (H00_pinv_H0.T @ LD_R[:, :-1].sum(axis=0)).reshape(
        #     n_classes, -1, order="F"
        # )
        for k in  range(n_classes - 1):
            LD_R_sum0[k] = 0
            for i in range(n_samples):
                LD_R_sum0[k] += LD_R[i, k]
        _gemv(
            ColMajor, Trans, n_classes - 1, n_classes * n_features,
            -1.0, &H00_pinv_H0[0, 0],
            n_classes - 1, &LD_R_sum0[0], 1, 1.0, &XtA[0, 0], 1,
        )
        free(LD_R_sum0)

    if alpha == 0:
        # ||X'R||_2^2
        # dual_norm_XtA = np.linalg.norm(XtA, ord="fro") ** 2
        dual_norm_XtA = _dot(n_classes * n_features, &XtA[0, 0], 1, &XtA[0, 0], 1)
        if beta == 0:
            # This is OLS, no dual gap available. Resort to first order condition
            #     X'R = 0
            #     gap = ||X'R||_2^2
            # Compare with stopping criterion of LSQR.
            gap = dual_norm_XtA
            return gap, dual_norm_XtA
        # This is Ridge regression with primal objective
        #     P(w) = 1/2 ||y - X w||_2^2 + beta/2 ||w||_2^2
        # We use the "alternative" dual formulation with alpha=0, see
        # https://github.com/scikit-learn/scikit-learn/issues/22836
        #     D(v) = -1/2 ||v||_2^2 - v'y - 1/(2 beta) |||X'v||_2^2
        # With v = Xw - y = -R (residual), the dual gap reads
        #     G = P(w) - D(-R)
        #       = ||R||_2^2  + beta/2 ||w||_2^2 - R'y + 1/(2 beta) ||X'R||_2^2
        dual = -0.5 * R_norm2 + Ry - 1 / (2 * beta) * dual_norm_XtA
    else:
        # XtA = X.T @ R - beta * w
        # XtA -= beta * w
        _axpy(n_classes * n_features, -beta, &w[0, 0], 1, &XtA[0, 0], 1)
        dual_norm_XtA = abs_max(n_classes * n_features, &XtA[0, 0])

        if dual_norm_XtA > alpha:
            scale = alpha / dual_norm_XtA
        else:
            scale = 1.0
        dual = -0.5 * (scale**2) * (R_norm2 + beta * w_l2_norm2)
        dual += scale * Ry

    gap = primal - dual
    if gap_smaller_eps and abs(gap) <= 2 * eps * primal:
        gap = 0.0
    return gap, dual_norm_XtA


cdef inline void update_LD_R(
    unsigned int k,
    unsigned int j,
    const floating w_kj,
    const floating[::1, :] X,
    bint X_is_sparse,
    floating[::1] X_data,
    int32_t[::1] X_indices,
    int32_t[::1] X_indptr,
    const floating[::1, :] sw_proba,
    const floating[::1, :] proba,
    const floating[::1, :] H00_pinv_H0,
    bint fit_intercept,
    floating[::1, :] LD_R,
    floating[::1] x_k,
    floating[::1, :] xx,
) noexcept nogil:
    """Update LD_R by X[:, j] * w_kj for class k and feature j.

    Note:
        - LDL = diag(proba) - proba proba', the last term being the outer product.
        - We multiply by sample weights because we want the full hessian,
          sw_proba = sample_weight * proba.
    """
    cdef unsigned int n_samples = LD_R.shape[0]
    cdef unsigned int n_classes = LD_R.shape[1]
    cdef int32_t i, l, i_ind
    cdef int32_t startptr
    cdef int32_t endptr

    if X_is_sparse:
        startptr = X_indptr[j]
        endptr = X_indptr[j + 1]

    # L sqrt(D) R += w * LDL[:, k] * X[:, j]
    if not fit_intercept:
        # numpy:
        #   x_k[:] = w_kj * X[:, j] * sw_proba[:, k]
        #   LD_R[:, k] += x_k
        #   LD_R -= proba * x_k[:, None]
        if not X_is_sparse:
            for i in range(n_samples):
                x_k[i] = w_kj * X[i, j] * sw_proba[i, k]
                LD_R[i, k] += (1 - proba[i, k]) * x_k[i]
        else:
            memset(&x_k[0], 0, n_samples * sizeof(floating))
            for i_ind in range(startptr, endptr):
                i = X_indices[i_ind]
                x_k[i] += X_data[i_ind] * w_kj * sw_proba[i, k]
                LD_R[i, k] += (1 - proba[i, k]) * x_k[i]

        for l in range(n_classes):
            if l != k:
                for i in range(n_samples):
                    LD_R[i, l] -= proba[i, l] * x_k[i]
    else:
        # numpy:
        #   xx[:, :-1] = -H00_pinv_H0[:, jn + k][None, :]
        #   xx[:, -1] = 0
        #   xx[:, k] += X[:, j]
        #   xx *= w_kj * sw_proba
        #   LD_R -= proba * np.sum(xx, axis=1)[:, None]
        #   LD_R += xx

        # x_k[:] = 0
        memset(&x_k[0], 0, n_samples * sizeof(floating))
        # xx[:, :] = 0
        memset(&xx[0, 0], 0, n_samples * n_classes * sizeof(floating))
        # xx[:, k] += X[:, j]
        if not X_is_sparse:
            _copy(n_samples, &X[0, j], 1, &xx[0, k], 1)
        else:
            for i_ind in range(startptr, endptr):
                i = X_indices[i_ind]
                xx[i] = X_data[i_ind]
        for l in range(n_classes - 1):
            for i in range(n_samples):
                xx[i, l] -= H00_pinv_H0[l, k + n_classes * j]
                xx[i, l] *= w_kj * sw_proba[i, l]
                LD_R[i, l] += xx[i, l]
                x_k[i] += xx[i, l]  # x_k = np.sum(xx, axis=1)
        if k == n_classes - 1:  # H00_pinv_H0[k, :] == 0
            l = n_classes - 1
            for i in range(n_samples):
                xx[i, l] *= w_kj * sw_proba[i, l]
                LD_R[i, l] += xx[i, l]
                x_k[i] += xx[i, l]  # x_k = np.sum(xx, axis=1)
        for l in range(n_classes):
            for i in range(n_samples):
                LD_R[i, l] -= proba[i, l] * x_k[i]


def enet_coordinate_descent_multinomial(
    W,
    float64_t alpha,
    float64_t beta,
    X,
    sample_weight,
    raw_prediction,
    grad_pointwise,
    proba,
    bint fit_intercept,
    unsigned int max_iter,
    float64_t tol,
    bint do_screening=True,
    bint early_stopping=True,
    int verbose=False,
):
    """Cython coordinate descent algorithm for Elastic-Net multinomial regression.

    See function enet_coordinate_descent.
    We minimize the primal

        P(w) = 1/2 ||y - X w||_2^2 + alpha ||w||_1 + beta/2 ||w||_2^2

    But we replace the variables y and X such that P(w) corresponds to the 2nd order
    approximation of the multinomial loss

        1/2 w' H w + (G' - coef' H) w + alpha ||w||_1 + beta/2 ||w||_2^2

    - w = W.ravel(order="F")
    - G = gradient = X.T @ (proba - Y)  with Y_k = (y==k)
    - H = X' * (diag(p) - pp') * X, the full multinomial hessian
      (* is kind of a Kronecker multiplication)

    We use the analytical LDL decomposition of diag(p) - pp' = L D L' of
    Tanabe & Sagae (1992) to replace

        X -> A = sqrt(D) L' * X
        y -> b = (L sqrt(D))^-1 (LDL' * X coef - g)
               = A coef - (L sqrt(D))^-1 g

    Which gives (up to a constant)

        P(w) = 1/2 ||b - A w||_2^2 + alpha ||w||_1 + beta/2 ||w||_2^2

    As the matrix A has n_samples * n_features * n_classes * n_classes entries, we
    avoid to explicitly build it.

    Note:
        - H = A'A
        - A'b = (coef' H - G') w

    Intercepts
    ----------
    If intercepts w0 of shape (n_classes,) are included, we have

        P(w) = 1/2 ||b - A w - A0 w0||_2^2 + alpha ||w||_1 + beta/2 ||w||_2^2

    with A0 = sqrt(D) L' 1 and 1 = 1_{n_samples, n_classes} for the n_classes intercept
    columns.

    Minimizing for w0 gives

        w0 = (A0' A0)^(-1) (A0' b - A0' A w)

        P(w) = 1/2 ||(I - A0 (A0' A0)^(-1) A0') (b - A w)||_2^2 + ...
             = 1/2 ||(b - A0 H00^(-1) q0) - (A - A0 H00^(-1) H0) w||_2^2 + ...

    We therefore replace

        A -> A - A0 H00^(-1) H0 = sqrt(D) L' (X - 1 H00^(-1) H0)

        b -> b - A0 H00^(-1) q0 = sqrt(D) L' (X - 1 H00^(-1) H0) coef
                                - (L sqrt(D))^-1 (I - LDL' 1 H00^(-1) 1') g

    Note:
        - A0 = sqrt(D) L' 1_{n_samples, n_classes}
        - H00 = A0 A0', i.e. the part of the hessian for the the intercept alone
        - H0 = A0' A, i.e. the part of the hessian mixing intercepts with the features:
          H[-n_classes, :-n_classes]
        - q0 = A0' b = 1' (LDL' X coef - g)

    Tracking the Residual
    ---------------------
    Usually, CD tracks the residual R = y - X w -> b - Aw. In the multinomial case, it
    is advantageous to track the rotated residual instead:

        LD_R = L sqrt(D) R

    This makes the coordinate update simple and fast involving just

        X[:, j] @ LD_R[:, k]

    The update of the rotated residual involves basically multiplying by the full
    LDL = diag(p) - pp', i.e.

        L sqrt(D) R -= (w_new_kj - w_old_kj) * LDL[:, k] * X[:, j]

    with LDL[:, k] begin the k-th column of the LDL matrix, having shape
    (n_sampels, n_classes).

    Returns
    -------
    W : ndarray of shape (n_classes, n_features + fit_intercept)
        ElasticNet coefficients.
    gap : float
        Achieved dual gap.
    tol : float
        Equals input `tol` times `np.dot(y, y)`. The tolerance used for the dual gap.
    n_iter : int
        Number of coordinate descent iterations.
    """
    dtype = X.dtype
    # cdef time_t time_total = time(NULL)
    # cdef time_t time_pre = 0
    # cdef time_t time_gap = 0
    # cdef time_t time_screen = 0
    # cdef time_t time_w = 0
    # cdef time_t time_r = 0
    # cdef time_t tic = time(NULL)

    # get the data information into easy vars
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]
    cdef unsigned int n_classes = W.shape[0]
    cdef bint X_is_sparse = sparse.issparse(X)

    if not W.flags["F_CONTIGUOUS"]:
        raise ValueError("Coefficient W must be F-contiguous.")

    if sample_weight is None:
        sw_sum = n_samples
        sw = np.full(fill_value=1 / sw_sum, shape=n_samples, dtype=dtype)
    else:
        sw_sum = np.sum(sample_weight)
        sw = sample_weight / sw_sum
    sqrt_sw = np.sqrt(sw)
    norm2_cols_X = np.empty((n_classes, n_features), dtype=dtype, order="C")
    # Note: Full hessian H of LinearModelLoss.gradient_hessian(coef=W, X=X, y=y)
    # divides by sw_sum. So each LDL.XXX_matmul computation must be multiplied by
    # sqrt(samples_weight / sw_sum) = sqrt_sw, i.e. as if sqrt(D) would have
    # been multiplied by sqrt_sw.
    LDL = Multinomial_LDL_Decomposition(proba=proba)
    XtA = np.zeros((n_classes, n_features), dtype=dtype, order="F")  # TODO: best order?
    R = np.empty((n_samples, n_classes), dtype=dtype, order="F")
    LD_R = np.empty_like(R, order="F")  # memory buffer for LDL.L_sqrt_D_matmul(R)
    sw_proba = sw[:, None] * proba  # precompute this one
    # memory buffers for temporaries
    x_k = np.empty(n_samples, dtype=dtype)
    xx = np.empty((n_samples, n_classes), dtype=dtype, order="F")
    t_k = np.empty(n_samples, dtype=dtype)
    if fit_intercept:
        tt = np.empty((n_samples, n_classes), dtype=dtype, order="F")

    cdef float64_t d_j
    cdef float64_t Xj_theta
    cdef float64_t tmp
    cdef float64_t w_kj
    cdef float64_t d_w_max
    cdef float64_t w_max
    cdef float64_t d_w_j
    cdef float64_t gap = tol + 1.0
    cdef float64_t d_w_tol = tol
    cdef float64_t dual_norm_XtA
    cdef uint32_t[::1] n_active = np.full(shape=n_classes, fill_value=n_features, dtype=np.uint32)
    cdef uint32_t[:, ::1] active_set
    # TODO: use binset instead of array of bools
    cdef uint8_t[:, ::1] excluded_set
    cdef unsigned int i, j, k, jn, l
    cdef int32_t endptr, startptr, i_ind  # indexing sparse X
    cdef unsigned int n_iter = 0
    cdef bint at_least_one_feature_updated

    if alpha <= 0:
        do_screening = False

    if do_screening:
        # active_set maps [k, :n_active[k]] -> j
        active_set = np.empty((n_classes, n_features), dtype=np.uint32)
        excluded_set = np.empty((n_classes, n_features), dtype=np.uint8)

    W_original = W
    H00_pinv_H0 = None
    if fit_intercept:
        # See LinearModelLoss.gradient_hessian and NewtonCDGramSolver.inner_solve.
        # We set the intercept of the last class to zero, loops and shapes often
        # have n_classes - 1 instead of n_classes. The missing entries are all zeros,
        # but we often do not add those zeros explicitly.
        W0 = W[:, -1]
        # W0[-1] = 0  # intercept of last class
        W0[n_classes - 1] = 0  # Cython does not like negative indices
        W = W[:, :-1]
        # H00 = H[-n_classes:-1, -n_classes:-1] = Hessian part of intercepts
        H00 = np.zeros((n_classes - 1, n_classes - 1), dtype=dtype)
        # H0 = H[-n_classes:-1, :-n_classes] = Hessian part mixing intercepts with
        # features
        H0 = np.zeros((n_classes - 1, n_classes * n_features), dtype=dtype)
        H0_coef = np.zeros((n_classes - 1,), dtype=dtype)  # 1' LDL raw_prediction
        h = x_k
        for k in range(n_classes - 1):
            # Diagonal terms h_kk.
            h[:] = proba[:, k] * (1 - proba[:, k]) * sw
            H00[k, k] = h.sum()
            H0[k, k::n_classes] = X.T @ h
            H0_coef[k] += (h * raw_prediction[:, k]).sum()
            # Off diagonal terms (in classes) hess_kl.
            for l in range(k + 1, n_classes - 1):
                # Upper triangle (in classes).
                h[:] = -proba[:, k] * proba[:, l] * sw
                H00[k, l] = h.sum()
                H00[l, k] = H00[k, l]
                H0[k, l::n_classes] = X.T @ h
                H0[l, k::n_classes] = H0[k, l::n_classes]
                H0_coef[k] += (h * raw_prediction[:, l]).sum()
                H0_coef[l] += (h * raw_prediction[:, k]).sum()
            # So far omitted term for last class where we still need it:
            l = n_classes - 1
            h[:] = -proba[:, k] * proba[:, l] * sw
            H0[k, l::n_classes] = X.T @ h
            H0_coef[k] += (h * raw_prediction[:, l]).sum()

        H00_pinv = np.linalg.pinv(H00, hermitian=True)
        # H0_coef is the same as:
        # H0_coef = H0 @ W.ravel(order="F") + H00 @ W0
        grad_p_sum = grad_pointwise[:, :-1].sum(axis=0)
        q0 = H0_coef - grad_p_sum
        # H00_pinv_H0 = H00_pinv @ H0  # shape (n_classes - 1, n_classes * n_features)
        H00_pinv_H0 = np.zeros(
            (n_classes - 1, n_classes * n_features), dtype=dtype, order="F"
        )
        # Double transpose to make the result F-contiguous.
        np.dot(H0.T, H00_pinv.T, out=H00_pinv_H0.T)

        # Centering X as
        #   X -= (H00_pinv_H0)[None, :]
        # is not an option because it would mean n_classes^2 copies of X.
        # Center raw_prediction: += -intercepts - 1 H00^(-1) H0) coef
        raw_prediction = raw_prediction - W0  # creates a copy
        raw_prediction[:, :-1] -= (H00_pinv_H0 @ W.ravel(order="F"))[None, :]
        # Center grad_pointwise: g -= LDL' 1 H00^(-1) 1' g
        t = H00_pinv @ grad_p_sum  # H00^(-1) 1' g
        for k in range(n_classes - 1):
            h = proba[:, k] * (1 - proba[:, k]) * sw
            grad_pointwise[:, k] -= h * t[k]
            for l in range(k + 1, n_classes - 1):
                h = -proba[:, k] * proba[:, l] * sw
                grad_pointwise[:, k] -= h * t[l]
                grad_pointwise[:, l] -= h * t[k]

    # A w = sqrt(D) L' X w = sqrt(D) L' raw = sqrt_D_Lt_raw
    sqrt_D_Lt_raw = LDL.sqrt_D_Lt_matmul(raw_prediction.copy(order="F")) * sqrt_sw[:, None]

    # residual R = b - A w = -(L sqrt(D))^-1 g
    # R[:, :] = LDL.inverse_L_sqrt_D_matmul(-grad_pointwise / sqrt_sw[:, None])
    np.divide(-grad_pointwise, sqrt_sw[:, None], out=R)
    LDL.inverse_L_sqrt_D_matmul(R)

    # b = A w - (L sqrt(D))^-1 g
    b = sqrt_D_Lt_raw + R  # shape (n_samples, n_classes)

    # tol = tol * linalg.norm(Y, ord='fro') ** 2
    tol = tol * np.linalg.norm(b, ord="fro") ** 2

    # Rotated residual: L sqrt(D) R sqrt(sw)
    # LD_R[:, :] = LDL.L_sqrt_D_matmul(R * sqrt_sw[:, None])
    np.multiply(R, sqrt_sw[:, None], out=LD_R)
    LDL.L_sqrt_D_matmul(LD_R)
    # IMPORTANT NOTE: With fit_intercept=True, np.sum(LD_R, axis=0) == 0.

    # Convention: memoryviews have a trailing underscore.
    cdef float64_t[::1, :] W_ = W
    cdef const float64_t[::1, :] X_
    cdef float64_t[::1] X_data
    cdef int32_t[::1] X_indices
    cdef int32_t[::1] X_indptr
    cdef float64_t[::1, :] b_ = b
    cdef float64_t[::1] sqrt_sw_ = sqrt_sw
    cdef float64_t[:, ::1] norm2_cols_X_ = norm2_cols_X
    cdef float64_t[::1, :] R_ = R
    cdef float64_t[::1, :] LD_R_ = LD_R
    cdef float64_t[::1, :] proba_ = proba
    cdef float64_t[::1, :] XtA_ = XtA
    cdef float64_t[::1, :] H00_pinv_H0_ = H00_pinv_H0
    cdef float64_t[::1, :] sw_proba_ = sw_proba
    cdef float64_t[::1, :] q_inv_ = LDL.q_inv
    cdef float64_t[::1, :] sqrt_d_ = LDL.sqrt_d
    cdef float64_t[::1] t_k_ = t_k
    cdef float64_t[::1] x_k_ = x_k
    cdef float64_t[::1, :] tt_
    cdef float64_t[::1, :] xx_ = xx
    if fit_intercept:
        tt_ = tt

    if X_is_sparse:
        if not (
            X.data.flags.c_contiguous
            and X.indices.flags.c_contiguous
            and X.indptr.flags.c_contiguous
        ):
            raise ValueError("Sparse X must have contiguous underlying arrays.")
        X_data = X.data
        X_indices = X.indices
        X_indptr = X.indptr
    else:
        X_ = X

    # time_pre += time(NULL) - tic
    # tic = time(NULL)

    # Check convergence before entering the main loop.
    gap, dual_norm_XtA = gap_enet_multinomial(
        w=W_,
        alpha=alpha,
        beta=beta,
        X=X_,
        X_is_sparse=X_is_sparse,
        X_data=X_data,
        X_indices=X_indices,
        X_indptr=X_indptr,
        y=b_,
        R=R_,
        LD_R=LD_R_,
        H00_pinv_H0=H00_pinv_H0_,
        fit_intercept=fit_intercept,
        XtA=XtA_,
        gap_smaller_eps=False,
    )
    if early_stopping and gap <= tol:
        if verbose:
            print(
                f"  inner coordinate descent solver stops early with gap={float(gap)}"
                f"  <= tol={float(tol)}"
            )
        return np.asarray(W), gap, tol, 0

    # time_gap += time(NULL) - tic
    # tic = time(NULL)

    # Compute squared norms of the columns of X.
    # Same as norm2_cols_X = np.square(X).sum(axis=0)
    # norm2_cols_X = np.einsum(
    #     "ij,ij->j", X, X, dtype=dtype, order="C"
    # )
    # X -> sqrt(D) L' X = A
    # sum_{i} X_{ij} X_{ij} -> sum_i X_{ij} LDL_i X_{ij}
    # These are just the diagonal elements of the full hessian H.
    # for k in range(n_classes):
    #     h = proba[:, k] * (1 - proba[:, k]) * sw
    #     norm2_cols_X[k, :] = np.einsum("ij, i, ij -> j", X, h, X)
    np.subtract(1, proba, out=xx)  # xx = 1 - proba
    xx *= sw_proba
    if not X_is_sparse:
        np.einsum("ij, ik, ij -> kj", X, xx, X, order="C", out=norm2_cols_X)
    else:
        memset(&norm2_cols_X_[0, 0], 0, n_classes * n_features * sizeof(float64_t))
        for j in range(n_features):
            startptr = X_indptr[j]
            endptr = X_indptr[j + 1]
            for i_ind in range(startptr, endptr):
                i = X_indices[i_ind]
                tmp = X_data[i_ind] ** 2
                for k in range(n_classes):
                    norm2_cols_X_[k, j] += tmp * xx_[i, k]
    if fit_intercept:
        # See Q_centered = Q - Q0.T @ Q00_inv @ Q0 in NewtonCDGramSolve.
        norm2_cols_X -= np.einsum("ji, ji -> i", H0, H00_pinv_H0).reshape(
            n_classes, -1, order="F"
        )
        # Guarantee norm2_cols_X >= 0:
        norm2_cols_X[norm2_cols_X < 0] = 0

    # time_pre += time(NULL) - tic
    # tic = time(NULL)

    # Gap Safe Screening Rules, see https://arxiv.org/abs/1802.07481, Eq. 11
    if do_screening:
        for k in range(n_classes):
            n_active[k] = 0
            for j in range(n_features):
                if norm2_cols_X_[k, j] == 0:
                    W_[k, j] = 0
                    excluded_set[k, j] = 1
                    continue
                Xj_theta = XtA_[k, j] / fmax(alpha, dual_norm_XtA)  # X[:,j] @ dual_theta
                d_j = (1 - fabs(Xj_theta)) / sqrt(norm2_cols_X_[k, j] + beta)
                if d_j <= sqrt(2 * gap) / alpha:
                    # include feature j of class k
                    active_set[k, n_active[k]] = j
                    excluded_set[k, j] = 0
                    n_active[k] += 1
                else:
                    if W_[k, j] != 0:
                        # R += w[j] * X[:,j] -> w[k, j] A[:, kj]
                        # LD_R += w[k, j] * LDL[:, k] * X[:, j]
                        update_LD_R(
                            k=k, j=j, w_kj=W_[k, j], X=X_,
                            X_is_sparse=X_is_sparse, X_data=X_data,
                            X_indices=X_indices, X_indptr=X_indptr,
                            sw_proba=sw_proba_,
                            proba=proba_, H00_pinv_H0=H00_pinv_H0_,
                            fit_intercept=fit_intercept, LD_R=LD_R_, x_k=x_k_, xx=xx_,
                        )
                        W_[k, j] = 0
                    excluded_set[k, j] = 1
        if np.sum(n_active) == 0:
            # We want to guarantee at least one CD step.
            # n_active[:] = n_features
            # active_set[:, :] = np.arange(n_features, dtype=np.uint32)[None, :]
            # excluded_set[:, :] = 0
            for k in range(n_classes):
                n_active[k] = n_features
                for j in range(n_features):
                    active_set[k, j] = j
                    excluded_set[k, j] = 0

    # time_screen += time(NULL) - tic

    with nogil:
        for n_iter in range(max_iter):
            w_max = 0.0
            d_w_max = 0.0
            for k in range(n_classes):  # Loop over coordinates
                # tic = time(NULL)
                at_least_one_feature_updated = False
                # t_k[:] = 0
                memset(&t_k_[0], 0, n_samples * sizeof(float64_t))
                if fit_intercept:
                    # tt[:, :] = 0
                    memset(&tt_[0, 0], 0, n_samples * n_classes * sizeof(float64_t))
                # time_r += time(NULL) - tic

                for j in range(n_active[k]):  # Loop over coordinates
                    # tic = time(NULL)
                    if do_screening:
                        j = active_set[k, j]

                    if norm2_cols_X_[k, j] == 0.0:
                        continue
                    w_kj = W_[k, j]  # Store previous value

                    if X_is_sparse:
                        startptr = X_indptr[j]
                        endptr = X_indptr[j + 1]

                    # tmp = X[:,j].T @ (R + w_j * X[:,j])
                    #    -> A[:,jk].T @ (R + w_kj * A[:,kj])
                    # With precomputed LD_R: A'R = X' L sqrt(D) R
                    # tmp = X[:, j] @ LD_R[:, k] + w_kj * norm2_cols_X[k, j]
                    if not X_is_sparse:
                        tmp = _dot(n_samples, &X_[0, j], 1, &LD_R_[0, k], 1)
                    else:
                        tmp = 0.0
                        for i_ind in range(startptr, endptr):
                            i = X_indices[i_ind]
                            tmp += LD_R_[i, k] * X_data[i_ind]
                    tmp += w_kj * norm2_cols_X_[k, j]
                    # Note: With fit_intercept=True, np.sum(LD_R, axis=0) == 0. So the
                    # additional term
                    #   tmp -= H00_pinv_H0[:, n_classes * j + k].T @ np.sum(LD_R, axis=0)
                    # is zero.

                    W_[k, j] = (
                        fsign(tmp) * fmax(fabs(tmp) - alpha, 0)
                        / (norm2_cols_X_[k, j] + beta)
                    )
                    # time_w += time(NULL) - tic
                    # tic = time(NULL)

                    if W_[k, j] != w_kj:
                        at_least_one_feature_updated = True
                        # Update residual
                        # R -= (w[j] - w_j) * X[:,j] -> (w[kj] - w_kj) * A[:,kj]
                        if not fit_intercept:
                            # Update rotated residual LD_R instead of R
                            # L sqrt(D) R -= (w[kj] - w_kj) * LDL[:, k] * X[:, j]
                            # numpy:
                            #   x_k[:] = (w_kj - W[k, j]) * X[:, j] * sw_proba[:, k]
                            #   LD_R[:, k] += (1 - proba[:, k]) * x_k  # diagonal LDL
                            #   for l in range(n_classes):
                            #       if l != k:
                            #           LD_R[:, l] -= proba[:, l] * x_k
                            # faster version:
                            #   LD_R[:, k] += x_k
                            #   LD_R -= proba * x_k[:, None]
                            if not X_is_sparse:
                                for i in range(n_samples):
                                    x_k_[i] = (w_kj - W_[k, j]) * X_[i, j] * sw_proba_[i, k]
                                    t_k_[i] += x_k_[i]  # accumulation for postponed update
                                    LD_R_[i, k] += (1 - proba_[i, k]) * x_k_[i]
                            else:
                                memset(&x_k_[0], 0, n_samples * sizeof(float64_t))
                                for i_ind in range(startptr, endptr):
                                    i = X_indices[i_ind]
                                    x_k_[i] = (w_kj - W_[k, j]) * X_data[i_ind] * sw_proba_[i, k]
                                    t_k_[i] += x_k_[i]
                                    LD_R_[i, k] += (1 - proba_[i, k]) * x_k_[i]
                            # Postpone updates of LD_R for classes != k.
                            # for l in range(n_classes):
                            #     if l != k:
                            #         for i in range(n_samples):
                            #             LD_R_[i, l] -= proba_[i, l] * x_k_[i]
                        else:
                            # Update rotated residual LD_R instead of R
                            jn = n_classes * j
                            # numpy:
                            #   xx[:, :-1] = -H00_pinv_H0[:, jn + k][None, :]
                            #   xx[:, -1] = 0
                            #   xx[:, k] += X[:, j]
                            #   xx *= (w_kj - W[k, j]) * sw_proba
                            # L sqrt(D) R += LDL xx
                            # numpy:
                            #   for l in range(n_classes):
                            #       LD_R[:, l] += (1 - proba[:, l]) * xx[:, l]
                            #       for m in range(l + 1, n_classes):
                            #           LD_R[:, l] -= proba[:, l] * xx[:, m]
                            #           LD_R[:, m] -= proba[:, m] * xx[:, l]
                            # faster version:
                            #   LD_R -= proba * np.sum(xx, axis=1)[:, None]
                            #   LD_R += xx

                            # x_k[:] = 0
                            memset(&x_k_[0], 0, n_samples * sizeof(float64_t))
                            # xx[:, :] = 0
                            memset(&xx_[0, 0], 0, n_samples * n_classes * sizeof(float64_t))
                            # xx[:, k] += X[:, j]
                            if not X_is_sparse:
                                _copy(n_samples, &X_[0, j], 1, &xx_[0, k], 1)
                            else:
                                for i_ind in range(startptr, endptr):
                                    i = X_indices[i_ind]
                                    xx_[i, k] = X_data[i_ind]
                            for l in range(n_classes - 1):
                                for i in range(n_samples):
                                    xx_[i, l] -= H00_pinv_H0_[l, jn + k]
                                    xx_[i, l] *= (w_kj - W_[k, j]) * sw_proba_[i, l]
                                    tt_[i, l] += xx_[i, l]
                                    # LD_R_[i, l] += xx_[i, l]  # postponed
                                    x_k_[i] += xx_[i, l]  # x_k = np.sum(xx, axis=1)
                            if k == n_classes - 1:  # H00_pinv_H0[k, :] == 0
                                l = n_classes - 1
                                for i in range(n_samples):
                                    xx_[i, l] *= (w_kj - W_[k, j]) * sw_proba_[i, l]
                                    tt_[i, l] += xx_[i, l]
                                    # LD_R_[i, l] += xx_[i, l]  # postponed
                                    x_k_[i] += xx_[i, l]  # x_k = np.sum(xx, axis=1)
                            # Postpone updates of LD_R for classes != k.
                            for i in range(n_samples):
                                LD_R_[i, k] += xx_[i, k] - proba_[i, k] * x_k_[i]
                                t_k_[i] += x_k_[i]
                            # The following is postponed.
                            # for l in range(n_classes):
                            #     for i in range(n_samples):
                            #         LD_R_[i, l] -= proba_[i, l] * x_k_[i]

                    # time_r += time(NULL) - tic

                    # update the maximum absolute coefficient update
                    d_w_j = fabs(W_[k, j] - w_kj)
                    d_w_max = fmax(d_w_max, d_w_j)
                    w_max = fmax(w_max, fabs(W_[k, j]))

                # tic = time(NULL)

                # Postponed updates of LD_R for classes != k .
                if at_least_one_feature_updated:
                    if not fit_intercept:
                        # We accumulated t_k_ = (X @ (w_k - W[k, :])) * sw_proba[:, k]
                        for l in range(n_classes):
                            if l != k:
                                for i in range(n_samples):
                                    LD_R_[i, l] -= proba_[i, l] * t_k_[i]
                    else:
                        for l in range(n_classes):
                            if l != k:
                                for i in range(n_samples):
                                    LD_R_[i, l] += tt_[i, l] - proba_[i, l] * t_k_[i]

                # time_r += time(NULL) - tic

            if w_max == 0.0 or d_w_max / w_max <= d_w_tol or n_iter == max_iter - 1:
                # tic = time(NULL)

                # The biggest coordinate update of this iteration was smaller than the
                # tolerance: check the duality gap as ultimate stopping criterion.
                for k in range(n_classes):
                    for i in range(n_samples):
                        R_[i, k] = LD_R_[i, k] / sqrt_sw_[i]
                # LDL.inverse_L_sqrt_D_matmul(R)
                inverse_L_sqrt_D_matmul(p=proba_, q_inv=q_inv_, sqrt_d=sqrt_d_, x=R_, buffer=x_k_)
                gap, dual_norm_XtA = gap_enet_multinomial(
                    w=W_,
                    alpha=alpha,
                    beta=beta,
                    X=X_,
                    X_is_sparse=X_is_sparse,
                    X_data=X_data,
                    X_indices=X_indices,
                    X_indptr=X_indptr,
                    y=b_,
                    R=R_,
                    LD_R=LD_R_,
                    H00_pinv_H0=H00_pinv_H0_,
                    fit_intercept=fit_intercept,
                    XtA=XtA_,
                    gap_smaller_eps=True,
                )
                if verbose:
                    with gil:
                        print(
                            f"  inner coordinate descent solver at iter {n_iter} "
                            f"checks for gap={float(gap)} <= tol={float(tol)} ? {gap <= tol}"
                        )
                # time_gap += time(NULL) - tic
                # tic = time(NULL)
                if gap <= tol:
                    # return if we reached desired tolerance
                    break

                # Gap Safe Screening Rules, see https://arxiv.org/abs/1802.07481, Eq. 11
                if do_screening:
                    for k in range(n_classes):
                        n_active[k] = 0
                        for j in range(n_features):
                            if excluded_set[k, j]:
                                continue
                            Xj_theta = XtA_[k, j] / fmax(alpha, dual_norm_XtA)  # X[:,j] @ dual_theta
                            d_j = (1 - fabs(Xj_theta)) / sqrt(norm2_cols_X_[k, j] + beta)
                            if d_j <= sqrt(2 * gap) / alpha:
                                # include feature j of class k
                                active_set[k, n_active[k]] = j
                                excluded_set[k, j] = 0
                                n_active[k] += 1
                            else:
                                if W_[k, j] != 0:
                                    # R += w[j] * X[:,j] -> w[k, j] A[:, kj]
                                    # LD_R += w[k, j] * LDL[:, k] * X[:, j]
                                    update_LD_R(
                                        k=k, j=j, w_kj=W_[k, j], X=X_,
                                        X_is_sparse=X_is_sparse, X_data=X_data,
                                        X_indices=X_indices, X_indptr=X_indptr,
                                        sw_proba=sw_proba_,
                                        proba=proba_, H00_pinv_H0=H00_pinv_H0_,
                                        fit_intercept=fit_intercept, LD_R=LD_R_, x_k=x_k_, xx=xx_,
                                    )
                                W_[k, j] = 0
                                excluded_set[k, j] = 1

                # time_screen += time(NULL) - tic
        else:
            with gil:
                # for/else, runs if for doesn't end with a `break`
                message = (
                    "Objective did not converge. You might want to increase "
                    "the number of iterations, check the scale of the "
                    "features or consider increasing regularisation."
                    f" Duality gap: {gap:.6e}, tolerance: {tol:.3e}"
                )
                warnings.warn(message, ConvergenceWarning)

    if fit_intercept:
        # W0[:-1] = H00_pinv @ (q0 - H0 @ W.ravel(order="F"))
        W0[:n_classes - 1] = H00_pinv @ (q0 - H0 @ W.ravel(order="F"))
        # W0[-1] = 0  # was already set at the beginning.

    # time_total = time(NULL) - time_total
    # print(
    #     f"{time_pre=:6.4f} {time_gap=:6.4f} {time_screen=:6.4f}\n"
    #     f"{time_w=:6.4f} {time_r=:6.4f} {time_total=:6.4f}"
    # )
    return np.asarray(W_original), gap, tol, n_iter + 1
