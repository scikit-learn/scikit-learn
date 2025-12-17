# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from libc.math cimport fabs, sqrt
import numpy as np

from cython cimport floating
import warnings
from sklearn.exceptions import ConvergenceWarning

from sklearn.utils._cython_blas cimport (
    _axpy, _dot, _asum, _gemv, _nrm2, _copy, _scal
)
from sklearn.utils._cython_blas cimport ColMajor, Trans, NoTrans
from sklearn.utils._typedefs cimport uint8_t, uint32_t
from sklearn.utils._random cimport our_rand_r


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
) noexcept nogil:
    """Compute dual gap for use in enet_coordinate_descent."""
    cdef floating gap = 0.0
    cdef floating dual_norm_XtA
    cdef floating R_norm2
    cdef floating w_norm2 = 0.0
    cdef floating l1_norm
    cdef floating A_norm2
    cdef floating const_

    # XtA = X.T @ R - beta * w
    _copy(n_features, &w[0], 1, &XtA[0], 1)
    _gemv(ColMajor, Trans, n_samples, n_features, 1.0, &X[0, 0],
          n_samples, &R[0], 1,
          -beta, &XtA[0], 1)

    if positive:
        dual_norm_XtA = max(n_features, &XtA[0])
    else:
        dual_norm_XtA = abs_max(n_features, &XtA[0])

    # R_norm2 = R @ R
    R_norm2 = _dot(n_samples, &R[0], 1, &R[0], 1)

    # w_norm2 = w @ w
    if beta > 0:
        w_norm2 = _dot(n_features, &w[0], 1, &w[0], 1)

    if (dual_norm_XtA > alpha):
        const_ = alpha / dual_norm_XtA
        A_norm2 = R_norm2 * (const_ ** 2)
        gap = 0.5 * (R_norm2 + A_norm2)
    else:
        const_ = 1.0
        gap = R_norm2

    l1_norm = _asum(n_features, &w[0], 1)

    gap += (
        alpha * l1_norm
        - const_ * _dot(n_samples, &R[0], 1, &y[0], 1)  # R @ y
        + 0.5 * beta * (1 + const_ ** 2) * w_norm2
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
):
    """
    Cython version of the coordinate descent algorithm for Elastic-Net regression.

    The algorithm mostly follows [Friedman 2010].
    We minimize the primal

        P(w) = 1/2 ||y - X w||_2^2 + alpha ||w||_1 + beta/2 ||w||_2^2

    The dual for beta = 0, see e.g. [Fercoq 2015] with v = alpha * theta, is

        D(v) = -1/2 ||v||_2^2 + y v

    with dual feasible condition ||X^T v||_inf <= alpha.
    For beta > 0, one uses extended versions of X and y by adding n_features rows

        X -> (           X)    y -> (y)
             (sqrt(beta) I)         (0)

    Note that the residual y - X w is an important ingredient for the estimation of a
    dual feasible point v.
    At optimum of primal w* and dual v*, one has

        v = y* - X w*

    The duality gap is

        G(w, v) = P(w) - D(v) <= P(w) - P(w*)

    The final stopping criterion is based on the duality gap

        tol ||y||_2^2 <= G(w, v)

    The tolerance here is multiplied by ||y||_2^2 to have an inequality that scales the
    same on both sides and because one has G(0, 0) = 1/2 ||y||_2^2.

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

    if alpha == 0 and beta == 0:
        warnings.warn("Coordinate descent with no regularization may lead to "
                      "unexpected results and is discouraged.")

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
            n_samples, n_features, w, alpha, beta, X, y, R, XtA, positive
        )
        if gap <= tol:
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
            ):
                # the biggest coordinate update of this iteration was smaller
                # than the tolerance: check the duality gap as ultimate
                # stopping criterion
                gap, dual_norm_XtA = gap_enet(
                    n_samples, n_features, w, alpha, beta, X, y, R, XtA, positive
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
                            # R += w[j] * X[:,j]
                            _axpy(n_samples, w[j], &X[0, j], 1, &R[0], 1)
                            w[j] = 0
                            excluded_set[j] = 1

        else:
            # for/else, runs if for doesn't end with a `break`
            with gil:
                message = (
                    message_conv +
                    f" Duality gap: {gap:.3e}, tolerance: {tol:.3e}"
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
) noexcept nogil:
    """Compute dual gap for use in sparse_enet_coordinate_descent."""
    cdef floating gap = 0.0
    cdef floating dual_norm_XtA
    cdef floating R_norm2
    cdef floating w_norm2 = 0.0
    cdef floating l1_norm
    cdef floating A_norm2
    cdef floating const_
    cdef unsigned int i, j

    # XtA = X.T @ R - beta * w
    # sparse X.T @ dense R
    for j in range(n_features):
        XtA[j] = 0.0
        for i in range(X_indptr[j], X_indptr[j + 1]):
            XtA[j] += X_data[i] * R[X_indices[i]]

        if center:
            XtA[j] -= X_mean[j] * R_sum
        XtA[j] -= beta * w[j]

    if positive:
        dual_norm_XtA = max(n_features, &XtA[0])
    else:
        dual_norm_XtA = abs_max(n_features, &XtA[0])

    # R_norm2 = R @ R
    if no_sample_weights:
        R_norm2 = _dot(n_samples, &R[0], 1, &R[0], 1)
    else:
        R_norm2 = 0.0
        for i in range(n_samples):
            # R is already multiplied by sample_weight
            if sample_weight[i] != 0:
                R_norm2 += (R[i] ** 2) / sample_weight[i]

    # w_norm2 = w @ w
    if beta > 0:
        w_norm2 = _dot(n_features, &w[0], 1, &w[0], 1)

    if (dual_norm_XtA > alpha):
        const_ = alpha / dual_norm_XtA
        A_norm2 = R_norm2 * const_**2
        gap = 0.5 * (R_norm2 + A_norm2)
    else:
        const_ = 1.0
        gap = R_norm2

    l1_norm = _asum(n_features, &w[0], 1)

    gap += (
        alpha * l1_norm
        - const_ * _dot(n_samples, &R[0], 1, &y[0], 1)  # R @ y
        + 0.5 * beta * (1 + const_ ** 2) * w_norm2
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
        )
        if gap <= tol:
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

            if w_max == 0.0 or d_w_max / w_max <= d_w_tol or n_iter == max_iter - 1:
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
                    f" Duality gap: {gap:.3e}, tolerance: {tol:.3e}"
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
) noexcept nogil:
    """Compute dual gap for use in enet_coordinate_descent."""
    cdef floating gap = 0.0
    cdef floating dual_norm_XtA
    cdef floating R_norm2
    cdef floating w_norm2 = 0.0
    cdef floating l1_norm
    cdef floating A_norm2
    cdef floating const_
    cdef floating q_dot_w
    cdef floating wQw
    cdef unsigned int j

    # q_dot_w = w @ q
    q_dot_w = _dot(n_features, &w[0], 1, &q[0], 1)

    # XtA = X.T @ R - beta * w = X.T @ y - X.T @ X @ w - beta * w
    for j in range(n_features):
        XtA[j] = q[j] - Qw[j] - beta * w[j]

    if positive:
        dual_norm_XtA = max(n_features, &XtA[0])
    else:
        dual_norm_XtA = abs_max(n_features, &XtA[0])

    # wQw = w @ Q @ w
    wQw = _dot(n_features, &w[0], 1, &Qw[0], 1)
    # R_norm2 = R @ R
    R_norm2 = y_norm2 + wQw - 2.0 * q_dot_w

    # w_norm2 = w @ w
    if beta > 0:
        w_norm2 = _dot(n_features, &w[0], 1, &w[0], 1)

    if (dual_norm_XtA > alpha):
        const_ = alpha / dual_norm_XtA
        A_norm2 = R_norm2 * (const_ ** 2)
        gap = 0.5 * (R_norm2 + A_norm2)
    else:
        const_ = 1.0
        gap = R_norm2

    l1_norm = _asum(n_features, &w[0], 1)

    gap += (
        alpha * l1_norm
        - const_ * (y_norm2 - q_dot_w)  # -const_ * R @ y
        + 0.5 * beta * (1 + const_ ** 2) * w_norm2
    )
    return gap, dual_norm_XtA


cdef inline uint32_t screen_features_enet_gram(
    const floating[:, ::1] Q,
    const floating[::1] XtA,
    floating[::1] w,
    floating[::1] Qw,
    uint32_t[::1] active_set,
    uint8_t[::1] excluded_set,
    floating alpha,
    floating beta,
    floating gap,
    floating dual_norm_XtA,
    uint32_t n_features,
) noexcept nogil:
    """Apply gap safe screening for all features within enet_coordinate_descent_gram"""
    cdef floating d_j
    cdef floating Xj_theta
    cdef uint32_t n_active = 0
    # Due to floating point issues, gap might be negative.
    cdef floating radius = sqrt(2 * fabs(gap)) / alpha

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
            # Qw -= w[j] * Q[j]  # Update Qw = Q @ w
            _axpy(n_features, -w[j], &Q[j, 0], 1, &Qw[0], 1)
            w[j] = 0
            excluded_set[j] = 1

    return n_active


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
        warnings.warn(
            "Coordinate descent without L1 regularization may "
            "lead to unexpected results and is discouraged. "
            "Set l1_ratio > 0 to add L1 regularization."
        )

    if do_screening:
        active_set = np.empty(n_features, dtype=np.uint32)  # map [:n_active] -> j
        excluded_set = np.empty(n_features, dtype=np.uint8)

    with nogil:
        tol *= y_norm2

        # Check convergence before entering the main loop.
        gap, dual_norm_XtA = gap_enet_gram(
            n_features, w, alpha, beta, Qw, q, y_norm2, XtA, positive
        )
        if 0 <= gap <= tol:
            # Only if gap >=0 as singular Q may cause dubious values of gap.
            with gil:
                return np.asarray(w), gap, tol, 0

        # Gap Safe Screening Rules, see https://arxiv.org/abs/1802.07481, Eq. 11
        if do_screening:
            n_active = screen_features_enet_gram(
                Q=Q,
                XtA=XtA,
                w=w,
                Qw=Qw,
                active_set=active_set,
                excluded_set=excluded_set,
                alpha=alpha,
                beta=beta,
                gap=gap,
                dual_norm_XtA=dual_norm_XtA,
                n_features=n_features,
            )

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

            if w_max == 0.0 or d_w_max / w_max <= d_w_tol or n_iter == max_iter - 1:
                # the biggest coordinate update of this iteration was smaller than
                # the tolerance: check the duality gap as ultimate stopping
                # criterion
                gap, dual_norm_XtA = gap_enet_gram(
                    n_features, w, alpha, beta, Qw, q, y_norm2, XtA, positive
                )

                if gap <= tol:
                    # return if we reached desired tolerance
                    break

                # Gap Safe Screening Rules, see https://arxiv.org/abs/1802.07481, Eq. 11
                if do_screening:
                    n_active = screen_features_enet_gram(
                        Q=Q,
                        XtA=XtA,
                        w=w,
                        Qw=Qw,
                        active_set=active_set,
                        excluded_set=excluded_set,
                        alpha=alpha,
                        beta=beta,
                        gap=gap,
                        dual_norm_XtA=dual_norm_XtA,
                        n_features=n_features,
                    )

        else:
            # for/else, runs if for doesn't end with a `break`
            with gil:
                message = (
                    message_conv +
                    f" Duality gap: {gap:.3e}, tolerance: {tol:.3e}"
                )
                warnings.warn(message, ConvergenceWarning)

    return np.asarray(w), gap, tol, n_iter + 1


cdef (floating, floating) gap_enet_multi_task(
    int n_samples,
    int n_features,
    int n_tasks,
    const floating[::1, :] W,  # in
    floating l1_reg,
    floating l2_reg,
    const floating[::1, :] X,  # in
    const floating[::1, :] Y,  # in
    const floating[::1, :] R,  # in
    floating[:, ::1] XtA,  # out
    floating[::1] XtA_row_norms,  # out
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
        Inplace calculated as XtA = X.T @ R - l2_reg * W.T
    XtA_row_norms : memoryview of shape n_features
        Inplace calculated as np.sqrt(np.sum(XtA ** 2, axis=1))
    """
    cdef floating gap = 0.0
    cdef floating dual_norm_XtA
    cdef floating R_norm2
    cdef floating w_norm2 = 0.0
    cdef floating l21_norm
    cdef floating A_norm2
    cdef floating const_
    cdef unsigned int t, j

    # XtA = X.T @ R - l2_reg * W.T
    for j in range(n_features):
        for t in range(n_tasks):
            XtA[j, t] = _dot(n_samples, &X[0, j], 1, &R[0, t], 1) - l2_reg * W[t, j]

    # dual_norm_XtA = np.max(np.sqrt(np.sum(XtA ** 2, axis=1)))
    dual_norm_XtA = 0.0
    for j in range(n_features):
        # np.sqrt(np.sum(XtA ** 2, axis=1))
        XtA_row_norms[j] = _nrm2(n_tasks, &XtA[j, 0], 1)
        if XtA_row_norms[j] > dual_norm_XtA:
            dual_norm_XtA = XtA_row_norms[j]

    # R_norm2 = linalg.norm(R, ord="fro") ** 2
    R_norm2 = _dot(n_samples * n_tasks, &R[0, 0], 1, &R[0, 0], 1)

    # w_norm2 = linalg.norm(W, ord="fro") ** 2
    if l2_reg > 0:
        w_norm2 = _dot(n_features * n_tasks, &W[0, 0], 1, &W[0, 0], 1)

    if (dual_norm_XtA > l1_reg):
        const_ = l1_reg / dual_norm_XtA
        A_norm2 = R_norm2 * (const_ ** 2)
        gap = 0.5 * (R_norm2 + A_norm2)
    else:
        const_ = 1.0
        gap = R_norm2

    # l21_norm = np.sqrt(np.sum(W ** 2, axis=0)).sum()
    l21_norm = 0.0
    for ii in range(n_features):
        l21_norm += _nrm2(n_tasks, &W[0, ii], 1)

    gap += (
        l1_reg * l21_norm
        - const_ * _dot(n_samples * n_tasks, &R[0, 0], 1, &Y[0, 0], 1)  # np.sum(R * Y)
        + 0.5 * l2_reg * (1 + const_ ** 2) * w_norm2
    )
    return gap, dual_norm_XtA


def enet_coordinate_descent_multi_task(
    floating[::1, :] W,
    floating l1_reg,
    floating l2_reg,
    const floating[::1, :] X,
    const floating[::1, :] Y,
    unsigned int max_iter,
    floating tol,
    object rng,
    bint random=0,
    bint do_screening=1,
):
    """Cython version of the coordinate descent algorithm
        for Elastic-Net multi-task regression

        We minimize

        0.5 * norm(Y - X W.T, 2)^2 + l1_reg ||W.T||_21 + 0.5 * l2_reg norm(W.T, 2)^2

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

    if l1_reg == 0:
        warnings.warn(
            "Coordinate descent with l1_reg=0 may lead to unexpected"
            " results and is discouraged."
        )

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
            n_samples, n_features, n_tasks, W, l1_reg, l2_reg, X, Y, R, XtA, XtA_row_norms
        )
        if gap <= tol:
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
                Xj_theta = XtA_row_norms[j] / fmax(l1_reg, dual_norm_XtA)
                d_j = (1 - Xj_theta) / sqrt(norm2_cols_X[j] + l2_reg)
                if d_j <= sqrt(2 * gap) / l1_reg:
                    # include feature j
                    active_set[n_active] = j
                    excluded_set[j] = 0
                    n_active += 1
                else:
                    # R += W[:, 1] * X[:, 1][:, None]
                    for t in range(n_tasks):
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

                # W[:, j] = tmp * fmax(1. - l1_reg / nn, 0) / (norm2_cols_X[j] + l2_reg)
                _copy(n_tasks, &tmp[0], 1, &W[0, j], 1)
                _scal(n_tasks, fmax(1. - l1_reg / nn, 0) / (norm2_cols_X[j] + l2_reg),
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

            if w_max == 0.0 or d_w_max / w_max <= d_w_tol or n_iter == max_iter - 1:
                # the biggest coordinate update of this iteration was smaller than
                # the tolerance: check the duality gap as ultimate stopping
                # criterion
                gap, dual_norm_XtA = gap_enet_multi_task(
                    n_samples, n_features, n_tasks, W, l1_reg, l2_reg, X, Y, R, XtA, XtA_row_norms
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
                        Xj_theta = XtA_row_norms[j] / fmax(l1_reg, dual_norm_XtA)
                        d_j = (1 - Xj_theta) / sqrt(norm2_cols_X[j] + l2_reg)
                        if d_j <= sqrt(2 * gap) / l1_reg:
                            # include feature j
                            active_set[n_active] = j
                            excluded_set[j] = 0
                            n_active += 1
                        else:
                            # R += W[:, 1] * X[:, 1][:, None]
                            for t in range(n_tasks):
                                _axpy(n_samples, W[t, j], &X[0, j], 1, &R[0, t], 1)
                                W[t, j] = 0
                            excluded_set[j] = 1

        else:
            # for/else, runs if for doesn't end with a `break`
            with gil:
                message = (
                    message_conv +
                    f" Duality gap: {gap:.3e}, tolerance: {tol:.3e}"
                )
                warnings.warn(message, ConvergenceWarning)

    return np.asarray(W), gap, tol, n_iter + 1
