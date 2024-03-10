# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#         Olivier Grisel <olivier.grisel@ensta.org>
#         Alexis Mignon <alexis.mignon@gmail.com>
#         Manoj Kumar <manojkumarsivaraj334@gmail.com>
#
# License: BSD 3 clause

from libc.math cimport fabs
cimport numpy as cnp
import numpy as np

from cython cimport floating
import warnings
from ..exceptions import ConvergenceWarning

from ..utils._cython_blas cimport (
    _axpy, _dot, _asum, _gemv, _nrm2, _copy, _scal
)
from ..utils._cython_blas cimport ColMajor, Trans, NoTrans
from ..utils._typedefs cimport uint32_t
from ..utils._random cimport our_rand_r


cnp.import_array()

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


cdef floating abs_max(int n, const floating* a) noexcept nogil:
    """np.max(np.abs(a))"""
    cdef int i
    cdef floating m = fabs(a[0])
    cdef floating d
    for i in range(1, n):
        d = fabs(a[i])
        if d > m:
            m = d
    return m


cdef floating max(int n, floating* a) noexcept nogil:
    """np.max(a)"""
    cdef int i
    cdef floating m = a[0]
    cdef floating d
    for i in range(1, n):
        d = a[i]
        if d > m:
            m = d
    return m


cdef floating diff_abs_max(int n, const floating* a, floating* b) noexcept nogil:
    """np.max(np.abs(a - b))"""
    cdef int i
    cdef floating m = fabs(a[0] - b[0])
    cdef floating d
    for i in range(1, n):
        d = fabs(a[i] - b[i])
        if d > m:
            m = d
    return m


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
    bint positive=0
):
    """Cython version of the coordinate descent algorithm
        for Elastic-Net regression

        We minimize

        (1/2) * norm(y - X w, 2)^2 + alpha norm(w, 1) + (beta/2) norm(w, 2)^2

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
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]

    # compute norms of the columns of X
    cdef floating[::1] norm_cols_X = np.square(X).sum(axis=0)

    # initial value of the residuals
    cdef floating[::1] R = np.empty(n_samples, dtype=dtype)
    cdef floating[::1] XtA = np.empty(n_features, dtype=dtype)

    cdef floating tmp
    cdef floating w_ii
    cdef floating d_w_max
    cdef floating w_max
    cdef floating d_w_ii
    cdef floating gap = tol + 1.0
    cdef floating d_w_tol = tol
    cdef floating dual_norm_XtA
    cdef floating R_norm2
    cdef floating w_norm2
    cdef floating l1_norm
    cdef floating const
    cdef floating A_norm2
    cdef unsigned int ii
    cdef unsigned int n_iter = 0
    cdef unsigned int f_iter
    cdef uint32_t rand_r_state_seed = rng.randint(0, RAND_R_MAX)
    cdef uint32_t* rand_r_state = &rand_r_state_seed

    if alpha == 0 and beta == 0:
        warnings.warn("Coordinate descent with no regularization may lead to "
                      "unexpected results and is discouraged.")

    with nogil:
        # R = y - np.dot(X, w)
        _copy(n_samples, &y[0], 1, &R[0], 1)
        _gemv(ColMajor, NoTrans, n_samples, n_features, -1.0, &X[0, 0],
              n_samples, &w[0], 1, 1.0, &R[0], 1)

        # tol *= np.dot(y, y)
        tol *= _dot(n_samples, &y[0], 1, &y[0], 1)

        for n_iter in range(max_iter):
            w_max = 0.0
            d_w_max = 0.0
            for f_iter in range(n_features):  # Loop over coordinates
                if random:
                    ii = rand_int(n_features, rand_r_state)
                else:
                    ii = f_iter

                if norm_cols_X[ii] == 0.0:
                    continue

                w_ii = w[ii]  # Store previous value

                if w_ii != 0.0:
                    # R += w_ii * X[:,ii]
                    _axpy(n_samples, w_ii, &X[0, ii], 1, &R[0], 1)

                # tmp = (X[:,ii]*R).sum()
                tmp = _dot(n_samples, &X[0, ii], 1, &R[0], 1)

                if positive and tmp < 0:
                    w[ii] = 0.0
                else:
                    w[ii] = (fsign(tmp) * fmax(fabs(tmp) - alpha, 0)
                             / (norm_cols_X[ii] + beta))

                if w[ii] != 0.0:
                    # R -=  w[ii] * X[:,ii] # Update residual
                    _axpy(n_samples, -w[ii], &X[0, ii], 1, &R[0], 1)

                # update the maximum absolute coefficient update
                d_w_ii = fabs(w[ii] - w_ii)
                d_w_max = fmax(d_w_max, d_w_ii)

                w_max = fmax(w_max, fabs(w[ii]))

            if (
                w_max == 0.0
                or d_w_max / w_max < d_w_tol
                or n_iter == max_iter - 1
            ):
                # the biggest coordinate update of this iteration was smaller
                # than the tolerance: check the duality gap as ultimate
                # stopping criterion

                # XtA = np.dot(X.T, R) - beta * w
                _copy(n_features, &w[0], 1, &XtA[0], 1)
                _gemv(ColMajor, Trans,
                      n_samples, n_features, 1.0, &X[0, 0], n_samples,
                      &R[0], 1,
                      -beta, &XtA[0], 1)

                if positive:
                    dual_norm_XtA = max(n_features, &XtA[0])
                else:
                    dual_norm_XtA = abs_max(n_features, &XtA[0])

                # R_norm2 = np.dot(R, R)
                R_norm2 = _dot(n_samples, &R[0], 1, &R[0], 1)

                # w_norm2 = np.dot(w, w)
                w_norm2 = _dot(n_features, &w[0], 1, &w[0], 1)

                if (dual_norm_XtA > alpha):
                    const = alpha / dual_norm_XtA
                    A_norm2 = R_norm2 * (const ** 2)
                    gap = 0.5 * (R_norm2 + A_norm2)
                else:
                    const = 1.0
                    gap = R_norm2

                l1_norm = _asum(n_features, &w[0], 1)

                # np.dot(R.T, y)
                gap += (alpha * l1_norm
                        - const * _dot(n_samples, &R[0], 1, &y[0], 1)
                        + 0.5 * beta * (1 + const ** 2) * (w_norm2))

                if gap < tol:
                    # return if we reached desired tolerance
                    break

        else:
            # for/else, runs if for doesn't end with a `break`
            with gil:
                message = (
                    "Objective did not converge. You might want to increase "
                    "the number of iterations, check the scale of the "
                    "features or consider increasing regularisation. "
                    f"Duality gap: {gap:.3e}, tolerance: {tol:.3e}"
                )
                if alpha < np.finfo(np.float64).eps:
                    message += (
                        " Linear regression models with null weight for the "
                        "l1 regularization term are more efficiently fitted "
                        "using one of the solvers implemented in "
                        "sklearn.linear_model.Ridge/RidgeCV instead."
                    )
                warnings.warn(message, ConvergenceWarning)

    return np.asarray(w), gap, tol, n_iter + 1


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
):
    """Cython version of the coordinate descent algorithm for Elastic-Net

    We minimize:

        1/2 * norm(y - Z w, 2)^2 + alpha * norm(w, 1) + (beta/2) * norm(w, 2)^2

    where Z = X - X_mean.
    With sample weights sw, this becomes

        1/2 * sum(sw * (y - Z w)^2, axis=0) + alpha * norm(w, 1)
        + (beta/2) * norm(w, 2)^2

    and X_mean is the weighted average of X (per column).

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
    #     yw = sample_weight
    #     R = sample_weight * residual
    #     norm_cols_X = np.sum(sample_weight * (X - X_mean)**2, axis=0)

    # get the data information into easy vars
    cdef unsigned int n_samples = y.shape[0]
    cdef unsigned int n_features = w.shape[0]

    # compute norms of the columns of X
    cdef unsigned int ii
    cdef floating[:] norm_cols_X

    cdef unsigned int startptr = X_indptr[0]
    cdef unsigned int endptr

    # initial value of the residuals
    # R = y - Zw, weighted version R = sample_weight * (y - Zw)
    cdef floating[::1] R
    cdef floating[::1] XtA
    cdef const floating[::1] yw

    if floating is float:
        dtype = np.float32
    else:
        dtype = np.float64

    norm_cols_X = np.zeros(n_features, dtype=dtype)
    XtA = np.zeros(n_features, dtype=dtype)

    cdef floating tmp
    cdef floating w_ii
    cdef floating d_w_max
    cdef floating w_max
    cdef floating d_w_ii
    cdef floating X_mean_ii
    cdef floating R_sum = 0.0
    cdef floating R_norm2
    cdef floating w_norm2
    cdef floating A_norm2
    cdef floating l1_norm
    cdef floating normalize_sum
    cdef floating gap = tol + 1.0
    cdef floating d_w_tol = tol
    cdef floating dual_norm_XtA
    cdef unsigned int jj
    cdef unsigned int n_iter = 0
    cdef unsigned int f_iter
    cdef uint32_t rand_r_state_seed = rng.randint(0, RAND_R_MAX)
    cdef uint32_t* rand_r_state = &rand_r_state_seed
    cdef bint center = False
    cdef bint no_sample_weights = sample_weight is None
    cdef int kk

    if no_sample_weights:
        yw = y
        R = y.copy()
    else:
        yw = np.multiply(sample_weight, y)
        R = yw.copy()

    with nogil:
        # center = (X_mean != 0).any()
        for ii in range(n_features):
            if X_mean[ii]:
                center = True
                break

        for ii in range(n_features):
            X_mean_ii = X_mean[ii]
            endptr = X_indptr[ii + 1]
            normalize_sum = 0.0
            w_ii = w[ii]

            if no_sample_weights:
                for jj in range(startptr, endptr):
                    normalize_sum += (X_data[jj] - X_mean_ii) ** 2
                    R[X_indices[jj]] -= X_data[jj] * w_ii
                norm_cols_X[ii] = normalize_sum + \
                    (n_samples - endptr + startptr) * X_mean_ii ** 2
                if center:
                    for jj in range(n_samples):
                        R[jj] += X_mean_ii * w_ii
            else:
                for jj in range(startptr, endptr):
                    tmp = sample_weight[X_indices[jj]]
                    # second term will be subtracted by loop over range(n_samples)
                    normalize_sum += (tmp * (X_data[jj] - X_mean_ii) ** 2
                                      - tmp * X_mean_ii ** 2)
                    R[X_indices[jj]] -= tmp * X_data[jj] * w_ii
                if center:
                    for jj in range(n_samples):
                        normalize_sum += sample_weight[jj] * X_mean_ii ** 2
                        R[jj] += sample_weight[jj] * X_mean_ii * w_ii
                norm_cols_X[ii] = normalize_sum
            startptr = endptr

        # tol *= np.dot(y, y)
        # with sample weights: tol *= y @ (sw * y)
        tol *= _dot(n_samples, &y[0], 1, &yw[0], 1)

        for n_iter in range(max_iter):

            w_max = 0.0
            d_w_max = 0.0

            for f_iter in range(n_features):  # Loop over coordinates
                if random:
                    ii = rand_int(n_features, rand_r_state)
                else:
                    ii = f_iter

                if norm_cols_X[ii] == 0.0:
                    continue

                startptr = X_indptr[ii]
                endptr = X_indptr[ii + 1]
                w_ii = w[ii]  # Store previous value
                X_mean_ii = X_mean[ii]

                if w_ii != 0.0:
                    # R += w_ii * X[:,ii]
                    if no_sample_weights:
                        for jj in range(startptr, endptr):
                            R[X_indices[jj]] += X_data[jj] * w_ii
                        if center:
                            for jj in range(n_samples):
                                R[jj] -= X_mean_ii * w_ii
                    else:
                        for jj in range(startptr, endptr):
                            tmp = sample_weight[X_indices[jj]]
                            R[X_indices[jj]] += tmp * X_data[jj] * w_ii
                        if center:
                            for jj in range(n_samples):
                                R[jj] -= sample_weight[jj] * X_mean_ii * w_ii

                # tmp = (X[:,ii] * R).sum()
                tmp = 0.0
                for jj in range(startptr, endptr):
                    tmp += R[X_indices[jj]] * X_data[jj]

                if center:
                    R_sum = 0.0
                    for jj in range(n_samples):
                        R_sum += R[jj]
                    tmp -= R_sum * X_mean_ii

                if positive and tmp < 0.0:
                    w[ii] = 0.0
                else:
                    w[ii] = fsign(tmp) * fmax(fabs(tmp) - alpha, 0) \
                            / (norm_cols_X[ii] + beta)

                if w[ii] != 0.0:
                    # R -=  w[ii] * X[:,ii] # Update residual
                    if no_sample_weights:
                        for jj in range(startptr, endptr):
                            R[X_indices[jj]] -= X_data[jj] * w[ii]
                        if center:
                            for jj in range(n_samples):
                                R[jj] += X_mean_ii * w[ii]
                    else:
                        for jj in range(startptr, endptr):
                            tmp = sample_weight[X_indices[jj]]
                            R[X_indices[jj]] -= tmp * X_data[jj] * w[ii]
                        if center:
                            for jj in range(n_samples):
                                R[jj] += sample_weight[jj] * X_mean_ii * w[ii]

                # update the maximum absolute coefficient update
                d_w_ii = fabs(w[ii] - w_ii)
                d_w_max = fmax(d_w_max, d_w_ii)

                w_max = fmax(w_max, fabs(w[ii]))

            if w_max == 0.0 or d_w_max / w_max < d_w_tol or n_iter == max_iter - 1:
                # the biggest coordinate update of this iteration was smaller than
                # the tolerance: check the duality gap as ultimate stopping
                # criterion

                # sparse X.T / dense R dot product
                if center:
                    R_sum = 0.0
                    for jj in range(n_samples):
                        R_sum += R[jj]

                # XtA = X.T @ R - beta * w
                for ii in range(n_features):
                    XtA[ii] = 0.0
                    for kk in range(X_indptr[ii], X_indptr[ii + 1]):
                        XtA[ii] += X_data[kk] * R[X_indices[kk]]

                    if center:
                        XtA[ii] -= X_mean[ii] * R_sum
                    XtA[ii] -= beta * w[ii]

                if positive:
                    dual_norm_XtA = max(n_features, &XtA[0])
                else:
                    dual_norm_XtA = abs_max(n_features, &XtA[0])

                # R_norm2 = np.dot(R, R)
                if no_sample_weights:
                    R_norm2 = _dot(n_samples, &R[0], 1, &R[0], 1)
                else:
                    R_norm2 = 0.0
                    for jj in range(n_samples):
                        # R is already multiplied by sample_weight
                        if sample_weight[jj] != 0:
                            R_norm2 += (R[jj] ** 2) / sample_weight[jj]

                # w_norm2 = np.dot(w, w)
                w_norm2 = _dot(n_features, &w[0], 1, &w[0], 1)
                if (dual_norm_XtA > alpha):
                    const = alpha / dual_norm_XtA
                    A_norm2 = R_norm2 * const**2
                    gap = 0.5 * (R_norm2 + A_norm2)
                else:
                    const = 1.0
                    gap = R_norm2

                l1_norm = _asum(n_features, &w[0], 1)

                gap += (alpha * l1_norm - const * _dot(
                            n_samples,
                            &R[0], 1,
                            &y[0], 1
                            )
                        + 0.5 * beta * (1 + const ** 2) * w_norm2)

                if gap < tol:
                    # return if we reached desired tolerance
                    break

        else:
            # for/else, runs if for doesn't end with a `break`
            with gil:
                warnings.warn("Objective did not converge. You might want to "
                              "increase the number of iterations. Duality "
                              "gap: {}, tolerance: {}".format(gap, tol),
                              ConvergenceWarning)

    return np.asarray(w), gap, tol, n_iter + 1


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
    bint positive=0
):
    """Cython version of the coordinate descent algorithm
        for Elastic-Net regression

        We minimize

        (1/2) * w^T Q w - q^T w + alpha norm(w, 1) + (beta/2) * norm(w, 2)^2

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
    cdef floating[:] H = np.dot(Q, w)

    cdef floating[:] XtA = np.zeros(n_features, dtype=dtype)
    cdef floating tmp
    cdef floating w_ii
    cdef floating d_w_max
    cdef floating w_max
    cdef floating d_w_ii
    cdef floating q_dot_w
    cdef floating w_norm2
    cdef floating gap = tol + 1.0
    cdef floating d_w_tol = tol
    cdef floating dual_norm_XtA
    cdef unsigned int ii
    cdef unsigned int n_iter = 0
    cdef unsigned int f_iter
    cdef uint32_t rand_r_state_seed = rng.randint(0, RAND_R_MAX)
    cdef uint32_t* rand_r_state = &rand_r_state_seed

    cdef floating y_norm2 = np.dot(y, y)
    cdef floating* w_ptr = &w[0]
    cdef const floating* Q_ptr = &Q[0, 0]
    cdef const floating* q_ptr = &q[0]
    cdef floating* H_ptr = &H[0]
    cdef floating* XtA_ptr = &XtA[0]
    tol = tol * y_norm2

    if alpha == 0:
        warnings.warn(
            "Coordinate descent without L1 regularization may "
            "lead to unexpected results and is discouraged. "
            "Set l1_ratio > 0 to add L1 regularization."
        )

    with nogil:
        for n_iter in range(max_iter):
            w_max = 0.0
            d_w_max = 0.0
            for f_iter in range(n_features):  # Loop over coordinates
                if random:
                    ii = rand_int(n_features, rand_r_state)
                else:
                    ii = f_iter

                if Q[ii, ii] == 0.0:
                    continue

                w_ii = w[ii]  # Store previous value

                if w_ii != 0.0:
                    # H -= w_ii * Q[ii]
                    _axpy(n_features, -w_ii, Q_ptr + ii * n_features, 1,
                          H_ptr, 1)

                tmp = q[ii] - H[ii]

                if positive and tmp < 0:
                    w[ii] = 0.0
                else:
                    w[ii] = fsign(tmp) * fmax(fabs(tmp) - alpha, 0) \
                        / (Q[ii, ii] + beta)

                if w[ii] != 0.0:
                    # H +=  w[ii] * Q[ii] # Update H = X.T X w
                    _axpy(n_features, w[ii], Q_ptr + ii * n_features, 1,
                          H_ptr, 1)

                # update the maximum absolute coefficient update
                d_w_ii = fabs(w[ii] - w_ii)
                if d_w_ii > d_w_max:
                    d_w_max = d_w_ii

                if fabs(w[ii]) > w_max:
                    w_max = fabs(w[ii])

            if w_max == 0.0 or d_w_max / w_max < d_w_tol or n_iter == max_iter - 1:
                # the biggest coordinate update of this iteration was smaller than
                # the tolerance: check the duality gap as ultimate stopping
                # criterion

                # q_dot_w = np.dot(w, q)
                q_dot_w = _dot(n_features, w_ptr, 1, q_ptr, 1)

                for ii in range(n_features):
                    XtA[ii] = q[ii] - H[ii] - beta * w[ii]
                if positive:
                    dual_norm_XtA = max(n_features, XtA_ptr)
                else:
                    dual_norm_XtA = abs_max(n_features, XtA_ptr)

                # temp = np.sum(w * H)
                tmp = 0.0
                for ii in range(n_features):
                    tmp += w[ii] * H[ii]
                R_norm2 = y_norm2 + tmp - 2.0 * q_dot_w

                # w_norm2 = np.dot(w, w)
                w_norm2 = _dot(n_features, &w[0], 1, &w[0], 1)

                if (dual_norm_XtA > alpha):
                    const = alpha / dual_norm_XtA
                    A_norm2 = R_norm2 * (const ** 2)
                    gap = 0.5 * (R_norm2 + A_norm2)
                else:
                    const = 1.0
                    gap = R_norm2

                # The call to asum is equivalent to the L1 norm of w
                gap += (
                    alpha * _asum(n_features, &w[0], 1)
                    - const * y_norm2
                    + const * q_dot_w
                    + 0.5 * beta * (1 + const ** 2) * w_norm2
                )

                if gap < tol:
                    # return if we reached desired tolerance
                    break

        else:
            # for/else, runs if for doesn't end with a `break`
            with gil:
                warnings.warn("Objective did not converge. You might want to "
                              "increase the number of iterations. Duality "
                              "gap: {}, tolerance: {}".format(gap, tol),
                              ConvergenceWarning)

    return np.asarray(w), gap, tol, n_iter + 1


def enet_coordinate_descent_multi_task(
    const floating[::1, :] W,
    floating l1_reg,
    floating l2_reg,
    const floating[::1, :] X,
    const floating[::1, :] Y,
    unsigned int max_iter,
    floating tol,
    object rng,
    bint random=0
):
    """Cython version of the coordinate descent algorithm
        for Elastic-Net mult-task regression

        We minimize

        0.5 * norm(Y - X W.T, 2)^2 + l1_reg ||W.T||_21 + 0.5 * l2_reg norm(W.T, 2)^2

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

    # to store XtA
    cdef floating[:, ::1] XtA = np.zeros((n_features, n_tasks), dtype=dtype)
    cdef floating XtA_axis1norm
    cdef floating dual_norm_XtA

    # initial value of the residuals
    cdef floating[::1, :] R = np.zeros((n_samples, n_tasks), dtype=dtype, order='F')

    cdef floating[::1] norm_cols_X = np.zeros(n_features, dtype=dtype)
    cdef floating[::1] tmp = np.zeros(n_tasks, dtype=dtype)
    cdef floating[::1] w_ii = np.zeros(n_tasks, dtype=dtype)
    cdef floating d_w_max
    cdef floating w_max
    cdef floating d_w_ii
    cdef floating nn
    cdef floating W_ii_abs_max
    cdef floating gap = tol + 1.0
    cdef floating d_w_tol = tol
    cdef floating R_norm
    cdef floating w_norm
    cdef floating ry_sum
    cdef floating l21_norm
    cdef unsigned int ii
    cdef unsigned int jj
    cdef unsigned int n_iter = 0
    cdef unsigned int f_iter
    cdef uint32_t rand_r_state_seed = rng.randint(0, RAND_R_MAX)
    cdef uint32_t* rand_r_state = &rand_r_state_seed

    cdef const floating* X_ptr = &X[0, 0]
    cdef const floating* Y_ptr = &Y[0, 0]

    if l1_reg == 0:
        warnings.warn(
            "Coordinate descent with l1_reg=0 may lead to unexpected"
            " results and is discouraged."
        )

    with nogil:
        # norm_cols_X = (np.asarray(X) ** 2).sum(axis=0)
        for ii in range(n_features):
            norm_cols_X[ii] = _nrm2(n_samples, X_ptr + ii * n_samples, 1) ** 2

        # R = Y - np.dot(X, W.T)
        _copy(n_samples * n_tasks, Y_ptr, 1, &R[0, 0], 1)
        for ii in range(n_features):
            for jj in range(n_tasks):
                if W[jj, ii] != 0:
                    _axpy(n_samples, -W[jj, ii], X_ptr + ii * n_samples, 1,
                          &R[0, jj], 1)

        # tol = tol * linalg.norm(Y, ord='fro') ** 2
        tol = tol * _nrm2(n_samples * n_tasks, Y_ptr, 1) ** 2

        for n_iter in range(max_iter):
            w_max = 0.0
            d_w_max = 0.0
            for f_iter in range(n_features):  # Loop over coordinates
                if random:
                    ii = rand_int(n_features, rand_r_state)
                else:
                    ii = f_iter

                if norm_cols_X[ii] == 0.0:
                    continue

                # w_ii = W[:, ii] # Store previous value
                _copy(n_tasks, &W[0, ii], 1, &w_ii[0], 1)

                # Using Numpy:
                # R += np.dot(X[:, ii][:, None], w_ii[None, :]) # rank 1 update
                # Using Blas Level2:
                # _ger(RowMajor, n_samples, n_tasks, 1.0,
                #      &X[0, ii], 1,
                #      &w_ii[0], 1, &R[0, 0], n_tasks)
                # Using Blas Level1 and for loop to avoid slower threads
                # for such small vectors
                for jj in range(n_tasks):
                    if w_ii[jj] != 0:
                        _axpy(n_samples, w_ii[jj], X_ptr + ii * n_samples, 1,
                              &R[0, jj], 1)

                # Using numpy:
                # tmp = np.dot(X[:, ii][None, :], R).ravel()
                # Using BLAS Level 2:
                # _gemv(RowMajor, Trans, n_samples, n_tasks, 1.0, &R[0, 0],
                #       n_tasks, &X[0, ii], 1, 0.0, &tmp[0], 1)
                # Using BLAS Level 1 (faster for small vectors like here):
                for jj in range(n_tasks):
                    tmp[jj] = _dot(n_samples, X_ptr + ii * n_samples, 1,
                                   &R[0, jj], 1)

                # nn = sqrt(np.sum(tmp ** 2))
                nn = _nrm2(n_tasks, &tmp[0], 1)

                # W[:, ii] = tmp * fmax(1. - l1_reg / nn, 0) / (norm_cols_X[ii] + l2_reg)
                _copy(n_tasks, &tmp[0], 1, &W[0, ii], 1)
                _scal(n_tasks, fmax(1. - l1_reg / nn, 0) / (norm_cols_X[ii] + l2_reg),
                      &W[0, ii], 1)

                # Using numpy:
                # R -= np.dot(X[:, ii][:, None], W[:, ii][None, :])
                # Using BLAS Level 2:
                # Update residual : rank 1 update
                # _ger(RowMajor, n_samples, n_tasks, -1.0,
                #      &X[0, ii], 1, &W[0, ii], 1,
                #      &R[0, 0], n_tasks)
                # Using BLAS Level 1 (faster for small vectors like here):
                for jj in range(n_tasks):
                    if W[jj, ii] != 0:
                        _axpy(n_samples, -W[jj, ii], X_ptr + ii * n_samples, 1,
                              &R[0, jj], 1)

                # update the maximum absolute coefficient update
                d_w_ii = diff_abs_max(n_tasks, &W[0, ii], &w_ii[0])

                if d_w_ii > d_w_max:
                    d_w_max = d_w_ii

                W_ii_abs_max = abs_max(n_tasks, &W[0, ii])
                if W_ii_abs_max > w_max:
                    w_max = W_ii_abs_max

            if w_max == 0.0 or d_w_max / w_max < d_w_tol or n_iter == max_iter - 1:
                # the biggest coordinate update of this iteration was smaller than
                # the tolerance: check the duality gap as ultimate stopping
                # criterion

                # XtA = np.dot(X.T, R) - l2_reg * W.T
                for ii in range(n_features):
                    for jj in range(n_tasks):
                        XtA[ii, jj] = _dot(
                            n_samples, X_ptr + ii * n_samples, 1, &R[0, jj], 1
                            ) - l2_reg * W[jj, ii]

                # dual_norm_XtA = np.max(np.sqrt(np.sum(XtA ** 2, axis=1)))
                dual_norm_XtA = 0.0
                for ii in range(n_features):
                    # np.sqrt(np.sum(XtA ** 2, axis=1))
                    XtA_axis1norm = _nrm2(n_tasks, &XtA[ii, 0], 1)
                    if XtA_axis1norm > dual_norm_XtA:
                        dual_norm_XtA = XtA_axis1norm

                # TODO: use squared L2 norm directly
                # R_norm = linalg.norm(R, ord='fro')
                # w_norm = linalg.norm(W, ord='fro')
                R_norm = _nrm2(n_samples * n_tasks, &R[0, 0], 1)
                w_norm = _nrm2(n_features * n_tasks, &W[0, 0], 1)
                if (dual_norm_XtA > l1_reg):
                    const = l1_reg / dual_norm_XtA
                    A_norm = R_norm * const
                    gap = 0.5 * (R_norm ** 2 + A_norm ** 2)
                else:
                    const = 1.0
                    gap = R_norm ** 2

                # ry_sum = np.sum(R * y)
                ry_sum = _dot(n_samples * n_tasks, &R[0, 0], 1, &Y[0, 0], 1)

                # l21_norm = np.sqrt(np.sum(W ** 2, axis=0)).sum()
                l21_norm = 0.0
                for ii in range(n_features):
                    l21_norm += _nrm2(n_tasks, &W[0, ii], 1)

                gap += (
                    l1_reg * l21_norm
                    - const * ry_sum
                    + 0.5 * l2_reg * (1 + const ** 2) * (w_norm ** 2)
                )

                if gap < tol:
                    # return if we reached desired tolerance
                    break
        else:
            # for/else, runs if for doesn't end with a `break`
            with gil:
                warnings.warn("Objective did not converge. You might want to "
                              "increase the number of iterations. Duality "
                              "gap: {}, tolerance: {}".format(gap, tol),
                              ConvergenceWarning)

    return np.asarray(W), gap, tol, n_iter + 1
