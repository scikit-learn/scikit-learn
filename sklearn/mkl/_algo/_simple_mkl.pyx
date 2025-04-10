"""
SimpleMKL algorithm.

Optimizes kernel weights for a Multiple Kernel Learning (MKL) model using
gradient-based methods.
"""

import numpy as np

from .._utils import kernel_generator

from ...utils._typedefs cimport float64_t, int32_t, int64_t, intp_t
from ._utils cimport BINARY, MULTICLASS, REGRESSION, ONECLASS
from ._utils cimport combine_kernels, fix_precision


def learn(
    object X,
    object y,
    object svm,
    object kernels,
    object kernels_scopes,
    object kernels_param_grids,
    const bint precomputed_kernels,
    const intp_t n_kernels,
    const intp_t n_samples,
    const double tol,
    const double numeric_tol=1e-8,
    const int verbose=0,
    const int max_iter=-1,
):
    classes = np.unique(y) if svm.__class__.__name__ == "SVC" else np.empty(0)

    cdef intp_t mu, M = n_kernels
    cdef double J, J_dagger, J_prev, gamma_max = 0.0, stopping_crit, \
        goldensearch_precision = 1e-1, max_goldensearch_precision = 1e-8
    cdef bint stopping_cond

    cdef str svm_name = svm.__class__.__name__
    cdef int n_classes = classes.shape[0]
    cdef int svm_type = BINARY if svm_name == "SVC" and n_classes == 2 else \
        MULTICLASS if svm_name == "SVC" else \
        REGRESSION if svm_name == "SVR" else \
        ONECLASS if svm_name == "OneClassSVM" \
        else -1

    cdef float64_t[:, ::1] alpha
    cdef int32_t[::1] alpha_lengths
    cdef float64_t[::1] delta_J, D, D_dagger, d_dagger
    cdef float64_t[::1] d = np.full(M, 1.0 / M, dtype=np.float64)

    if verbose:
        print(f"[SimpleMKL] Starting weights optimization for {M} kernels.")

    # Reduced gradient descent
    cdef int iterations = 1
    while True:
        old_d = np.copy(d)

        J = _objective_value(
            svm,
            d,
            kernel_generator(
                X,
                kernels,
                kernels_scopes,
                kernels_param_grids,
                precomputed_kernels,
            ),
            n_samples,
            y,
        )
        alpha, alpha_lengths = np.copy(svm.alpha_raw_), np.copy(svm.alpha_raw_lengths_)
        delta_J = _gradient(
            svm_type,
            kernel_generator(
                X,
                kernels,
                kernels_scopes,
                kernels_param_grids,
                precomputed_kernels,
            ),
            n_kernels,
            n_samples,
            y,
            classes,
            alpha,
            alpha_lengths,
        )

        mu = np.argmax(d)
        D = _gradient_direction(d, delta_J, mu, M, numeric_tol)
        J_prev, J_dagger, d_dagger, D_dagger = J, float('-inf'), d, D

        while J_dagger < J:
            d, D = d_dagger, D_dagger

            if J_dagger != float('-inf'):
                J = J_dagger

            gamma_max = _compute_gamma_max(d, D, M)
            d_dagger = np.add(d, np.multiply(gamma_max, D))
            J_dagger = _objective_value(
                svm,
                d_dagger,
                kernel_generator(
                    X,
                    kernels,
                    kernels_scopes,
                    kernels_param_grids,
                    precomputed_kernels,
                ),
                n_samples,
                y,
            )
            if J_dagger < J:
                D_dagger = _update_gradient_direction(D, d_dagger, mu, M, numeric_tol)

        # Line search for the optimal step size
        d, J, alpha, alpha_lengths = _gamma_linesearch(
            svm,
            0.0,
            gamma_max,
            gamma_max,
            J,
            J_dagger,
            d,
            D,
            X,
            kernels,
            kernels_scopes,
            kernels_param_grids,
            precomputed_kernels,
            y,
            n_samples,
            J_prev,
            alpha,
            alpha_lengths,
            goldensearch_precision,
        )
        d = fix_precision(d, numeric_tol)
        delta_J = _gradient(
            svm_type,
            kernel_generator(
                X,
                kernels,
                kernels_scopes,
                kernels_param_grids,
                precomputed_kernels,
            ),
            n_kernels,
            n_samples,
            y,
            classes,
            alpha,
            alpha_lengths,
        )

        # Precision enhancement for line search
        if (
            max(abs(d - old_d)) < numeric_tol
            and goldensearch_precision > max_goldensearch_precision
        ):
            goldensearch_precision = goldensearch_precision/10

        # Stopping criterion (DualGap or KKT Constraint for Multi-Classification)
        stopping_crit, stopping_cond = _stopping_criterion(
            svm_type, J, d, delta_J, y, alpha, M, getattr(svm, "epsilon", None), tol
        )

        # Verbose output
        if verbose:
            if iterations % 20 == 1:
                if iterations == 1:
                    print("┌──────┬────────────┬───────────┬─────────────┐")
                else:
                    print("├──────┼────────────┼───────────┼─────────────┤")
                print("│ Iter │    Obj.    │ DiffBetas │ Stop. crit. │")
                print("├──────┼────────────┼───────────┼─────────────┤")
            print("│ {:4d} │ {:10.3e} │ {:9.6f} │ {:11.6f} │".format(
                iterations, J, np.max(np.abs(np.subtract(d, old_d))), stopping_crit
            ))

        # Stopping conditions
        if iterations == max_iter or stopping_cond:
            # It fits the final SVM
            J = _objective_value(
                svm,
                d,
                kernel_generator(
                    X,
                    kernels,
                    kernels_scopes,
                    kernels_param_grids,
                    precomputed_kernels,
                ),
                n_samples,
                y,
            )

            # Final verbose output
            if verbose:
                print("└──────┴────────────┴───────────┴─────────────┘")
                if iterations == max_iter:
                    print("Maximum number of iterations reached.")
                else:
                    print("Converged: Stopping criterion reached.")
                print("Final objective value:", J)
                print("Optimal weights:", d.base)

            return d.base, svm, iterations

        iterations += 1


cdef double _objective_value(
    object svm,
    const float64_t[::1] d,
    object kernels,
    const intp_t n_samples,
    object y,
):
    if hasattr(svm, "alpha_raw_"):
        svm.alpha_init_ = svm.alpha_raw_
    svm.fit(combine_kernels(d, kernels, n_samples), y)

    return svm.objective_val_


cdef float64_t[::1] _gradient_direction(
    const float64_t[::1] d,
    const float64_t[::1] delta_J,
    const intp_t mu,
    const intp_t M,
    const double numeric_tol,
):
    cdef delta_J_norm = np.divide(delta_J, np.linalg.norm(delta_J))
    cdef float64_t[::1] reduced_grad = np.subtract(delta_J_norm, delta_J_norm[mu])

    cdef float64_t[::1] D = np.empty(M, dtype=np.float64)
    for m in range(M):
        if d[m] < numeric_tol and reduced_grad[m] > 0.0:
            D[m] = 0.0
        else:
            D[m] = -reduced_grad[m]
    D[mu] = -np.sum(D)

    return D


cdef double _compute_gamma_max(
    const float64_t[::1] d,
    const float64_t[::1] D,
    const intp_t M,
):
    cdef double gamma_max = float('inf')
    for m in range(M):
        f = -d[m]/D[m]
        if D[m] < 0.0 and f < gamma_max:
            gamma_max = f

    if gamma_max == float('inf'):
        return 0.0
    return gamma_max


cdef float64_t[::1] _update_gradient_direction(
    float64_t[::1] D,
    const float64_t[::1] d,
    const intp_t mu,
    const intp_t M,
    const double numeric_tol,
):
    cdef intp_t m

    for m in range(M):
        if D[m] < 0.0 and d[m] < numeric_tol:
            D[m] = 0.0
    D[mu] = -np.sum(D[:mu]) - np.sum(D[mu+1:])

    return D


def _gamma_linesearch(
    object svm,
    const double gamma_min,
    const double gamma_max,
    const double delta_max,
    const double cost_min,
    const double cost_max,
    const float64_t[::1] d,
    const float64_t[::1] D,
    object X,
    object kernels,
    object kernels_scopes,
    object kernels_param_grids,
    const bint precomputed_kernels,
    object y,
    const intp_t n_samples,
    const double J_prev,
    float64_t[:, ::1] alpha,
    int32_t[::1] alpha_lengths,
    const double goldensearch_precision_factor,
):
    cdef double gamma_medr, gamma_medl, cost_medr, cost_medl
    cdef float64_t[::1] tmp_d_r, tmp_d_l
    cdef float64_t[:, ::1] alpha_r, alpha_l
    cdef int32_t[::1] alpha_lengths_r, alpha_lengths_l

    cdef double gold_ratio = (np.sqrt(5) + 1) / 2
    cdef float64_t[::1] gamma_arr = np.array([gamma_min, gamma_max], dtype=np.float64)
    cdef float64_t[::1] cost_arr = np.array([cost_min, cost_max], dtype=np.float64)
    cdef intp_t coord = np.argmin(cost_arr)

    while (
        (gamma_max - gamma_min) > goldensearch_precision_factor * abs(delta_max)
        and gamma_max > np.finfo(float).eps
    ):
        gamma_medr = gamma_min + (gamma_max - gamma_min) / gold_ratio
        gamma_medl = gamma_min + (gamma_medr - gamma_min) / gold_ratio

        tmp_d_r = np.add(d, np.multiply(gamma_medr, D))
        cost_medr = _objective_value(
            svm,
            tmp_d_r,
            kernel_generator(
                X,
                kernels,
                kernels_scopes,
                kernels_param_grids,
                precomputed_kernels,
            ),
            n_samples,
            y,
        )
        alpha_r = np.copy(svm.alpha_raw_)
        alpha_lengths_r = np.copy(svm.alpha_raw_lengths_)

        tmp_d_l = np.add(d, np.multiply(gamma_medl, D))
        cost_medl = _objective_value(
            svm,
            tmp_d_l,
            kernel_generator(
                X,
                kernels,
                kernels_scopes,
                kernels_param_grids,
                precomputed_kernels,
            ),
            n_samples,
            y,
        )
        alpha_l = np.copy(svm.alpha_raw_)
        alpha_lengths_l = np.copy(svm.alpha_raw_lengths_)

        cost_arr = np.array([cost_min, cost_medl, cost_medr, cost_max])
        gamma_arr = np.array([gamma_min, gamma_medl, gamma_medr, gamma_max])
        coord = np.argmin(cost_arr)

        if coord == 0:
            gamma_max = gamma_medl
            cost_max = cost_medl
            alpha = alpha_l
            alpha_lengths = alpha_lengths_l
        elif coord == 1:
            gamma_max = gamma_medr
            cost_max = cost_medr
            alpha = alpha_r
            alpha_lengths = alpha_lengths_r
        elif coord == 2:
            gamma_min = gamma_medl
            cost_min = cost_medl
            alpha = alpha_l
            alpha_lengths = alpha_lengths_l
        elif coord == 3:
            gamma_min = gamma_medr
            cost_min = cost_medr
            alpha = alpha_r
            alpha_lengths = alpha_lengths_r

    if cost_arr[coord] < J_prev:
        d = np.add(d, np.multiply(gamma_arr[coord], D))
        return d, cost_arr[coord], alpha, alpha_lengths
    else:
        d = np.add(d, np.multiply(gamma_min, D))
        return d, cost_min, alpha, alpha_lengths


cpdef tuple _stopping_criterion(
    const int svm_type,
    const double J,
    const float64_t[::1] d,
    const float64_t[::1] delta_J,
    object y,
    const float64_t[:, ::1] alpha,
    const intp_t M,
    const double epsilon,
    const double tol,
):
    if svm_type == MULTICLASS:
        return _kkt_constraint(d, delta_J, M, tol)
    else:
        return _dual_gap(svm_type, J, delta_J, y, alpha, epsilon, tol)


cpdef tuple _dual_gap(
    const int svm_type,
    const double J,
    const float64_t[::1] delta_J,
    object y,
    const float64_t[:, ::1] alpha,
    const double epsilon,
    const double dual_gap_epsilon,
):
    cdef double duality_gap

    if svm_type == BINARY:
        # (J(d) + max(-∂J/∂dₘ) - Σᵢαᵢ) / J(d)
        duality_gap = (J - np.min(delta_J) - np.sum(alpha)) / J
    elif svm_type == REGRESSION:
        # (J(d) + max(-∂J/∂dₘ) - Σᵢ(βᵢ - αᵢ)·yᵢ - ε·Σᵢ(βᵢ + αᵢ)) / J(d)
        # Here, alpha = [β, α]
        vec = np.concatenate([np.subtract(y, epsilon), -np.add(y, epsilon)])
        duality_gap = (J - np.min(delta_J) - np.dot(alpha, vec).item()) / J
    elif svm_type == ONECLASS:
        # (J(d) + max(-∂J/∂dₘ)) / J(d)
        duality_gap = np.abs((J - np.min(delta_J)) / J)
    else:
        raise ValueError("Unknown SVM type.")

    return duality_gap, duality_gap < dual_gap_epsilon


cpdef tuple _kkt_constraint(
    const float64_t[::1] d,
    const float64_t[::1] delta_J,
    const intp_t M,
    const double kkt_epsilon,
):
    cdef intp_t m
    cdef double min_dJ_eq_0, min_dJ_gt_0, max_dJ_gt_0

    # We retrieve the extrema values of the gradient when dₘ = 0 and dₘ > 0
    min_dJ_eq_0, min_dJ_gt_0, max_dJ_gt_0 = float('inf'), float('inf'), float('-inf')
    for m in range(M):
        if d[m] > 0.0:
            if delta_J[m] < min_dJ_gt_0:
                min_dJ_gt_0 = delta_J[m]
            if delta_J[m] > max_dJ_gt_0:
                max_dJ_gt_0 = delta_J[m]
        else:
            if delta_J[m] < min_dJ_eq_0:
                min_dJ_eq_0 = delta_J[m]

    cdef double kkt_const = np.abs((min_dJ_gt_0 - max_dJ_gt_0) / min_dJ_gt_0)
    return kkt_const, kkt_const <= kkt_epsilon and min_dJ_eq_0 >= max_dJ_gt_0


cdef float64_t[::1] _gradient(
    object svm_type,
    object kernels,
    const intp_t n_kernels,
    const intp_t n_samples,
    object y,
    object classes,
    const float64_t[:, ::1] alpha,
    const int32_t[::1] alpha_lengths,
):
    if svm_type == BINARY or svm_type == MULTICLASS:
        return _svc_gradient(kernels, n_kernels, y, classes, alpha, alpha_lengths)
    elif svm_type == REGRESSION:
        return _svr_gradient(kernels, n_kernels, n_samples, alpha)
    elif svm_type == ONECLASS:
        return _one_class_gradient(kernels, n_kernels, n_samples, alpha)
    else:
        raise ValueError("Unknown SVM type.")


cdef float64_t[::1] _svc_gradient(
    object kernels,
    const intp_t n_kernels,
    object y,
    object classes,
    const float64_t[:, ::1] alpha,
    const int32_t[::1] alpha_lengths,
):
    cdef intp_t c1, c2, i, j, p, l, l_c1, idx_i, idx_j, k
    cdef double s
    cdef int yi_x_yj
    cdef const float64_t[:, ::1] kernel

    # LibSVM orders its x, y, and consequently alpha, by class
    # We need to retrieve the real indices in order to compute the kernel
    cdef int size = len(classes) * (len(classes) - 1) // 2
    cdef int64_t[:, ::1] indices = np.empty((size, np.max(alpha_lengths)),
                                            dtype=np.int64)
    cdef int64_t[::1] left_class_sizes = np.empty(size, dtype=np.int64)
    p, i = 0, 0
    for c1 in range(len(classes)):
        for c2 in range(c1 + 1, len(classes)):
            for i, idx_i in enumerate(np.where(np.equal(y, classes[c1]))[0]):
                indices[p, i] = idx_i
            left_class_sizes[p] = i + 1
            for j, idx_j in enumerate(np.where(np.equal(y, classes[c2]))[0], i+1):
                indices[p, j] = idx_j
            p += 1

    # ∂J/∂dₘ = -1/2·ΣₚΣᵢΣⱼαₚᵢ·αₚⱼ·yᵢ·yⱼ·Kₘ(xᵢ, xⱼ)
    cdef float64_t[::1] delta_J = np.empty(n_kernels, dtype=np.float64)
    for k, kernel in enumerate(kernels):
        p = 0
        s = 0.
        for c1 in range(len(classes)):
            for c2 in range(c1 + 1, len(classes)):
                l = alpha_lengths[p]
                l_c1 = left_class_sizes[p]
                for i in range(l):
                    for j in range(l):
                        yi_x_yj = 1 if (i < l_c1) == (j < l_c1) else -1
                        s += (alpha[p, i] * alpha[p, j] * yi_x_yj *
                              kernel[indices[p, i], indices[p, j]])
                p += 1
        delta_J[k] = -0.5 * s

    return delta_J


cdef float64_t[::1] _svr_gradient(
    object kernels,
    const intp_t n_kernels,
    const intp_t n_samples,
    const float64_t[:, ::1] alpha_raw_,
):
    cdef intp_t i, j, k
    cdef double s
    cdef const float64_t[:, ::1] kernel

    # ∂J/∂dₘ = -1/2·ΣᵢΣⱼ(βᵢ - αᵢ)·(βⱼ - αⱼ)·Kₘ(xᵢ, xⱼ)
    cdef float64_t[::1] delta_J = np.empty(n_kernels, dtype=np.float64)
    for k, kernel in enumerate(kernels):
        s = 0.
        for i in range(n_samples):
            for j in range(n_samples):
                s += ((alpha_raw_[0, i] - alpha_raw_[0, i+n_samples]) *
                      (alpha_raw_[0, j] - alpha_raw_[0, j+n_samples]) *
                      kernel[i, j])
        delta_J[k] = -0.5 * s

    return delta_J


cdef float64_t[::1] _one_class_gradient(
    object kernels,
    const intp_t n_kernels,
    const intp_t n_samples,
    const float64_t[:, ::1] alpha_raw_,
):
    cdef intp_t i, j, k
    cdef double s
    cdef const float64_t[:, ::1] kernel

    # ∂J/∂dₘ = -1/2·ΣᵢΣⱼαᵢ·αⱼ·Kₘ(xᵢ, xⱼ)
    cdef float64_t[::1] delta_J = np.empty(n_kernels, dtype=np.float64)
    for k, kernel in enumerate(kernels):
        s = 0.
        for i in range(n_samples):
            for j in range(n_samples):
                s += (alpha_raw_[0, i] * alpha_raw_[0, j] * kernel[i, j])
        delta_J[k] = -0.5 * s

    return delta_J
