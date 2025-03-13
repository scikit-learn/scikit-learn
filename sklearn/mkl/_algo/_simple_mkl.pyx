import numpy as np

from ...utils._typedefs cimport float64_t, int32_t, int64_t, intp_t
from .._utils cimport BINARY, MULTICLASS, REGRESSION, ONECLASS
from .._utils cimport combine_kernels, fix_precision


def learn(
    object svm,
    const float64_t[:, :, ::1] kernels,
    object y,
    const double epsilon,
    const double tol=1e-8,
    const int max_iter=-1,
    const int verbose=1,
):
    cdef intp_t mu
    cdef double J, J_dagger, J_prev, gamma_max = 0.0
    cdef double goldensearch_precision = 1e-1, max_goldensearch_precision = 1e-08
    cdef str svm_name = svm.__class__.__name__
    cdef int64_t[::1] classes = np.unique(y) if svm_name == "SVC" else \
        np.empty(0, dtype=np.int64)
    cdef int n_classes = classes.shape[0]
    cdef int svm_type = BINARY if svm_name == "SVC" and n_classes == 2 else \
        MULTICLASS if svm_name == "SVC" else \
        REGRESSION if svm_name == "SVR" else \
        ONECLASS if svm_name == "OneClassSVM" \
        else -1
    cdef int32_t[::1] alpha_lengths
    cdef float64_t[:, ::1] alpha
    cdef float64_t[::1] delta_J, D, D_dagger, d_dagger
    cdef intp_t M = kernels.shape[0]
    cdef float64_t[::1] d = np.full(M, 1.0 / M, dtype=np.float64)

    if verbose:
        print(f"[SimpleMKL] Starting weights optimization for {M} kernels.")

    # Reduced gradient descent
    cdef int iterations = 1
    while True:
        old_d = np.copy(d)

        J = objective_value(svm, d, kernels, y)
        alpha, alpha_lengths = np.copy(svm.alpha_raw_), np.copy(svm.alpha_raw_lengths_)
        delta_J = gradient(svm_type, kernels, y, classes, alpha, alpha_lengths)

        mu = np.argmax(d)
        D = gradient_direction(d, delta_J, mu, M, tol)
        J_prev, J_dagger, d_dagger, D_dagger = J, float('-inf'), d, D

        while J_dagger < J:
            d, D = d_dagger, D_dagger

            if J_dagger != float('-inf'):
                J = J_dagger

            gamma_max = compute_gamma_max(d, D, M)
            d_dagger = np.add(d, np.multiply(gamma_max, D))
            J_dagger = objective_value(svm, d_dagger, kernels, y)
            if J_dagger < J:
                D_dagger = update_gradient_direction(D, d_dagger, mu, M, tol)

        # Line search for the optimal step size
        d, J, alpha, alpha_lengths = gamma_linesearch(
            svm, 0.0, gamma_max, gamma_max, J, J_dagger, d, D, kernels,
            J_prev, y, alpha, alpha_lengths, goldensearch_precision
        )
        d = fix_precision(d, tol)
        delta_J = gradient(svm_type, kernels, y, classes, alpha, alpha_lengths)

        # Precision enhancement for line search
        if max(abs(d-old_d))<tol and goldensearch_precision>max_goldensearch_precision:
            goldensearch_precision = goldensearch_precision/10

        # Stopping criterion (DualGap or KKT Constraint for Multi-Classification)
        stopping_crit, stopping_cond = stopping_criterion(
            svm_type, J, d, delta_J, y, alpha, M, getattr(svm, "epsilon", None), epsilon
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
            print("│ {:4d} │ {:10.4f} │ {:9.6f} │ {:11.6f} │".format(
                iterations, J, np.max(np.abs(np.subtract(d, old_d))), stopping_crit
            ))

        # Stopping conditions
        if iterations == max_iter or stopping_cond:
            J = objective_value(svm, d, kernels, y)  # It fits the final SVM

            if verbose:
                print("└──────┴────────────┴───────────┴─────────────┘")
                if iterations == max_iter:
                    print("Maximum number of iterations reached.")
                else:
                    print("Converged: stopping criterion reached.")
                print("Optimal weights:", d.base)

            return d.base, svm

        iterations += 1


def objective_value(
    object svm,
    const float64_t[::1] d,
    const float64_t[:, :, ::1] kernels,
    object y,
):
    if hasattr(svm, "alpha_raw_"):
        svm.alpha_init_ = svm.alpha_raw_
    svm.fit(combine_kernels(d, kernels), y)

    return svm.objective_val_


def gradient_direction(
    const float64_t[::1] d,
    const float64_t[::1] delta_J,
    const intp_t mu,
    const intp_t M,
    const double tol,
):
    cdef delta_J_norm = np.divide(delta_J, np.linalg.norm(delta_J))
    cdef float64_t[::1] reduced_grad = np.subtract(delta_J_norm, delta_J_norm[mu])

    cdef float64_t[::1] D = np.empty(M, dtype=np.float64)
    for m in range(M):
        if d[m] < tol and reduced_grad[m] > 0.0:
            D[m] = 0.0
        else:
            D[m] = -reduced_grad[m]
    D[mu] = -np.sum(D)

    return D


def compute_gamma_max(
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


def update_gradient_direction(
    float64_t[::1] D,
    const float64_t[::1] d,
    const intp_t mu,
    const intp_t M,
    const double tol,
):
    cdef intp_t m

    for m in range(M):
        if D[m] < 0.0 and d[m] < tol:
            D[m] = 0.0
    D[mu] = -np.sum(D[:mu]) - np.sum(D[mu+1:])

    return D


def gamma_linesearch(
    object svm,
    const double gamma_min,
    const double gamma_max,
    const double delta_max,
    const double cost_min,
    const double cost_max,
    const float64_t[::1] d,
    const float64_t[::1] D,
    const float64_t[:, :, ::1] kernels,
    const double J_prev,
    object y,
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
        cost_medr = objective_value(svm, tmp_d_r, kernels, y)
        alpha_r = np.copy(svm.alpha_raw_)
        alpha_lengths_r = np.copy(svm.alpha_raw_lengths_)

        tmp_d_l = np.add(d, np.multiply(gamma_medl, D))
        cost_medl = objective_value(svm, tmp_d_l, kernels, y)
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


def stopping_criterion(
    const int svm_type,
    const double J,
    const float64_t[::1] d,
    const float64_t[::1] delta_J,
    object y,
    const float64_t[:, ::1] alpha,
    const intp_t M,
    const double epsilon,
    const double cond_epsilon,
):
    if svm_type == MULTICLASS:
        return kkt_constraint(d, delta_J, M, cond_epsilon)
    else:
        return dual_gap(svm_type, J, delta_J, y, alpha, epsilon, cond_epsilon)


def dual_gap(
    const int svm_type,
    const double J,
    const float64_t[::1] delta_J,
    object y,
    const float64_t[:, ::1] alpha,
    const double epsilon,
    const double dual_gap_epsilon = 0.01,
):
    cdef double duality_gap

    if svm_type == BINARY:
        # (J(d) + max(-∂J/∂dₘ) - Σᵢαᵢ) / J(d)
        duality_gap = (J - np.min(delta_J) - np.sum(alpha)) / J
    elif svm_type == REGRESSION:
        # (J(d) + max(-∂J/∂dₘ) - Σᵢ(βᵢ - αᵢ)·yᵢ - ε·Σᵢ(βᵢ + αᵢ)) / J(d)
        # Here, alpha = [β, α]
        vec = np.concatenate([np.subtract(y, epsilon), np.subtract(-y, epsilon)])
        duality_gap = (J - np.min(delta_J) - np.dot(alpha, vec)) / J
    elif svm_type == ONECLASS:
        # (J(d) + max(-∂J/∂dₘ)) / J(d)
        duality_gap = (J - np.min(delta_J)) / J
    else:
        raise ValueError("Unknown SVM type.")

    return duality_gap, duality_gap < dual_gap_epsilon


def kkt_constraint(
    const float64_t[::1] d,
    const float64_t[::1] delta_J,
    const intp_t M,
    const double kkt_epsilon = 0.1,
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

    cdef double kkt_cons = np.abs((min_dJ_gt_0 - max_dJ_gt_0) / min_dJ_gt_0)
    return kkt_cons, kkt_cons <= kkt_epsilon and min_dJ_eq_0 >= max_dJ_gt_0


def gradient(
    object svm_type,
    const float64_t[:, :, ::1] kernels,
    object y,
    const int64_t[::1] classes,
    const float64_t[:, ::1] alpha,
    const int32_t[::1] alpha_lengths,
):
    if svm_type == BINARY or svm_type == MULTICLASS:
        return svc_gradient(kernels, y, classes, alpha, alpha_lengths)
    elif svm_type == REGRESSION:
        return svr_gradient(kernels, alpha)
    elif svm_type == ONECLASS:
        return one_class_gradient(kernels, alpha)
    else:
        raise ValueError("Unknown SVM type.")


def svc_gradient(
    const float64_t[:, :, ::1] kernels,
    const int64_t[::1] y,
    const int64_t[::1] classes,
    const float64_t[:, ::1] alpha,
    const int32_t[::1] alpha_lengths,
):
    cdef intp_t c1, c2, i, j, p, l, l_c1, idx_i, idx_j
    cdef double s
    cdef int yi_x_yj

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
    cdef float64_t[::1] delta_J = np.empty(kernels.shape[0], dtype=np.float64)
    for k in range(kernels.shape[0]):
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
                              kernels[k, indices[p, i], indices[p, j]])
                p += 1
        delta_J[k] = -0.5 * s

    return delta_J


def svr_gradient(
    const float64_t[:, :, ::1] kernels,
    const float64_t[:, ::1] alpha_raw_,
):
    cdef intp_t i, j, l
    cdef double s

    # ∂J/∂dₘ = -1/2·ΣᵢΣⱼ(βᵢ - αᵢ)·(βⱼ - αⱼ)·Kₘ(xᵢ, xⱼ)
    cdef float64_t[::1] delta_J = np.empty(kernels.shape[0], dtype=np.float64)
    l = kernels.shape[1]
    for k in range(kernels.shape[0]):
        s = 0.
        for i in range(l):
            for j in range(l):
                s += ((alpha_raw_[0, i] - alpha_raw_[0, i+l]) *
                      (alpha_raw_[0, j] - alpha_raw_[0, j+l]) *
                      kernels[k, i, j])
        delta_J[k] = -0.5 * s

    return delta_J


def one_class_gradient(
    const float64_t[:, :, ::1] kernels,
    const float64_t[:, ::1] alpha_raw_,
):
    cdef intp_t i, j, l
    cdef double s

    # ∂J/∂dₘ = -1/2·ΣᵢΣⱼαᵢ·αⱼ·Kₘ(xᵢ, xⱼ)
    cdef float64_t[::1] delta_J = np.empty(kernels.shape[0], dtype=np.float64)
    l = kernels.shape[1]
    for k in range(kernels.shape[0]):
        s = 0.
        for i in range(l):
            for j in range(l):
                s += (alpha_raw_[0, i] * alpha_raw_[0, j] * kernels[k, i, j])
        delta_J[k] = -0.5 * s

    return delta_J
