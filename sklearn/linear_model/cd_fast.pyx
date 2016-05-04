# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#         Olivier Grisel <olivier.grisel@ensta.org>
#         Alexis Mignon <alexis.mignon@gmail.com>
#         Manoj Kumar <manojkumarsivaraj334@gmail.com>
#         Olivier Fercoq <olivier.fercoq@telecom-paristech.fr>
#         Joseph Salmon <joseph.salmon@telecom-paristech.fr>
#
# Licence: BSD 3 clause

from libc.math cimport fabs, sqrt
cimport numpy as np
import numpy as np
import numpy.linalg as linalg

cimport cython
from cpython cimport bool
import warnings

ctypedef np.float64_t DOUBLE
ctypedef np.uint32_t UINT32_t

np.import_array()

# The following two functions are shamelessly copied from the tree code.

cdef enum:
    # Max value for our rand_r replacement (near the bottom).
    # We don't use RAND_MAX because it's different across platforms and
    # particularly tiny on Windows/MSVC.
    RAND_R_MAX = 0x7FFFFFFF


cdef inline UINT32_t our_rand_r(UINT32_t* seed) nogil:
    seed[0] ^= <UINT32_t>(seed[0] << 13)
    seed[0] ^= <UINT32_t>(seed[0] >> 17)
    seed[0] ^= <UINT32_t>(seed[0] << 5)

    return seed[0] % (<UINT32_t>RAND_R_MAX + 1)


cdef inline UINT32_t rand_int(UINT32_t end, UINT32_t* random_state) nogil:
    """Generate a random integer in [0; end)."""
    return our_rand_r(random_state) % end


cdef inline double fmax(double x, double y) nogil:
    if x > y:
        return x
    return y


cdef inline double fmin(double x, double y) nogil:
    if x < y:
        return x
    return y


cdef inline double fsign(double f) nogil:
    if f == 0:
        return 0
    elif f > 0:
        return 1.0
    else:
        return -1.0


cdef double abs_max(int n, double* a) nogil:
    """np.max(np.abs(a))"""
    cdef int i
    cdef double m = fabs(a[0])
    cdef double d
    for i in range(1, n):
        d = fabs(a[i])
        if d > m:
            m = d
    return m


cdef double max(int n, double* a) nogil:
    """np.max(a)"""
    cdef int i
    cdef double m = a[0]
    cdef double d
    for i in range(1, n):
        d = a[i]
        if d > m:
            m = d
    return m


cdef double diff_abs_max(int n, double* a, double* b) nogil:
    """np.max(np.abs(a - b))"""
    cdef int i
    cdef double m = fabs(a[0] - b[0])
    cdef double d
    for i in range(1, n):
        d = fabs(a[i] - b[i])
        if d > m:
            m = d
    return m


cdef extern from "cblas.h":
    enum CBLAS_ORDER:
        CblasRowMajor=101
        CblasColMajor=102
    enum CBLAS_TRANSPOSE:
        CblasNoTrans=111
        CblasTrans=112
        CblasConjTrans=113
        AtlasConj=114

    void daxpy "cblas_daxpy"(int N, double alpha, double *X, int incX,
                             double *Y, int incY) nogil
    double ddot "cblas_ddot"(int N, double *X, int incX, double *Y, int incY
                             ) nogil
    double dasum "cblas_dasum"(int N, double *X, int incX) nogil
    void dger "cblas_dger"(CBLAS_ORDER Order, int M, int N, double alpha,
                double *X, int incX, double *Y, int incY, double *A, int lda) nogil
    void dgemv "cblas_dgemv"(CBLAS_ORDER Order,
                      CBLAS_TRANSPOSE TransA, int M, int N,
                      double alpha, double *A, int lda,
                      double *X, int incX, double beta,
                      double *Y, int incY) nogil
    double dnrm2 "cblas_dnrm2"(int N, double *X, int incX) nogil
    void dcopy "cblas_dcopy"(int N, double *X, int incX, double *Y, int incY) nogil
    void dscal "cblas_dscal"(int N, double alpha, double *X, int incX) nogil



# Function to compute the duality gap
cdef double enet_duality_gap( # Data
                             unsigned int n_samples,
                             unsigned int n_features,
                             unsigned int n_tasks,
                             DOUBLE * X_data,
                             DOUBLE * y_data,
                             DOUBLE * R_data,
                             DOUBLE * w_data,
                             # Variables intended to be modified
                             DOUBLE * XtA_data,
                             double * dual_scaling,
                             # Parameters
                             double alpha,
                             double beta,
                             bint positive,
                             # active set for improved performance
                             np.int32_t* disabled_data = NULL,
                             unsigned int n_disabled = 0
                             ) nogil:

    cdef double R_norm2
    cdef double w_norm2
    cdef double y_norm2
    cdef double l1_norm
    cdef double const
    cdef double gap
    cdef unsigned int i
    cdef double yTA
    cdef double dual_norm_XtA

    # XtA = np.dot(X.T, R) - beta * w
    for i in range(n_features):
        if (n_disabled == 0 or disabled_data[i] == 0):
            XtA_data[i] = ddot(
                n_samples,
                X_data + i * n_samples,
                1, R_data, 1) - beta * w_data[i]
        else:
            # We only need the infinity norm of XtA and by KKT, we know that
            # in the present case, XtA[i] will not reach the maximum
            XtA_data[i] = 0

    yTA = ddot(n_samples, y_data, 1, R_data, 1)
    # R_norm2 = np.dot(R, R)
    R_norm2 = ddot(n_samples, R_data, 1, R_data, 1)

    if alpha == 0:
        dual_scaling[0] = 1
    else:
        if positive:
            dual_norm_XtA = max(n_features, XtA_data)
        else:
            dual_norm_XtA = abs_max(n_features, XtA_data)

        if dual_norm_XtA <= 0:
            if R_norm2 == 0:
                dual_scaling[0] = 1. / alpha
            else:
                dual_scaling[0] = yTA / R_norm2 / alpha
        elif positive:
            dual_scaling[0] = fmin(yTA / (alpha * R_norm2),
                                 1. / dual_norm_XtA)
        else:
            dual_scaling[0] = fmin(fmax(yTA / (alpha * R_norm2), -1. / dual_norm_XtA),
                               1. / dual_norm_XtA)

    dual_scaling[0] = 1. / dual_scaling[0]

    # w_norm2 = np.dot(w, w)
    w_norm2 = ddot(n_features, w_data, 1, w_data, 1)

    const = alpha / dual_scaling[0]
    A_norm2 = R_norm2 * (const ** 2)
    gap = 0.5 * (R_norm2 + A_norm2)

    l1_norm = dasum(n_features, w_data, 1)

    # np.dot(R.T, y)
    gap += (alpha * l1_norm - const * yTA
            + 0.5 * beta * (1. + const ** 2) * (w_norm2))

    return gap


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def enet_coordinate_descent(np.ndarray[DOUBLE, ndim=1] w,
                            double alpha, double beta,
                            np.ndarray[DOUBLE, ndim=2, mode='fortran'] X,
                            np.ndarray[DOUBLE, ndim=1, mode='c'] y,
                            int max_iter, double tol,
                            object rng, bint random=0,
                            bint positive=0,
                            int screening=10
                            ):
    """Cython version of the coordinate descent algorithm
        for Elastic-Net regression

        We minimize

        (1/2) * norm(y - X w, 2)^2 + alpha norm(w, 1) + (beta/2) norm(w, 2)^2

    """
    # get the data information into easy vars
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]

    # get the number of tasks indirectly, using strides
    cdef unsigned int n_tasks = y.strides[0] / sizeof(DOUBLE)

    # compute norms of the columns of X
    cdef np.ndarray[DOUBLE, ndim=1] norm2_cols_X = (X**2).sum(axis=0)
    cdef np.ndarray[DOUBLE, ndim=1] norm_cols_X = np.sqrt(norm2_cols_X)

    # initial value of the residuals
    cdef np.ndarray[DOUBLE, ndim=1] R = np.empty(n_samples)

    cdef np.ndarray[DOUBLE, ndim=1] XtA = np.empty(n_features)
    cdef np.ndarray[DOUBLE, ndim=1] Xty = np.empty(n_features)
    cdef np.ndarray[np.int32_t, ndim=1] active_set = np.empty(n_features, dtype=np.intc)
    cdef np.ndarray[np.int32_t, ndim=1] disabled = np.zeros(n_features, dtype=np.intc)
    #
    cdef double tmp
    cdef double w_ii
    cdef double d_w_max
    cdef double w_max
    cdef double d_w_ii
    cdef double gap = tol + 1.0
    cdef double d_w_tol = tol
    cdef np.ndarray[DOUBLE, ndim=1] dual_scaling = np.empty(1)
    cdef double r_screening
    cdef double tmp_XtA_dual_scaling
    cdef double Xtymax
    cdef unsigned int n_active = n_features
    cdef unsigned int ii
    cdef unsigned int i
    cdef unsigned int n_iter = 0
    cdef unsigned int f_iter
    cdef unsigned int n_features_dual_gap
    cdef UINT32_t rand_r_state_seed = rng.randint(0, RAND_R_MAX)
    cdef UINT32_t* rand_r_state = &rand_r_state_seed
    cdef bint do_gap = 0

    if alpha == 0:
        warnings.warn("Coordinate descent with alpha=0 may lead to unexpected"
                      " results and is discouraged.")

    y_norm2 = ddot(n_samples, <DOUBLE*>y.data, 1, <DOUBLE*>y.data, 1)
    cdef unsigned int done = False

    with nogil:
        if screening == 0:
            n_active = 0
            for ii in range(n_features):
                if norm2_cols_X[ii] != 0.0:
                    active_set[n_active] = ii
                    n_active += 1
                else:
                    disabled[ii] = 1

        # R = y - np.dot(X, w)
        for j in range(n_samples):
            R[j] = y[j]
        for i in range(n_features):
            # R -=  w[ii] * X[:,ii] # Update residual
            daxpy(n_samples, -w[i],
                  <DOUBLE*>(X.data + i * n_samples * sizeof(DOUBLE)),
                  1, <DOUBLE*>R.data, 1)

        ### scaling tolerance
        # tol *= np.dot(y, y)
        tol *= y_norm2

        ### main loop
        for n_iter in range(max_iter):
            ####################
            # variable screening

            do_gap = False
            if (n_iter and w_max == 0.0  # heuristic termination criterion
                    or d_w_max / w_max < d_w_tol):
                # the biggest coordinate update of this iteration was smaller
                # than the tolerance: check the duality gap as ultimate
                # stopping criterion
                do_gap = True

            do_gap = (do_gap or (screening > 0 and n_iter % screening == 0)
                      or (n_iter == max_iter - 1)) or (n_iter == 0)

            if do_gap:  # Screening
                if screening > 0:
                    n_features_dual_gap = n_features - n_active
                else:
                    n_features_dual_gap = 0

                gap = enet_duality_gap(n_samples, n_features, n_tasks,
                                       <DOUBLE*>X.data, <DOUBLE*>y.data, <DOUBLE*>R.data,
                                       <DOUBLE*>w.data, <DOUBLE*>XtA.data,
                                       <DOUBLE*>dual_scaling.data, alpha, beta, positive,
                                       <np.int32_t*>disabled.data, n_features_dual_gap)

                # break if we reached desired tolerance
                if gap < tol:
                    break

            if do_gap and screening > 0:  # Screening

                r_screening = sqrt(2. * gap) / alpha
                n_active = 0

                for ii in range(n_features):  # Loop over coordinates

                    if disabled[ii] == 1:
                         continue

                    tmp = XtA[ii] / dual_scaling[0]  # ours

                    if not positive:
                        tmp = fabs(tmp)

                    if tmp >= (1. - r_screening * norm_cols_X[ii]) or w[ii] != 0.:
                        active_set[n_active] = ii
                        n_active += 1
                    else:
                        disabled[ii] = 1

            ###################
            # Coordinate descent
            w_max = 0.0
            d_w_max = 0.0
            for f_iter in range(n_active):  # Loop over coordinates
                if random:
                    ii = active_set[rand_int(n_active, rand_r_state)]
                else:
                    ii = active_set[f_iter]

                if norm2_cols_X[ii] == 0.0:
                    continue

                w_ii = w[ii]  # Store previous value

                # tmp = (X[:,ii]*R).sum() + norm2_cols_X[ii] * w_ii
                tmp = ddot(n_samples,
                           <DOUBLE*>(X.data + ii * n_samples * sizeof(DOUBLE)),
                           1, <DOUBLE*>R.data, 1) + norm2_cols_X[ii] * w_ii

                if positive and tmp < 0:
                    w[ii] = 0.0
                else:
                    w[ii] = (fsign(tmp) * fmax(fabs(tmp) - alpha, 0)
                             / (norm2_cols_X[ii] + beta))

                d_w_ii = w[ii] - w_ii
                if w[ii] != w_ii:
                    # R -=  (w[ii]-w_ii) * X[:,ii] # Update residual
                    daxpy(n_samples, -d_w_ii,
                          <DOUBLE*>(X.data + ii * n_samples * sizeof(DOUBLE)),
                          1, <DOUBLE*>R.data, 1)

                # update the maximum absolute coefficient update
                d_w_ii = fabs(d_w_ii)
                if d_w_ii > d_w_max:
                    d_w_max = d_w_ii

                if fabs(w[ii]) > w_max:
                    w_max = fabs(w[ii])

    return w, gap, tol, n_iter + 1


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef double sparse_enet_duality_gap(unsigned int n_samples,
                       unsigned int n_features,
                       double[:] X_data,
                       int[:] X_indices,
                       int[:] X_indptr,
                       double[:] X_mean,
                       double[:] y,
                       double y_sum,
                       double[:] R,
                       double R_shift,
                       double[:] w,
                       # Variables intended to be modified
                       double[:] XtA,
                       double[:] X_T_R,
                       double[:] dual_scaling,
                       # Parameters
                       double alpha,
                       double beta,
                       bint positive,
                       bint center,
                       # active set for improved performance
                       int[:] disabled,
                       unsigned int n_disabled=0
                       ) nogil:

    cdef double R_norm2
    cdef double R_sum = 0.
    cdef double w_norm2
    cdef double y_norm2
    cdef double l1_norm
    cdef double const
    cdef double gap
    cdef unsigned int ii
    cdef unsigned int jj
    cdef double yTA
    cdef double dual_norm_XtA

    # sparse X.T * R product for non-disabled features
    if center:
        # Rtot_sum = R_sum - n_samples * R_shift  is a constant when columns are centered
        R_sum = y_sum + n_samples * R_shift

    for ii in range(n_features):
        if (n_disabled == 0 or disabled[ii] == 0):
            X_T_R[ii] = 0.0
            for jj in range(X_indptr[ii], X_indptr[ii + 1]):
                X_T_R[ii] += X_data[jj] * R[X_indices[jj]]
            if center:
                X_T_R[ii] -= X_mean[ii] * R_sum
            XtA[ii] = X_T_R[ii] - beta * w[ii]
        else:
            XtA[ii] = 0.

    if positive:
        dual_norm_XtA = max(n_features, &XtA[0])
    else:
        dual_norm_XtA = abs_max(n_features, &XtA[0])

    yTA = ddot(n_samples, <DOUBLE*>&y[0], 1, <DOUBLE*>&R[0], 1)

    if center:
        yTA -= y_sum * R_shift
    # R_norm2 = np.dot(R, R)
    R_norm2 = ddot(n_samples, <DOUBLE*>&R[0], 1, <DOUBLE*>&R[0], 1)
    if center:
        R_norm2 += n_samples * R_shift**2 - 2 * R_sum * R_shift

    if alpha == 0:
        dual_scaling[0] = 1.
    else:
        if dual_norm_XtA <= 0:
            if R_norm2 == 0:
                dual_scaling[0] = 1. / alpha
            else:
                dual_scaling[0] = yTA / R_norm2 / alpha
        elif positive:
            dual_scaling[0] = fmin(yTA / (alpha * R_norm2),
                                 1. / dual_norm_XtA)
        else:
            dual_scaling[0] = fmin(fmax(yTA / (alpha * R_norm2),
                              -1. / dual_norm_XtA), 1. / dual_norm_XtA)

    dual_scaling[0] = 1. / dual_scaling[0]

    # w_norm2 = np.dot(w, w)
    if beta > 0:
        w_norm2 = ddot(n_features, <DOUBLE*>&w[0], 1, <DOUBLE*>&w[0], 1)
    else:
        w_norm2 = 0.

    const = alpha / dual_scaling[0]
    A_norm2 = R_norm2 * const**2
    l1_norm = dasum(n_features, <DOUBLE*>&w[0], 1)
    gap = 0.5 * (R_norm2 + A_norm2)
    gap += (alpha * l1_norm - const * yTA
            + 0.5 * beta * (1. + const ** 2) * (w_norm2))

    return gap


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def sparse_enet_coordinate_descent(double[:] w,
                            double alpha, double beta,
                            double[:] X_data, int[:] X_indices,
                            int[:] X_indptr, double[:] y,
                            double[:] X_mean, int max_iter,
                            double tol, object rng, bint random=0,
                            bint positive=0, int screening=0
                            ):
    """Cython version of the coordinate descent algorithm for Elastic-Net

    We minimize:

        (1/2) * norm(y - X w, 2)^2 + alpha norm(w, 1) + (beta/2) * norm(w, 2)^2

    """

    # get the data information into easy vars
    cdef unsigned int n_samples = y.shape[0]
    cdef unsigned int n_features = w.shape[0]

    # compute norms of the columns of X
    cdef unsigned int ii
    cdef double[:] norm2_cols_X = np.zeros(n_features, np.float64)
    cdef double[:] norm_cols_X = np.zeros(n_features, np.float64)

    cdef unsigned int startptr = X_indptr[0]
    cdef unsigned int endptr

    # get the number of tasks indirectly, using strides
    cdef unsigned int n_tasks = y.strides[0] / sizeof(DOUBLE)

    # initial value of the residuals
    # The total residuals are: Rtot = R - R_shift
    cdef double[:] R = y.copy()
    cdef double R_shift = 0

    cdef double[:] X_T_R = np.zeros(n_features)
    cdef double[:] Xty = np.zeros(n_features)
    cdef double[:] XtA = np.zeros(n_features)
    cdef double[:] Xty_div_alpha = np.zeros(n_features)

    cdef int[:] active_set = np.empty(n_features, dtype=np.intc)
    cdef int[:] disabled = np.zeros(n_features, dtype=np.intc)

    cdef double d_w_max
    cdef double w_max
    cdef double tmp
    cdef double w_ii
    cdef double d_w_ii
    cdef double X_mean_ii

    cdef double y_sum
    cdef double R_sum
    cdef double Rtot_sum

    cdef double normalize_sum
    cdef double gap = tol + 1.0
    cdef double d_w_tol = tol
    cdef double Xtymax
    cdef double y_norm2

    cdef double r_small2
    cdef double r_large
    cdef double tmp_XtA_dual_scaling

    cdef double[:] dual_scaling = np.zeros(1)
    cdef unsigned int jj
    cdef unsigned int n_iter = 0
    cdef unsigned int f_iter
    cdef unsigned int n_features_dual_gap
    cdef UINT32_t rand_r_state_seed = rng.randint(0, RAND_R_MAX)
    cdef UINT32_t* rand_r_state = &rand_r_state_seed
    cdef bint center = False
    cdef unsigned int n_active = 0
    cdef bint do_gap = False

    cdef double[:] Xistar = np.zeros(n_samples)
    cdef double[:] XtXistar = np.zeros(n_features)

    if alpha == 0:
        warnings.warn("Coordinate descent with alpha=0 may lead to unexpected"
                      " results and is discouraged.")

    with nogil:

        # center = (X_mean != 0).any()
        for ii in range(n_features):
            if X_mean[ii]:
                center = True
                break

        if center:
            y_sum = 0.0
            for jj in range(n_samples):
                y_sum += y[jj]
            Rtot_sum = y_sum

        # R = y - np.dot(X, w)
        for ii in range(n_features):
            X_mean_ii = X_mean[ii]
            endptr = X_indptr[ii + 1]
            normalize_sum = 0.0
            w_ii = w[ii]

            for jj in range(startptr, endptr):
                normalize_sum += (X_data[jj] - X_mean_ii) ** 2
                R[X_indices[jj]] -= X_data[jj] * w_ii

            norm2_cols_X[ii] = normalize_sum + \
                (n_samples - endptr + startptr) * X_mean_ii ** 2

            if center:
                #for jj in range(n_samples):
                #    R[jj] += X_mean_ii * w_ii
                R_shift -= X_mean_ii * w_ii
            startptr = endptr

        #norm of columns
        for ii in range(n_features):
            norm_cols_X[ii] = sqrt(norm2_cols_X[ii])

        if screening == 0:
            n_active = 0
            for ii in range(n_features):
                if norm2_cols_X[ii] != 0.0:
                    active_set[n_active] = ii
                    n_active += 1
                else:
                    disabled[ii] = 1

        ### scaling tolerance
        # tol *= np.dot(y, y)
        y_norm2 = ddot(n_samples, <DOUBLE*>&y[0], 1, <DOUBLE*>&y[0], 1)
        tol *= y_norm2

        for n_iter in range(max_iter):

            ####################
            # variable screening
            do_gap = False
            if (n_iter and w_max == 0.0  # heuristic termination criterion
                    or d_w_max / w_max < d_w_tol):
                # the biggest coordinate update of this iteration was smaller
                # than the tolerance: check the duality gap as ultimate
                # stopping criterion
                do_gap = True

            do_gap = (do_gap or (screening > 0 and n_iter % screening == 0)
                      or (n_iter == max_iter - 1)) or (n_iter == 0)

            if do_gap and screening > 0:  # Screening
                if screening > 0:
                    n_features_dual_gap = n_features - n_active
                else:
                    n_features_dual_gap = 0

                gap = sparse_enet_duality_gap(n_samples, n_features, X_data,
                                              X_indices, X_indptr, X_mean,
                                              y, y_sum, R, R_shift, w,
                                              XtA, X_T_R, dual_scaling,
                                              alpha, beta, positive, center,
                                              disabled, n_features_dual_gap)

                # return if we reached desired tolerance
                if gap < tol:
                    break

            if do_gap and screening > 0:  # Screening

                r_screening = sqrt(2. * gap) / alpha
                n_active = 0

                for ii in range(n_features):  # Loop over coordinates

                    if disabled[ii] == 1:
                         continue

                    tmp = XtA[ii] / dual_scaling[0]

                    if not positive:
                        tmp = fabs(tmp)

                    if tmp >= (1. - r_screening * norm_cols_X[ii]) or w[ii] != 0.:
                        active_set[n_active] = ii
                        n_active += 1
                    else:
                        disabled[ii] = 1


            ###################
            # Coordinate descent

            w_max = 0.0
            d_w_max = 0.0

            for f_iter in range(n_active):  # Loop over coordinates
                if random:
                    ii = active_set[rand_int(n_active, rand_r_state)]
                else:
                    ii = active_set[f_iter]

                startptr = X_indptr[ii]
                endptr = X_indptr[ii + 1]
                w_ii = w[ii]  # Store previous value
                X_mean_ii = X_mean[ii]

                # tmp = (X[:,ii] * R).sum() - X_mean_ii * R_sum + norm2_cols_X[ii] * w[ii]
                tmp = 0.0
                # tmp = (X[:,ii] * R).sum()
                for jj in range(startptr, endptr):
                    tmp += R[X_indices[jj]] * X_data[jj]

                if center:
                    # tmp -= X_mean_ii * R_sum
                    # Note that Rtot_sum is a constant when the columns are centered
                    tmp -= (Rtot_sum + n_samples*R_shift) * X_mean_ii

                tmp += norm2_cols_X[ii] * w[ii] # end tmp computation

                if positive and tmp < 0.0:
                    w[ii] = 0.0
                else:
                    w[ii] = fsign(tmp) * fmax(fabs(tmp) - alpha, 0) \
                            / (norm2_cols_X[ii] + beta)

                d_w_ii = w[ii] - w_ii

                if w[ii] != w_ii:
                    # R +=  d_w_ii * X[:,ii] # Update residual
                    for jj in range(startptr, endptr):
                        R[X_indices[jj]] -=  X_data[jj] * d_w_ii

                    if center:
                        R_shift -= X_mean_ii * d_w_ii

                # update the maximum absolute coefficient update
                d_w_ii = fabs(d_w_ii)
                if d_w_ii > d_w_max:
                    d_w_max = d_w_ii

                if fabs(w[ii]) > w_max:
                    w_max = fabs(w[ii])

    return w, gap, tol, n_iter + 1


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def enet_coordinate_descent_gram(double[:] w, double alpha, double beta,
                                 np.ndarray[double, ndim=2, mode='c'] Q,
                                 np.ndarray[double, ndim=1, mode='c'] q,
                                 np.ndarray[double, ndim=1] y,
                                 int max_iter, double tol, object rng,
                                 bint random=0, bint positive=0):
    """Cython version of the coordinate descent algorithm
        for Elastic-Net regression

        We minimize

        (1/2) * w^T Q w - q^T w + alpha norm(w, 1) + (beta/2) * norm(w, 2)^2

        which amount to the Elastic-Net problem when:
        Q = X^T X (Gram matrix)
        q = X^T y
    """

    # get the data information into easy vars
    cdef unsigned int n_samples = y.shape[0]
    cdef unsigned int n_features = Q.shape[0]

    # initial value "Q w" which will be kept of up to date in the iterations
    cdef double[:] H = np.dot(Q, w)

    cdef double[:] XtA = np.zeros(n_features)
    cdef double tmp
    cdef double w_ii
    cdef double d_w_max
    cdef double w_max
    cdef double d_w_ii
    cdef double gap = tol + 1.0
    cdef double d_w_tol = tol
    cdef double dual_norm_XtA
    cdef unsigned int ii
    cdef unsigned int n_iter = 0
    cdef unsigned int f_iter
    cdef UINT32_t rand_r_state_seed = rng.randint(0, RAND_R_MAX)
    cdef UINT32_t* rand_r_state = &rand_r_state_seed

    cdef double y_norm2 = np.dot(y, y)
    cdef double* w_ptr = <double*>&w[0]
    cdef double* Q_ptr = &Q[0, 0]
    cdef double* q_ptr = <double*>q.data
    cdef double* H_ptr = &H[0]
    cdef double* XtA_ptr = &XtA[0]
    tol = tol * y_norm2

    if alpha == 0:
        warnings.warn("Coordinate descent with alpha=0 may lead to unexpected"
            " results and is discouraged.")

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
                    daxpy(n_features, -w_ii, Q_ptr + ii * n_features, 1,
                          H_ptr, 1)

                tmp = q[ii] - H[ii]

                if positive and tmp < 0:
                    w[ii] = 0.0
                else:
                    w[ii] = fsign(tmp) * fmax(fabs(tmp) - alpha, 0) \
                        / (Q[ii, ii] + beta)

                if w[ii] != 0.0:
                    # H +=  w[ii] * Q[ii] # Update H = X.T X w
                    daxpy(n_features, w[ii], Q_ptr + ii * n_features, 1,
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
                q_dot_w = ddot(n_features, w_ptr, 1, q_ptr, 1)

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
                w_norm2 = ddot(n_features, &w[0], 1, &w[0], 1)

                if (dual_norm_XtA > alpha):
                    const = alpha / dual_norm_XtA
                    A_norm2 = R_norm2 * (const ** 2)
                    gap = 0.5 * (R_norm2 + A_norm2)
                else:
                    const = 1.0
                    gap = R_norm2

                # The call to dasum is equivalent to the L1 norm of w
                gap += (alpha * dasum(n_features, &w[0], 1) -
                        const * y_norm2 +  const * q_dot_w +
                        0.5 * beta * (1 + const ** 2) * w_norm2)

                if gap < tol:
                    # return if we reached desired tolerance
                    break

    return np.asarray(w), gap, tol, n_iter + 1




@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def enet_coordinate_descent_multi_task(double[::1, :] W, double l1_reg,
                                       double l2_reg,
                                       np.ndarray[double, ndim=2, mode='fortran'] X,
                                       np.ndarray[double, ndim=2] Y,
                                       int max_iter, double tol, object rng,
                                       bint random=0):
    """Cython version of the coordinate descent algorithm
        for Elastic-Net mult-task regression

        We minimize

        (1/2) * norm(y - X w, 2)^2 + l1_reg ||w||_21 + (1/2) * l2_reg norm(w, 2)^2

    """
    # get the data information into easy vars
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]
    cdef unsigned int n_tasks = Y.shape[1]

    # to store XtA
    cdef double[:, ::1] XtA = np.zeros((n_features, n_tasks))
    cdef double XtA_axis1norm
    cdef double dual_norm_XtA

    # initial value of the residuals
    cdef double[:, ::1] R = np.zeros((n_samples, n_tasks))

    cdef double[:] norm2_cols_X = np.zeros(n_features)
    cdef double[::1] tmp = np.zeros(n_tasks, dtype=np.float)
    cdef double[:] w_ii = np.zeros(n_tasks, dtype=np.float)
    cdef double d_w_max
    cdef double w_max
    cdef double d_w_ii
    cdef double nn
    cdef double W_ii_abs_max
    cdef double gap = tol + 1.0
    cdef double d_w_tol = tol
    cdef double ry_sum
    cdef double l21_norm
    cdef unsigned int ii
    cdef unsigned int jj
    cdef unsigned int n_iter = 0
    cdef unsigned int f_iter
    cdef UINT32_t rand_r_state_seed = rng.randint(0, RAND_R_MAX)
    cdef UINT32_t* rand_r_state = &rand_r_state_seed

    cdef double* X_ptr = &X[0, 0]
    cdef double* W_ptr = &W[0, 0]
    cdef double* Y_ptr = &Y[0, 0]
    cdef double* wii_ptr = &w_ii[0]

    if l1_reg == 0:
        warnings.warn("Coordinate descent with l1_reg=0 may lead to unexpected"
            " results and is discouraged.")

    with nogil:
        # norm2_cols_X = (np.asarray(X) ** 2).sum(axis=0)
        for ii in range(n_features):
            for jj in range(n_samples):
                norm2_cols_X[ii] += X[jj, ii] ** 2

        # R = Y - np.dot(X, W.T)
        for ii in range(n_samples):
            for jj in range(n_tasks):
                R[ii, jj] = Y[ii, jj] - (
                    ddot(n_features, X_ptr + ii, n_samples, W_ptr + jj, n_tasks)
                    )

        # tol = tol * linalg.norm(Y, ord='fro') ** 2
        tol = tol * dnrm2(n_samples * n_tasks, Y_ptr, 1) ** 2

        for n_iter in range(max_iter):
            w_max = 0.0
            d_w_max = 0.0
            for f_iter in range(n_features):  # Loop over coordinates
                if random:
                    ii = rand_int(n_features, rand_r_state)
                else:
                    ii = f_iter

                if norm2_cols_X[ii] == 0.0:
                    continue

                # w_ii = W[:, ii] # Store previous value
                dcopy(n_tasks, W_ptr + ii * n_tasks, 1, wii_ptr, 1)

                # if np.sum(w_ii ** 2) != 0.0:  # can do better
                if dnrm2(n_tasks, wii_ptr, 1) != 0.0:
                    # R += np.dot(X[:, ii][:, None], w_ii[None, :]) # rank 1 update
                    dger(CblasRowMajor, n_samples, n_tasks, 1.0,
                         X_ptr + ii * n_samples, 1,
                         wii_ptr, 1, &R[0, 0], n_tasks)

                # tmp = np.dot(X[:, ii][None, :], R).ravel()
                dgemv(CblasRowMajor, CblasTrans,
                      n_samples, n_tasks, 1.0, &R[0, 0], n_tasks,
                      X_ptr + ii * n_samples, 1, 0.0, &tmp[0], 1)

                # nn = sqrt(np.sum(tmp ** 2))
                nn = dnrm2(n_tasks, &tmp[0], 1)

                # W[:, ii] = tmp * fmax(1. - l1_reg / nn, 0) / (norm2_cols_X[ii] + l2_reg)
                dcopy(n_tasks, &tmp[0], 1, W_ptr + ii * n_tasks, 1)
                dscal(n_tasks, fmax(1. - l1_reg / nn, 0) / (norm2_cols_X[ii] + l2_reg),
                          W_ptr + ii * n_tasks, 1)

                # if np.sum(W[:, ii] ** 2) != 0.0:  # can do better
                if dnrm2(n_tasks, W_ptr + ii * n_tasks, 1) != 0.0:
                    # R -= np.dot(X[:, ii][:, None], W[:, ii][None, :])
                    # Update residual : rank 1 update
                    dger(CblasRowMajor, n_samples, n_tasks, -1.0,
                         X_ptr + ii * n_samples, 1, W_ptr + ii * n_tasks, 1,
                         &R[0, 0], n_tasks)

                # update the maximum absolute coefficient update
                d_w_ii = diff_abs_max(n_tasks, W_ptr + ii * n_tasks, wii_ptr)

                if d_w_ii > d_w_max:
                    d_w_max = d_w_ii

                W_ii_abs_max = abs_max(n_tasks, W_ptr + ii * n_tasks)
                if W_ii_abs_max > w_max:
                    w_max = W_ii_abs_max

            if w_max == 0.0 or d_w_max / w_max < d_w_tol or n_iter == max_iter - 1:
                # the biggest coordinate update of this iteration was smaller than
                # the tolerance: check the duality gap as ultimate stopping
                # criterion

                # XtA = np.dot(X.T, R) - l2_reg * W.T
                for ii in range(n_features):
                    for jj in range(n_tasks):
                        XtA[ii, jj] = ddot(
                            n_samples, X_ptr + ii * n_samples, 1,
                            &R[0, 0] + jj, n_tasks
                            ) - l2_reg * W[jj, ii]

                # dual_norm_XtA = np.max(np.sqrt(np.sum(XtA ** 2, axis=1)))
                dual_norm_XtA = 0.0
                for ii in range(n_features):
                    # np.sqrt(np.sum(XtA ** 2, axis=1))
                    XtA_axis1norm = dnrm2(n_tasks, &XtA[0, 0] + ii * n_tasks, 1)
                    if XtA_axis1norm > dual_norm_XtA:
                        dual_norm_XtA = XtA_axis1norm

                # TODO: use squared L2 norm directly
                # R_norm = linalg.norm(R, ord='fro')
                # w_norm = linalg.norm(W, ord='fro')
                R_norm = dnrm2(n_samples * n_tasks, &R[0, 0], 1)
                w_norm = dnrm2(n_features * n_tasks, W_ptr, 1)
                if (dual_norm_XtA > l1_reg):
                    const =  l1_reg / dual_norm_XtA
                    A_norm = R_norm * const
                    gap = 0.5 * (R_norm ** 2 + A_norm ** 2)
                else:
                    const = 1.0
                    gap = R_norm ** 2

                # ry_sum = np.sum(R * y)
                ry_sum = 0.0
                for ii in range(n_samples):
                    for jj in range(n_tasks):
                        ry_sum += R[ii, jj] * Y[ii, jj]

                # l21_norm = np.sqrt(np.sum(W ** 2, axis=0)).sum()
                l21_norm = 0.0
                for ii in range(n_features):
                    # np.sqrt(np.sum(W ** 2, axis=0))
                    l21_norm += dnrm2(n_tasks, W_ptr + n_tasks * ii, 1)

                gap += l1_reg * l21_norm - const * ry_sum + \
                     0.5 * l2_reg * (1 + const ** 2) * (w_norm ** 2)

                if gap < tol:
                    # return if we reached desired tolerance
                    break

    return np.asarray(W), gap, tol, n_iter + 1
