# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#         Olivier Grisel <olivier.grisel@ensta.org>
#         Alexis Mignon <alexis.mignon@gmail.com>
#
# License: BSD Style.

cimport numpy as np
import numpy as np
import numpy.linalg as linalg
cimport cython
from cpython cimport bool
import warnings

cdef extern from "math.h":
    double fabs(double f)
    double sqrt(double f)

cdef inline double fmax(double x, double y):
    if x > y:
        return x
    return y

cdef inline double fsign(double f):
    if f == 0:
        return 0
    elif f > 0:
        return 1.0
    else:
        return -1.0

cdef extern from "cblas.h":
    void daxpy "cblas_daxpy"(int N, double alpha, double *X, int incX,
                             double *Y, int incY)
    double ddot "cblas_ddot"(int N, double *X, int incX, double *Y, int incY)


ctypedef np.float64_t DOUBLE
ctypedef np.int32_t INTEGER


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def sparse_std(unsigned int n_samples,
               unsigned int n_features,
               np.ndarray[DOUBLE, ndim=1] X_data,
               np.ndarray[INTEGER, ndim=1] X_indices,
               np.ndarray[INTEGER, ndim=1] X_indptr,
               np.ndarray[DOUBLE, ndim=1] X_mean=None):
    cdef unsigned int ii
    cdef unsigned int jj
    cdef unsigned int nnz_ii
    cdef double X_sum_ii
    cdef double X_mean_ii
    cdef double diff
    cdef double X_std_ii

    cdef np.ndarray[DOUBLE, ndim = 1] X_std = np.zeros(n_features, np.float64)

    if X_mean is None:
        X_mean = np.zeros(n_features, np.float64)

        for ii in xrange(n_features):
            # Computes the mean
            X_sum_ii = 0.0
            for jj in xrange(X_indptr[ii], X_indptr[ii + 1]):
                X_sum_ii += X_data[jj]
            X_mean[ii] = X_sum_ii / n_samples

    for ii in xrange(n_features):
        X_mean_ii = X_mean[ii]
        X_sum_ii = 0.0
        nnz_ii = 0
        for jj in xrange(X_indptr[ii], X_indptr[ii + 1]):
            diff = X_data[jj] - X_mean_ii
            X_sum_ii += diff * diff
            nnz_ii += 1

        X_std[ii] = (X_sum_ii + (n_samples - nnz_ii) * X_mean_ii * X_mean_ii)
    return np.sqrt(X_std)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def enet_coordinate_descent(np.ndarray[DOUBLE, ndim=1] w,
                            double alpha, double beta,
                            np.ndarray[DOUBLE, ndim=2] X,
                            np.ndarray[DOUBLE, ndim=1] y,
                            int max_iter, double tol, bool positive=False):
    """Cython version of the coordinate descent algorithm
        for Elastic-Net regression

        We minimize

        1 norm(y - X w, 2)^2 + alpha norm(w, 1) + beta norm(w, 2)^2
        -                                         ----
        2                                           2

    """

    # get the data information into easy vars
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]

    # compute norms of the columns of X
    cdef np.ndarray[DOUBLE, ndim=1] norm_cols_X = (X**2).sum(axis=0)

    # initial value of the residuals
    cdef np.ndarray[DOUBLE, ndim=1] R

    cdef double tmp
    cdef double w_ii
    cdef double d_w_max
    cdef double w_max
    cdef double d_w_ii
    cdef double gap = tol + 1.0
    cdef double d_w_tol = tol
    cdef unsigned int ii
    cdef unsigned int n_iter

    if alpha == 0:
        warnings.warn("Coordinate descent with alpha=0 may lead to unexpected"
            " results and is discouraged.")

    R = y - np.dot(X, w)

    tol = tol * linalg.norm(y) ** 2

    for n_iter in range(max_iter):
        w_max = 0.0
        d_w_max = 0.0
        for ii in xrange(n_features):  # Loop over coordinates
            if norm_cols_X[ii] == 0.0:
                continue

            w_ii = w[ii]  # Store previous value

            if w_ii != 0.0:
                # R += w_ii * X[:,ii]
                daxpy(n_samples, w_ii,
                      <DOUBLE*>(X.data + ii * n_samples * sizeof(DOUBLE)), 1,
                      <DOUBLE*>R.data, 1)

            # tmp = (X[:,ii]*R).sum()
            tmp = ddot(n_samples,
                       <DOUBLE*>(X.data + ii * n_samples * sizeof(DOUBLE)), 1,
                       <DOUBLE*>R.data, 1)

            if positive and tmp < 0:
                w[ii] = 0.0
            else:
                w[ii] = fsign(tmp) * fmax(fabs(tmp) - alpha, 0) \
                    / (norm_cols_X[ii] + beta)

            if w[ii] != 0.0:
                # R -=  w[ii] * X[:,ii] # Update residual
                daxpy(n_samples, -w[ii],
                      <DOUBLE*>(X.data + ii * n_samples * sizeof(DOUBLE)), 1,
                      <DOUBLE*>R.data, 1)

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

            XtA = np.dot(X.T, R) - beta * w
            if positive:
                dual_norm_XtA = np.max(XtA)
            else:
                dual_norm_XtA = linalg.norm(XtA, np.inf)

            # TODO: use squared L2 norm directly
            R_norm = linalg.norm(R)
            w_norm = linalg.norm(w, 2)
            if (dual_norm_XtA > alpha):
                const = alpha / dual_norm_XtA
                A_norm = R_norm * const
                gap = 0.5 * (R_norm ** 2 + A_norm ** 2)
            else:
                const = 1.0
                gap = R_norm ** 2

            gap += alpha * linalg.norm(w, 1) - const * np.dot(R.T, y) + \
                  0.5 * beta * (1 + const ** 2) * (w_norm ** 2)

            if gap < tol:
                # return if we reached desired tolerance
                break

    return w, gap, tol


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def sparse_enet_coordinate_descent(np.ndarray[DOUBLE, ndim=1] w,
                            double alpha, double beta,
                            np.ndarray[DOUBLE, ndim=1] X_data,
                            np.ndarray[INTEGER, ndim=1] X_indices,
                            np.ndarray[INTEGER, ndim=1] X_indptr,
                            np.ndarray[DOUBLE, ndim=1] y,
                            np.ndarray[DOUBLE, ndim=1] X_mean,
                            int max_iter, double tol, bint positive=False):
    """Cython version of the coordinate descent algorithm for Elastic-Net

    We minimize:

        1 norm(y - X w, 2)^2 + alpha norm(w, 1) + beta norm(w, 2)^2
        -                                         ----
        2                                           2

    """

    # get the data information into easy vars
    cdef unsigned int n_samples = y.shape[0]
    cdef unsigned int n_features = w.shape[0]

    # compute norms of the columns of X
    cdef unsigned int ii
    cdef np.ndarray[DOUBLE, ndim = 1] norm_cols_X = np.zeros(n_features,
                                                           np.float64)
    for ii in xrange(n_features):
        norm_cols_X[ii] = ((X_data[X_indptr[ii]:X_indptr[ii + 1]] - \
        X_mean[ii]) ** 2).sum() + \
        (n_samples - X_indptr[ii + 1] + X_indptr[ii]) * X_mean[ii] ** 2

    # initial value of the residuals
    cdef np.ndarray[DOUBLE, ndim = 1] R

    cdef double tmp
    cdef double w_ii
    cdef double d_w_max
    cdef double w_max
    cdef double d_w_ii
    cdef double X_mean_ii
    cdef double R_sum
    cdef double gap = tol + 1.0
    cdef double d_w_tol = tol
    cdef unsigned int jj
    cdef unsigned int n_iter
    cdef bint center = (X_mean != 0).any()

    # initialize the residuals
    R = y.copy()

    for ii in xrange(n_features):
        # sparse X column / dense w dot product
        for jj in xrange(X_indptr[ii], X_indptr[ii + 1]):
            R[X_indices[jj]] -= X_data[jj] * w[ii]
        if center:
            R += X_mean[ii] * w[ii]

    tol = tol * linalg.norm(y) ** 2

    for n_iter in range(max_iter):

        w_max = 0.0
        d_w_max = 0.0

        for ii in xrange(n_features):  # Loop over coordinates

            if norm_cols_X[ii] == 0.0:
                continue

            w_ii = w[ii]  # Store previous value
            X_mean_ii = X_mean[ii]

            if w_ii != 0.0:
                # R += w_ii * X[:,ii]
                for jj in xrange(X_indptr[ii], X_indptr[ii + 1]):
                    R[X_indices[jj]] += X_data[jj] * w_ii
                if center:
                    for jj in xrange(n_samples):
                        R[jj] -= X_mean_ii * w[ii]

            # tmp = (X[:,ii] * R).sum()
            tmp = 0.0
            for jj in xrange(X_indptr[ii], X_indptr[ii + 1]):
                tmp += R[X_indices[jj]] * X_data[jj]

            if center:
                R_sum = 0.0
                for jj in xrange(n_samples):
                    R_sum += R[jj]
                tmp -= R_sum * X_mean_ii

            if positive and tmp < 0.0:
                w[ii] = 0.0
            else:
                w[ii] = fsign(tmp) * fmax(fabs(tmp) - alpha, 0) \
                        / (norm_cols_X[ii] + beta)

            if w[ii] != 0.0:
                # R -=  w[ii] * X[:,ii] # Update residual
                for jj in xrange(X_indptr[ii], X_indptr[ii + 1]):
                    R[X_indices[jj]] -= X_data[jj] * w[ii]

                if center:
                    for jj in xrange(n_samples):
                        R[jj] += X_mean_ii * w[ii]

            # update the maximum absolute coefficient update
            d_w_ii = fabs(w[ii] - w_ii)
            if d_w_ii > d_w_max:
                d_w_max = d_w_ii

            if w[ii] > w_max:
                w_max = w[ii]

        if w_max == 0.0 or d_w_max / w_max < d_w_tol or n_iter == max_iter - 1:
            # the biggest coordinate update of this iteration was smaller than
            # the tolerance: check the duality gap as ultimate stopping
            # criterion

            # sparse X.T / dense R dot product
            X_T_R = np.zeros(n_features)
            for ii in xrange(n_features):
                for jj in xrange(X_indptr[ii], X_indptr[ii + 1]):
                    X_T_R[ii] += X_data[jj] * R[X_indices[jj]]
                X_T_R[ii] -= X_mean[ii] * R.sum()

            XtA = X_T_R - beta * w
            if positive:
                dual_norm_XtA = np.max(XtA)
            else:
                dual_norm_XtA = linalg.norm(XtA, np.inf)

            # TODO: use squared L2 norm directly
            R_norm = linalg.norm(R)
            w_norm = linalg.norm(w, 2)
            if (dual_norm_XtA > alpha):
                const = alpha / dual_norm_XtA
                A_norm = R_norm * const
                gap = 0.5 * (R_norm ** 2 + A_norm ** 2)
            else:
                const = 1.0
                gap = R_norm ** 2

            gap += alpha * linalg.norm(w, 1) - const * np.dot(R.T, y) + \
                  0.5 * beta * (1 + const ** 2) * (w_norm ** 2)

            if gap < tol:
                # return if we reached desired tolerance
                break

    return w, gap, tol


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def enet_coordinate_descent_gram(np.ndarray[DOUBLE, ndim=1] w,
                            double alpha, double beta,
                            np.ndarray[DOUBLE, ndim=2] Q,
                            np.ndarray[DOUBLE, ndim=1] q,
                            np.ndarray[DOUBLE, ndim=1] y,
                            int max_iter, double tol, bool positive=False):
    """Cython version of the coordinate descent algorithm
        for Elastic-Net regression

        We minimize

        1 w^T Q w - q^T w + alpha norm(w, 1) + beta norm(w, 2)^2
        -                                      ----
        2                                        2

        which amount to the Elastic-Net problem when:
        Q = X^T X (Gram matrix)
        q = X^T y
    """

    # get the data information into easy vars
    cdef unsigned int n_samples = y.shape[0]
    cdef unsigned int n_features = Q.shape[0]

    # initial value "Q w" which will be kept of up to date in the iterations
    cdef np.ndarray[DOUBLE, ndim=1] H = np.dot(Q, w)

    cdef double tmp
    cdef double w_ii
    cdef double d_w_max
    cdef double w_max
    cdef double d_w_ii
    cdef double gap = tol + 1.0
    cdef double d_w_tol = tol
    cdef unsigned int ii
    cdef unsigned int n_iter

    cdef double y_norm2 = linalg.norm(y) ** 2
    tol = tol * y_norm2

    if alpha == 0:
        warnings.warn("Coordinate descent with alpha=0 may lead to unexpected"
            " results and is discouraged.")

    for n_iter in range(max_iter):
        w_max = 0.0
        d_w_max = 0.0
        for ii in xrange(n_features):  # Loop over coordinates
            if Q[ii, ii] == 0.0:
                continue

            w_ii = w[ii]  # Store previous value

            if w_ii != 0.0:
                # H -= w_ii * Q[ii]
                daxpy(n_features, -w_ii,
                      <DOUBLE*>(Q.data + ii * n_features * sizeof(DOUBLE)), 1,
                      <DOUBLE*>H.data, 1)

            tmp = q[ii] - H[ii]

            if positive and tmp < 0:
                w[ii] = 0.0
            else:
                w[ii] = fsign(tmp) * fmax(fabs(tmp) - alpha, 0) \
                    / (Q[ii, ii] + beta)

            if w[ii] != 0.0:
                # H +=  w[ii] * Q[ii] # Update H = X.T X w
                daxpy(n_features, w[ii],
                      <DOUBLE*>(Q.data + ii * n_features * sizeof(DOUBLE)), 1,
                      <DOUBLE*>H.data, 1)

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

            q_dot_w = np.dot(w, q)

            XtA = q - H - beta * w
            if positive:
                dual_norm_XtA = np.max(XtA)
            else:
                dual_norm_XtA = linalg.norm(XtA, np.inf)

            R_norm2 = y_norm2 + np.sum(w * H) - 2.0 * q_dot_w
            w_norm = linalg.norm(w, 2)
            if (dual_norm_XtA > alpha):
                const = alpha / dual_norm_XtA
                A_norm2 = R_norm2 * (const ** 2)
                gap = 0.5 * (R_norm2 + A_norm2)
            else:
                const = 1.0
                gap = R_norm2

            gap += alpha * linalg.norm(w, 1) \
                   - const * y_norm2 \
                   + const * q_dot_w + \
                  0.5 * beta * (1 + const ** 2) * (w_norm ** 2)

            if gap < tol:
                # return if we reached desired tolerance
                break

    return w, gap, tol
