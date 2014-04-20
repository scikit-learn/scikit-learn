# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#         Olivier Grisel <olivier.grisel@ensta.org>
#         Alexis Mignon <alexis.mignon@gmail.com>
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


np.import_array()


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
    enum CBLAS_ORDER:
        CblasRowMajor=101
        CblasColMajor=102
    enum CBLAS_TRANSPOSE:
        CblasNoTrans=111
        CblasTrans=112
        CblasConjTrans=113
        AtlasConj=114

    void daxpy "cblas_daxpy"(int N, double alpha, double *X, int incX,
                             double *Y, int incY)
    double ddot "cblas_ddot"(int N, double *X, int incX, double *Y, int incY)
    void dger "cblas_dger"(CBLAS_ORDER Order, int M, int N, double alpha,
                double *X, int incX, double *Y, int incY, double *A, int lda)
    void dgemv "cblas_dgemv"(CBLAS_ORDER Order,
                      CBLAS_TRANSPOSE TransA, int M, int N,
                      double alpha, double *A, int lda,
                      double *X, int incX, double beta,
                      double *Y, int incY)
    double dnrm2 "cblas_dnrm2"(int N, double *X, int incX)
    void dcopy "cblas_dcopy"(int N, double *X, int incX, double *Y, int incY)
    void dscal "cblas_dscal"(int N, double alpha, double *X, int incX)


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

    tol = tol * np.dot(y, y)

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

            R_norm2 = np.dot(R, R)
            w_norm2 = np.dot(w, w)
            if (dual_norm_XtA > alpha):
                const = alpha / dual_norm_XtA
                A_norm2 = R_norm2 * (const**2)
                gap = 0.5 * (R_norm2 + A_norm2)
            else:
                const = 1.0
                gap = R_norm2

            gap += alpha * linalg.norm(w, 1) - const * np.dot(R.T, y) + \
                  0.5 * beta * (1 + const ** 2) * (w_norm2)

            if gap < tol:
                # return if we reached desired tolerance
                break

    return w, gap, tol


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def sparse_enet_coordinate_descent(double[:] w,
                            double alpha, double beta,
                            double[:] X_data, int[:] X_indices,
                            int[:] X_indptr, double[:] y,
                            double[:] X_mean, int max_iter,
                            double tol, bint positive=False):
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
    cdef double[:] norm_cols_X = np.zeros(n_features, np.float64)

    cdef unsigned int startptr = X_indptr[0]
    cdef unsigned int endptr

    # initial value of the residuals
    cdef double[:] R = y.copy()

    cdef double[:] X_T_R = np.zeros(n_features)
    cdef double[:] XtA = np.zeros(n_features)

    cdef double tmp
    cdef double w_ii
    cdef double d_w_max
    cdef double w_max
    cdef double d_w_ii
    cdef double X_mean_ii
    cdef double R_sum
    cdef double normalize_sum
    cdef double gap = tol + 1.0
    cdef double d_w_tol = tol
    cdef unsigned int jj
    cdef unsigned int n_iter
    cdef bint center = False

    # center = (X_mean != 0).any()
    for ii in xrange(n_samples):
        if X_mean[ii]:
           center = True
           break

    for ii in xrange(n_features):
        X_mean_ii = X_mean[ii]
        endptr = X_indptr[ii + 1]
        normalize_sum = 0.0
        w_ii = w[ii]

        for jj in xrange(startptr, endptr):
            normalize_sum += (X_data[jj] - X_mean_ii) ** 2
            R[X_indices[jj]] -= X_data[jj] * w_ii
        norm_cols_X[ii] = normalize_sum + \
            (n_samples - endptr + startptr) * X_mean_ii ** 2

        if center:
            for jj in xrange(n_samples):
                R[jj] += X_mean_ii * w_ii
        startptr = endptr

    #tol *= np.dot(y, y)
    tol *= ddot(n_samples, <DOUBLE*>&y[0], 1, <DOUBLE*>&y[0], 1)

    for n_iter in range(max_iter):

        w_max = 0.0
        d_w_max = 0.0
        startptr = X_indptr[0]

        for ii in xrange(n_features):  # Loop over coordinates

            if norm_cols_X[ii] == 0.0:
                continue

            endptr = X_indptr[ii + 1]
            w_ii = w[ii]  # Store previous value
            X_mean_ii = X_mean[ii]

            if w_ii != 0.0:
                # R += w_ii * X[:,ii]
                for jj in xrange(startptr, endptr):
                    R[X_indices[jj]] += X_data[jj] * w_ii
                if center:
                    for jj in xrange(n_samples):
                        R[jj] -= X_mean_ii * w_ii

            # tmp = (X[:,ii] * R).sum()
            tmp = 0.0
            for jj in xrange(startptr, endptr):
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
                for jj in xrange(startptr, endptr):
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
            startptr = endptr
        if w_max == 0.0 or d_w_max / w_max < d_w_tol or n_iter == max_iter - 1:
            # the biggest coordinate update of this iteration was smaller than
            # the tolerance: check the duality gap as ultimate stopping
            # criterion

            # sparse X.T / dense R dot product
            for ii in xrange(n_features):
                for jj in xrange(X_indptr[ii], X_indptr[ii + 1]):
                    X_T_R[ii] += X_data[jj] * R[X_indices[jj]]
                R_sum = 0.0
                for jj in xrange(n_samples):
                    R_sum += R[jj]
                X_T_R[ii] -= X_mean[ii] * R_sum

            for jj in xrange(n_features):
                XtA[jj] = X_T_R[jj] - beta * w[jj]
            if positive:
                dual_norm_XtA = np.max(XtA)
            else:
                dual_norm_XtA = linalg.norm(XtA, np.inf)

            #R_norm2 = np.dot(R, R)
            R_norm2 = ddot(n_samples, <DOUBLE*>&R[0], 1, <DOUBLE*>&R[0], 1)

            #w_norm2 = np.dot(w, w)
            w_norm2 = ddot(n_features, <DOUBLE*>&w[0], 1, <DOUBLE*>&w[0], 1)
            if (dual_norm_XtA > alpha):
                const = alpha / dual_norm_XtA
                A_norm2 = R_norm2 * const**2
                gap = 0.5 * (R_norm2 + A_norm2)
            else:
                const = 1.0
                gap = R_norm2

            gap += alpha * linalg.norm(w, 1) - const * np.dot(R.T, y) + \
                  0.5 * beta * (1 + const ** 2) * (w_norm2)

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

    cdef double y_norm2 = np.dot(y, y)
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
            w_norm2 = np.dot(w, w)
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
                  0.5 * beta * (1 + const ** 2) * (w_norm2)

            if gap < tol:
                # return if we reached desired tolerance
                break

    return w, gap, tol


cdef double abs_max(int n, double* a):
    """np.max(np.abs(a))"""
    cdef int i
    cdef double m = fabs(a[0])
    cdef double d
    for i in xrange(1, n):
        d = fabs(a[i])
        if d > m:
            m = d
    return m


cdef double diff_abs_max(int n, double* a, double* b):
    """np.max(np.abs(a - b))"""
    cdef int i
    cdef double m = fabs(a[0] - b[0])
    cdef double d
    for i in xrange(1, n):
        d = fabs(a[i] - b[i])
        if d > m:
            m = d
    return m


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def enet_coordinate_descent_multi_task(np.ndarray[DOUBLE, ndim=2, mode='fortran'] W,
                                       double l1_reg, double l2_reg,
                                       np.ndarray[DOUBLE, ndim=2, mode='fortran'] X,
                                       np.ndarray[DOUBLE, ndim=2] Y,
                                       int max_iter, double tol):
    """Cython version of the coordinate descent algorithm
        for Elastic-Net mult-task regression

        We minimize

        1 norm(y - X w, 2)^2 + l1_reg ||w||_21 + l2_reg norm(w, 2)^2
        -                                       ----
        2                                        2

    """
    # get the data information into easy vars
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]
    cdef unsigned int n_tasks = Y.shape[1]

    # compute norms of the columns of X
    cdef np.ndarray[DOUBLE, ndim=1] norm_cols_X = (X ** 2).sum(axis=0)

    # initial value of the residuals
    cdef np.ndarray[DOUBLE, ndim=2, mode='c'] R

    cdef np.ndarray[DOUBLE, ndim=1, mode='c'] tmp = np.zeros(n_tasks, dtype=np.float)
    cdef np.ndarray[DOUBLE, ndim=1] w_ii = np.zeros(n_tasks, dtype=np.float)
    cdef double d_w_max
    cdef double w_max
    cdef double d_w_ii
    cdef double nn
    cdef double W_ii_abs_max
    cdef double gap = tol + 1.0
    cdef double d_w_tol = tol
    cdef unsigned int ii
    cdef unsigned int n_iter

    if l1_reg == 0:
        warnings.warn("Coordinate descent with l1_reg=0 may lead to unexpected"
            " results and is discouraged.")

    R = Y - np.dot(X, W.T)
    R = np.asarray(R, order='C')

    # tol = tol * linalg.norm(Y, ord='fro') ** 2
    tol = tol * dnrm2(n_samples * n_tasks, <DOUBLE*>Y.data, 1) ** 2

    for n_iter in range(max_iter):
        w_max = 0.0
        d_w_max = 0.0
        for ii in xrange(n_features): # Loop over coordinates
            if norm_cols_X[ii] == 0.0:
                continue

            # w_ii = W[:, ii] # Store previous value
            dcopy(n_tasks, <DOUBLE*>(W.data + ii * n_tasks * sizeof(DOUBLE)),
                  1, <DOUBLE*>w_ii.data, 1)

            # if np.sum(w_ii ** 2) != 0.0:  # can do better
            if dnrm2(n_tasks, <DOUBLE*>w_ii.data, 1) != 0.0:
                # R += np.dot(X[:, ii][:, None], w_ii[None, :]) # rank 1 update
                dger(CblasRowMajor, n_samples, n_tasks, 1.0,
                    <DOUBLE*>(X.data + ii * n_samples * sizeof(DOUBLE)), 1,
                    <DOUBLE*>w_ii.data, 1,
                    <DOUBLE*>R.data, n_tasks)

            # tmp = np.dot(X[:, ii][None, :], R).ravel()
            dgemv(CblasRowMajor, CblasTrans,
                  n_samples, n_tasks, 1.0, <DOUBLE*>R.data,
                  n_tasks, <DOUBLE*>(X.data + ii * n_samples * sizeof(DOUBLE)),
                  1, 0.0, <DOUBLE*>tmp.data, 1)

            # nn = sqrt(np.sum(tmp ** 2))
            nn = dnrm2(n_tasks, <DOUBLE*>tmp.data, 1)

            # W[:, ii] = tmp * fmax(1. - l1_reg / nn, 0) / (norm_cols_X[ii] + l2_reg)
            dcopy(n_tasks, <DOUBLE*>tmp.data,
                  1, <DOUBLE*>(W.data + ii * n_tasks * sizeof(DOUBLE)), 1)
            dscal(n_tasks, fmax(1. - l1_reg / nn, 0) / (norm_cols_X[ii] + l2_reg),
                  <DOUBLE*>(W.data + ii * n_tasks * sizeof(DOUBLE)), 1)

            # if np.sum(W[:, ii] ** 2) != 0.0:  # can do better
            if dnrm2(n_tasks, <DOUBLE*>(W.data + ii * n_tasks * sizeof(DOUBLE)), 1) != 0.0:
                # R -= np.dot(X[:, ii][:, None], W[:, ii][None, :]) # Update residual : rank 1 update
                dger(CblasRowMajor, n_samples, n_tasks, -1.0,
                    <DOUBLE*>(X.data + ii * n_samples * sizeof(DOUBLE)), 1,
                    <DOUBLE*>(W.data + ii * n_tasks * sizeof(DOUBLE)), 1,
                    <DOUBLE*>R.data, n_tasks)

            # update the maximum absolute coefficient update
            d_w_ii = diff_abs_max(n_tasks,
                                  <DOUBLE*>(W.data + ii * n_tasks * sizeof(DOUBLE)),
                                  <DOUBLE*>w_ii.data)
            if d_w_ii > d_w_max:
                d_w_max = d_w_ii

            W_ii_abs_max = abs_max(n_tasks,
                                   <DOUBLE*>(W.data + ii * n_tasks * sizeof(DOUBLE)))
            if W_ii_abs_max > w_max:
                w_max = W_ii_abs_max

        if w_max == 0.0 or d_w_max / w_max < d_w_tol or n_iter == max_iter - 1:
            # the biggest coordinate update of this iteration was smaller than
            # the tolerance: check the duality gap as ultimate stopping
            # criterion

            XtA = np.dot(X.T, R) - l2_reg * W.T
            dual_norm_XtA = np.max(np.sqrt(np.sum(XtA ** 2, axis=1)))

            # TODO: use squared L2 norm directly
            # R_norm = linalg.norm(R, ord='fro')
            # w_norm = linalg.norm(W, ord='fro')
            R_norm = dnrm2(n_samples * n_tasks, <DOUBLE*>R.data, 1)
            w_norm = dnrm2(n_features * n_tasks, <DOUBLE*>W.data, 1)
            if (dual_norm_XtA > l1_reg):
                const =  l1_reg / dual_norm_XtA
                A_norm = R_norm * const
                gap = 0.5 * (R_norm ** 2 + A_norm ** 2)
            else:
                const = 1.0
                gap = R_norm ** 2

            gap += l1_reg * np.sqrt(np.sum(W ** 2, axis=0)).sum() - const * np.sum(R * Y) + \
                  0.5 * l2_reg * (1 + const ** 2) * (w_norm ** 2)

            if gap < tol:
                # return if we reached desired tolerance
                break

    return W, gap, tol
