# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#         Olivier Grisel <olivier.grisel@ensta.org>
#
# License: BSD Style.

cimport numpy as np
import numpy as np
import numpy.linalg as linalg
cimport cython

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
def sparse_normalize(np.ndarray[DOUBLE, ndim=1] X_data,
                    np.ndarray[INTEGER, ndim=1] X_indices,
                    np.ndarray[INTEGER, ndim=1] X_indptr,
                    np.ndarray[DOUBLE, ndim=1] X_std):

    cdef unsigned int n_features = X_std.shape[0]
    cdef unsigned int ii
    cdef unsigned int jj

    for ii in xrange(n_features):
        # Computes the mean
        for jj in xrange(X_indptr[ii], X_indptr[ii + 1]):
            X_data[jj] /= X_std[ii]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def enet_coordinate_descent(np.ndarray[DOUBLE, ndim=1] w,
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
    cdef bint center = (X_mean!=0).any()

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
