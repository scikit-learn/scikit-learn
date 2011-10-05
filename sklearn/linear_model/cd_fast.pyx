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
    if x > y: return x
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

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def enet_coordinate_descent(np.ndarray[DOUBLE, ndim=1] w,
                            double alpha, double beta,
                            np.ndarray[DOUBLE, ndim=2] X,
                            np.ndarray[DOUBLE, ndim=1] y,
                            int max_iter, double tol):
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

    R = y - np.dot(X, w)

    tol = tol * linalg.norm(y) ** 2

    for n_iter in range(max_iter):
        w_max = 0.0
        d_w_max = 0.0
        for ii in xrange(n_features): # Loop over coordinates
            if norm_cols_X[ii] == 0.0:
                continue

            w_ii = w[ii] # Store previous value

            if w_ii != 0.0:
                # R += w_ii * X[:,ii]
                daxpy(n_samples, w_ii,
                      <DOUBLE*>(X.data + ii * n_samples * sizeof(DOUBLE)), 1,
                      <DOUBLE*>R.data, 1)

            # tmp = (X[:,ii]*R).sum()
            tmp = ddot(n_samples,
                       <DOUBLE*>(X.data + ii * n_samples * sizeof(DOUBLE)), 1,
                       <DOUBLE*>R.data, 1)

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

            dual_norm_XtA = linalg.norm(np.dot(X.T, R) - beta * w, np.inf)
            # TODO: use squared L2 norm directly
            R_norm = linalg.norm(R)
            w_norm = linalg.norm(w, 2)
            if (dual_norm_XtA > alpha):
                const =  alpha / dual_norm_XtA
                A_norm = R_norm * const
                gap = 0.5 * (R_norm**2 + A_norm**2)
            else:
                const = 1.0
                gap = R_norm**2

            gap += alpha * linalg.norm(w, 1) - const * np.dot(R.T, y) + \
                  0.5 * beta * (1 + const**2) * (w_norm**2)

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
                            int max_iter, double tol):
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

    for n_iter in range(max_iter):
        w_max = 0.0
        d_w_max = 0.0
        for ii in xrange(n_features): # Loop over coordinates
            if Q[ii,ii] == 0.0:
                continue

            w_ii = w[ii] # Store previous value

            if w_ii != 0.0:
                # H -= w_ii * Q[ii]
                daxpy(n_features, -w_ii,
                      <DOUBLE*>(Q.data + ii * n_features * sizeof(DOUBLE)), 1,
                      <DOUBLE*>H.data, 1)

            tmp = q[ii] - H[ii]

            w[ii] = fsign(tmp) * fmax(fabs(tmp) - alpha, 0) \
                    / (Q[ii,ii] + beta)

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
            dual_norm_XtA = linalg.norm(q - H - beta * w, np.inf)
            R_norm2 = y_norm2 + np.sum(w * H) - 2.0 * q_dot_w
            w_norm = linalg.norm(w, 2)
            if (dual_norm_XtA > alpha):
                const =  alpha / dual_norm_XtA
                A_norm2 = R_norm2 * (const**2)
                gap = 0.5 * (R_norm2 + A_norm2)
            else:
                const = 1.0
                gap = R_norm2

            gap += alpha * linalg.norm(w, 1) \
                   - const * y_norm2 \
                   + const * q_dot_w + \
                  0.5 * beta * (1 + const**2) * (w_norm**2)

            if gap < tol:
                # return if we reached desired tolerance
                break

    return w, gap, tol
