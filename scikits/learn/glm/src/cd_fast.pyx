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
    double exp(double f)
    double rand()

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
def enet_coordinate_descent(np.ndarray[DOUBLE, ndim=1] w,
                            double alpha, double beta,
                            np.ndarray[DOUBLE, ndim=2] X,
                            np.ndarray[DOUBLE, ndim=1] y,
                            int maxit, double tol):
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
    cdef double gab = tol + 1.0
    cdef unsigned int ii
    cdef unsigned int n_iter

    R = y - np.dot(X, w)
    tol = tol * linalg.norm(y) ** 2

    for n_iter in range(maxit):
        w_max = 0.0
        d_w_max = 0.0
        for ii in xrange(n_features): # Loop over coordinates
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

            if w[ii] > w_max:
                w_max = w[ii]

        if w_max == 0.0 or d_w_max / w_max < tol or n_iter == maxit - 1:
            # the biggest coordinate update of this iteration was smaller than
            # the tolerance: check the duality gap as ultimate stopping
            # criterion

            dual_norm_XtA = linalg.norm(np.dot(X.T, R) - beta * w, np.inf)
            # TODO: use C sqrt
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

