# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

# $Id$

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

ctypedef np.float64_t DOUBLE

# @cython.boundscheck(False)
# @cython.wraparound(False)
def lasso_coordinate_descent(np.ndarray[DOUBLE, ndim=1] w,
                             double alpha,
                             np.ndarray[DOUBLE, ndim=2] X,
                             np.ndarray[DOUBLE, ndim=1] y,
                             int maxit, double tol):
    """Cython version of the coordinate descent algorithm
        for Lasso regression
    """

    # get the data information into easy vars

    cdef unsigned int nsamples = X.shape[0]
    cdef unsigned int nfeatures = X.shape[1]
    cdef unsigned int nclasses = w.shape[1]

    cdef np.ndarray[DOUBLE, ndim=1] norm_cols_X = (X**2).sum(axis=0) # Compute norms of the columns of X
    cdef np.ndarray[DOUBLE, ndim=1] R = y - np.dot(X, w) # Init residual

    cdef double tmp, w_ii, gap
    cdef unsigned int ii
    cdef unsigned int jj
    cdef unsigned int n_iter
    cdef double dual_norm_XtA

    goon = True
    for n_iter in range(maxit):
        for ii in xrange(nfeatures): # Loop over coordinates
            w_ii = w[ii] # Store previous value

            if w_ii != 0.0:
                # R += w_ii * X[:,ii]
                for jj in range(nsamples):
                    R[jj] += w_ii * X[jj, ii]

            # tmp = (X[:,ii]*R).sum()
            tmp = 0.0
            for jj in range(nsamples):
                tmp += R[jj] * X[jj, ii]

            w[ii] = fsign(tmp) * fmax(fabs(tmp) - alpha, 0) / norm_cols_X[ii]

            if w[ii] != 0.0:
                # R -=  w[ii] * X[:,ii] # Update residual
                for jj in range(nsamples):
                    R[jj] -=  w[ii] * X[jj, ii] # Update residual


        dual_norm_XtA = linalg.norm(np.dot(X.T, R), ord=np.inf)
        gap = 0.0
        const = 1.0
        if (dual_norm_XtA > alpha):
            R_norm = linalg.norm(R)
            const =  np.abs(alpha / dual_norm_XtA)
            A_norm = R_norm * const
            gap = 0.5 * (R_norm**2 - A_norm**2)

        gap = gap + alpha * linalg.norm(w, ord=1) + const * np.dot(R.T, y)

        if gap < tol:
            # TODO: something about an exp that alex told me
            break

    return w, gap
