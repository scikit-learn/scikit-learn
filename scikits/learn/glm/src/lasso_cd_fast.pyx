# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

# $Id$

cimport numpy as np
import numpy as np
import scipy.linalg as linalg
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
                             unsigned int maxit, double tol):
    """Cython version of the coordinate descent algorithm
        for Lasso regression
    """

    # get the data information into easy vars

    cdef unsigned int nsamples = X.shape[0]
    cdef unsigned int nfeatures = X.shape[1]
    cdef unsigned int nclasses = w.shape[1]

    cdef np.ndarray[DOUBLE, ndim=1] norm_cols_X = (X**2).sum(axis=0) # Compute norms of the columns of X
    cdef np.ndarray[DOUBLE, ndim=1] R = y - np.dot(X, w) # Init residual

    cdef double tmp
    cdef double w_ii
    cdef unsigned int ii
    cdef unsigned int jj
    cdef unsigned int n_iter

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

        
        A = R
        XtA = np.dot(X.T, A)
        dual_norm_XtA = np.max(np.abs(XtA))
        A_norm = linalg.norm(A)
        R_norm = A_norm
        if (dual_norm_XtA > alpha):
            A_norm *= np.abs(alpha / dual_norm_XtA)
        pobj = 0.5 * R_norm**2 + alpha * np.abs(w).sum();
        dobj = - 0.5 * A_norm**2 + np.dot(A.T, y)
        gap = pobj - dobj

        if gap < tol:
            # TODO: 
            break

    return w
