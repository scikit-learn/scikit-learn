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
    double fmax(double f1, double f2)
    double rand()

cdef inline double fsign(double f):
    if f == 0:
        return 0
    elif f > 0:
        return 1.0
    else:
        return -1.0

ctypedef np.float64_t DOUBLE

@cython.boundscheck(False)
@cython.wraparound(False)
def enet_coordinate_descent(model,
                            np.ndarray[DOUBLE, ndim=2] X,
                            np.ndarray[DOUBLE, ndim=1] y,
                            unsigned int maxit):
    """Cython version of the coordinate descent algorithm
        for Elastic-Net regression
    """

    # get the data information into easy vars
    cdef float alpha = model.alpha
    cdef float beta = model.beta
    cdef np.ndarray[DOUBLE, ndim=1] w = model.coef_
    callbacks = model.callbacks

    cdef unsigned int nsamples = X.shape[0]
    cdef unsigned int nfeatures = X.shape[1]
    cdef unsigned int nclasses = w.shape[1]

    # compute norms of the columns of X
    cdef np.ndarray[DOUBLE, ndim=1] norm_cols_X = (X**2).sum(axis=0)

    # initial value of the residuals
    cdef np.ndarray[DOUBLE, ndim=1] R = np.empty(nfeatures + nsamples)
    R[:nsamples] = y - np.dot(X, w)
    R[nsamples:] = - sqrt(beta) * w

    cdef float tmp
    cdef float w_ii
    cdef unsigned int ii
    cdef unsigned int jj
    cdef unsigned int n_iter

    for callback in callbacks:
        callback(0) # Init callback

    goon = True
    for n_iter in range(maxit):
        for ii in xrange(nfeatures): # Loop over coordinates
            w_ii = w[ii] # Store previous value

            if w_ii != 0.0:
                # R += w_ii * X[:,ii]
                for jj in range(nsamples):
                    R[jj] += w_ii * X[jj,ii]
                R[nsamples+ii] += w_ii * sqrt(beta)

            # tmp = (X[:,ii]*R).sum()
            tmp = 0.0
            for jj in range(nsamples):
                tmp += R[jj] * X[jj,ii]

            w[ii] = fsign(tmp) * fmax(fabs(tmp) - alpha, 0) \
                    / (norm_cols_X[ii] + beta)

            if w[ii] != 0.0:
                # R -=  w[ii] * X[:,ii] # Update residual
                for jj in range(nsamples):
                    R[jj] -=  w[ii] * X[jj,ii] # Update residual
                R[nsamples+ii] -= w[ii] * sqrt(beta)

        for callback in callbacks:
            if not callback(n_iter, X=X, y=y, w=w, alpha=alpha, beta=beta, R=R):
                goon *= False

        if not goon:
            break

    return w
