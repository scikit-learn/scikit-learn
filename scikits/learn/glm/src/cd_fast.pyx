# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
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

ctypedef np.float64_t DOUBLE

@cython.boundscheck(False)
@cython.wraparound(False)
def lasso_coordinate_descent(np.ndarray[DOUBLE, ndim=1] w,
                             double alpha,
                             np.ndarray[DOUBLE, ndim=2] X,
                             np.ndarray[DOUBLE, ndim=1] y,
                             int maxit, int gap_step,
                             double tol):
    """Cython version of the coordinate descent algorithm
        for Lasso regression

    gap_step : int
        calculate the duality each `gap_step` interations.

    Notes
    -----
    tolerance (tol) is scaled by ||y||^2
    """

    # get the data information into easy vars

    cdef unsigned int nsamples = X.shape[0]
    cdef unsigned int nfeatures = X.shape[1]
    cdef unsigned int nclasses = w.shape[1]

    cdef np.ndarray[DOUBLE, ndim=1] norm_cols_X
    cdef np.ndarray[DOUBLE, ndim=1] R = y - np.dot(X, w) # Init residual

    cdef double tmp, w_ii, gap
    cdef unsigned int ii
    cdef unsigned int jj
    cdef unsigned int n_iter
    cdef double dual_norm_XtA

    norm_cols_X = (X**2).sum(axis=0) # Compute norms of the columns of X

    tol = tol * linalg.norm(y)**2

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


        # calculate the dual gap each gap_step iterations
        if n_iter % gap_step == 0: continue
        dual_norm_XtA = linalg.norm(np.dot(X.T, R), np.inf)
        R_norm = linalg.norm(R)
        if (dual_norm_XtA > alpha):
            const =  alpha / dual_norm_XtA
            A_norm = R_norm * (const)
            gap = 0.5 * (R_norm**2 + A_norm**2)
        else:
            const = 1.0
            gap = R_norm**2

        gap = gap + alpha * linalg.norm(w, 1) - const * np.dot(R.T, y)

        if gap < tol:
            # return if we reached desired tolerance
            break

    return w, gap, tol



def enet_coordinate_descent(np.ndarray[DOUBLE, ndim=1] w,
                            double alpha, double beta,
                            np.ndarray[DOUBLE, ndim=2] X,
                            np.ndarray[DOUBLE, ndim=1] y,
                            int maxit, int gap_step, double tol):
    """Cython version of the coordinate descent algorithm
        for Elastic-Net regression

        We minimize

        1 norm(y - X w, 2)^2 + alpha norm(w, 1) + beta norm(w, 2)^2
        -                                         ----  
        2                                           2  

    """

    # get the data information into easy vars
    cdef unsigned int nsamples = X.shape[0]
    cdef unsigned int nfeatures = X.shape[1]
    cdef unsigned int nclasses = w.shape[1]

    # compute norms of the columns of X
    cdef np.ndarray[DOUBLE, ndim=1] norm_cols_X = (X**2).sum(axis=0)

    # initial value of the residuals
    cdef np.ndarray[DOUBLE, ndim=1] R

    cdef double tmp
    cdef double w_ii
    cdef unsigned int ii
    cdef unsigned int jj
    cdef unsigned int n_iter

    R = y - np.dot(X, w)

    tol = tol * linalg.norm(y)**2

    for n_iter in range(maxit):
        for ii in xrange(nfeatures): # Loop over coordinates
            w_ii = w[ii] # Store previous value

            if w_ii != 0.0:
                # R += w_ii * X[:,ii]
                for jj in range(nsamples):
                    R[jj] += w_ii * X[jj,ii]


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


        # calculate the dual gap each gap_step iterations
        # if n_iter % gap_step == 0: continue

        dual_norm_XtA = linalg.norm(np.dot(X.T, R) - np.sqrt(beta) * w, np.inf)
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
