# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr> 
# License: BSD Style.

# $Id$

import numpy as np
import scipy.linalg as linalg

def lasso_coordinate_descent(X, y, alpha, w, maxit=10):
    """coordinate descent for Lasso model
    """
    E = []
    norm_cols_X = np.sum(X**2, axis=0) # Compute norms of the columns of X
    R = y.copy() # Init residual
    nsamples, nfeatures = X.shape

    for _ in xrange(maxit):
        for ii in xrange(nfeatures): # Loop over coordinates
            w_ii = w[ii] # Store previous value
            R += w_ii * X[:,ii]
            tmp = (X[:,ii]*R).sum()
            w[ii] = np.sign(tmp) * np.maximum(abs(tmp) - alpha,0) / norm_cols_X[ii]
            R -= w[ii] * X[:,ii] # Update residual

        E.append(0.5*linalg.norm(R)**2 + alpha*np.abs(w).sum())

    return w, E
