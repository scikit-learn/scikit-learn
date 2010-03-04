# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

# $Id$

from math import sqrt
import numpy as np
import scipy.linalg as linalg

def enet_coordinate_descent(X, y, alpha, beta, w, maxit=10, callback=None):
    """coordinate descent for Elastic-Net model
    """
    E = []
    norm_cols_X = np.sum(X ** 2, axis=0) # Compute norms of the columns of X
    nsamples, nfeatures = X.shape

    # Init residual
    R = np.empty(nfeatures+nsamples)
    R[:nsamples] = y - np.dot(X,w)
    R[nsamples:] = - sqrt(beta) * w

    for iter in xrange(maxit):
        for ii in xrange(nfeatures): # Loop over coordinates
            w_ii = w[ii] # Store previous value
            R[:nsamples] += w_ii * X[:, ii]
            R[nsamples + ii] += w_ii * sqrt(beta)
            tmp = (X[:, ii] * R[:nsamples]).sum()
            tmp += sqrt(beta) * R[nsamples + ii]
            w[ii] = np.sign(tmp) * np.maximum(abs(tmp) - alpha, 0) \
                    / (norm_cols_X[ii] + beta)
            R[:nsamples] -= w[ii] * X[:, ii] # Update residual
            R[nsamples + ii] -= w[ii] * sqrt(beta)

        E.append(0.5 * linalg.norm(R) ** 2 + alpha * np.abs(w).sum() +
                 0.5 * beta * linalg.norm(w) ** 2)
        if (callback is not None and not callback(X, y, alpha, w, iter)):
            break

    return w, E

