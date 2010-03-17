# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

# $Id$

from math import sqrt
import numpy as np
import scipy.linalg as linalg

def enet_coordinate_descent(model, X, y, maxit):
    """coordinate descent for Elastic-Net model"""
    norm_cols_X = np.sum(X ** 2, axis=0) # Compute norms of the columns of X
    n_samples, n_features = X.shape

    alpha = model.alpha
    beta = model.beta
    callbacks = model.callbacks
    w = model.w

    # Init residual
    R = np.empty(n_features + n_samples)
    R[:n_samples] = y - np.dot(X,w)
    R[n_samples:] = - sqrt(beta) * w

    for callback in callbacks:
        callback(0) # Init callback

    goon = True
    for n_iter in range(maxit):
        for ii in xrange(n_features): # Loop over coordinates
            w_ii = w[ii] # Store previous value
            if w_ii != 0.0:
                R[:n_samples] += w_ii * X[:, ii]
                R[n_samples + ii] += w_ii * sqrt(beta)
            tmp = (X[:, ii] * R[:n_samples]).sum()
            tmp += sqrt(beta) * R[n_samples + ii]
            w[ii] = np.sign(tmp) * np.maximum(abs(tmp) - alpha, 0) \
                    / (norm_cols_X[ii] + beta)

            if w[ii] != 0.0:
                R[:n_samples] -= w[ii] * X[:, ii] # Update residual
                R[n_samples + ii] -= w[ii] * sqrt(beta)

        for callback in callbacks:
            if not callback(n_iter, X=X, y=y, w=w, alpha=alpha, beta=beta, R=R):
                goon *= False

        if not goon:
            break

    return w

