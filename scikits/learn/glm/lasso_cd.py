# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

# $Id$

import numpy as np
import scipy.linalg as linalg

def lasso_coordinate_descent(model, X, y, maxit):
    """Coordinate descent for Lasso model"""
    norm_cols_X = np.sum(X**2, axis=0) # Compute norms of the columns of X
    n_samples, n_features = X.shape

    alpha = model.alpha
    callbacks = model.callbacks
    w = model.w

    R = y - np.dot(X,w) # Init residual

    for callback in callbacks:
        callback(0) # Init callback

    goon = True
    for n_iter in range(maxit):
        for ii in xrange(n_features): # Loop over coordinates
            w_ii = w[ii] # Store previous value
            if w_ii != 0.0:
                R += w_ii * X[:, ii]
            tmp = (X[:, ii] * R).sum()
            w[ii] = np.sign(tmp) * np.maximum(abs(tmp) - alpha, 0) \
                    / norm_cols_X[ii]
            if w[ii] != 0.0:
                R -= w[ii] * X[:, ii] # Update residual

        for callback in callbacks:
            if not callback(n_iter, X=X, y=y, w=w, alpha=alpha, R=R):
                goon *= False

        if not goon:
            break

    return w
