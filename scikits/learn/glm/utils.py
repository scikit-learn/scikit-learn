# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

# $Id$

import numpy as np
import scipy.linalg as linalg

def lasso_objective(X, y, w, alpha, **kwargs):
    """Compute objective for Lasso model

    Returns
    -------
    obj : the objective value

    """
    if kwargs.has_key('R'):
        R = kwargs['R']
    else:
        R = y - np.dot(X, w)

    cost = 0.5 * linalg.norm(R) ** 2 + alpha * np.abs(w).sum()
    return cost

def enet_objective(X, y, w, alpha, beta, **kwargs):
    """Compute objective for Elastic-Net model

    Returns
    -------
    obj : the objective value

    """
    cost = lasso_objective(X, y, w, alpha, **kwargs)
    cost += 0.5 * beta * linalg.norm(w) ** 2
    return cost

def density(w, **kwargs):
    """Compute density of a sparse vector
        Return a value between 0 and 1
    """
    d = 0 if w is None else float((w != 0).sum()) / w.size
    return d
