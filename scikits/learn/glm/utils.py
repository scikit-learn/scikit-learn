# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

# $Id$

import numpy as np
import scipy.linalg as linalg

def enet_dual_gap(X, y, w, alpha, beta=0):
    """Compute dual gap for Elastic-Net model to check KKT optimality conditions

    Returns
    -------
    gap : the difference  primal_objective - dual_objective (should be positive)
        A value less that 1e-5 means convergence in practice
    primal_objective : the value of the objective function of the primal problem
    dual_objective : the value of the objective function of the dual problem

    """
    Xw = np.dot(X, w)
    A = (y - Xw)
    if beta > 0:
        B = - np.sqrt(beta) * w
    XtA = np.dot(X.T, A)
    if beta > 0:
        XtA += np.sqrt(beta) * B
    dual_norm_XtA = np.max(XtA)
    if (dual_norm_XtA > alpha):
        A *= alpha / dual_norm_XtA
        if beta > 0:
            B *= alpha / dual_norm_XtA
    pobj = 0.5 * linalg.norm(y - Xw)**2 + alpha * np.abs(w).sum() \
           + 0.5 * beta * linalg.norm(w)**2
    dobj = - 0.5 * linalg.norm(A)**2 + np.dot(A.T, y)
    if beta > 0:
        dobj += - 0.5 * linalg.norm(B)**2
    gap = pobj - dobj
    return gap, pobj, dobj

def lasso_dual_gap(X, y, w, alpha):
    """Compute dual gap for Lasso model to check KKT optimality conditions

    Returns
    -------
    gap : the difference  primal_objective - dual_objective (should be positive)
        A value less that 1e-5 means convergence in practice
    primal_objective : the value of the objective function of the primal problem
    dual_objective : the value of the objective function of the dual problem

    """
    return enet_dual_gap(X, y, w, alpha, beta=0)

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
