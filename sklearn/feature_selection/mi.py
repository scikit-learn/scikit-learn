"""
Methods for calculating Mutual Information in an embarrassingly parallel way.
"""

# Author: Daniel Homola <dani.homola@gmail.com>
# License: BSD 3 clause


import numpy as np
from scipy.special import gamma, psi
from ..neighbors import NearestNeighbors
from joblib import Parallel, delayed
import multiprocessing

num_cores = multiprocessing.cpu_count()


def get_mi_vector(MI_FS, F, s):
    """
    Calculates the Mututal Information between each feature in F and s.

    This function is for when |S| > 1. s is the previously selected feature.
    We exploite the fact that this step is embarrassingly parallel.
    """

    MIs = Parallel(n_jobs=num_cores)(delayed(_get_mi)(f, s, MI_FS) for f in F)
    return MIs


def _get_mi(f, s, MI_FS):
    n, p = MI_FS.X.shape
    if MI_FS.method in ['JMI', 'JMIM']:
        # JMI & JMIM
        joint = MI_FS.X[:, (s, f)]
        if MI_FS.categorical:
            MI = _mi_dc(joint, MI_FS.y, MI_FS.k)
        else:
            vars = (joint, MI_FS.y)
            MI = _mi_cc(vars, MI_FS.k)
    else:
        # MRMR
        vars = (MI_FS.X[:, s].reshape(n, 1), MI_FS.X[:, f].reshape(n, 1))
        MI = _mi_cc(vars, MI_FS.k)

    # MI must be non-negative
    if MI > 0:
        return MI
    else:
        return np.nan


def get_first_mi_vector(MI_FS, k):
    """
    Calculates the Mututal Information between each feature in X and y.

    This function is for when |S| = 0. We select the first feautre in S.
    """
    n, p = MI_FS.X.shape
    MIs = Parallel(
        n_jobs=num_cores)(
        delayed(_get_first_mi)(
            i, k, MI_FS) for i in xrange(p))
    return MIs


def _get_first_mi(i, k, MI_FS):
    n, p = MI_FS.X.shape
    if MI_FS.categorical:
        x = MI_FS.X[:, i].reshape((n, 1))
        MI = _mi_dc(x, MI_FS.y, k)
    else:
        vars = (MI_FS.X[:, i].reshape((n, 1)), MI_FS.y)
        MI = _mi_cc(vars, k)

    # MI must be non-negative
    if MI > 0:
        return MI
    else:
        return np.nan


def _mi_dc(x, y, k):
    """
    Calculates the mututal information between a continuous vector x and a
    disrete class vector y.

    This implementation can calculate the MI between the joint distribution of
    one or more continuous variables (X[:, 1:3]) with a discrete variable (y).

    Thanks to Adam Pocock, the author of the FEAST package for the idea.

    Brian C. Ross, 2014, PLOS ONE
    Mutual Information between Discrete and Continuous Data Sets
    """

    y = y.flatten()
    n = x.shape[0]
    classes = np.unique(y)
    knn = NearestNeighbors(n_neighbors=k)
    # distance to kth in-class neighbour
    d2k = np.empty(n)
    # number of points within each point's class
    Nx = []
    for yi in y:
        Nx.append(np.sum(y == yi))

    # find the distance of the kth in-class point
    for c in classes:
        mask = np.where(y == c)[0]
        knn.fit(x[mask, :])
        d2k[mask] = knn.kneighbors()[0][:, -1]

    # find the number of points within the distance of the kth in-class point
    knn.fit(x)
    m = knn.radius_neighbors(radius=d2k, return_distance=False)
    m = [i.shape[0] for i in m]

    # calculate MI based on Equation 2 in Ross 2014
    MI = psi(n) - np.mean(psi(Nx)) + psi(k) - np.mean(psi(m))
    return MI


def _mi_cc(variables, k=1):
    """
    Returns the mutual information between any number of variables.

    Here it is used to estimate MI between continuous X(s) and y.
    Written by Gael Varoquaux:
    https://gist.github.com/GaelVaroquaux/ead9898bd3c973c40429
    """

    all_vars = np.hstack(variables)
    return (sum([_entropy(X, k=k) for X in variables]) -
            _entropy(all_vars, k=k))


def _nearest_distances(X, k=1):
    """
    Returns the distance to the kth nearest neighbor for every point in X
    """

    knn = NearestNeighbors(n_neighbors=k, metric='chebyshev')
    knn.fit(X)
    # the first nearest neighbor is itself
    d, _ = knn.kneighbors(X)
    # returns the distance to the kth nearest neighbor
    return d[:, -1]


def _entropy(X, k=1):
    """
    Returns the entropy of the X.

    Written by Gael Varoquaux:
    https://gist.github.com/GaelVaroquaux/ead9898bd3c973c40429

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)
        The data the entropy of which is computed
    k : int, optional
        number of nearest neighbors for density estimation

    References
    ----------
    Kozachenko, L. F. & Leonenko, N. N. 1987 Sample estimate of entropy
    of a random vector. Probl. Inf. Transm. 23, 95-101.
    See also: Evans, D. 2008 A computationally efficient estimator for
    mutual information, Proc. R. Soc. A 464 (2093), 1203-1215.
    and:

    Kraskov A, Stogbauer H, Grassberger P. (2004). Estimating mutual
    information. Phys Rev E 69(6 Pt 2):066138.

    F. Perez-Cruz, (2008). Estimation of Information Theoretic Measures
    for Continuous Random Variables. Advances in Neural Information
    Processing Systems 21 (NIPS). Vancouver (Canada), December.
    return d*mean(log(r))+log(volume_unit_ball)+log(n-1)-log(k)

    """

    # Distance to kth nearest neighbor
    r = _nearest_distances(X, k)
    n, d = X.shape
    volume_unit_ball = (np.pi ** (.5 * d)) / gamma(.5 * d + 1)
    return (d * np.mean(np.log(r + np.finfo(X.dtype).eps)) +
            np.log(volume_unit_ball) + psi(n) - psi(k))
