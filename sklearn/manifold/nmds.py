"""
Non-metric Multdimensional Scaling
"""

import numpy as np

from ..base import BaseEstimator
from ..metrics import euclidean_distances


# XXX FIXME: should we use proximities or similarities ??
def PAV(distances, similarities, copy=False, verbose=False):
    """
    Pool adjancent violators

    Parameters
    ----------
        distances:

        similarities:

        copy: boolean, optional
    """
    # FIXME ties ?
    indxs = similarities.argsort()
    dis = distances.copy()

    sort = True
    while sort:
        block = []
        # FIXME: deal with floating points errors. Else, it stays stuck in the
        # loop
        el = dis[indxs][:-1] <= dis[indxs][1:] + np.finfo(np.float).resolution
        print "number of wrongly sorted elements: %d" % el.sum()
        sort = False
        for i, element in enumerate(el):
            if not element:
                sort = True
                block.append(i)
            if element and block:
                # Average on the block
                block.append(i)
                dis[indxs[block]] = dis[indxs][block].mean()
                block = []
    return dis


def guttman_image_ranking(distances, similarities, verbose=False):
    """
    Guttman image ranking

    """
    sim_arg = similarities.argsort()
    dis_arg = distances.argsort()
    disparities = distances.copy()
    for i, j in zip(sim_arg, dis_arg):
        disparities[i] = distances[j]
    return disparities


def smacof(similarities, metric=True, p=2, init=None,
           max_iter=300, verbose=0, eps=1e-3):
    """
    Computes multidimensional scaling using SMACOF algorithm

    Parameters
    ----------
    similarities: symmetric ndarray, shape [n * n]
        similarities between the points

    metric: boolean, optional, default: True
        compute metric or nonmetric SMACOF algorithm

    p: int, optional, default: 2
        number of dimension in which to immerse the similarities
        overridden if initial array is provided.

    init: {None or ndarray}
        if None, randomly chooses the initial configuration
        if ndarray, initialize the SMACOF algorithm with this array

    max_iter: int, optional, default: 300
        Maximum number of iterations of the SMACOF algorithm for a single run

    verbose: int, optional, default: 0
        level of verbosity

    eps: float, optional, default: 1e-6
        relative tolerance w.r.t stress to declare converge

    Returns
    -------
    X: ndarray (n, p)
        coordinates of the n points in a p-space
    """
    n = similarities.shape[0]

    if similarities.shape[0] != similarities.shape[1]:
        raise ValueError("similarities must be a square array (shape=%d)" % \
                            n)

    if np.any(similarities != similarities.T):
        raise ValueError("similarities must be symmetric")

    sim_flat = ((1 - np.tri(n)) * similarities).flatten()
    sim_flat_w = sim_flat[sim_flat != 0]
    if init is None:
        # Randomly choose initial configuration
        X = np.random.random(size=(n, p))
    else:
        # overrides the parameter p
        p = init.shape[1]
        if n != init.shape[0]:
            raise ValueError("init matrix should be of shape (%d, %d)" % \
                                 (n, p))
        X = init

    old_stress = None
    for it in range(max_iter):
        # Compute distance and monotonic regression
        dis = euclidean_distances(X)

        if metric:
            disparities = similarities
        else:
            dis_flat = dis.flatten()
            # similarities with 0 are considered as missing values
            dis_flat_w = dis_flat[sim_flat != 0]

            # Compute the disparities using a monotonic regression
            disparities_flat = guttman_image_ranking(dis_flat_w,
                                                     sim_flat_w)
            disparities = dis_flat.copy()
            disparities[sim_flat != 0] = disparities_flat
            disparities = disparities.reshape((n, n))
            disparities *= np.sqrt((n * (n - 1) / 2) / \
                           (disparities ** 2).sum())

        # Compute stress
        stress = ((dis.flatten() - \
                    disparities.flatten()) ** 2).sum() / 2

        # Update X using the Guttman transform
        ratio = disparities / dis
        ratio[np.isinf(ratio) | np.isnan(ratio)] = 0
        B = - ratio + np.diag(ratio.sum(axis=1))
        X = 1. / n * np.dot(B, X)
        if verbose == 2:
            print 'it: %d, stress %s' % (it, stress)
        if old_stress is not None:
            if(old_stress - stress) < eps:
                if verbose:
                    print 'breaking at iteration %d with stress %s' % (it,
                                                                       stress)
                break
        old_stress = stress

    return X


class MDS(BaseEstimator):
    """
    Multidimensional scaling

    Parameters
    ----------

    Notes
    -----
    """
    def __init__(self, p=2, init=None, max_iter=300, eps=1e-3):
        self.p = p
        self.init = init
        self.max_iter = max_iter
        self.eps = eps

    def fit():
        """
        """
