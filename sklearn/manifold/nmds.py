"""
Non-metric Multdimensional Scaling
"""

import numpy as np
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
    indxs = similarities.argsort()[::-1]
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
    # FIXME this isn't working
    sim_arg = similarities.argsort()
    dis_arg = distances.argsort()
    disparities = distances.copy()
    for i, j in zip(sim_arg, dis_arg):
        disparities[i] = distances[j]
    return disparities


def smacof(similarities, metric=True, p=3, init=None,
           max_iter=300, verbose=False, eps=1e-3):
    """
    """
    n = similarities.shape[0]

    if similarities.shape[0] != similarities.shape[1]:
        raise ValueError("similarities must be a square array (shape=%r)" % \
                            similarities.shape)

    if np.any(similarities != similarities.T):
        raise ValueError("similarities must be symmetric")

    sim_flat = ((1 - np.tri(n)) * similarities).flatten()
    sim_flat_w = sim_flat[sim_flat != 0]
    if init is None:
        # Randomly choose initial configuration
        X = np.random.random(size=(n, p))
    else:
        if n != init.shape[0] or p != init.shape[1]:
            raise ValueError("init matrix should be of shape (%d, %d)" % \
                                 (n, p))
        X = init

    old_stress = None
    for it in range(max_iter):
        # Compute distance and monotonic regression
        dis = euclidean_distances(X)
        dis_flat = dis.flatten()
        # similarities with 0 are considered as missing values
        dis_flat_w = dis_flat[sim_flat != 0]

        if metric:
            disparities = similarities
        else:
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
                    (disparities + disparities.T).flatten()) ** 2).sum() / 2

        # Update X using the Guttman transform
        ratio = (disparities + disparities.T) / dis
        ratio[np.isinf(ratio) | np.isnan(ratio)] = 0
        B = - ratio + np.diag(ratio.sum(axis=1))
        X = 1. / n * np.dot(B, X)
        if verbose:
            print 'it: %d, stress %s' % (it, stress)
        if old_stress is not None:
            if(old_stress - stress) < eps:
                if verbose:
                    print 'breaking at iteration %d with stress %s' % (it,
                                                                       stress)
        old_stress = stress

    return X
