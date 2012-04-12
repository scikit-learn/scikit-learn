"""
Non-metric Multdimensional Scaling
"""

import numpy as np


def PVA(distances, similarities, copy=False):
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
