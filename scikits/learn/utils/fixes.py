"""
Fixes for older version of numpy and scipy.
"""
# Authors: Emmanuelle Gouillart <emmanuelle.gouillart@normalesup.org>
#          Gael Varoquaux <gael.varoquaux@normalesup.org>
# License: BSD

import numpy as np

def _unique(ar, return_index=False, return_inverse=False):
    """ A replacement for np.unique that appeared in numpy 1.4.
    """
    try:
        ar = ar.flatten()
    except AttributeError:
        if not return_inverse and not return_index:
            items = sorted(set(ar))
            return np.asarray(items)
        else:
            ar = np.asanyarray(ar).flatten()

    if ar.size == 0:
        if return_inverse and return_index:
            return ar, np.empty(0, np.bool), np.empty(0, np.bool)
        elif return_inverse or return_index:
            return ar, np.empty(0, np.bool)
        else:
            return ar

    if return_inverse or return_index:
        perm = ar.argsort()
        aux = ar[perm]
        flag = np.concatenate(([True], aux[1:] != aux[:-1]))
        if return_inverse:
            iflag = np.cumsum(flag) - 1
            iperm = perm.argsort()
            if return_index:
                return aux[flag], perm[flag], iflag[iperm]
            else:
                return aux[flag], iflag[iperm]
        else:
            return aux[flag], perm[flag]

    else:
        ar.sort()
        flag = np.concatenate(([True], ar[1:] != ar[:-1]))
        return ar[flag]


def _copysign (x1, x2):
    """
    (slow) Replacement for np.copysign, which was introduced in numpy 1.4
    """
    return np.abs(x1) * np.sign(x2)

def _in1d(ar1, ar2, assume_unique=False):
    """ Replacement for in1d that is provided for numpy >= 1.4
    """
    if not assume_unique:
        ar1, rev_idx = unique(ar1, return_inverse=True)
        ar2 = np.unique(ar2)
    ar = np.concatenate( (ar1, ar2) )
    # We need this to be a stable sort, so always use 'mergesort'
    # here. The values from the first array should always come before
    # the values from the second array.
    order = ar.argsort(kind='mergesort')
    sar = ar[order]
    equal_adj = (sar[1:] == sar[:-1])
    flag = np.concatenate( (equal_adj, [False] ) )
    indx = order.argsort(kind='mergesort')[:len( ar1 )]

    if assume_unique:
        return flag[indx]
    else:
        return flag[indx][rev_idx]


def qr_economic(A, **kwargs):
    """
    Scipy 0.9 changed the keyword econ=True to mode='economic'
    """
    import scipy.linalg
    # trick: triangular solve has introduced in 0.9
    if hasattr(scipy.linalg, 'triangular_solve'):
        return scipy.linalg.qr(A, mode='economic', **kwargs)
    else:
        return scipy.linalg.qr(A, econ=True, **kwargs)

def arpack_eigsh(A, **kwargs):
    """
    Scipy 0.9 renamed eigen_symmetric to eigsh in
    scipy.sparse.linalg.eigen.arpack
    """
    from scipy.sparse.linalg.eigen import arpack
    if hasattr(arpack, 'eigsh'):
        return arpack.eigsh(A, **kwargs)
    else:
        return arpack.eigen_symmetric(A, **kwargs)

# export fixes for np < 1.4
if not hasattr(np, 'copysign'):
    in1d = _in1d
    copysign = _copysign
    unique = _unique
else:
    from numpy import in1d, copysign, unique




