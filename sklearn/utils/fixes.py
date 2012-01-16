"""Compatibility fixes for older version of python, numpy and scipy

If you add content to this file, please give the version of the package
at which the fixe is no longer needed.
"""
# Authors: Emmanuelle Gouillart <emmanuelle.gouillart@normalesup.org>
#          Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Fabian Pedregosa <fpedregosa@acm.org>
#          Lars Buitinck <L.J.Buitinck@uva.nl>
# License: BSD

import collections
import numpy as np
from operator import itemgetter


try:
    Counter = collections.Counter
except AttributeError:
# Partial replacement for Python 2.7 collections.Counter
    class Counter(collections.defaultdict):
        def __init__(self, iterable=(), **kwargs):
            super(Counter, self).__init__(int, **kwargs)
            self.update(iterable)

        def most_common(self):
            return sorted(self.iteritems(), key=itemgetter(1), reverse=True)

        def update(self, other):
            """Adds counts for elements in other"""
            if isinstance(other, self.__class__):
                for x, n in other.iteritems():
                    self[x] += n
            else:
                for x in other:
                    self[x] += 1


def _unique(ar, return_index=False, return_inverse=False):
    """A replacement for the np.unique that appeared in numpy 1.4.

    While np.unique existed long before, keyword return_inverse was
    only added in 1.4.
    """
    try:
        ar = ar.flatten()
    except AttributeError:
        if not return_inverse and not return_index:
            items = sorted(set(ar))
            return np.asarray(items)
        else:
            ar = np.asarray(ar).flatten()

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

np_version = np.__version__.split('.')
if int(np_version[0]) < 2 and int(np_version[1]) < 5:
    unique = _unique
else:
    unique = np.unique


def _copysign(x1, x2):
    """Slow replacement for np.copysign, which was introduced in numpy 1.4"""
    return np.abs(x1) * np.sign(x2)

if not hasattr(np, 'copysign'):
    copysign = _copysign
else:
    copysign = np.copysign


def _in1d(ar1, ar2, assume_unique=False):
    """Replacement for in1d that is provided for numpy >= 1.4"""
    if not assume_unique:
        ar1, rev_idx = unique(ar1, return_inverse=True)
        ar2 = np.unique(ar2)
    ar = np.concatenate((ar1, ar2))
    # We need this to be a stable sort, so always use 'mergesort'
    # here. The values from the first array should always come before
    # the values from the second array.
    order = ar.argsort(kind='mergesort')
    sar = ar[order]
    equal_adj = (sar[1:] == sar[:-1])
    flag = np.concatenate((equal_adj, [False]))
    indx = order.argsort(kind='mergesort')[:len(ar1)]

    if assume_unique:
        return flag[indx]
    else:
        return flag[indx][rev_idx]

if not hasattr(np, 'in1d'):
    in1d = _in1d
else:
    in1d = np.in1d


def qr_economic(A, **kwargs):
    """Compat function for the QR-decomposition in economic mode

    Scipy 0.9 changed the keyword econ=True to mode='economic'
    """
    import scipy.linalg
    # trick: triangular solve has introduced in 0.9
    if hasattr(scipy.linalg, 'solve_triangular'):
        return scipy.linalg.qr(A, mode='economic', **kwargs)
    else:
        return scipy.linalg.qr(A, econ=True, **kwargs)


def savemat(file_name, mdict, oned_as="column", **kwargs):
    """MATLAB-format output routine that is compatible with SciPy 0.7's.

    0.7.2 (or .1?) added the oned_as keyword arg with 'column' as the default
    value. It issues a warning if this is not provided, stating that "This will
    change to 'row' in future versions."
    """
    import scipy.io
    try:
        return scipy.io.savemat(file_name, mdict, oned_as=oned_as, **kwargs)
    except TypeError:
        return scipy.io.savemat(file_name, mdict, **kwargs)

try:
    from numpy import count_nonzero
except ImportError:
    def count_nonzero(X):
        return len(np.flatnonzero(X))
