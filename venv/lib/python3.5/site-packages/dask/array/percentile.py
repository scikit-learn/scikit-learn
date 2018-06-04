from __future__ import absolute_import, division, print_function

from functools import wraps
from collections import Iterator
from numbers import Number

import numpy as np
from toolz import merge, merge_sorted

from .core import Array
from ..base import tokenize
from .. import sharedict


@wraps(np.percentile)
def _percentile(a, q, interpolation='linear'):
    n = len(a)
    if not len(a):
        return None, n
    if isinstance(q, Iterator):
        q = list(q)
    if a.dtype.name == 'category':
        result = np.percentile(a.codes, q, interpolation=interpolation)
        import pandas as pd
        return pd.Categorical.from_codes(result, a.categories, a.ordered), n
    if np.issubdtype(a.dtype, np.datetime64):
        a2 = a.astype('i8')
        result = np.percentile(a2, q, interpolation=interpolation)
        return result.astype(a.dtype), n
    if not np.issubdtype(a.dtype, np.number):
        interpolation = 'nearest'
    return np.percentile(a, q, interpolation=interpolation), n


def percentile(a, q, interpolation='linear'):
    """ Approximate percentile of 1-D array

    See :func:`numpy.percentile` for more information
    """
    if not a.ndim == 1:
        raise NotImplementedError(
            "Percentiles only implemented for 1-d arrays")
    if isinstance(q, Number):
        q = [q]
    q = np.array(q)
    token = tokenize(a, list(q), interpolation)
    name = 'percentile_chunk-' + token
    dsk = dict(((name, i), (_percentile, (key), q, interpolation))
               for i, key in enumerate(a.__dask_keys__()))

    name2 = 'percentile-' + token
    dsk2 = {(name2, 0): (merge_percentiles, q, [q] * len(a.chunks[0]),
                         sorted(dsk), interpolation)}

    dtype = a.dtype
    if np.issubdtype(dtype, np.integer):
        dtype = (np.array([], dtype=dtype) / 0.5).dtype

    dsk = merge(dsk, dsk2)
    dsk = sharedict.merge(a.dask, (name2, dsk))
    return Array(dsk, name2, chunks=((len(q),),), dtype=dtype)


def merge_percentiles(finalq, qs, vals, interpolation='lower', Ns=None):
    """ Combine several percentile calculations of different data.

    Parameters
    ----------

    finalq : numpy.array
        Percentiles to compute (must use same scale as ``qs``).
    qs : sequence of :class:`numpy.array`s
        Percentiles calculated on different sets of data.
    vals : sequence of :class:`numpy.array`s
        Resulting values associated with percentiles ``qs``.
    Ns : sequence of integers
        The number of data elements associated with each data set.
    interpolation : {'linear', 'lower', 'higher', 'midpoint', 'nearest'}
        Specify the type of interpolation to use to calculate final
        percentiles.  For more information, see :func:`numpy.percentile`.

    Examples
    --------

    >>> finalq = [10, 20, 30, 40, 50, 60, 70, 80]
    >>> qs = [[20, 40, 60, 80], [20, 40, 60, 80]]
    >>> vals = [np.array([1, 2, 3, 4]), np.array([10, 11, 12, 13])]
    >>> Ns = [100, 100]  # Both original arrays had 100 elements

    >>> merge_percentiles(finalq, qs, vals, Ns=Ns)
    array([ 1,  2,  3,  4, 10, 11, 12, 13])
    """
    if isinstance(finalq, Iterator):
        finalq = list(finalq)
    finalq = np.array(finalq)
    qs = list(map(list, qs))
    vals = list(vals)
    if Ns is None:
        vals, Ns = zip(*vals)
    Ns = list(Ns)

    L = list(zip(*[(q, val, N) for q, val, N in zip(qs, vals, Ns) if N]))
    if not L:
        raise ValueError("No non-trivial arrays found")
    qs, vals, Ns = L

    # TODO: Perform this check above in percentile once dtype checking is easy
    #       Here we silently change meaning
    if vals[0].dtype.name == 'category':
        result = merge_percentiles(finalq, qs, [v.codes for v in vals], interpolation, Ns)
        import pandas as pd
        return pd.Categorical.from_codes(result, vals[0].categories, vals[0].ordered)
    if not np.issubdtype(vals[0].dtype, np.number):
        interpolation = 'nearest'

    if len(vals) != len(qs) or len(Ns) != len(qs):
        raise ValueError('qs, vals, and Ns parameters must be the same length')

    # transform qs and Ns into number of observations between percentiles
    counts = []
    for q, N in zip(qs, Ns):
        count = np.empty(len(q))
        count[1:] = np.diff(q)
        count[0] = q[0]
        count *= N
        counts.append(count)

    # Sort by calculated percentile values, then number of observations.
    # >95% of the time in this function is spent in `merge_sorted` below.
    # An alternative that uses numpy sort is shown.  It is sometimes
    # comparable to, but typically slower than, `merge_sorted`.
    #
    # >>> A = np.concatenate(map(np.array, map(zip, vals, counts)))
    # >>> A.sort(0, kind='mergesort')

    combined_vals_counts = merge_sorted(*map(zip, vals, counts))
    combined_vals, combined_counts = zip(*combined_vals_counts)

    combined_vals = np.array(combined_vals)
    combined_counts = np.array(combined_counts)

    # percentile-like, but scaled by total number of observations
    combined_q = np.cumsum(combined_counts)

    # rescale finalq percentiles to match combined_q
    desired_q = finalq * sum(Ns)

    # the behavior of different interpolation methods should be
    # investigated further.
    if interpolation == 'linear':
        rv = np.interp(desired_q, combined_q, combined_vals)
    else:
        left = np.searchsorted(combined_q, desired_q, side='left')
        right = np.searchsorted(combined_q, desired_q, side='right') - 1
        np.minimum(left, len(combined_vals) - 1, left) # don't exceed max index
        lower = np.minimum(left, right)
        upper = np.maximum(left, right)
        if interpolation == 'lower':
            rv = combined_vals[lower]
        elif interpolation == 'higher':
            rv = combined_vals[upper]
        elif interpolation == 'midpoint':
            rv = 0.5 * (combined_vals[lower] + combined_vals[upper])
        elif interpolation == 'nearest':
            lower_residual = np.abs(combined_q[lower] - desired_q)
            upper_residual = np.abs(combined_q[upper] - desired_q)
            mask = lower_residual > upper_residual
            index = lower  # alias; we no longer need lower
            index[mask] = upper[mask]
            rv = combined_vals[index]
        else:
            raise ValueError("interpolation can only be 'linear', 'lower', "
                             "'higher', 'midpoint', or 'nearest'")
    return rv
