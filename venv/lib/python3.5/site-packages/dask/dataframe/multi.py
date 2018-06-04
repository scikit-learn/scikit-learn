"""
Algorithms that Involve Multiple DataFrames
===========================================

The pandas operations ``concat``, ``join``, and ``merge`` combine multiple
DataFrames.  This module contains analogous algorithms in the parallel case.

There are two important cases:

1.  We combine along a partitioned index
2.  We combine along an unpartitioned index or other column

In the first case we know which partitions of each dataframe interact with
which others.  This lets uss be significantly more clever and efficient.

In the second case each partition from one dataset interacts with all
partitions from the other.  We handle this through a shuffle operation.

Partitioned Joins
-----------------

In the first case where we join along a partitioned index we proceed in the
following stages.

1.  Align the partitions of all inputs to be the same.  This involves a call
    to ``dd.repartition`` which will split up and concat existing partitions as
    necessary.  After this step all inputs have partitions that align with
    each other.  This step is relatively cheap.
    See the function ``align_partitions``.
2.  Remove unnecessary partitions based on the type of join we perform (left,
    right, inner, outer).  We can do this at the partition level before any
    computation happens.  We'll do it again on each partition when we call the
    in-memory function.  See the function ``require``.
3.  Embarrassingly parallel calls to ``pd.concat``, ``pd.join``, or
    ``pd.merge``.  Now that the data is aligned and unnecessary blocks have
    been removed we can rely on the fast in-memory Pandas join machinery to
    execute joins per-partition.  We know that all intersecting records exist
    within the same partition


Hash Joins via Shuffle
----------------------

When we join along an unpartitioned index or along an arbitrary column any
partition from one input might interact with any partition in another.  In
this case we perform a hash-join by shuffling data in each input by that
column.  This results in new inputs with the same partition structure cleanly
separated along that column.

We proceed with hash joins in the following stages:

1.  Shuffle each input on the specified column.  See the function
    ``dask.dataframe.shuffle.shuffle``.
2.  Perform embarrassingly parallel join across shuffled inputs.
"""
from __future__ import absolute_import, division, print_function

from functools import wraps, partial
from warnings import warn

from toolz import merge_sorted, unique, first
import toolz
import pandas as pd

from ..base import tokenize
from ..compatibility import apply
from .core import (_Frame, DataFrame, Series, map_partitions, Index,
                   _maybe_from_pandas, new_dd_object, is_broadcastable)
from .io import from_pandas
from . import methods
from .shuffle import shuffle, rearrange_by_divisions
from .utils import strip_unknown_categories


def align_partitions(*dfs):
    """ Mutually partition and align DataFrame blocks

    This serves as precursor to multi-dataframe operations like join, concat,
    or merge.

    Parameters
    ----------
    dfs: sequence of dd.DataFrame, dd.Series and dd.base.Scalar
        Sequence of dataframes to be aligned on their index

    Returns
    -------
    dfs: sequence of dd.DataFrame, dd.Series and dd.base.Scalar
        These must have consistent divisions with each other
    divisions: tuple
        Full divisions sequence of the entire result
    result: list
        A list of lists of keys that show which data exist on which
        divisions
    """
    _is_broadcastable = partial(is_broadcastable, dfs)
    dfs1 = [df for df in dfs
            if isinstance(df, _Frame) and
            not _is_broadcastable(df)]
    if len(dfs) == 0:
        raise ValueError("dfs contains no DataFrame and Series")
    if not all(df.known_divisions for df in dfs1):
        raise ValueError("Not all divisions are known, can't align "
                         "partitions. Please use `set_index` "
                         "to set the index.")

    divisions = list(unique(merge_sorted(*[df.divisions for df in dfs1])))
    if len(divisions) == 1:  # single value for index
        divisions = (divisions[0], divisions[0])
    dfs2 = [df.repartition(divisions, force=True)
            if isinstance(df, _Frame) else df for df in dfs]

    result = list()
    inds = [0 for df in dfs]
    for d in divisions[:-1]:
        L = list()
        for i, df in enumerate(dfs2):
            if isinstance(df, _Frame):
                j = inds[i]
                divs = df.divisions
                if j < len(divs) - 1 and divs[j] == d:
                    L.append((df._name, inds[i]))
                    inds[i] += 1
                else:
                    L.append(None)
            else:    # Scalar has no divisions
                L.append(None)
        result.append(L)
    return dfs2, tuple(divisions), result


def _maybe_align_partitions(args):
    """Align DataFrame blocks if divisions are different.

    Note that if all divisions are unknown, but have equal npartitions, then
    they will be passed through unchanged. This is different than
    `align_partitions`, which will fail if divisions aren't all known"""
    _is_broadcastable = partial(is_broadcastable, args)
    dfs = [df for df in args
           if isinstance(df, _Frame) and
           not _is_broadcastable(df)]
    if not dfs:
        return args

    divisions = dfs[0].divisions
    if not all(df.divisions == divisions for df in dfs):
        dfs2 = iter(align_partitions(*dfs)[0])
        return [a if not isinstance(a, _Frame) else next(dfs2) for a in args]
    return args


def require(divisions, parts, required=None):
    """ Clear out divisions where required components are not present

    In left, right, or inner joins we exclude portions of the dataset if one
    side or the other is not present.  We can achieve this at the partition
    level as well

    >>> divisions = [1, 3, 5, 7, 9]
    >>> parts = [(('a', 0), None),
    ...          (('a', 1), ('b', 0)),
    ...          (('a', 2), ('b', 1)),
    ...          (None, ('b', 2))]

    >>> divisions2, parts2 = require(divisions, parts, required=[0])
    >>> divisions2
    (1, 3, 5, 7)
    >>> parts2  # doctest: +NORMALIZE_WHITESPACE
    ((('a', 0), None),
     (('a', 1), ('b', 0)),
     (('a', 2), ('b', 1)))

    >>> divisions2, parts2 = require(divisions, parts, required=[1])
    >>> divisions2
    (3, 5, 7, 9)
    >>> parts2  # doctest: +NORMALIZE_WHITESPACE
    ((('a', 1), ('b', 0)),
     (('a', 2), ('b', 1)),
     (None, ('b', 2)))

    >>> divisions2, parts2 = require(divisions, parts, required=[0, 1])
    >>> divisions2
    (3, 5, 7)
    >>> parts2  # doctest: +NORMALIZE_WHITESPACE
    ((('a', 1), ('b', 0)),
     (('a', 2), ('b', 1)))
    """
    if not required:
        return divisions, parts
    for i in required:
        present = [j for j, p in enumerate(parts) if p[i] is not None]
        divisions = tuple(divisions[min(present): max(present) + 2])
        parts = tuple(parts[min(present): max(present) + 1])
    return divisions, parts


###############################################################
# Join / Merge
###############################################################


required = {'left': [0], 'right': [1], 'inner': [0, 1], 'outer': []}


def merge_indexed_dataframes(lhs, rhs, how='left', lsuffix='', rsuffix='',
                             indicator=False, left_on=None, right_on=None,
                             left_index=True, right_index=True):
    """ Join two partitioned dataframes along their index """

    (lhs, rhs), divisions, parts = align_partitions(lhs, rhs)
    divisions, parts = require(divisions, parts, required[how])

    left_empty = lhs._meta
    right_empty = rhs._meta

    name = 'join-indexed-' + tokenize(lhs, rhs, how, lsuffix, rsuffix,
                                      indicator, left_on, right_on,
                                      left_index, right_index)

    dsk = dict()
    for i, (a, b) in enumerate(parts):
        if a is None and how in ('right', 'outer'):
            a = left_empty
        if b is None and how in ('left', 'outer'):
            b = right_empty

        dsk[(name, i)] = (methods.merge, a, b, how, left_on, right_on,
                          left_index, right_index,
                          indicator, (lsuffix, rsuffix), left_empty,
                          right_empty)

    meta = pd.merge(lhs._meta_nonempty, rhs._meta_nonempty, how=how,
                    left_index=left_index, right_index=right_index,
                    left_on=left_on, right_on=right_on,
                    suffixes=(lsuffix, rsuffix), indicator=indicator)
    return new_dd_object(toolz.merge(lhs.dask, rhs.dask, dsk),
                         name, meta, divisions)


shuffle_func = shuffle  # name sometimes conflicts with keyword argument


def hash_join(lhs, left_on, rhs, right_on, how='inner',
              npartitions=None, suffixes=('_x', '_y'), shuffle=None,
              indicator=False):
    """ Join two DataFrames on particular columns with hash join

    This shuffles both datasets on the joined column and then performs an
    embarrassingly parallel join partition-by-partition

    >>> hash_join(a, 'id', rhs, 'id', how='left', npartitions=10)  # doctest: +SKIP
    """
    if npartitions is None:
        npartitions = max(lhs.npartitions, rhs.npartitions)

    lhs2 = shuffle_func(lhs, left_on, npartitions=npartitions, shuffle=shuffle)
    rhs2 = shuffle_func(rhs, right_on, npartitions=npartitions, shuffle=shuffle)

    if isinstance(left_on, Index):
        left_on = None
        left_index = True
    else:
        left_index = False

    if isinstance(right_on, Index):
        right_on = None
        right_index = True
    else:
        right_index = False

    # dummy result
    meta = pd.merge(lhs._meta_nonempty, rhs._meta_nonempty, how=how,
                    left_on=left_on, right_on=right_on,
                    left_index=left_index, right_index=right_index,
                    suffixes=suffixes, indicator=indicator)

    if isinstance(left_on, list):
        left_on = (list, tuple(left_on))
    if isinstance(right_on, list):
        right_on = (list, tuple(right_on))

    token = tokenize(lhs2, left_on, rhs2, right_on, left_index, right_index,
                     how, npartitions, suffixes, shuffle, indicator)
    name = 'hash-join-' + token

    dsk = {(name, i): (methods.merge,
                       (lhs2._name, i), (rhs2._name, i),
                       how, left_on, right_on,
                       left_index, right_index, indicator,
                       suffixes, lhs._meta, rhs._meta)
           for i in range(npartitions)}

    divisions = [None] * (npartitions + 1)
    return new_dd_object(toolz.merge(lhs2.dask, rhs2.dask, dsk),
                         name, meta, divisions)


def single_partition_join(left, right, **kwargs):
    # if the merge is perfomed on_index, divisions can be kept, otherwise the
    # new index will not necessarily correspond the current divisions

    meta = pd.merge(left._meta_nonempty, right._meta_nonempty, **kwargs)
    name = 'merge-' + tokenize(left, right, **kwargs)
    if left.npartitions == 1:
        left_key = first(left.__dask_keys__())
        dsk = {(name, i): (apply, pd.merge, [left_key, right_key], kwargs)
               for i, right_key in enumerate(right.__dask_keys__())}

        if kwargs.get('right_index') or right._contains_index_name(
                kwargs.get('right_on')):
            divisions = right.divisions
        else:
            divisions = [None for _ in right.divisions]

    elif right.npartitions == 1:
        right_key = first(right.__dask_keys__())
        dsk = {(name, i): (apply, pd.merge, [left_key, right_key], kwargs)
               for i, left_key in enumerate(left.__dask_keys__())}

        if kwargs.get('left_index') or left._contains_index_name(
                kwargs.get('left_on')):
            divisions = left.divisions
        else:
            divisions = [None for _ in left.divisions]

    return new_dd_object(toolz.merge(dsk, left.dask, right.dask), name,
                         meta, divisions)


@wraps(pd.merge)
def merge(left, right, how='inner', on=None, left_on=None, right_on=None,
          left_index=False, right_index=False, suffixes=('_x', '_y'),
          indicator=False, npartitions=None, shuffle=None, max_branch=None):
    for o in [on, left_on, right_on]:
        if isinstance(o, _Frame):
            raise NotImplementedError(
                "Dask collections not currently allowed in merge columns")
    if not on and not left_on and not right_on and not left_index and not right_index:
        on = [c for c in left.columns if c in right.columns]
        if not on:
            left_index = right_index = True

    if on and not left_on and not right_on:
        left_on = right_on = on
        on = None

    if (isinstance(left, (pd.Series, pd.DataFrame)) and
            isinstance(right, (pd.Series, pd.DataFrame))):
        return pd.merge(left, right, how=how, on=on, left_on=left_on,
                        right_on=right_on, left_index=left_index,
                        right_index=right_index, suffixes=suffixes,
                        indicator=indicator)

    # Transform pandas objects into dask.dataframe objects
    if isinstance(left, (pd.Series, pd.DataFrame)):
        if right_index and left_on:  # change to join on index
            left = left.set_index(left[left_on])
            left_on = False
            left_index = True
        left = from_pandas(left, npartitions=1)  # turn into DataFrame

    if isinstance(right, (pd.Series, pd.DataFrame)):
        if left_index and right_on:  # change to join on index
            right = right.set_index(right[right_on])
            right_on = False
            right_index = True
        right = from_pandas(right, npartitions=1)  # turn into DataFrame

    # Both sides are now dd.DataFrame or dd.Series objects
    merge_indexed_left = (left_index or left._contains_index_name(
        left_on)) and left.known_divisions

    merge_indexed_right = (right_index or right._contains_index_name(
        right_on)) and right.known_divisions

    # Both sides indexed
    if merge_indexed_left and merge_indexed_right:  # Do indexed join
        return merge_indexed_dataframes(left, right, how=how,
                                        lsuffix=suffixes[0],
                                        rsuffix=suffixes[1],
                                        indicator=indicator,
                                        left_on=left_on,
                                        right_on=right_on,
                                        left_index=left_index,
                                        right_index=right_index)

    # Single partition on one side
    elif (left.npartitions == 1 and how in ('inner', 'right') or
          right.npartitions == 1 and how in ('inner', 'left')):
        return single_partition_join(left, right, how=how, right_on=right_on,
                                     left_on=left_on, left_index=left_index,
                                     right_index=right_index,
                                     suffixes=suffixes, indicator=indicator)

    # One side is indexed, the other not
    elif (left_index and left.known_divisions and not right_index or
          right_index and right.known_divisions and not left_index):
        left_empty = left._meta_nonempty
        right_empty = right._meta_nonempty
        meta = pd.merge(left_empty, right_empty, how=how, on=on,
                        left_on=left_on, right_on=right_on,
                        left_index=left_index, right_index=right_index,
                        suffixes=suffixes, indicator=indicator)
        if merge_indexed_left and left.known_divisions:
            right = rearrange_by_divisions(right, right_on, left.divisions,
                                           max_branch, shuffle=shuffle)
            left = left.clear_divisions()
        elif merge_indexed_right and right.known_divisions:
            left = rearrange_by_divisions(left, left_on, right.divisions,
                                          max_branch, shuffle=shuffle)
            right = right.clear_divisions()
        return map_partitions(pd.merge, left, right, meta=meta, how=how, on=on,
                              left_on=left_on, right_on=right_on,
                              left_index=left_index, right_index=right_index,
                              suffixes=suffixes, indicator=indicator)
    # Catch all hash join
    else:
        return hash_join(left, left.index if left_index else left_on,
                         right, right.index if right_index else right_on,
                         how, npartitions, suffixes, shuffle=shuffle,
                         indicator=indicator)


###############################################################
# Concat
###############################################################

def concat_and_check(dfs):
    if len(set(map(len, dfs))) != 1:
        raise ValueError("Concatenated DataFrames of different lengths")
    return pd.concat(dfs, axis=1)


def concat_unindexed_dataframes(dfs):
    name = 'concat-' + tokenize(*dfs)

    dsk = {(name, i): (concat_and_check, [(df._name, i) for df in dfs])
           for i in range(dfs[0].npartitions)}

    meta = pd.concat([df._meta for df in dfs], axis=1)

    return new_dd_object(toolz.merge(dsk, *[df.dask for df in dfs]),
                         name, meta, dfs[0].divisions)


def concat_indexed_dataframes(dfs, axis=0, join='outer'):
    """ Concatenate indexed dataframes together along the index """
    meta = methods.concat([df._meta for df in dfs], axis=axis, join=join)
    empties = [strip_unknown_categories(df._meta) for df in dfs]

    dfs2, divisions, parts = align_partitions(*dfs)

    name = 'concat-indexed-' + tokenize(join, *dfs)

    parts2 = [[df if df is not None else empty
               for df, empty in zip(part, empties)]
              for part in parts]

    dsk = dict(((name, i), (methods.concat, part, axis, join))
               for i, part in enumerate(parts2))
    for df in dfs2:
        dsk.update(df.dask)

    return new_dd_object(dsk, name, meta, divisions)


def stack_partitions(dfs, divisions, join='outer'):
    """Concatenate partitions on axis=0 by doing a simple stack"""
    meta = methods.concat([df._meta for df in dfs], join=join)
    empty = strip_unknown_categories(meta)

    name = 'concat-{0}'.format(tokenize(*dfs))
    dsk = {}
    i = 0
    for df in dfs:
        dsk.update(df.dask)
        # An error will be raised if the schemas or categories don't match. In
        # this case we need to pass along the meta object to transform each
        # partition, so they're all equivalent.
        try:
            df._meta == meta
            match = True
        except (ValueError, TypeError):
            match = False

        for key in df.__dask_keys__():
            if match:
                dsk[(name, i)] = key
            else:
                dsk[(name, i)] = (methods.concat, [empty, key], 0, join)
            i += 1

    return new_dd_object(dsk, name, meta, divisions)


def concat(dfs, axis=0, join='outer', interleave_partitions=False):
    """ Concatenate DataFrames along rows.

    - When axis=0 (default), concatenate DataFrames row-wise:

      - If all divisions are known and ordered, concatenate DataFrames keeping
        divisions. When divisions are not ordered, specifying
        interleave_partition=True allows concatenate divisions each by each.

      - If any of division is unknown, concatenate DataFrames resetting its
        division to unknown (None)

    - When axis=1, concatenate DataFrames column-wise:

      - Allowed if all divisions are known.

      - If any of division is unknown, it raises ValueError.

    Parameters
    ----------

    dfs : list
        List of dask.DataFrames to be concatenated
    axis : {0, 1, 'index', 'columns'}, default 0
        The axis to concatenate along
    join : {'inner', 'outer'}, default 'outer'
        How to handle indexes on other axis
    interleave_partitions : bool, default False
        Whether to concatenate DataFrames ignoring its order. If True, every
        divisions are concatenated each by each.

    Notes
    -----
    This differs in from ``pd.concat`` in the when concatenating Categoricals
    with different categories. Pandas currently coerces those to objects
    before concatenating. Coercing to objects is very expensive for large
    arrays, so dask preserves the Categoricals by taking the union of
    the categories.

    Examples
    --------

    If all divisions are known and ordered, divisions are kept.

    >>> a                                               # doctest: +SKIP
    dd.DataFrame<x, divisions=(1, 3, 5)>
    >>> b                                               # doctest: +SKIP
    dd.DataFrame<y, divisions=(6, 8, 10)>
    >>> dd.concat([a, b])                               # doctest: +SKIP
    dd.DataFrame<concat-..., divisions=(1, 3, 6, 8, 10)>

    Unable to concatenate if divisions are not ordered.

    >>> a                                               # doctest: +SKIP
    dd.DataFrame<x, divisions=(1, 3, 5)>
    >>> b                                               # doctest: +SKIP
    dd.DataFrame<y, divisions=(2, 3, 6)>
    >>> dd.concat([a, b])                               # doctest: +SKIP
    ValueError: All inputs have known divisions which cannot be concatenated
    in order. Specify interleave_partitions=True to ignore order

    Specify interleave_partitions=True to ignore the division order.

    >>> dd.concat([a, b], interleave_partitions=True)   # doctest: +SKIP
    dd.DataFrame<concat-..., divisions=(1, 2, 3, 5, 6)>

    If any of division is unknown, the result division will be unknown

    >>> a                                               # doctest: +SKIP
    dd.DataFrame<x, divisions=(None, None)>
    >>> b                                               # doctest: +SKIP
    dd.DataFrame<y, divisions=(1, 4, 10)>
    >>> dd.concat([a, b])                               # doctest: +SKIP
    dd.DataFrame<concat-..., divisions=(None, None, None, None)>

    Different categoricals are unioned

    >> dd.concat([                                     # doctest: +SKIP
    ...     dd.from_pandas(pd.Series(['a', 'b'], dtype='category'), 1),
    ...     dd.from_pandas(pd.Series(['a', 'c'], dtype='category'), 1),
    ... ], interleave_partitions=True).dtype
    CategoricalDtype(categories=['a', 'b', 'c'], ordered=False)
    """
    if not isinstance(dfs, list):
        raise TypeError("dfs must be a list of DataFrames/Series objects")
    if len(dfs) == 0:
        raise ValueError('No objects to concatenate')
    if len(dfs) == 1:
        if axis == 1 and isinstance(dfs[0], Series):
            return dfs[0].to_frame()
        else:
            return dfs[0]

    if join not in ('inner', 'outer'):
        raise ValueError("'join' must be 'inner' or 'outer'")

    axis = DataFrame._validate_axis(axis)
    dasks = [df for df in dfs if isinstance(df, _Frame)]
    dfs = _maybe_from_pandas(dfs)

    if axis == 1:
        if all(df.known_divisions for df in dasks):
            return concat_indexed_dataframes(dfs, axis=axis, join=join)
        elif (len(dasks) == len(dfs) and
              all(not df.known_divisions for df in dfs) and
              len({df.npartitions for df in dasks}) == 1):
            warn("Concatenating dataframes with unknown divisions.\n"
                 "We're assuming that the indexes of each dataframes are \n"
                 "aligned. This assumption is not generally safe.")
            return concat_unindexed_dataframes(dfs)
        else:
            raise ValueError('Unable to concatenate DataFrame with unknown '
                             'division specifying axis=1')
    else:
        if all(df.known_divisions for df in dasks):
            # each DataFrame's division must be greater than previous one
            if all(dfs[i].divisions[-1] < dfs[i + 1].divisions[0]
                   for i in range(len(dfs) - 1)):
                divisions = []
                for df in dfs[:-1]:
                    # remove last to concatenate with next
                    divisions += df.divisions[:-1]
                divisions += dfs[-1].divisions
                return stack_partitions(dfs, divisions, join=join)
            elif interleave_partitions:
                return concat_indexed_dataframes(dfs, join=join)
            else:
                raise ValueError('All inputs have known divisions which '
                                 'cannot be concatenated in order. Specify '
                                 'interleave_partitions=True to ignore order')
        else:
            divisions = [None] * (sum([df.npartitions for df in dfs]) + 1)
            return stack_partitions(dfs, divisions, join=join)
