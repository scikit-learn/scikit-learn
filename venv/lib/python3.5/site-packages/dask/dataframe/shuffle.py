from __future__ import absolute_import, division, print_function

import math
from operator import getitem
import uuid

import numpy as np
import pandas as pd
from toolz import merge

from .methods import drop_columns
from .core import DataFrame, Series, _Frame, _concat, map_partitions
from .hashing import hash_pandas_object
from .utils import PANDAS_VERSION

from .. import base
from ..base import tokenize, compute, compute_as_if_collection
from ..context import _globals
from ..delayed import delayed
from ..sizeof import sizeof
from ..utils import digit, insert, M

if PANDAS_VERSION >= '0.20.0':
    from pandas._libs.algos import groupsort_indexer
else:
    from pandas.algos import groupsort_indexer


def set_index(df, index, npartitions=None, shuffle=None, compute=False,
              drop=True, upsample=1.0, divisions=None,
              partition_size=128e6, **kwargs):
    """ See _Frame.set_index for docstring """
    if (isinstance(index, Series) and index._name == df.index._name):
        return df
    if isinstance(index, (DataFrame, tuple, list)):
        raise NotImplementedError(
            "Dask dataframe does not yet support multi-indexes.\n"
            "You tried to index with this index: %s\n"
            "Indexes must be single columns only." % str(index))

    if npartitions == 'auto':
        repartition = True
        npartitions = max(100, df.npartitions)
    else:
        if npartitions is None:
            npartitions = df.npartitions
        repartition = False

    if not isinstance(index, Series):
        index2 = df[index]
    else:
        index2 = index

    if divisions is None:
        divisions = index2._repartition_quantiles(npartitions, upsample=upsample)
        if repartition:
            parts = df.to_delayed()
            sizes = [delayed(sizeof)(part) for part in parts]
        else:
            sizes = []
        iparts = index2.to_delayed()
        mins = [ipart.min() for ipart in iparts]
        maxes = [ipart.max() for ipart in iparts]
        divisions, sizes, mins, maxes = base.compute(divisions, sizes, mins, maxes)
        divisions = divisions.tolist()

        empty_dataframe_detected = pd.isnull(divisions).all()
        if repartition or empty_dataframe_detected:
            total = sum(sizes)
            npartitions = max(math.ceil(total / partition_size), 1)
            npartitions = min(npartitions, df.npartitions)
            n = len(divisions)
            try:
                divisions = np.interp(x=np.linspace(0, n - 1, npartitions + 1),
                                      xp=np.linspace(0, n - 1, n),
                                      fp=divisions).tolist()
            except (TypeError, ValueError):  # str type
                indexes = np.linspace(0, n - 1, npartitions + 1).astype(int)
                divisions = [divisions[i] for i in indexes]

        mins = remove_nans(mins)
        maxes = remove_nans(maxes)

        if (mins == sorted(mins) and maxes == sorted(maxes) and
                all(mx < mn for mx, mn in zip(maxes[:-1], mins[1:]))):
            divisions = mins + [maxes[-1]]
            result = set_sorted_index(df, index, drop=drop, divisions=divisions)
            # There are cases where this still may not be sorted
            # so sort_index to be sure. https://github.com/dask/dask/issues/2288
            return result.map_partitions(M.sort_index)

    return set_partition(df, index, divisions, shuffle=shuffle, drop=drop,
                         compute=compute, **kwargs)


def remove_nans(divisions):
    """ Remove nans from divisions

    These sometime pop up when we call min/max on an empty partition

    Examples
    --------
    >>> remove_nans((np.nan, 1, 2))
    [1, 1, 2]
    >>> remove_nans((1, np.nan, 2))
    [1, 2, 2]
    >>> remove_nans((1, 2, np.nan))
    [1, 2, 2]
    """
    divisions = list(divisions)

    for i in range(len(divisions) - 2, -1, -1):
        if pd.isnull(divisions[i]):
            divisions[i] = divisions[i + 1]

    for i in range(len(divisions) - 1, -1, -1):
        if not pd.isnull(divisions[i]):
            for j in range(i + 1, len(divisions)):
                divisions[j] = divisions[i]
            break

    return divisions


def set_partition(df, index, divisions, max_branch=32, drop=True, shuffle=None,
                  compute=None):
    """ Group DataFrame by index

    Sets a new index and partitions data along that index according to
    divisions.  Divisions are often found by computing approximate quantiles.
    The function ``set_index`` will do both of these steps.

    Parameters
    ----------
    df: DataFrame/Series
        Data that we want to re-partition
    index: string or Series
        Column to become the new index
    divisions: list
        Values to form new divisions between partitions
    drop: bool, default True
        Whether to delete columns to be used as the new index
    shuffle: str (optional)
        Either 'disk' for an on-disk shuffle or 'tasks' to use the task
        scheduling framework.  Use 'disk' if you are on a single machine
        and 'tasks' if you are on a distributed cluster.
    max_branch: int (optional)
        If using the task-based shuffle, the amount of splitting each
        partition undergoes.  Increase this for fewer copies but more
        scheduler overhead.

    See Also
    --------
    set_index
    shuffle
    partd
    """
    if np.isscalar(index):
        partitions = df[index].map_partitions(set_partitions_pre,
                                              divisions=divisions,
                                              meta=pd.Series([0]))
        df2 = df.assign(_partitions=partitions)
    else:
        partitions = index.map_partitions(set_partitions_pre,
                                          divisions=divisions,
                                          meta=pd.Series([0]))
        df2 = df.assign(_partitions=partitions, _index=index)

    df3 = rearrange_by_column(df2, '_partitions', max_branch=max_branch,
                              npartitions=len(divisions) - 1, shuffle=shuffle,
                              compute=compute)

    if np.isscalar(index):
        df4 = df3.map_partitions(set_index_post_scalar, index_name=index,
                                 drop=drop, column_dtype=df.columns.dtype)
    else:
        df4 = df3.map_partitions(set_index_post_series, index_name=index.name,
                                 drop=drop, column_dtype=df.columns.dtype)

    df4.divisions = divisions

    return df4.map_partitions(M.sort_index)


def shuffle(df, index, shuffle=None, npartitions=None, max_branch=32,
            compute=None):
    """ Group DataFrame by index

    Hash grouping of elements. After this operation all elements that have
    the same index will be in the same partition. Note that this requires
    full dataset read, serialization and shuffle. This is expensive. If
    possible you should avoid shuffles.

    This does not preserve a meaningful index/partitioning scheme. This is not
    deterministic if done in parallel.

    See Also
    --------
    set_index
    set_partition
    shuffle_disk
    shuffle_tasks
    """
    if not isinstance(index, _Frame):
        index = df._select_columns_or_index(index)

    partitions = index.map_partitions(partitioning_index,
                                      npartitions=npartitions or df.npartitions,
                                      meta=pd.Series([0]))
    df2 = df.assign(_partitions=partitions)
    df3 = rearrange_by_column(df2, '_partitions', npartitions=npartitions,
                              max_branch=max_branch, shuffle=shuffle,
                              compute=compute)
    df4 = df3.map_partitions(drop_columns, '_partitions', df.columns.dtype)
    return df4


def rearrange_by_divisions(df, column, divisions, max_branch=None, shuffle=None):
    """ Shuffle dataframe so that column separates along divisions """
    partitions = df[column].map_partitions(set_partitions_pre,
                                           divisions=divisions,
                                           meta=pd.Series([0]))
    df2 = df.assign(_partitions=partitions)
    df3 = rearrange_by_column(df2, '_partitions', max_branch=max_branch,
                              npartitions=len(divisions) - 1, shuffle=shuffle)
    df4 = df3.map_partitions(drop_columns, '_partitions', df.columns.dtype)
    return df4


def rearrange_by_column(df, col, npartitions=None, max_branch=None,
                        shuffle=None, compute=None):
    shuffle = shuffle or _globals.get('shuffle', 'disk')
    if shuffle == 'disk':
        return rearrange_by_column_disk(df, col, npartitions, compute=compute)
    elif shuffle == 'tasks':
        return rearrange_by_column_tasks(df, col, max_branch, npartitions)
    else:
        raise NotImplementedError("Unknown shuffle method %s" % shuffle)


class maybe_buffered_partd(object):
    """If serialized, will return non-buffered partd. Otherwise returns a
    buffered partd"""
    def __init__(self, buffer=True, tempdir=None):
        self.tempdir = tempdir or _globals.get('temporary_directory')
        self.buffer = buffer

    def __reduce__(self):
        if self.tempdir:
            return (maybe_buffered_partd, (False, self.tempdir))
        else:
            return (maybe_buffered_partd, (False,))

    def __call__(self, *args, **kwargs):
        import partd
        if self.tempdir:
            file = partd.File(dir=self.tempdir)
        else:
            file = partd.File()
        if self.buffer:
            return partd.PandasBlocks(partd.Buffer(partd.Dict(), file))
        else:
            return partd.PandasBlocks(file)


def rearrange_by_column_disk(df, column, npartitions=None, compute=False):
    """ Shuffle using local disk """
    if npartitions is None:
        npartitions = df.npartitions

    token = tokenize(df, column, npartitions)
    always_new_token = uuid.uuid1().hex

    p = ('zpartd-' + always_new_token,)
    dsk1 = {p: (maybe_buffered_partd(),)}

    # Partition data on disk
    name = 'shuffle-partition-' + always_new_token
    dsk2 = {(name, i): (shuffle_group_3, key, column, npartitions, p)
            for i, key in enumerate(df.__dask_keys__())}

    dsk = merge(df.dask, dsk1, dsk2)
    if compute:
        keys = [p, sorted(dsk2)]
        pp, values = compute_as_if_collection(DataFrame, dsk, keys)
        dsk1 = {p: pp}
        dsk = dict(zip(sorted(dsk2), values))

    # Barrier
    barrier_token = 'barrier-' + always_new_token
    dsk3 = {barrier_token: (barrier, list(dsk2))}

    # Collect groups
    name = 'shuffle-collect-' + token
    dsk4 = {(name, i): (collect, p, i, df._meta, barrier_token)
            for i in range(npartitions)}

    divisions = (None,) * (npartitions + 1)

    dsk = merge(dsk, dsk1, dsk3, dsk4)

    return DataFrame(dsk, name, df._meta, divisions)


def rearrange_by_column_tasks(df, column, max_branch=32, npartitions=None):
    """ Order divisions of DataFrame so that all values within column align

    This enacts a task-based shuffle

    See also:
        rearrange_by_column_disk
        set_partitions_tasks
        shuffle_tasks
    """
    max_branch = max_branch or 32
    n = df.npartitions

    stages = int(math.ceil(math.log(n) / math.log(max_branch)))
    if stages > 1:
        k = int(math.ceil(n ** (1 / stages)))
    else:
        k = n

    groups = []
    splits = []
    joins = []

    inputs = [tuple(digit(i, j, k) for j in range(stages))
              for i in range(k**stages)]

    token = tokenize(df, column, max_branch)

    start = dict((('shuffle-join-' + token, 0, inp),
                  (df._name, i) if i < df.npartitions else df._meta)
                 for i, inp in enumerate(inputs))

    for stage in range(1, stages + 1):
        group = dict((('shuffle-group-' + token, stage, inp),
                      (shuffle_group, ('shuffle-join-' + token, stage - 1, inp),
                       column, stage - 1, k, n))
                     for inp in inputs)

        split = dict((('shuffle-split-' + token, stage, i, inp),
                      (getitem, ('shuffle-group-' + token, stage, inp), i))
                     for i in range(k)
                     for inp in inputs)

        join = dict((('shuffle-join-' + token, stage, inp),
                     (_concat,
                      [('shuffle-split-' + token, stage, inp[stage - 1],
                       insert(inp, stage - 1, j)) for j in range(k)]))
                    for inp in inputs)
        groups.append(group)
        splits.append(split)
        joins.append(join)

    end = dict((('shuffle-' + token, i),
                ('shuffle-join-' + token, stages, inp))
               for i, inp in enumerate(inputs))

    dsk = merge(df.dask, start, end, *(groups + splits + joins))
    df2 = DataFrame(dsk, 'shuffle-' + token, df, df.divisions)

    if npartitions is not None and npartitions != df.npartitions:
        parts = [i % df.npartitions for i in range(npartitions)]
        token = tokenize(df2, npartitions)

        dsk = {('repartition-group-' + token, i): (shuffle_group_2, k, column)
               for i, k in enumerate(df2.__dask_keys__())}
        for p in range(npartitions):
            dsk[('repartition-get-' + token, p)] = \
                (shuffle_group_get, ('repartition-group-' + token, parts[p]), p)

        df3 = DataFrame(merge(df2.dask, dsk), 'repartition-get-' + token, df2,
                        [None] * (npartitions + 1))
    else:
        df3 = df2
        df3.divisions = (None,) * (df.npartitions + 1)

    return df3


########################################################
# Various convenience functions to be run by the above #
########################################################


def partitioning_index(df, npartitions):
    """
    Computes a deterministic index mapping each record to a partition.

    Identical rows are mapped to the same partition.

    Parameters
    ----------
    df : DataFrame/Series/Index
    npartitions : int
        The number of partitions to group into.

    Returns
    -------
    partitions : ndarray
        An array of int64 values mapping each record to a partition.
    """
    return hash_pandas_object(df, index=False) % int(npartitions)


def barrier(args):
    list(args)
    return 0


def collect(p, part, meta, barrier_token):
    """ Collect partitions from partd, yield dataframes """
    res = p.get(part)
    return res if len(res) > 0 else meta


def set_partitions_pre(s, divisions):
    partitions = pd.Series(divisions).searchsorted(s, side='right') - 1
    partitions[(s >= divisions[-1]).values] = len(divisions) - 2
    return partitions


def shuffle_group_2(df, col):
    if not len(df):
        return {}, df
    ind = df[col]._values.astype(np.int64)
    n = ind.max() + 1
    indexer, locations = groupsort_indexer(ind.view(np.int64), n)
    df2 = df.take(indexer)
    locations = locations.cumsum()
    parts = [df2.iloc[a:b] for a, b in zip(locations[:-1], locations[1:])]
    result2 = dict(zip(range(n), parts))
    return result2, df.iloc[:0]


def shuffle_group_get(g_head, i):
    g, head = g_head
    if i in g:
        return g[i]
    else:
        return head


def shuffle_group(df, col, stage, k, npartitions):
    """ Splits dataframe into groups

    The group is determined by their final partition, and which stage we are in
    in the shuffle
    """
    if col == '_partitions':
        ind = df[col]
    else:
        ind = hash_pandas_object(df[col], index=False)

    c = ind._values
    typ = np.min_scalar_type(npartitions * 2)

    npartitions, k, stage = [np.array(x, dtype=np.min_scalar_type(x))[()]
                             for x in [npartitions, k, stage]]

    c = np.mod(c, npartitions).astype(typ, copy=False)
    c = np.floor_divide(c, k ** stage, out=c)
    c = np.mod(c, k, out=c)

    indexer, locations = groupsort_indexer(c.astype(np.int64), k)
    df2 = df.take(indexer)
    locations = locations.cumsum()
    parts = [df2.iloc[a:b] for a, b in zip(locations[:-1], locations[1:])]

    return dict(zip(range(k), parts))


def shuffle_group_3(df, col, npartitions, p):
    g = df.groupby(col)
    d = {i: g.get_group(i) for i in g.groups}
    p.append(d, fsync=True)


def set_index_post_scalar(df, index_name, drop, column_dtype):
    df2 = df.drop('_partitions', axis=1).set_index(index_name, drop=drop)
    df2.columns = df2.columns.astype(column_dtype)
    return df2


def set_index_post_series(df, index_name, drop, column_dtype):
    df2 = df.drop('_partitions', axis=1).set_index('_index', drop=True)
    df2.index.name = index_name
    df2.columns = df2.columns.astype(column_dtype)
    return df2


def set_sorted_index(df, index, drop=True, divisions=None, **kwargs):
    if not isinstance(index, Series):
        meta = df._meta.set_index(index, drop=drop)
    else:
        meta = df._meta.set_index(index._meta, drop=drop)

    result = map_partitions(M.set_index, df, index, drop=drop, meta=meta)

    if not divisions:
        divisions = compute_divisions(result, **kwargs)
    elif len(divisions) != len(df.divisions):
        msg = ("When doing `df.set_index(col, sorted=True, divisions=...)`, "
               "divisions indicates known splits in the index column. In this "
               "case divisions must be the same length as the existing "
               "divisions in `df`\n\n"
               "If the intent is to repartition into new divisions after "
               "setting the index, you probably want:\n\n"
               "`df.set_index(col, sorted=True).repartition(divisions=divisions)`")
        raise ValueError(msg)

    result.divisions = tuple(divisions)
    return result


def compute_divisions(df, **kwargs):
    mins = df.index.map_partitions(M.min, meta=df.index)
    maxes = df.index.map_partitions(M.max, meta=df.index)
    mins, maxes = compute(mins, maxes, **kwargs)

    if (sorted(mins) != list(mins) or
            sorted(maxes) != list(maxes) or
            any(a > b for a, b in zip(mins, maxes))):
        raise ValueError("Partitions must be sorted ascending with the index",
                         mins, maxes)

    divisions = tuple(mins) + (list(maxes)[-1],)
    return divisions
