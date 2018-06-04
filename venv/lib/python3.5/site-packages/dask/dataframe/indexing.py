from __future__ import absolute_import, division, print_function

from datetime import datetime
from collections import defaultdict

from toolz import merge
import bisect
import numpy as np
import pandas as pd

from .core import new_dd_object, Series
from . import methods
from ..base import tokenize


class _LocIndexer(object):
    """ Helper class for the .loc accessor """

    def __init__(self, obj):
        self.obj = obj

    @property
    def _name(self):
        return self.obj._name

    def _make_meta(self, iindexer, cindexer):
        """
        get metadata
        """
        if cindexer is None:
            return self.obj
        else:
            return self.obj._meta.loc[:, cindexer]

    def __getitem__(self, key):

        if isinstance(key, tuple):
            # multi-dimensional selection
            if len(key) > self.obj.ndim:
                # raise from pandas
                msg = 'Too many indexers'
                raise pd.core.indexing.IndexingError(msg)

            iindexer = key[0]
            cindexer = key[1]
        else:
            # if self.obj is Series, cindexer is always None
            iindexer = key
            cindexer = None
        return self._loc(iindexer, cindexer)

    def _loc(self, iindexer, cindexer):
        """ Helper function for the .loc accessor """
        if isinstance(iindexer, Series):
            return self._loc_series(iindexer, cindexer)

        if self.obj.known_divisions:
            iindexer = self._maybe_partial_time_string(iindexer)

            if isinstance(iindexer, slice):
                return self._loc_slice(iindexer, cindexer)
            elif isinstance(iindexer, (list, np.ndarray)):
                return self._loc_list(iindexer, cindexer)
            else:
                # element should raise KeyError
                return self._loc_element(iindexer, cindexer)
        else:
            if isinstance(iindexer, (list, np.ndarray)):
                # applying map_pattition to each partitions
                # results in duplicated NaN rows
                msg = 'Cannot index with list against unknown division'
                raise KeyError(msg)
            elif not isinstance(iindexer, slice):
                iindexer = slice(iindexer, iindexer)

            meta = self._make_meta(iindexer, cindexer)
            return self.obj.map_partitions(methods.try_loc, iindexer, cindexer,
                                           meta=meta)

    def _maybe_partial_time_string(self, iindexer):
        """
        Convert index-indexer for partial time string slicing
        if obj.index is DatetimeIndex / PeriodIndex
        """
        iindexer = _maybe_partial_time_string(self.obj._meta_nonempty.index,
                                              iindexer, kind='loc')
        return iindexer

    def _loc_series(self, iindexer, cindexer):
        meta = self._make_meta(iindexer, cindexer)
        return self.obj.map_partitions(methods.loc, iindexer, cindexer,
                                       token='loc-series', meta=meta)

    def _loc_list(self, iindexer, cindexer):
        name = 'loc-%s' % tokenize(iindexer, self.obj)
        parts = self._get_partitions(iindexer)
        meta = self._make_meta(iindexer, cindexer)

        if len(iindexer):
            dsk = {}
            divisions = []
            items = sorted(parts.items())
            for i, (div, indexer) in enumerate(items):
                dsk[name, i] = (methods.loc, (self._name, div),
                                indexer, cindexer)
                # append minimum value as division
                divisions.append(sorted(indexer)[0])
            # append maximum value of the last division
            divisions.append(sorted(items[-1][1])[-1])
        else:
            divisions = [None, None]
            dsk = {(name, 0): meta.head(0)}
        return new_dd_object(merge(self.obj.dask, dsk), name,
                             meta=meta, divisions=divisions)

    def _loc_element(self, iindexer, cindexer):
        name = 'loc-%s' % tokenize(iindexer, self.obj)
        part = self._get_partitions(iindexer)

        if iindexer < self.obj.divisions[0] or iindexer > self.obj.divisions[-1]:
            raise KeyError('the label [%s] is not in the index' % str(iindexer))

        dsk = {(name, 0): (methods.loc, (self._name, part),
                           slice(iindexer, iindexer), cindexer)}

        meta = self._make_meta(iindexer, cindexer)
        return new_dd_object(merge(self.obj.dask, dsk), name,
                             meta=meta, divisions=[iindexer, iindexer])

    def _get_partitions(self, keys):
        if isinstance(keys, (list, np.ndarray)):
            return _partitions_of_index_values(self.obj.divisions, keys)
        else:
            # element
            return _partition_of_index_value(self.obj.divisions, keys)

    def _coerce_loc_index(self, key):
        return _coerce_loc_index(self.obj.divisions, key)

    def _loc_slice(self, iindexer, cindexer):
        name = 'loc-%s' % tokenize(iindexer, cindexer, self)

        assert isinstance(iindexer, slice)
        assert iindexer.step in (None, 1)

        if iindexer.start is not None:
            start = self._get_partitions(iindexer.start)
        else:
            start = 0
        if iindexer.stop is not None:
            stop = self._get_partitions(iindexer.stop)
        else:
            stop = self.obj.npartitions - 1

        if iindexer.start is None and self.obj.known_divisions:
            istart = self.obj.divisions[0]
        else:
            istart = self._coerce_loc_index(iindexer.start)
        if iindexer.stop is None and self.obj.known_divisions:
            istop = self.obj.divisions[-1]
        else:
            istop = self._coerce_loc_index(iindexer.stop)

        if stop == start:
            dsk = {(name, 0): (methods.loc, (self._name, start),
                               slice(iindexer.start, iindexer.stop), cindexer)}
            divisions = [istart, istop]
        else:
            dsk = {(name, 0): (methods.loc, (self._name, start),
                               slice(iindexer.start, None), cindexer)}
            for i in range(1, stop - start):
                if cindexer is None:
                    dsk[name, i] = (self._name, start + i)
                else:
                    dsk[name, i] = (methods.loc, (self._name, start + i),
                                    slice(None, None), cindexer)

            dsk[name, stop - start] = (methods.loc, (self._name, stop),
                                       slice(None, iindexer.stop), cindexer)

            if iindexer.start is None:
                div_start = self.obj.divisions[0]
            else:
                div_start = max(istart, self.obj.divisions[start])

            if iindexer.stop is None:
                div_stop = self.obj.divisions[-1]
            else:
                div_stop = min(istop, self.obj.divisions[stop + 1])

            divisions = ((div_start, ) +
                         self.obj.divisions[start + 1:stop + 1] +
                         (div_stop, ))

        assert len(divisions) == len(dsk) + 1

        meta = self._make_meta(iindexer, cindexer)
        return new_dd_object(merge(self.obj.dask, dsk), name,
                             meta=meta, divisions=divisions)


def _partition_of_index_value(divisions, val):
    """ In which partition does this value lie?

    >>> _partition_of_index_value([0, 5, 10], 3)
    0
    >>> _partition_of_index_value([0, 5, 10], 8)
    1
    >>> _partition_of_index_value([0, 5, 10], 100)
    1
    >>> _partition_of_index_value([0, 5, 10], 5)  # left-inclusive divisions
    1
    """
    if divisions[0] is None:
        msg = "Can not use loc on DataFrame without known divisions"
        raise ValueError(msg)
    val = _coerce_loc_index(divisions, val)
    i = bisect.bisect_right(divisions, val)
    return min(len(divisions) - 2, max(0, i - 1))


def _partitions_of_index_values(divisions, values):
    """ Return defaultdict of division and values pairs
    Each key corresponds to the division which values are index values belong
    to the division.

    >>> sorted(_partitions_of_index_values([0, 5, 10], [3]).items())
    [(0, [3])]
    >>> sorted(_partitions_of_index_values([0, 5, 10], [3, 8, 5]).items())
    [(0, [3]), (1, [8, 5])]
    """
    if divisions[0] is None:
        msg = "Can not use loc on DataFrame without known divisions"
        raise ValueError(msg)

    results = defaultdict(list)
    values = pd.Index(values, dtype=object)
    for val in values:
        i = bisect.bisect_right(divisions, val)
        div = min(len(divisions) - 2, max(0, i - 1))
        results[div].append(val)
    return results


def _coerce_loc_index(divisions, o):
    """ Transform values to be comparable against divisions

    This is particularly valuable to use with pandas datetimes
    """
    if divisions and isinstance(divisions[0], datetime):
        return pd.Timestamp(o)
    if divisions and isinstance(divisions[0], np.datetime64):
        return np.datetime64(o).astype(divisions[0].dtype)
    return o


def _maybe_partial_time_string(index, indexer, kind):
    """
    Convert indexer for partial string selection
    if data has DatetimeIndex/PeriodIndex
    """
    # do not pass dd.Index
    assert isinstance(index, pd.Index)

    if not isinstance(index, (pd.DatetimeIndex, pd.PeriodIndex)):
        return indexer

    if isinstance(indexer, slice):
        if isinstance(indexer.start, pd.compat.string_types):
            start = index._maybe_cast_slice_bound(indexer.start, 'left', kind)
        else:
            start = indexer.start

        if isinstance(indexer.stop, pd.compat.string_types):
            stop = index._maybe_cast_slice_bound(indexer.stop, 'right', kind)
        else:
            stop = indexer.stop
        return slice(start, stop)

    elif isinstance(indexer, pd.compat.string_types):
        start = index._maybe_cast_slice_bound(indexer, 'left', 'loc')
        stop = index._maybe_cast_slice_bound(indexer, 'right', 'loc')
        return slice(min(start, stop), max(start, stop))

    return indexer
