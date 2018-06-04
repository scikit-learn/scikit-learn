from __future__ import absolute_import, division, print_function

import datetime

import pandas as pd
from pandas.core.window import Rolling as pd_Rolling

from ..base import tokenize
from ..utils import M, funcname, derived_from
from .core import _emulate
from .utils import make_meta, PANDAS_VERSION


def overlap_chunk(func, prev_part, current_part, next_part, before, after,
                  args, kwargs):

    msg = ("Partition size is less than overlapping "
           "window size. Try using ``df.repartition`` "
           "to increase the partition size.")

    if prev_part is not None and isinstance(before, int):
        if prev_part.shape[0] != before:
            raise NotImplementedError(msg)

    if next_part is not None and isinstance(after, int):
        if next_part.shape[0] != after:
            raise NotImplementedError(msg)
    # We validate that the window isn't too large for tiemdeltas in map_overlap

    parts = [p for p in (prev_part, current_part, next_part) if p is not None]
    combined = pd.concat(parts)
    out = func(combined, *args, **kwargs)
    if prev_part is None:
        before = None
    if isinstance(before, datetime.timedelta):
        before = len(prev_part)

    if next_part is None:
        return out.iloc[before:]
    if isinstance(after, datetime.timedelta):
        after = len(next_part)
    return out.iloc[before:-after]


def map_overlap(func, df, before, after, *args, **kwargs):
    """Apply a function to each partition, sharing rows with adjacent partitions.

    Parameters
    ----------
    func : function
        Function applied to each partition.
    df : dd.DataFrame, dd.Series
    before : int or timedelta
        The rows to prepend to partition ``i`` from the end of
        partition ``i - 1``.
    after : int or timedelta
        The rows to append to partition ``i`` from the beginning
        of partition ``i + 1``.
    args, kwargs :
        Arguments and keywords to pass to the function. The partition will
        be the first argument, and these will be passed *after*.

    See Also
    --------
    dd.DataFrame.map_overlap
    """
    if (isinstance(before, datetime.timedelta) or isinstance(after, datetime.timedelta)):
        if not df.index._meta_nonempty.is_all_dates:
            raise TypeError("Must have a `DatetimeIndex` when using string offset "
                            "for `before` and `after`")
    else:
        if not (isinstance(before, int) and before >= 0 and
                isinstance(after, int) and after >= 0):
            raise ValueError("before and after must be positive integers")

    if 'token' in kwargs:
        func_name = kwargs.pop('token')
        token = tokenize(df, before, after, *args, **kwargs)
    else:
        func_name = 'overlap-' + funcname(func)
        token = tokenize(func, df, before, after, *args, **kwargs)

    if 'meta' in kwargs:
        meta = kwargs.pop('meta')
    else:
        meta = _emulate(func, df, *args, **kwargs)
    meta = make_meta(meta)

    name = '{0}-{1}'.format(func_name, token)
    name_a = 'overlap-prepend-' + tokenize(df, before)
    name_b = 'overlap-append-' + tokenize(df, after)
    df_name = df._name

    dsk = df.dask.copy()

    # Have to do the checks for too large windows in the time-delta case
    # here instead of in `overlap_chunk`, since we can't rely on fix-frequency
    # index

    timedelta_partition_message = (
        "Partition size is less than specified window. "
        "Try using ``df.repartition`` to increase the partition size"
    )

    if before and isinstance(before, int):
        dsk.update({(name_a, i): (M.tail, (df_name, i), before)
                    for i in range(df.npartitions - 1)})
        prevs = [None] + [(name_a, i) for i in range(df.npartitions - 1)]
    elif isinstance(before, datetime.timedelta):
        # Assumes monotonic (increasing?) index
        deltas = pd.Series(df.divisions).diff().iloc[1:-1]
        if (before > deltas).any():
            raise ValueError(timedelta_partition_message)
        dsk.update({(name_a, i): (_tail_timedelta, (df_name, i), (df_name, i + 1), before)
                    for i in range(df.npartitions - 1)})
        prevs = [None] + [(name_a, i) for i in range(df.npartitions - 1)]
    else:
        prevs = [None] * df.npartitions

    if after and isinstance(after, int):
        dsk.update({(name_b, i): (M.head, (df_name, i), after)
                    for i in range(1, df.npartitions)})
        nexts = [(name_b, i) for i in range(1, df.npartitions)] + [None]
    elif isinstance(after, datetime.timedelta):
        # TODO: Do we have a use-case for this? Pandas doesn't allow negative rolling windows
        deltas = pd.Series(df.divisions).diff().iloc[1:-1]
        if (after > deltas).any():
            raise ValueError(timedelta_partition_message)

        dsk.update({(name_b, i): (_head_timedelta, (df_name, i - 0), (df_name, i), after)
                    for i in range(1, df.npartitions)})
        nexts = [(name_b, i) for i in range(1, df.npartitions)] + [None]
    else:
        nexts = [None] * df.npartitions

    for i, (prev, current, next) in enumerate(zip(prevs, df.__dask_keys__(), nexts)):
        dsk[(name, i)] = (overlap_chunk, func, prev, current, next, before,
                          after, args, kwargs)

    return df._constructor(dsk, name, meta, df.divisions)


def _head_timedelta(current, next_, after):
    """Return rows of ``next_`` whose index is before the last
    observation in ``current`` + ``after``.

    Parameters
    ----------
    current : DataFrame
    next_ : DataFrame
    after : timedelta

    Returns
    -------
    overlapped : DataFrame
    """
    return next_[next_.index < (current.index.max() + after)]


def _tail_timedelta(prev, current, before):
    """Return rows of ``prev`` whose index is after the first
    observation in ``current`` - ``before``.

    Parameters
    ----------
    current : DataFrame
    next_ : DataFrame
    before : timedelta

    Returns
    -------
    overlapped : DataFrame
    """
    return prev[prev.index > (current.index.min() - before)]


def pandas_rolling_method(df, rolling_kwargs, name, *args, **kwargs):
    rolling = df.rolling(**rolling_kwargs)
    return getattr(rolling, name)(*args, **kwargs)


class Rolling(object):
    """Provides rolling window calculations."""

    def __init__(self, obj, window=None, min_periods=None, freq=None,
                 center=False, win_type=None, axis=0):
        if freq is not None:
            msg = 'The deprecated freq argument is not supported.'
            raise NotImplementedError(msg)

        self.obj = obj     # dataframe or series
        self.window = window
        self.min_periods = min_periods
        self.center = center
        self.axis = axis
        self.win_type = win_type
        # Allow pandas to raise if appropriate
        pd_roll = obj._meta.rolling(**self._rolling_kwargs())
        # Using .rolling(window='2s'), pandas will convert the
        # offset str to a window in nanoseconds. But pandas doesn't
        # accept the integer window with win_type='freq', so we store
        # that information here.
        # See https://github.com/pandas-dev/pandas/issues/15969
        self._window = pd_roll.window
        self._win_type = pd_roll.win_type
        self._min_periods = pd_roll.min_periods

    def _rolling_kwargs(self):
        return {'window': self.window,
                'min_periods': self.min_periods,
                'center': self.center,
                'win_type': self.win_type,
                'axis': self.axis}

    @property
    def _has_single_partition(self):
        """
        Indicator for whether the object has a single partition (True)
        or multiple (False).
        """
        return (self.axis in (1, 'columns') or
                (isinstance(self.window, int) and self.window <= 1) or
                self.obj.npartitions == 1)

    def _call_method(self, method_name, *args, **kwargs):
        rolling_kwargs = self._rolling_kwargs()
        meta = pandas_rolling_method(self.obj._meta_nonempty, rolling_kwargs,
                                     method_name, *args, **kwargs)

        if self._has_single_partition:
            # There's no overlap just use map_partitions
            return self.obj.map_partitions(pandas_rolling_method,
                                           rolling_kwargs, method_name,
                                           *args, token=method_name, meta=meta,
                                           **kwargs)
        # Convert window to overlap
        if self.center:
            before = self.window // 2
            after = self.window - before - 1
        elif self._win_type == 'freq':
            before = pd.Timedelta(self.window)
            after = 0
        else:
            before = self.window - 1
            after = 0
        return map_overlap(pandas_rolling_method, self.obj, before, after,
                           rolling_kwargs, method_name, *args,
                           token=method_name, meta=meta, **kwargs)

    @derived_from(pd_Rolling)
    def count(self):
        return self._call_method('count')

    @derived_from(pd_Rolling)
    def sum(self):
        return self._call_method('sum')

    @derived_from(pd_Rolling)
    def mean(self):
        return self._call_method('mean')

    @derived_from(pd_Rolling)
    def median(self):
        return self._call_method('median')

    @derived_from(pd_Rolling)
    def min(self):
        return self._call_method('min')

    @derived_from(pd_Rolling)
    def max(self):
        return self._call_method('max')

    @derived_from(pd_Rolling)
    def std(self, ddof=1):
        return self._call_method('std', ddof=1)

    @derived_from(pd_Rolling)
    def var(self, ddof=1):
        return self._call_method('var', ddof=1)

    @derived_from(pd_Rolling)
    def skew(self):
        return self._call_method('skew')

    @derived_from(pd_Rolling)
    def kurt(self):
        return self._call_method('kurt')

    @derived_from(pd_Rolling)
    def quantile(self, quantile):
        return self._call_method('quantile', quantile)

    @derived_from(pd_Rolling)
    def apply(self, func, args=(), kwargs={}, **kwds):
        # TODO: In a future version of pandas this will change to
        # raw=False. Think about inspecting the function signature and setting
        # to that?
        if PANDAS_VERSION >= '0.23.0':
            kwds.setdefault("raw", None)
        else:
            if kwargs:
                msg = ("Invalid argument to 'apply'. Keyword arguments "
                       "should be given as a dict to the 'kwargs' arugment. ")
                raise TypeError(msg)
        return self._call_method('apply', func, args=args,
                                 kwargs=kwargs, **kwds)

    def __repr__(self):

        def order(item):
            k, v = item
            _order = {'window': 0, 'min_periods': 1, 'center': 2,
                      'win_type': 3, 'axis': 4}
            return _order[k]

        rolling_kwargs = self._rolling_kwargs()
        # pandas translates the '2S' offset to nanoseconds
        rolling_kwargs['window'] = self._window
        rolling_kwargs['win_type'] = self._win_type
        return 'Rolling [{}]'.format(','.join(
            '{}={}'.format(k, v)
            for k, v in sorted(rolling_kwargs.items(), key=order)
            if v is not None))
