from datetime import timedelta
import numpy as np
import warnings
import copy
from textwrap import dedent

import pandas as pd
from pandas.core.base import GroupByMixin

from pandas.core.groupby.groupby import (
    BinGrouper, Grouper, _GroupBy, GroupBy, SeriesGroupBy, groupby,
    PanelGroupBy, _pipe_template
)

from pandas.tseries.frequencies import to_offset, is_subperiod, is_superperiod
from pandas.core.indexes.datetimes import DatetimeIndex, date_range
from pandas.core.indexes.timedeltas import TimedeltaIndex
from pandas.tseries.offsets import DateOffset, Tick, Day, delta_to_nanoseconds
from pandas.core.indexes.period import PeriodIndex
from pandas.errors import AbstractMethodError
import pandas.core.algorithms as algos
from pandas.core.dtypes.generic import ABCDataFrame, ABCSeries

import pandas.compat as compat
from pandas.compat.numpy import function as nv

from pandas._libs import lib, tslib
from pandas._libs.tslib import Timestamp
from pandas._libs.tslibs.period import IncompatibleFrequency

from pandas.util._decorators import Appender, Substitution
from pandas.core.generic import _shared_docs
_shared_docs_kwargs = dict()


class Resampler(_GroupBy):

    """
    Class for resampling datetimelike data, a groupby-like operation.
    See aggregate, transform, and apply functions on this object.

    It's easiest to use obj.resample(...) to use Resampler.

    Parameters
    ----------
    obj : pandas object
    groupby : a TimeGrouper object
    axis : int, default 0
    kind : str or None
        'period', 'timestamp' to override default index treatement

    Notes
    -----
    After resampling, see aggregate, apply, and transform functions.

    Returns
    -------
    a Resampler of the appropriate type
    """

    # to the groupby descriptor
    _attributes = ['freq', 'axis', 'closed', 'label', 'convention',
                   'loffset', 'base', 'kind']

    def __init__(self, obj, groupby=None, axis=0, kind=None, **kwargs):
        self.groupby = groupby
        self.keys = None
        self.sort = True
        self.axis = axis
        self.kind = kind
        self.squeeze = False
        self.group_keys = True
        self.as_index = True
        self.exclusions = set()
        self.binner = None
        self.grouper = None

        if self.groupby is not None:
            self.groupby._set_grouper(self._convert_obj(obj), sort=True)

    def __unicode__(self):
        """ provide a nice str repr of our rolling object """
        attrs = ["{k}={v}".format(k=k, v=getattr(self.groupby, k))
                 for k in self._attributes if
                 getattr(self.groupby, k, None) is not None]
        return "{klass} [{attrs}]".format(klass=self.__class__.__name__,
                                          attrs=', '.join(attrs))

    def __getattr__(self, attr):
        if attr in self._internal_names_set:
            return object.__getattribute__(self, attr)
        if attr in self._attributes:
            return getattr(self.groupby, attr)
        if attr in self.obj:
            return self[attr]

        return object.__getattribute__(self, attr)

    @property
    def obj(self):
        return self.groupby.obj

    @property
    def ax(self):
        return self.groupby.ax

    @property
    def _typ(self):
        """ masquerade for compat as a Series or a DataFrame """
        if isinstance(self._selected_obj, pd.Series):
            return 'series'
        return 'dataframe'

    @property
    def _from_selection(self):
        """ is the resampling from a DataFrame column or MultiIndex level """
        # upsampling and PeriodIndex resampling do not work
        # with selection, this state used to catch and raise an error
        return (self.groupby is not None and
                (self.groupby.key is not None or
                 self.groupby.level is not None))

    def _convert_obj(self, obj):
        """
        provide any conversions for the object in order to correctly handle

        Parameters
        ----------
        obj : the object to be resampled

        Returns
        -------
        obj : converted object
        """
        obj = obj._consolidate()
        return obj

    def _get_binner_for_time(self):
        raise AbstractMethodError(self)

    def _set_binner(self):
        """
        setup our binners
        cache these as we are an immutable object
        """

        if self.binner is None:
            self.binner, self.grouper = self._get_binner()

    def _get_binner(self):
        """
        create the BinGrouper, assume that self.set_grouper(obj)
        has already been called
        """

        binner, bins, binlabels = self._get_binner_for_time()
        bin_grouper = BinGrouper(bins, binlabels, indexer=self.groupby.indexer)
        return binner, bin_grouper

    def _assure_grouper(self):
        """ make sure that we are creating our binner & grouper """
        self._set_binner()

    @Substitution(klass='Resampler',
                  versionadded='.. versionadded:: 0.23.0',
                  examples="""
>>> df = pd.DataFrame({'A': [1, 2, 3, 4]},
...                   index=pd.date_range('2012-08-02', periods=4))
>>> df
            A
2012-08-02  1
2012-08-03  2
2012-08-04  3
2012-08-05  4

To get the difference between each 2-day period's maximum and minimum value in
one pass, you can do

>>> df.resample('2D').pipe(lambda x: x.max() - x.min())
            A
2012-08-02  1
2012-08-04  1""")
    @Appender(_pipe_template)
    def pipe(self, func, *args, **kwargs):
        return super(Resampler, self).pipe(func, *args, **kwargs)

    _agg_doc = dedent("""

    Examples
    --------
    >>> s = Series([1,2,3,4,5],
                    index=pd.date_range('20130101',
                                        periods=5,freq='s'))
    2013-01-01 00:00:00    1
    2013-01-01 00:00:01    2
    2013-01-01 00:00:02    3
    2013-01-01 00:00:03    4
    2013-01-01 00:00:04    5
    Freq: S, dtype: int64

    >>> r = s.resample('2s')
    DatetimeIndexResampler [freq=<2 * Seconds>, axis=0, closed=left,
                            label=left, convention=start, base=0]

    >>> r.agg(np.sum)
    2013-01-01 00:00:00    3
    2013-01-01 00:00:02    7
    2013-01-01 00:00:04    5
    Freq: 2S, dtype: int64

    >>> r.agg(['sum','mean','max'])
                         sum  mean  max
    2013-01-01 00:00:00    3   1.5    2
    2013-01-01 00:00:02    7   3.5    4
    2013-01-01 00:00:04    5   5.0    5

    >>> r.agg({'result' : lambda x: x.mean() / x.std(),
               'total' : np.sum})
                         total    result
    2013-01-01 00:00:00      3  2.121320
    2013-01-01 00:00:02      7  4.949747
    2013-01-01 00:00:04      5       NaN

    See also
    --------
    pandas.DataFrame.groupby.aggregate
    pandas.DataFrame.resample.transform
    pandas.DataFrame.aggregate

    """)

    @Appender(_agg_doc)
    @Appender(_shared_docs['aggregate'] % dict(
        klass='DataFrame',
        versionadded='',
        axis=''))
    def aggregate(self, arg, *args, **kwargs):

        self._set_binner()
        result, how = self._aggregate(arg, *args, **kwargs)
        if result is None:
            result = self._groupby_and_aggregate(arg,
                                                 *args,
                                                 **kwargs)

        result = self._apply_loffset(result)
        return result

    agg = aggregate
    apply = aggregate

    def transform(self, arg, *args, **kwargs):
        """
        Call function producing a like-indexed Series on each group and return
        a Series with the transformed values

        Parameters
        ----------
        func : function
            To apply to each group. Should return a Series with the same index

        Examples
        --------
        >>> resampled.transform(lambda x: (x - x.mean()) / x.std())

        Returns
        -------
        transformed : Series
        """
        return self._selected_obj.groupby(self.groupby).transform(
            arg, *args, **kwargs)

    def _downsample(self, f):
        raise AbstractMethodError(self)

    def _upsample(self, f, limit=None, fill_value=None):
        raise AbstractMethodError(self)

    def _gotitem(self, key, ndim, subset=None):
        """
        sub-classes to define
        return a sliced object

        Parameters
        ----------
        key : string / list of selections
        ndim : 1,2
            requested ndim of result
        subset : object, default None
            subset to act on
        """
        self._set_binner()
        grouper = self.grouper
        if subset is None:
            subset = self.obj
        grouped = groupby(subset, by=None, grouper=grouper, axis=self.axis)

        # try the key selection
        try:
            return grouped[key]
        except KeyError:
            return grouped

    def _groupby_and_aggregate(self, how, grouper=None, *args, **kwargs):
        """ re-evaluate the obj with a groupby aggregation """

        if grouper is None:
            self._set_binner()
            grouper = self.grouper

        obj = self._selected_obj

        try:
            grouped = groupby(obj, by=None, grouper=grouper, axis=self.axis)
        except TypeError:

            # panel grouper
            grouped = PanelGroupBy(obj, grouper=grouper, axis=self.axis)

        try:
            if isinstance(obj, ABCDataFrame) and compat.callable(how):
                # Check if the function is reducing or not.
                result = grouped._aggregate_item_by_item(how, *args, **kwargs)
            else:
                result = grouped.aggregate(how, *args, **kwargs)
        except Exception:

            # we have a non-reducing function
            # try to evaluate
            result = grouped.apply(how, *args, **kwargs)

        result = self._apply_loffset(result)
        return self._wrap_result(result)

    def _apply_loffset(self, result):
        """
        if loffset is set, offset the result index

        This is NOT an idempotent routine, it will be applied
        exactly once to the result.

        Parameters
        ----------
        result : Series or DataFrame
            the result of resample
        """

        needs_offset = (
            isinstance(self.loffset, (DateOffset, timedelta)) and
            isinstance(result.index, DatetimeIndex) and
            len(result.index) > 0
        )

        if needs_offset:
            result.index = result.index + self.loffset

        self.loffset = None
        return result

    def _get_resampler_for_grouping(self, groupby, **kwargs):
        """ return the correct class for resampling with groupby """
        return self._resampler_for_grouping(self, groupby=groupby, **kwargs)

    def _wrap_result(self, result):
        """ potentially wrap any results """
        if isinstance(result, ABCSeries) and self._selection is not None:
            result.name = self._selection

        if isinstance(result, ABCSeries) and result.empty:
            obj = self.obj
            result.index = obj.index._shallow_copy(freq=to_offset(self.freq))
            result.name = getattr(obj, 'name', None)

        return result

    def pad(self, limit=None):
        """
        Forward fill the values

        Parameters
        ----------
        limit : integer, optional
            limit of how many values to fill

        Returns
        -------
        an upsampled Series

        See Also
        --------
        Series.fillna
        DataFrame.fillna
        """
        return self._upsample('pad', limit=limit)
    ffill = pad

    def nearest(self, limit=None):
        """
        Fill values with nearest neighbor starting from center

        Parameters
        ----------
        limit : integer, optional
            limit of how many values to fill

            .. versionadded:: 0.21.0

        Returns
        -------
        an upsampled Series

        See Also
        --------
        Series.fillna
        DataFrame.fillna
        """
        return self._upsample('nearest', limit=limit)

    def backfill(self, limit=None):
        """
        Backward fill the new missing values in the resampled data.

        In statistics, imputation is the process of replacing missing data with
        substituted values [1]_. When resampling data, missing values may
        appear (e.g., when the resampling frequency is higher than the original
        frequency). The backward fill will replace NaN values that appeared in
        the resampled data with the next value in the original sequence.
        Missing values that existed in the orginal data will not be modified.

        Parameters
        ----------
        limit : integer, optional
            Limit of how many values to fill.

        Returns
        -------
        Series, DataFrame
            An upsampled Series or DataFrame with backward filled NaN values.

        See Also
        --------
        bfill : Alias of backfill.
        fillna : Fill NaN values using the specified method, which can be
            'backfill'.
        nearest : Fill NaN values with nearest neighbor starting from center.
        pad : Forward fill NaN values.
        pandas.Series.fillna : Fill NaN values in the Series using the
            specified method, which can be 'backfill'.
        pandas.DataFrame.fillna : Fill NaN values in the DataFrame using the
            specified method, which can be 'backfill'.

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Imputation_(statistics)

        Examples
        --------

        Resampling a Series:

        >>> s = pd.Series([1, 2, 3],
        ...               index=pd.date_range('20180101', periods=3, freq='h'))
        >>> s
        2018-01-01 00:00:00    1
        2018-01-01 01:00:00    2
        2018-01-01 02:00:00    3
        Freq: H, dtype: int64

        >>> s.resample('30min').backfill()
        2018-01-01 00:00:00    1
        2018-01-01 00:30:00    2
        2018-01-01 01:00:00    2
        2018-01-01 01:30:00    3
        2018-01-01 02:00:00    3
        Freq: 30T, dtype: int64

        >>> s.resample('15min').backfill(limit=2)
        2018-01-01 00:00:00    1.0
        2018-01-01 00:15:00    NaN
        2018-01-01 00:30:00    2.0
        2018-01-01 00:45:00    2.0
        2018-01-01 01:00:00    2.0
        2018-01-01 01:15:00    NaN
        2018-01-01 01:30:00    3.0
        2018-01-01 01:45:00    3.0
        2018-01-01 02:00:00    3.0
        Freq: 15T, dtype: float64

        Resampling a DataFrame that has missing values:

        >>> df = pd.DataFrame({'a': [2, np.nan, 6], 'b': [1, 3, 5]},
        ...                   index=pd.date_range('20180101', periods=3,
        ...                                       freq='h'))
        >>> df
                               a  b
        2018-01-01 00:00:00  2.0  1
        2018-01-01 01:00:00  NaN  3
        2018-01-01 02:00:00  6.0  5

        >>> df.resample('30min').backfill()
                               a  b
        2018-01-01 00:00:00  2.0  1
        2018-01-01 00:30:00  NaN  3
        2018-01-01 01:00:00  NaN  3
        2018-01-01 01:30:00  6.0  5
        2018-01-01 02:00:00  6.0  5

        >>> df.resample('15min').backfill(limit=2)
                               a    b
        2018-01-01 00:00:00  2.0  1.0
        2018-01-01 00:15:00  NaN  NaN
        2018-01-01 00:30:00  NaN  3.0
        2018-01-01 00:45:00  NaN  3.0
        2018-01-01 01:00:00  NaN  3.0
        2018-01-01 01:15:00  NaN  NaN
        2018-01-01 01:30:00  6.0  5.0
        2018-01-01 01:45:00  6.0  5.0
        2018-01-01 02:00:00  6.0  5.0
        """
        return self._upsample('backfill', limit=limit)
    bfill = backfill

    def fillna(self, method, limit=None):
        """
        Fill missing values introduced by upsampling.

        In statistics, imputation is the process of replacing missing data with
        substituted values [1]_. When resampling data, missing values may
        appear (e.g., when the resampling frequency is higher than the original
        frequency).

        Missing values that existed in the orginal data will
        not be modified.

        Parameters
        ----------
        method : {'pad', 'backfill', 'ffill', 'bfill', 'nearest'}
            Method to use for filling holes in resampled data

            * 'pad' or 'ffill': use previous valid observation to fill gap
              (forward fill).
            * 'backfill' or 'bfill': use next valid observation to fill gap.
            * 'nearest': use nearest valid observation to fill gap.

        limit : integer, optional
            Limit of how many consecutive missing values to fill.

        Returns
        -------
        Series or DataFrame
            An upsampled Series or DataFrame with missing values filled.

        See Also
        --------
        backfill : Backward fill NaN values in the resampled data.
        pad : Forward fill NaN values in the resampled data.
        nearest : Fill NaN values in the resampled data
            with nearest neighbor starting from center.
        interpolate : Fill NaN values using interpolation.
        pandas.Series.fillna : Fill NaN values in the Series using the
            specified method, which can be 'bfill' and 'ffill'.
        pandas.DataFrame.fillna : Fill NaN values in the DataFrame using the
            specified method, which can be 'bfill' and 'ffill'.

        Examples
        --------
        Resampling a Series:

        >>> s = pd.Series([1, 2, 3],
        ...               index=pd.date_range('20180101', periods=3, freq='h'))
        >>> s
        2018-01-01 00:00:00    1
        2018-01-01 01:00:00    2
        2018-01-01 02:00:00    3
        Freq: H, dtype: int64

        Without filling the missing values you get:

        >>> s.resample("30min").asfreq()
        2018-01-01 00:00:00    1.0
        2018-01-01 00:30:00    NaN
        2018-01-01 01:00:00    2.0
        2018-01-01 01:30:00    NaN
        2018-01-01 02:00:00    3.0
        Freq: 30T, dtype: float64

        >>> s.resample('30min').fillna("backfill")
        2018-01-01 00:00:00    1
        2018-01-01 00:30:00    2
        2018-01-01 01:00:00    2
        2018-01-01 01:30:00    3
        2018-01-01 02:00:00    3
        Freq: 30T, dtype: int64

        >>> s.resample('15min').fillna("backfill", limit=2)
        2018-01-01 00:00:00    1.0
        2018-01-01 00:15:00    NaN
        2018-01-01 00:30:00    2.0
        2018-01-01 00:45:00    2.0
        2018-01-01 01:00:00    2.0
        2018-01-01 01:15:00    NaN
        2018-01-01 01:30:00    3.0
        2018-01-01 01:45:00    3.0
        2018-01-01 02:00:00    3.0
        Freq: 15T, dtype: float64

        >>> s.resample('30min').fillna("pad")
        2018-01-01 00:00:00    1
        2018-01-01 00:30:00    1
        2018-01-01 01:00:00    2
        2018-01-01 01:30:00    2
        2018-01-01 02:00:00    3
        Freq: 30T, dtype: int64

        >>> s.resample('30min').fillna("nearest")
        2018-01-01 00:00:00    1
        2018-01-01 00:30:00    2
        2018-01-01 01:00:00    2
        2018-01-01 01:30:00    3
        2018-01-01 02:00:00    3
        Freq: 30T, dtype: int64

        Missing values present before the upsampling are not affected.

        >>> sm = pd.Series([1, None, 3],
        ...               index=pd.date_range('20180101', periods=3, freq='h'))
        >>> sm
        2018-01-01 00:00:00    1.0
        2018-01-01 01:00:00    NaN
        2018-01-01 02:00:00    3.0
        Freq: H, dtype: float64

        >>> sm.resample('30min').fillna('backfill')
        2018-01-01 00:00:00    1.0
        2018-01-01 00:30:00    NaN
        2018-01-01 01:00:00    NaN
        2018-01-01 01:30:00    3.0
        2018-01-01 02:00:00    3.0
        Freq: 30T, dtype: float64

        >>> sm.resample('30min').fillna('pad')
        2018-01-01 00:00:00    1.0
        2018-01-01 00:30:00    1.0
        2018-01-01 01:00:00    NaN
        2018-01-01 01:30:00    NaN
        2018-01-01 02:00:00    3.0
        Freq: 30T, dtype: float64

        >>> sm.resample('30min').fillna('nearest')
        2018-01-01 00:00:00    1.0
        2018-01-01 00:30:00    NaN
        2018-01-01 01:00:00    NaN
        2018-01-01 01:30:00    3.0
        2018-01-01 02:00:00    3.0
        Freq: 30T, dtype: float64

        DataFrame resampling is done column-wise. All the same options are
        available.

        >>> df = pd.DataFrame({'a': [2, np.nan, 6], 'b': [1, 3, 5]},
        ...                   index=pd.date_range('20180101', periods=3,
        ...                                       freq='h'))
        >>> df
                               a  b
        2018-01-01 00:00:00  2.0  1
        2018-01-01 01:00:00  NaN  3
        2018-01-01 02:00:00  6.0  5

        >>> df.resample('30min').fillna("bfill")
                               a  b
        2018-01-01 00:00:00  2.0  1
        2018-01-01 00:30:00  NaN  3
        2018-01-01 01:00:00  NaN  3
        2018-01-01 01:30:00  6.0  5
        2018-01-01 02:00:00  6.0  5

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Imputation_(statistics)
        """
        return self._upsample(method, limit=limit)

    @Appender(_shared_docs['interpolate'] % _shared_docs_kwargs)
    def interpolate(self, method='linear', axis=0, limit=None, inplace=False,
                    limit_direction='forward', limit_area=None,
                    downcast=None, **kwargs):
        """
        Interpolate values according to different methods.

        .. versionadded:: 0.18.1
        """
        result = self._upsample(None)
        return result.interpolate(method=method, axis=axis, limit=limit,
                                  inplace=inplace,
                                  limit_direction=limit_direction,
                                  limit_area=limit_area,
                                  downcast=downcast, **kwargs)

    def asfreq(self, fill_value=None):
        """
        return the values at the new freq,
        essentially a reindex

        Parameters
        ----------
        fill_value: scalar, optional
            Value to use for missing values, applied during upsampling (note
            this does not fill NaNs that already were present).

            .. versionadded:: 0.20.0

        See Also
        --------
        Series.asfreq
        DataFrame.asfreq
        """
        return self._upsample('asfreq', fill_value=fill_value)

    def std(self, ddof=1, *args, **kwargs):
        """
        Compute standard deviation of groups, excluding missing values

        Parameters
        ----------
        ddof : integer, default 1
        degrees of freedom
        """
        nv.validate_resampler_func('std', args, kwargs)
        return self._downsample('std', ddof=ddof)

    def var(self, ddof=1, *args, **kwargs):
        """
        Compute variance of groups, excluding missing values

        Parameters
        ----------
        ddof : integer, default 1
        degrees of freedom
        """
        nv.validate_resampler_func('var', args, kwargs)
        return self._downsample('var', ddof=ddof)

    @Appender(GroupBy.size.__doc__)
    def size(self):
        # It's a special case as higher level does return
        # a copy of 0-len objects. GH14962
        result = self._downsample('size')
        if not len(self.ax) and isinstance(self._selected_obj, ABCDataFrame):
            result = pd.Series([], index=result.index, dtype='int64')
        return result


# downsample methods
for method in ['sum', 'prod']:

    def f(self, _method=method, min_count=0, *args, **kwargs):
        nv.validate_resampler_func(_method, args, kwargs)
        return self._downsample(_method, min_count=min_count)
    f.__doc__ = getattr(GroupBy, method).__doc__
    setattr(Resampler, method, f)


# downsample methods
for method in ['min', 'max', 'first', 'last', 'mean', 'sem',
               'median', 'ohlc']:

    def f(self, _method=method, *args, **kwargs):
        nv.validate_resampler_func(_method, args, kwargs)
        return self._downsample(_method)
    f.__doc__ = getattr(GroupBy, method).__doc__
    setattr(Resampler, method, f)

# groupby & aggregate methods
for method in ['count']:
    def f(self, _method=method):
        return self._downsample(_method)
    f.__doc__ = getattr(GroupBy, method).__doc__
    setattr(Resampler, method, f)

# series only methods
for method in ['nunique']:
    def f(self, _method=method):
        return self._downsample(_method)
    f.__doc__ = getattr(SeriesGroupBy, method).__doc__
    setattr(Resampler, method, f)


def _maybe_process_deprecations(r, how=None, fill_method=None, limit=None):
    """ potentially we might have a deprecation warning, show it
    but call the appropriate methods anyhow """

    if how is not None:

        # .resample(..., how='sum')
        if isinstance(how, compat.string_types):
            method = "{0}()".format(how)

            # .resample(..., how=lambda x: ....)
        else:
            method = ".apply(<func>)"

        # if we have both a how and fill_method, then show
        # the following warning
        if fill_method is None:
            warnings.warn("how in .resample() is deprecated\n"
                          "the new syntax is "
                          ".resample(...).{method}".format(
                              method=method),
                          FutureWarning, stacklevel=3)
        r = r.aggregate(how)

    if fill_method is not None:

        # show the prior function call
        method = '.' + method if how is not None else ''

        args = "limit={0}".format(limit) if limit is not None else ""
        warnings.warn("fill_method is deprecated to .resample()\n"
                      "the new syntax is .resample(...){method}"
                      ".{fill_method}({args})".format(
                          method=method,
                          fill_method=fill_method,
                          args=args),
                      FutureWarning, stacklevel=3)

        if how is not None:
            r = getattr(r, fill_method)(limit=limit)
        else:
            r = r.aggregate(fill_method, limit=limit)

    return r


class _GroupByMixin(GroupByMixin):
    """ provide the groupby facilities """

    def __init__(self, obj, *args, **kwargs):

        parent = kwargs.pop('parent', None)
        groupby = kwargs.pop('groupby', None)
        if parent is None:
            parent = obj

        # initialize our GroupByMixin object with
        # the resampler attributes
        for attr in self._attributes:
            setattr(self, attr, kwargs.get(attr, getattr(parent, attr)))

        super(_GroupByMixin, self).__init__(None)
        self._groupby = groupby
        self._groupby.mutated = True
        self._groupby.grouper.mutated = True
        self.groupby = copy.copy(parent.groupby)

    def _apply(self, f, **kwargs):
        """
        dispatch to _upsample; we are stripping all of the _upsample kwargs and
        performing the original function call on the grouped object
        """

        def func(x):
            x = self._shallow_copy(x, groupby=self.groupby)

            if isinstance(f, compat.string_types):
                return getattr(x, f)(**kwargs)

            return x.apply(f, **kwargs)

        result = self._groupby.apply(func)
        return self._wrap_result(result)

    _upsample = _apply
    _downsample = _apply
    _groupby_and_aggregate = _apply


class DatetimeIndexResampler(Resampler):

    @property
    def _resampler_for_grouping(self):
        return DatetimeIndexResamplerGroupby

    def _get_binner_for_time(self):

        # this is how we are actually creating the bins
        if self.kind == 'period':
            return self.groupby._get_time_period_bins(self.ax)
        return self.groupby._get_time_bins(self.ax)

    def _downsample(self, how, **kwargs):
        """
        Downsample the cython defined function

        Parameters
        ----------
        how : string / cython mapped function
        **kwargs : kw args passed to how function
        """
        self._set_binner()
        how = self._is_cython_func(how) or how
        ax = self.ax
        obj = self._selected_obj

        if not len(ax):
            # reset to the new freq
            obj = obj.copy()
            obj.index.freq = self.freq
            return obj

        # do we have a regular frequency
        if ax.freq is not None or ax.inferred_freq is not None:

            if len(self.grouper.binlabels) > len(ax) and how is None:

                # let's do an asfreq
                return self.asfreq()

        # we are downsampling
        # we want to call the actual grouper method here
        result = obj.groupby(
            self.grouper, axis=self.axis).aggregate(how, **kwargs)

        result = self._apply_loffset(result)
        return self._wrap_result(result)

    def _adjust_binner_for_upsample(self, binner):
        """ adjust our binner when upsampling """
        if self.closed == 'right':
            binner = binner[1:]
        else:
            binner = binner[:-1]
        return binner

    def _upsample(self, method, limit=None, fill_value=None):
        """
        method : string {'backfill', 'bfill', 'pad',
            'ffill', 'asfreq'} method for upsampling
        limit : int, default None
            Maximum size gap to fill when reindexing
        fill_value : scalar, default None
            Value to use for missing values

        See also
        --------
        .fillna

        """
        self._set_binner()
        if self.axis:
            raise AssertionError('axis must be 0')
        if self._from_selection:
            raise ValueError("Upsampling from level= or on= selection"
                             " is not supported, use .set_index(...)"
                             " to explicitly set index to"
                             " datetime-like")

        ax = self.ax
        obj = self._selected_obj
        binner = self.binner
        res_index = self._adjust_binner_for_upsample(binner)

        # if we have the same frequency as our axis, then we are equal sampling
        if limit is None and to_offset(ax.inferred_freq) == self.freq:
            result = obj.copy()
            result.index = res_index
        else:
            result = obj.reindex(res_index, method=method,
                                 limit=limit, fill_value=fill_value)

        result = self._apply_loffset(result)
        return self._wrap_result(result)

    def _wrap_result(self, result):
        result = super(DatetimeIndexResampler, self)._wrap_result(result)

        # we may have a different kind that we were asked originally
        # convert if needed
        if self.kind == 'period' and not isinstance(result.index, PeriodIndex):
            result.index = result.index.to_period(self.freq)
        return result


class DatetimeIndexResamplerGroupby(_GroupByMixin, DatetimeIndexResampler):
    """
    Provides a resample of a groupby implementation

    .. versionadded:: 0.18.1

    """
    @property
    def _constructor(self):
        return DatetimeIndexResampler


class PeriodIndexResampler(DatetimeIndexResampler):

    @property
    def _resampler_for_grouping(self):
        return PeriodIndexResamplerGroupby

    def _get_binner_for_time(self):
        if self.kind == 'timestamp':
            return super(PeriodIndexResampler, self)._get_binner_for_time()
        return self.groupby._get_period_bins(self.ax)

    def _convert_obj(self, obj):
        obj = super(PeriodIndexResampler, self)._convert_obj(obj)

        if self._from_selection:
            # see GH 14008, GH 12871
            msg = ("Resampling from level= or on= selection"
                   " with a PeriodIndex is not currently supported,"
                   " use .set_index(...) to explicitly set index")
            raise NotImplementedError(msg)

        if self.loffset is not None:
            # Cannot apply loffset/timedelta to PeriodIndex -> convert to
            # timestamps
            self.kind = 'timestamp'

        # convert to timestamp
        if self.kind == 'timestamp':
            obj = obj.to_timestamp(how=self.convention)

        return obj

    def _downsample(self, how, **kwargs):
        """
        Downsample the cython defined function

        Parameters
        ----------
        how : string / cython mapped function
        **kwargs : kw args passed to how function
        """

        # we may need to actually resample as if we are timestamps
        if self.kind == 'timestamp':
            return super(PeriodIndexResampler, self)._downsample(how, **kwargs)

        how = self._is_cython_func(how) or how
        ax = self.ax

        if is_subperiod(ax.freq, self.freq):
            # Downsampling
            return self._groupby_and_aggregate(how, grouper=self.grouper)
        elif is_superperiod(ax.freq, self.freq):
            if how == 'ohlc':
                # GH #13083
                # upsampling to subperiods is handled as an asfreq, which works
                # for pure aggregating/reducing methods
                # OHLC reduces along the time dimension, but creates multiple
                # values for each period -> handle by _groupby_and_aggregate()
                return self._groupby_and_aggregate(how, grouper=self.grouper)
            return self.asfreq()
        elif ax.freq == self.freq:
            return self.asfreq()

        raise IncompatibleFrequency(
            'Frequency {} cannot be resampled to {}, as they are not '
            'sub or super periods'.format(ax.freq, self.freq))

    def _upsample(self, method, limit=None, fill_value=None):
        """
        method : string {'backfill', 'bfill', 'pad', 'ffill'}
            method for upsampling
        limit : int, default None
            Maximum size gap to fill when reindexing
        fill_value : scalar, default None
            Value to use for missing values

        See also
        --------
        .fillna

        """

        # we may need to actually resample as if we are timestamps
        if self.kind == 'timestamp':
            return super(PeriodIndexResampler, self)._upsample(
                method, limit=limit, fill_value=fill_value)

        self._set_binner()
        ax = self.ax
        obj = self.obj
        new_index = self.binner

        # Start vs. end of period
        memb = ax.asfreq(self.freq, how=self.convention)

        # Get the fill indexer
        indexer = memb.get_indexer(new_index, method=method, limit=limit)
        return self._wrap_result(_take_new_index(
            obj, indexer, new_index, axis=self.axis))


class PeriodIndexResamplerGroupby(_GroupByMixin, PeriodIndexResampler):
    """
    Provides a resample of a groupby implementation

    .. versionadded:: 0.18.1

    """
    @property
    def _constructor(self):
        return PeriodIndexResampler


class TimedeltaIndexResampler(DatetimeIndexResampler):

    @property
    def _resampler_for_grouping(self):
        return TimedeltaIndexResamplerGroupby

    def _get_binner_for_time(self):
        return self.groupby._get_time_delta_bins(self.ax)

    def _adjust_binner_for_upsample(self, binner):
        """ adjust our binner when upsampling """
        ax = self.ax

        if is_subperiod(ax.freq, self.freq):
            # We are actually downsampling
            # but are in the asfreq path
            # GH 12926
            if self.closed == 'right':
                binner = binner[1:]
            else:
                binner = binner[:-1]
        return binner


class TimedeltaIndexResamplerGroupby(_GroupByMixin, TimedeltaIndexResampler):
    """
    Provides a resample of a groupby implementation

    .. versionadded:: 0.18.1

    """
    @property
    def _constructor(self):
        return TimedeltaIndexResampler


def resample(obj, kind=None, **kwds):
    """ create a TimeGrouper and return our resampler """
    tg = TimeGrouper(**kwds)
    return tg._get_resampler(obj, kind=kind)


resample.__doc__ = Resampler.__doc__


def get_resampler_for_grouping(groupby, rule, how=None, fill_method=None,
                               limit=None, kind=None, **kwargs):
    """ return our appropriate resampler when grouping as well """

    # .resample uses 'on' similar to how .groupby uses 'key'
    kwargs['key'] = kwargs.pop('on', None)

    tg = TimeGrouper(freq=rule, **kwargs)
    resampler = tg._get_resampler(groupby.obj, kind=kind)
    r = resampler._get_resampler_for_grouping(groupby=groupby)
    return _maybe_process_deprecations(r,
                                       how=how,
                                       fill_method=fill_method,
                                       limit=limit)


class TimeGrouper(Grouper):
    """
    Custom groupby class for time-interval grouping

    Parameters
    ----------
    freq : pandas date offset or offset alias for identifying bin edges
    closed : closed end of interval; 'left' or 'right'
    label : interval boundary to use for labeling; 'left' or 'right'
    convention : {'start', 'end', 'e', 's'}
        If axis is PeriodIndex
    """
    _attributes = Grouper._attributes + ('closed', 'label', 'how',
                                         'loffset', 'kind', 'convention',
                                         'base')

    def __init__(self, freq='Min', closed=None, label=None, how='mean',
                 axis=0, fill_method=None, limit=None, loffset=None,
                 kind=None, convention=None, base=0, **kwargs):
        # Check for correctness of the keyword arguments which would
        # otherwise silently use the default if misspelled
        if label not in {None, 'left', 'right'}:
            raise ValueError('Unsupported value {} for `label`'.format(label))
        if closed not in {None, 'left', 'right'}:
            raise ValueError('Unsupported value {} for `closed`'.format(
                closed))
        if convention not in {None, 'start', 'end', 'e', 's'}:
            raise ValueError('Unsupported value {} for `convention`'
                             .format(convention))

        freq = to_offset(freq)

        end_types = set(['M', 'A', 'Q', 'BM', 'BA', 'BQ', 'W'])
        rule = freq.rule_code
        if (rule in end_types or
                ('-' in rule and rule[:rule.find('-')] in end_types)):
            if closed is None:
                closed = 'right'
            if label is None:
                label = 'right'
        else:
            if closed is None:
                closed = 'left'
            if label is None:
                label = 'left'

        self.closed = closed
        self.label = label
        self.kind = kind

        self.convention = convention or 'E'
        self.convention = self.convention.lower()

        if isinstance(loffset, compat.string_types):
            loffset = to_offset(loffset)
        self.loffset = loffset

        self.how = how
        self.fill_method = fill_method
        self.limit = limit
        self.base = base

        # always sort time groupers
        kwargs['sort'] = True

        super(TimeGrouper, self).__init__(freq=freq, axis=axis, **kwargs)

    def _get_resampler(self, obj, kind=None):
        """
        return my resampler or raise if we have an invalid axis

        Parameters
        ----------
        obj : input object
        kind : string, optional
            'period','timestamp','timedelta' are valid

        Returns
        -------
        a Resampler

        Raises
        ------
        TypeError if incompatible axis

        """
        self._set_grouper(obj)

        ax = self.ax
        if isinstance(ax, DatetimeIndex):
            return DatetimeIndexResampler(obj,
                                          groupby=self,
                                          kind=kind,
                                          axis=self.axis)
        elif isinstance(ax, PeriodIndex) or kind == 'period':
            return PeriodIndexResampler(obj,
                                        groupby=self,
                                        kind=kind,
                                        axis=self.axis)
        elif isinstance(ax, TimedeltaIndex):
            return TimedeltaIndexResampler(obj,
                                           groupby=self,
                                           axis=self.axis)

        raise TypeError("Only valid with DatetimeIndex, "
                        "TimedeltaIndex or PeriodIndex, "
                        "but got an instance of %r" % type(ax).__name__)

    def _get_grouper(self, obj, validate=True):
        # create the resampler and return our binner
        r = self._get_resampler(obj)
        r._set_binner()
        return r.binner, r.grouper, r.obj

    def _get_time_bins(self, ax):
        if not isinstance(ax, DatetimeIndex):
            raise TypeError('axis must be a DatetimeIndex, but got '
                            'an instance of %r' % type(ax).__name__)

        if len(ax) == 0:
            binner = labels = DatetimeIndex(
                data=[], freq=self.freq, name=ax.name)
            return binner, [], labels

        first, last = ax.min(), ax.max()
        first, last = _get_range_edges(first, last, self.freq,
                                       closed=self.closed,
                                       base=self.base)
        tz = ax.tz
        # GH #12037
        # use first/last directly instead of call replace() on them
        # because replace() will swallow the nanosecond part
        # thus last bin maybe slightly before the end if the end contains
        # nanosecond part and lead to `Values falls after last bin` error
        binner = labels = DatetimeIndex(freq=self.freq,
                                        start=first,
                                        end=last,
                                        tz=tz,
                                        name=ax.name)

        # GH 15549
        # In edge case of tz-aware resapmling binner last index can be
        # less than the last variable in data object, this happens because of
        # DST time change
        if len(binner) > 1 and binner[-1] < last:
            extra_date_range = pd.date_range(binner[-1], last + self.freq,
                                             freq=self.freq, tz=tz,
                                             name=ax.name)
            binner = labels = binner.append(extra_date_range[1:])

        # a little hack
        trimmed = False
        if (len(binner) > 2 and binner[-2] == last and
                self.closed == 'right'):

            binner = binner[:-1]
            trimmed = True

        ax_values = ax.asi8
        binner, bin_edges = self._adjust_bin_edges(binner, ax_values)

        # general version, knowing nothing about relative frequencies
        bins = lib.generate_bins_dt64(
            ax_values, bin_edges, self.closed, hasnans=ax.hasnans)

        if self.closed == 'right':
            labels = binner
            if self.label == 'right':
                labels = labels[1:]
            elif not trimmed:
                labels = labels[:-1]
        else:
            if self.label == 'right':
                labels = labels[1:]
            elif not trimmed:
                labels = labels[:-1]

        if ax.hasnans:
            binner = binner.insert(0, tslib.NaT)
            labels = labels.insert(0, tslib.NaT)

        # if we end up with more labels than bins
        # adjust the labels
        # GH4076
        if len(bins) < len(labels):
            labels = labels[:len(bins)]

        return binner, bins, labels

    def _adjust_bin_edges(self, binner, ax_values):
        # Some hacks for > daily data, see #1471, #1458, #1483

        bin_edges = binner.asi8

        if self.freq != 'D' and is_superperiod(self.freq, 'D'):
            day_nanos = delta_to_nanoseconds(timedelta(1))
            if self.closed == 'right':
                bin_edges = bin_edges + day_nanos - 1

            # intraday values on last day
            if bin_edges[-2] > ax_values.max():
                bin_edges = bin_edges[:-1]
                binner = binner[:-1]

        return binner, bin_edges

    def _get_time_delta_bins(self, ax):
        if not isinstance(ax, TimedeltaIndex):
            raise TypeError('axis must be a TimedeltaIndex, but got '
                            'an instance of %r' % type(ax).__name__)

        if not len(ax):
            binner = labels = TimedeltaIndex(
                data=[], freq=self.freq, name=ax.name)
            return binner, [], labels

        start = ax[0]
        end = ax[-1]
        labels = binner = TimedeltaIndex(start=start,
                                         end=end,
                                         freq=self.freq,
                                         name=ax.name)

        end_stamps = labels + 1
        bins = ax.searchsorted(end_stamps, side='left')

        # Addresses GH #10530
        if self.base > 0:
            labels += type(self.freq)(self.base)

        return binner, bins, labels

    def _get_time_period_bins(self, ax):
        if not isinstance(ax, DatetimeIndex):
            raise TypeError('axis must be a DatetimeIndex, but got '
                            'an instance of %r' % type(ax).__name__)

        if not len(ax):
            binner = labels = PeriodIndex(
                data=[], freq=self.freq, name=ax.name)
            return binner, [], labels

        labels = binner = PeriodIndex(start=ax[0],
                                      end=ax[-1],
                                      freq=self.freq,
                                      name=ax.name)

        end_stamps = (labels + 1).asfreq(self.freq, 's').to_timestamp()
        if ax.tzinfo:
            end_stamps = end_stamps.tz_localize(ax.tzinfo)
        bins = ax.searchsorted(end_stamps, side='left')

        return binner, bins, labels

    def _get_period_bins(self, ax):
        if not isinstance(ax, PeriodIndex):
            raise TypeError('axis must be a PeriodIndex, but got '
                            'an instance of %r' % type(ax).__name__)

        memb = ax.asfreq(self.freq, how=self.convention)

        # NaT handling as in pandas._lib.lib.generate_bins_dt64()
        nat_count = 0
        if memb.hasnans:
            nat_count = np.sum(memb._isnan)
            memb = memb[~memb._isnan]

        # if index contains no valid (non-NaT) values, return empty index
        if not len(memb):
            binner = labels = PeriodIndex(
                data=[], freq=self.freq, name=ax.name)
            return binner, [], labels

        start = ax.min().asfreq(self.freq, how=self.convention)
        end = ax.max().asfreq(self.freq, how='end')

        labels = binner = PeriodIndex(start=start, end=end,
                                      freq=self.freq, name=ax.name)

        i8 = memb.asi8
        freq_mult = self.freq.n

        # when upsampling to subperiods, we need to generate enough bins
        expected_bins_count = len(binner) * freq_mult
        i8_extend = expected_bins_count - (i8[-1] - i8[0])
        rng = np.arange(i8[0], i8[-1] + i8_extend, freq_mult)
        rng += freq_mult
        bins = memb.searchsorted(rng, side='left')

        if nat_count > 0:
            # NaT handling as in pandas._lib.lib.generate_bins_dt64()
            # shift bins by the number of NaT
            bins += nat_count
            bins = np.insert(bins, 0, nat_count)
            binner = binner.insert(0, tslib.NaT)
            labels = labels.insert(0, tslib.NaT)

        return binner, bins, labels


def _take_new_index(obj, indexer, new_index, axis=0):
    from pandas.core.api import Series, DataFrame

    if isinstance(obj, Series):
        new_values = algos.take_1d(obj.values, indexer)
        return Series(new_values, index=new_index, name=obj.name)
    elif isinstance(obj, DataFrame):
        if axis == 1:
            raise NotImplementedError("axis 1 is not supported")
        return DataFrame(obj._data.reindex_indexer(
            new_axis=new_index, indexer=indexer, axis=1))
    else:
        raise ValueError("'obj' should be either a Series or a DataFrame")


def _get_range_edges(first, last, offset, closed='left', base=0):
    if isinstance(offset, compat.string_types):
        offset = to_offset(offset)

    if isinstance(offset, Tick):
        is_day = isinstance(offset, Day)
        day_nanos = delta_to_nanoseconds(timedelta(1))

        # #1165
        if (is_day and day_nanos % offset.nanos == 0) or not is_day:
            return _adjust_dates_anchored(first, last, offset,
                                          closed=closed, base=base)

    if not isinstance(offset, Tick):  # and first.time() != last.time():
        # hack!
        first = first.normalize()
        last = last.normalize()

    if closed == 'left':
        first = Timestamp(offset.rollback(first))
    else:
        first = Timestamp(first - offset)

    last = Timestamp(last + offset)

    return first, last


def _adjust_dates_anchored(first, last, offset, closed='right', base=0):
    # First and last offsets should be calculated from the start day to fix an
    # error cause by resampling across multiple days when a one day period is
    # not a multiple of the frequency.
    #
    # See https://github.com/pandas-dev/pandas/issues/8683

    # 14682 - Since we need to drop the TZ information to perform
    # the adjustment in the presence of a DST change,
    # save TZ Info and the DST state of the first and last parameters
    # so that we can accurately rebuild them at the end.
    first_tzinfo = first.tzinfo
    last_tzinfo = last.tzinfo
    first_dst = bool(first.dst())
    last_dst = bool(last.dst())

    first = first.tz_localize(None)
    last = last.tz_localize(None)

    start_day_nanos = first.normalize().value

    base_nanos = (base % offset.n) * offset.nanos // offset.n
    start_day_nanos += base_nanos

    foffset = (first.value - start_day_nanos) % offset.nanos
    loffset = (last.value - start_day_nanos) % offset.nanos

    if closed == 'right':
        if foffset > 0:
            # roll back
            fresult = first.value - foffset
        else:
            fresult = first.value - offset.nanos

        if loffset > 0:
            # roll forward
            lresult = last.value + (offset.nanos - loffset)
        else:
            # already the end of the road
            lresult = last.value
    else:  # closed == 'left'
        if foffset > 0:
            fresult = first.value - foffset
        else:
            # start of the road
            fresult = first.value

        if loffset > 0:
            # roll forward
            lresult = last.value + (offset.nanos - loffset)
        else:
            lresult = last.value + offset.nanos

    return (Timestamp(fresult).tz_localize(first_tzinfo, ambiguous=first_dst),
            Timestamp(lresult).tz_localize(last_tzinfo, ambiguous=last_dst))


def asfreq(obj, freq, method=None, how=None, normalize=False, fill_value=None):
    """
    Utility frequency conversion method for Series/DataFrame
    """
    if isinstance(obj.index, PeriodIndex):
        if method is not None:
            raise NotImplementedError("'method' argument is not supported")

        if how is None:
            how = 'E'

        new_obj = obj.copy()
        new_obj.index = obj.index.asfreq(freq, how=how)

    elif len(obj.index) == 0:
        new_obj = obj.copy()
        new_obj.index = obj.index._shallow_copy(freq=to_offset(freq))

    else:
        dti = date_range(obj.index[0], obj.index[-1], freq=freq)
        dti.name = obj.index.name
        new_obj = obj.reindex(dti, method=method, fill_value=fill_value)
        if normalize:
            new_obj.index = new_obj.index.normalize()

    return new_obj
