""" implement the TimedeltaIndex """

from datetime import timedelta

import numpy as np
from pandas.core.dtypes.common import (
    _TD_DTYPE,
    is_integer,
    is_float,
    is_bool_dtype,
    is_list_like,
    is_scalar,
    is_timedelta64_dtype,
    is_timedelta64_ns_dtype,
    pandas_dtype,
    _ensure_int64)
from pandas.core.dtypes.missing import isna
from pandas.core.dtypes.generic import ABCSeries

from pandas.core.indexes.base import Index
from pandas.core.indexes.numeric import Int64Index
import pandas.compat as compat
from pandas.compat import u
from pandas.tseries.frequencies import to_offset
from pandas.core.algorithms import checked_add_with_arr
from pandas.core.base import _shared_docs
from pandas.core.indexes.base import _index_shared_docs
import pandas.core.common as com
import pandas.core.dtypes.concat as _concat
from pandas.util._decorators import Appender, Substitution, deprecate_kwarg
from pandas.core.indexes.datetimelike import (
    TimelikeOps, DatetimeIndexOpsMixin)
from pandas.core.tools.timedeltas import (
    to_timedelta, _coerce_scalar_to_timedelta_type)
from pandas.tseries.offsets import Tick, DateOffset
from pandas._libs import (lib, index as libindex, tslib as libts,
                          join as libjoin, Timedelta, NaT, iNaT)
from pandas._libs.tslibs.timedeltas import array_to_timedelta64
from pandas._libs.tslibs.fields import get_timedelta_field


def _field_accessor(name, alias, docstring=None):
    def f(self):
        values = self.asi8
        result = get_timedelta_field(values, alias)
        if self.hasnans:
            result = self._maybe_mask_results(result, convert='float64')

        return Index(result, name=self.name)

    f.__name__ = name
    f.__doc__ = docstring
    return property(f)


def _td_index_cmp(opname, cls):
    """
    Wrap comparison operations to convert timedelta-like to timedelta64
    """
    nat_result = True if opname == '__ne__' else False

    def wrapper(self, other):
        msg = "cannot compare a {cls} with type {typ}"
        func = getattr(super(TimedeltaIndex, self), opname)
        if _is_convertible_to_td(other) or other is NaT:
            try:
                other = _to_m8(other)
            except ValueError:
                # failed to parse as timedelta
                raise TypeError(msg.format(cls=type(self).__name__,
                                           typ=type(other).__name__))
            result = func(other)
            if isna(other):
                result.fill(nat_result)

        elif not is_list_like(other):
            raise TypeError(msg.format(cls=type(self).__name__,
                                       typ=type(other).__name__))
        else:
            other = TimedeltaIndex(other).values
            result = func(other)
            result = com._values_from_object(result)

            o_mask = np.array(isna(other))
            if o_mask.any():
                result[o_mask] = nat_result

        if self.hasnans:
            result[self._isnan] = nat_result

        # support of bool dtype indexers
        if is_bool_dtype(result):
            return result
        return Index(result)

    return compat.set_function_name(wrapper, opname, cls)


class TimedeltaIndex(DatetimeIndexOpsMixin, TimelikeOps, Int64Index):
    """
    Immutable ndarray of timedelta64 data, represented internally as int64, and
    which can be boxed to timedelta objects

    Parameters
    ----------
    data  : array-like (1-dimensional), optional
        Optional timedelta-like data to construct index with
    unit: unit of the arg (D,h,m,s,ms,us,ns) denote the unit, optional
        which is an integer/float number
    freq: a frequency for the index, optional
    copy  : bool
        Make a copy of input ndarray
    start : starting value, timedelta-like, optional
        If data is None, start is used as the start point in generating regular
        timedelta data.
    periods  : int, optional, > 0
        Number of periods to generate, if generating index. Takes precedence
        over end argument
    end   : end time, timedelta-like, optional
        If periods is none, generated index will extend to first conforming
        time on or just past end argument
    closed : string or None, default None
        Make the interval closed with respect to the given frequency to
        the 'left', 'right', or both sides (None)
    name : object
        Name to be stored in the index

    Notes
    -----

    To learn more about the frequency strings, please see `this link
    <http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases>`__.

    See Also
    ---------
    Index : The base pandas Index type
    Timedelta : Represents a duration between two dates or times.
    DatetimeIndex : Index of datetime64 data
    PeriodIndex : Index of Period data

    Attributes
    ----------
    days
    seconds
    microseconds
    nanoseconds
    components
    inferred_freq

    Methods
    -------
    to_pytimedelta
    to_series
    round
    floor
    ceil
    to_frame
    """

    _typ = 'timedeltaindex'
    _join_precedence = 10

    def _join_i8_wrapper(joinf, **kwargs):
        return DatetimeIndexOpsMixin._join_i8_wrapper(
            joinf, dtype='m8[ns]', **kwargs)

    _inner_indexer = _join_i8_wrapper(libjoin.inner_join_indexer_int64)
    _outer_indexer = _join_i8_wrapper(libjoin.outer_join_indexer_int64)
    _left_indexer = _join_i8_wrapper(libjoin.left_join_indexer_int64)
    _left_indexer_unique = _join_i8_wrapper(
        libjoin.left_join_indexer_unique_int64, with_indexers=False)

    # define my properties & methods for delegation
    _other_ops = []
    _bool_ops = []
    _object_ops = ['freq']
    _field_ops = ['days', 'seconds', 'microseconds', 'nanoseconds']
    _datetimelike_ops = _field_ops + _object_ops + _bool_ops
    _datetimelike_methods = ["to_pytimedelta", "total_seconds",
                             "round", "floor", "ceil"]

    @classmethod
    def _add_comparison_methods(cls):
        """ add in comparison methods """
        cls.__eq__ = _td_index_cmp('__eq__', cls)
        cls.__ne__ = _td_index_cmp('__ne__', cls)
        cls.__lt__ = _td_index_cmp('__lt__', cls)
        cls.__gt__ = _td_index_cmp('__gt__', cls)
        cls.__le__ = _td_index_cmp('__le__', cls)
        cls.__ge__ = _td_index_cmp('__ge__', cls)

    _engine_type = libindex.TimedeltaEngine

    _comparables = ['name', 'freq']
    _attributes = ['name', 'freq']
    _is_numeric_dtype = True
    _infer_as_myclass = True

    _freq = None

    def __new__(cls, data=None, unit=None, freq=None, start=None, end=None,
                periods=None, closed=None, dtype=None, copy=False,
                name=None, verify_integrity=True):

        if isinstance(data, TimedeltaIndex) and freq is None and name is None:
            if copy:
                return data.copy()
            else:
                return data._shallow_copy()

        freq_infer = False
        if not isinstance(freq, DateOffset):

            # if a passed freq is None, don't infer automatically
            if freq != 'infer':
                freq = to_offset(freq)
            else:
                freq_infer = True
                freq = None

        if periods is not None:
            if is_float(periods):
                periods = int(periods)
            elif not is_integer(periods):
                msg = 'periods must be a number, got {periods}'
                raise TypeError(msg.format(periods=periods))

        if data is None:
            if freq is None and com._any_none(periods, start, end):
                msg = 'Must provide freq argument if no data is supplied'
                raise ValueError(msg)
            else:
                return cls._generate(start, end, periods, name, freq,
                                     closed=closed)

        if unit is not None:
            data = to_timedelta(data, unit=unit, box=False)

        if not isinstance(data, (np.ndarray, Index, ABCSeries)):
            if is_scalar(data):
                raise ValueError('TimedeltaIndex() must be called with a '
                                 'collection of some kind, %s was passed'
                                 % repr(data))

        # convert if not already
        if getattr(data, 'dtype', None) != _TD_DTYPE:
            data = to_timedelta(data, unit=unit, box=False)
        elif copy:
            data = np.array(data, copy=True)

        # check that we are matching freqs
        if verify_integrity and len(data) > 0:
            if freq is not None and not freq_infer:
                index = cls._simple_new(data, name=name)
                cls._validate_frequency(index, freq)
                index.freq = freq
                return index

        if freq_infer:
            index = cls._simple_new(data, name=name)
            inferred = index.inferred_freq
            if inferred:
                index.freq = to_offset(inferred)
            return index

        return cls._simple_new(data, name=name, freq=freq)

    @classmethod
    def _generate(cls, start, end, periods, name, freq, closed=None):
        if com._count_not_none(start, end, periods, freq) != 3:
            raise ValueError('Of the four parameters: start, end, periods, '
                             'and freq, exactly three must be specified')

        if start is not None:
            start = Timedelta(start)

        if end is not None:
            end = Timedelta(end)

        left_closed = False
        right_closed = False

        if start is None and end is None:
            if closed is not None:
                raise ValueError("Closed has to be None if not both of start"
                                 "and end are defined")

        if closed is None:
            left_closed = True
            right_closed = True
        elif closed == "left":
            left_closed = True
        elif closed == "right":
            right_closed = True
        else:
            raise ValueError("Closed has to be either 'left', 'right' or None")

        if freq is not None:
            index = _generate_regular_range(start, end, periods, freq)
            index = cls._simple_new(index, name=name, freq=freq)
        else:
            index = to_timedelta(np.linspace(start.value, end.value, periods))

        if not left_closed:
            index = index[1:]
        if not right_closed:
            index = index[:-1]

        return index

    @property
    def _box_func(self):
        return lambda x: Timedelta(x, unit='ns')

    @classmethod
    def _simple_new(cls, values, name=None, freq=None, **kwargs):
        values = np.array(values, copy=False)
        if values.dtype == np.object_:
            values = array_to_timedelta64(values)
        if values.dtype != _TD_DTYPE:
            values = _ensure_int64(values).view(_TD_DTYPE)

        result = object.__new__(cls)
        result._data = values
        result.name = name
        result._freq = freq
        result._reset_identity()
        return result

    @property
    def _formatter_func(self):
        from pandas.io.formats.format import _get_format_timedelta64
        return _get_format_timedelta64(self, box=True)

    def __setstate__(self, state):
        """Necessary for making this object picklable"""
        if isinstance(state, dict):
            super(TimedeltaIndex, self).__setstate__(state)
        else:
            raise Exception("invalid pickle state")
    _unpickle_compat = __setstate__

    def _maybe_update_attributes(self, attrs):
        """ Update Index attributes (e.g. freq) depending on op """
        freq = attrs.get('freq', None)
        if freq is not None:
            # no need to infer if freq is None
            attrs['freq'] = 'infer'
        return attrs

    def _add_offset(self, other):
        assert not isinstance(other, Tick)
        raise TypeError("cannot add the type {typ} to a {cls}"
                        .format(typ=type(other).__name__,
                                cls=type(self).__name__))

    def _add_delta(self, delta):
        """
        Add a timedelta-like, Tick, or TimedeltaIndex-like object
        to self.

        Parameters
        ----------
        delta : {timedelta, np.timedelta64, Tick, TimedeltaIndex}

        Returns
        -------
        result : TimedeltaIndex

        Notes
        -----
        The result's name is set outside of _add_delta by the calling
        method (__add__ or __sub__)
        """
        if isinstance(delta, (Tick, timedelta, np.timedelta64)):
            new_values = self._add_delta_td(delta)
        elif isinstance(delta, TimedeltaIndex):
            new_values = self._add_delta_tdi(delta)
        elif is_timedelta64_dtype(delta):
            # ndarray[timedelta64] --> wrap in TimedeltaIndex
            delta = TimedeltaIndex(delta)
            new_values = self._add_delta_tdi(delta)
        else:
            raise TypeError("cannot add the type {0} to a TimedeltaIndex"
                            .format(type(delta)))

        return TimedeltaIndex(new_values, freq='infer')

    def _evaluate_with_timedelta_like(self, other, op):
        if isinstance(other, ABCSeries):
            # GH#19042
            return NotImplemented

        opstr = '__{opname}__'.format(opname=op.__name__).replace('__r', '__')
        # allow division by a timedelta
        if opstr in ['__div__', '__truediv__', '__floordiv__']:
            if _is_convertible_to_td(other):
                other = Timedelta(other)
                if isna(other):
                    raise NotImplementedError(
                        "division by pd.NaT not implemented")

                i8 = self.asi8
                left, right = i8, other.value

                if opstr in ['__floordiv__']:
                    result = op(left, right)
                else:
                    result = op(left, np.float64(right))
                result = self._maybe_mask_results(result, convert='float64')
                return Index(result, name=self.name, copy=False)

        return NotImplemented

    def _add_datelike(self, other):
        # adding a timedeltaindex to a datetimelike
        from pandas import Timestamp, DatetimeIndex
        if isinstance(other, (DatetimeIndex, np.ndarray)):
            # if other is an ndarray, we assume it is datetime64-dtype
            # defer to implementation in DatetimeIndex
            other = DatetimeIndex(other)
            return other + self
        else:
            assert other is not NaT
            other = Timestamp(other)
            i8 = self.asi8
            result = checked_add_with_arr(i8, other.value,
                                          arr_mask=self._isnan)
            result = self._maybe_mask_results(result, fill_value=iNaT)
            return DatetimeIndex(result)

    def _sub_datelike(self, other):
        assert other is not NaT
        raise TypeError("cannot subtract a datelike from a {cls}"
                        .format(cls=type(self).__name__))

    def _addsub_offset_array(self, other, op):
        # Add or subtract Array-like of DateOffset objects
        try:
            # TimedeltaIndex can only operate with a subset of DateOffset
            # subclasses.  Incompatible classes will raise AttributeError,
            # which we re-raise as TypeError
            return DatetimeIndexOpsMixin._addsub_offset_array(self, other, op)
        except AttributeError:
            raise TypeError("Cannot add/subtract non-tick DateOffset to {cls}"
                            .format(cls=type(self).__name__))

    def _format_native_types(self, na_rep=u('NaT'),
                             date_format=None, **kwargs):
        from pandas.io.formats.format import Timedelta64Formatter
        return Timedelta64Formatter(values=self,
                                    nat_rep=na_rep,
                                    justify='all').get_result()

    days = _field_accessor("days", "days",
                           " Number of days for each element. ")
    seconds = _field_accessor("seconds", "seconds",
                              " Number of seconds (>= 0 and less than 1 day) "
                              "for each element. ")
    microseconds = _field_accessor("microseconds", "microseconds",
                                   "\nNumber of microseconds (>= 0 and less "
                                   "than 1 second) for each\nelement. ")
    nanoseconds = _field_accessor("nanoseconds", "nanoseconds",
                                  "\nNumber of nanoseconds (>= 0 and less "
                                  "than 1 microsecond) for each\nelement.\n")

    @property
    def components(self):
        """
        Return a dataframe of the components (days, hours, minutes,
        seconds, milliseconds, microseconds, nanoseconds) of the Timedeltas.

        Returns
        -------
        a DataFrame
        """
        from pandas import DataFrame

        columns = ['days', 'hours', 'minutes', 'seconds',
                   'milliseconds', 'microseconds', 'nanoseconds']
        hasnans = self.hasnans
        if hasnans:
            def f(x):
                if isna(x):
                    return [np.nan] * len(columns)
                return x.components
        else:
            def f(x):
                return x.components

        result = DataFrame([f(x) for x in self])
        result.columns = columns
        if not hasnans:
            result = result.astype('int64')
        return result

    def total_seconds(self):
        """
        Return total duration of each element expressed in seconds.

        This method is available directly on TimedeltaIndex and on Series
        containing timedelta values under the ``.dt`` namespace.

        Returns
        -------
        seconds : Float64Index or Series
            When the calling object is a TimedeltaIndex, the return type is a
            Float64Index. When the calling object is a Series, the return type
            is Series of type `float64` whose index is the same as the
            original.

        See Also
        --------
        datetime.timedelta.total_seconds : Standard library version
            of this method.
        TimedeltaIndex.components : Return a DataFrame with components of
            each Timedelta.

        Examples
        --------
        **Series**

        >>> s = pd.Series(pd.to_timedelta(np.arange(5), unit='d'))
        >>> s
        0   0 days
        1   1 days
        2   2 days
        3   3 days
        4   4 days
        dtype: timedelta64[ns]

        >>> s.dt.total_seconds()
        0         0.0
        1     86400.0
        2    172800.0
        3    259200.0
        4    345600.0
        dtype: float64

        **TimedeltaIndex**

        >>> idx = pd.to_timedelta(np.arange(5), unit='d')
        >>> idx
        TimedeltaIndex(['0 days', '1 days', '2 days', '3 days', '4 days'],
                       dtype='timedelta64[ns]', freq=None)

        >>> idx.total_seconds()
        Float64Index([0.0, 86400.0, 172800.0, 259200.00000000003, 345600.0],
                     dtype='float64')
        """
        return Index(self._maybe_mask_results(1e-9 * self.asi8),
                     name=self.name)

    def to_pytimedelta(self):
        """
        Return TimedeltaIndex as object ndarray of datetime.timedelta objects

        Returns
        -------
        datetimes : ndarray
        """
        return libts.ints_to_pytimedelta(self.asi8)

    @Appender(_index_shared_docs['astype'])
    def astype(self, dtype, copy=True):
        dtype = pandas_dtype(dtype)
        if is_timedelta64_dtype(dtype) and not is_timedelta64_ns_dtype(dtype):
            # return an index (essentially this is division)
            result = self.values.astype(dtype, copy=copy)
            if self.hasnans:
                values = self._maybe_mask_results(result, convert='float64')
                return Index(values, name=self.name)
            return Index(result.astype('i8'), name=self.name)
        return super(TimedeltaIndex, self).astype(dtype, copy=copy)

    def union(self, other):
        """
        Specialized union for TimedeltaIndex objects. If combine
        overlapping ranges with the same DateOffset, will be much
        faster than Index.union

        Parameters
        ----------
        other : TimedeltaIndex or array-like

        Returns
        -------
        y : Index or TimedeltaIndex
        """
        self._assert_can_do_setop(other)
        if not isinstance(other, TimedeltaIndex):
            try:
                other = TimedeltaIndex(other)
            except (TypeError, ValueError):
                pass
        this, other = self, other

        if this._can_fast_union(other):
            return this._fast_union(other)
        else:
            result = Index.union(this, other)
            if isinstance(result, TimedeltaIndex):
                if result.freq is None:
                    result.freq = to_offset(result.inferred_freq)
            return result

    def join(self, other, how='left', level=None, return_indexers=False,
             sort=False):
        """
        See Index.join
        """
        if _is_convertible_to_index(other):
            try:
                other = TimedeltaIndex(other)
            except (TypeError, ValueError):
                pass

        return Index.join(self, other, how=how, level=level,
                          return_indexers=return_indexers,
                          sort=sort)

    def _wrap_joined_index(self, joined, other):
        name = self.name if self.name == other.name else None
        if (isinstance(other, TimedeltaIndex) and self.freq == other.freq and
                self._can_fast_union(other)):
            joined = self._shallow_copy(joined, name=name)
            return joined
        else:
            return self._simple_new(joined, name)

    def _can_fast_union(self, other):
        if not isinstance(other, TimedeltaIndex):
            return False

        freq = self.freq

        if freq is None or freq != other.freq:
            return False

        if not self.is_monotonic or not other.is_monotonic:
            return False

        if len(self) == 0 or len(other) == 0:
            return True

        # to make our life easier, "sort" the two ranges
        if self[0] <= other[0]:
            left, right = self, other
        else:
            left, right = other, self

        right_start = right[0]
        left_end = left[-1]

        # Only need to "adjoin", not overlap
        return (right_start == left_end + freq) or right_start in left

    def _fast_union(self, other):
        if len(other) == 0:
            return self.view(type(self))

        if len(self) == 0:
            return other.view(type(self))

        # to make our life easier, "sort" the two ranges
        if self[0] <= other[0]:
            left, right = self, other
        else:
            left, right = other, self

        left_end = left[-1]
        right_end = right[-1]

        # concatenate
        if left_end < right_end:
            loc = right.searchsorted(left_end, side='right')
            right_chunk = right.values[loc:]
            dates = _concat._concat_compat((left.values, right_chunk))
            return self._shallow_copy(dates)
        else:
            return left

    def _wrap_union_result(self, other, result):
        name = self.name if self.name == other.name else None
        return self._simple_new(result, name=name, freq=None)

    def intersection(self, other):
        """
        Specialized intersection for TimedeltaIndex objects. May be much faster
        than Index.intersection

        Parameters
        ----------
        other : TimedeltaIndex or array-like

        Returns
        -------
        y : Index or TimedeltaIndex
        """
        self._assert_can_do_setop(other)
        if not isinstance(other, TimedeltaIndex):
            try:
                other = TimedeltaIndex(other)
            except (TypeError, ValueError):
                pass
            result = Index.intersection(self, other)
            return result

        if len(self) == 0:
            return self
        if len(other) == 0:
            return other
        # to make our life easier, "sort" the two ranges
        if self[0] <= other[0]:
            left, right = self, other
        else:
            left, right = other, self

        end = min(left[-1], right[-1])
        start = right[0]

        if end < start:
            return type(self)(data=[])
        else:
            lslice = slice(*left.slice_locs(start, end))
            left_chunk = left.values[lslice]
            return self._shallow_copy(left_chunk)

    def _maybe_promote(self, other):
        if other.inferred_type == 'timedelta':
            other = TimedeltaIndex(other)
        return self, other

    def get_value(self, series, key):
        """
        Fast lookup of value from 1-dimensional ndarray. Only use this if you
        know what you're doing
        """

        if _is_convertible_to_td(key):
            key = Timedelta(key)
            return self.get_value_maybe_box(series, key)

        try:
            return com._maybe_box(self, Index.get_value(self, series, key),
                                  series, key)
        except KeyError:
            try:
                loc = self._get_string_slice(key)
                return series[loc]
            except (TypeError, ValueError, KeyError):
                pass

            try:
                return self.get_value_maybe_box(series, key)
            except (TypeError, ValueError, KeyError):
                raise KeyError(key)

    def get_value_maybe_box(self, series, key):
        if not isinstance(key, Timedelta):
            key = Timedelta(key)
        values = self._engine.get_value(com._values_from_object(series), key)
        return com._maybe_box(self, values, series, key)

    def get_loc(self, key, method=None, tolerance=None):
        """
        Get integer location for requested label

        Returns
        -------
        loc : int
        """
        if is_list_like(key):
            raise TypeError

        if isna(key):
            key = NaT

        if tolerance is not None:
            # try converting tolerance now, so errors don't get swallowed by
            # the try/except clauses below
            tolerance = self._convert_tolerance(tolerance, np.asarray(key))

        if _is_convertible_to_td(key):
            key = Timedelta(key)
            return Index.get_loc(self, key, method, tolerance)

        try:
            return Index.get_loc(self, key, method, tolerance)
        except (KeyError, ValueError, TypeError):
            try:
                return self._get_string_slice(key)
            except (TypeError, KeyError, ValueError):
                pass

            try:
                stamp = Timedelta(key)
                return Index.get_loc(self, stamp, method, tolerance)
            except (KeyError, ValueError):
                raise KeyError(key)

    def _maybe_cast_slice_bound(self, label, side, kind):
        """
        If label is a string, cast it to timedelta according to resolution.


        Parameters
        ----------
        label : object
        side : {'left', 'right'}
        kind : {'ix', 'loc', 'getitem'}

        Returns
        -------
        label :  object

        """
        assert kind in ['ix', 'loc', 'getitem', None]

        if isinstance(label, compat.string_types):
            parsed = _coerce_scalar_to_timedelta_type(label, box=True)
            lbound = parsed.round(parsed.resolution)
            if side == 'left':
                return lbound
            else:
                return (lbound + to_offset(parsed.resolution) -
                        Timedelta(1, 'ns'))
        elif ((is_integer(label) or is_float(label)) and
              not is_timedelta64_dtype(label)):
            self._invalid_indexer('slice', label)

        return label

    def _get_string_slice(self, key, use_lhs=True, use_rhs=True):
        freq = getattr(self, 'freqstr',
                       getattr(self, 'inferred_freq', None))
        if is_integer(key) or is_float(key) or key is NaT:
            self._invalid_indexer('slice', key)
        loc = self._partial_td_slice(key, freq, use_lhs=use_lhs,
                                     use_rhs=use_rhs)
        return loc

    def _partial_td_slice(self, key, freq, use_lhs=True, use_rhs=True):

        # given a key, try to figure out a location for a partial slice
        if not isinstance(key, compat.string_types):
            return key

        raise NotImplementedError

        # TODO(wesm): dead code
        # parsed = _coerce_scalar_to_timedelta_type(key, box=True)

        # is_monotonic = self.is_monotonic

        # # figure out the resolution of the passed td
        # # and round to it

        # # t1 = parsed.round(reso)

        # t2 = t1 + to_offset(parsed.resolution) - Timedelta(1, 'ns')

        # stamps = self.asi8

        # if is_monotonic:

        #     # we are out of range
        #     if (len(stamps) and ((use_lhs and t1.value < stamps[0] and
        #                           t2.value < stamps[0]) or
        #                          ((use_rhs and t1.value > stamps[-1] and
        #                            t2.value > stamps[-1])))):
        #         raise KeyError

        #     # a monotonic (sorted) series can be sliced
        #     left = (stamps.searchsorted(t1.value, side='left')
        #             if use_lhs else None)
        #     right = (stamps.searchsorted(t2.value, side='right')
        #              if use_rhs else None)

        #     return slice(left, right)

        # lhs_mask = (stamps >= t1.value) if use_lhs else True
        # rhs_mask = (stamps <= t2.value) if use_rhs else True

        # # try to find a the dates
        # return (lhs_mask & rhs_mask).nonzero()[0]

    @Substitution(klass='TimedeltaIndex')
    @Appender(_shared_docs['searchsorted'])
    @deprecate_kwarg(old_arg_name='key', new_arg_name='value')
    def searchsorted(self, value, side='left', sorter=None):
        if isinstance(value, (np.ndarray, Index)):
            value = np.array(value, dtype=_TD_DTYPE, copy=False)
        else:
            value = _to_m8(value)

        return self.values.searchsorted(value, side=side, sorter=sorter)

    def is_type_compatible(self, typ):
        return typ == self.inferred_type or typ == 'timedelta'

    @property
    def inferred_type(self):
        return 'timedelta64'

    @property
    def dtype(self):
        return _TD_DTYPE

    @property
    def is_all_dates(self):
        return True

    def insert(self, loc, item):
        """
        Make new Index inserting new item at location

        Parameters
        ----------
        loc : int
        item : object
            if not either a Python datetime or a numpy integer-like, returned
            Index dtype will be object rather than datetime.

        Returns
        -------
        new_index : Index
        """
        # try to convert if possible
        if _is_convertible_to_td(item):
            try:
                item = Timedelta(item)
            except Exception:
                pass
        elif is_scalar(item) and isna(item):
            # GH 18295
            item = self._na_value

        freq = None
        if isinstance(item, Timedelta) or (is_scalar(item) and isna(item)):

            # check freq can be preserved on edge cases
            if self.freq is not None:
                if ((loc == 0 or loc == -len(self)) and
                        item + self.freq == self[0]):
                    freq = self.freq
                elif (loc == len(self)) and item - self.freq == self[-1]:
                    freq = self.freq
            item = _to_m8(item)

        try:
            new_tds = np.concatenate((self[:loc].asi8, [item.view(np.int64)],
                                      self[loc:].asi8))
            return TimedeltaIndex(new_tds, name=self.name, freq=freq)

        except (AttributeError, TypeError):

            # fall back to object index
            if isinstance(item, compat.string_types):
                return self.astype(object).insert(loc, item)
            raise TypeError(
                "cannot insert TimedeltaIndex with incompatible label")

    def delete(self, loc):
        """
        Make a new TimedeltaIndex with passed location(s) deleted.

        Parameters
        ----------
        loc: int, slice or array of ints
            Indicate which sub-arrays to remove.

        Returns
        -------
        new_index : TimedeltaIndex
        """
        new_tds = np.delete(self.asi8, loc)

        freq = 'infer'
        if is_integer(loc):
            if loc in (0, -len(self), -1, len(self) - 1):
                freq = self.freq
        else:
            if is_list_like(loc):
                loc = lib.maybe_indices_to_slice(
                    _ensure_int64(np.array(loc)), len(self))
            if isinstance(loc, slice) and loc.step in (1, None):
                if (loc.start in (0, None) or loc.stop in (len(self), None)):
                    freq = self.freq

        return TimedeltaIndex(new_tds, name=self.name, freq=freq)


TimedeltaIndex._add_comparison_methods()
TimedeltaIndex._add_numeric_methods()
TimedeltaIndex._add_logical_methods_disabled()
TimedeltaIndex._add_datetimelike_methods()


def _is_convertible_to_index(other):
    """
    return a boolean whether I can attempt conversion to a TimedeltaIndex
    """
    if isinstance(other, TimedeltaIndex):
        return True
    elif (len(other) > 0 and
          other.inferred_type not in ('floating', 'mixed-integer', 'integer',
                                      'mixed-integer-float', 'mixed')):
        return True
    return False


def _is_convertible_to_td(key):
    return isinstance(key, (Tick, timedelta,
                            np.timedelta64, compat.string_types))


def _to_m8(key):
    """
    Timedelta-like => dt64
    """
    if not isinstance(key, Timedelta):
        # this also converts strings
        key = Timedelta(key)

    # return an type that can be compared
    return np.int64(key.value).view(_TD_DTYPE)


def _generate_regular_range(start, end, periods, offset):
    stride = offset.nanos
    if periods is None:
        b = Timedelta(start).value
        e = Timedelta(end).value
        e += stride - e % stride
    elif start is not None:
        b = Timedelta(start).value
        e = b + periods * stride
    elif end is not None:
        e = Timedelta(end).value + stride
        b = e - periods * stride
    else:
        raise ValueError("at least 'start' or 'end' should be specified "
                         "if a 'period' is given.")

    data = np.arange(b, e, stride, dtype=np.int64)
    data = TimedeltaIndex._simple_new(data, None)

    return data


def timedelta_range(start=None, end=None, periods=None, freq=None,
                    name=None, closed=None):
    """
    Return a fixed frequency TimedeltaIndex, with day as the default
    frequency

    Parameters
    ----------
    start : string or timedelta-like, default None
        Left bound for generating timedeltas
    end : string or timedelta-like, default None
        Right bound for generating timedeltas
    periods : integer, default None
        Number of periods to generate
    freq : string or DateOffset, default 'D' (calendar daily)
        Frequency strings can have multiples, e.g. '5H'
    name : string, default None
        Name of the resulting TimedeltaIndex
    closed : string, default None
        Make the interval closed with respect to the given frequency to
        the 'left', 'right', or both sides (None)

    Returns
    -------
    rng : TimedeltaIndex

    Notes
    -----
    Of the four parameters ``start``, ``end``, ``periods``, and ``freq``,
    exactly three must be specified. If ``freq`` is omitted, the resulting
    ``TimedeltaIndex`` will have ``periods`` linearly spaced elements between
    ``start`` and ``end`` (closed on both sides).

    To learn more about the frequency strings, please see `this link
    <http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases>`__.

    Examples
    --------

    >>> pd.timedelta_range(start='1 day', periods=4)
    TimedeltaIndex(['1 days', '2 days', '3 days', '4 days'],
                   dtype='timedelta64[ns]', freq='D')

    The ``closed`` parameter specifies which endpoint is included.  The default
    behavior is to include both endpoints.

    >>> pd.timedelta_range(start='1 day', periods=4, closed='right')
    TimedeltaIndex(['2 days', '3 days', '4 days'],
                   dtype='timedelta64[ns]', freq='D')

    The ``freq`` parameter specifies the frequency of the TimedeltaIndex.
    Only fixed frequencies can be passed, non-fixed frequencies such as
    'M' (month end) will raise.

    >>> pd.timedelta_range(start='1 day', end='2 days', freq='6H')
    TimedeltaIndex(['1 days 00:00:00', '1 days 06:00:00', '1 days 12:00:00',
                    '1 days 18:00:00', '2 days 00:00:00'],
                   dtype='timedelta64[ns]', freq='6H')

    Specify ``start``, ``end``, and ``periods``; the frequency is generated
    automatically (linearly spaced).

    >>> pd.timedelta_range(start='1 day', end='5 days', periods=4)
    TimedeltaIndex(['1 days 00:00:00', '2 days 08:00:00', '3 days 16:00:00',
                '5 days 00:00:00'],
               dtype='timedelta64[ns]', freq=None)
    """
    if freq is None and com._any_none(periods, start, end):
        freq = 'D'

    return TimedeltaIndex(start=start, end=end, periods=periods,
                          freq=freq, name=name, closed=closed)
