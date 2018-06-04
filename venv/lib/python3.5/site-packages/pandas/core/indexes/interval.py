""" define the IntervalIndex """

import numpy as np
import warnings

from pandas.core.dtypes.missing import notna, isna
from pandas.core.dtypes.generic import ABCDatetimeIndex, ABCPeriodIndex
from pandas.core.dtypes.dtypes import IntervalDtype
from pandas.core.dtypes.cast import (
    maybe_convert_platform, find_common_type, maybe_downcast_to_dtype)
from pandas.core.dtypes.common import (
    _ensure_platform_int,
    is_list_like,
    is_datetime_or_timedelta_dtype,
    is_datetime64tz_dtype,
    is_categorical_dtype,
    is_string_dtype,
    is_integer_dtype,
    is_float_dtype,
    is_interval_dtype,
    is_object_dtype,
    is_scalar,
    is_float,
    is_number,
    is_integer,
    pandas_dtype)
from pandas.core.indexes.base import (
    Index, _ensure_index,
    default_pprint, _index_shared_docs)

from pandas._libs import Timestamp, Timedelta
from pandas._libs.interval import (
    Interval, IntervalMixin, IntervalTree,
    intervals_to_interval_bounds)

from pandas.core.indexes.datetimes import date_range
from pandas.core.indexes.timedeltas import timedelta_range
from pandas.core.indexes.multi import MultiIndex
from pandas.compat.numpy import function as nv
import pandas.core.common as com
from pandas.util._decorators import cache_readonly, Appender
from pandas.core.config import get_option
from pandas.tseries.frequencies import to_offset
from pandas.tseries.offsets import DateOffset

import pandas.core.indexes.base as ibase
_index_doc_kwargs = dict(ibase._index_doc_kwargs)
_index_doc_kwargs.update(
    dict(klass='IntervalIndex',
         target_klass='IntervalIndex or list of Intervals'))


_VALID_CLOSED = set(['left', 'right', 'both', 'neither'])


def _get_next_label(label):
    dtype = getattr(label, 'dtype', type(label))
    if isinstance(label, (Timestamp, Timedelta)):
        dtype = 'datetime64'
    if is_datetime_or_timedelta_dtype(dtype) or is_datetime64tz_dtype(dtype):
        return label + np.timedelta64(1, 'ns')
    elif is_integer_dtype(dtype):
        return label + 1
    elif is_float_dtype(dtype):
        return np.nextafter(label, np.infty)
    else:
        raise TypeError('cannot determine next label for type {typ!r}'
                        .format(typ=type(label)))


def _get_prev_label(label):
    dtype = getattr(label, 'dtype', type(label))
    if isinstance(label, (Timestamp, Timedelta)):
        dtype = 'datetime64'
    if is_datetime_or_timedelta_dtype(dtype) or is_datetime64tz_dtype(dtype):
        return label - np.timedelta64(1, 'ns')
    elif is_integer_dtype(dtype):
        return label - 1
    elif is_float_dtype(dtype):
        return np.nextafter(label, -np.infty)
    else:
        raise TypeError('cannot determine next label for type {typ!r}'
                        .format(typ=type(label)))


def _get_interval_closed_bounds(interval):
    """
    Given an Interval or IntervalIndex, return the corresponding interval with
    closed bounds.
    """
    left, right = interval.left, interval.right
    if interval.open_left:
        left = _get_next_label(left)
    if interval.open_right:
        right = _get_prev_label(right)
    return left, right


def maybe_convert_platform_interval(values):
    """
    Try to do platform conversion, with special casing for IntervalIndex.
    Wrapper around maybe_convert_platform that alters the default return
    dtype in certain cases to be compatible with IntervalIndex.  For example,
    empty lists return with integer dtype instead of object dtype, which is
    prohibited for IntervalIndex.

    Parameters
    ----------
    values : array-like

    Returns
    -------
    array
    """
    if isinstance(values, (list, tuple)) and len(values) == 0:
        # GH 19016
        # empty lists/tuples get object dtype by default, but this is not
        # prohibited for IntervalIndex, so coerce to integer instead
        return np.array([], dtype=np.int64)
    return maybe_convert_platform(values)


def _new_IntervalIndex(cls, d):
    """
    This is called upon unpickling, rather than the default which doesn't have
    arguments and breaks __new__
    """
    return cls.from_arrays(**d)


class IntervalIndex(IntervalMixin, Index):
    """
    Immutable Index implementing an ordered, sliceable set. IntervalIndex
    represents an Index of Interval objects that are all closed on the same
    side.

    .. versionadded:: 0.20.0

    .. warning::

       The indexing behaviors are provisional and may change in
       a future version of pandas.

    Parameters
    ----------
    data : array-like (1-dimensional)
        Array-like containing Interval objects from which to build the
        IntervalIndex
    closed : {'left', 'right', 'both', 'neither'}, default 'right'
        Whether the intervals are closed on the left-side, right-side, both or
        neither.
    name : object, optional
        Name to be stored in the index.
    copy : boolean, default False
        Copy the meta-data
    dtype : dtype or None, default None
        If None, dtype will be inferred

        ..versionadded:: 0.23.0

    Attributes
    ----------
    closed
    is_non_overlapping_monotonic
    left
    length
    mid
    right
    values

    Methods
    -------
    contains
    from_arrays
    from_breaks
    from_tuples
    get_indexer
    get_loc

    Examples
    ---------
    A new ``IntervalIndex`` is typically constructed using
    :func:`interval_range`:

    >>> pd.interval_range(start=0, end=5)
    IntervalIndex([(0, 1], (1, 2], (2, 3], (3, 4], (4, 5]]
                  closed='right', dtype='interval[int64]')

    It may also be constructed using one of the constructor
    methods: :meth:`IntervalIndex.from_arrays`,
    :meth:`IntervalIndex.from_breaks`, and :meth:`IntervalIndex.from_tuples`.

    See further examples in the doc strings of ``interval_range`` and the
    mentioned constructor methods.

    Notes
    ------
    See the `user guide
    <http://pandas.pydata.org/pandas-docs/stable/advanced.html#intervalindex>`_
    for more.

    See Also
    --------
    Index : The base pandas Index type
    Interval : A bounded slice-like interval; the elements of an IntervalIndex
    interval_range : Function to create a fixed frequency IntervalIndex
    cut, qcut : Convert arrays of continuous data into Categoricals/Series of
                Intervals
    """
    _typ = 'intervalindex'
    _comparables = ['name']
    _attributes = ['name', 'closed']

    # we would like our indexing holder to defer to us
    _defer_to_indexing = True

    _mask = None

    def __new__(cls, data, closed=None, dtype=None, copy=False,
                name=None, fastpath=False, verify_integrity=True):

        if fastpath:
            return cls._simple_new(data.left, data.right, closed, name,
                                   copy=copy, verify_integrity=False)

        if name is None and hasattr(data, 'name'):
            name = data.name

        if isinstance(data, IntervalIndex):
            left = data.left
            right = data.right
            closed = data.closed
        else:

            # don't allow scalars
            if is_scalar(data):
                cls._scalar_data_error(data)

            data = maybe_convert_platform_interval(data)
            left, right, infer_closed = intervals_to_interval_bounds(data)

            if (com._all_not_none(closed, infer_closed) and
                    closed != infer_closed):
                # GH 18421
                msg = ("conflicting values for closed: constructor got "
                       "'{closed}', inferred from data '{infer_closed}'"
                       .format(closed=closed, infer_closed=infer_closed))
                raise ValueError(msg)

            closed = closed or infer_closed

        return cls._simple_new(left, right, closed, name, copy=copy,
                               dtype=dtype, verify_integrity=verify_integrity)

    @classmethod
    def _simple_new(cls, left, right, closed=None, name=None, copy=False,
                    dtype=None, verify_integrity=True):
        result = IntervalMixin.__new__(cls)

        closed = closed or 'right'
        left = _ensure_index(left, copy=copy)
        right = _ensure_index(right, copy=copy)

        if dtype is not None:
            # GH 19262: dtype must be an IntervalDtype to override inferred
            dtype = pandas_dtype(dtype)
            if not is_interval_dtype(dtype):
                msg = 'dtype must be an IntervalDtype, got {dtype}'
                raise TypeError(msg.format(dtype=dtype))
            elif dtype.subtype is not None:
                left = left.astype(dtype.subtype)
                right = right.astype(dtype.subtype)

        # coerce dtypes to match if needed
        if is_float_dtype(left) and is_integer_dtype(right):
            right = right.astype(left.dtype)
        elif is_float_dtype(right) and is_integer_dtype(left):
            left = left.astype(right.dtype)

        if type(left) != type(right):
            msg = ('must not have differing left [{ltype}] and right '
                   '[{rtype}] types')
            raise ValueError(msg.format(ltype=type(left).__name__,
                                        rtype=type(right).__name__))
        elif is_categorical_dtype(left.dtype) or is_string_dtype(left.dtype):
            # GH 19016
            msg = ('category, object, and string subtypes are not supported '
                   'for IntervalIndex')
            raise TypeError(msg)
        elif isinstance(left, ABCPeriodIndex):
            msg = 'Period dtypes are not supported, use a PeriodIndex instead'
            raise ValueError(msg)
        elif (isinstance(left, ABCDatetimeIndex) and
                str(left.tz) != str(right.tz)):
            msg = ("left and right must have the same time zone, got "
                   "'{left_tz}' and '{right_tz}'")
            raise ValueError(msg.format(left_tz=left.tz, right_tz=right.tz))

        result._left = left
        result._right = right
        result._closed = closed
        result.name = name
        if verify_integrity:
            result._validate()
        result._reset_identity()
        return result

    @Appender(_index_shared_docs['_shallow_copy'])
    def _shallow_copy(self, left=None, right=None, **kwargs):
        if left is None:

            # no values passed
            left, right = self.left, self.right

        elif right is None:

            # only single value passed, could be an IntervalIndex
            # or array of Intervals
            if not isinstance(left, IntervalIndex):
                left = self._constructor(left)

            left, right = left.left, left.right
        else:

            # both left and right are values
            pass

        attributes = self._get_attributes_dict()
        attributes.update(kwargs)
        attributes['verify_integrity'] = False
        return self._simple_new(left, right, **attributes)

    def _validate(self):
        """
        Verify that the IntervalIndex is valid.
        """
        if self.closed not in _VALID_CLOSED:
            raise ValueError("invalid option for 'closed': {closed}"
                             .format(closed=self.closed))
        if len(self.left) != len(self.right):
            raise ValueError('left and right must have the same length')
        left_mask = notna(self.left)
        right_mask = notna(self.right)
        if not (left_mask == right_mask).all():
            raise ValueError('missing values must be missing in the same '
                             'location both left and right sides')
        if not (self.left[left_mask] <= self.right[left_mask]).all():
            raise ValueError('left side of interval must be <= right side')
        self._mask = ~left_mask

    @cache_readonly
    def hasnans(self):
        """
        Return if the IntervalIndex has any nans; enables various performance
        speedups
        """
        return self._isnan.any()

    @cache_readonly
    def _isnan(self):
        """Return a mask indicating if each value is NA"""
        if self._mask is None:
            self._mask = isna(self.left)
        return self._mask

    @cache_readonly
    def _engine(self):
        return IntervalTree(self.left, self.right, closed=self.closed)

    @property
    def _constructor(self):
        return type(self)

    def __contains__(self, key):
        """
        return a boolean if this key is IN the index
        We *only* accept an Interval

        Parameters
        ----------
        key : Interval

        Returns
        -------
        boolean
        """
        if not isinstance(key, Interval):
            return False

        try:
            self.get_loc(key)
            return True
        except KeyError:
            return False

    def contains(self, key):
        """
        Return a boolean indicating if the key is IN the index

        We accept / allow keys to be not *just* actual
        objects.

        Parameters
        ----------
        key : int, float, Interval

        Returns
        -------
        boolean
        """
        try:
            self.get_loc(key)
            return True
        except KeyError:
            return False

    @classmethod
    def from_breaks(cls, breaks, closed='right', name=None, copy=False,
                    dtype=None):
        """
        Construct an IntervalIndex from an array of splits

        Parameters
        ----------
        breaks : array-like (1-dimensional)
            Left and right bounds for each interval.
        closed : {'left', 'right', 'both', 'neither'}, default 'right'
            Whether the intervals are closed on the left-side, right-side, both
            or neither.
        name : object, optional
            Name to be stored in the index.
        copy : boolean, default False
            copy the data
        dtype : dtype or None, default None
            If None, dtype will be inferred

            ..versionadded:: 0.23.0

        Examples
        --------
        >>> pd.IntervalIndex.from_breaks([0, 1, 2, 3])
        IntervalIndex([(0, 1], (1, 2], (2, 3]]
                      closed='right',
                      dtype='interval[int64]')

        See Also
        --------
        interval_range : Function to create a fixed frequency IntervalIndex
        IntervalIndex.from_arrays : Construct an IntervalIndex from a left and
                                    right array
        IntervalIndex.from_tuples : Construct an IntervalIndex from a
                                    list/array of tuples
        """
        breaks = maybe_convert_platform_interval(breaks)

        return cls.from_arrays(breaks[:-1], breaks[1:], closed,
                               name=name, copy=copy, dtype=dtype)

    @classmethod
    def from_arrays(cls, left, right, closed='right', name=None, copy=False,
                    dtype=None):
        """
        Construct from two arrays defining the left and right bounds.

        Parameters
        ----------
        left : array-like (1-dimensional)
            Left bounds for each interval.
        right : array-like (1-dimensional)
            Right bounds for each interval.
        closed : {'left', 'right', 'both', 'neither'}, default 'right'
            Whether the intervals are closed on the left-side, right-side, both
            or neither.
        name : object, optional
            Name to be stored in the index.
        copy : boolean, default False
            Copy the data.
        dtype : dtype, optional
            If None, dtype will be inferred.

            .. versionadded:: 0.23.0

        Returns
        -------
        index : IntervalIndex

        Notes
        -----
        Each element of `left` must be less than or equal to the `right`
        element at the same position. If an element is missing, it must be
        missing in both `left` and `right`. A TypeError is raised when
        using an unsupported type for `left` or `right`. At the moment,
        'category', 'object', and 'string' subtypes are not supported.

        Raises
        ------
        ValueError
            When a value is missing in only one of `left` or `right`.
            When a value in `left` is greater than the corresponding value
            in `right`.

        See Also
        --------
        interval_range : Function to create a fixed frequency IntervalIndex.
        IntervalIndex.from_breaks : Construct an IntervalIndex from an array of
            splits.
        IntervalIndex.from_tuples : Construct an IntervalIndex from a
            list/array of tuples.

        Examples
        --------
        >>> pd.IntervalIndex.from_arrays([0, 1, 2], [1, 2, 3])
        IntervalIndex([(0, 1], (1, 2], (2, 3]]
                      closed='right',
                      dtype='interval[int64]')

        If you want to segment different groups of people based on
        ages, you can apply the method as follows:

        >>> ages = pd.IntervalIndex.from_arrays([0, 2, 13],
        ...                                     [2, 13, 19], closed='left')
        >>> ages
        IntervalIndex([[0, 2), [2, 13), [13, 19)]
                      closed='left',
                      dtype='interval[int64]')
        >>> s = pd.Series(['baby', 'kid', 'teen'], ages)
        >>> s
        [0, 2)      baby
        [2, 13)      kid
        [13, 19)    teen
        dtype: object

        Values may be missing, but they must be missing in both arrays.

        >>> pd.IntervalIndex.from_arrays([0, np.nan, 13],
        ...                              [2, np.nan, 19])
        IntervalIndex([(0.0, 2.0], nan, (13.0, 19.0]]
                      closed='right',
                      dtype='interval[float64]')
        """
        left = maybe_convert_platform_interval(left)
        right = maybe_convert_platform_interval(right)

        return cls._simple_new(left, right, closed, name=name, copy=copy,
                               dtype=dtype, verify_integrity=True)

    @classmethod
    def from_intervals(cls, data, closed=None, name=None, copy=False,
                       dtype=None):
        """
        Construct an IntervalIndex from a 1d array of Interval objects

        .. deprecated:: 0.23.0

        Parameters
        ----------
        data : array-like (1-dimensional)
            Array of Interval objects. All intervals must be closed on the same
            sides.
        name : object, optional
            Name to be stored in the index.
        copy : boolean, default False
            by-default copy the data, this is compat only and ignored
        dtype : dtype or None, default None
            If None, dtype will be inferred

            ..versionadded:: 0.23.0

        Examples
        --------
        >>> pd.IntervalIndex.from_intervals([pd.Interval(0, 1),
        ...                                  pd.Interval(1, 2)])
        IntervalIndex([(0, 1], (1, 2]]
                      closed='right', dtype='interval[int64]')

        The generic Index constructor work identically when it infers an array
        of all intervals:

        >>> pd.Index([pd.Interval(0, 1), pd.Interval(1, 2)])
        IntervalIndex([(0, 1], (1, 2]]
                      closed='right', dtype='interval[int64]')

        See Also
        --------
        interval_range : Function to create a fixed frequency IntervalIndex
        IntervalIndex.from_arrays : Construct an IntervalIndex from a left and
                                    right array
        IntervalIndex.from_breaks : Construct an IntervalIndex from an array of
                                    splits
        IntervalIndex.from_tuples : Construct an IntervalIndex from a
                                    list/array of tuples
        """
        msg = ('IntervalIndex.from_intervals is deprecated and will be '
               'removed in a future version; use IntervalIndex(...) instead')
        warnings.warn(msg, FutureWarning, stacklevel=2)
        return cls(data, closed=closed, name=name, copy=copy, dtype=dtype)

    @classmethod
    def from_tuples(cls, data, closed='right', name=None, copy=False,
                    dtype=None):
        """
        Construct an IntervalIndex from a list/array of tuples

        Parameters
        ----------
        data : array-like (1-dimensional)
            Array of tuples
        closed : {'left', 'right', 'both', 'neither'}, default 'right'
            Whether the intervals are closed on the left-side, right-side, both
            or neither.
        name : object, optional
            Name to be stored in the index.
        copy : boolean, default False
            by-default copy the data, this is compat only and ignored
        dtype : dtype or None, default None
            If None, dtype will be inferred

            ..versionadded:: 0.23.0

        Examples
        --------
        >>>  pd.IntervalIndex.from_tuples([(0, 1), (1, 2)])
        IntervalIndex([(0, 1], (1, 2]],
                      closed='right', dtype='interval[int64]')

        See Also
        --------
        interval_range : Function to create a fixed frequency IntervalIndex
        IntervalIndex.from_arrays : Construct an IntervalIndex from a left and
                                    right array
        IntervalIndex.from_breaks : Construct an IntervalIndex from an array of
                                    splits
        """
        if len(data):
            left, right = [], []
        else:
            left = right = data

        for d in data:
            if isna(d):
                lhs = rhs = np.nan
            else:
                try:
                    # need list of length 2 tuples, e.g. [(0, 1), (1, 2), ...]
                    lhs, rhs = d
                except ValueError:
                    msg = ('IntervalIndex.from_tuples requires tuples of '
                           'length 2, got {tpl}').format(tpl=d)
                    raise ValueError(msg)
                except TypeError:
                    msg = ('IntervalIndex.from_tuples received an invalid '
                           'item, {tpl}').format(tpl=d)
                    raise TypeError(msg)
            left.append(lhs)
            right.append(rhs)

        return cls.from_arrays(left, right, closed, name=name, copy=False,
                               dtype=dtype)

    def to_tuples(self, na_tuple=True):
        """
        Return an Index of tuples of the form (left, right)

        Parameters
        ----------
        na_tuple : boolean, default True
            Returns NA as a tuple if True, ``(nan, nan)``, or just as the NA
            value itself if False, ``nan``.

            ..versionadded:: 0.23.0

        Examples
        --------
        >>>  idx = pd.IntervalIndex.from_arrays([0, np.nan, 2], [1, np.nan, 3])
        >>>  idx.to_tuples()
        Index([(0.0, 1.0), (nan, nan), (2.0, 3.0)], dtype='object')
        >>>  idx.to_tuples(na_tuple=False)
        Index([(0.0, 1.0), nan, (2.0, 3.0)], dtype='object')
        """
        tuples = com._asarray_tuplesafe(zip(self.left, self.right))
        if not na_tuple:
            # GH 18756
            tuples = np.where(~self._isnan, tuples, np.nan)
        return Index(tuples)

    @cache_readonly
    def _multiindex(self):
        return MultiIndex.from_arrays([self.left, self.right],
                                      names=['left', 'right'])

    @property
    def left(self):
        """
        Return the left endpoints of each Interval in the IntervalIndex as
        an Index
        """
        return self._left

    @property
    def right(self):
        """
        Return the right endpoints of each Interval in the IntervalIndex as
        an Index
        """
        return self._right

    @property
    def closed(self):
        """
        Whether the intervals are closed on the left-side, right-side, both or
        neither
        """
        return self._closed

    @property
    def length(self):
        """
        Return an Index with entries denoting the length of each Interval in
        the IntervalIndex
        """
        try:
            return self.right - self.left
        except TypeError:
            # length not defined for some types, e.g. string
            msg = ('IntervalIndex contains Intervals without defined length, '
                   'e.g. Intervals with string endpoints')
            raise TypeError(msg)

    @property
    def size(self):
        # Avoid materializing self.values
        return self.left.size

    @property
    def shape(self):
        # Avoid materializing self.values
        return self.left.shape

    def __len__(self):
        return len(self.left)

    @cache_readonly
    def values(self):
        """
        Return the IntervalIndex's data as a numpy array of Interval
        objects (with dtype='object')
        """
        left = self.left
        right = self.right
        mask = self._isnan
        closed = self._closed

        result = np.empty(len(left), dtype=object)
        for i in range(len(left)):
            if mask[i]:
                result[i] = np.nan
            else:
                result[i] = Interval(left[i], right[i], closed)
        return result

    def __array__(self, result=None):
        """ the array interface, return my values """
        return self.values

    def __array_wrap__(self, result, context=None):
        # we don't want the superclass implementation
        return result

    def _array_values(self):
        return self.values

    def __reduce__(self):
        d = dict(left=self.left,
                 right=self.right)
        d.update(self._get_attributes_dict())
        return _new_IntervalIndex, (self.__class__, d), None

    @Appender(_index_shared_docs['copy'])
    def copy(self, deep=False, name=None):
        left = self.left.copy(deep=True) if deep else self.left
        right = self.right.copy(deep=True) if deep else self.right
        name = name if name is not None else self.name
        closed = self.closed
        return type(self).from_arrays(left, right, closed=closed, name=name)

    @Appender(_index_shared_docs['astype'])
    def astype(self, dtype, copy=True):
        dtype = pandas_dtype(dtype)
        if is_interval_dtype(dtype) and dtype != self.dtype:
            try:
                new_left = self.left.astype(dtype.subtype)
                new_right = self.right.astype(dtype.subtype)
            except TypeError:
                msg = ('Cannot convert {dtype} to {new_dtype}; subtypes are '
                       'incompatible')
                raise TypeError(msg.format(dtype=self.dtype, new_dtype=dtype))
            return self._shallow_copy(new_left, new_right)
        return super(IntervalIndex, self).astype(dtype, copy=copy)

    @cache_readonly
    def dtype(self):
        """Return the dtype object of the underlying data"""
        return IntervalDtype.construct_from_string(str(self.left.dtype))

    @property
    def inferred_type(self):
        """Return a string of the type inferred from the values"""
        return 'interval'

    @Appender(Index.memory_usage.__doc__)
    def memory_usage(self, deep=False):
        # we don't use an explicit engine
        # so return the bytes here
        return (self.left.memory_usage(deep=deep) +
                self.right.memory_usage(deep=deep))

    @cache_readonly
    def mid(self):
        """
        Return the midpoint of each Interval in the IntervalIndex as an Index
        """
        try:
            return 0.5 * (self.left + self.right)
        except TypeError:
            # datetime safe version
            return self.left + 0.5 * self.length

    @cache_readonly
    def is_monotonic(self):
        """
        Return True if the IntervalIndex is monotonic increasing (only equal or
        increasing values), else False
        """
        return self._multiindex.is_monotonic

    @cache_readonly
    def is_monotonic_increasing(self):
        """
        Return True if the IntervalIndex is monotonic increasing (only equal or
        increasing values), else False
        """
        return self._multiindex.is_monotonic_increasing

    @cache_readonly
    def is_monotonic_decreasing(self):
        """
        Return True if the IntervalIndex is monotonic decreasing (only equal or
        decreasing values), else False
        """
        return self._multiindex.is_monotonic_decreasing

    @cache_readonly
    def is_unique(self):
        """
        Return True if the IntervalIndex contains unique elements, else False
        """
        return self._multiindex.is_unique

    @cache_readonly
    def is_non_overlapping_monotonic(self):
        """
        Return True if the IntervalIndex is non-overlapping (no Intervals share
        points) and is either monotonic increasing or monotonic decreasing,
        else False
        """
        # must be increasing  (e.g., [0, 1), [1, 2), [2, 3), ... )
        # or decreasing (e.g., [-1, 0), [-2, -1), [-3, -2), ...)
        # we already require left <= right

        # strict inequality for closed == 'both'; equality implies overlapping
        # at a point when both sides of intervals are included
        if self.closed == 'both':
            return bool((self.right[:-1] < self.left[1:]).all() or
                        (self.left[:-1] > self.right[1:]).all())

        # non-strict inequality when closed != 'both'; at least one side is
        # not included in the intervals, so equality does not imply overlapping
        return bool((self.right[:-1] <= self.left[1:]).all() or
                    (self.left[:-1] >= self.right[1:]).all())

    @Appender(_index_shared_docs['_convert_scalar_indexer'])
    def _convert_scalar_indexer(self, key, kind=None):
        if kind == 'iloc':
            return super(IntervalIndex, self)._convert_scalar_indexer(
                key, kind=kind)
        return key

    def _maybe_cast_slice_bound(self, label, side, kind):
        return getattr(self, side)._maybe_cast_slice_bound(label, side, kind)

    @Appender(_index_shared_docs['_convert_list_indexer'])
    def _convert_list_indexer(self, keyarr, kind=None):
        """
        we are passed a list-like indexer. Return the
        indexer for matching intervals.
        """
        locs = self.get_indexer_for(keyarr)

        # we have missing values
        if (locs == -1).any():
            raise KeyError

        return locs

    def _maybe_cast_indexed(self, key):
        """
        we need to cast the key, which could be a scalar
        or an array-like to the type of our subtype
        """
        if isinstance(key, IntervalIndex):
            return key

        subtype = self.dtype.subtype
        if is_float_dtype(subtype):
            if is_integer(key):
                key = float(key)
            elif isinstance(key, (np.ndarray, Index)):
                key = key.astype('float64')
        elif is_integer_dtype(subtype):
            if is_integer(key):
                key = int(key)

        return key

    def _check_method(self, method):
        if method is None:
            return

        if method in ['bfill', 'backfill', 'pad', 'ffill', 'nearest']:
            msg = 'method {method} not yet implemented for IntervalIndex'
            raise NotImplementedError(msg.format(method=method))

        raise ValueError("Invalid fill method")

    def _searchsorted_monotonic(self, label, side, exclude_label=False):
        if not self.is_non_overlapping_monotonic:
            raise KeyError('can only get slices from an IntervalIndex if '
                           'bounds are non-overlapping and all monotonic '
                           'increasing or decreasing')

        if isinstance(label, IntervalMixin):
            raise NotImplementedError

        # GH 20921: "not is_monotonic_increasing" for the second condition
        # instead of "is_monotonic_decreasing" to account for single element
        # indexes being both increasing and decreasing
        if ((side == 'left' and self.left.is_monotonic_increasing) or
                (side == 'right' and not self.left.is_monotonic_increasing)):
            sub_idx = self.right
            if self.open_right or exclude_label:
                label = _get_next_label(label)
        else:
            sub_idx = self.left
            if self.open_left or exclude_label:
                label = _get_prev_label(label)

        return sub_idx._searchsorted_monotonic(label, side)

    def _get_loc_only_exact_matches(self, key):
        if isinstance(key, Interval):

            if not self.is_unique:
                raise ValueError("cannot index with a slice Interval"
                                 " and a non-unique index")

            # TODO: this expands to a tuple index, see if we can
            # do better
            return Index(self._multiindex.values).get_loc(key)
        raise KeyError

    def _find_non_overlapping_monotonic_bounds(self, key):
        if isinstance(key, IntervalMixin):
            start = self._searchsorted_monotonic(
                key.left, 'left', exclude_label=key.open_left)
            stop = self._searchsorted_monotonic(
                key.right, 'right', exclude_label=key.open_right)
        elif isinstance(key, slice):
            # slice
            start, stop = key.start, key.stop
            if (key.step or 1) != 1:
                raise NotImplementedError("cannot slice with a slice step")
            if start is None:
                start = 0
            else:
                start = self._searchsorted_monotonic(start, 'left')
            if stop is None:
                stop = len(self)
            else:
                stop = self._searchsorted_monotonic(stop, 'right')
        else:
            # scalar or index-like

            start = self._searchsorted_monotonic(key, 'left')
            stop = self._searchsorted_monotonic(key, 'right')
        return start, stop

    def get_loc(self, key, method=None):
        """Get integer location, slice or boolean mask for requested label.

        Parameters
        ----------
        key : label
        method : {None}, optional
            * default: matches where the label is within an interval only.

        Returns
        -------
        loc : int if unique index, slice if monotonic index, else mask

        Examples
        ---------
        >>> i1, i2 = pd.Interval(0, 1), pd.Interval(1, 2)
        >>> index = pd.IntervalIndex([i1, i2])
        >>> index.get_loc(1)
        0

        You can also supply an interval or an location for a point inside an
        interval.

        >>> index.get_loc(pd.Interval(0, 2))
        array([0, 1], dtype=int64)
        >>> index.get_loc(1.5)
        1

        If a label is in several intervals, you get the locations of all the
        relevant intervals.

        >>> i3 = pd.Interval(0, 2)
        >>> overlapping_index = pd.IntervalIndex([i2, i3])
        >>> overlapping_index.get_loc(1.5)
        array([0, 1], dtype=int64)
        """
        self._check_method(method)

        original_key = key
        key = self._maybe_cast_indexed(key)

        if self.is_non_overlapping_monotonic:
            if isinstance(key, Interval):
                left = self._maybe_cast_slice_bound(key.left, 'left', None)
                right = self._maybe_cast_slice_bound(key.right, 'right', None)
                key = Interval(left, right, key.closed)
            else:
                key = self._maybe_cast_slice_bound(key, 'left', None)

            start, stop = self._find_non_overlapping_monotonic_bounds(key)

            if start is None or stop is None:
                return slice(start, stop)
            elif start + 1 == stop:
                return start
            elif start < stop:
                return slice(start, stop)
            else:
                raise KeyError(original_key)

        else:
            # use the interval tree
            if isinstance(key, Interval):
                left, right = _get_interval_closed_bounds(key)
                return self._engine.get_loc_interval(left, right)
            else:
                return self._engine.get_loc(key)

    def get_value(self, series, key):
        if com.is_bool_indexer(key):
            loc = key
        elif is_list_like(key):
            loc = self.get_indexer(key)
        elif isinstance(key, slice):

            if not (key.step is None or key.step == 1):
                raise ValueError("cannot support not-default step in a slice")

            try:
                loc = self.get_loc(key)
            except TypeError:
                # we didn't find exact intervals or are non-unique
                msg = "unable to slice with this key: {key}".format(key=key)
                raise ValueError(msg)

        else:
            loc = self.get_loc(key)
        return series.iloc[loc]

    @Appender(_index_shared_docs['get_indexer'] % _index_doc_kwargs)
    def get_indexer(self, target, method=None, limit=None, tolerance=None):

        self._check_method(method)
        target = _ensure_index(target)
        target = self._maybe_cast_indexed(target)

        if self.equals(target):
            return np.arange(len(self), dtype='intp')

        if self.is_non_overlapping_monotonic:
            start, stop = self._find_non_overlapping_monotonic_bounds(target)

            start_plus_one = start + 1
            if not ((start_plus_one < stop).any()):
                return np.where(start_plus_one == stop, start, -1)

        if not self.is_unique:
            raise ValueError("cannot handle non-unique indices")

        # IntervalIndex
        if isinstance(target, IntervalIndex):
            indexer = self._get_reindexer(target)

        # non IntervalIndex
        else:
            indexer = np.concatenate([self.get_loc(i) for i in target])

        return _ensure_platform_int(indexer)

    def _get_reindexer(self, target):
        """
        Return an indexer for a target IntervalIndex with self
        """

        # find the left and right indexers
        lindexer = self._engine.get_indexer(target.left.values)
        rindexer = self._engine.get_indexer(target.right.values)

        # we want to return an indexer on the intervals
        # however, our keys could provide overlapping of multiple
        # intervals, so we iterate thru the indexers and construct
        # a set of indexers

        indexer = []
        n = len(self)

        for i, (lhs, rhs) in enumerate(zip(lindexer, rindexer)):

            target_value = target[i]

            # matching on the lhs bound
            if (lhs != -1 and
                    self.closed == 'right' and
                    target_value.left == self[lhs].right):
                lhs += 1

            # matching on the lhs bound
            if (rhs != -1 and
                    self.closed == 'left' and
                    target_value.right == self[rhs].left):
                rhs -= 1

            # not found
            if lhs == -1 and rhs == -1:
                indexer.append(np.array([-1]))

            elif rhs == -1:

                indexer.append(np.arange(lhs, n))

            elif lhs == -1:

                # care about left/right closed here
                value = self[i]

                # target.closed same as self.closed
                if self.closed == target.closed:
                    if target_value.left < value.left:
                        indexer.append(np.array([-1]))
                        continue

                # target.closed == 'left'
                elif self.closed == 'right':
                    if target_value.left <= value.left:
                        indexer.append(np.array([-1]))
                        continue

                # target.closed == 'right'
                elif self.closed == 'left':
                    if target_value.left <= value.left:
                        indexer.append(np.array([-1]))
                        continue

                indexer.append(np.arange(0, rhs + 1))

            else:
                indexer.append(np.arange(lhs, rhs + 1))

        return np.concatenate(indexer)

    @Appender(_index_shared_docs['get_indexer_non_unique'] % _index_doc_kwargs)
    def get_indexer_non_unique(self, target):
        target = self._maybe_cast_indexed(_ensure_index(target))
        return super(IntervalIndex, self).get_indexer_non_unique(target)

    @Appender(_index_shared_docs['where'])
    def where(self, cond, other=None):
        if other is None:
            other = self._na_value
        values = np.where(cond, self.values, other)
        return self._shallow_copy(values)

    def delete(self, loc):
        """
        Return a new IntervalIndex with passed location(-s) deleted

        Returns
        -------
        new_index : IntervalIndex
        """
        new_left = self.left.delete(loc)
        new_right = self.right.delete(loc)
        return self._shallow_copy(new_left, new_right)

    def insert(self, loc, item):
        """
        Return a new IntervalIndex inserting new item at location. Follows
        Python list.append semantics for negative values.  Only Interval
        objects and NA can be inserted into an IntervalIndex

        Parameters
        ----------
        loc : int
        item : object

        Returns
        -------
        new_index : IntervalIndex
        """
        if isinstance(item, Interval):
            if item.closed != self.closed:
                raise ValueError('inserted item must be closed on the same '
                                 'side as the index')
            left_insert = item.left
            right_insert = item.right
        elif is_scalar(item) and isna(item):
            # GH 18295
            left_insert = right_insert = item
        else:
            raise ValueError('can only insert Interval objects and NA into '
                             'an IntervalIndex')

        new_left = self.left.insert(loc, left_insert)
        new_right = self.right.insert(loc, right_insert)
        return self._shallow_copy(new_left, new_right)

    def _as_like_interval_index(self, other):
        self._assert_can_do_setop(other)
        other = _ensure_index(other)
        if not isinstance(other, IntervalIndex):
            msg = ('the other index needs to be an IntervalIndex too, but '
                   'was type {}').format(other.__class__.__name__)
            raise TypeError(msg)
        elif self.closed != other.closed:
            msg = ('can only do set operations between two IntervalIndex '
                   'objects that are closed on the same side')
            raise ValueError(msg)
        return other

    def _concat_same_dtype(self, to_concat, name):
        """
        assert that we all have the same .closed
        we allow a 0-len index here as well
        """
        if not len({i.closed for i in to_concat if len(i)}) == 1:
            msg = ('can only append two IntervalIndex objects '
                   'that are closed on the same side')
            raise ValueError(msg)
        return super(IntervalIndex, self)._concat_same_dtype(to_concat, name)

    @Appender(_index_shared_docs['take'] % _index_doc_kwargs)
    def take(self, indices, axis=0, allow_fill=True,
             fill_value=None, **kwargs):
        nv.validate_take(tuple(), kwargs)
        indices = _ensure_platform_int(indices)
        left, right = self.left, self.right

        if fill_value is None:
            fill_value = self._na_value
        mask = indices == -1

        if not mask.any():
            # we won't change dtype here in this case
            # if we don't need
            allow_fill = False

        taker = lambda x: x.take(indices, allow_fill=allow_fill,
                                 fill_value=fill_value)

        try:
            new_left = taker(left)
            new_right = taker(right)
        except ValueError:

            # we need to coerce; migth have NA's in an
            # integer dtype
            new_left = taker(left.astype(float))
            new_right = taker(right.astype(float))

        return self._shallow_copy(new_left, new_right)

    def __getitem__(self, value):
        mask = self._isnan[value]
        if is_scalar(mask) and mask:
            return self._na_value

        left = self.left[value]
        right = self.right[value]

        # scalar
        if not isinstance(left, Index):
            return Interval(left, right, self.closed)

        return self._shallow_copy(left, right)

    # __repr__ associated methods are based on MultiIndex

    def _format_with_header(self, header, **kwargs):
        return header + list(self._format_native_types(**kwargs))

    def _format_native_types(self, na_rep='', quoting=None, **kwargs):
        """ actually format my specific types """
        from pandas.io.formats.format import IntervalArrayFormatter
        return IntervalArrayFormatter(values=self,
                                      na_rep=na_rep,
                                      justify='all').get_result()

    def _format_data(self, name=None):

        # TODO: integrate with categorical and make generic
        # name argument is unused here; just for compat with base / categorical
        n = len(self)
        max_seq_items = min((get_option(
            'display.max_seq_items') or n) // 10, 10)

        formatter = str

        if n == 0:
            summary = '[]'
        elif n == 1:
            first = formatter(self[0])
            summary = '[{first}]'.format(first=first)
        elif n == 2:
            first = formatter(self[0])
            last = formatter(self[-1])
            summary = '[{first}, {last}]'.format(first=first, last=last)
        else:

            if n > max_seq_items:
                n = min(max_seq_items // 2, 10)
                head = [formatter(x) for x in self[:n]]
                tail = [formatter(x) for x in self[-n:]]
                summary = '[{head} ... {tail}]'.format(
                    head=', '.join(head), tail=', '.join(tail))
            else:
                head = []
                tail = [formatter(x) for x in self]
                summary = '[{tail}]'.format(tail=', '.join(tail))

        return summary + self._format_space()

    def _format_attrs(self):
        attrs = [('closed', repr(self.closed))]
        if self.name is not None:
            attrs.append(('name', default_pprint(self.name)))
        attrs.append(('dtype', "'{dtype}'".format(dtype=self.dtype)))
        return attrs

    def _format_space(self):
        space = ' ' * (len(self.__class__.__name__) + 1)
        return "\n{space}".format(space=space)

    def argsort(self, *args, **kwargs):
        return np.lexsort((self.right, self.left))

    def equals(self, other):
        """
        Determines if two IntervalIndex objects contain the same elements
        """
        if self.is_(other):
            return True

        # if we can coerce to an II
        # then we can compare
        if not isinstance(other, IntervalIndex):
            if not is_interval_dtype(other):
                return False
            other = Index(getattr(other, '.values', other))

        return (self.left.equals(other.left) and
                self.right.equals(other.right) and
                self.closed == other.closed)

    def _setop(op_name):
        def func(self, other):
            other = self._as_like_interval_index(other)

            # GH 19016: ensure set op will not return a prohibited dtype
            subtypes = [self.dtype.subtype, other.dtype.subtype]
            common_subtype = find_common_type(subtypes)
            if is_object_dtype(common_subtype):
                msg = ('can only do {op} between two IntervalIndex '
                       'objects that have compatible dtypes')
                raise TypeError(msg.format(op=op_name))

            result = getattr(self._multiindex, op_name)(other._multiindex)
            result_name = self.name if self.name == other.name else None

            # GH 19101: ensure empty results have correct dtype
            if result.empty:
                result = result.values.astype(self.dtype.subtype)
            else:
                result = result.values

            return type(self).from_tuples(result, closed=self.closed,
                                          name=result_name)
        return func

    union = _setop('union')
    intersection = _setop('intersection')
    difference = _setop('difference')
    symmetric_difference = _setop('symmetric_difference')

    # TODO: arithmetic operations


IntervalIndex._add_logical_methods_disabled()


def _is_valid_endpoint(endpoint):
    """helper for interval_range to check if start/end are valid types"""
    return any([is_number(endpoint),
                isinstance(endpoint, Timestamp),
                isinstance(endpoint, Timedelta),
                endpoint is None])


def _is_type_compatible(a, b):
    """helper for interval_range to check type compat of start/end/freq"""
    is_ts_compat = lambda x: isinstance(x, (Timestamp, DateOffset))
    is_td_compat = lambda x: isinstance(x, (Timedelta, DateOffset))
    return ((is_number(a) and is_number(b)) or
            (is_ts_compat(a) and is_ts_compat(b)) or
            (is_td_compat(a) and is_td_compat(b)) or
            com._any_none(a, b))


def interval_range(start=None, end=None, periods=None, freq=None,
                   name=None, closed='right'):
    """
    Return a fixed frequency IntervalIndex

    Parameters
    ----------
    start : numeric or datetime-like, default None
        Left bound for generating intervals
    end : numeric or datetime-like, default None
        Right bound for generating intervals
    periods : integer, default None
        Number of periods to generate
    freq : numeric, string, or DateOffset, default None
        The length of each interval. Must be consistent with the type of start
        and end, e.g. 2 for numeric, or '5H' for datetime-like.  Default is 1
        for numeric and 'D' (calendar daily) for datetime-like.
    name : string, default None
        Name of the resulting IntervalIndex
    closed : {'left', 'right', 'both', 'neither'}, default 'right'
        Whether the intervals are closed on the left-side, right-side, both
        or neither.

    Notes
    -----
    Of the four parameters ``start``, ``end``, ``periods``, and ``freq``,
    exactly three must be specified. If ``freq`` is omitted, the resulting
    ``IntervalIndex`` will have ``periods`` linearly spaced elements between
    ``start`` and ``end``, inclusively.

    To learn more about datetime-like frequency strings, please see `this link
    <http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases>`__.

    Returns
    -------
    rng : IntervalIndex

    Examples
    --------
    Numeric ``start`` and  ``end`` is supported.

    >>> pd.interval_range(start=0, end=5)
    IntervalIndex([(0, 1], (1, 2], (2, 3], (3, 4], (4, 5]]
                  closed='right', dtype='interval[int64]')

    Additionally, datetime-like input is also supported.

    >>> pd.interval_range(start=pd.Timestamp('2017-01-01'),
                          end=pd.Timestamp('2017-01-04'))
    IntervalIndex([(2017-01-01, 2017-01-02], (2017-01-02, 2017-01-03],
                   (2017-01-03, 2017-01-04]]
                  closed='right', dtype='interval[datetime64[ns]]')

    The ``freq`` parameter specifies the frequency between the left and right.
    endpoints of the individual intervals within the ``IntervalIndex``.  For
    numeric ``start`` and ``end``, the frequency must also be numeric.

    >>> pd.interval_range(start=0, periods=4, freq=1.5)
    IntervalIndex([(0.0, 1.5], (1.5, 3.0], (3.0, 4.5], (4.5, 6.0]]
                  closed='right', dtype='interval[float64]')

    Similarly, for datetime-like ``start`` and ``end``, the frequency must be
    convertible to a DateOffset.

    >>> pd.interval_range(start=pd.Timestamp('2017-01-01'),
                          periods=3, freq='MS')
    IntervalIndex([(2017-01-01, 2017-02-01], (2017-02-01, 2017-03-01],
                   (2017-03-01, 2017-04-01]]
                  closed='right', dtype='interval[datetime64[ns]]')

    Specify ``start``, ``end``, and ``periods``; the frequency is generated
    automatically (linearly spaced).

    >>> pd.interval_range(start=0, end=6, periods=4)
    IntervalIndex([(0.0, 1.5], (1.5, 3.0], (3.0, 4.5], (4.5, 6.0]]
              closed='right',
              dtype='interval[float64]')

    The ``closed`` parameter specifies which endpoints of the individual
    intervals within the ``IntervalIndex`` are closed.

    >>> pd.interval_range(end=5, periods=4, closed='both')
    IntervalIndex([[1, 2], [2, 3], [3, 4], [4, 5]]
                  closed='both', dtype='interval[int64]')

    See Also
    --------
    IntervalIndex : an Index of intervals that are all closed on the same side.
    """
    start = com._maybe_box_datetimelike(start)
    end = com._maybe_box_datetimelike(end)
    endpoint = start if start is not None else end

    if freq is None and com._any_none(periods, start, end):
        freq = 1 if is_number(endpoint) else 'D'

    if com._count_not_none(start, end, periods, freq) != 3:
        raise ValueError('Of the four parameters: start, end, periods, and '
                         'freq, exactly three must be specified')

    if not _is_valid_endpoint(start):
        msg = 'start must be numeric or datetime-like, got {start}'
        raise ValueError(msg.format(start=start))
    elif not _is_valid_endpoint(end):
        msg = 'end must be numeric or datetime-like, got {end}'
        raise ValueError(msg.format(end=end))

    if is_float(periods):
        periods = int(periods)
    elif not is_integer(periods) and periods is not None:
        msg = 'periods must be a number, got {periods}'
        raise TypeError(msg.format(periods=periods))

    if freq is not None and not is_number(freq):
        try:
            freq = to_offset(freq)
        except ValueError:
            raise ValueError('freq must be numeric or convertible to '
                             'DateOffset, got {freq}'.format(freq=freq))

    # verify type compatibility
    if not all([_is_type_compatible(start, end),
                _is_type_compatible(start, freq),
                _is_type_compatible(end, freq)]):
        raise TypeError("start, end, freq need to be type compatible")

    # +1 to convert interval count to breaks count (n breaks = n-1 intervals)
    if periods is not None:
        periods += 1

    if is_number(endpoint):
        # compute the period/start/end if unspecified (at most one)
        if periods is None:
            periods = int((end - start) // freq) + 1
        elif start is None:
            start = end - (periods - 1) * freq
        elif end is None:
            end = start + (periods - 1) * freq

        # force end to be consistent with freq (lower if freq skips end)
        if freq is not None:
            end -= end % freq

        breaks = np.linspace(start, end, periods)
        if all(is_integer(x) for x in com._not_none(start, end, freq)):
            # np.linspace always produces float output
            breaks = maybe_downcast_to_dtype(breaks, 'int64')
    else:
        # delegate to the appropriate range function
        if isinstance(endpoint, Timestamp):
            range_func = date_range
        else:
            range_func = timedelta_range

        breaks = range_func(start=start, end=end, periods=periods, freq=freq)

    return IntervalIndex.from_breaks(breaks, name=name, closed=closed)
