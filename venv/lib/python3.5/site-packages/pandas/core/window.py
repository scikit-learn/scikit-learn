"""

provide a generic structure to support window functions,
similar to how we have a Groupby object


"""
from __future__ import division

import warnings
import numpy as np
from collections import defaultdict
from datetime import timedelta

from pandas.core.dtypes.generic import (
    ABCSeries,
    ABCDataFrame,
    ABCDatetimeIndex,
    ABCTimedeltaIndex,
    ABCPeriodIndex,
    ABCDateOffset)
from pandas.core.dtypes.common import (
    is_integer,
    is_bool,
    is_float_dtype,
    is_integer_dtype,
    needs_i8_conversion,
    is_timedelta64_dtype,
    is_list_like,
    _ensure_float64,
    is_scalar)

from pandas.core.base import (PandasObject, SelectionMixin,
                              GroupByMixin)
import pandas.core.common as com
import pandas._libs.window as _window

from pandas import compat
from pandas.compat.numpy import function as nv
from pandas.util._decorators import (Substitution, Appender,
                                     cache_readonly)
from pandas.core.generic import _shared_docs
from textwrap import dedent


_shared_docs = dict(**_shared_docs)
_doc_template = """

Returns
-------
same type as input

See also
--------
pandas.Series.%(name)s
pandas.DataFrame.%(name)s
"""


class _Window(PandasObject, SelectionMixin):
    _attributes = ['window', 'min_periods', 'center', 'win_type',
                   'axis', 'on', 'closed']
    exclusions = set()

    def __init__(self, obj, window=None, min_periods=None,
                 center=False, win_type=None, axis=0, on=None, closed=None,
                 **kwargs):

        self.__dict__.update(kwargs)
        self.blocks = []
        self.obj = obj
        self.on = on
        self.closed = closed
        self.window = window
        self.min_periods = min_periods
        self.center = center
        self.win_type = win_type
        self.win_freq = None
        self.axis = obj._get_axis_number(axis) if axis is not None else None
        self.validate()

    @property
    def _constructor(self):
        return Window

    @property
    def is_datetimelike(self):
        return None

    @property
    def _on(self):
        return None

    @property
    def is_freq_type(self):
        return self.win_type == 'freq'

    def validate(self):
        if self.center is not None and not is_bool(self.center):
            raise ValueError("center must be a boolean")
        if self.min_periods is not None and not \
           is_integer(self.min_periods):
            raise ValueError("min_periods must be an integer")
        if self.closed is not None and self.closed not in \
           ['right', 'both', 'left', 'neither']:
            raise ValueError("closed must be 'right', 'left', 'both' or "
                             "'neither'")

    def _convert_freq(self):
        """ resample according to the how, return a new object """

        obj = self._selected_obj
        index = None
        return obj, index

    def _create_blocks(self):
        """ split data into blocks & return conformed data """

        obj, index = self._convert_freq()
        if index is not None:
            index = self._on

        # filter out the on from the object
        if self.on is not None:
            if obj.ndim == 2:
                obj = obj.reindex(columns=obj.columns.difference([self.on]),
                                  copy=False)
        blocks = obj._to_dict_of_blocks(copy=False).values()

        return blocks, obj, index

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

        # create a new object to prevent aliasing
        if subset is None:
            subset = self.obj
        self = self._shallow_copy(subset)
        self._reset_cache()
        if subset.ndim == 2:
            if is_scalar(key) and key in subset or is_list_like(key):
                self._selection = key
        return self

    def __getattr__(self, attr):
        if attr in self._internal_names_set:
            return object.__getattribute__(self, attr)
        if attr in self.obj:
            return self[attr]

        raise AttributeError("%r object has no attribute %r" %
                             (type(self).__name__, attr))

    def _dir_additions(self):
        return self.obj._dir_additions()

    def _get_window(self, other=None):
        return self.window

    @property
    def _window_type(self):
        return self.__class__.__name__

    def __unicode__(self):
        """ provide a nice str repr of our rolling object """

        attrs = ["{k}={v}".format(k=k, v=getattr(self, k))
                 for k in self._attributes
                 if getattr(self, k, None) is not None]
        return "{klass} [{attrs}]".format(klass=self._window_type,
                                          attrs=','.join(attrs))

    def __iter__(self):
        url = 'https://github.com/pandas-dev/pandas/issues/11704'
        raise NotImplementedError('See issue #11704 {url}'.format(url=url))

    def _get_index(self, index=None):
        """
        Return index as ndarrays

        Returns
        -------
        tuple of (index, index_as_ndarray)
        """

        if self.is_freq_type:
            if index is None:
                index = self._on
            return index, index.asi8
        return index, index

    def _prep_values(self, values=None, kill_inf=True):

        if values is None:
            values = getattr(self._selected_obj, 'values', self._selected_obj)

        # GH #12373 : rolling functions error on float32 data
        # make sure the data is coerced to float64
        if is_float_dtype(values.dtype):
            values = _ensure_float64(values)
        elif is_integer_dtype(values.dtype):
            values = _ensure_float64(values)
        elif needs_i8_conversion(values.dtype):
            raise NotImplementedError("ops for {action} for this "
                                      "dtype {dtype} are not "
                                      "implemented".format(
                                          action=self._window_type,
                                          dtype=values.dtype))
        else:
            try:
                values = _ensure_float64(values)
            except (ValueError, TypeError):
                raise TypeError("cannot handle this type -> {0}"
                                "".format(values.dtype))

        if kill_inf:
            values = values.copy()
            values[np.isinf(values)] = np.NaN

        return values

    def _wrap_result(self, result, block=None, obj=None):
        """ wrap a single result """

        if obj is None:
            obj = self._selected_obj
        index = obj.index

        if isinstance(result, np.ndarray):

            # coerce if necessary
            if block is not None:
                if is_timedelta64_dtype(block.values.dtype):
                    from pandas import to_timedelta
                    result = to_timedelta(
                        result.ravel(), unit='ns').values.reshape(result.shape)

            if result.ndim == 1:
                from pandas import Series
                return Series(result, index, name=obj.name)

            return type(obj)(result, index=index, columns=block.columns)
        return result

    def _wrap_results(self, results, blocks, obj):
        """
        wrap the results

        Parameters
        ----------
        results : list of ndarrays
        blocks : list of blocks
        obj : conformed data (may be resampled)
        """

        from pandas import Series, concat
        from pandas.core.index import _ensure_index

        final = []
        for result, block in zip(results, blocks):

            result = self._wrap_result(result, block=block, obj=obj)
            if result.ndim == 1:
                return result
            final.append(result)

        # if we have an 'on' column
        # we want to put it back into the results
        # in the same location
        columns = self._selected_obj.columns
        if self.on is not None and not self._on.equals(obj.index):

            name = self._on.name
            final.append(Series(self._on, index=obj.index, name=name))

            if self._selection is not None:

                selection = _ensure_index(self._selection)

                # need to reorder to include original location of
                # the on column (if its not already there)
                if name not in selection:
                    columns = self.obj.columns
                    indexer = columns.get_indexer(selection.tolist() + [name])
                    columns = columns.take(sorted(indexer))

        if not len(final):
            return obj.astype('float64')
        return concat(final, axis=1).reindex(columns=columns, copy=False)

    def _center_window(self, result, window):
        """ center the result in the window """
        if self.axis > result.ndim - 1:
            raise ValueError("Requested axis is larger then no. of argument "
                             "dimensions")

        offset = _offset(window, True)
        if offset > 0:
            if isinstance(result, (ABCSeries, ABCDataFrame)):
                result = result.slice_shift(-offset, axis=self.axis)
            else:
                lead_indexer = [slice(None)] * result.ndim
                lead_indexer[self.axis] = slice(offset, None)
                result = np.copy(result[tuple(lead_indexer)])
        return result

    def aggregate(self, arg, *args, **kwargs):
        result, how = self._aggregate(arg, *args, **kwargs)
        if result is None:
            return self.apply(arg, raw=False, args=args, kwargs=kwargs)
        return result

    agg = aggregate

    _shared_docs['sum'] = dedent("""
    Calculate %(name)s sum of given DataFrame or Series.

    Parameters
    ----------
    *args, **kwargs
        For compatibility with other %(name)s methods. Has no effect
        on the computed value.

    Returns
    -------
    Series or DataFrame
        Same type as the input, with the same index, containing the
        %(name)s sum.

    See Also
    --------
    Series.sum : Reducing sum for Series.
    DataFrame.sum : Reducing sum for DataFrame.

    Examples
    --------
    >>> s = pd.Series([1, 2, 3, 4, 5])
    >>> s
    0    1
    1    2
    2    3
    3    4
    4    5
    dtype: int64

    >>> s.rolling(3).sum()
    0     NaN
    1     NaN
    2     6.0
    3     9.0
    4    12.0
    dtype: float64

    >>> s.expanding(3).sum()
    0     NaN
    1     NaN
    2     6.0
    3    10.0
    4    15.0
    dtype: float64

    >>> s.rolling(3, center=True).sum()
    0     NaN
    1     6.0
    2     9.0
    3    12.0
    4     NaN
    dtype: float64

    For DataFrame, each %(name)s sum is computed column-wise.

    >>> df = pd.DataFrame({"A": s, "B": s ** 2})
    >>> df
       A   B
    0  1   1
    1  2   4
    2  3   9
    3  4  16
    4  5  25

    >>> df.rolling(3).sum()
          A     B
    0   NaN   NaN
    1   NaN   NaN
    2   6.0  14.0
    3   9.0  29.0
    4  12.0  50.0
    """)

    _shared_docs['mean'] = dedent("""
    Calculate the %(name)s mean of the values.

    Parameters
    ----------
    *args
        Under Review.
    **kwargs
        Under Review.

    Returns
    -------
    Series or DataFrame
        Returned object type is determined by the caller of the %(name)s
        calculation.

    See Also
    --------
    Series.%(name)s : Calling object with Series data
    DataFrame.%(name)s : Calling object with DataFrames
    Series.mean : Equivalent method for Series
    DataFrame.mean : Equivalent method for DataFrame

    Examples
    --------
    The below examples will show rolling mean calculations with window sizes of
    two and three, respectively.

    >>> s = pd.Series([1, 2, 3, 4])
    >>> s.rolling(2).mean()
    0    NaN
    1    1.5
    2    2.5
    3    3.5
    dtype: float64

    >>> s.rolling(3).mean()
    0    NaN
    1    NaN
    2    2.0
    3    3.0
    dtype: float64
    """)


class Window(_Window):
    """
    Provides rolling window calculations.

    .. versionadded:: 0.18.0

    Parameters
    ----------
    window : int, or offset
        Size of the moving window. This is the number of observations used for
        calculating the statistic. Each window will be a fixed size.

        If its an offset then this will be the time period of each window. Each
        window will be a variable sized based on the observations included in
        the time-period. This is only valid for datetimelike indexes. This is
        new in 0.19.0
    min_periods : int, default None
        Minimum number of observations in window required to have a value
        (otherwise result is NA). For a window that is specified by an offset,
        this will default to 1.
    center : boolean, default False
        Set the labels at the center of the window.
    win_type : string, default None
        Provide a window type. If ``None``, all points are evenly weighted.
        See the notes below for further information.
    on : string, optional
        For a DataFrame, column on which to calculate
        the rolling window, rather than the index
    closed : string, default None
        Make the interval closed on the 'right', 'left', 'both' or
        'neither' endpoints.
        For offset-based windows, it defaults to 'right'.
        For fixed windows, defaults to 'both'. Remaining cases not implemented
        for fixed windows.

        .. versionadded:: 0.20.0

    axis : int or string, default 0

    Returns
    -------
    a Window or Rolling sub-classed for the particular operation

    Examples
    --------

    >>> df = pd.DataFrame({'B': [0, 1, 2, np.nan, 4]})
    >>> df
         B
    0  0.0
    1  1.0
    2  2.0
    3  NaN
    4  4.0

    Rolling sum with a window length of 2, using the 'triang'
    window type.

    >>> df.rolling(2, win_type='triang').sum()
         B
    0  NaN
    1  1.0
    2  2.5
    3  NaN
    4  NaN

    Rolling sum with a window length of 2, min_periods defaults
    to the window length.

    >>> df.rolling(2).sum()
         B
    0  NaN
    1  1.0
    2  3.0
    3  NaN
    4  NaN

    Same as above, but explicitly set the min_periods

    >>> df.rolling(2, min_periods=1).sum()
         B
    0  0.0
    1  1.0
    2  3.0
    3  2.0
    4  4.0

    A ragged (meaning not-a-regular frequency), time-indexed DataFrame

    >>> df = pd.DataFrame({'B': [0, 1, 2, np.nan, 4]},
    ...                   index = [pd.Timestamp('20130101 09:00:00'),
    ...                            pd.Timestamp('20130101 09:00:02'),
    ...                            pd.Timestamp('20130101 09:00:03'),
    ...                            pd.Timestamp('20130101 09:00:05'),
    ...                            pd.Timestamp('20130101 09:00:06')])

    >>> df
                           B
    2013-01-01 09:00:00  0.0
    2013-01-01 09:00:02  1.0
    2013-01-01 09:00:03  2.0
    2013-01-01 09:00:05  NaN
    2013-01-01 09:00:06  4.0


    Contrasting to an integer rolling window, this will roll a variable
    length window corresponding to the time period.
    The default for min_periods is 1.

    >>> df.rolling('2s').sum()
                           B
    2013-01-01 09:00:00  0.0
    2013-01-01 09:00:02  1.0
    2013-01-01 09:00:03  3.0
    2013-01-01 09:00:05  NaN
    2013-01-01 09:00:06  4.0

    Notes
    -----
    By default, the result is set to the right edge of the window. This can be
    changed to the center of the window by setting ``center=True``.

    To learn more about the offsets & frequency strings, please see `this link
    <http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases>`__.

    The recognized win_types are:

    * ``boxcar``
    * ``triang``
    * ``blackman``
    * ``hamming``
    * ``bartlett``
    * ``parzen``
    * ``bohman``
    * ``blackmanharris``
    * ``nuttall``
    * ``barthann``
    * ``kaiser`` (needs beta)
    * ``gaussian`` (needs std)
    * ``general_gaussian`` (needs power, width)
    * ``slepian`` (needs width).

    If ``win_type=None`` all points are evenly weighted. To learn more about
    different window types see `scipy.signal window functions
    <https://docs.scipy.org/doc/scipy/reference/signal.html#window-functions>`__.

    See Also
    --------
    expanding : Provides expanding transformations.
    ewm : Provides exponential weighted functions
    """

    def validate(self):
        super(Window, self).validate()

        window = self.window
        if isinstance(window, (list, tuple, np.ndarray)):
            pass
        elif is_integer(window):
            if window < 0:
                raise ValueError("window must be non-negative")
            try:
                import scipy.signal as sig
            except ImportError:
                raise ImportError('Please install scipy to generate window '
                                  'weight')

            if not isinstance(self.win_type, compat.string_types):
                raise ValueError('Invalid win_type {0}'.format(self.win_type))
            if getattr(sig, self.win_type, None) is None:
                raise ValueError('Invalid win_type {0}'.format(self.win_type))
        else:
            raise ValueError('Invalid window {0}'.format(window))

    def _prep_window(self, **kwargs):
        """
        provide validation for our window type, return the window
        we have already been validated
        """

        window = self._get_window()
        if isinstance(window, (list, tuple, np.ndarray)):
            return com._asarray_tuplesafe(window).astype(float)
        elif is_integer(window):
            import scipy.signal as sig

            # the below may pop from kwargs
            def _validate_win_type(win_type, kwargs):
                arg_map = {'kaiser': ['beta'],
                           'gaussian': ['std'],
                           'general_gaussian': ['power', 'width'],
                           'slepian': ['width']}
                if win_type in arg_map:
                    return tuple([win_type] + _pop_args(win_type,
                                                        arg_map[win_type],
                                                        kwargs))
                return win_type

            def _pop_args(win_type, arg_names, kwargs):
                msg = '%s window requires %%s' % win_type
                all_args = []
                for n in arg_names:
                    if n not in kwargs:
                        raise ValueError(msg % n)
                    all_args.append(kwargs.pop(n))
                return all_args

            win_type = _validate_win_type(self.win_type, kwargs)
            # GH #15662. `False` makes symmetric window, rather than periodic.
            return sig.get_window(win_type, window, False).astype(float)

    def _apply_window(self, mean=True, **kwargs):
        """
        Applies a moving window of type ``window_type`` on the data.

        Parameters
        ----------
        mean : boolean, default True
            If True computes weighted mean, else weighted sum

        Returns
        -------
        y : type of input argument

        """
        window = self._prep_window(**kwargs)
        center = self.center

        blocks, obj, index = self._create_blocks()
        results = []
        for b in blocks:
            try:
                values = self._prep_values(b.values)
            except TypeError:
                results.append(b.values.copy())
                continue

            if values.size == 0:
                results.append(values.copy())
                continue

            offset = _offset(window, center)
            additional_nans = np.array([np.NaN] * offset)

            def f(arg, *args, **kwargs):
                minp = _use_window(self.min_periods, len(window))
                return _window.roll_window(np.concatenate((arg,
                                                           additional_nans))
                                           if center else arg, window, minp,
                                           avg=mean)

            result = np.apply_along_axis(f, self.axis, values)

            if center:
                result = self._center_window(result, window)
            results.append(result)

        return self._wrap_results(results, blocks, obj)

    _agg_doc = dedent("""
    Examples
    --------

    >>> df = pd.DataFrame(np.random.randn(10, 3), columns=['A', 'B', 'C'])
    >>> df
              A         B         C
    0 -2.385977 -0.102758  0.438822
    1 -1.004295  0.905829 -0.954544
    2  0.735167 -0.165272 -1.619346
    3 -0.702657 -1.340923 -0.706334
    4 -0.246845  0.211596 -0.901819
    5  2.463718  3.157577 -1.380906
    6 -1.142255  2.340594 -0.039875
    7  1.396598 -1.647453  1.677227
    8 -0.543425  1.761277 -0.220481
    9 -0.640505  0.289374 -1.550670

    >>> df.rolling(3, win_type='boxcar').agg('mean')
              A         B         C
    0       NaN       NaN       NaN
    1       NaN       NaN       NaN
    2 -0.885035  0.212600 -0.711689
    3 -0.323928 -0.200122 -1.093408
    4 -0.071445 -0.431533 -1.075833
    5  0.504739  0.676083 -0.996353
    6  0.358206  1.903256 -0.774200
    7  0.906020  1.283573  0.085482
    8 -0.096361  0.818139  0.472290
    9  0.070889  0.134399 -0.031308

    See also
    --------
    pandas.DataFrame.rolling.aggregate
    pandas.DataFrame.aggregate

    """)

    @Appender(_agg_doc)
    @Appender(_shared_docs['aggregate'] % dict(
        versionadded='',
        klass='Series/DataFrame',
        axis=''))
    def aggregate(self, arg, *args, **kwargs):
        result, how = self._aggregate(arg, *args, **kwargs)
        if result is None:

            # these must apply directly
            result = arg(self)

        return result

    agg = aggregate

    @Substitution(name='window')
    @Appender(_shared_docs['sum'])
    def sum(self, *args, **kwargs):
        nv.validate_window_func('sum', args, kwargs)
        return self._apply_window(mean=False, **kwargs)

    @Substitution(name='window')
    @Appender(_shared_docs['mean'])
    def mean(self, *args, **kwargs):
        nv.validate_window_func('mean', args, kwargs)
        return self._apply_window(mean=True, **kwargs)


class _GroupByMixin(GroupByMixin):
    """ provide the groupby facilities """

    def __init__(self, obj, *args, **kwargs):
        parent = kwargs.pop('parent', None)  # noqa
        groupby = kwargs.pop('groupby', None)
        if groupby is None:
            groupby, obj = obj, obj.obj
        self._groupby = groupby
        self._groupby.mutated = True
        self._groupby.grouper.mutated = True
        super(GroupByMixin, self).__init__(obj, *args, **kwargs)

    count = GroupByMixin._dispatch('count')
    corr = GroupByMixin._dispatch('corr', other=None, pairwise=None)
    cov = GroupByMixin._dispatch('cov', other=None, pairwise=None)

    def _apply(self, func, name, window=None, center=None,
               check_minp=None, **kwargs):
        """
        dispatch to apply; we are stripping all of the _apply kwargs and
        performing the original function call on the grouped object
        """

        def f(x, name=name, *args):
            x = self._shallow_copy(x)

            if isinstance(name, compat.string_types):
                return getattr(x, name)(*args, **kwargs)

            return x.apply(name, *args, **kwargs)

        return self._groupby.apply(f)


class _Rolling(_Window):

    @property
    def _constructor(self):
        return Rolling

    def _apply(self, func, name=None, window=None, center=None,
               check_minp=None, **kwargs):
        """
        Rolling statistical measure using supplied function. Designed to be
        used with passed-in Cython array-based functions.

        Parameters
        ----------
        func : string/callable to apply
        name : string, optional
           name of this function
        window : int/array, default to _get_window()
        center : boolean, default to self.center
        check_minp : function, default to _use_window

        Returns
        -------
        y : type of input
        """
        if center is None:
            center = self.center
        if window is None:
            window = self._get_window()

        if check_minp is None:
            check_minp = _use_window

        blocks, obj, index = self._create_blocks()
        index, indexi = self._get_index(index=index)
        results = []
        for b in blocks:
            values = self._prep_values(b.values)

            if values.size == 0:
                results.append(values.copy())
                continue

            # if we have a string function name, wrap it
            if isinstance(func, compat.string_types):
                cfunc = getattr(_window, func, None)
                if cfunc is None:
                    raise ValueError("we do not support this function "
                                     "in _window.{0}".format(func))

                def func(arg, window, min_periods=None, closed=None):
                    minp = check_minp(min_periods, window)
                    # ensure we are only rolling on floats
                    arg = _ensure_float64(arg)
                    return cfunc(arg,
                                 window, minp, indexi, closed, **kwargs)

            # calculation function
            if center:
                offset = _offset(window, center)
                additional_nans = np.array([np.NaN] * offset)

                def calc(x):
                    return func(np.concatenate((x, additional_nans)),
                                window, min_periods=self.min_periods,
                                closed=self.closed)
            else:

                def calc(x):
                    return func(x, window, min_periods=self.min_periods,
                                closed=self.closed)

            with np.errstate(all='ignore'):
                if values.ndim > 1:
                    result = np.apply_along_axis(calc, self.axis, values)
                else:
                    result = calc(values)

            if center:
                result = self._center_window(result, window)

            results.append(result)

        return self._wrap_results(results, blocks, obj)


class _Rolling_and_Expanding(_Rolling):

    _shared_docs['count'] = dedent(r"""
    The %(name)s count of any non-NaN observations inside the window.

    Returns
    -------
    Series or DataFrame
        Returned object type is determined by the caller of the %(name)s
        calculation.

    See Also
    --------
    pandas.Series.%(name)s : Calling object with Series data
    pandas.DataFrame.%(name)s : Calling object with DataFrames
    pandas.DataFrame.count : Count of the full DataFrame

    Examples
    --------
    >>> s = pd.Series([2, 3, np.nan, 10])
    >>> s.rolling(2).count()
    0    1.0
    1    2.0
    2    1.0
    3    1.0
    dtype: float64
    >>> s.rolling(3).count()
    0    1.0
    1    2.0
    2    2.0
    3    2.0
    dtype: float64
    >>> s.rolling(4).count()
    0    1.0
    1    2.0
    2    2.0
    3    3.0
    dtype: float64
    """)

    def count(self):

        blocks, obj, index = self._create_blocks()
        index, indexi = self._get_index(index=index)

        window = self._get_window()
        window = min(window, len(obj)) if not self.center else window

        results = []
        for b in blocks:
            result = b.notna().astype(int)
            result = self._constructor(result, window=window, min_periods=0,
                                       center=self.center,
                                       closed=self.closed).sum()
            results.append(result)

        return self._wrap_results(results, blocks, obj)

    _shared_docs['apply'] = dedent(r"""
    %(name)s function apply

    Parameters
    ----------
    func : function
        Must produce a single value from an ndarray input if ``raw=True``
        or a Series if ``raw=False``
    raw : bool, default None
        * ``False`` : passes each row or column as a Series to the
          function.
        * ``True`` or ``None`` : the passed function will receive ndarray
          objects instead.
          If you are just applying a NumPy reduction function this will
          achieve much better performance.

        The `raw` parameter is required and will show a FutureWarning if
        not passed. In the future `raw` will default to False.

        .. versionadded:: 0.23.0

    \*args and \*\*kwargs are passed to the function""")

    def apply(self, func, raw=None, args=(), kwargs={}):
        from pandas import Series

        # TODO: _level is unused?
        _level = kwargs.pop('_level', None)  # noqa
        window = self._get_window()
        offset = _offset(window, self.center)
        index, indexi = self._get_index()

        # TODO: default is for backward compat
        # change to False in the future
        if raw is None:
            warnings.warn(
                "Currently, 'apply' passes the values as ndarrays to the "
                "applied function. In the future, this will change to passing "
                "it as Series objects. You need to specify 'raw=True' to keep "
                "the current behaviour, and you can pass 'raw=False' to "
                "silence this warning", FutureWarning, stacklevel=3)
            raw = True

        def f(arg, window, min_periods, closed):
            minp = _use_window(min_periods, window)
            if not raw:
                arg = Series(arg, index=self.obj.index)
            return _window.roll_generic(
                arg, window, minp, indexi,
                closed, offset, func, raw, args, kwargs)

        return self._apply(f, func, args=args, kwargs=kwargs,
                           center=False, raw=raw)

    def sum(self, *args, **kwargs):
        nv.validate_window_func('sum', args, kwargs)
        return self._apply('roll_sum', 'sum', **kwargs)

    _shared_docs['max'] = dedent("""
    %(name)s maximum
    """)

    def max(self, *args, **kwargs):
        nv.validate_window_func('max', args, kwargs)
        return self._apply('roll_max', 'max', **kwargs)

    _shared_docs['min'] = dedent("""
    Calculate the %(name)s minimum.

    Parameters
    ----------
    **kwargs
        Under Review.

    Returns
    -------
    Series or DataFrame
        Returned object type is determined by the caller of the %(name)s
        calculation.

    See Also
    --------
    Series.%(name)s : Calling object with a Series
    DataFrame.%(name)s : Calling object with a DataFrame
    Series.min : Similar method for Series
    DataFrame.min : Similar method for DataFrame

    Examples
    --------
    Performing a rolling minimum with a window size of 3.

    >>> s = pd.Series([4, 3, 5, 2, 6])
    >>> s.rolling(3).min()
    0    NaN
    1    NaN
    2    3.0
    3    2.0
    4    2.0
    dtype: float64
    """)

    def min(self, *args, **kwargs):
        nv.validate_window_func('min', args, kwargs)
        return self._apply('roll_min', 'min', **kwargs)

    def mean(self, *args, **kwargs):
        nv.validate_window_func('mean', args, kwargs)
        return self._apply('roll_mean', 'mean', **kwargs)

    _shared_docs['median'] = dedent("""
    Calculate the %(name)s median.

    Parameters
    ----------
    **kwargs
        For compatibility with other %(name)s methods. Has no effect
        on the computed median.

    Returns
    -------
    Series or DataFrame
        Returned type is the same as the original object.

    See Also
    --------
    Series.%(name)s : Calling object with Series data
    DataFrame.%(name)s : Calling object with DataFrames
    Series.median : Equivalent method for Series
    DataFrame.median : Equivalent method for DataFrame

    Examples
    --------
    Compute the rolling median of a series with a window size of 3.

    >>> s = pd.Series([0, 1, 2, 3, 4])
    >>> s.rolling(3).median()
    0    NaN
    1    NaN
    2    1.0
    3    2.0
    4    3.0
    dtype: float64
    """)

    def median(self, **kwargs):
        return self._apply('roll_median_c', 'median', **kwargs)

    _shared_docs['std'] = dedent("""
    Calculate %(name)s standard deviation.

    Normalized by N-1 by default. This can be changed using the `ddof`
    argument.

    Parameters
    ----------
    ddof : int, default 1
        Delta Degrees of Freedom.  The divisor used in calculations
        is ``N - ddof``, where ``N`` represents the number of elements.
    *args, **kwargs
        For NumPy compatibility. No additional arguments are used.

    Returns
    -------
    Series or DataFrame
        Returns the same object type as the caller of the %(name)s calculation.

    See Also
    --------
    Series.%(name)s : Calling object with Series data
    DataFrame.%(name)s : Calling object with DataFrames
    Series.std : Equivalent method for Series
    DataFrame.std : Equivalent method for DataFrame
    numpy.std : Equivalent method for Numpy array

    Notes
    -----
    The default `ddof` of 1 used in Series.std is different than the default
    `ddof` of 0 in numpy.std.

    A minimum of one period is required for the rolling calculation.

    Examples
    --------
    >>> s = pd.Series([5, 5, 6, 7, 5, 5, 5])
    >>> s.rolling(3).std()
    0         NaN
    1         NaN
    2    0.577350
    3    1.000000
    4    1.000000
    5    1.154701
    6    0.000000
    dtype: float64

    >>> s.expanding(3).std()
    0         NaN
    1         NaN
    2    0.577350
    3    0.957427
    4    0.894427
    5    0.836660
    6    0.786796
    dtype: float64
    """)

    def std(self, ddof=1, *args, **kwargs):
        nv.validate_window_func('std', args, kwargs)
        window = self._get_window()
        index, indexi = self._get_index()

        def f(arg, *args, **kwargs):
            minp = _require_min_periods(1)(self.min_periods, window)
            return _zsqrt(_window.roll_var(arg, window, minp, indexi,
                                           self.closed, ddof))

        return self._apply(f, 'std', check_minp=_require_min_periods(1),
                           ddof=ddof, **kwargs)

    _shared_docs['var'] = dedent("""
    Calculate unbiased %(name)s variance.

    Normalized by N-1 by default. This can be changed using the `ddof`
    argument.

    Parameters
    ----------
    ddof : int, default 1
        Delta Degrees of Freedom.  The divisor used in calculations
        is ``N - ddof``, where ``N`` represents the number of elements.
    *args, **kwargs
        For NumPy compatibility. No additional arguments are used.

    Returns
    -------
    Series or DataFrame
        Returns the same object type as the caller of the %(name)s calculation.

    See Also
    --------
    Series.%(name)s : Calling object with Series data
    DataFrame.%(name)s : Calling object with DataFrames
    Series.var : Equivalent method for Series
    DataFrame.var : Equivalent method for DataFrame
    numpy.var : Equivalent method for Numpy array

    Notes
    -----
    The default `ddof` of 1 used in :meth:`Series.var` is different than the
    default `ddof` of 0 in :func:`numpy.var`.

    A minimum of 1 period is required for the rolling calculation.

    Examples
    --------
    >>> s = pd.Series([5, 5, 6, 7, 5, 5, 5])
    >>> s.rolling(3).var()
    0         NaN
    1         NaN
    2    0.333333
    3    1.000000
    4    1.000000
    5    1.333333
    6    0.000000
    dtype: float64

    >>> s.expanding(3).var()
    0         NaN
    1         NaN
    2    0.333333
    3    0.916667
    4    0.800000
    5    0.700000
    6    0.619048
    dtype: float64
    """)

    def var(self, ddof=1, *args, **kwargs):
        nv.validate_window_func('var', args, kwargs)
        return self._apply('roll_var', 'var',
                           check_minp=_require_min_periods(1), ddof=ddof,
                           **kwargs)

    _shared_docs['skew'] = """Unbiased %(name)s skewness"""

    def skew(self, **kwargs):
        return self._apply('roll_skew', 'skew',
                           check_minp=_require_min_periods(3), **kwargs)

    _shared_docs['kurt'] = dedent("""
    Calculate unbiased %(name)s kurtosis.

    This function uses Fisher's definition of kurtosis without bias.

    Parameters
    ----------
    **kwargs
        Under Review.

    Returns
    -------
    Series or DataFrame
        Returned object type is determined by the caller of the %(name)s
        calculation

    See Also
    --------
    Series.%(name)s : Calling object with Series data
    DataFrame.%(name)s : Calling object with DataFrames
    Series.kurt : Equivalent method for Series
    DataFrame.kurt : Equivalent method for DataFrame
    scipy.stats.skew : Third moment of a probability density
    scipy.stats.kurtosis : Reference SciPy method

    Notes
    -----
    A minimum of 4 periods is required for the %(name)s calculation.
    """)

    def kurt(self, **kwargs):
        return self._apply('roll_kurt', 'kurt',
                           check_minp=_require_min_periods(4), **kwargs)

    _shared_docs['quantile'] = dedent("""
    %(name)s quantile.

    Parameters
    ----------
    quantile : float
        Quantile to compute. 0 <= quantile <= 1.
    interpolation : {'linear', 'lower', 'higher', 'midpoint', 'nearest'}
        .. versionadded:: 0.23.0

        This optional parameter specifies the interpolation method to use,
        when the desired quantile lies between two data points `i` and `j`:

            * linear: `i + (j - i) * fraction`, where `fraction` is the
              fractional part of the index surrounded by `i` and `j`.
            * lower: `i`.
            * higher: `j`.
            * nearest: `i` or `j` whichever is nearest.
            * midpoint: (`i` + `j`) / 2.
    **kwargs:
        For compatibility with other %(name)s methods. Has no effect on
        the result.

    Returns
    -------
    Series or DataFrame
        Returned object type is determined by the caller of the %(name)s
        calculation.

    Examples
    --------
    >>> s = pd.Series([1, 2, 3, 4])
    >>> s.rolling(2).quantile(.4, interpolation='lower')
    0    NaN
    1    1.0
    2    2.0
    3    3.0
    dtype: float64

    >>> s.rolling(2).quantile(.4, interpolation='midpoint')
    0    NaN
    1    1.5
    2    2.5
    3    3.5
    dtype: float64

    See Also
    --------
    pandas.Series.quantile : Computes value at the given quantile over all data
        in Series.
    pandas.DataFrame.quantile : Computes values at the given quantile over
        requested axis in DataFrame.
    """)

    def quantile(self, quantile, interpolation='linear', **kwargs):
        window = self._get_window()
        index, indexi = self._get_index()

        def f(arg, *args, **kwargs):
            minp = _use_window(self.min_periods, window)
            if quantile == 1.0:
                return _window.roll_max(arg, window, minp, indexi,
                                        self.closed)
            elif quantile == 0.0:
                return _window.roll_min(arg, window, minp, indexi,
                                        self.closed)
            else:
                return _window.roll_quantile(arg, window, minp, indexi,
                                             self.closed, quantile,
                                             interpolation)

        return self._apply(f, 'quantile', quantile=quantile,
                           **kwargs)

    _shared_docs['cov'] = dedent("""
    %(name)s sample covariance

    Parameters
    ----------
    other : Series, DataFrame, or ndarray, optional
        if not supplied then will default to self and produce pairwise output
    pairwise : bool, default None
        If False then only matching columns between self and other will be used
        and the output will be a DataFrame.
        If True then all pairwise combinations will be calculated and the
        output will be a MultiIndexed DataFrame in the case of DataFrame
        inputs. In the case of missing elements, only complete pairwise
        observations will be used.
    ddof : int, default 1
        Delta Degrees of Freedom.  The divisor used in calculations
        is ``N - ddof``, where ``N`` represents the number of elements.""")

    def cov(self, other=None, pairwise=None, ddof=1, **kwargs):
        if other is None:
            other = self._selected_obj
            # only default unset
            pairwise = True if pairwise is None else pairwise
        other = self._shallow_copy(other)

        # GH 16058: offset window
        if self.is_freq_type:
            window = self.win_freq
        else:
            window = self._get_window(other)

        def _get_cov(X, Y):
            # GH #12373 : rolling functions error on float32 data
            # to avoid potential overflow, cast the data to float64
            X = X.astype('float64')
            Y = Y.astype('float64')
            mean = lambda x: x.rolling(window, self.min_periods,
                                       center=self.center).mean(**kwargs)
            count = (X + Y).rolling(window=window,
                                    center=self.center).count(**kwargs)
            bias_adj = count / (count - ddof)
            return (mean(X * Y) - mean(X) * mean(Y)) * bias_adj

        return _flex_binary_moment(self._selected_obj, other._selected_obj,
                                   _get_cov, pairwise=bool(pairwise))

    _shared_docs['corr'] = dedent("""
    %(name)s sample correlation

    Parameters
    ----------
    other : Series, DataFrame, or ndarray, optional
        if not supplied then will default to self and produce pairwise output
    pairwise : bool, default None
        If False then only matching columns between self and other will be
        used and the output will be a DataFrame.
        If True then all pairwise combinations will be calculated and the
        output will be a MultiIndex DataFrame in the case of DataFrame inputs.
        In the case of missing elements, only complete pairwise observations
        will be used.""")

    def corr(self, other=None, pairwise=None, **kwargs):
        if other is None:
            other = self._selected_obj
            # only default unset
            pairwise = True if pairwise is None else pairwise
        other = self._shallow_copy(other)
        window = self._get_window(other)

        def _get_corr(a, b):
            a = a.rolling(window=window, min_periods=self.min_periods,
                          center=self.center)
            b = b.rolling(window=window, min_periods=self.min_periods,
                          center=self.center)

            return a.cov(b, **kwargs) / (a.std(**kwargs) * b.std(**kwargs))

        return _flex_binary_moment(self._selected_obj, other._selected_obj,
                                   _get_corr, pairwise=bool(pairwise))


class Rolling(_Rolling_and_Expanding):

    @cache_readonly
    def is_datetimelike(self):
        return isinstance(self._on,
                          (ABCDatetimeIndex,
                           ABCTimedeltaIndex,
                           ABCPeriodIndex))

    @cache_readonly
    def _on(self):

        if self.on is None:
            return self.obj.index
        elif (isinstance(self.obj, ABCDataFrame) and
              self.on in self.obj.columns):
            from pandas import Index
            return Index(self.obj[self.on])
        else:
            raise ValueError("invalid on specified as {0}, "
                             "must be a column (if DataFrame) "
                             "or None".format(self.on))

    def validate(self):
        super(Rolling, self).validate()

        # we allow rolling on a datetimelike index
        if ((self.obj.empty or self.is_datetimelike) and
                isinstance(self.window, (compat.string_types, ABCDateOffset,
                                         timedelta))):

            self._validate_monotonic()
            freq = self._validate_freq()

            # we don't allow center
            if self.center:
                raise NotImplementedError("center is not implemented "
                                          "for datetimelike and offset "
                                          "based windows")

            # this will raise ValueError on non-fixed freqs
            self.win_freq = self.window
            self.window = freq.nanos
            self.win_type = 'freq'

            # min_periods must be an integer
            if self.min_periods is None:
                self.min_periods = 1

        elif not is_integer(self.window):
            raise ValueError("window must be an integer")
        elif self.window < 0:
            raise ValueError("window must be non-negative")

        if not self.is_datetimelike and self.closed is not None:
            raise ValueError("closed only implemented for datetimelike "
                             "and offset based windows")

    def _validate_monotonic(self):
        """ validate on is monotonic """
        if not self._on.is_monotonic:
            formatted = self.on or 'index'
            raise ValueError("{0} must be "
                             "monotonic".format(formatted))

    def _validate_freq(self):
        """ validate & return window frequency """
        from pandas.tseries.frequencies import to_offset
        try:
            return to_offset(self.window)
        except (TypeError, ValueError):
            raise ValueError("passed window {0} is not "
                             "compatible with a datetimelike "
                             "index".format(self.window))

    _agg_doc = dedent("""
    Examples
    --------

    >>> df = pd.DataFrame(np.random.randn(10, 3), columns=['A', 'B', 'C'])
    >>> df
              A         B         C
    0 -2.385977 -0.102758  0.438822
    1 -1.004295  0.905829 -0.954544
    2  0.735167 -0.165272 -1.619346
    3 -0.702657 -1.340923 -0.706334
    4 -0.246845  0.211596 -0.901819
    5  2.463718  3.157577 -1.380906
    6 -1.142255  2.340594 -0.039875
    7  1.396598 -1.647453  1.677227
    8 -0.543425  1.761277 -0.220481
    9 -0.640505  0.289374 -1.550670

    >>> df.rolling(3).sum()
              A         B         C
    0       NaN       NaN       NaN
    1       NaN       NaN       NaN
    2 -2.655105  0.637799 -2.135068
    3 -0.971785 -0.600366 -3.280224
    4 -0.214334 -1.294599 -3.227500
    5  1.514216  2.028250 -2.989060
    6  1.074618  5.709767 -2.322600
    7  2.718061  3.850718  0.256446
    8 -0.289082  2.454418  1.416871
    9  0.212668  0.403198 -0.093924


    >>> df.rolling(3).agg({'A':'sum', 'B':'min'})
              A         B
    0       NaN       NaN
    1       NaN       NaN
    2 -2.655105 -0.165272
    3 -0.971785 -1.340923
    4 -0.214334 -1.340923
    5  1.514216 -1.340923
    6  1.074618  0.211596
    7  2.718061 -1.647453
    8 -0.289082 -1.647453
    9  0.212668 -1.647453

    See also
    --------
    pandas.Series.rolling
    pandas.DataFrame.rolling

    """)

    @Appender(_agg_doc)
    @Appender(_shared_docs['aggregate'] % dict(
        versionadded='',
        klass='Series/DataFrame',
        axis=''))
    def aggregate(self, arg, *args, **kwargs):
        return super(Rolling, self).aggregate(arg, *args, **kwargs)

    agg = aggregate

    @Substitution(name='rolling')
    @Appender(_shared_docs['count'])
    def count(self):

        # different impl for freq counting
        if self.is_freq_type:
            return self._apply('roll_count', 'count')

        return super(Rolling, self).count()

    @Substitution(name='rolling')
    @Appender(_doc_template)
    @Appender(_shared_docs['apply'])
    def apply(self, func, raw=None, args=(), kwargs={}):
        return super(Rolling, self).apply(
            func, raw=raw, args=args, kwargs=kwargs)

    @Substitution(name='rolling')
    @Appender(_shared_docs['sum'])
    def sum(self, *args, **kwargs):
        nv.validate_rolling_func('sum', args, kwargs)
        return super(Rolling, self).sum(*args, **kwargs)

    @Substitution(name='rolling')
    @Appender(_doc_template)
    @Appender(_shared_docs['max'])
    def max(self, *args, **kwargs):
        nv.validate_rolling_func('max', args, kwargs)
        return super(Rolling, self).max(*args, **kwargs)

    @Substitution(name='rolling')
    @Appender(_shared_docs['min'])
    def min(self, *args, **kwargs):
        nv.validate_rolling_func('min', args, kwargs)
        return super(Rolling, self).min(*args, **kwargs)

    @Substitution(name='rolling')
    @Appender(_shared_docs['mean'])
    def mean(self, *args, **kwargs):
        nv.validate_rolling_func('mean', args, kwargs)
        return super(Rolling, self).mean(*args, **kwargs)

    @Substitution(name='rolling')
    @Appender(_shared_docs['median'])
    def median(self, **kwargs):
        return super(Rolling, self).median(**kwargs)

    @Substitution(name='rolling')
    @Appender(_shared_docs['std'])
    def std(self, ddof=1, *args, **kwargs):
        nv.validate_rolling_func('std', args, kwargs)
        return super(Rolling, self).std(ddof=ddof, **kwargs)

    @Substitution(name='rolling')
    @Appender(_shared_docs['var'])
    def var(self, ddof=1, *args, **kwargs):
        nv.validate_rolling_func('var', args, kwargs)
        return super(Rolling, self).var(ddof=ddof, **kwargs)

    @Substitution(name='rolling')
    @Appender(_doc_template)
    @Appender(_shared_docs['skew'])
    def skew(self, **kwargs):
        return super(Rolling, self).skew(**kwargs)

    _agg_doc = dedent("""
    Examples
    --------

    The example below will show a rolling calculation with a window size of
    four matching the equivalent function call using `scipy.stats`.

    >>> arr = [1, 2, 3, 4, 999]
    >>> fmt = "{0:.6f}"  # limit the printed precision to 6 digits
    >>> import scipy.stats
    >>> print(fmt.format(scipy.stats.kurtosis(arr[:-1], bias=False)))
    -1.200000
    >>> print(fmt.format(scipy.stats.kurtosis(arr[1:], bias=False)))
    3.999946
    >>> s = pd.Series(arr)
    >>> s.rolling(4).kurt()
    0         NaN
    1         NaN
    2         NaN
    3   -1.200000
    4    3.999946
    dtype: float64
    """)

    @Appender(_agg_doc)
    @Substitution(name='rolling')
    @Appender(_shared_docs['kurt'])
    def kurt(self, **kwargs):
        return super(Rolling, self).kurt(**kwargs)

    @Substitution(name='rolling')
    @Appender(_shared_docs['quantile'])
    def quantile(self, quantile, interpolation='linear', **kwargs):
        return super(Rolling, self).quantile(quantile=quantile,
                                             interpolation=interpolation,
                                             **kwargs)

    @Substitution(name='rolling')
    @Appender(_doc_template)
    @Appender(_shared_docs['cov'])
    def cov(self, other=None, pairwise=None, ddof=1, **kwargs):
        return super(Rolling, self).cov(other=other, pairwise=pairwise,
                                        ddof=ddof, **kwargs)

    @Substitution(name='rolling')
    @Appender(_doc_template)
    @Appender(_shared_docs['corr'])
    def corr(self, other=None, pairwise=None, **kwargs):
        return super(Rolling, self).corr(other=other, pairwise=pairwise,
                                         **kwargs)


class RollingGroupby(_GroupByMixin, Rolling):
    """
    Provides a rolling groupby implementation

    .. versionadded:: 0.18.1

    """
    @property
    def _constructor(self):
        return Rolling

    def _gotitem(self, key, ndim, subset=None):

        # we are setting the index on the actual object
        # here so our index is carried thru to the selected obj
        # when we do the splitting for the groupby
        if self.on is not None:
            self._groupby.obj = self._groupby.obj.set_index(self._on)
            self.on = None
        return super(RollingGroupby, self)._gotitem(key, ndim, subset=subset)

    def _validate_monotonic(self):
        """
        validate that on is monotonic;
        we don't care for groupby.rolling
        because we have already validated at a higher
        level
        """
        pass


class Expanding(_Rolling_and_Expanding):
    """
    Provides expanding transformations.

    .. versionadded:: 0.18.0

    Parameters
    ----------
    min_periods : int, default 1
        Minimum number of observations in window required to have a value
        (otherwise result is NA).
    center : boolean, default False
        Set the labels at the center of the window.
    axis : int or string, default 0

    Returns
    -------
    a Window sub-classed for the particular operation

    Examples
    --------

    >>> df = DataFrame({'B': [0, 1, 2, np.nan, 4]})
         B
    0  0.0
    1  1.0
    2  2.0
    3  NaN
    4  4.0

    >>> df.expanding(2).sum()
         B
    0  NaN
    1  1.0
    2  3.0
    3  3.0
    4  7.0

    Notes
    -----
    By default, the result is set to the right edge of the window. This can be
    changed to the center of the window by setting ``center=True``.

    See Also
    --------
    rolling : Provides rolling window calculations
    ewm : Provides exponential weighted functions
    """

    _attributes = ['min_periods', 'center', 'axis']

    def __init__(self, obj, min_periods=1, center=False, axis=0,
                 **kwargs):
        super(Expanding, self).__init__(obj=obj, min_periods=min_periods,
                                        center=center, axis=axis)

    @property
    def _constructor(self):
        return Expanding

    def _get_window(self, other=None):
        obj = self._selected_obj
        if other is None:
            return (max(len(obj), self.min_periods) if self.min_periods
                    else len(obj))
        return (max((len(obj) + len(obj)), self.min_periods)
                if self.min_periods else (len(obj) + len(obj)))

    _agg_doc = dedent("""
    Examples
    --------

    >>> df = pd.DataFrame(np.random.randn(10, 3), columns=['A', 'B', 'C'])
    >>> df
              A         B         C
    0 -2.385977 -0.102758  0.438822
    1 -1.004295  0.905829 -0.954544
    2  0.735167 -0.165272 -1.619346
    3 -0.702657 -1.340923 -0.706334
    4 -0.246845  0.211596 -0.901819
    5  2.463718  3.157577 -1.380906
    6 -1.142255  2.340594 -0.039875
    7  1.396598 -1.647453  1.677227
    8 -0.543425  1.761277 -0.220481
    9 -0.640505  0.289374 -1.550670

    >>> df.ewm(alpha=0.5).mean()
              A         B         C
    0 -2.385977 -0.102758  0.438822
    1 -1.464856  0.569633 -0.490089
    2 -0.207700  0.149687 -1.135379
    3 -0.471677 -0.645305 -0.906555
    4 -0.355635 -0.203033 -0.904111
    5  1.076417  1.503943 -1.146293
    6 -0.041654  1.925562 -0.588728
    7  0.680292  0.132049  0.548693
    8  0.067236  0.948257  0.163353
    9 -0.286980  0.618493 -0.694496

    See also
    --------
    pandas.DataFrame.expanding.aggregate
    pandas.DataFrame.rolling.aggregate
    pandas.DataFrame.aggregate

    """)

    @Appender(_agg_doc)
    @Appender(_shared_docs['aggregate'] % dict(
        versionadded='',
        klass='Series/DataFrame',
        axis=''))
    def aggregate(self, arg, *args, **kwargs):
        return super(Expanding, self).aggregate(arg, *args, **kwargs)

    agg = aggregate

    @Substitution(name='expanding')
    @Appender(_shared_docs['count'])
    def count(self, **kwargs):
        return super(Expanding, self).count(**kwargs)

    @Substitution(name='expanding')
    @Appender(_doc_template)
    @Appender(_shared_docs['apply'])
    def apply(self, func, raw=None, args=(), kwargs={}):
        return super(Expanding, self).apply(
            func, raw=raw, args=args, kwargs=kwargs)

    @Substitution(name='expanding')
    @Appender(_shared_docs['sum'])
    def sum(self, *args, **kwargs):
        nv.validate_expanding_func('sum', args, kwargs)
        return super(Expanding, self).sum(*args, **kwargs)

    @Substitution(name='expanding')
    @Appender(_doc_template)
    @Appender(_shared_docs['max'])
    def max(self, *args, **kwargs):
        nv.validate_expanding_func('max', args, kwargs)
        return super(Expanding, self).max(*args, **kwargs)

    @Substitution(name='expanding')
    @Appender(_shared_docs['min'])
    def min(self, *args, **kwargs):
        nv.validate_expanding_func('min', args, kwargs)
        return super(Expanding, self).min(*args, **kwargs)

    @Substitution(name='expanding')
    @Appender(_shared_docs['mean'])
    def mean(self, *args, **kwargs):
        nv.validate_expanding_func('mean', args, kwargs)
        return super(Expanding, self).mean(*args, **kwargs)

    @Substitution(name='expanding')
    @Appender(_shared_docs['median'])
    def median(self, **kwargs):
        return super(Expanding, self).median(**kwargs)

    @Substitution(name='expanding')
    @Appender(_shared_docs['std'])
    def std(self, ddof=1, *args, **kwargs):
        nv.validate_expanding_func('std', args, kwargs)
        return super(Expanding, self).std(ddof=ddof, **kwargs)

    @Substitution(name='expanding')
    @Appender(_shared_docs['var'])
    def var(self, ddof=1, *args, **kwargs):
        nv.validate_expanding_func('var', args, kwargs)
        return super(Expanding, self).var(ddof=ddof, **kwargs)

    @Substitution(name='expanding')
    @Appender(_doc_template)
    @Appender(_shared_docs['skew'])
    def skew(self, **kwargs):
        return super(Expanding, self).skew(**kwargs)

    _agg_doc = dedent("""
    Examples
    --------

    The example below will show an expanding calculation with a window size of
    four matching the equivalent function call using `scipy.stats`.

    >>> arr = [1, 2, 3, 4, 999]
    >>> import scipy.stats
    >>> fmt = "{0:.6f}"  # limit the printed precision to 6 digits
    >>> print(fmt.format(scipy.stats.kurtosis(arr[:-1], bias=False)))
    -1.200000
    >>> print(fmt.format(scipy.stats.kurtosis(arr, bias=False)))
    4.999874
    >>> s = pd.Series(arr)
    >>> s.expanding(4).kurt()
    0         NaN
    1         NaN
    2         NaN
    3   -1.200000
    4    4.999874
    dtype: float64
    """)

    @Appender(_agg_doc)
    @Substitution(name='expanding')
    @Appender(_shared_docs['kurt'])
    def kurt(self, **kwargs):
        return super(Expanding, self).kurt(**kwargs)

    @Substitution(name='expanding')
    @Appender(_shared_docs['quantile'])
    def quantile(self, quantile, interpolation='linear', **kwargs):
        return super(Expanding, self).quantile(quantile=quantile,
                                               interpolation=interpolation,
                                               **kwargs)

    @Substitution(name='expanding')
    @Appender(_doc_template)
    @Appender(_shared_docs['cov'])
    def cov(self, other=None, pairwise=None, ddof=1, **kwargs):
        return super(Expanding, self).cov(other=other, pairwise=pairwise,
                                          ddof=ddof, **kwargs)

    @Substitution(name='expanding')
    @Appender(_doc_template)
    @Appender(_shared_docs['corr'])
    def corr(self, other=None, pairwise=None, **kwargs):
        return super(Expanding, self).corr(other=other, pairwise=pairwise,
                                           **kwargs)


class ExpandingGroupby(_GroupByMixin, Expanding):
    """
    Provides a expanding groupby implementation

    .. versionadded:: 0.18.1

    """
    @property
    def _constructor(self):
        return Expanding


_bias_template = """

Parameters
----------
bias : boolean, default False
    Use a standard estimation bias correction
"""

_pairwise_template = """

Parameters
----------
other : Series, DataFrame, or ndarray, optional
    if not supplied then will default to self and produce pairwise output
pairwise : bool, default None
    If False then only matching columns between self and other will be used and
    the output will be a DataFrame.
    If True then all pairwise combinations will be calculated and the output
    will be a MultiIndex DataFrame in the case of DataFrame inputs.
    In the case of missing elements, only complete pairwise observations will
    be used.
bias : boolean, default False
   Use a standard estimation bias correction
"""


class EWM(_Rolling):
    r"""
    Provides exponential weighted functions

    .. versionadded:: 0.18.0

    Parameters
    ----------
    com : float, optional
        Specify decay in terms of center of mass,
        :math:`\alpha = 1 / (1 + com),\text{ for } com \geq 0`
    span : float, optional
        Specify decay in terms of span,
        :math:`\alpha = 2 / (span + 1),\text{ for } span \geq 1`
    halflife : float, optional
        Specify decay in terms of half-life,
        :math:`\alpha = 1 - exp(log(0.5) / halflife),\text{ for } halflife > 0`
    alpha : float, optional
        Specify smoothing factor :math:`\alpha` directly,
        :math:`0 < \alpha \leq 1`

        .. versionadded:: 0.18.0

    min_periods : int, default 0
        Minimum number of observations in window required to have a value
        (otherwise result is NA).
    adjust : boolean, default True
        Divide by decaying adjustment factor in beginning periods to account
        for imbalance in relative weightings (viewing EWMA as a moving average)
    ignore_na : boolean, default False
        Ignore missing values when calculating weights;
        specify True to reproduce pre-0.15.0 behavior

    Returns
    -------
    a Window sub-classed for the particular operation

    Examples
    --------

    >>> df = DataFrame({'B': [0, 1, 2, np.nan, 4]})
         B
    0  0.0
    1  1.0
    2  2.0
    3  NaN
    4  4.0

    >>> df.ewm(com=0.5).mean()
              B
    0  0.000000
    1  0.750000
    2  1.615385
    3  1.615385
    4  3.670213

    Notes
    -----
    Exactly one of center of mass, span, half-life, and alpha must be provided.
    Allowed values and relationship between the parameters are specified in the
    parameter descriptions above; see the link at the end of this section for
    a detailed explanation.

    When adjust is True (default), weighted averages are calculated using
    weights (1-alpha)**(n-1), (1-alpha)**(n-2), ..., 1-alpha, 1.

    When adjust is False, weighted averages are calculated recursively as:
       weighted_average[0] = arg[0];
       weighted_average[i] = (1-alpha)*weighted_average[i-1] + alpha*arg[i].

    When ignore_na is False (default), weights are based on absolute positions.
    For example, the weights of x and y used in calculating the final weighted
    average of [x, None, y] are (1-alpha)**2 and 1 (if adjust is True), and
    (1-alpha)**2 and alpha (if adjust is False).

    When ignore_na is True (reproducing pre-0.15.0 behavior), weights are based
    on relative positions. For example, the weights of x and y used in
    calculating the final weighted average of [x, None, y] are 1-alpha and 1
    (if adjust is True), and 1-alpha and alpha (if adjust is False).

    More details can be found at
    http://pandas.pydata.org/pandas-docs/stable/computation.html#exponentially-weighted-windows

    See Also
    --------
    rolling : Provides rolling window calculations
    expanding : Provides expanding transformations.
    """
    _attributes = ['com', 'min_periods', 'adjust', 'ignore_na', 'axis']

    def __init__(self, obj, com=None, span=None, halflife=None, alpha=None,
                 min_periods=0, adjust=True, ignore_na=False,
                 axis=0):
        self.obj = obj
        self.com = _get_center_of_mass(com, span, halflife, alpha)
        self.min_periods = min_periods
        self.adjust = adjust
        self.ignore_na = ignore_na
        self.axis = axis
        self.on = None

    @property
    def _constructor(self):
        return EWM

    _agg_doc = dedent("""
    Examples
    --------

    >>> df = pd.DataFrame(np.random.randn(10, 3), columns=['A', 'B', 'C'])
    >>> df
              A         B         C
    0 -2.385977 -0.102758  0.438822
    1 -1.004295  0.905829 -0.954544
    2  0.735167 -0.165272 -1.619346
    3 -0.702657 -1.340923 -0.706334
    4 -0.246845  0.211596 -0.901819
    5  2.463718  3.157577 -1.380906
    6 -1.142255  2.340594 -0.039875
    7  1.396598 -1.647453  1.677227
    8 -0.543425  1.761277 -0.220481
    9 -0.640505  0.289374 -1.550670

    >>> df.ewm(alpha=0.5).mean()
              A         B         C
    0 -2.385977 -0.102758  0.438822
    1 -1.464856  0.569633 -0.490089
    2 -0.207700  0.149687 -1.135379
    3 -0.471677 -0.645305 -0.906555
    4 -0.355635 -0.203033 -0.904111
    5  1.076417  1.503943 -1.146293
    6 -0.041654  1.925562 -0.588728
    7  0.680292  0.132049  0.548693
    8  0.067236  0.948257  0.163353
    9 -0.286980  0.618493 -0.694496

    See also
    --------
    pandas.DataFrame.rolling.aggregate

    """)

    @Appender(_agg_doc)
    @Appender(_shared_docs['aggregate'] % dict(
        versionadded='',
        klass='Series/DataFrame',
        axis=''))
    def aggregate(self, arg, *args, **kwargs):
        return super(EWM, self).aggregate(arg, *args, **kwargs)

    agg = aggregate

    def _apply(self, func, **kwargs):
        """Rolling statistical measure using supplied function. Designed to be
        used with passed-in Cython array-based functions.

        Parameters
        ----------
        func : string/callable to apply

        Returns
        -------
        y : type of input argument

        """
        blocks, obj, index = self._create_blocks()
        results = []
        for b in blocks:
            try:
                values = self._prep_values(b.values)
            except TypeError:
                results.append(b.values.copy())
                continue

            if values.size == 0:
                results.append(values.copy())
                continue

            # if we have a string function name, wrap it
            if isinstance(func, compat.string_types):
                cfunc = getattr(_window, func, None)
                if cfunc is None:
                    raise ValueError("we do not support this function "
                                     "in _window.{0}".format(func))

                def func(arg):
                    return cfunc(arg, self.com, int(self.adjust),
                                 int(self.ignore_na), int(self.min_periods))

            results.append(np.apply_along_axis(func, self.axis, values))

        return self._wrap_results(results, blocks, obj)

    @Substitution(name='ewm')
    @Appender(_doc_template)
    def mean(self, *args, **kwargs):
        """exponential weighted moving average"""
        nv.validate_window_func('mean', args, kwargs)
        return self._apply('ewma', **kwargs)

    @Substitution(name='ewm')
    @Appender(_doc_template)
    @Appender(_bias_template)
    def std(self, bias=False, *args, **kwargs):
        """exponential weighted moving stddev"""
        nv.validate_window_func('std', args, kwargs)
        return _zsqrt(self.var(bias=bias, **kwargs))

    vol = std

    @Substitution(name='ewm')
    @Appender(_doc_template)
    @Appender(_bias_template)
    def var(self, bias=False, *args, **kwargs):
        """exponential weighted moving variance"""
        nv.validate_window_func('var', args, kwargs)

        def f(arg):
            return _window.ewmcov(arg, arg, self.com, int(self.adjust),
                                  int(self.ignore_na), int(self.min_periods),
                                  int(bias))

        return self._apply(f, **kwargs)

    @Substitution(name='ewm')
    @Appender(_doc_template)
    @Appender(_pairwise_template)
    def cov(self, other=None, pairwise=None, bias=False, **kwargs):
        """exponential weighted sample covariance"""
        if other is None:
            other = self._selected_obj
            # only default unset
            pairwise = True if pairwise is None else pairwise
        other = self._shallow_copy(other)

        def _get_cov(X, Y):
            X = self._shallow_copy(X)
            Y = self._shallow_copy(Y)
            cov = _window.ewmcov(X._prep_values(), Y._prep_values(), self.com,
                                 int(self.adjust), int(self.ignore_na),
                                 int(self.min_periods), int(bias))
            return X._wrap_result(cov)

        return _flex_binary_moment(self._selected_obj, other._selected_obj,
                                   _get_cov, pairwise=bool(pairwise))

    @Substitution(name='ewm')
    @Appender(_doc_template)
    @Appender(_pairwise_template)
    def corr(self, other=None, pairwise=None, **kwargs):
        """exponential weighted sample correlation"""
        if other is None:
            other = self._selected_obj
            # only default unset
            pairwise = True if pairwise is None else pairwise
        other = self._shallow_copy(other)

        def _get_corr(X, Y):
            X = self._shallow_copy(X)
            Y = self._shallow_copy(Y)

            def _cov(x, y):
                return _window.ewmcov(x, y, self.com, int(self.adjust),
                                      int(self.ignore_na),
                                      int(self.min_periods),
                                      1)

            x_values = X._prep_values()
            y_values = Y._prep_values()
            with np.errstate(all='ignore'):
                cov = _cov(x_values, y_values)
                x_var = _cov(x_values, x_values)
                y_var = _cov(y_values, y_values)
                corr = cov / _zsqrt(x_var * y_var)
            return X._wrap_result(corr)

        return _flex_binary_moment(self._selected_obj, other._selected_obj,
                                   _get_corr, pairwise=bool(pairwise))

# Helper Funcs


def _flex_binary_moment(arg1, arg2, f, pairwise=False):

    if not (isinstance(arg1, (np.ndarray, ABCSeries, ABCDataFrame)) and
            isinstance(arg2, (np.ndarray, ABCSeries, ABCDataFrame))):
        raise TypeError("arguments to moment function must be of type "
                        "np.ndarray/Series/DataFrame")

    if (isinstance(arg1, (np.ndarray, ABCSeries)) and
            isinstance(arg2, (np.ndarray, ABCSeries))):
        X, Y = _prep_binary(arg1, arg2)
        return f(X, Y)

    elif isinstance(arg1, ABCDataFrame):
        from pandas import DataFrame

        def dataframe_from_int_dict(data, frame_template):
            result = DataFrame(data, index=frame_template.index)
            if len(result.columns) > 0:
                result.columns = frame_template.columns[result.columns]
            return result

        results = {}
        if isinstance(arg2, ABCDataFrame):
            if pairwise is False:
                if arg1 is arg2:
                    # special case in order to handle duplicate column names
                    for i, col in enumerate(arg1.columns):
                        results[i] = f(arg1.iloc[:, i], arg2.iloc[:, i])
                    return dataframe_from_int_dict(results, arg1)
                else:
                    if not arg1.columns.is_unique:
                        raise ValueError("'arg1' columns are not unique")
                    if not arg2.columns.is_unique:
                        raise ValueError("'arg2' columns are not unique")
                    with warnings.catch_warnings(record=True):
                        X, Y = arg1.align(arg2, join='outer')
                    X = X + 0 * Y
                    Y = Y + 0 * X

                    with warnings.catch_warnings(record=True):
                        res_columns = arg1.columns.union(arg2.columns)
                    for col in res_columns:
                        if col in X and col in Y:
                            results[col] = f(X[col], Y[col])
                    return DataFrame(results, index=X.index,
                                     columns=res_columns)
            elif pairwise is True:
                results = defaultdict(dict)
                for i, k1 in enumerate(arg1.columns):
                    for j, k2 in enumerate(arg2.columns):
                        if j < i and arg2 is arg1:
                            # Symmetric case
                            results[i][j] = results[j][i]
                        else:
                            results[i][j] = f(*_prep_binary(arg1.iloc[:, i],
                                                            arg2.iloc[:, j]))

                from pandas import MultiIndex, concat

                result_index = arg1.index.union(arg2.index)
                if len(result_index):

                    # construct result frame
                    result = concat(
                        [concat([results[i][j]
                                 for j, c in enumerate(arg2.columns)],
                                ignore_index=True)
                         for i, c in enumerate(arg1.columns)],
                        ignore_index=True,
                        axis=1)
                    result.columns = arg1.columns

                    # set the index and reorder
                    if arg2.columns.nlevels > 1:
                        result.index = MultiIndex.from_product(
                            arg2.columns.levels + [result_index])
                        result = result.reorder_levels([2, 0, 1]).sort_index()
                    else:
                        result.index = MultiIndex.from_product(
                            [range(len(arg2.columns)),
                             range(len(result_index))])
                        result = result.swaplevel(1, 0).sort_index()
                        result.index = MultiIndex.from_product(
                            [result_index] + [arg2.columns])
                else:

                    # empty result
                    result = DataFrame(
                        index=MultiIndex(levels=[arg1.index, arg2.columns],
                                         labels=[[], []]),
                        columns=arg2.columns,
                        dtype='float64')

                # reset our index names to arg1 names
                # reset our column names to arg2 names
                # careful not to mutate the original names
                result.columns = result.columns.set_names(
                    arg1.columns.names)
                result.index = result.index.set_names(
                    result_index.names + arg2.columns.names)

                return result

            else:
                raise ValueError("'pairwise' is not True/False")
        else:
            results = {}
            for i, col in enumerate(arg1.columns):
                results[i] = f(*_prep_binary(arg1.iloc[:, i], arg2))
            return dataframe_from_int_dict(results, arg1)

    else:
        return _flex_binary_moment(arg2, arg1, f)


def _get_center_of_mass(comass, span, halflife, alpha):
    valid_count = com._count_not_none(comass, span, halflife, alpha)
    if valid_count > 1:
        raise ValueError("comass, span, halflife, and alpha "
                         "are mutually exclusive")

    # Convert to center of mass; domain checks ensure 0 < alpha <= 1
    if comass is not None:
        if comass < 0:
            raise ValueError("comass must satisfy: comass >= 0")
    elif span is not None:
        if span < 1:
            raise ValueError("span must satisfy: span >= 1")
        comass = (span - 1) / 2.
    elif halflife is not None:
        if halflife <= 0:
            raise ValueError("halflife must satisfy: halflife > 0")
        decay = 1 - np.exp(np.log(0.5) / halflife)
        comass = 1 / decay - 1
    elif alpha is not None:
        if alpha <= 0 or alpha > 1:
            raise ValueError("alpha must satisfy: 0 < alpha <= 1")
        comass = (1.0 - alpha) / alpha
    else:
        raise ValueError("Must pass one of comass, span, halflife, or alpha")

    return float(comass)


def _offset(window, center):
    if not is_integer(window):
        window = len(window)
    offset = (window - 1) / 2. if center else 0
    try:
        return int(offset)
    except:
        return offset.astype(int)


def _require_min_periods(p):
    def _check_func(minp, window):
        if minp is None:
            return window
        else:
            return max(p, minp)

    return _check_func


def _use_window(minp, window):
    if minp is None:
        return window
    else:
        return minp


def _zsqrt(x):
    with np.errstate(all='ignore'):
        result = np.sqrt(x)
        mask = x < 0

    if isinstance(x, ABCDataFrame):
        if mask.values.any():
            result[mask] = 0
    else:
        if mask.any():
            result[mask] = 0

    return result


def _prep_binary(arg1, arg2):
    if not isinstance(arg2, type(arg1)):
        raise Exception('Input arrays must be of the same type!')

    # mask out values, this also makes a common index...
    X = arg1 + 0 * arg2
    Y = arg2 + 0 * arg1

    return X, Y


# Top-level exports


def rolling(obj, win_type=None, **kwds):
    if not isinstance(obj, (ABCSeries, ABCDataFrame)):
        raise TypeError('invalid type: %s' % type(obj))

    if win_type is not None:
        return Window(obj, win_type=win_type, **kwds)

    return Rolling(obj, **kwds)


rolling.__doc__ = Window.__doc__


def expanding(obj, **kwds):
    if not isinstance(obj, (ABCSeries, ABCDataFrame)):
        raise TypeError('invalid type: %s' % type(obj))

    return Expanding(obj, **kwds)


expanding.__doc__ = Expanding.__doc__


def ewm(obj, **kwds):
    if not isinstance(obj, (ABCSeries, ABCDataFrame)):
        raise TypeError('invalid type: %s' % type(obj))

    return EWM(obj, **kwds)


ewm.__doc__ = EWM.__doc__
