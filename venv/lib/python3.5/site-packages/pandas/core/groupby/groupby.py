import types
from functools import wraps, partial
import numpy as np
import datetime
import collections
import warnings
import copy
from textwrap import dedent
from contextlib import contextmanager

from pandas.compat import (
    zip, range, lzip,
    callable, map
)

from pandas import compat
from pandas.compat.numpy import function as nv
from pandas.compat import set_function_name

from pandas.core.dtypes.common import (
    is_numeric_dtype,
    is_timedelta64_dtype, is_datetime64_dtype,
    is_categorical_dtype,
    is_interval_dtype,
    is_datetimelike,
    is_datetime64_any_dtype,
    is_bool, is_integer_dtype,
    is_complex_dtype,
    is_bool_dtype,
    is_scalar,
    is_list_like,
    is_hashable,
    needs_i8_conversion,
    _ensure_float64,
    _ensure_platform_int,
    _ensure_int64,
    _ensure_object,
    _ensure_categorical,
    _ensure_float)
from pandas.core.dtypes.cast import maybe_downcast_to_dtype
from pandas.core.dtypes.generic import ABCSeries
from pandas.core.dtypes.missing import isna, isnull, notna, _maybe_fill

from pandas.core.base import (PandasObject, SelectionMixin, GroupByError,
                              DataError, SpecificationError)
from pandas.core.index import (Index, MultiIndex,
                               CategoricalIndex, _ensure_index)
from pandas.core.arrays import ExtensionArray, Categorical
from pandas.core.frame import DataFrame
from pandas.core.generic import NDFrame, _shared_docs
from pandas.core.internals import BlockManager, make_block
from pandas.core.series import Series
from pandas.core.panel import Panel
from pandas.core.sorting import (get_group_index_sorter, get_group_index,
                                 compress_group_index, get_flattened_iterator,
                                 decons_obs_group_ids, get_indexer_dict)
from pandas.util._decorators import (cache_readonly, Substitution,
                                     Appender, make_signature)
from pandas.io.formats.printing import pprint_thing
from pandas.util._validators import validate_kwargs

import pandas.core.common as com
import pandas.core.algorithms as algorithms
from pandas.core.config import option_context

from pandas.plotting._core import boxplot_frame_groupby

from pandas._libs import (lib, reduction,
                          groupby as libgroupby,
                          Timestamp, NaT, iNaT)
from pandas._libs.lib import count_level_2d

_doc_template = """

        See also
        --------
        pandas.Series.%(name)s
        pandas.DataFrame.%(name)s
        pandas.Panel.%(name)s
"""

_apply_docs = dict(
    template="""
    Apply function ``func``  group-wise and combine the results together.

    The function passed to ``apply`` must take a {input} as its first
    argument and return a dataframe, a series or a scalar. ``apply`` will
    then take care of combining the results back together into a single
    dataframe or series. ``apply`` is therefore a highly flexible
    grouping method.

    While ``apply`` is a very flexible method, its downside is that
    using it can be quite a bit slower than using more specific methods.
    Pandas offers a wide range of method that will be much faster
    than using ``apply`` for their specific purposes, so try to use them
    before reaching for ``apply``.

    Parameters
    ----------
    func : function
        A callable that takes a {input} as its first argument, and
        returns a dataframe, a series or a scalar. In addition the
        callable may take positional and keyword arguments
    args, kwargs : tuple and dict
        Optional positional and keyword arguments to pass to ``func``

    Returns
    -------
    applied : Series or DataFrame

    Notes
    -----
    In the current implementation ``apply`` calls func twice on the
    first group to decide whether it can take a fast or slow code
    path. This can lead to unexpected behavior if func has
    side-effects, as they will take effect twice for the first
    group.

    Examples
    --------
    {examples}

    See also
    --------
    pipe : Apply function to the full GroupBy object instead of to each
        group.
    aggregate, transform
    """,
    dataframe_examples="""
    >>> df = pd.DataFrame({'A': 'a a b'.split(), 'B': [1,2,3], 'C': [4,6, 5]})
    >>> g = df.groupby('A')

    From ``df`` above we can see that ``g`` has two groups, ``a``, ``b``.
    Calling ``apply`` in various ways, we can get different grouping results:

    Example 1: below the function passed to ``apply`` takes a dataframe as
    its argument and returns a dataframe. ``apply`` combines the result for
    each group together into a new dataframe:

    >>> g.apply(lambda x: x / x.sum())
              B    C
    0  0.333333  0.4
    1  0.666667  0.6
    2  1.000000  1.0

    Example 2: The function passed to ``apply`` takes a dataframe as
    its argument and returns a series.  ``apply`` combines the result for
    each group together into a new dataframe:

    >>> g.apply(lambda x: x.max() - x.min())
       B  C
    A
    a  1  2
    b  0  0

    Example 3: The function passed to ``apply`` takes a dataframe as
    its argument and returns a scalar. ``apply`` combines the result for
    each group together into a series, including setting the index as
    appropriate:

    >>> g.apply(lambda x: x.C.max() - x.B.min())
    A
    a    5
    b    2
    dtype: int64
    """,
    series_examples="""
    >>> ser = pd.Series([0, 1, 2], index='a a b'.split())
    >>> g = ser.groupby(ser.index)

    From ``ser`` above we can see that ``g`` has two groups, ``a``, ``b``.
    Calling ``apply`` in various ways, we can get different grouping results:

    Example 1: The function passed to ``apply`` takes a series as
    its argument and returns a series.  ``apply`` combines the result for
    each group together into a new series:

    >>> g.apply(lambda x:  x*2 if x.name == 'b' else x/2)
    0    0.0
    1    0.5
    2    4.0
    dtype: float64

    Example 2: The function passed to ``apply`` takes a series as
    its argument and returns a scalar. ``apply`` combines the result for
    each group together into a series, including setting the index as
    appropriate:

    >>> g.apply(lambda x: x.max() - x.min())
    a    1
    b    0
    dtype: int64
    """)

_pipe_template = """\
Apply a function ``func`` with arguments to this %(klass)s object and return
the function's result.

%(versionadded)s

Use ``.pipe`` when you want to improve readability by chaining together
functions that expect Series, DataFrames, GroupBy or Resampler objects.
Instead of writing

>>> h(g(f(df.groupby('group')), arg1=a), arg2=b, arg3=c)

You can write

>>> (df.groupby('group')
...    .pipe(f)
...    .pipe(g, arg1=a)
...    .pipe(h, arg2=b, arg3=c))

which is much more readable.

Parameters
----------
func : callable or tuple of (callable, string)
    Function to apply to this %(klass)s object or, alternatively,
    a ``(callable, data_keyword)`` tuple where ``data_keyword`` is a
    string indicating the keyword of ``callable`` that expects the
    %(klass)s object.
args : iterable, optional
       positional arguments passed into ``func``.
kwargs : dict, optional
         a dictionary of keyword arguments passed into ``func``.

Returns
-------
object : the return type of ``func``.

Notes
-----
See more `here
<http://pandas.pydata.org/pandas-docs/stable/groupby.html#piping-function-calls>`_

Examples
--------
%(examples)s

See Also
--------
pandas.Series.pipe : Apply a function with arguments to a series
pandas.DataFrame.pipe: Apply a function with arguments to a dataframe
apply : Apply function to each group instead of to the
    full %(klass)s object.
"""

_transform_template = """
Call function producing a like-indexed %(klass)s on each group and
return a %(klass)s having the same indexes as the original object
filled with the transformed values

Parameters
----------
f : function
    Function to apply to each group

Notes
-----
Each group is endowed the attribute 'name' in case you need to know
which group you are working on.

The current implementation imposes three requirements on f:

* f must return a value that either has the same shape as the input
  subframe or can be broadcast to the shape of the input subframe.
  For example, f returns a scalar it will be broadcast to have the
  same shape as the input subframe.
* if this is a DataFrame, f must support application column-by-column
  in the subframe. If f also supports application to the entire subframe,
  then a fast path is used starting from the second chunk.
* f must not mutate groups. Mutation is not supported and may
  produce unexpected results.

Returns
-------
%(klass)s

See also
--------
aggregate, transform

Examples
--------

# Same shape
>>> df = pd.DataFrame({'A' : ['foo', 'bar', 'foo', 'bar',
...                           'foo', 'bar'],
...                    'B' : ['one', 'one', 'two', 'three',
...                          'two', 'two'],
...                    'C' : [1, 5, 5, 2, 5, 5],
...                    'D' : [2.0, 5., 8., 1., 2., 9.]})
>>> grouped = df.groupby('A')
>>> grouped.transform(lambda x: (x - x.mean()) / x.std())
          C         D
0 -1.154701 -0.577350
1  0.577350  0.000000
2  0.577350  1.154701
3 -1.154701 -1.000000
4  0.577350 -0.577350
5  0.577350  1.000000

# Broadcastable
>>> grouped.transform(lambda x: x.max() - x.min())
   C    D
0  4  6.0
1  3  8.0
2  4  6.0
3  3  8.0
4  4  6.0
5  3  8.0

"""


# special case to prevent duplicate plots when catching exceptions when
# forwarding methods from NDFrames
_plotting_methods = frozenset(['plot', 'boxplot', 'hist'])

_common_apply_whitelist = frozenset([
    'last', 'first',
    'head', 'tail', 'median',
    'mean', 'sum', 'min', 'max',
    'cumcount', 'ngroup',
    'resample',
    'rank', 'quantile',
    'fillna',
    'mad',
    'any', 'all',
    'take',
    'idxmax', 'idxmin',
    'shift', 'tshift',
    'ffill', 'bfill',
    'pct_change', 'skew',
    'corr', 'cov', 'diff',
]) | _plotting_methods

_series_apply_whitelist = ((_common_apply_whitelist |
                            {'nlargest', 'nsmallest',
                             'is_monotonic_increasing',
                             'is_monotonic_decreasing'}) -
                           {'boxplot'}) | frozenset(['dtype', 'unique'])

_dataframe_apply_whitelist = ((_common_apply_whitelist |
                              frozenset(['dtypes', 'corrwith'])) -
                              {'boxplot'})

_cython_transforms = frozenset(['cumprod', 'cumsum', 'shift',
                                'cummin', 'cummax'])

_cython_cast_blacklist = frozenset(['rank', 'count', 'size'])


class Grouper(object):
    """
    A Grouper allows the user to specify a groupby instruction for a target
    object

    This specification will select a column via the key parameter, or if the
    level and/or axis parameters are given, a level of the index of the target
    object.

    These are local specifications and will override 'global' settings,
    that is the parameters axis and level which are passed to the groupby
    itself.

    Parameters
    ----------
    key : string, defaults to None
        groupby key, which selects the grouping column of the target
    level : name/number, defaults to None
        the level for the target index
    freq : string / frequency object, defaults to None
        This will groupby the specified frequency if the target selection
        (via key or level) is a datetime-like object. For full specification
        of available frequencies, please see `here
        <http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases>`_.
    axis : number/name of the axis, defaults to 0
    sort : boolean, default to False
        whether to sort the resulting labels

    additional kwargs to control time-like groupers (when ``freq`` is passed)

    closed : closed end of interval; 'left' or 'right'
    label : interval boundary to use for labeling; 'left' or 'right'
    convention : {'start', 'end', 'e', 's'}
        If grouper is PeriodIndex
    base, loffset

    Returns
    -------
    A specification for a groupby instruction

    Examples
    --------

    Syntactic sugar for ``df.groupby('A')``

    >>> df.groupby(Grouper(key='A'))

    Specify a resample operation on the column 'date'

    >>> df.groupby(Grouper(key='date', freq='60s'))

    Specify a resample operation on the level 'date' on the columns axis
    with a frequency of 60s

    >>> df.groupby(Grouper(level='date', freq='60s', axis=1))
    """
    _attributes = ('key', 'level', 'freq', 'axis', 'sort')

    def __new__(cls, *args, **kwargs):
        if kwargs.get('freq') is not None:
            from pandas.core.resample import TimeGrouper
            cls = TimeGrouper
        return super(Grouper, cls).__new__(cls)

    def __init__(self, key=None, level=None, freq=None, axis=0, sort=False):
        self.key = key
        self.level = level
        self.freq = freq
        self.axis = axis
        self.sort = sort

        self.grouper = None
        self.obj = None
        self.indexer = None
        self.binner = None
        self._grouper = None

    @property
    def ax(self):
        return self.grouper

    def _get_grouper(self, obj, validate=True):
        """
        Parameters
        ----------
        obj : the subject object
        validate : boolean, default True
            if True, validate the grouper

        Returns
        -------
        a tuple of binner, grouper, obj (possibly sorted)
        """

        self._set_grouper(obj)
        self.grouper, exclusions, self.obj = _get_grouper(self.obj, [self.key],
                                                          axis=self.axis,
                                                          level=self.level,
                                                          sort=self.sort,
                                                          validate=validate)
        return self.binner, self.grouper, self.obj

    def _set_grouper(self, obj, sort=False):
        """
        given an object and the specifications, setup the internal grouper
        for this particular specification

        Parameters
        ----------
        obj : the subject object
        sort : bool, default False
            whether the resulting grouper should be sorted
        """

        if self.key is not None and self.level is not None:
            raise ValueError(
                "The Grouper cannot specify both a key and a level!")

        # Keep self.grouper value before overriding
        if self._grouper is None:
            self._grouper = self.grouper

        # the key must be a valid info item
        if self.key is not None:
            key = self.key
            # The 'on' is already defined
            if getattr(self.grouper, 'name', None) == key and \
                    isinstance(obj, ABCSeries):
                ax = self._grouper.take(obj.index)
            else:
                if key not in obj._info_axis:
                    raise KeyError(
                        "The grouper name {0} is not found".format(key))
                ax = Index(obj[key], name=key)

        else:
            ax = obj._get_axis(self.axis)
            if self.level is not None:
                level = self.level

                # if a level is given it must be a mi level or
                # equivalent to the axis name
                if isinstance(ax, MultiIndex):
                    level = ax._get_level_number(level)
                    ax = Index(ax._get_level_values(level),
                               name=ax.names[level])

                else:
                    if level not in (0, ax.name):
                        raise ValueError(
                            "The level {0} is not valid".format(level))

        # possibly sort
        if (self.sort or sort) and not ax.is_monotonic:
            # use stable sort to support first, last, nth
            indexer = self.indexer = ax.argsort(kind='mergesort')
            ax = ax.take(indexer)
            obj = obj._take(indexer, axis=self.axis, is_copy=False)

        self.obj = obj
        self.grouper = ax
        return self.grouper

    @property
    def groups(self):
        return self.grouper.groups

    def __repr__(self):
        attrs_list = ["{}={!r}".format(attr_name, getattr(self, attr_name))
                      for attr_name in self._attributes
                      if getattr(self, attr_name) is not None]
        attrs = ", ".join(attrs_list)
        cls_name = self.__class__.__name__
        return "{}({})".format(cls_name, attrs)


class GroupByPlot(PandasObject):
    """
    Class implementing the .plot attribute for groupby objects
    """

    def __init__(self, groupby):
        self._groupby = groupby

    def __call__(self, *args, **kwargs):
        def f(self):
            return self.plot(*args, **kwargs)
        f.__name__ = 'plot'
        return self._groupby.apply(f)

    def __getattr__(self, name):
        def attr(*args, **kwargs):
            def f(self):
                return getattr(self.plot, name)(*args, **kwargs)
            return self._groupby.apply(f)
        return attr


@contextmanager
def _group_selection_context(groupby):
    """
    set / reset the _group_selection_context
    """
    groupby._set_group_selection()
    yield groupby
    groupby._reset_group_selection()


class _GroupBy(PandasObject, SelectionMixin):
    _group_selection = None
    _apply_whitelist = frozenset([])

    def __init__(self, obj, keys=None, axis=0, level=None,
                 grouper=None, exclusions=None, selection=None, as_index=True,
                 sort=True, group_keys=True, squeeze=False,
                 observed=False, **kwargs):

        self._selection = selection

        if isinstance(obj, NDFrame):
            obj._consolidate_inplace()

        self.level = level

        if not as_index:
            if not isinstance(obj, DataFrame):
                raise TypeError('as_index=False only valid with DataFrame')
            if axis != 0:
                raise ValueError('as_index=False only valid for axis=0')

        self.as_index = as_index
        self.keys = keys
        self.sort = sort
        self.group_keys = group_keys
        self.squeeze = squeeze
        self.observed = observed
        self.mutated = kwargs.pop('mutated', False)

        if grouper is None:
            grouper, exclusions, obj = _get_grouper(obj, keys,
                                                    axis=axis,
                                                    level=level,
                                                    sort=sort,
                                                    observed=observed,
                                                    mutated=self.mutated)

        self.obj = obj
        self.axis = obj._get_axis_number(axis)
        self.grouper = grouper
        self.exclusions = set(exclusions) if exclusions else set()

        # we accept no other args
        validate_kwargs('group', kwargs, {})

    def __len__(self):
        return len(self.groups)

    def __unicode__(self):
        # TODO: Better unicode/repr for GroupBy object
        return object.__repr__(self)

    def _assure_grouper(self):
        """
        we create the grouper on instantiation
        sub-classes may have a different policy
        """
        pass

    @property
    def groups(self):
        """ dict {group name -> group labels} """
        self._assure_grouper()
        return self.grouper.groups

    @property
    def ngroups(self):
        self._assure_grouper()
        return self.grouper.ngroups

    @property
    def indices(self):
        """ dict {group name -> group indices} """
        self._assure_grouper()
        return self.grouper.indices

    def _get_indices(self, names):
        """
        safe get multiple indices, translate keys for
        datelike to underlying repr
        """

        def get_converter(s):
            # possibly convert to the actual key types
            # in the indices, could be a Timestamp or a np.datetime64
            if isinstance(s, (Timestamp, datetime.datetime)):
                return lambda key: Timestamp(key)
            elif isinstance(s, np.datetime64):
                return lambda key: Timestamp(key).asm8
            else:
                return lambda key: key

        if len(names) == 0:
            return []

        if len(self.indices) > 0:
            index_sample = next(iter(self.indices))
        else:
            index_sample = None     # Dummy sample

        name_sample = names[0]
        if isinstance(index_sample, tuple):
            if not isinstance(name_sample, tuple):
                msg = ("must supply a tuple to get_group with multiple"
                       " grouping keys")
                raise ValueError(msg)
            if not len(name_sample) == len(index_sample):
                try:
                    # If the original grouper was a tuple
                    return [self.indices[name] for name in names]
                except KeyError:
                    # turns out it wasn't a tuple
                    msg = ("must supply a a same-length tuple to get_group"
                           " with multiple grouping keys")
                    raise ValueError(msg)

            converters = [get_converter(s) for s in index_sample]
            names = [tuple(f(n) for f, n in zip(converters, name))
                     for name in names]

        else:
            converter = get_converter(index_sample)
            names = [converter(name) for name in names]

        return [self.indices.get(name, []) for name in names]

    def _get_index(self, name):
        """ safe get index, translate keys for datelike to underlying repr """
        return self._get_indices([name])[0]

    @cache_readonly
    def _selected_obj(self):

        if self._selection is None or isinstance(self.obj, Series):
            if self._group_selection is not None:
                return self.obj[self._group_selection]
            return self.obj
        else:
            return self.obj[self._selection]

    def _reset_group_selection(self):
        """
        Clear group based selection. Used for methods needing to return info on
        each group regardless of whether a group selection was previously set.
        """
        if self._group_selection is not None:
            # GH12839 clear cached selection too when changing group selection
            self._group_selection = None
            self._reset_cache('_selected_obj')

    def _set_group_selection(self):
        """
        Create group based selection. Used when selection is not passed
        directly but instead via a grouper.

        NOTE: this should be paired with a call to _reset_group_selection
        """
        grp = self.grouper
        if not (self.as_index and
                getattr(grp, 'groupings', None) is not None and
                self.obj.ndim > 1 and
                self._group_selection is None):
            return

        ax = self.obj._info_axis
        groupers = [g.name for g in grp.groupings
                    if g.level is None and g.in_axis]

        if len(groupers):
            # GH12839 clear selected obj cache when group selection changes
            self._group_selection = ax.difference(Index(groupers)).tolist()
            self._reset_cache('_selected_obj')

    def _set_result_index_ordered(self, result):
        # set the result index on the passed values object and
        # return the new object, xref 8046

        # the values/counts are repeated according to the group index
        # shortcut if we have an already ordered grouper
        if not self.grouper.is_monotonic:
            index = Index(np.concatenate(
                self._get_indices(self.grouper.result_index)))
            result.set_axis(index, axis=self.axis, inplace=True)
            result = result.sort_index(axis=self.axis)

        result.set_axis(self.obj._get_axis(self.axis), axis=self.axis,
                        inplace=True)
        return result

    def _dir_additions(self):
        return self.obj._dir_additions() | self._apply_whitelist

    def __getattr__(self, attr):
        if attr in self._internal_names_set:
            return object.__getattribute__(self, attr)
        if attr in self.obj:
            return self[attr]
        if hasattr(self.obj, attr):
            return self._make_wrapper(attr)

        raise AttributeError("%r object has no attribute %r" %
                             (type(self).__name__, attr))

    @Substitution(klass='GroupBy',
                  versionadded='.. versionadded:: 0.21.0',
                  examples="""\
>>> df = pd.DataFrame({'A': 'a b a b'.split(), 'B': [1, 2, 3, 4]})
>>> df
   A  B
0  a  1
1  b  2
2  a  3
3  b  4

To get the difference between each groups maximum and minimum value in one
pass, you can do

>>> df.groupby('A').pipe(lambda x: x.max() - x.min())
   B
A
a  2
b  2""")
    @Appender(_pipe_template)
    def pipe(self, func, *args, **kwargs):
        return com._pipe(self, func, *args, **kwargs)

    plot = property(GroupByPlot)

    def _make_wrapper(self, name):
        if name not in self._apply_whitelist:
            is_callable = callable(getattr(self._selected_obj, name, None))
            kind = ' callable ' if is_callable else ' '
            msg = ("Cannot access{0}attribute {1!r} of {2!r} objects, try "
                   "using the 'apply' method".format(kind, name,
                                                     type(self).__name__))
            raise AttributeError(msg)

        self._set_group_selection()

        # need to setup the selection
        # as are not passed directly but in the grouper
        f = getattr(self._selected_obj, name)
        if not isinstance(f, types.MethodType):
            return self.apply(lambda self: getattr(self, name))

        f = getattr(type(self._selected_obj), name)

        def wrapper(*args, **kwargs):
            # a little trickery for aggregation functions that need an axis
            # argument
            kwargs_with_axis = kwargs.copy()
            if 'axis' not in kwargs_with_axis or \
               kwargs_with_axis['axis'] is None:
                kwargs_with_axis['axis'] = self.axis

            def curried_with_axis(x):
                return f(x, *args, **kwargs_with_axis)

            def curried(x):
                return f(x, *args, **kwargs)

            # preserve the name so we can detect it when calling plot methods,
            # to avoid duplicates
            curried.__name__ = curried_with_axis.__name__ = name

            # special case otherwise extra plots are created when catching the
            # exception below
            if name in _plotting_methods:
                return self.apply(curried)

            try:
                return self.apply(curried_with_axis)
            except Exception:
                try:
                    return self.apply(curried)
                except Exception:

                    # related to : GH3688
                    # try item-by-item
                    # this can be called recursively, so need to raise
                    # ValueError
                    # if we don't have this method to indicated to aggregate to
                    # mark this column as an error
                    try:
                        return self._aggregate_item_by_item(name,
                                                            *args, **kwargs)
                    except (AttributeError):
                        raise ValueError

        return wrapper

    def get_group(self, name, obj=None):
        """
        Constructs NDFrame from group with provided name

        Parameters
        ----------
        name : object
            the name of the group to get as a DataFrame
        obj : NDFrame, default None
            the NDFrame to take the DataFrame out of.  If
            it is None, the object groupby was called on will
            be used

        Returns
        -------
        group : type of obj
        """
        if obj is None:
            obj = self._selected_obj

        inds = self._get_index(name)
        if not len(inds):
            raise KeyError(name)

        return obj._take(inds, axis=self.axis)

    def __iter__(self):
        """
        Groupby iterator

        Returns
        -------
        Generator yielding sequence of (name, subsetted object)
        for each group
        """
        return self.grouper.get_iterator(self.obj, axis=self.axis)

    @Appender(_apply_docs['template']
              .format(input="dataframe",
                      examples=_apply_docs['dataframe_examples']))
    def apply(self, func, *args, **kwargs):

        func = self._is_builtin_func(func)

        # this is needed so we don't try and wrap strings. If we could
        # resolve functions to their callable functions prior, this
        # wouldn't be needed
        if args or kwargs:
            if callable(func):

                @wraps(func)
                def f(g):
                    with np.errstate(all='ignore'):
                        return func(g, *args, **kwargs)
            else:
                raise ValueError('func must be a callable if args or '
                                 'kwargs are supplied')
        else:
            f = func

        # ignore SettingWithCopy here in case the user mutates
        with option_context('mode.chained_assignment', None):
            try:
                result = self._python_apply_general(f)
            except Exception:

                # gh-20949
                # try again, with .apply acting as a filtering
                # operation, by excluding the grouping column
                # This would normally not be triggered
                # except if the udf is trying an operation that
                # fails on *some* columns, e.g. a numeric operation
                # on a string grouper column

                with _group_selection_context(self):
                    return self._python_apply_general(f)

        return result

    def _python_apply_general(self, f):
        keys, values, mutated = self.grouper.apply(f, self._selected_obj,
                                                   self.axis)

        return self._wrap_applied_output(
            keys,
            values,
            not_indexed_same=mutated or self.mutated)

    def _iterate_slices(self):
        yield self._selection_name, self._selected_obj

    def transform(self, func, *args, **kwargs):
        raise com.AbstractMethodError(self)

    def _cumcount_array(self, ascending=True):
        """
        Parameters
        ----------
        ascending : bool, default True
            If False, number in reverse, from length of group - 1 to 0.

        Notes
        -----
        this is currently implementing sort=False
        (though the default is sort=True) for groupby in general
        """
        ids, _, ngroups = self.grouper.group_info
        sorter = get_group_index_sorter(ids, ngroups)
        ids, count = ids[sorter], len(ids)

        if count == 0:
            return np.empty(0, dtype=np.int64)

        run = np.r_[True, ids[:-1] != ids[1:]]
        rep = np.diff(np.r_[np.nonzero(run)[0], count])
        out = (~run).cumsum()

        if ascending:
            out -= np.repeat(out[run], rep)
        else:
            out = np.repeat(out[np.r_[run[1:], True]], rep) - out

        rev = np.empty(count, dtype=np.intp)
        rev[sorter] = np.arange(count, dtype=np.intp)
        return out[rev].astype(np.int64, copy=False)

    def _try_cast(self, result, obj, numeric_only=False):
        """
        try to cast the result to our obj original type,
        we may have roundtripped thru object in the mean-time

        if numeric_only is True, then only try to cast numerics
        and not datetimelikes

        """
        if obj.ndim > 1:
            dtype = obj.values.dtype
        else:
            dtype = obj.dtype

        if not is_scalar(result):
            if numeric_only and is_numeric_dtype(dtype) or not numeric_only:
                result = maybe_downcast_to_dtype(result, dtype)

        return result

    def _transform_should_cast(self, func_nm):
        """
        Parameters:
        -----------
        func_nm: str
            The name of the aggregation function being performed

        Returns:
        --------
        bool
            Whether transform should attempt to cast the result of aggregation
        """
        return (self.size().fillna(0) > 0).any() and (func_nm not in
                                                      _cython_cast_blacklist)

    def _cython_transform(self, how, numeric_only=True, **kwargs):
        output = collections.OrderedDict()
        for name, obj in self._iterate_slices():
            is_numeric = is_numeric_dtype(obj.dtype)
            if numeric_only and not is_numeric:
                continue

            try:
                result, names = self.grouper.transform(obj.values, how,
                                                       **kwargs)
            except NotImplementedError:
                continue
            except AssertionError as e:
                raise GroupByError(str(e))
            if self._transform_should_cast(how):
                output[name] = self._try_cast(result, obj)
            else:
                output[name] = result

        if len(output) == 0:
            raise DataError('No numeric types to aggregate')

        return self._wrap_transformed_output(output, names)

    def _cython_agg_general(self, how, alt=None, numeric_only=True,
                            min_count=-1):
        output = {}
        for name, obj in self._iterate_slices():
            is_numeric = is_numeric_dtype(obj.dtype)
            if numeric_only and not is_numeric:
                continue

            try:
                result, names = self.grouper.aggregate(obj.values, how,
                                                       min_count=min_count)
            except AssertionError as e:
                raise GroupByError(str(e))
            output[name] = self._try_cast(result, obj)

        if len(output) == 0:
            raise DataError('No numeric types to aggregate')

        return self._wrap_aggregated_output(output, names)

    def _python_agg_general(self, func, *args, **kwargs):
        func = self._is_builtin_func(func)
        f = lambda x: func(x, *args, **kwargs)

        # iterate through "columns" ex exclusions to populate output dict
        output = {}
        for name, obj in self._iterate_slices():
            try:
                result, counts = self.grouper.agg_series(obj, f)
                output[name] = self._try_cast(result, obj, numeric_only=True)
            except TypeError:
                continue

        if len(output) == 0:
            return self._python_apply_general(f)

        if self.grouper._filter_empty_groups:

            mask = counts.ravel() > 0
            for name, result in compat.iteritems(output):

                # since we are masking, make sure that we have a float object
                values = result
                if is_numeric_dtype(values.dtype):
                    values = _ensure_float(values)

                output[name] = self._try_cast(values[mask], result)

        return self._wrap_aggregated_output(output)

    def _wrap_applied_output(self, *args, **kwargs):
        raise com.AbstractMethodError(self)

    def _concat_objects(self, keys, values, not_indexed_same=False):
        from pandas.core.reshape.concat import concat

        def reset_identity(values):
            # reset the identities of the components
            # of the values to prevent aliasing
            for v in com._not_none(*values):
                ax = v._get_axis(self.axis)
                ax._reset_identity()
            return values

        if not not_indexed_same:
            result = concat(values, axis=self.axis)
            ax = self._selected_obj._get_axis(self.axis)

            if isinstance(result, Series):
                result = result.reindex(ax)
            else:

                # this is a very unfortunate situation
                # we have a multi-index that is NOT lexsorted
                # and we have a result which is duplicated
                # we can't reindex, so we resort to this
                # GH 14776
                if isinstance(ax, MultiIndex) and not ax.is_unique:
                    indexer = algorithms.unique1d(
                        result.index.get_indexer_for(ax.values))
                    result = result.take(indexer, axis=self.axis)
                else:
                    result = result.reindex(ax, axis=self.axis)

        elif self.group_keys:

            values = reset_identity(values)
            if self.as_index:

                # possible MI return case
                group_keys = keys
                group_levels = self.grouper.levels
                group_names = self.grouper.names

                result = concat(values, axis=self.axis, keys=group_keys,
                                levels=group_levels, names=group_names,
                                sort=False)
            else:

                # GH5610, returns a MI, with the first level being a
                # range index
                keys = list(range(len(values)))
                result = concat(values, axis=self.axis, keys=keys)
        else:
            values = reset_identity(values)
            result = concat(values, axis=self.axis)

        if (isinstance(result, Series) and
                getattr(self, '_selection_name', None) is not None):

            result.name = self._selection_name

        return result

    def _apply_filter(self, indices, dropna):
        if len(indices) == 0:
            indices = np.array([], dtype='int64')
        else:
            indices = np.sort(np.concatenate(indices))
        if dropna:
            filtered = self._selected_obj.take(indices, axis=self.axis)
        else:
            mask = np.empty(len(self._selected_obj.index), dtype=bool)
            mask.fill(False)
            mask[indices.astype(int)] = True
            # mask fails to broadcast when passed to where; broadcast manually.
            mask = np.tile(mask, list(self._selected_obj.shape[1:]) + [1]).T
            filtered = self._selected_obj.where(mask)  # Fill with NaNs.
        return filtered


class GroupBy(_GroupBy):

    """
    Class for grouping and aggregating relational data. See aggregate,
    transform, and apply functions on this object.

    It's easiest to use obj.groupby(...) to use GroupBy, but you can also do:

    ::

        grouped = groupby(obj, ...)

    Parameters
    ----------
    obj : pandas object
    axis : int, default 0
    level : int, default None
        Level of MultiIndex
    groupings : list of Grouping objects
        Most users should ignore this
    exclusions : array-like, optional
        List of columns to exclude
    name : string
        Most users should ignore this

    Notes
    -----
    After grouping, see aggregate, apply, and transform functions. Here are
    some other brief notes about usage. When grouping by multiple groups, the
    result index will be a MultiIndex (hierarchical) by default.

    Iteration produces (key, group) tuples, i.e. chunking the data by group. So
    you can write code like:

    ::

        grouped = obj.groupby(keys, axis=axis)
        for key, group in grouped:
            # do something with the data

    Function calls on GroupBy, if not specially implemented, "dispatch" to the
    grouped data. So if you group a DataFrame and wish to invoke the std()
    method on each group, you can simply do:

    ::

        df.groupby(mapper).std()

    rather than

    ::

        df.groupby(mapper).aggregate(np.std)

    You can pass arguments to these "wrapped" functions, too.

    See the online documentation for full exposition on these topics and much
    more

    Returns
    -------
    **Attributes**
    groups : dict
        {group name -> group labels}
    len(grouped) : int
        Number of groups
    """
    _apply_whitelist = _common_apply_whitelist

    def _bool_agg(self, val_test, skipna):
        """Shared func to call any / all Cython GroupBy implementations"""

        def objs_to_bool(vals):
            try:
                vals = vals.astype(np.bool)
            except ValueError:  # for objects
                vals = np.array([bool(x) for x in vals])

            return vals.view(np.uint8)

        def result_to_bool(result):
            return result.astype(np.bool, copy=False)

        return self._get_cythonized_result('group_any_all', self.grouper,
                                           aggregate=True,
                                           cython_dtype=np.uint8,
                                           needs_values=True,
                                           needs_mask=True,
                                           pre_processing=objs_to_bool,
                                           post_processing=result_to_bool,
                                           val_test=val_test, skipna=skipna)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def any(self, skipna=True):
        """
        Returns True if any value in the group is truthful, else False

        Parameters
        ----------
        skipna : bool, default True
            Flag to ignore nan values during truth testing
        """
        return self._bool_agg('any', skipna)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def all(self, skipna=True):
        """Returns True if all values in the group are truthful, else False

        Parameters
        ----------
        skipna : bool, default True
            Flag to ignore nan values during truth testing
        """
        return self._bool_agg('all', skipna)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def count(self):
        """Compute count of group, excluding missing values"""

        # defined here for API doc
        raise NotImplementedError

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def mean(self, *args, **kwargs):
        """
        Compute mean of groups, excluding missing values

        For multiple groupings, the result index will be a MultiIndex
        """
        nv.validate_groupby_func('mean', args, kwargs, ['numeric_only'])
        try:
            return self._cython_agg_general('mean', **kwargs)
        except GroupByError:
            raise
        except Exception:  # pragma: no cover
            with _group_selection_context(self):
                f = lambda x: x.mean(axis=self.axis, **kwargs)
                return self._python_agg_general(f)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def median(self, **kwargs):
        """
        Compute median of groups, excluding missing values

        For multiple groupings, the result index will be a MultiIndex
        """
        try:
            return self._cython_agg_general('median', **kwargs)
        except GroupByError:
            raise
        except Exception:  # pragma: no cover

            def f(x):
                if isinstance(x, np.ndarray):
                    x = Series(x)
                return x.median(axis=self.axis, **kwargs)
            with _group_selection_context(self):
                return self._python_agg_general(f)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def std(self, ddof=1, *args, **kwargs):
        """
        Compute standard deviation of groups, excluding missing values

        For multiple groupings, the result index will be a MultiIndex

        Parameters
        ----------
        ddof : integer, default 1
            degrees of freedom
        """

        # TODO: implement at Cython level?
        nv.validate_groupby_func('std', args, kwargs)
        return np.sqrt(self.var(ddof=ddof, **kwargs))

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def var(self, ddof=1, *args, **kwargs):
        """
        Compute variance of groups, excluding missing values

        For multiple groupings, the result index will be a MultiIndex

        Parameters
        ----------
        ddof : integer, default 1
            degrees of freedom
        """
        nv.validate_groupby_func('var', args, kwargs)
        if ddof == 1:
            return self._cython_agg_general('var', **kwargs)
        else:
            f = lambda x: x.var(ddof=ddof, **kwargs)
            with _group_selection_context(self):
                return self._python_agg_general(f)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def sem(self, ddof=1):
        """
        Compute standard error of the mean of groups, excluding missing values

        For multiple groupings, the result index will be a MultiIndex

        Parameters
        ----------
        ddof : integer, default 1
            degrees of freedom
        """

        return self.std(ddof=ddof) / np.sqrt(self.count())

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def size(self):
        """Compute group sizes"""
        result = self.grouper.size()

        if isinstance(self.obj, Series):
            result.name = getattr(self.obj, 'name', None)
        return result

    @classmethod
    def _add_numeric_operations(cls):
        """ add numeric operations to the GroupBy generically """

        def groupby_function(name, alias, npfunc,
                             numeric_only=True, _convert=False,
                             min_count=-1):

            _local_template = "Compute %(f)s of group values"

            @Substitution(name='groupby', f=name)
            @Appender(_doc_template)
            @Appender(_local_template)
            def f(self, **kwargs):
                if 'numeric_only' not in kwargs:
                    kwargs['numeric_only'] = numeric_only
                if 'min_count' not in kwargs:
                    kwargs['min_count'] = min_count

                self._set_group_selection()
                try:
                    return self._cython_agg_general(
                        alias, alt=npfunc, **kwargs)
                except AssertionError as e:
                    raise SpecificationError(str(e))
                except Exception:
                    result = self.aggregate(
                        lambda x: npfunc(x, axis=self.axis))
                    if _convert:
                        result = result._convert(datetime=True)
                    return result

            set_function_name(f, name, cls)

            return f

        def first_compat(x, axis=0):

            def first(x):

                x = np.asarray(x)
                x = x[notna(x)]
                if len(x) == 0:
                    return np.nan
                return x[0]

            if isinstance(x, DataFrame):
                return x.apply(first, axis=axis)
            else:
                return first(x)

        def last_compat(x, axis=0):

            def last(x):

                x = np.asarray(x)
                x = x[notna(x)]
                if len(x) == 0:
                    return np.nan
                return x[-1]

            if isinstance(x, DataFrame):
                return x.apply(last, axis=axis)
            else:
                return last(x)

        cls.sum = groupby_function('sum', 'add', np.sum, min_count=0)
        cls.prod = groupby_function('prod', 'prod', np.prod, min_count=0)
        cls.min = groupby_function('min', 'min', np.min, numeric_only=False)
        cls.max = groupby_function('max', 'max', np.max, numeric_only=False)
        cls.first = groupby_function('first', 'first', first_compat,
                                     numeric_only=False)
        cls.last = groupby_function('last', 'last', last_compat,
                                    numeric_only=False)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def ohlc(self):
        """
        Compute sum of values, excluding missing values
        For multiple groupings, the result index will be a MultiIndex
        """

        return self._apply_to_column_groupbys(
            lambda x: x._cython_agg_general('ohlc'))

    @Appender(DataFrame.describe.__doc__)
    def describe(self, **kwargs):
        with _group_selection_context(self):
            result = self.apply(lambda x: x.describe(**kwargs))
            if self.axis == 1:
                return result.T
            return result.unstack()

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def resample(self, rule, *args, **kwargs):
        """
        Provide resampling when using a TimeGrouper
        Return a new grouper with our resampler appended
        """
        from pandas.core.resample import get_resampler_for_grouping
        return get_resampler_for_grouping(self, rule, *args, **kwargs)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def rolling(self, *args, **kwargs):
        """
        Return a rolling grouper, providing rolling
        functionality per group

        """
        from pandas.core.window import RollingGroupby
        return RollingGroupby(self, *args, **kwargs)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def expanding(self, *args, **kwargs):
        """
        Return an expanding grouper, providing expanding
        functionality per group

        """
        from pandas.core.window import ExpandingGroupby
        return ExpandingGroupby(self, *args, **kwargs)

    def _fill(self, direction, limit=None):
        """Shared function for `pad` and `backfill` to call Cython method

        Parameters
        ----------
        direction : {'ffill', 'bfill'}
            Direction passed to underlying Cython function. `bfill` will cause
            values to be filled backwards. `ffill` and any other values will
            default to a forward fill
        limit : int, default None
            Maximum number of consecutive values to fill. If `None`, this
            method will convert to -1 prior to passing to Cython

        Returns
        -------
        `Series` or `DataFrame` with filled values

        See Also
        --------
        pad
        backfill
        """
        # Need int value for Cython
        if limit is None:
            limit = -1

        return self._get_cythonized_result('group_fillna_indexer',
                                           self.grouper, needs_mask=True,
                                           cython_dtype=np.int64,
                                           result_is_index=True,
                                           direction=direction, limit=limit)

    @Substitution(name='groupby')
    def pad(self, limit=None):
        """
        Forward fill the values

        Parameters
        ----------
        limit : integer, optional
            limit of how many values to fill

        See Also
        --------
        Series.pad
        DataFrame.pad
        Series.fillna
        DataFrame.fillna
        """
        return self._fill('ffill', limit=limit)
    ffill = pad

    @Substitution(name='groupby')
    def backfill(self, limit=None):
        """
        Backward fill the values

        Parameters
        ----------
        limit : integer, optional
            limit of how many values to fill

        See Also
        --------
        Series.backfill
        DataFrame.backfill
        Series.fillna
        DataFrame.fillna
        """
        return self._fill('bfill', limit=limit)
    bfill = backfill

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def nth(self, n, dropna=None):
        """
        Take the nth row from each group if n is an int, or a subset of rows
        if n is a list of ints.

        If dropna, will take the nth non-null row, dropna is either
        Truthy (if a Series) or 'all', 'any' (if a DataFrame);
        this is equivalent to calling dropna(how=dropna) before the
        groupby.

        Parameters
        ----------
        n : int or list of ints
            a single nth value for the row or a list of nth values
        dropna : None or str, optional
            apply the specified dropna operation before counting which row is
            the nth row. Needs to be None, 'any' or 'all'

        Examples
        --------

        >>> df = pd.DataFrame({'A': [1, 1, 2, 1, 2],
        ...                    'B': [np.nan, 2, 3, 4, 5]}, columns=['A', 'B'])
        >>> g = df.groupby('A')
        >>> g.nth(0)
             B
        A
        1  NaN
        2  3.0
        >>> g.nth(1)
             B
        A
        1  2.0
        2  5.0
        >>> g.nth(-1)
             B
        A
        1  4.0
        2  5.0
        >>> g.nth([0, 1])
             B
        A
        1  NaN
        1  2.0
        2  3.0
        2  5.0

        Specifying ``dropna`` allows count ignoring NaN

        >>> g.nth(0, dropna='any')
             B
        A
        1  2.0
        2  3.0

        NaNs denote group exhausted when using dropna

        >>> g.nth(3, dropna='any')
            B
        A
        1 NaN
        2 NaN

        Specifying ``as_index=False`` in ``groupby`` keeps the original index.

        >>> df.groupby('A', as_index=False).nth(1)
           A    B
        1  1  2.0
        4  2  5.0
        """

        if isinstance(n, int):
            nth_values = [n]
        elif isinstance(n, (set, list, tuple)):
            nth_values = list(set(n))
            if dropna is not None:
                raise ValueError(
                    "dropna option with a list of nth values is not supported")
        else:
            raise TypeError("n needs to be an int or a list/set/tuple of ints")

        nth_values = np.array(nth_values, dtype=np.intp)
        self._set_group_selection()

        if not dropna:
            mask = np.in1d(self._cumcount_array(), nth_values) | \
                np.in1d(self._cumcount_array(ascending=False) + 1, -nth_values)

            out = self._selected_obj[mask]
            if not self.as_index:
                return out

            ids, _, _ = self.grouper.group_info
            out.index = self.grouper.result_index[ids[mask]]

            return out.sort_index() if self.sort else out

        if dropna not in ['any', 'all']:
            if isinstance(self._selected_obj, Series) and dropna is True:
                warnings.warn("the dropna={dropna} keyword is deprecated,"
                              "use dropna='all' instead. "
                              "For a Series groupby, dropna must be "
                              "either None, 'any' or 'all'.".format(
                                  dropna=dropna),
                              FutureWarning,
                              stacklevel=2)
                dropna = 'all'
            else:
                # Note: when agg-ing picker doesn't raise this,
                # just returns NaN
                raise ValueError("For a DataFrame groupby, dropna must be "
                                 "either None, 'any' or 'all', "
                                 "(was passed %s)." % (dropna),)

        # old behaviour, but with all and any support for DataFrames.
        # modified in GH 7559 to have better perf
        max_len = n if n >= 0 else - 1 - n
        dropped = self.obj.dropna(how=dropna, axis=self.axis)

        # get a new grouper for our dropped obj
        if self.keys is None and self.level is None:

            # we don't have the grouper info available
            # (e.g. we have selected out
            # a column that is not in the current object)
            axis = self.grouper.axis
            grouper = axis[axis.isin(dropped.index)]

        else:

            # create a grouper with the original parameters, but on the dropped
            # object
            grouper, _, _ = _get_grouper(dropped, key=self.keys,
                                         axis=self.axis, level=self.level,
                                         sort=self.sort,
                                         mutated=self.mutated)

        grb = dropped.groupby(grouper, as_index=self.as_index, sort=self.sort)
        sizes, result = grb.size(), grb.nth(n)
        mask = (sizes < max_len).values

        # set the results which don't meet the criteria
        if len(result) and mask.any():
            result.loc[mask] = np.nan

        # reset/reindex to the original groups
        if len(self.obj) == len(dropped) or \
           len(result) == len(self.grouper.result_index):
            result.index = self.grouper.result_index
        else:
            result = result.reindex(self.grouper.result_index)

        return result

    @Substitution(name='groupby')
    def ngroup(self, ascending=True):
        """
        Number each group from 0 to the number of groups - 1.

        This is the enumerative complement of cumcount.  Note that the
        numbers given to the groups match the order in which the groups
        would be seen when iterating over the groupby object, not the
        order they are first observed.

        .. versionadded:: 0.20.2

        Parameters
        ----------
        ascending : bool, default True
            If False, number in reverse, from number of group - 1 to 0.

        Examples
        --------

        >>> df = pd.DataFrame({"A": list("aaabba")})
        >>> df
           A
        0  a
        1  a
        2  a
        3  b
        4  b
        5  a
        >>> df.groupby('A').ngroup()
        0    0
        1    0
        2    0
        3    1
        4    1
        5    0
        dtype: int64
        >>> df.groupby('A').ngroup(ascending=False)
        0    1
        1    1
        2    1
        3    0
        4    0
        5    1
        dtype: int64
        >>> df.groupby(["A", [1,1,2,3,2,1]]).ngroup()
        0    0
        1    0
        2    1
        3    3
        4    2
        5    0
        dtype: int64

        See also
        --------
        .cumcount : Number the rows in each group.
        """

        with _group_selection_context(self):
            index = self._selected_obj.index
            result = Series(self.grouper.group_info[0], index)
            if not ascending:
                result = self.ngroups - 1 - result
            return result

    @Substitution(name='groupby')
    def cumcount(self, ascending=True):
        """
        Number each item in each group from 0 to the length of that group - 1.

        Essentially this is equivalent to

        >>> self.apply(lambda x: Series(np.arange(len(x)), x.index))

        Parameters
        ----------
        ascending : bool, default True
            If False, number in reverse, from length of group - 1 to 0.

        Examples
        --------

        >>> df = pd.DataFrame([['a'], ['a'], ['a'], ['b'], ['b'], ['a']],
        ...                   columns=['A'])
        >>> df
           A
        0  a
        1  a
        2  a
        3  b
        4  b
        5  a
        >>> df.groupby('A').cumcount()
        0    0
        1    1
        2    2
        3    0
        4    1
        5    3
        dtype: int64
        >>> df.groupby('A').cumcount(ascending=False)
        0    3
        1    2
        2    1
        3    1
        4    0
        5    0
        dtype: int64

        See also
        --------
        .ngroup : Number the groups themselves.
        """

        with _group_selection_context(self):
            index = self._selected_obj.index
            cumcounts = self._cumcount_array(ascending=ascending)
            return Series(cumcounts, index)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def rank(self, method='average', ascending=True, na_option='keep',
             pct=False, axis=0):
        """
        Provides the rank of values within each group.

        Parameters
        ----------
        method : {'average', 'min', 'max', 'first', 'dense'}, default 'average'
            * average: average rank of group
            * min: lowest rank in group
            * max: highest rank in group
            * first: ranks assigned in order they appear in the array
            * dense: like 'min', but rank always increases by 1 between groups
        ascending : boolean, default True
            False for ranks by high (1) to low (N)
        na_option :  {'keep', 'top', 'bottom'}, default 'keep'
            * keep: leave NA values where they are
            * top: smallest rank if ascending
            * bottom: smallest rank if descending
        pct : boolean, default False
            Compute percentage rank of data within each group
        axis : int, default 0
            The axis of the object over which to compute the rank.

        Returns
        -----
        DataFrame with ranking of values within each group
        """
        return self._cython_transform('rank', numeric_only=False,
                                      ties_method=method, ascending=ascending,
                                      na_option=na_option, pct=pct, axis=axis)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def cumprod(self, axis=0, *args, **kwargs):
        """Cumulative product for each group"""
        nv.validate_groupby_func('cumprod', args, kwargs,
                                 ['numeric_only', 'skipna'])
        if axis != 0:
            return self.apply(lambda x: x.cumprod(axis=axis, **kwargs))

        return self._cython_transform('cumprod', **kwargs)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def cumsum(self, axis=0, *args, **kwargs):
        """Cumulative sum for each group"""
        nv.validate_groupby_func('cumsum', args, kwargs,
                                 ['numeric_only', 'skipna'])
        if axis != 0:
            return self.apply(lambda x: x.cumsum(axis=axis, **kwargs))

        return self._cython_transform('cumsum', **kwargs)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def cummin(self, axis=0, **kwargs):
        """Cumulative min for each group"""
        if axis != 0:
            return self.apply(lambda x: np.minimum.accumulate(x, axis))

        return self._cython_transform('cummin', numeric_only=False)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def cummax(self, axis=0, **kwargs):
        """Cumulative max for each group"""
        if axis != 0:
            return self.apply(lambda x: np.maximum.accumulate(x, axis))

        return self._cython_transform('cummax', numeric_only=False)

    def _get_cythonized_result(self, how, grouper, aggregate=False,
                               cython_dtype=None, needs_values=False,
                               needs_mask=False, needs_ngroups=False,
                               result_is_index=False,
                               pre_processing=None, post_processing=None,
                               **kwargs):
        """Get result for Cythonized functions

        Parameters
        ----------
        how : str, Cythonized function name to be called
        grouper : Grouper object containing pertinent group info
        aggregate : bool, default False
            Whether the result should be aggregated to match the number of
            groups
        cython_dtype : default None
            Type of the array that will be modified by the Cython call. If
            `None`, the type will be inferred from the values of each slice
        needs_values : bool, default False
            Whether the values should be a part of the Cython call
            signature
        needs_mask : bool, default False
            Whether boolean mask needs to be part of the Cython call
            signature
        needs_ngroups : bool, default False
            Whether number of groups is part of the Cython call signature
        result_is_index : bool, default False
            Whether the result of the Cython operation is an index of
            values to be retrieved, instead of the actual values themselves
        pre_processing : function, default None
            Function to be applied to `values` prior to passing to Cython
            Raises if `needs_values` is False
        post_processing : function, default None
            Function to be applied to result of Cython function
        **kwargs : dict
            Extra arguments to be passed back to Cython funcs

        Returns
        -------
        `Series` or `DataFrame`  with filled values
        """
        if result_is_index and aggregate:
            raise ValueError("'result_is_index' and 'aggregate' cannot both "
                             "be True!")
        if post_processing:
            if not callable(pre_processing):
                raise ValueError("'post_processing' must be a callable!")
        if pre_processing:
            if not callable(pre_processing):
                raise ValueError("'pre_processing' must be a callable!")
            if not needs_values:
                raise ValueError("Cannot use 'pre_processing' without "
                                 "specifying 'needs_values'!")

        labels, _, ngroups = grouper.group_info
        output = collections.OrderedDict()
        base_func = getattr(libgroupby, how)

        for name, obj in self._iterate_slices():
            if aggregate:
                result_sz = ngroups
            else:
                result_sz = len(obj.values)

            if not cython_dtype:
                cython_dtype = obj.values.dtype

            result = np.zeros(result_sz, dtype=cython_dtype)
            func = partial(base_func, result, labels)
            if needs_values:
                vals = obj.values
                if pre_processing:
                    vals = pre_processing(vals)
                func = partial(func, vals)

            if needs_mask:
                mask = isnull(obj.values).view(np.uint8)
                func = partial(func, mask)

            if needs_ngroups:
                func = partial(func, ngroups)

            func(**kwargs)  # Call func to modify indexer values in place

            if result_is_index:
                result = algorithms.take_nd(obj.values, result)

            if post_processing:
                result = post_processing(result)

            output[name] = result

        if aggregate:
            return self._wrap_aggregated_output(output)
        else:
            return self._wrap_transformed_output(output)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def shift(self, periods=1, freq=None, axis=0):
        """
        Shift each group by periods observations

        Parameters
        ----------
        periods : integer, default 1
            number of periods to shift
        freq : frequency string
        axis : axis to shift, default 0
        """

        if freq is not None or axis != 0:
            return self.apply(lambda x: x.shift(periods, freq, axis))

        return self._get_cythonized_result('group_shift_indexer',
                                           self.grouper, cython_dtype=np.int64,
                                           needs_ngroups=True,
                                           result_is_index=True,
                                           periods=periods)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def pct_change(self, periods=1, fill_method='pad', limit=None, freq=None,
                   axis=0):
        """Calcuate pct_change of each value to previous entry in group"""
        if freq is not None or axis != 0:
            return self.apply(lambda x: x.pct_change(periods=periods,
                                                     fill_method=fill_method,
                                                     limit=limit, freq=freq,
                                                     axis=axis))

        filled = getattr(self, fill_method)(limit=limit).drop(
            self.grouper.names, axis=1)
        shifted = filled.shift(periods=periods, freq=freq)

        return (filled / shifted) - 1

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def head(self, n=5):
        """
        Returns first n rows of each group.

        Essentially equivalent to ``.apply(lambda x: x.head(n))``,
        except ignores as_index flag.

        Examples
        --------

        >>> df = DataFrame([[1, 2], [1, 4], [5, 6]],
                           columns=['A', 'B'])
        >>> df.groupby('A', as_index=False).head(1)
           A  B
        0  1  2
        2  5  6
        >>> df.groupby('A').head(1)
           A  B
        0  1  2
        2  5  6
        """
        self._reset_group_selection()
        mask = self._cumcount_array() < n
        return self._selected_obj[mask]

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def tail(self, n=5):
        """
        Returns last n rows of each group

        Essentially equivalent to ``.apply(lambda x: x.tail(n))``,
        except ignores as_index flag.

        Examples
        --------

        >>> df = DataFrame([['a', 1], ['a', 2], ['b', 1], ['b', 2]],
                           columns=['A', 'B'])
        >>> df.groupby('A').tail(1)
           A  B
        1  a  2
        3  b  2
        >>> df.groupby('A').head(1)
           A  B
        0  a  1
        2  b  1
        """
        self._reset_group_selection()
        mask = self._cumcount_array(ascending=False) < n
        return self._selected_obj[mask]


GroupBy._add_numeric_operations()


@Appender(GroupBy.__doc__)
def groupby(obj, by, **kwds):
    if isinstance(obj, Series):
        klass = SeriesGroupBy
    elif isinstance(obj, DataFrame):
        klass = DataFrameGroupBy
    else:  # pragma: no cover
        raise TypeError('invalid type: %s' % type(obj))

    return klass(obj, by, **kwds)


def _get_axes(group):
    if isinstance(group, Series):
        return [group.index]
    else:
        return group.axes


def _is_indexed_like(obj, axes):
    if isinstance(obj, Series):
        if len(axes) > 1:
            return False
        return obj.index.equals(axes[0])
    elif isinstance(obj, DataFrame):
        return obj.index.equals(axes[0])

    return False


class BaseGrouper(object):
    """
    This is an internal Grouper class, which actually holds
    the generated groups

    Parameters
    ----------
    axis : int
        the axis to group
    groupings : array of grouping
        all the grouping instances to handle in this grouper
        for example for grouper list to groupby, need to pass the list
    sort : boolean, default True
        whether this grouper will give sorted result or not
    group_keys : boolean, default True
    mutated : boolean, default False
    indexer : intp array, optional
        the indexer created by Grouper
        some groupers (TimeGrouper) will sort its axis and its
        group_info is also sorted, so need the indexer to reorder

    """

    def __init__(self, axis, groupings, sort=True, group_keys=True,
                 mutated=False, indexer=None):
        self._filter_empty_groups = self.compressed = len(groupings) != 1
        self.axis = axis
        self.groupings = groupings
        self.sort = sort
        self.group_keys = group_keys
        self.mutated = mutated
        self.indexer = indexer

    @property
    def shape(self):
        return tuple(ping.ngroups for ping in self.groupings)

    def __iter__(self):
        return iter(self.indices)

    @property
    def nkeys(self):
        return len(self.groupings)

    def get_iterator(self, data, axis=0):
        """
        Groupby iterator

        Returns
        -------
        Generator yielding sequence of (name, subsetted object)
        for each group
        """
        splitter = self._get_splitter(data, axis=axis)
        keys = self._get_group_keys()
        for key, (i, group) in zip(keys, splitter):
            yield key, group

    def _get_splitter(self, data, axis=0):
        comp_ids, _, ngroups = self.group_info
        return get_splitter(data, comp_ids, ngroups, axis=axis)

    def _get_group_keys(self):
        if len(self.groupings) == 1:
            return self.levels[0]
        else:
            comp_ids, _, ngroups = self.group_info

            # provide "flattened" iterator for multi-group setting
            return get_flattened_iterator(comp_ids,
                                          ngroups,
                                          self.levels,
                                          self.labels)

    def apply(self, f, data, axis=0):
        mutated = self.mutated
        splitter = self._get_splitter(data, axis=axis)
        group_keys = self._get_group_keys()

        # oh boy
        f_name = com._get_callable_name(f)
        if (f_name not in _plotting_methods and
                hasattr(splitter, 'fast_apply') and axis == 0):
            try:
                values, mutated = splitter.fast_apply(f, group_keys)
                return group_keys, values, mutated
            except reduction.InvalidApply:
                # we detect a mutation of some kind
                # so take slow path
                pass
            except Exception:
                # raise this error to the caller
                pass

        result_values = []
        for key, (i, group) in zip(group_keys, splitter):
            object.__setattr__(group, 'name', key)

            # group might be modified
            group_axes = _get_axes(group)
            res = f(group)
            if not _is_indexed_like(res, group_axes):
                mutated = True
            result_values.append(res)

        return group_keys, result_values, mutated

    @cache_readonly
    def indices(self):
        """ dict {group name -> group indices} """
        if len(self.groupings) == 1:
            return self.groupings[0].indices
        else:
            label_list = [ping.labels for ping in self.groupings]
            keys = [com._values_from_object(ping.group_index)
                    for ping in self.groupings]
            return get_indexer_dict(label_list, keys)

    @property
    def labels(self):
        return [ping.labels for ping in self.groupings]

    @property
    def levels(self):
        return [ping.group_index for ping in self.groupings]

    @property
    def names(self):
        return [ping.name for ping in self.groupings]

    def size(self):
        """
        Compute group sizes

        """
        ids, _, ngroup = self.group_info
        ids = _ensure_platform_int(ids)
        if ngroup:
            out = np.bincount(ids[ids != -1], minlength=ngroup)
        else:
            out = ids
        return Series(out,
                      index=self.result_index,
                      dtype='int64')

    @cache_readonly
    def groups(self):
        """ dict {group name -> group labels} """
        if len(self.groupings) == 1:
            return self.groupings[0].groups
        else:
            to_groupby = lzip(*(ping.grouper for ping in self.groupings))
            to_groupby = Index(to_groupby)
            return self.axis.groupby(to_groupby)

    @cache_readonly
    def is_monotonic(self):
        # return if my group orderings are monotonic
        return Index(self.group_info[0]).is_monotonic

    @cache_readonly
    def group_info(self):
        comp_ids, obs_group_ids = self._get_compressed_labels()

        ngroups = len(obs_group_ids)
        comp_ids = _ensure_int64(comp_ids)
        return comp_ids, obs_group_ids, ngroups

    @cache_readonly
    def label_info(self):
        # return the labels of items in original grouped axis
        labels, _, _ = self.group_info
        if self.indexer is not None:
            sorter = np.lexsort((labels, self.indexer))
            labels = labels[sorter]
        return labels

    def _get_compressed_labels(self):
        all_labels = [ping.labels for ping in self.groupings]
        if len(all_labels) > 1:
            group_index = get_group_index(all_labels, self.shape,
                                          sort=True, xnull=True)
            return compress_group_index(group_index, sort=self.sort)

        ping = self.groupings[0]
        return ping.labels, np.arange(len(ping.group_index))

    @cache_readonly
    def ngroups(self):
        return len(self.result_index)

    @property
    def recons_labels(self):
        comp_ids, obs_ids, _ = self.group_info
        labels = (ping.labels for ping in self.groupings)
        return decons_obs_group_ids(
            comp_ids, obs_ids, self.shape, labels, xnull=True)

    @cache_readonly
    def result_index(self):
        if not self.compressed and len(self.groupings) == 1:
            return self.groupings[0].result_index.rename(self.names[0])

        labels = self.recons_labels
        levels = [ping.result_index for ping in self.groupings]
        result = MultiIndex(levels=levels,
                            labels=labels,
                            verify_integrity=False,
                            names=self.names)
        return result

    def get_group_levels(self):
        if not self.compressed and len(self.groupings) == 1:
            return [self.groupings[0].result_index]

        name_list = []
        for ping, labels in zip(self.groupings, self.recons_labels):
            labels = _ensure_platform_int(labels)
            levels = ping.result_index.take(labels)

            name_list.append(levels)

        return name_list

    # ------------------------------------------------------------
    # Aggregation functions

    _cython_functions = {
        'aggregate': {
            'add': 'group_add',
            'prod': 'group_prod',
            'min': 'group_min',
            'max': 'group_max',
            'mean': 'group_mean',
            'median': {
                'name': 'group_median'
            },
            'var': 'group_var',
            'first': {
                'name': 'group_nth',
                'f': lambda func, a, b, c, d, e: func(a, b, c, d, 1, -1)
            },
            'last': 'group_last',
            'ohlc': 'group_ohlc',
        },

        'transform': {
            'cumprod': 'group_cumprod',
            'cumsum': 'group_cumsum',
            'cummin': 'group_cummin',
            'cummax': 'group_cummax',
            'rank': {
                'name': 'group_rank',
                'f': lambda func, a, b, c, d, **kwargs: func(
                    a, b, c, d,
                    kwargs.get('ties_method', 'average'),
                    kwargs.get('ascending', True),
                    kwargs.get('pct', False),
                    kwargs.get('na_option', 'keep')
                )
            }
        }
    }

    _cython_arity = {
        'ohlc': 4,  # OHLC
    }

    _name_functions = {
        'ohlc': lambda *args: ['open', 'high', 'low', 'close']
    }

    def _is_builtin_func(self, arg):
        """
        if we define an builtin function for this argument, return it,
        otherwise return the arg
        """
        return SelectionMixin._builtin_table.get(arg, arg)

    def _get_cython_function(self, kind, how, values, is_numeric):

        dtype_str = values.dtype.name

        def get_func(fname):
            # see if there is a fused-type version of function
            # only valid for numeric
            f = getattr(libgroupby, fname, None)
            if f is not None and is_numeric:
                return f

            # otherwise find dtype-specific version, falling back to object
            for dt in [dtype_str, 'object']:
                f = getattr(libgroupby, "%s_%s" % (fname, dtype_str), None)
                if f is not None:
                    return f

        ftype = self._cython_functions[kind][how]

        if isinstance(ftype, dict):
            func = afunc = get_func(ftype['name'])

            # a sub-function
            f = ftype.get('f')
            if f is not None:

                def wrapper(*args, **kwargs):
                    return f(afunc, *args, **kwargs)

                # need to curry our sub-function
                func = wrapper

        else:
            func = get_func(ftype)

        if func is None:
            raise NotImplementedError("function is not implemented for this"
                                      "dtype: [how->%s,dtype->%s]" %
                                      (how, dtype_str))
        return func

    def _cython_operation(self, kind, values, how, axis, min_count=-1,
                          **kwargs):
        assert kind in ['transform', 'aggregate']

        # can we do this operation with our cython functions
        # if not raise NotImplementedError

        # we raise NotImplemented if this is an invalid operation
        # entirely, e.g. adding datetimes

        # categoricals are only 1d, so we
        # are not setup for dim transforming
        if is_categorical_dtype(values):
            raise NotImplementedError(
                "categoricals are not support in cython ops ATM")
        elif is_datetime64_any_dtype(values):
            if how in ['add', 'prod', 'cumsum', 'cumprod']:
                raise NotImplementedError(
                    "datetime64 type does not support {} "
                    "operations".format(how))
        elif is_timedelta64_dtype(values):
            if how in ['prod', 'cumprod']:
                raise NotImplementedError(
                    "timedelta64 type does not support {} "
                    "operations".format(how))

        arity = self._cython_arity.get(how, 1)

        vdim = values.ndim
        swapped = False
        if vdim == 1:
            values = values[:, None]
            out_shape = (self.ngroups, arity)
        else:
            if axis > 0:
                swapped = True
                values = values.swapaxes(0, axis)
            if arity > 1:
                raise NotImplementedError("arity of more than 1 is not "
                                          "supported for the 'how' argument")
            out_shape = (self.ngroups,) + values.shape[1:]

        is_datetimelike = needs_i8_conversion(values.dtype)
        is_numeric = is_numeric_dtype(values.dtype)

        if is_datetimelike:
            values = values.view('int64')
            is_numeric = True
        elif is_bool_dtype(values.dtype):
            values = _ensure_float64(values)
        elif is_integer_dtype(values):
            # we use iNaT for the missing value on ints
            # so pre-convert to guard this condition
            if (values == iNaT).any():
                values = _ensure_float64(values)
            else:
                values = values.astype('int64', copy=False)
        elif is_numeric and not is_complex_dtype(values):
            values = _ensure_float64(values)
        else:
            values = values.astype(object)

        try:
            func = self._get_cython_function(
                kind, how, values, is_numeric)
        except NotImplementedError:
            if is_numeric:
                values = _ensure_float64(values)
                func = self._get_cython_function(
                    kind, how, values, is_numeric)
            else:
                raise

        if how == 'rank':
            out_dtype = 'float'
        else:
            if is_numeric:
                out_dtype = '%s%d' % (values.dtype.kind, values.dtype.itemsize)
            else:
                out_dtype = 'object'

        labels, _, _ = self.group_info

        if kind == 'aggregate':
            result = _maybe_fill(np.empty(out_shape, dtype=out_dtype),
                                 fill_value=np.nan)
            counts = np.zeros(self.ngroups, dtype=np.int64)
            result = self._aggregate(
                result, counts, values, labels, func, is_numeric,
                is_datetimelike, min_count)
        elif kind == 'transform':
            result = _maybe_fill(np.empty_like(values, dtype=out_dtype),
                                 fill_value=np.nan)

            # TODO: min_count
            result = self._transform(
                result, values, labels, func, is_numeric, is_datetimelike,
                **kwargs)

        if is_integer_dtype(result) and not is_datetimelike:
            mask = result == iNaT
            if mask.any():
                result = result.astype('float64')
                result[mask] = np.nan

        if kind == 'aggregate' and \
           self._filter_empty_groups and not counts.all():
            if result.ndim == 2:
                try:
                    result = lib.row_bool_subset(
                        result, (counts > 0).view(np.uint8))
                except ValueError:
                    result = lib.row_bool_subset_object(
                        _ensure_object(result),
                        (counts > 0).view(np.uint8))
            else:
                result = result[counts > 0]

        if vdim == 1 and arity == 1:
            result = result[:, 0]

        if how in self._name_functions:
            # TODO
            names = self._name_functions[how]()
        else:
            names = None

        if swapped:
            result = result.swapaxes(0, axis)

        return result, names

    def aggregate(self, values, how, axis=0, min_count=-1):
        return self._cython_operation('aggregate', values, how, axis,
                                      min_count=min_count)

    def transform(self, values, how, axis=0, **kwargs):
        return self._cython_operation('transform', values, how, axis, **kwargs)

    def _aggregate(self, result, counts, values, comp_ids, agg_func,
                   is_numeric, is_datetimelike, min_count=-1):
        if values.ndim > 3:
            # punting for now
            raise NotImplementedError("number of dimensions is currently "
                                      "limited to 3")
        elif values.ndim > 2:
            for i, chunk in enumerate(values.transpose(2, 0, 1)):

                chunk = chunk.squeeze()
                agg_func(result[:, :, i], counts, chunk, comp_ids,
                         min_count)
        else:
            agg_func(result, counts, values, comp_ids, min_count)

        return result

    def _transform(self, result, values, comp_ids, transform_func,
                   is_numeric, is_datetimelike, **kwargs):

        comp_ids, _, ngroups = self.group_info
        if values.ndim > 3:
            # punting for now
            raise NotImplementedError("number of dimensions is currently "
                                      "limited to 3")
        elif values.ndim > 2:
            for i, chunk in enumerate(values.transpose(2, 0, 1)):

                chunk = chunk.squeeze()
                transform_func(result[:, :, i], values,
                               comp_ids, is_datetimelike, **kwargs)
        else:
            transform_func(result, values, comp_ids, is_datetimelike, **kwargs)

        return result

    def agg_series(self, obj, func):
        try:
            return self._aggregate_series_fast(obj, func)
        except Exception:
            return self._aggregate_series_pure_python(obj, func)

    def _aggregate_series_fast(self, obj, func):
        func = self._is_builtin_func(func)

        if obj.index._has_complex_internals:
            raise TypeError('Incompatible index for Cython grouper')

        group_index, _, ngroups = self.group_info

        # avoids object / Series creation overhead
        dummy = obj._get_values(slice(None, 0)).to_dense()
        indexer = get_group_index_sorter(group_index, ngroups)
        obj = obj._take(indexer).to_dense()
        group_index = algorithms.take_nd(
            group_index, indexer, allow_fill=False)
        grouper = reduction.SeriesGrouper(obj, func, group_index, ngroups,
                                          dummy)
        result, counts = grouper.get_result()
        return result, counts

    def _aggregate_series_pure_python(self, obj, func):

        group_index, _, ngroups = self.group_info

        counts = np.zeros(ngroups, dtype=int)
        result = None

        splitter = get_splitter(obj, group_index, ngroups, axis=self.axis)

        for label, group in splitter:
            res = func(group)
            if result is None:
                if (isinstance(res, (Series, Index, np.ndarray))):
                    raise ValueError('Function does not reduce')
                result = np.empty(ngroups, dtype='O')

            counts[label] = group.shape[0]
            result[label] = res

        result = lib.maybe_convert_objects(result, try_float=0)
        return result, counts


def generate_bins_generic(values, binner, closed):
    """
    Generate bin edge offsets and bin labels for one array using another array
    which has bin edge values. Both arrays must be sorted.

    Parameters
    ----------
    values : array of values
    binner : a comparable array of values representing bins into which to bin
        the first array. Note, 'values' end-points must fall within 'binner'
        end-points.
    closed : which end of bin is closed; left (default), right

    Returns
    -------
    bins : array of offsets (into 'values' argument) of bins.
        Zero and last edge are excluded in result, so for instance the first
        bin is values[0:bin[0]] and the last is values[bin[-1]:]
    """
    lenidx = len(values)
    lenbin = len(binner)

    if lenidx <= 0 or lenbin <= 0:
        raise ValueError("Invalid length for values or for binner")

    # check binner fits data
    if values[0] < binner[0]:
        raise ValueError("Values falls before first bin")

    if values[lenidx - 1] > binner[lenbin - 1]:
        raise ValueError("Values falls after last bin")

    bins = np.empty(lenbin - 1, dtype=np.int64)

    j = 0  # index into values
    bc = 0  # bin count

    # linear scan, presume nothing about values/binner except that it fits ok
    for i in range(0, lenbin - 1):
        r_bin = binner[i + 1]

        # count values in current bin, advance to next bin
        while j < lenidx and (values[j] < r_bin or
                              (closed == 'right' and values[j] == r_bin)):
            j += 1

        bins[bc] = j
        bc += 1

    return bins


class BinGrouper(BaseGrouper):

    """
    This is an internal Grouper class

    Parameters
    ----------
    bins : the split index of binlabels to group the item of axis
    binlabels : the label list
    filter_empty : boolean, default False
    mutated : boolean, default False
    indexer : a intp array

    Examples
    --------
    bins: [2, 4, 6, 8, 10]
    binlabels: DatetimeIndex(['2005-01-01', '2005-01-03',
        '2005-01-05', '2005-01-07', '2005-01-09'],
        dtype='datetime64[ns]', freq='2D')

    the group_info, which contains the label of each item in grouped
    axis, the index of label in label list, group number, is

    (array([0, 0, 1, 1, 2, 2, 3, 3, 4, 4]), array([0, 1, 2, 3, 4]), 5)

    means that, the grouped axis has 10 items, can be grouped into 5
    labels, the first and second items belong to the first label, the
    third and forth items belong to the second label, and so on

    """

    def __init__(self, bins, binlabels, filter_empty=False, mutated=False,
                 indexer=None):
        self.bins = _ensure_int64(bins)
        self.binlabels = _ensure_index(binlabels)
        self._filter_empty_groups = filter_empty
        self.mutated = mutated
        self.indexer = indexer

    @cache_readonly
    def groups(self):
        """ dict {group name -> group labels} """

        # this is mainly for compat
        # GH 3881
        result = {}
        for key, value in zip(self.binlabels, self.bins):
            if key is not NaT:
                result[key] = value
        return result

    @property
    def nkeys(self):
        return 1

    def get_iterator(self, data, axis=0):
        """
        Groupby iterator

        Returns
        -------
        Generator yielding sequence of (name, subsetted object)
        for each group
        """
        if isinstance(data, NDFrame):
            slicer = lambda start, edge: data._slice(
                slice(start, edge), axis=axis)
            length = len(data.axes[axis])
        else:
            slicer = lambda start, edge: data[slice(start, edge)]
            length = len(data)

        start = 0
        for edge, label in zip(self.bins, self.binlabels):
            if label is not NaT:
                yield label, slicer(start, edge)
            start = edge

        if start < length:
            yield self.binlabels[-1], slicer(start, None)

    @cache_readonly
    def indices(self):
        indices = collections.defaultdict(list)

        i = 0
        for label, bin in zip(self.binlabels, self.bins):
            if i < bin:
                if label is not NaT:
                    indices[label] = list(range(i, bin))
                i = bin
        return indices

    @cache_readonly
    def group_info(self):
        ngroups = self.ngroups
        obs_group_ids = np.arange(ngroups)
        rep = np.diff(np.r_[0, self.bins])

        rep = _ensure_platform_int(rep)
        if ngroups == len(self.bins):
            comp_ids = np.repeat(np.arange(ngroups), rep)
        else:
            comp_ids = np.repeat(np.r_[-1, np.arange(ngroups)], rep)

        return comp_ids.astype('int64', copy=False), \
            obs_group_ids.astype('int64', copy=False), ngroups

    @cache_readonly
    def ngroups(self):
        return len(self.result_index)

    @cache_readonly
    def result_index(self):
        if len(self.binlabels) != 0 and isna(self.binlabels[0]):
            return self.binlabels[1:]

        return self.binlabels

    @property
    def levels(self):
        return [self.binlabels]

    @property
    def names(self):
        return [self.binlabels.name]

    @property
    def groupings(self):
        return [Grouping(lvl, lvl, in_axis=False, level=None, name=name)
                for lvl, name in zip(self.levels, self.names)]

    def agg_series(self, obj, func):
        dummy = obj[:0]
        grouper = reduction.SeriesBinGrouper(obj, func, self.bins, dummy)
        return grouper.get_result()

    # ----------------------------------------------------------------------
    # cython aggregation

    _cython_functions = copy.deepcopy(BaseGrouper._cython_functions)


class Grouping(object):

    """
    Holds the grouping information for a single key

    Parameters
    ----------
    index : Index
    grouper :
    obj :
    name :
    level :
    observed : boolean, default False
        If we are a Categorical, use the observed values
    in_axis : if the Grouping is a column in self.obj and hence among
        Groupby.exclusions list

    Returns
    -------
    **Attributes**:
      * indices : dict of {group -> index_list}
      * labels : ndarray, group labels
      * ids : mapping of label -> group
      * counts : array of group counts
      * group_index : unique groups
      * groups : dict of {group -> label_list}
    """

    def __init__(self, index, grouper=None, obj=None, name=None, level=None,
                 sort=True, observed=False, in_axis=False):

        self.name = name
        self.level = level
        self.grouper = _convert_grouper(index, grouper)
        self.all_grouper = None
        self.index = index
        self.sort = sort
        self.obj = obj
        self.observed = observed
        self.in_axis = in_axis

        # right place for this?
        if isinstance(grouper, (Series, Index)) and name is None:
            self.name = grouper.name

        if isinstance(grouper, MultiIndex):
            self.grouper = grouper.values

        # we have a single grouper which may be a myriad of things,
        # some of which are dependent on the passing in level

        if level is not None:
            if not isinstance(level, int):
                if level not in index.names:
                    raise AssertionError('Level %s not in index' % str(level))
                level = index.names.index(level)

            if self.name is None:
                self.name = index.names[level]

            self.grouper, self._labels, self._group_index = \
                index._get_grouper_for_level(self.grouper, level)

        # a passed Grouper like, directly get the grouper in the same way
        # as single grouper groupby, use the group_info to get labels
        elif isinstance(self.grouper, Grouper):
            # get the new grouper; we already have disambiguated
            # what key/level refer to exactly, don't need to
            # check again as we have by this point converted these
            # to an actual value (rather than a pd.Grouper)
            _, grouper, _ = self.grouper._get_grouper(self.obj, validate=False)
            if self.name is None:
                self.name = grouper.result_index.name
            self.obj = self.grouper.obj
            self.grouper = grouper

        else:
            if self.grouper is None and self.name is not None:
                self.grouper = self.obj[self.name]

            elif isinstance(self.grouper, (list, tuple)):
                self.grouper = com._asarray_tuplesafe(self.grouper)

            # a passed Categorical
            elif is_categorical_dtype(self.grouper):

                self.all_grouper = self.grouper
                self.grouper = self.grouper._codes_for_groupby(
                    self.sort, observed)
                categories = self.grouper.categories

                # we make a CategoricalIndex out of the cat grouper
                # preserving the categories / ordered attributes
                self._labels = self.grouper.codes
                if observed:
                    codes = algorithms.unique1d(self.grouper.codes)
                else:
                    codes = np.arange(len(categories))

                self._group_index = CategoricalIndex(
                    Categorical.from_codes(
                        codes=codes,
                        categories=categories,
                        ordered=self.grouper.ordered))

            # we are done
            if isinstance(self.grouper, Grouping):
                self.grouper = self.grouper.grouper

            # no level passed
            elif not isinstance(self.grouper,
                                (Series, Index, ExtensionArray, np.ndarray)):
                if getattr(self.grouper, 'ndim', 1) != 1:
                    t = self.name or str(type(self.grouper))
                    raise ValueError("Grouper for '%s' not 1-dimensional" % t)
                self.grouper = self.index.map(self.grouper)
                if not (hasattr(self.grouper, "__len__") and
                        len(self.grouper) == len(self.index)):
                    errmsg = ('Grouper result violates len(labels) == '
                              'len(data)\nresult: %s' %
                              pprint_thing(self.grouper))
                    self.grouper = None  # Try for sanity
                    raise AssertionError(errmsg)

        # if we have a date/time-like grouper, make sure that we have
        # Timestamps like
        if getattr(self.grouper, 'dtype', None) is not None:
            if is_datetime64_dtype(self.grouper):
                from pandas import to_datetime
                self.grouper = to_datetime(self.grouper)
            elif is_timedelta64_dtype(self.grouper):
                from pandas import to_timedelta
                self.grouper = to_timedelta(self.grouper)

    def __repr__(self):
        return 'Grouping({0})'.format(self.name)

    def __iter__(self):
        return iter(self.indices)

    _labels = None
    _group_index = None

    @property
    def ngroups(self):
        return len(self.group_index)

    @cache_readonly
    def indices(self):
        # we have a list of groupers
        if isinstance(self.grouper, BaseGrouper):
            return self.grouper.indices

        values = _ensure_categorical(self.grouper)
        return values._reverse_indexer()

    @property
    def labels(self):
        if self._labels is None:
            self._make_labels()
        return self._labels

    @cache_readonly
    def result_index(self):
        if self.all_grouper is not None:
            all_categories = self.all_grouper.categories

            # we re-order to the original category orderings
            if self.sort:
                return self.group_index.set_categories(all_categories)

            # we are not sorting, so add unobserved to the end
            categories = self.group_index.categories
            return self.group_index.add_categories(
                all_categories[~all_categories.isin(categories)])

        return self.group_index

    @property
    def group_index(self):
        if self._group_index is None:
            self._make_labels()
        return self._group_index

    def _make_labels(self):
        if self._labels is None or self._group_index is None:
            # we have a list of groupers
            if isinstance(self.grouper, BaseGrouper):
                labels = self.grouper.label_info
                uniques = self.grouper.result_index
            else:
                labels, uniques = algorithms.factorize(
                    self.grouper, sort=self.sort)
                uniques = Index(uniques, name=self.name)
            self._labels = labels
            self._group_index = uniques

    @cache_readonly
    def groups(self):
        return self.index.groupby(Categorical.from_codes(self.labels,
                                                         self.group_index))


def _get_grouper(obj, key=None, axis=0, level=None, sort=True,
                 observed=False, mutated=False, validate=True):
    """
    create and return a BaseGrouper, which is an internal
    mapping of how to create the grouper indexers.
    This may be composed of multiple Grouping objects, indicating
    multiple groupers

    Groupers are ultimately index mappings. They can originate as:
    index mappings, keys to columns, functions, or Groupers

    Groupers enable local references to axis,level,sort, while
    the passed in axis, level, and sort are 'global'.

    This routine tries to figure out what the passing in references
    are and then creates a Grouping for each one, combined into
    a BaseGrouper.

    If observed & we have a categorical grouper, only show the observed
    values

    If validate, then check for key/level overlaps

    """
    group_axis = obj._get_axis(axis)

    # validate that the passed single level is compatible with the passed
    # axis of the object
    if level is not None:
        # TODO: These if-block and else-block are almost same.
        # MultiIndex instance check is removable, but it seems that there are
        # some processes only for non-MultiIndex in else-block,
        # eg. `obj.index.name != level`. We have to consider carefully whether
        # these are applicable for MultiIndex. Even if these are applicable,
        # we need to check if it makes no side effect to subsequent processes
        # on the outside of this condition.
        # (GH 17621)
        if isinstance(group_axis, MultiIndex):
            if is_list_like(level) and len(level) == 1:
                level = level[0]

            if key is None and is_scalar(level):
                # Get the level values from group_axis
                key = group_axis.get_level_values(level)
                level = None

        else:
            # allow level to be a length-one list-like object
            # (e.g., level=[0])
            # GH 13901
            if is_list_like(level):
                nlevels = len(level)
                if nlevels == 1:
                    level = level[0]
                elif nlevels == 0:
                    raise ValueError('No group keys passed!')
                else:
                    raise ValueError('multiple levels only valid with '
                                     'MultiIndex')

            if isinstance(level, compat.string_types):
                if obj.index.name != level:
                    raise ValueError('level name %s is not the name of the '
                                     'index' % level)
            elif level > 0 or level < -1:
                raise ValueError('level > 0 or level < -1 only valid with '
                                 ' MultiIndex')

            # NOTE: `group_axis` and `group_axis.get_level_values(level)`
            # are same in this section.
            level = None
            key = group_axis

    # a passed-in Grouper, directly convert
    if isinstance(key, Grouper):
        binner, grouper, obj = key._get_grouper(obj, validate=False)
        if key.key is None:
            return grouper, [], obj
        else:
            return grouper, set([key.key]), obj

    # already have a BaseGrouper, just return it
    elif isinstance(key, BaseGrouper):
        return key, [], obj

    # In the future, a tuple key will always mean an actual key,
    # not an iterable of keys. In the meantime, we attempt to provide
    # a warning. We can assume that the user wanted a list of keys when
    # the key is not in the index. We just have to be careful with
    # unhashble elements of `key`. Any unhashable elements implies that
    # they wanted a list of keys.
    # https://github.com/pandas-dev/pandas/issues/18314
    is_tuple = isinstance(key, tuple)
    all_hashable = is_tuple and is_hashable(key)

    if is_tuple:
        if ((all_hashable and key not in obj and set(key).issubset(obj))
                or not all_hashable):
            # column names ('a', 'b') -> ['a', 'b']
            # arrays like (a, b) -> [a, b]
            msg = ("Interpreting tuple 'by' as a list of keys, rather than "
                   "a single key. Use 'by=[...]' instead of 'by=(...)'. In "
                   "the future, a tuple will always mean a single key.")
            warnings.warn(msg, FutureWarning, stacklevel=5)
            key = list(key)

    if not isinstance(key, list):
        keys = [key]
        match_axis_length = False
    else:
        keys = key
        match_axis_length = len(keys) == len(group_axis)

    # what are we after, exactly?
    any_callable = any(callable(g) or isinstance(g, dict) for g in keys)
    any_groupers = any(isinstance(g, Grouper) for g in keys)
    any_arraylike = any(isinstance(g, (list, tuple, Series, Index, np.ndarray))
                        for g in keys)

    try:
        if isinstance(obj, DataFrame):
            all_in_columns_index = all(g in obj.columns or g in obj.index.names
                                       for g in keys)
        else:
            all_in_columns_index = False
    except Exception:
        all_in_columns_index = False

    if not any_callable and not all_in_columns_index and \
       not any_arraylike and not any_groupers and \
       match_axis_length and level is None:
        keys = [com._asarray_tuplesafe(keys)]

    if isinstance(level, (tuple, list)):
        if key is None:
            keys = [None] * len(level)
        levels = level
    else:
        levels = [level] * len(keys)

    groupings = []
    exclusions = []

    # if the actual grouper should be obj[key]
    def is_in_axis(key):
        if not _is_label_like(key):
            try:
                obj._data.items.get_loc(key)
            except Exception:
                return False

        return True

    # if the grouper is obj[name]
    def is_in_obj(gpr):
        try:
            return id(gpr) == id(obj[gpr.name])
        except Exception:
            return False

    for i, (gpr, level) in enumerate(zip(keys, levels)):

        if is_in_obj(gpr):  # df.groupby(df['name'])
            in_axis, name = True, gpr.name
            exclusions.append(name)

        elif is_in_axis(gpr):  # df.groupby('name')
            if gpr in obj:
                if validate:
                    stacklevel = 5  # Number of stack levels from df.groupby
                    obj._check_label_or_level_ambiguity(
                        gpr, stacklevel=stacklevel)
                in_axis, name, gpr = True, gpr, obj[gpr]
                exclusions.append(name)
            elif obj._is_level_reference(gpr):
                in_axis, name, level, gpr = False, None, gpr, None
            else:
                raise KeyError(gpr)
        elif isinstance(gpr, Grouper) and gpr.key is not None:
            # Add key to exclusions
            exclusions.append(gpr.key)
            in_axis, name = False, None
        else:
            in_axis, name = False, None

        if is_categorical_dtype(gpr) and len(gpr) != obj.shape[axis]:
            raise ValueError(
                ("Length of grouper ({len_gpr}) and axis ({len_axis})"
                 " must be same length"
                 .format(len_gpr=len(gpr), len_axis=obj.shape[axis])))

        # create the Grouping
        # allow us to passing the actual Grouping as the gpr
        ping = Grouping(group_axis,
                        gpr,
                        obj=obj,
                        name=name,
                        level=level,
                        sort=sort,
                        observed=observed,
                        in_axis=in_axis) \
            if not isinstance(gpr, Grouping) else gpr

        groupings.append(ping)

    if len(groupings) == 0:
        raise ValueError('No group keys passed!')

    # create the internals grouper
    grouper = BaseGrouper(group_axis, groupings, sort=sort, mutated=mutated)
    return grouper, exclusions, obj


def _is_label_like(val):
    return (isinstance(val, (compat.string_types, tuple)) or
            (val is not None and is_scalar(val)))


def _convert_grouper(axis, grouper):
    if isinstance(grouper, dict):
        return grouper.get
    elif isinstance(grouper, Series):
        if grouper.index.equals(axis):
            return grouper._values
        else:
            return grouper.reindex(axis)._values
    elif isinstance(grouper, (list, Series, Index, np.ndarray)):
        if len(grouper) != len(axis):
            raise ValueError('Grouper and axis must be same length')
        return grouper
    else:
        return grouper


def _whitelist_method_generator(klass, whitelist):
    """
    Yields all GroupBy member defs for DataFrame/Series names in _whitelist.

    Parameters
    ----------
    klass - class where members are defined.  Should be Series or DataFrame

    whitelist - list of names of klass methods to be constructed

    Returns
    -------
    The generator yields a sequence of strings, each suitable for exec'ing,
    that define implementations of the named methods for DataFrameGroupBy
    or SeriesGroupBy.

    Since we don't want to override methods explicitly defined in the
    base class, any such name is skipped.
    """

    method_wrapper_template = \
        """def %(name)s(%(sig)s) :
    \"""
    %(doc)s
    \"""
    f = %(self)s.__getattr__('%(name)s')
    return f(%(args)s)"""
    property_wrapper_template = \
        """@property
def %(name)s(self) :
    \"""
    %(doc)s
    \"""
    return self.__getattr__('%(name)s')"""
    for name in whitelist:
        # don't override anything that was explicitly defined
        # in the base class
        if hasattr(GroupBy, name):
            continue
        # ugly, but we need the name string itself in the method.
        f = getattr(klass, name)
        doc = f.__doc__
        doc = doc if type(doc) == str else ''
        if isinstance(f, types.MethodType):
            wrapper_template = method_wrapper_template
            decl, args = make_signature(f)
            # pass args by name to f because otherwise
            # GroupBy._make_wrapper won't know whether
            # we passed in an axis parameter.
            args_by_name = ['{0}={0}'.format(arg) for arg in args[1:]]
            params = {'name': name,
                      'doc': doc,
                      'sig': ','.join(decl),
                      'self': args[0],
                      'args': ','.join(args_by_name)}
        else:
            wrapper_template = property_wrapper_template
            params = {'name': name, 'doc': doc}
        yield wrapper_template % params


class SeriesGroupBy(GroupBy):
    #
    # Make class defs of attributes on SeriesGroupBy whitelist
    _apply_whitelist = _series_apply_whitelist
    for _def_str in _whitelist_method_generator(Series,
                                                _series_apply_whitelist):
        exec(_def_str)

    @property
    def _selection_name(self):
        """
        since we are a series, we by definition only have
        a single name, but may be the result of a selection or
        the name of our object
        """
        if self._selection is None:
            return self.obj.name
        else:
            return self._selection

    _agg_doc = dedent("""
    Examples
    --------

    >>> s = Series([1, 2, 3, 4])

    >>> s
    0    1
    1    2
    2    3
    3    4
    dtype: int64

    >>> s.groupby([1, 1, 2, 2]).min()
    1    1
    2    3
    dtype: int64

    >>> s.groupby([1, 1, 2, 2]).agg('min')
    1    1
    2    3
    dtype: int64

    >>> s.groupby([1, 1, 2, 2]).agg(['min', 'max'])
       min  max
    1    1    2
    2    3    4

    See also
    --------
    pandas.Series.groupby.apply
    pandas.Series.groupby.transform
    pandas.Series.aggregate

    """)

    @Appender(_apply_docs['template']
              .format(input='series',
                      examples=_apply_docs['series_examples']))
    def apply(self, func, *args, **kwargs):
        return super(SeriesGroupBy, self).apply(func, *args, **kwargs)

    @Appender(_agg_doc)
    @Appender(_shared_docs['aggregate'] % dict(
        klass='Series',
        versionadded='',
        axis=''))
    def aggregate(self, func_or_funcs, *args, **kwargs):
        _level = kwargs.pop('_level', None)
        if isinstance(func_or_funcs, compat.string_types):
            return getattr(self, func_or_funcs)(*args, **kwargs)

        if isinstance(func_or_funcs, collections.Iterable):
            # Catch instances of lists / tuples
            # but not the class list / tuple itself.
            ret = self._aggregate_multiple_funcs(func_or_funcs,
                                                 (_level or 0) + 1)
        else:
            cyfunc = self._is_cython_func(func_or_funcs)
            if cyfunc and not args and not kwargs:
                return getattr(self, cyfunc)()

            if self.grouper.nkeys > 1:
                return self._python_agg_general(func_or_funcs, *args, **kwargs)

            try:
                return self._python_agg_general(func_or_funcs, *args, **kwargs)
            except Exception:
                result = self._aggregate_named(func_or_funcs, *args, **kwargs)

            index = Index(sorted(result), name=self.grouper.names[0])
            ret = Series(result, index=index)

        if not self.as_index:  # pragma: no cover
            print('Warning, ignoring as_index=True')

        # _level handled at higher
        if not _level and isinstance(ret, dict):
            from pandas import concat
            ret = concat(ret, axis=1)
        return ret

    agg = aggregate

    def _aggregate_multiple_funcs(self, arg, _level):
        if isinstance(arg, dict):

            # show the deprecation, but only if we
            # have not shown a higher level one
            # GH 15931
            if isinstance(self._selected_obj, Series) and _level <= 1:
                warnings.warn(
                    ("using a dict on a Series for aggregation\n"
                     "is deprecated and will be removed in a future "
                     "version"),
                    FutureWarning, stacklevel=3)

            columns = list(arg.keys())
            arg = list(arg.items())
        elif any(isinstance(x, (tuple, list)) for x in arg):
            arg = [(x, x) if not isinstance(x, (tuple, list)) else x
                   for x in arg]

            # indicated column order
            columns = lzip(*arg)[0]
        else:
            # list of functions / function names
            columns = []
            for f in arg:
                if isinstance(f, compat.string_types):
                    columns.append(f)
                else:
                    # protect against callables without names
                    columns.append(com._get_callable_name(f))
            arg = lzip(columns, arg)

        results = {}
        for name, func in arg:
            obj = self
            if name in results:
                raise SpecificationError('Function names must be unique, '
                                         'found multiple named %s' % name)

            # reset the cache so that we
            # only include the named selection
            if name in self._selected_obj:
                obj = copy.copy(obj)
                obj._reset_cache()
                obj._selection = name
            results[name] = obj.aggregate(func)

        if isinstance(list(compat.itervalues(results))[0],
                      DataFrame):

            # let higher level handle
            if _level:
                return results
            return list(compat.itervalues(results))[0]
        return DataFrame(results, columns=columns)

    def _wrap_output(self, output, index, names=None):
        """ common agg/transform wrapping logic """
        output = output[self._selection_name]

        if names is not None:
            return DataFrame(output, index=index, columns=names)
        else:
            name = self._selection_name
            if name is None:
                name = self._selected_obj.name
            return Series(output, index=index, name=name)

    def _wrap_aggregated_output(self, output, names=None):
        return self._wrap_output(output=output,
                                 index=self.grouper.result_index,
                                 names=names)

    def _wrap_transformed_output(self, output, names=None):
        return self._wrap_output(output=output,
                                 index=self.obj.index,
                                 names=names)

    def _wrap_applied_output(self, keys, values, not_indexed_same=False):
        if len(keys) == 0:
            # GH #6265
            return Series([], name=self._selection_name, index=keys)

        def _get_index():
            if self.grouper.nkeys > 1:
                index = MultiIndex.from_tuples(keys, names=self.grouper.names)
            else:
                index = Index(keys, name=self.grouper.names[0])
            return index

        if isinstance(values[0], dict):
            # GH #823
            index = _get_index()
            result = DataFrame(values, index=index).stack()
            result.name = self._selection_name
            return result

        if isinstance(values[0], (Series, dict)):
            return self._concat_objects(keys, values,
                                        not_indexed_same=not_indexed_same)
        elif isinstance(values[0], DataFrame):
            # possible that Series -> DataFrame by applied function
            return self._concat_objects(keys, values,
                                        not_indexed_same=not_indexed_same)
        else:
            # GH #6265
            return Series(values, index=_get_index(),
                          name=self._selection_name)

    def _aggregate_named(self, func, *args, **kwargs):
        result = {}

        for name, group in self:
            group.name = name
            output = func(group, *args, **kwargs)
            if isinstance(output, (Series, Index, np.ndarray)):
                raise Exception('Must produce aggregated value')
            result[name] = self._try_cast(output, group)

        return result

    @Substitution(klass='Series', selected='A.')
    @Appender(_transform_template)
    def transform(self, func, *args, **kwargs):
        func = self._is_cython_func(func) or func

        # if string function
        if isinstance(func, compat.string_types):
            if func in _cython_transforms:
                # cythonized transform
                return getattr(self, func)(*args, **kwargs)
            else:
                # cythonized aggregation and merge
                return self._transform_fast(
                    lambda: getattr(self, func)(*args, **kwargs), func)

        # reg transform
        klass = self._selected_obj.__class__
        results = []
        wrapper = lambda x: func(x, *args, **kwargs)
        for name, group in self:
            object.__setattr__(group, 'name', name)
            res = wrapper(group)

            if hasattr(res, 'values'):
                res = res.values

            indexer = self._get_index(name)
            s = klass(res, indexer)
            results.append(s)

        from pandas.core.reshape.concat import concat
        result = concat(results).sort_index()

        # we will only try to coerce the result type if
        # we have a numeric dtype, as these are *always* udfs
        # the cython take a different path (and casting)
        dtype = self._selected_obj.dtype
        if is_numeric_dtype(dtype):
            result = maybe_downcast_to_dtype(result, dtype)

        result.name = self._selected_obj.name
        result.index = self._selected_obj.index
        return result

    def _transform_fast(self, func, func_nm):
        """
        fast version of transform, only applicable to
        builtin/cythonizable functions
        """
        if isinstance(func, compat.string_types):
            func = getattr(self, func)

        ids, _, ngroup = self.grouper.group_info
        cast = self._transform_should_cast(func_nm)
        out = algorithms.take_1d(func().values, ids)
        if cast:
            out = self._try_cast(out, self.obj)
        return Series(out, index=self.obj.index, name=self.obj.name)

    def filter(self, func, dropna=True, *args, **kwargs):  # noqa
        """
        Return a copy of a Series excluding elements from groups that
        do not satisfy the boolean criterion specified by func.

        Parameters
        ----------
        func : function
            To apply to each group. Should return True or False.
        dropna : Drop groups that do not pass the filter. True by default;
            if False, groups that evaluate False are filled with NaNs.

        Examples
        --------
        >>> import pandas as pd
        >>> df = pd.DataFrame({'A' : ['foo', 'bar', 'foo', 'bar',
        ...                           'foo', 'bar'],
        ...                    'B' : [1, 2, 3, 4, 5, 6],
        ...                    'C' : [2.0, 5., 8., 1., 2., 9.]})
        >>> grouped = df.groupby('A')
        >>> df.groupby('A').B.filter(lambda x: x.mean() > 3.)
        1    2
        3    4
        5    6
        Name: B, dtype: int64

        Returns
        -------
        filtered : Series
        """
        if isinstance(func, compat.string_types):
            wrapper = lambda x: getattr(x, func)(*args, **kwargs)
        else:
            wrapper = lambda x: func(x, *args, **kwargs)

        # Interpret np.nan as False.
        def true_and_notna(x, *args, **kwargs):
            b = wrapper(x, *args, **kwargs)
            return b and notna(b)

        try:
            indices = [self._get_index(name) for name, group in self
                       if true_and_notna(group)]
        except ValueError:
            raise TypeError("the filter must return a boolean result")
        except TypeError:
            raise TypeError("the filter must return a boolean result")

        filtered = self._apply_filter(indices, dropna)
        return filtered

    def nunique(self, dropna=True):
        """ Returns number of unique elements in the group """
        ids, _, _ = self.grouper.group_info

        val = self.obj.get_values()

        try:
            sorter = np.lexsort((val, ids))
        except TypeError:  # catches object dtypes
            assert val.dtype == object, \
                'val.dtype must be object, got %s' % val.dtype
            val, _ = algorithms.factorize(val, sort=False)
            sorter = np.lexsort((val, ids))
            _isna = lambda a: a == -1
        else:
            _isna = isna

        ids, val = ids[sorter], val[sorter]

        # group boundaries are where group ids change
        # unique observations are where sorted values change
        idx = np.r_[0, 1 + np.nonzero(ids[1:] != ids[:-1])[0]]
        inc = np.r_[1, val[1:] != val[:-1]]

        # 1st item of each group is a new unique observation
        mask = _isna(val)
        if dropna:
            inc[idx] = 1
            inc[mask] = 0
        else:
            inc[mask & np.r_[False, mask[:-1]]] = 0
            inc[idx] = 1

        out = np.add.reduceat(inc, idx).astype('int64', copy=False)
        if len(ids):
            # NaN/NaT group exists if the head of ids is -1,
            # so remove it from res and exclude its index from idx
            if ids[0] == -1:
                res = out[1:]
                idx = idx[np.flatnonzero(idx)]
            else:
                res = out
        else:
            res = out[1:]
        ri = self.grouper.result_index

        # we might have duplications among the bins
        if len(res) != len(ri):
            res, out = np.zeros(len(ri), dtype=out.dtype), res
            res[ids[idx]] = out

        return Series(res,
                      index=ri,
                      name=self._selection_name)

    @Appender(Series.describe.__doc__)
    def describe(self, **kwargs):
        result = self.apply(lambda x: x.describe(**kwargs))
        if self.axis == 1:
            return result.T
        return result.unstack()

    def value_counts(self, normalize=False, sort=True, ascending=False,
                     bins=None, dropna=True):

        from pandas.core.reshape.tile import cut
        from pandas.core.reshape.merge import _get_join_indexers

        if bins is not None and not np.iterable(bins):
            # scalar bins cannot be done at top level
            # in a backward compatible way
            return self.apply(Series.value_counts,
                              normalize=normalize,
                              sort=sort,
                              ascending=ascending,
                              bins=bins)

        ids, _, _ = self.grouper.group_info
        val = self.obj.get_values()

        # groupby removes null keys from groupings
        mask = ids != -1
        ids, val = ids[mask], val[mask]

        if bins is None:
            lab, lev = algorithms.factorize(val, sort=True)
            llab = lambda lab, inc: lab[inc]
        else:

            # lab is a Categorical with categories an IntervalIndex
            lab = cut(Series(val), bins, include_lowest=True)
            lev = lab.cat.categories
            lab = lev.take(lab.cat.codes)
            llab = lambda lab, inc: lab[inc]._multiindex.labels[-1]

        if is_interval_dtype(lab):
            # TODO: should we do this inside II?
            sorter = np.lexsort((lab.left, lab.right, ids))
        else:
            sorter = np.lexsort((lab, ids))

        ids, lab = ids[sorter], lab[sorter]

        # group boundaries are where group ids change
        idx = np.r_[0, 1 + np.nonzero(ids[1:] != ids[:-1])[0]]

        # new values are where sorted labels change
        lchanges = llab(lab, slice(1, None)) != llab(lab, slice(None, -1))
        inc = np.r_[True, lchanges]
        inc[idx] = True  # group boundaries are also new values
        out = np.diff(np.nonzero(np.r_[inc, True])[0])  # value counts

        # num. of times each group should be repeated
        rep = partial(np.repeat, repeats=np.add.reduceat(inc, idx))

        # multi-index components
        labels = list(map(rep, self.grouper.recons_labels)) + [llab(lab, inc)]
        levels = [ping.group_index for ping in self.grouper.groupings] + [lev]
        names = self.grouper.names + [self._selection_name]

        if dropna:
            mask = labels[-1] != -1
            if mask.all():
                dropna = False
            else:
                out, labels = out[mask], [label[mask] for label in labels]

        if normalize:
            out = out.astype('float')
            d = np.diff(np.r_[idx, len(ids)])
            if dropna:
                m = ids[lab == -1]
                np.add.at(d, m, -1)
                acc = rep(d)[mask]
            else:
                acc = rep(d)
            out /= acc

        if sort and bins is None:
            cat = ids[inc][mask] if dropna else ids[inc]
            sorter = np.lexsort((out if ascending else -out, cat))
            out, labels[-1] = out[sorter], labels[-1][sorter]

        if bins is None:
            mi = MultiIndex(levels=levels, labels=labels, names=names,
                            verify_integrity=False)

            if is_integer_dtype(out):
                out = _ensure_int64(out)
            return Series(out, index=mi, name=self._selection_name)

        # for compat. with libgroupby.value_counts need to ensure every
        # bin is present at every index level, null filled with zeros
        diff = np.zeros(len(out), dtype='bool')
        for lab in labels[:-1]:
            diff |= np.r_[True, lab[1:] != lab[:-1]]

        ncat, nbin = diff.sum(), len(levels[-1])

        left = [np.repeat(np.arange(ncat), nbin),
                np.tile(np.arange(nbin), ncat)]

        right = [diff.cumsum() - 1, labels[-1]]

        _, idx = _get_join_indexers(left, right, sort=False, how='left')
        out = np.where(idx != -1, out[idx], 0)

        if sort:
            sorter = np.lexsort((out if ascending else -out, left[0]))
            out, left[-1] = out[sorter], left[-1][sorter]

        # build the multi-index w/ full levels
        labels = list(map(lambda lab: np.repeat(lab[diff], nbin), labels[:-1]))
        labels.append(left[-1])

        mi = MultiIndex(levels=levels, labels=labels, names=names,
                        verify_integrity=False)

        if is_integer_dtype(out):
            out = _ensure_int64(out)
        return Series(out, index=mi, name=self._selection_name)

    def count(self):
        """ Compute count of group, excluding missing values """
        ids, _, ngroups = self.grouper.group_info
        val = self.obj.get_values()

        mask = (ids != -1) & ~isna(val)
        ids = _ensure_platform_int(ids)
        out = np.bincount(ids[mask], minlength=ngroups or 0)

        return Series(out,
                      index=self.grouper.result_index,
                      name=self._selection_name,
                      dtype='int64')

    def _apply_to_column_groupbys(self, func):
        """ return a pass thru """
        return func(self)

    def pct_change(self, periods=1, fill_method='pad', limit=None, freq=None):
        """Calculate percent change of each value to previous entry in group"""
        filled = getattr(self, fill_method)(limit=limit)
        shifted = filled.shift(periods=periods, freq=freq)

        return (filled / shifted) - 1


class NDFrameGroupBy(GroupBy):

    def _iterate_slices(self):
        if self.axis == 0:
            # kludge
            if self._selection is None:
                slice_axis = self.obj.columns
            else:
                slice_axis = self._selection_list
            slicer = lambda x: self.obj[x]
        else:
            slice_axis = self.obj.index
            slicer = self.obj.xs

        for val in slice_axis:
            if val in self.exclusions:
                continue
            yield val, slicer(val)

    def _cython_agg_general(self, how, alt=None, numeric_only=True,
                            min_count=-1):
        new_items, new_blocks = self._cython_agg_blocks(
            how, alt=alt, numeric_only=numeric_only, min_count=min_count)
        return self._wrap_agged_blocks(new_items, new_blocks)

    def _wrap_agged_blocks(self, items, blocks):
        obj = self._obj_with_exclusions

        new_axes = list(obj._data.axes)

        # more kludge
        if self.axis == 0:
            new_axes[0], new_axes[1] = new_axes[1], self.grouper.result_index
        else:
            new_axes[self.axis] = self.grouper.result_index

        # Make sure block manager integrity check passes.
        assert new_axes[0].equals(items)
        new_axes[0] = items

        mgr = BlockManager(blocks, new_axes)

        new_obj = type(obj)(mgr)

        return self._post_process_cython_aggregate(new_obj)

    _block_agg_axis = 0

    def _cython_agg_blocks(self, how, alt=None, numeric_only=True,
                           min_count=-1):
        # TODO: the actual managing of mgr_locs is a PITA
        # here, it should happen via BlockManager.combine

        data, agg_axis = self._get_data_to_aggregate()

        if numeric_only:
            data = data.get_numeric_data(copy=False)

        new_blocks = []
        new_items = []
        deleted_items = []
        for block in data.blocks:

            locs = block.mgr_locs.as_array
            try:
                result, _ = self.grouper.aggregate(
                    block.values, how, axis=agg_axis, min_count=min_count)
            except NotImplementedError:
                # generally if we have numeric_only=False
                # and non-applicable functions
                # try to python agg

                if alt is None:
                    # we cannot perform the operation
                    # in an alternate way, exclude the block
                    deleted_items.append(locs)
                    continue

                # call our grouper again with only this block
                obj = self.obj[data.items[locs]]
                s = groupby(obj, self.grouper)
                result = s.aggregate(lambda x: alt(x, axis=self.axis))
                newb = result._data.blocks[0]

            finally:

                # see if we can cast the block back to the original dtype
                result = block._try_coerce_and_cast_result(result)
                newb = block.make_block(result)

            new_items.append(locs)
            new_blocks.append(newb)

        if len(new_blocks) == 0:
            raise DataError('No numeric types to aggregate')

        # reset the locs in the blocks to correspond to our
        # current ordering
        indexer = np.concatenate(new_items)
        new_items = data.items.take(np.sort(indexer))

        if len(deleted_items):

            # we need to adjust the indexer to account for the
            # items we have removed
            # really should be done in internals :<

            deleted = np.concatenate(deleted_items)
            ai = np.arange(len(data))
            mask = np.zeros(len(data))
            mask[deleted] = 1
            indexer = (ai - mask.cumsum())[indexer]

        offset = 0
        for b in new_blocks:
            loc = len(b.mgr_locs)
            b.mgr_locs = indexer[offset:(offset + loc)]
            offset += loc

        return new_items, new_blocks

    def _get_data_to_aggregate(self):
        obj = self._obj_with_exclusions
        if self.axis == 0:
            return obj.swapaxes(0, 1)._data, 1
        else:
            return obj._data, self.axis

    def _post_process_cython_aggregate(self, obj):
        # undoing kludge from below
        if self.axis == 0:
            obj = obj.swapaxes(0, 1)
        return obj

    def aggregate(self, arg, *args, **kwargs):

        _level = kwargs.pop('_level', None)
        result, how = self._aggregate(arg, _level=_level, *args, **kwargs)
        if how is None:
            return result

        if result is None:

            # grouper specific aggregations
            if self.grouper.nkeys > 1:
                return self._python_agg_general(arg, *args, **kwargs)
            else:

                # try to treat as if we are passing a list
                try:
                    assert not args and not kwargs
                    result = self._aggregate_multiple_funcs(
                        [arg], _level=_level, _axis=self.axis)
                    result.columns = Index(
                        result.columns.levels[0],
                        name=self._selected_obj.columns.name)
                except Exception:
                    result = self._aggregate_generic(arg, *args, **kwargs)

        if not self.as_index:
            self._insert_inaxis_grouper_inplace(result)
            result.index = np.arange(len(result))

        return result._convert(datetime=True)

    agg = aggregate

    def _aggregate_generic(self, func, *args, **kwargs):
        if self.grouper.nkeys != 1:
            raise AssertionError('Number of keys must be 1')

        axis = self.axis
        obj = self._obj_with_exclusions

        result = {}
        if axis != obj._info_axis_number:
            try:
                for name, data in self:
                    result[name] = self._try_cast(func(data, *args, **kwargs),
                                                  data)
            except Exception:
                return self._aggregate_item_by_item(func, *args, **kwargs)
        else:
            for name in self.indices:
                try:
                    data = self.get_group(name, obj=obj)
                    result[name] = self._try_cast(func(data, *args, **kwargs),
                                                  data)
                except Exception:
                    wrapper = lambda x: func(x, *args, **kwargs)
                    result[name] = data.apply(wrapper, axis=axis)

        return self._wrap_generic_output(result, obj)

    def _wrap_aggregated_output(self, output, names=None):
        raise com.AbstractMethodError(self)

    def _aggregate_item_by_item(self, func, *args, **kwargs):
        # only for axis==0

        obj = self._obj_with_exclusions
        result = {}
        cannot_agg = []
        errors = None
        for item in obj:
            try:
                data = obj[item]
                colg = SeriesGroupBy(data, selection=item,
                                     grouper=self.grouper)
                result[item] = self._try_cast(
                    colg.aggregate(func, *args, **kwargs), data)
            except ValueError:
                cannot_agg.append(item)
                continue
            except TypeError as e:
                cannot_agg.append(item)
                errors = e
                continue

        result_columns = obj.columns
        if cannot_agg:
            result_columns = result_columns.drop(cannot_agg)

            # GH6337
            if not len(result_columns) and errors is not None:
                raise errors

        return DataFrame(result, columns=result_columns)

    def _decide_output_index(self, output, labels):
        if len(output) == len(labels):
            output_keys = labels
        else:
            output_keys = sorted(output)
            try:
                output_keys.sort()
            except Exception:  # pragma: no cover
                pass

            if isinstance(labels, MultiIndex):
                output_keys = MultiIndex.from_tuples(output_keys,
                                                     names=labels.names)

        return output_keys

    def _wrap_applied_output(self, keys, values, not_indexed_same=False):
        from pandas.core.index import _all_indexes_same
        from pandas.core.tools.numeric import to_numeric

        if len(keys) == 0:
            return DataFrame(index=keys)

        key_names = self.grouper.names

        # GH12824.
        def first_not_none(values):
            try:
                return next(com._not_none(*values))
            except StopIteration:
                return None

        v = first_not_none(values)

        if v is None:
            # GH9684. If all values are None, then this will throw an error.
            # We'd prefer it return an empty dataframe.
            return DataFrame()
        elif isinstance(v, DataFrame):
            return self._concat_objects(keys, values,
                                        not_indexed_same=not_indexed_same)
        elif self.grouper.groupings is not None:
            if len(self.grouper.groupings) > 1:
                key_index = self.grouper.result_index

            else:
                ping = self.grouper.groupings[0]
                if len(keys) == ping.ngroups:
                    key_index = ping.group_index
                    key_index.name = key_names[0]

                    key_lookup = Index(keys)
                    indexer = key_lookup.get_indexer(key_index)

                    # reorder the values
                    values = [values[i] for i in indexer]
                else:

                    key_index = Index(keys, name=key_names[0])

                # don't use the key indexer
                if not self.as_index:
                    key_index = None

            # make Nones an empty object
            v = first_not_none(values)
            if v is None:
                return DataFrame()
            elif isinstance(v, NDFrame):
                values = [
                    x if x is not None else
                    v._constructor(**v._construct_axes_dict())
                    for x in values
                ]

            v = values[0]

            if isinstance(v, (np.ndarray, Index, Series)):
                if isinstance(v, Series):
                    applied_index = self._selected_obj._get_axis(self.axis)
                    all_indexed_same = _all_indexes_same([
                        x.index for x in values
                    ])
                    singular_series = (len(values) == 1 and
                                       applied_index.nlevels == 1)

                    # GH3596
                    # provide a reduction (Frame -> Series) if groups are
                    # unique
                    if self.squeeze:

                        # assign the name to this series
                        if singular_series:
                            values[0].name = keys[0]

                            # GH2893
                            # we have series in the values array, we want to
                            # produce a series:
                            # if any of the sub-series are not indexed the same
                            # OR we don't have a multi-index and we have only a
                            # single values
                            return self._concat_objects(
                                keys, values, not_indexed_same=not_indexed_same
                            )

                        # still a series
                        # path added as of GH 5545
                        elif all_indexed_same:
                            from pandas.core.reshape.concat import concat
                            return concat(values)

                    if not all_indexed_same:
                        # GH 8467
                        return self._concat_objects(
                            keys, values, not_indexed_same=True,
                        )

                try:
                    if self.axis == 0:
                        # GH6124 if the list of Series have a consistent name,
                        # then propagate that name to the result.
                        index = v.index.copy()
                        if index.name is None:
                            # Only propagate the series name to the result
                            # if all series have a consistent name.  If the
                            # series do not have a consistent name, do
                            # nothing.
                            names = {v.name for v in values}
                            if len(names) == 1:
                                index.name = list(names)[0]

                        # normally use vstack as its faster than concat
                        # and if we have mi-columns
                        if (isinstance(v.index, MultiIndex) or
                                key_index is None or
                                isinstance(key_index, MultiIndex)):
                            stacked_values = np.vstack(map(np.asarray, values))
                            result = DataFrame(stacked_values, index=key_index,
                                               columns=index)
                        else:
                            # GH5788 instead of stacking; concat gets the
                            # dtypes correct
                            from pandas.core.reshape.concat import concat
                            result = concat(values, keys=key_index,
                                            names=key_index.names,
                                            axis=self.axis).unstack()
                            result.columns = index
                    else:
                        stacked_values = np.vstack(map(np.asarray, values))
                        result = DataFrame(stacked_values.T, index=v.index,
                                           columns=key_index)

                except (ValueError, AttributeError):
                    # GH1738: values is list of arrays of unequal lengths fall
                    # through to the outer else caluse
                    return Series(values, index=key_index,
                                  name=self._selection_name)

                # if we have date/time like in the original, then coerce dates
                # as we are stacking can easily have object dtypes here
                so = self._selected_obj
                if (so.ndim == 2 and so.dtypes.apply(is_datetimelike).any()):
                    result = result.apply(
                        lambda x: to_numeric(x, errors='ignore'))
                    date_cols = self._selected_obj.select_dtypes(
                        include=['datetime', 'timedelta']).columns
                    date_cols = date_cols.intersection(result.columns)
                    result[date_cols] = (result[date_cols]
                                         ._convert(datetime=True,
                                                   coerce=True))
                else:
                    result = result._convert(datetime=True)

                return self._reindex_output(result)

            # values are not series or array-like but scalars
            else:
                # only coerce dates if we find at least 1 datetime
                coerce = any(isinstance(x, Timestamp) for x in values)
                # self._selection_name not passed through to Series as the
                # result should not take the name of original selection
                # of columns
                return (Series(values, index=key_index)
                        ._convert(datetime=True,
                                  coerce=coerce))

        else:
            # Handle cases like BinGrouper
            return self._concat_objects(keys, values,
                                        not_indexed_same=not_indexed_same)

    def _transform_general(self, func, *args, **kwargs):
        from pandas.core.reshape.concat import concat

        applied = []
        obj = self._obj_with_exclusions
        gen = self.grouper.get_iterator(obj, axis=self.axis)
        fast_path, slow_path = self._define_paths(func, *args, **kwargs)

        path = None
        for name, group in gen:
            object.__setattr__(group, 'name', name)

            if path is None:
                # Try slow path and fast path.
                try:
                    path, res = self._choose_path(fast_path, slow_path, group)
                except TypeError:
                    return self._transform_item_by_item(obj, fast_path)
                except ValueError:
                    msg = 'transform must return a scalar value for each group'
                    raise ValueError(msg)
            else:
                res = path(group)

            if isinstance(res, Series):

                # we need to broadcast across the
                # other dimension; this will preserve dtypes
                # GH14457
                if not np.prod(group.shape):
                    continue
                elif res.index.is_(obj.index):
                    r = concat([res] * len(group.columns), axis=1)
                    r.columns = group.columns
                    r.index = group.index
                else:
                    r = DataFrame(
                        np.concatenate([res.values] * len(group.index)
                                       ).reshape(group.shape),
                        columns=group.columns, index=group.index)

                applied.append(r)
            else:
                applied.append(res)

        concat_index = obj.columns if self.axis == 0 else obj.index
        concatenated = concat(applied, join_axes=[concat_index],
                              axis=self.axis, verify_integrity=False)
        return self._set_result_index_ordered(concatenated)

    @Substitution(klass='DataFrame', selected='')
    @Appender(_transform_template)
    def transform(self, func, *args, **kwargs):

        # optimized transforms
        func = self._is_cython_func(func) or func
        if isinstance(func, compat.string_types):
            if func in _cython_transforms:
                # cythonized transform
                return getattr(self, func)(*args, **kwargs)
            else:
                # cythonized aggregation and merge
                result = getattr(self, func)(*args, **kwargs)
        else:
            return self._transform_general(func, *args, **kwargs)

        # a reduction transform
        if not isinstance(result, DataFrame):
            return self._transform_general(func, *args, **kwargs)

        obj = self._obj_with_exclusions

        # nuiscance columns
        if not result.columns.equals(obj.columns):
            return self._transform_general(func, *args, **kwargs)

        return self._transform_fast(result, obj, func)

    def _transform_fast(self, result, obj, func_nm):
        """
        Fast transform path for aggregations
        """
        # if there were groups with no observations (Categorical only?)
        # try casting data to original dtype
        cast = self._transform_should_cast(func_nm)

        # for each col, reshape to to size of original frame
        # by take operation
        ids, _, ngroup = self.grouper.group_info
        output = []
        for i, _ in enumerate(result.columns):
            res = algorithms.take_1d(result.iloc[:, i].values, ids)
            if cast:
                res = self._try_cast(res, obj.iloc[:, i])
            output.append(res)

        return DataFrame._from_arrays(output, columns=result.columns,
                                      index=obj.index)

    def _define_paths(self, func, *args, **kwargs):
        if isinstance(func, compat.string_types):
            fast_path = lambda group: getattr(group, func)(*args, **kwargs)
            slow_path = lambda group: group.apply(
                lambda x: getattr(x, func)(*args, **kwargs), axis=self.axis)
        else:
            fast_path = lambda group: func(group, *args, **kwargs)
            slow_path = lambda group: group.apply(
                lambda x: func(x, *args, **kwargs), axis=self.axis)
        return fast_path, slow_path

    def _choose_path(self, fast_path, slow_path, group):
        path = slow_path
        res = slow_path(group)

        # if we make it here, test if we can use the fast path
        try:
            res_fast = fast_path(group)

            # compare that we get the same results
            if res.shape == res_fast.shape:
                res_r = res.values.ravel()
                res_fast_r = res_fast.values.ravel()
                mask = notna(res_r)
            if (res_r[mask] == res_fast_r[mask]).all():
                path = fast_path

        except Exception:
            pass
        return path, res

    def _transform_item_by_item(self, obj, wrapper):
        # iterate through columns
        output = {}
        inds = []
        for i, col in enumerate(obj):
            try:
                output[col] = self[col].transform(wrapper)
                inds.append(i)
            except Exception:
                pass

        if len(output) == 0:  # pragma: no cover
            raise TypeError('Transform function invalid for data types')

        columns = obj.columns
        if len(output) < len(obj.columns):
            columns = columns.take(inds)

        return DataFrame(output, index=obj.index, columns=columns)

    def filter(self, func, dropna=True, *args, **kwargs):  # noqa
        """
        Return a copy of a DataFrame excluding elements from groups that
        do not satisfy the boolean criterion specified by func.

        Parameters
        ----------
        f : function
            Function to apply to each subframe. Should return True or False.
        dropna : Drop groups that do not pass the filter. True by default;
            if False, groups that evaluate False are filled with NaNs.

        Notes
        -----
        Each subframe is endowed the attribute 'name' in case you need to know
        which group you are working on.

        Examples
        --------
        >>> import pandas as pd
        >>> df = pd.DataFrame({'A' : ['foo', 'bar', 'foo', 'bar',
        ...                           'foo', 'bar'],
        ...                    'B' : [1, 2, 3, 4, 5, 6],
        ...                    'C' : [2.0, 5., 8., 1., 2., 9.]})
        >>> grouped = df.groupby('A')
        >>> grouped.filter(lambda x: x['B'].mean() > 3.)
             A  B    C
        1  bar  2  5.0
        3  bar  4  1.0
        5  bar  6  9.0

        Returns
        -------
        filtered : DataFrame
        """

        indices = []

        obj = self._selected_obj
        gen = self.grouper.get_iterator(obj, axis=self.axis)

        for name, group in gen:
            object.__setattr__(group, 'name', name)

            res = func(group, *args, **kwargs)

            try:
                res = res.squeeze()
            except AttributeError:  # allow e.g., scalars and frames to pass
                pass

            # interpret the result of the filter
            if is_bool(res) or (is_scalar(res) and isna(res)):
                if res and notna(res):
                    indices.append(self._get_index(name))
            else:
                # non scalars aren't allowed
                raise TypeError("filter function returned a %s, "
                                "but expected a scalar bool" %
                                type(res).__name__)

        return self._apply_filter(indices, dropna)


class DataFrameGroupBy(NDFrameGroupBy):
    _apply_whitelist = _dataframe_apply_whitelist
    #
    # Make class defs of attributes on DataFrameGroupBy whitelist.
    for _def_str in _whitelist_method_generator(DataFrame, _apply_whitelist):
        exec(_def_str)

    _block_agg_axis = 1

    _agg_doc = dedent("""
    Examples
    --------

    >>> df = pd.DataFrame({'A': [1, 1, 2, 2],
    ...                    'B': [1, 2, 3, 4],
    ...                    'C': np.random.randn(4)})

    >>> df
       A  B         C
    0  1  1  0.362838
    1  1  2  0.227877
    2  2  3  1.267767
    3  2  4 -0.562860

    The aggregation is for each column.

    >>> df.groupby('A').agg('min')
       B         C
    A
    1  1  0.227877
    2  3 -0.562860

    Multiple aggregations

    >>> df.groupby('A').agg(['min', 'max'])
        B             C
      min max       min       max
    A
    1   1   2  0.227877  0.362838
    2   3   4 -0.562860  1.267767

    Select a column for aggregation

    >>> df.groupby('A').B.agg(['min', 'max'])
       min  max
    A
    1    1    2
    2    3    4

    Different aggregations per column

    >>> df.groupby('A').agg({'B': ['min', 'max'], 'C': 'sum'})
        B             C
      min max       sum
    A
    1   1   2  0.590716
    2   3   4  0.704907

    See also
    --------
    pandas.DataFrame.groupby.apply
    pandas.DataFrame.groupby.transform
    pandas.DataFrame.aggregate

    """)

    @Appender(_agg_doc)
    @Appender(_shared_docs['aggregate'] % dict(
        klass='DataFrame',
        versionadded='',
        axis=''))
    def aggregate(self, arg, *args, **kwargs):
        return super(DataFrameGroupBy, self).aggregate(arg, *args, **kwargs)

    agg = aggregate

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

        if ndim == 2:
            if subset is None:
                subset = self.obj
            return DataFrameGroupBy(subset, self.grouper, selection=key,
                                    grouper=self.grouper,
                                    exclusions=self.exclusions,
                                    as_index=self.as_index)
        elif ndim == 1:
            if subset is None:
                subset = self.obj[key]
            return SeriesGroupBy(subset, selection=key,
                                 grouper=self.grouper)

        raise AssertionError("invalid ndim for _gotitem")

    def _wrap_generic_output(self, result, obj):
        result_index = self.grouper.levels[0]

        if self.axis == 0:
            return DataFrame(result, index=obj.columns,
                             columns=result_index).T
        else:
            return DataFrame(result, index=obj.index,
                             columns=result_index)

    def _get_data_to_aggregate(self):
        obj = self._obj_with_exclusions
        if self.axis == 1:
            return obj.T._data, 1
        else:
            return obj._data, 1

    def _insert_inaxis_grouper_inplace(self, result):
        # zip in reverse so we can always insert at loc 0
        izip = zip(* map(reversed, (
            self.grouper.names,
            self.grouper.get_group_levels(),
            [grp.in_axis for grp in self.grouper.groupings])))

        for name, lev, in_axis in izip:
            if in_axis:
                result.insert(0, name, lev)

    def _wrap_aggregated_output(self, output, names=None):
        agg_axis = 0 if self.axis == 1 else 1
        agg_labels = self._obj_with_exclusions._get_axis(agg_axis)

        output_keys = self._decide_output_index(output, agg_labels)

        if not self.as_index:
            result = DataFrame(output, columns=output_keys)
            self._insert_inaxis_grouper_inplace(result)
            result = result._consolidate()
        else:
            index = self.grouper.result_index
            result = DataFrame(output, index=index, columns=output_keys)

        if self.axis == 1:
            result = result.T

        return self._reindex_output(result)._convert(datetime=True)

    def _wrap_transformed_output(self, output, names=None):
        return DataFrame(output, index=self.obj.index)

    def _wrap_agged_blocks(self, items, blocks):
        if not self.as_index:
            index = np.arange(blocks[0].values.shape[1])
            mgr = BlockManager(blocks, [items, index])
            result = DataFrame(mgr)

            self._insert_inaxis_grouper_inplace(result)
            result = result._consolidate()
        else:
            index = self.grouper.result_index
            mgr = BlockManager(blocks, [items, index])
            result = DataFrame(mgr)

        if self.axis == 1:
            result = result.T

        return self._reindex_output(result)._convert(datetime=True)

    def _reindex_output(self, result):
        """
        If we have categorical groupers, then we want to make sure that
        we have a fully reindex-output to the levels. These may have not
        participated in the groupings (e.g. may have all been
        nan groups);

        This can re-expand the output space
        """

        # we need to re-expand the output space to accomodate all values
        # whether observed or not in the cartesian product of our groupes
        groupings = self.grouper.groupings
        if groupings is None:
            return result
        elif len(groupings) == 1:
            return result

        # if we only care about the observed values
        # we are done
        elif self.observed:
            return result

        # reindexing only applies to a Categorical grouper
        elif not any(isinstance(ping.grouper, (Categorical, CategoricalIndex))
                     for ping in groupings):
            return result

        levels_list = [ping.group_index for ping in groupings]
        index, _ = MultiIndex.from_product(
            levels_list, names=self.grouper.names).sortlevel()

        if self.as_index:
            d = {self.obj._get_axis_name(self.axis): index, 'copy': False}
            return result.reindex(**d)

        # GH 13204
        # Here, the categorical in-axis groupers, which need to be fully
        # expanded, are columns in `result`. An idea is to do:
        # result = result.set_index(self.grouper.names)
        #                .reindex(index).reset_index()
        # but special care has to be taken because of possible not-in-axis
        # groupers.
        # So, we manually select and drop the in-axis grouper columns,
        # reindex `result`, and then reset the in-axis grouper columns.

        # Select in-axis groupers
        in_axis_grps = [(i, ping.name) for (i, ping)
                        in enumerate(groupings) if ping.in_axis]
        g_nums, g_names = zip(*in_axis_grps)

        result = result.drop(labels=list(g_names), axis=1)

        # Set a temp index and reindex (possibly expanding)
        result = result.set_index(self.grouper.result_index
                                  ).reindex(index, copy=False)

        # Reset in-axis grouper columns
        # (using level numbers `g_nums` because level names may not be unique)
        result = result.reset_index(level=g_nums)

        return result.reset_index(drop=True)

    def _iterate_column_groupbys(self):
        for i, colname in enumerate(self._selected_obj.columns):
            yield colname, SeriesGroupBy(self._selected_obj.iloc[:, i],
                                         selection=colname,
                                         grouper=self.grouper,
                                         exclusions=self.exclusions)

    def _apply_to_column_groupbys(self, func):
        from pandas.core.reshape.concat import concat
        return concat(
            (func(col_groupby) for _, col_groupby
             in self._iterate_column_groupbys()),
            keys=self._selected_obj.columns, axis=1)

    def _fill(self, direction, limit=None):
        """Overridden method to join grouped columns in output"""
        res = super(DataFrameGroupBy, self)._fill(direction, limit=limit)
        output = collections.OrderedDict(
            (grp.name, grp.grouper) for grp in self.grouper.groupings)

        from pandas import concat
        return concat((self._wrap_transformed_output(output), res), axis=1)

    def count(self):
        """ Compute count of group, excluding missing values """
        from pandas.core.dtypes.missing import _isna_ndarraylike as isna

        data, _ = self._get_data_to_aggregate()
        ids, _, ngroups = self.grouper.group_info
        mask = ids != -1

        val = ((mask & ~isna(np.atleast_2d(blk.get_values())))
               for blk in data.blocks)
        loc = (blk.mgr_locs for blk in data.blocks)

        counter = partial(count_level_2d, labels=ids, max_bin=ngroups, axis=1)
        blk = map(make_block, map(counter, val), loc)

        return self._wrap_agged_blocks(data.items, list(blk))

    def nunique(self, dropna=True):
        """
        Return DataFrame with number of distinct observations per group for
        each column.

        .. versionadded:: 0.20.0

        Parameters
        ----------
        dropna : boolean, default True
            Don't include NaN in the counts.

        Returns
        -------
        nunique: DataFrame

        Examples
        --------
        >>> df = pd.DataFrame({'id': ['spam', 'egg', 'egg', 'spam',
        ...                           'ham', 'ham'],
        ...                    'value1': [1, 5, 5, 2, 5, 5],
        ...                    'value2': list('abbaxy')})
        >>> df
             id  value1 value2
        0  spam       1      a
        1   egg       5      b
        2   egg       5      b
        3  spam       2      a
        4   ham       5      x
        5   ham       5      y

        >>> df.groupby('id').nunique()
            id  value1  value2
        id
        egg    1       1       1
        ham    1       1       2
        spam   1       2       1

        # check for rows with the same id but conflicting values
        >>> df.groupby('id').filter(lambda g: (g.nunique() > 1).any())
             id  value1 value2
        0  spam       1      a
        3  spam       2      a
        4   ham       5      x
        5   ham       5      y
        """

        obj = self._selected_obj

        def groupby_series(obj, col=None):
            return SeriesGroupBy(obj,
                                 selection=col,
                                 grouper=self.grouper).nunique(dropna=dropna)

        if isinstance(obj, Series):
            results = groupby_series(obj)
        else:
            from pandas.core.reshape.concat import concat
            results = [groupby_series(obj[col], col) for col in obj.columns]
            results = concat(results, axis=1)

        if not self.as_index:
            results.index = com._default_index(len(results))
        return results

    boxplot = boxplot_frame_groupby


class PanelGroupBy(NDFrameGroupBy):

    def aggregate(self, arg, *args, **kwargs):
        return super(PanelGroupBy, self).aggregate(arg, *args, **kwargs)

    agg = aggregate

    def _iterate_slices(self):
        if self.axis == 0:
            # kludge
            if self._selection is None:
                slice_axis = self._selected_obj.items
            else:
                slice_axis = self._selection_list
            slicer = lambda x: self._selected_obj[x]
        else:
            raise NotImplementedError("axis other than 0 is not supported")

        for val in slice_axis:
            if val in self.exclusions:
                continue

            yield val, slicer(val)

    def aggregate(self, arg, *args, **kwargs):
        """
        Aggregate using input function or dict of {column -> function}

        Parameters
        ----------
        arg : function or dict
            Function to use for aggregating groups. If a function, must either
            work when passed a Panel or when passed to Panel.apply. If
            pass a dict, the keys must be DataFrame column names

        Returns
        -------
        aggregated : Panel
        """
        if isinstance(arg, compat.string_types):
            return getattr(self, arg)(*args, **kwargs)

        return self._aggregate_generic(arg, *args, **kwargs)

    def _wrap_generic_output(self, result, obj):
        if self.axis == 0:
            new_axes = list(obj.axes)
            new_axes[0] = self.grouper.result_index
        elif self.axis == 1:
            x, y, z = obj.axes
            new_axes = [self.grouper.result_index, z, x]
        else:
            x, y, z = obj.axes
            new_axes = [self.grouper.result_index, y, x]

        result = Panel._from_axes(result, new_axes)

        if self.axis == 1:
            result = result.swapaxes(0, 1).swapaxes(0, 2)
        elif self.axis == 2:
            result = result.swapaxes(0, 2)

        return result

    def _aggregate_item_by_item(self, func, *args, **kwargs):
        obj = self._obj_with_exclusions
        result = {}

        if self.axis > 0:
            for item in obj:
                try:
                    itemg = DataFrameGroupBy(obj[item],
                                             axis=self.axis - 1,
                                             grouper=self.grouper)
                    result[item] = itemg.aggregate(func, *args, **kwargs)
                except (ValueError, TypeError):
                    raise
            new_axes = list(obj.axes)
            new_axes[self.axis] = self.grouper.result_index
            return Panel._from_axes(result, new_axes)
        else:
            raise ValueError("axis value must be greater than 0")

    def _wrap_aggregated_output(self, output, names=None):
        raise com.AbstractMethodError(self)


# ----------------------------------------------------------------------
# Splitting / application


class DataSplitter(object):

    def __init__(self, data, labels, ngroups, axis=0):
        self.data = data
        self.labels = _ensure_int64(labels)
        self.ngroups = ngroups

        self.axis = axis

    @cache_readonly
    def slabels(self):
        # Sorted labels
        return algorithms.take_nd(self.labels, self.sort_idx, allow_fill=False)

    @cache_readonly
    def sort_idx(self):
        # Counting sort indexer
        return get_group_index_sorter(self.labels, self.ngroups)

    def __iter__(self):
        sdata = self._get_sorted_data()

        if self.ngroups == 0:
            # we are inside a generator, rather than raise StopIteration
            # we merely return signal the end
            return

        starts, ends = lib.generate_slices(self.slabels, self.ngroups)

        for i, (start, end) in enumerate(zip(starts, ends)):
            # Since I'm now compressing the group ids, it's now not "possible"
            # to produce empty slices because such groups would not be observed
            # in the data
            # if start >= end:
            #     raise AssertionError('Start %s must be less than end %s'
            #                          % (str(start), str(end)))
            yield i, self._chop(sdata, slice(start, end))

    def _get_sorted_data(self):
        return self.data._take(self.sort_idx, axis=self.axis)

    def _chop(self, sdata, slice_obj):
        return sdata.iloc[slice_obj]

    def apply(self, f):
        raise com.AbstractMethodError(self)


class SeriesSplitter(DataSplitter):

    def _chop(self, sdata, slice_obj):
        return sdata._get_values(slice_obj).to_dense()


class FrameSplitter(DataSplitter):

    def __init__(self, data, labels, ngroups, axis=0):
        super(FrameSplitter, self).__init__(data, labels, ngroups, axis=axis)

    def fast_apply(self, f, names):
        # must return keys::list, values::list, mutated::bool
        try:
            starts, ends = lib.generate_slices(self.slabels, self.ngroups)
        except Exception:
            # fails when all -1
            return [], True

        sdata = self._get_sorted_data()
        results, mutated = reduction.apply_frame_axis0(sdata, f, names,
                                                       starts, ends)

        return results, mutated

    def _chop(self, sdata, slice_obj):
        if self.axis == 0:
            return sdata.iloc[slice_obj]
        else:
            return sdata._slice(slice_obj, axis=1)  # .loc[:, slice_obj]


class NDFrameSplitter(DataSplitter):

    def __init__(self, data, labels, ngroups, axis=0):
        super(NDFrameSplitter, self).__init__(data, labels, ngroups, axis=axis)

        self.factory = data._constructor

    def _get_sorted_data(self):
        # this is the BlockManager
        data = self.data._data

        # this is sort of wasteful but...
        sorted_axis = data.axes[self.axis].take(self.sort_idx)
        sorted_data = data.reindex_axis(sorted_axis, axis=self.axis)

        return sorted_data

    def _chop(self, sdata, slice_obj):
        return self.factory(sdata.get_slice(slice_obj, axis=self.axis))


def get_splitter(data, *args, **kwargs):
    if isinstance(data, Series):
        klass = SeriesSplitter
    elif isinstance(data, DataFrame):
        klass = FrameSplitter
    else:
        klass = NDFrameSplitter

    return klass(data, *args, **kwargs)
