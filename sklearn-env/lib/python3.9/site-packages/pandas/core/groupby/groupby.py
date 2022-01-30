"""
Provide the groupby split-apply-combine paradigm. Define the GroupBy
class providing the base-class of operations.

The SeriesGroupBy and DataFrameGroupBy sub-class
(defined in pandas.core.groupby.generic)
expose these user-facing objects to provide specific functionality.
"""
from __future__ import annotations

from contextlib import contextmanager
import datetime
from functools import (
    partial,
    wraps,
)
import inspect
from textwrap import dedent
import types
from typing import (
    Callable,
    Hashable,
    Iterable,
    Iterator,
    List,
    Literal,
    Mapping,
    Sequence,
    TypeVar,
    Union,
    cast,
    final,
)
import warnings

import numpy as np

from pandas._config.config import option_context

from pandas._libs import (
    Timestamp,
    lib,
)
import pandas._libs.groupby as libgroupby
from pandas._typing import (
    ArrayLike,
    IndexLabel,
    NDFrameT,
    PositionalIndexer,
    RandomState,
    Scalar,
    T,
    npt,
)
from pandas.compat.numpy import function as nv
from pandas.errors import AbstractMethodError
from pandas.util._decorators import (
    Appender,
    Substitution,
    cache_readonly,
    doc,
)
from pandas.util._exceptions import find_stack_level

from pandas.core.dtypes.common import (
    is_bool_dtype,
    is_datetime64_dtype,
    is_float_dtype,
    is_integer,
    is_integer_dtype,
    is_numeric_dtype,
    is_object_dtype,
    is_scalar,
    is_timedelta64_dtype,
)
from pandas.core.dtypes.missing import (
    isna,
    notna,
)

from pandas.core import nanops
from pandas.core._numba import executor
import pandas.core.algorithms as algorithms
from pandas.core.arrays import (
    BaseMaskedArray,
    BooleanArray,
    Categorical,
    ExtensionArray,
)
from pandas.core.base import (
    DataError,
    PandasObject,
    SelectionMixin,
)
import pandas.core.common as com
from pandas.core.frame import DataFrame
from pandas.core.generic import NDFrame
from pandas.core.groupby import (
    base,
    numba_,
    ops,
)
from pandas.core.groupby.indexing import (
    GroupByIndexingMixin,
    GroupByNthSelector,
)
from pandas.core.indexes.api import (
    CategoricalIndex,
    Index,
    MultiIndex,
)
from pandas.core.internals.blocks import ensure_block_shape
import pandas.core.sample as sample
from pandas.core.series import Series
from pandas.core.sorting import get_group_index_sorter
from pandas.core.util.numba_ import (
    NUMBA_FUNC_CACHE,
    maybe_use_numba,
)

_common_see_also = """
        See Also
        --------
        Series.%(name)s : Apply a function %(name)s to a Series.
        DataFrame.%(name)s : Apply a function %(name)s
            to each row or column of a DataFrame.
"""

_apply_docs = {
    "template": """
    Apply function ``func`` group-wise and combine the results together.

    The function passed to ``apply`` must take a {input} as its first
    argument and return a DataFrame, Series or scalar. ``apply`` will
    then take care of combining the results back together into a single
    dataframe or series. ``apply`` is therefore a highly flexible
    grouping method.

    While ``apply`` is a very flexible method, its downside is that
    using it can be quite a bit slower than using more specific methods
    like ``agg`` or ``transform``. Pandas offers a wide range of method that will
    be much faster than using ``apply`` for their specific purposes, so try to
    use them before reaching for ``apply``.

    Parameters
    ----------
    func : callable
        A callable that takes a {input} as its first argument, and
        returns a dataframe, a series or a scalar. In addition the
        callable may take positional and keyword arguments.
    args, kwargs : tuple and dict
        Optional positional and keyword arguments to pass to ``func``.

    Returns
    -------
    applied : Series or DataFrame

    See Also
    --------
    pipe : Apply function to the full GroupBy object instead of to each
        group.
    aggregate : Apply aggregate function to the GroupBy object.
    transform : Apply function column-by-column to the GroupBy object.
    Series.apply : Apply a function to a Series.
    DataFrame.apply : Apply a function to each row or column of a DataFrame.

    Notes
    -----

    .. versionchanged:: 1.3.0

        The resulting dtype will reflect the return value of the passed ``func``,
        see the examples below.

    Functions that mutate the passed object can produce unexpected
    behavior or errors and are not supported. See :ref:`gotchas.udf-mutation`
    for more details.

    Examples
    --------
    {examples}
    """,
    "dataframe_examples": """
    >>> df = pd.DataFrame({'A': 'a a b'.split(),
    ...                    'B': [1,2,3],
    ...                    'C': [4,6,5]})
    >>> g = df.groupby('A')

    Notice that ``g`` has two groups, ``a`` and ``b``.
    Calling `apply` in various ways, we can get different grouping results:

    Example 1: below the function passed to `apply` takes a DataFrame as
    its argument and returns a DataFrame. `apply` combines the result for
    each group together into a new DataFrame:

    >>> g[['B', 'C']].apply(lambda x: x / x.sum())
              B    C
    0  0.333333  0.4
    1  0.666667  0.6
    2  1.000000  1.0

    Example 2: The function passed to `apply` takes a DataFrame as
    its argument and returns a Series.  `apply` combines the result for
    each group together into a new DataFrame.

    .. versionchanged:: 1.3.0

        The resulting dtype will reflect the return value of the passed ``func``.

    >>> g[['B', 'C']].apply(lambda x: x.astype(float).max() - x.min())
         B    C
    A
    a  1.0  2.0
    b  0.0  0.0

    Example 3: The function passed to `apply` takes a DataFrame as
    its argument and returns a scalar. `apply` combines the result for
    each group together into a Series, including setting the index as
    appropriate:

    >>> g.apply(lambda x: x.C.max() - x.B.min())
    A
    a    5
    b    2
    dtype: int64""",
    "series_examples": """
    >>> s = pd.Series([0, 1, 2], index='a a b'.split())
    >>> g = s.groupby(s.index)

    From ``s`` above we can see that ``g`` has two groups, ``a`` and ``b``.
    Calling `apply` in various ways, we can get different grouping results:

    Example 1: The function passed to `apply` takes a Series as
    its argument and returns a Series.  `apply` combines the result for
    each group together into a new Series.

    .. versionchanged:: 1.3.0

        The resulting dtype will reflect the return value of the passed ``func``.

    >>> g.apply(lambda x: x*2 if x.name == 'a' else x/2)
    a    0.0
    a    2.0
    b    1.0
    dtype: float64

    Example 2: The function passed to `apply` takes a Series as
    its argument and returns a scalar. `apply` combines the result for
    each group together into a Series, including setting the index as
    appropriate:

    >>> g.apply(lambda x: x.max() - x.min())
    a    1
    b    0
    dtype: int64""",
}

_groupby_agg_method_template = """
Compute {fname} of group values.

Parameters
----------
numeric_only : bool, default {no}
    Include only float, int, boolean columns. If None, will attempt to use
    everything, then use only numeric data.
min_count : int, default {mc}
    The required number of valid values to perform the operation. If fewer
    than ``min_count`` non-NA values are present the result will be NA.

Returns
-------
Series or DataFrame
    Computed {fname} of values within each group.
"""

_pipe_template = """
Apply a function `func` with arguments to this %(klass)s object and return
the function's result.

Use `.pipe` when you want to improve readability by chaining together
functions that expect Series, DataFrames, GroupBy or Resampler objects.
Instead of writing

>>> h(g(f(df.groupby('group')), arg1=a), arg2=b, arg3=c)  # doctest: +SKIP

You can write

>>> (df.groupby('group')
...    .pipe(f)
...    .pipe(g, arg1=a)
...    .pipe(h, arg2=b, arg3=c))  # doctest: +SKIP

which is much more readable.

Parameters
----------
func : callable or tuple of (callable, str)
    Function to apply to this %(klass)s object or, alternatively,
    a `(callable, data_keyword)` tuple where `data_keyword` is a
    string indicating the keyword of `callable` that expects the
    %(klass)s object.
args : iterable, optional
       Positional arguments passed into `func`.
kwargs : dict, optional
         A dictionary of keyword arguments passed into `func`.

Returns
-------
object : the return type of `func`.

See Also
--------
Series.pipe : Apply a function with arguments to a series.
DataFrame.pipe: Apply a function with arguments to a dataframe.
apply : Apply function to each group instead of to the
    full %(klass)s object.

Notes
-----
See more `here
<https://pandas.pydata.org/pandas-docs/stable/user_guide/groupby.html#piping-function-calls>`_

Examples
--------
%(examples)s
"""

_transform_template = """
Call function producing a like-indexed %(klass)s on each group and
return a %(klass)s having the same indexes as the original object
filled with the transformed values.

Parameters
----------
f : function
    Function to apply to each group.

    Can also accept a Numba JIT function with
    ``engine='numba'`` specified.

    If the ``'numba'`` engine is chosen, the function must be
    a user defined function with ``values`` and ``index`` as the
    first and second arguments respectively in the function signature.
    Each group's index will be passed to the user defined function
    and optionally available for use.

    .. versionchanged:: 1.1.0
*args
    Positional arguments to pass to func.
engine : str, default None
    * ``'cython'`` : Runs the function through C-extensions from cython.
    * ``'numba'`` : Runs the function through JIT compiled code from numba.
    * ``None`` : Defaults to ``'cython'`` or the global setting ``compute.use_numba``

    .. versionadded:: 1.1.0
engine_kwargs : dict, default None
    * For ``'cython'`` engine, there are no accepted ``engine_kwargs``
    * For ``'numba'`` engine, the engine can accept ``nopython``, ``nogil``
      and ``parallel`` dictionary keys. The values must either be ``True`` or
      ``False``. The default ``engine_kwargs`` for the ``'numba'`` engine is
      ``{'nopython': True, 'nogil': False, 'parallel': False}`` and will be
      applied to the function

    .. versionadded:: 1.1.0
**kwargs
    Keyword arguments to be passed into func.

Returns
-------
%(klass)s

See Also
--------
%(klass)s.groupby.apply : Apply function ``func`` group-wise and combine
    the results together.
%(klass)s.groupby.aggregate : Aggregate using one or more
    operations over the specified axis.
%(klass)s.transform : Call ``func`` on self producing a %(klass)s with the
    same axis shape as self.

Notes
-----
Each group is endowed the attribute 'name' in case you need to know
which group you are working on.

The current implementation imposes three requirements on f:

* f must return a value that either has the same shape as the input
  subframe or can be broadcast to the shape of the input subframe.
  For example, if `f` returns a scalar it will be broadcast to have the
  same shape as the input subframe.
* if this is a DataFrame, f must support application column-by-column
  in the subframe. If f also supports application to the entire subframe,
  then a fast path is used starting from the second chunk.
* f must not mutate groups. Mutation is not supported and may
  produce unexpected results. See :ref:`gotchas.udf-mutation` for more details.

When using ``engine='numba'``, there will be no "fall back" behavior internally.
The group data and group index will be passed as numpy arrays to the JITed
user defined function, and no alternative execution attempts will be tried.

.. versionchanged:: 1.3.0

    The resulting dtype will reflect the return value of the passed ``func``,
    see the examples below.

Examples
--------

>>> df = pd.DataFrame({'A' : ['foo', 'bar', 'foo', 'bar',
...                           'foo', 'bar'],
...                    'B' : ['one', 'one', 'two', 'three',
...                           'two', 'two'],
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

Broadcast result of the transformation

>>> grouped.transform(lambda x: x.max() - x.min())
   C    D
0  4  6.0
1  3  8.0
2  4  6.0
3  3  8.0
4  4  6.0
5  3  8.0

.. versionchanged:: 1.3.0

    The resulting dtype will reflect the return value of the passed ``func``,
    for example:

>>> grouped[['C', 'D']].transform(lambda x: x.astype(int).max())
   C  D
0  5  8
1  5  9
2  5  8
3  5  9
4  5  8
5  5  9
"""

_agg_template = """
Aggregate using one or more operations over the specified axis.

Parameters
----------
func : function, str, list or dict
    Function to use for aggregating the data. If a function, must either
    work when passed a {klass} or when passed to {klass}.apply.

    Accepted combinations are:

    - function
    - string function name
    - list of functions and/or function names, e.g. ``[np.sum, 'mean']``
    - dict of axis labels -> functions, function names or list of such.

    Can also accept a Numba JIT function with
    ``engine='numba'`` specified. Only passing a single function is supported
    with this engine.

    If the ``'numba'`` engine is chosen, the function must be
    a user defined function with ``values`` and ``index`` as the
    first and second arguments respectively in the function signature.
    Each group's index will be passed to the user defined function
    and optionally available for use.

    .. versionchanged:: 1.1.0
*args
    Positional arguments to pass to func.
engine : str, default None
    * ``'cython'`` : Runs the function through C-extensions from cython.
    * ``'numba'`` : Runs the function through JIT compiled code from numba.
    * ``None`` : Defaults to ``'cython'`` or globally setting ``compute.use_numba``

    .. versionadded:: 1.1.0
engine_kwargs : dict, default None
    * For ``'cython'`` engine, there are no accepted ``engine_kwargs``
    * For ``'numba'`` engine, the engine can accept ``nopython``, ``nogil``
      and ``parallel`` dictionary keys. The values must either be ``True`` or
      ``False``. The default ``engine_kwargs`` for the ``'numba'`` engine is
      ``{{'nopython': True, 'nogil': False, 'parallel': False}}`` and will be
      applied to the function

    .. versionadded:: 1.1.0
**kwargs
    Keyword arguments to be passed into func.

Returns
-------
{klass}

See Also
--------
{klass}.groupby.apply : Apply function func group-wise
    and combine the results together.
{klass}.groupby.transform : Aggregate using one or more
    operations over the specified axis.
{klass}.aggregate : Transforms the Series on each group
    based on the given function.

Notes
-----
When using ``engine='numba'``, there will be no "fall back" behavior internally.
The group data and group index will be passed as numpy arrays to the JITed
user defined function, and no alternative execution attempts will be tried.

Functions that mutate the passed object can produce unexpected
behavior or errors and are not supported. See :ref:`gotchas.udf-mutation`
for more details.

.. versionchanged:: 1.3.0

    The resulting dtype will reflect the return value of the passed ``func``,
    see the examples below.
{examples}"""


@final
class GroupByPlot(PandasObject):
    """
    Class implementing the .plot attribute for groupby objects.
    """

    def __init__(self, groupby: GroupBy):
        self._groupby = groupby

    def __call__(self, *args, **kwargs):
        def f(self):
            return self.plot(*args, **kwargs)

        f.__name__ = "plot"
        return self._groupby.apply(f)

    def __getattr__(self, name: str):
        def attr(*args, **kwargs):
            def f(self):
                return getattr(self.plot, name)(*args, **kwargs)

            return self._groupby.apply(f)

        return attr


_KeysArgType = Union[
    Hashable,
    List[Hashable],
    Callable[[Hashable], Hashable],
    List[Callable[[Hashable], Hashable]],
    Mapping[Hashable, Hashable],
]


class BaseGroupBy(PandasObject, SelectionMixin[NDFrameT], GroupByIndexingMixin):
    _group_selection: IndexLabel | None = None
    _apply_allowlist: frozenset[str] = frozenset()
    _hidden_attrs = PandasObject._hidden_attrs | {
        "as_index",
        "axis",
        "dropna",
        "exclusions",
        "grouper",
        "group_keys",
        "keys",
        "level",
        "mutated",
        "obj",
        "observed",
        "sort",
        "squeeze",
    }

    axis: int
    grouper: ops.BaseGrouper
    group_keys: bool

    @final
    def __len__(self) -> int:
        return len(self.groups)

    @final
    def __repr__(self) -> str:
        # TODO: Better repr for GroupBy object
        return object.__repr__(self)

    @final
    @property
    def groups(self) -> dict[Hashable, np.ndarray]:
        """
        Dict {group name -> group labels}.
        """
        return self.grouper.groups

    @final
    @property
    def ngroups(self) -> int:
        return self.grouper.ngroups

    @final
    @property
    def indices(self):
        """
        Dict {group name -> group indices}.
        """
        return self.grouper.indices

    @final
    def _get_indices(self, names):
        """
        Safe get multiple indices, translate keys for
        datelike to underlying repr.
        """

        def get_converter(s):
            # possibly convert to the actual key types
            # in the indices, could be a Timestamp or a np.datetime64
            if isinstance(s, datetime.datetime):
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
            index_sample = None  # Dummy sample

        name_sample = names[0]
        if isinstance(index_sample, tuple):
            if not isinstance(name_sample, tuple):
                msg = "must supply a tuple to get_group with multiple grouping keys"
                raise ValueError(msg)
            if not len(name_sample) == len(index_sample):
                try:
                    # If the original grouper was a tuple
                    return [self.indices[name] for name in names]
                except KeyError as err:
                    # turns out it wasn't a tuple
                    msg = (
                        "must supply a same-length tuple to get_group "
                        "with multiple grouping keys"
                    )
                    raise ValueError(msg) from err

            converters = [get_converter(s) for s in index_sample]
            names = (tuple(f(n) for f, n in zip(converters, name)) for name in names)

        else:
            converter = get_converter(index_sample)
            names = (converter(name) for name in names)

        return [self.indices.get(name, []) for name in names]

    @final
    def _get_index(self, name):
        """
        Safe get index, translate keys for datelike to underlying repr.
        """
        return self._get_indices([name])[0]

    @final
    @cache_readonly
    def _selected_obj(self):
        # Note: _selected_obj is always just `self.obj` for SeriesGroupBy

        if self._selection is None or isinstance(self.obj, Series):
            if self._group_selection is not None:
                return self.obj[self._group_selection]
            return self.obj
        else:
            return self.obj[self._selection]

    @final
    def _dir_additions(self) -> set[str]:
        return self.obj._dir_additions() | self._apply_allowlist

    @Substitution(
        klass="GroupBy",
        examples=dedent(
            """\
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
        b  2"""
        ),
    )
    @Appender(_pipe_template)
    def pipe(
        self,
        func: Callable[..., T] | tuple[Callable[..., T], str],
        *args,
        **kwargs,
    ) -> T:
        return com.pipe(self, func, *args, **kwargs)

    plot = property(GroupByPlot)

    @final
    def get_group(self, name, obj=None) -> DataFrame | Series:
        """
        Construct DataFrame from group with provided name.

        Parameters
        ----------
        name : object
            The name of the group to get as a DataFrame.
        obj : DataFrame, default None
            The DataFrame to take the DataFrame out of.  If
            it is None, the object groupby was called on will
            be used.

        Returns
        -------
        group : same type as obj
        """
        if obj is None:
            obj = self._selected_obj

        inds = self._get_index(name)
        if not len(inds):
            raise KeyError(name)

        return obj._take_with_is_copy(inds, axis=self.axis)

    @final
    def __iter__(self) -> Iterator[tuple[Hashable, NDFrameT]]:
        """
        Groupby iterator.

        Returns
        -------
        Generator yielding sequence of (name, subsetted object)
        for each group
        """
        return self.grouper.get_iterator(self._selected_obj, axis=self.axis)


# To track operations that expand dimensions, like ohlc
OutputFrameOrSeries = TypeVar("OutputFrameOrSeries", bound=NDFrame)


class GroupBy(BaseGroupBy[NDFrameT]):
    """
    Class for grouping and aggregating relational data.

    See aggregate, transform, and apply functions on this object.

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
    name : str
        Most users should ignore this

    Returns
    -------
    **Attributes**
    groups : dict
        {group name -> group labels}
    len(grouped) : int
        Number of groups

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
    """

    grouper: ops.BaseGrouper
    as_index: bool

    @final
    def __init__(
        self,
        obj: NDFrameT,
        keys: _KeysArgType | None = None,
        axis: int = 0,
        level: IndexLabel | None = None,
        grouper: ops.BaseGrouper | None = None,
        exclusions: frozenset[Hashable] | None = None,
        selection: IndexLabel | None = None,
        as_index: bool = True,
        sort: bool = True,
        group_keys: bool = True,
        squeeze: bool = False,
        observed: bool = False,
        mutated: bool = False,
        dropna: bool = True,
    ):

        self._selection = selection

        assert isinstance(obj, NDFrame), type(obj)

        self.level = level

        if not as_index:
            if not isinstance(obj, DataFrame):
                raise TypeError("as_index=False only valid with DataFrame")
            if axis != 0:
                raise ValueError("as_index=False only valid for axis=0")

        self.as_index = as_index
        self.keys = keys
        self.sort = sort
        self.group_keys = group_keys
        self.squeeze = squeeze
        self.observed = observed
        self.mutated = mutated
        self.dropna = dropna

        if grouper is None:
            from pandas.core.groupby.grouper import get_grouper

            grouper, exclusions, obj = get_grouper(
                obj,
                keys,
                axis=axis,
                level=level,
                sort=sort,
                observed=observed,
                mutated=self.mutated,
                dropna=self.dropna,
            )

        self.obj = obj
        self.axis = obj._get_axis_number(axis)
        self.grouper = grouper
        self.exclusions = frozenset(exclusions) if exclusions else frozenset()

    def __getattr__(self, attr: str):
        if attr in self._internal_names_set:
            return object.__getattribute__(self, attr)
        if attr in self.obj:
            return self[attr]

        raise AttributeError(
            f"'{type(self).__name__}' object has no attribute '{attr}'"
        )

    def __getattribute__(self, attr: str):
        # Intercept nth to allow both call and index
        if attr == "nth":
            return GroupByNthSelector(self)
        elif attr == "nth_actual":
            return super().__getattribute__("nth")
        else:
            return super().__getattribute__(attr)

    @final
    def _make_wrapper(self, name: str) -> Callable:
        assert name in self._apply_allowlist

        with self._group_selection_context():
            # need to setup the selection
            # as are not passed directly but in the grouper
            f = getattr(self._obj_with_exclusions, name)
            if not isinstance(f, types.MethodType):
                return self.apply(lambda self: getattr(self, name))

        f = getattr(type(self._obj_with_exclusions), name)
        sig = inspect.signature(f)

        def wrapper(*args, **kwargs):
            # a little trickery for aggregation functions that need an axis
            # argument
            if "axis" in sig.parameters:
                if kwargs.get("axis", None) is None:
                    kwargs["axis"] = self.axis

            def curried(x):
                return f(x, *args, **kwargs)

            # preserve the name so we can detect it when calling plot methods,
            # to avoid duplicates
            curried.__name__ = name

            # special case otherwise extra plots are created when catching the
            # exception below
            if name in base.plotting_methods:
                return self.apply(curried)

            return self._python_apply_general(curried, self._obj_with_exclusions)

        wrapper.__name__ = name
        return wrapper

    # -----------------------------------------------------------------
    # Selection

    @final
    def _set_group_selection(self) -> None:
        """
        Create group based selection.

        Used when selection is not passed directly but instead via a grouper.

        NOTE: this should be paired with a call to _reset_group_selection
        """
        # This is a no-op for SeriesGroupBy
        grp = self.grouper
        if not (
            self.as_index
            and grp.groupings is not None
            and self.obj.ndim > 1
            and self._group_selection is None
        ):
            return

        groupers = [g.name for g in grp.groupings if g.level is None and g.in_axis]

        if len(groupers):
            # GH12839 clear selected obj cache when group selection changes
            ax = self.obj._info_axis
            self._group_selection = ax.difference(Index(groupers), sort=False).tolist()
            self._reset_cache("_selected_obj")

    @final
    def _reset_group_selection(self) -> None:
        """
        Clear group based selection.

        Used for methods needing to return info on each group regardless of
        whether a group selection was previously set.
        """
        if self._group_selection is not None:
            # GH12839 clear cached selection too when changing group selection
            self._group_selection = None
            self._reset_cache("_selected_obj")

    @contextmanager
    def _group_selection_context(self) -> Iterator[GroupBy]:
        """
        Set / reset the _group_selection_context.
        """
        self._set_group_selection()
        try:
            yield self
        finally:
            self._reset_group_selection()

    def _iterate_slices(self) -> Iterable[Series]:
        raise AbstractMethodError(self)

    # -----------------------------------------------------------------
    # Dispatch/Wrapping

    @final
    def _concat_objects(self, values, not_indexed_same: bool = False):
        from pandas.core.reshape.concat import concat

        def reset_identity(values):
            # reset the identities of the components
            # of the values to prevent aliasing
            for v in com.not_none(*values):
                ax = v._get_axis(self.axis)
                ax._reset_identity()
            return values

        if not not_indexed_same:
            result = concat(values, axis=self.axis)

            ax = self._selected_obj._get_axis(self.axis)
            if self.dropna:
                labels = self.grouper.group_info[0]
                mask = labels != -1
                ax = ax[mask]

            # this is a very unfortunate situation
            # we can't use reindex to restore the original order
            # when the ax has duplicates
            # so we resort to this
            # GH 14776, 30667
            if ax.has_duplicates and not result.axes[self.axis].equals(ax):
                indexer, _ = result.index.get_indexer_non_unique(ax._values)
                indexer = algorithms.unique1d(indexer)
                result = result.take(indexer, axis=self.axis)
            else:
                result = result.reindex(ax, axis=self.axis, copy=False)

        elif self.group_keys:

            values = reset_identity(values)
            if self.as_index:

                # possible MI return case
                group_keys = self.grouper.result_index
                group_levels = self.grouper.levels
                group_names = self.grouper.names

                result = concat(
                    values,
                    axis=self.axis,
                    keys=group_keys,
                    levels=group_levels,
                    names=group_names,
                    sort=False,
                )
            else:

                # GH5610, returns a MI, with the first level being a
                # range index
                keys = list(range(len(values)))
                result = concat(values, axis=self.axis, keys=keys)
        else:
            values = reset_identity(values)
            result = concat(values, axis=self.axis)

        name = self.obj.name if self.obj.ndim == 1 else self._selection
        if isinstance(result, Series) and name is not None:

            result.name = name

        return result

    @final
    def _set_result_index_ordered(
        self, result: OutputFrameOrSeries
    ) -> OutputFrameOrSeries:
        # set the result index on the passed values object and
        # return the new object, xref 8046

        if self.grouper.is_monotonic:
            # shortcut if we have an already ordered grouper
            result.set_axis(self.obj._get_axis(self.axis), axis=self.axis, inplace=True)
            return result

        # row order is scrambled => sort the rows by position in original index
        original_positions = Index(
            np.concatenate(self._get_indices(self.grouper.result_index))
        )
        result.set_axis(original_positions, axis=self.axis, inplace=True)
        result = result.sort_index(axis=self.axis)

        dropped_rows = len(result.index) < len(self.obj.index)

        if dropped_rows:
            # get index by slicing original index according to original positions
            # slice drops attrs => use set_axis when no rows were dropped
            sorted_indexer = result.index
            result.index = self._selected_obj.index[sorted_indexer]
        else:
            result.set_axis(self.obj._get_axis(self.axis), axis=self.axis, inplace=True)

        return result

    def _indexed_output_to_ndframe(
        self, result: Mapping[base.OutputKey, ArrayLike]
    ) -> Series | DataFrame:
        raise AbstractMethodError(self)

    @final
    def _wrap_aggregated_output(
        self,
        output: Series | DataFrame | Mapping[base.OutputKey, ArrayLike],
        qs: npt.NDArray[np.float64] | None = None,
    ):
        """
        Wraps the output of GroupBy aggregations into the expected result.

        Parameters
        ----------
        output : Series, DataFrame, or Mapping[base.OutputKey, ArrayLike]
           Data to wrap.

        Returns
        -------
        Series or DataFrame
        """

        if isinstance(output, (Series, DataFrame)):
            # We get here (for DataFrameGroupBy) if we used Manager.grouped_reduce,
            #  in which case our columns are already set correctly.
            # ATM we do not get here for SeriesGroupBy; when we do, we will
            #  need to require that result.name already match self.obj.name
            result = output
        else:
            result = self._indexed_output_to_ndframe(output)

        if not self.as_index:
            # `not self.as_index` is only relevant for DataFrameGroupBy,
            #   enforced in __init__
            self._insert_inaxis_grouper_inplace(result)
            result = result._consolidate()
            index = Index(range(self.grouper.ngroups))

        else:
            index = self.grouper.result_index

        if qs is not None:
            # We get here with len(qs) != 1 and not self.as_index
            #  in test_pass_args_kwargs
            index = _insert_quantile_level(index, qs)

        result.index = index

        if self.axis == 1:
            # Only relevant for DataFrameGroupBy, no-op for SeriesGroupBy
            result = result.T
            if result.index.equals(self.obj.index):
                # Retain e.g. DatetimeIndex/TimedeltaIndex freq
                result.index = self.obj.index.copy()
                # TODO: Do this more systematically

        return self._reindex_output(result, qs=qs)

    @final
    def _wrap_transformed_output(
        self, output: Mapping[base.OutputKey, ArrayLike]
    ) -> Series | DataFrame:
        """
        Wraps the output of GroupBy transformations into the expected result.

        Parameters
        ----------
        output : Mapping[base.OutputKey, ArrayLike]
            Data to wrap.

        Returns
        -------
        Series or DataFrame
            Series for SeriesGroupBy, DataFrame for DataFrameGroupBy
        """
        if isinstance(output, (Series, DataFrame)):
            result = output
        else:
            result = self._indexed_output_to_ndframe(output)

        if self.axis == 1:
            # Only relevant for DataFrameGroupBy
            result = result.T
            result.columns = self.obj.columns

        result.index = self.obj.index
        return result

    def _wrap_applied_output(self, data, values: list, not_indexed_same: bool = False):
        raise AbstractMethodError(self)

    def _resolve_numeric_only(self, numeric_only: bool | lib.NoDefault) -> bool:
        """
        Determine subclass-specific default value for 'numeric_only'.

        For SeriesGroupBy we want the default to be False (to match Series behavior).
        For DataFrameGroupBy we want it to be True (for backwards-compat).

        Parameters
        ----------
        numeric_only : bool or lib.no_default

        Returns
        -------
        bool
        """
        # GH#41291
        if numeric_only is lib.no_default:
            # i.e. not explicitly passed by user
            if self.obj.ndim == 2:
                # i.e. DataFrameGroupBy
                numeric_only = True
                # GH#42395 GH#43108 GH#43154
                # Regression from 1.2.5 to 1.3 caused object columns to be dropped
                if self.axis:
                    obj = self._obj_with_exclusions.T
                else:
                    obj = self._obj_with_exclusions
                check = obj._get_numeric_data()
                if len(obj.columns) and not len(check.columns) and not obj.empty:
                    numeric_only = False
                    # TODO: v1.4+ Add FutureWarning

            else:
                numeric_only = False

        # error: Incompatible return value type (got "Union[bool, NoDefault]",
        # expected "bool")
        return numeric_only  # type: ignore[return-value]

    # -----------------------------------------------------------------
    # numba

    @final
    def _numba_prep(self, func, data):
        if not callable(func):
            raise NotImplementedError(
                "Numba engine can only be used with a single function."
            )
        ids, _, ngroups = self.grouper.group_info
        sorted_index = get_group_index_sorter(ids, ngroups)
        sorted_ids = algorithms.take_nd(ids, sorted_index, allow_fill=False)

        sorted_data = data.take(sorted_index, axis=self.axis).to_numpy()
        sorted_index_data = data.index.take(sorted_index).to_numpy()

        starts, ends = lib.generate_slices(sorted_ids, ngroups)
        return (
            starts,
            ends,
            sorted_index_data,
            sorted_data,
        )

    def _numba_agg_general(
        self,
        func: Callable,
        engine_kwargs: dict[str, bool] | None,
        numba_cache_key_str: str,
        *aggregator_args,
    ):
        """
        Perform groupby with a standard numerical aggregation function (e.g. mean)
        with Numba.
        """
        if not self.as_index:
            raise NotImplementedError(
                "as_index=False is not supported. Use .reset_index() instead."
            )
        if self.axis == 1:
            raise NotImplementedError("axis=1 is not supported.")

        with self._group_selection_context():
            data = self._selected_obj
        df = data if data.ndim == 2 else data.to_frame()
        starts, ends, sorted_index, sorted_data = self._numba_prep(func, df)
        aggregator = executor.generate_shared_aggregator(
            func, engine_kwargs, numba_cache_key_str
        )
        result = aggregator(sorted_data, starts, ends, 0, *aggregator_args)

        cache_key = (func, numba_cache_key_str)
        if cache_key not in NUMBA_FUNC_CACHE:
            NUMBA_FUNC_CACHE[cache_key] = aggregator

        index = self.grouper.result_index
        if data.ndim == 1:
            result_kwargs = {"name": data.name}
            result = result.ravel()
        else:
            result_kwargs = {"columns": data.columns}
        return data._constructor(result, index=index, **result_kwargs)

    @final
    def _transform_with_numba(self, data, func, *args, engine_kwargs=None, **kwargs):
        """
        Perform groupby transform routine with the numba engine.

        This routine mimics the data splitting routine of the DataSplitter class
        to generate the indices of each group in the sorted data and then passes the
        data and indices into a Numba jitted function.
        """
        starts, ends, sorted_index, sorted_data = self._numba_prep(func, data)

        numba_transform_func = numba_.generate_numba_transform_func(
            kwargs, func, engine_kwargs
        )
        result = numba_transform_func(
            sorted_data,
            sorted_index,
            starts,
            ends,
            len(data.columns),
            *args,
        )

        cache_key = (func, "groupby_transform")
        if cache_key not in NUMBA_FUNC_CACHE:
            NUMBA_FUNC_CACHE[cache_key] = numba_transform_func

        # result values needs to be resorted to their original positions since we
        # evaluated the data sorted by group
        return result.take(np.argsort(sorted_index), axis=0)

    @final
    def _aggregate_with_numba(self, data, func, *args, engine_kwargs=None, **kwargs):
        """
        Perform groupby aggregation routine with the numba engine.

        This routine mimics the data splitting routine of the DataSplitter class
        to generate the indices of each group in the sorted data and then passes the
        data and indices into a Numba jitted function.
        """
        starts, ends, sorted_index, sorted_data = self._numba_prep(func, data)

        numba_agg_func = numba_.generate_numba_agg_func(kwargs, func, engine_kwargs)
        result = numba_agg_func(
            sorted_data,
            sorted_index,
            starts,
            ends,
            len(data.columns),
            *args,
        )

        cache_key = (func, "groupby_agg")
        if cache_key not in NUMBA_FUNC_CACHE:
            NUMBA_FUNC_CACHE[cache_key] = numba_agg_func

        return result

    # -----------------------------------------------------------------
    # apply/agg/transform

    @Appender(
        _apply_docs["template"].format(
            input="dataframe", examples=_apply_docs["dataframe_examples"]
        )
    )
    def apply(self, func, *args, **kwargs):

        func = com.is_builtin_func(func)

        # this is needed so we don't try and wrap strings. If we could
        # resolve functions to their callable functions prior, this
        # wouldn't be needed
        if args or kwargs:
            if callable(func):

                @wraps(func)
                def f(g):
                    with np.errstate(all="ignore"):
                        return func(g, *args, **kwargs)

            elif hasattr(nanops, "nan" + func):
                # TODO: should we wrap this in to e.g. _is_builtin_func?
                f = getattr(nanops, "nan" + func)

            else:
                raise ValueError(
                    "func must be a callable if args or kwargs are supplied"
                )
        elif isinstance(func, str):
            if hasattr(self, func):
                res = getattr(self, func)
                if callable(res):
                    return res()
                return res

            else:
                raise TypeError(f"apply func should be callable, not '{func}'")
        else:

            f = func

        # ignore SettingWithCopy here in case the user mutates
        with option_context("mode.chained_assignment", None):
            try:
                result = self._python_apply_general(f, self._selected_obj)
            except TypeError:
                # gh-20949
                # try again, with .apply acting as a filtering
                # operation, by excluding the grouping column
                # This would normally not be triggered
                # except if the udf is trying an operation that
                # fails on *some* columns, e.g. a numeric operation
                # on a string grouper column

                with self._group_selection_context():
                    return self._python_apply_general(f, self._selected_obj)

        return result

    @final
    def _python_apply_general(
        self,
        f: Callable,
        data: DataFrame | Series,
        not_indexed_same: bool | None = None,
    ) -> DataFrame | Series:
        """
        Apply function f in python space

        Parameters
        ----------
        f : callable
            Function to apply
        data : Series or DataFrame
            Data to apply f to
        not_indexed_same: bool, optional
            When specified, overrides the value of not_indexed_same. Apply behaves
            differently when the result index is equal to the input index, but
            this can be coincidental leading to value-dependent behavior.

        Returns
        -------
        Series or DataFrame
            data after applying f
        """
        values, mutated = self.grouper.apply(f, data, self.axis)

        if not_indexed_same is None:
            not_indexed_same = mutated or self.mutated

        return self._wrap_applied_output(
            data, values, not_indexed_same=not_indexed_same
        )

    @final
    def _python_agg_general(self, func, *args, **kwargs):
        func = com.is_builtin_func(func)
        f = lambda x: func(x, *args, **kwargs)

        # iterate through "columns" ex exclusions to populate output dict
        output: dict[base.OutputKey, ArrayLike] = {}

        if self.ngroups == 0:
            # agg_series below assumes ngroups > 0
            return self._python_apply_general(f, self._selected_obj)

        for idx, obj in enumerate(self._iterate_slices()):
            name = obj.name

            try:
                # if this function is invalid for this dtype, we will ignore it.
                result = self.grouper.agg_series(obj, f)
            except TypeError:
                warn_dropping_nuisance_columns_deprecated(type(self), "agg")
                continue

            key = base.OutputKey(label=name, position=idx)
            output[key] = result

        if not output:
            return self._python_apply_general(f, self._selected_obj)

        return self._wrap_aggregated_output(output)

    @final
    def _agg_general(
        self,
        numeric_only: bool = True,
        min_count: int = -1,
        *,
        alias: str,
        npfunc: Callable,
    ):

        with self._group_selection_context():
            # try a cython aggregation if we can
            result = self._cython_agg_general(
                how=alias,
                alt=npfunc,
                numeric_only=numeric_only,
                min_count=min_count,
            )
            return result.__finalize__(self.obj, method="groupby")

    def _agg_py_fallback(
        self, values: ArrayLike, ndim: int, alt: Callable
    ) -> ArrayLike:
        """
        Fallback to pure-python aggregation if _cython_operation raises
        NotImplementedError.
        """
        # We get here with a) EADtypes and b) object dtype

        if values.ndim == 1:
            # For DataFrameGroupBy we only get here with ExtensionArray
            ser = Series(values)
        else:
            # We only get here with values.dtype == object
            # TODO: special case not needed with ArrayManager
            df = DataFrame(values.T)
            # bc we split object blocks in grouped_reduce, we have only 1 col
            # otherwise we'd have to worry about block-splitting GH#39329
            assert df.shape[1] == 1
            # Avoid call to self.values that can occur in DataFrame
            #  reductions; see GH#28949
            ser = df.iloc[:, 0]

        # We do not get here with UDFs, so we know that our dtype
        #  should always be preserved by the implemented aggregations
        # TODO: Is this exactly right; see WrappedCythonOp get_result_dtype?
        res_values = self.grouper.agg_series(ser, alt, preserve_dtype=True)

        if isinstance(values, Categorical):
            # Because we only get here with known dtype-preserving
            #  reductions, we cast back to Categorical.
            # TODO: if we ever get "rank" working, exclude it here.
            res_values = type(values)._from_sequence(res_values, dtype=values.dtype)

        # If we are DataFrameGroupBy and went through a SeriesGroupByPath
        # then we need to reshape
        # GH#32223 includes case with IntegerArray values, ndarray res_values
        # test_groupby_duplicate_columns with object dtype values
        return ensure_block_shape(res_values, ndim=ndim)

    @final
    def _cython_agg_general(
        self, how: str, alt: Callable, numeric_only: bool, min_count: int = -1
    ):
        # Note: we never get here with how="ohlc" for DataFrameGroupBy;
        #  that goes through SeriesGroupBy

        data = self._get_data_to_aggregate()
        is_ser = data.ndim == 1

        if numeric_only:
            if is_ser and not is_numeric_dtype(self._selected_obj.dtype):
                # GH#41291 match Series behavior
                kwd_name = "numeric_only"
                if how in ["any", "all"]:
                    kwd_name = "bool_only"
                raise NotImplementedError(
                    f"{type(self).__name__}.{how} does not implement {kwd_name}."
                )
            elif not is_ser:
                data = data.get_numeric_data(copy=False)

        def array_func(values: ArrayLike) -> ArrayLike:
            try:
                result = self.grouper._cython_operation(
                    "aggregate", values, how, axis=data.ndim - 1, min_count=min_count
                )
            except NotImplementedError:
                # generally if we have numeric_only=False
                # and non-applicable functions
                # try to python agg
                # TODO: shouldn't min_count matter?
                result = self._agg_py_fallback(values, ndim=data.ndim, alt=alt)

            return result

        # TypeError -> we may have an exception in trying to aggregate
        #  continue and exclude the block
        new_mgr = data.grouped_reduce(array_func, ignore_failures=True)

        if not is_ser and len(new_mgr) < len(data):
            warn_dropping_nuisance_columns_deprecated(type(self), how)

        res = self._wrap_agged_manager(new_mgr)
        if is_ser:
            res.index = self.grouper.result_index
            return self._reindex_output(res)
        else:
            return res

    def _cython_transform(
        self, how: str, numeric_only: bool = True, axis: int = 0, **kwargs
    ):
        raise AbstractMethodError(self)

    @final
    def _transform(self, func, *args, engine=None, engine_kwargs=None, **kwargs):

        if maybe_use_numba(engine):
            # TODO: tests with self._selected_obj.ndim == 1 on DataFrameGroupBy
            with self._group_selection_context():
                data = self._selected_obj
            df = data if data.ndim == 2 else data.to_frame()
            result = self._transform_with_numba(
                df, func, *args, engine_kwargs=engine_kwargs, **kwargs
            )
            if self.obj.ndim == 2:
                return cast(DataFrame, self.obj)._constructor(
                    result, index=data.index, columns=data.columns
                )
            else:
                return cast(Series, self.obj)._constructor(
                    result.ravel(), index=data.index, name=data.name
                )

        # optimized transforms
        func = com.get_cython_func(func) or func

        if not isinstance(func, str):
            return self._transform_general(func, *args, **kwargs)

        elif func not in base.transform_kernel_allowlist:
            msg = f"'{func}' is not a valid function name for transform(name)"
            raise ValueError(msg)
        elif func in base.cythonized_kernels or func in base.transformation_kernels:
            # cythonized transform or canned "agg+broadcast"
            return getattr(self, func)(*args, **kwargs)

        else:
            # i.e. func in base.reduction_kernels

            # GH#30918 Use _transform_fast only when we know func is an aggregation
            # If func is a reduction, we need to broadcast the
            # result to the whole group. Compute func result
            # and deal with possible broadcasting below.
            # Temporarily set observed for dealing with categoricals.
            with com.temp_setattr(self, "observed", True):
                result = getattr(self, func)(*args, **kwargs)

            if self._can_use_transform_fast(result):
                return self._wrap_transform_fast_result(result)

            # only reached for DataFrameGroupBy
            return self._transform_general(func, *args, **kwargs)

    @final
    def _wrap_transform_fast_result(self, result: NDFrameT) -> NDFrameT:
        """
        Fast transform path for aggregations.
        """
        obj = self._obj_with_exclusions

        # for each col, reshape to size of original frame by take operation
        ids, _, _ = self.grouper.group_info
        result = result.reindex(self.grouper.result_index, copy=False)

        if self.obj.ndim == 1:
            # i.e. SeriesGroupBy
            out = algorithms.take_nd(result._values, ids)
            output = obj._constructor(out, index=obj.index, name=obj.name)
        else:
            output = result.take(ids, axis=0)
            output.index = obj.index
        return output

    # -----------------------------------------------------------------
    # Utilities

    @final
    def _apply_filter(self, indices, dropna):
        if len(indices) == 0:
            indices = np.array([], dtype="int64")
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

    @final
    def _cumcount_array(self, ascending: bool = True) -> np.ndarray:
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

    # -----------------------------------------------------------------

    @final
    @property
    def _obj_1d_constructor(self) -> type[Series]:
        # GH28330 preserve subclassed Series/DataFrames
        if isinstance(self.obj, DataFrame):
            return self.obj._constructor_sliced
        assert isinstance(self.obj, Series)
        return self.obj._constructor

    @final
    def _bool_agg(self, val_test: Literal["any", "all"], skipna: bool):
        """
        Shared func to call any / all Cython GroupBy implementations.
        """

        def objs_to_bool(vals: ArrayLike) -> tuple[np.ndarray, type]:
            if is_object_dtype(vals.dtype):
                # GH#37501: don't raise on pd.NA when skipna=True
                if skipna:
                    func = np.vectorize(
                        lambda x: bool(x) if not isna(x) else True, otypes=[bool]
                    )
                    vals = func(vals)
                else:
                    vals = vals.astype(bool, copy=False)

                vals = cast(np.ndarray, vals)
            elif isinstance(vals, BaseMaskedArray):
                vals = vals._data.astype(bool, copy=False)
            else:
                vals = vals.astype(bool, copy=False)

            return vals.view(np.int8), bool

        def result_to_bool(
            result: np.ndarray,
            inference: type,
            nullable: bool = False,
        ) -> ArrayLike:
            if nullable:
                return BooleanArray(result.astype(bool, copy=False), result == -1)
            else:
                return result.astype(inference, copy=False)

        return self._get_cythonized_result(
            libgroupby.group_any_all,
            numeric_only=False,
            cython_dtype=np.dtype(np.int8),
            needs_mask=True,
            needs_nullable=True,
            pre_processing=objs_to_bool,
            post_processing=result_to_bool,
            val_test=val_test,
            skipna=skipna,
        )

    @final
    @Substitution(name="groupby")
    @Appender(_common_see_also)
    def any(self, skipna: bool = True):
        """
        Return True if any value in the group is truthful, else False.

        Parameters
        ----------
        skipna : bool, default True
            Flag to ignore nan values during truth testing.

        Returns
        -------
        Series or DataFrame
            DataFrame or Series of boolean values, where a value is True if any element
            is True within its respective group, False otherwise.
        """
        return self._bool_agg("any", skipna)

    @final
    @Substitution(name="groupby")
    @Appender(_common_see_also)
    def all(self, skipna: bool = True):
        """
        Return True if all values in the group are truthful, else False.

        Parameters
        ----------
        skipna : bool, default True
            Flag to ignore nan values during truth testing.

        Returns
        -------
        Series or DataFrame
            DataFrame or Series of boolean values, where a value is True if all elements
            are True within its respective group, False otherwise.
        """
        return self._bool_agg("all", skipna)

    @final
    @Substitution(name="groupby")
    @Appender(_common_see_also)
    def count(self) -> Series | DataFrame:
        """
        Compute count of group, excluding missing values.

        Returns
        -------
        Series or DataFrame
            Count of values within each group.
        """
        data = self._get_data_to_aggregate()
        ids, _, ngroups = self.grouper.group_info
        mask = ids != -1

        is_series = data.ndim == 1

        def hfunc(bvalues: ArrayLike) -> ArrayLike:
            # TODO(EA2D): reshape would not be necessary with 2D EAs
            if bvalues.ndim == 1:
                # EA
                masked = mask & ~isna(bvalues).reshape(1, -1)
            else:
                masked = mask & ~isna(bvalues)

            counted = lib.count_level_2d(masked, labels=ids, max_bin=ngroups, axis=1)
            if is_series:
                assert counted.ndim == 2
                assert counted.shape[0] == 1
                return counted[0]
            return counted

        new_mgr = data.grouped_reduce(hfunc)

        # If we are grouping on categoricals we want unobserved categories to
        # return zero, rather than the default of NaN which the reindexing in
        # _wrap_agged_manager() returns. GH 35028
        with com.temp_setattr(self, "observed", True):
            result = self._wrap_agged_manager(new_mgr)

        if result.ndim == 1:
            result.index = self.grouper.result_index

        return self._reindex_output(result, fill_value=0)

    @final
    @Substitution(name="groupby")
    @Substitution(see_also=_common_see_also)
    def mean(
        self,
        numeric_only: bool | lib.NoDefault = lib.no_default,
        engine: str = "cython",
        engine_kwargs: dict[str, bool] | None = None,
    ):
        """
        Compute mean of groups, excluding missing values.

        Parameters
        ----------
        numeric_only : bool, default True
            Include only float, int, boolean columns. If None, will attempt to use
            everything, then use only numeric data.

        engine : str, default None
            * ``'cython'`` : Runs the operation through C-extensions from cython.
            * ``'numba'`` : Runs the operation through JIT compiled code from numba.
            * ``None`` : Defaults to ``'cython'`` or globally setting
              ``compute.use_numba``

            .. versionadded:: 1.4.0

        engine_kwargs : dict, default None
            * For ``'cython'`` engine, there are no accepted ``engine_kwargs``
            * For ``'numba'`` engine, the engine can accept ``nopython``, ``nogil``
              and ``parallel`` dictionary keys. The values must either be ``True`` or
              ``False``. The default ``engine_kwargs`` for the ``'numba'`` engine is
              ``{{'nopython': True, 'nogil': False, 'parallel': False}}``

            .. versionadded:: 1.4.0

        Returns
        -------
        pandas.Series or pandas.DataFrame
        %(see_also)s
        Examples
        --------
        >>> df = pd.DataFrame({'A': [1, 1, 2, 1, 2],
        ...                    'B': [np.nan, 2, 3, 4, 5],
        ...                    'C': [1, 2, 1, 1, 2]}, columns=['A', 'B', 'C'])

        Groupby one column and return the mean of the remaining columns in
        each group.

        >>> df.groupby('A').mean()
             B         C
        A
        1  3.0  1.333333
        2  4.0  1.500000

        Groupby two columns and return the mean of the remaining column.

        >>> df.groupby(['A', 'B']).mean()
                 C
        A B
        1 2.0  2.0
          4.0  1.0
        2 3.0  1.0
          5.0  2.0

        Groupby one column and return the mean of only particular column in
        the group.

        >>> df.groupby('A')['B'].mean()
        A
        1    3.0
        2    4.0
        Name: B, dtype: float64
        """
        numeric_only_bool = self._resolve_numeric_only(numeric_only)

        if maybe_use_numba(engine):
            from pandas.core._numba.kernels import sliding_mean

            return self._numba_agg_general(sliding_mean, engine_kwargs, "groupby_mean")
        else:
            result = self._cython_agg_general(
                "mean",
                alt=lambda x: Series(x).mean(numeric_only=numeric_only_bool),
                numeric_only=numeric_only_bool,
            )
            return result.__finalize__(self.obj, method="groupby")

    @final
    @Substitution(name="groupby")
    @Appender(_common_see_also)
    def median(self, numeric_only: bool | lib.NoDefault = lib.no_default):
        """
        Compute median of groups, excluding missing values.

        For multiple groupings, the result index will be a MultiIndex

        Parameters
        ----------
        numeric_only : bool, default True
            Include only float, int, boolean columns. If None, will attempt to use
            everything, then use only numeric data.

        Returns
        -------
        Series or DataFrame
            Median of values within each group.
        """
        numeric_only_bool = self._resolve_numeric_only(numeric_only)

        result = self._cython_agg_general(
            "median",
            alt=lambda x: Series(x).median(numeric_only=numeric_only_bool),
            numeric_only=numeric_only_bool,
        )
        return result.__finalize__(self.obj, method="groupby")

    @final
    @Substitution(name="groupby")
    @Appender(_common_see_also)
    def std(
        self,
        ddof: int = 1,
        engine: str | None = None,
        engine_kwargs: dict[str, bool] | None = None,
    ):
        """
        Compute standard deviation of groups, excluding missing values.

        For multiple groupings, the result index will be a MultiIndex.

        Parameters
        ----------
        ddof : int, default 1
            Degrees of freedom.

        engine : str, default None
            * ``'cython'`` : Runs the operation through C-extensions from cython.
            * ``'numba'`` : Runs the operation through JIT compiled code from numba.
            * ``None`` : Defaults to ``'cython'`` or globally setting
              ``compute.use_numba``

            .. versionadded:: 1.4.0

        engine_kwargs : dict, default None
            * For ``'cython'`` engine, there are no accepted ``engine_kwargs``
            * For ``'numba'`` engine, the engine can accept ``nopython``, ``nogil``
              and ``parallel`` dictionary keys. The values must either be ``True`` or
              ``False``. The default ``engine_kwargs`` for the ``'numba'`` engine is
              ``{{'nopython': True, 'nogil': False, 'parallel': False}}``

            .. versionadded:: 1.4.0

        Returns
        -------
        Series or DataFrame
            Standard deviation of values within each group.
        """
        if maybe_use_numba(engine):
            from pandas.core._numba.kernels import sliding_var

            return np.sqrt(
                self._numba_agg_general(sliding_var, engine_kwargs, "groupby_std", ddof)
            )
        else:
            return self._get_cythonized_result(
                libgroupby.group_var,
                needs_counts=True,
                cython_dtype=np.dtype(np.float64),
                post_processing=lambda vals, inference: np.sqrt(vals),
                ddof=ddof,
            )

    @final
    @Substitution(name="groupby")
    @Appender(_common_see_also)
    def var(
        self,
        ddof: int = 1,
        engine: str | None = None,
        engine_kwargs: dict[str, bool] | None = None,
    ):
        """
        Compute variance of groups, excluding missing values.

        For multiple groupings, the result index will be a MultiIndex.

        Parameters
        ----------
        ddof : int, default 1
            Degrees of freedom.

        engine : str, default None
            * ``'cython'`` : Runs the operation through C-extensions from cython.
            * ``'numba'`` : Runs the operation through JIT compiled code from numba.
            * ``None`` : Defaults to ``'cython'`` or globally setting
              ``compute.use_numba``

            .. versionadded:: 1.4.0

        engine_kwargs : dict, default None
            * For ``'cython'`` engine, there are no accepted ``engine_kwargs``
            * For ``'numba'`` engine, the engine can accept ``nopython``, ``nogil``
              and ``parallel`` dictionary keys. The values must either be ``True`` or
              ``False``. The default ``engine_kwargs`` for the ``'numba'`` engine is
              ``{{'nopython': True, 'nogil': False, 'parallel': False}}``

            .. versionadded:: 1.4.0

        Returns
        -------
        Series or DataFrame
            Variance of values within each group.
        """
        if maybe_use_numba(engine):
            from pandas.core._numba.kernels import sliding_var

            return self._numba_agg_general(
                sliding_var, engine_kwargs, "groupby_var", ddof
            )
        else:
            if ddof == 1:
                numeric_only = self._resolve_numeric_only(lib.no_default)
                return self._cython_agg_general(
                    "var",
                    alt=lambda x: Series(x).var(ddof=ddof),
                    numeric_only=numeric_only,
                )
            else:
                func = lambda x: x.var(ddof=ddof)
                with self._group_selection_context():
                    return self._python_agg_general(func)

    @final
    @Substitution(name="groupby")
    @Appender(_common_see_also)
    def sem(self, ddof: int = 1):
        """
        Compute standard error of the mean of groups, excluding missing values.

        For multiple groupings, the result index will be a MultiIndex.

        Parameters
        ----------
        ddof : int, default 1
            Degrees of freedom.

        Returns
        -------
        Series or DataFrame
            Standard error of the mean of values within each group.
        """
        result = self.std(ddof=ddof)
        if result.ndim == 1:
            result /= np.sqrt(self.count())
        else:
            cols = result.columns.difference(self.exclusions).unique()
            counts = self.count()
            result_ilocs = result.columns.get_indexer_for(cols)
            count_ilocs = counts.columns.get_indexer_for(cols)
            result.iloc[:, result_ilocs] /= np.sqrt(counts.iloc[:, count_ilocs])
        return result

    @final
    @Substitution(name="groupby")
    @Appender(_common_see_also)
    def size(self) -> DataFrame | Series:
        """
        Compute group sizes.

        Returns
        -------
        DataFrame or Series
            Number of rows in each group as a Series if as_index is True
            or a DataFrame if as_index is False.
        """
        result = self.grouper.size()

        # GH28330 preserve subclassed Series/DataFrames through calls
        if issubclass(self.obj._constructor, Series):
            result = self._obj_1d_constructor(result, name=self.obj.name)
        else:
            result = self._obj_1d_constructor(result)

        if not self.as_index:
            # Item "None" of "Optional[Series]" has no attribute "reset_index"
            result = result.rename("size").reset_index()  # type: ignore[union-attr]

        return self._reindex_output(result, fill_value=0)

    @final
    @doc(_groupby_agg_method_template, fname="sum", no=True, mc=0)
    def sum(
        self,
        numeric_only: bool | lib.NoDefault = lib.no_default,
        min_count: int = 0,
        engine: str | None = None,
        engine_kwargs: dict[str, bool] | None = None,
    ):
        if maybe_use_numba(engine):
            from pandas.core._numba.kernels import sliding_sum

            return self._numba_agg_general(
                sliding_sum,
                engine_kwargs,
                "groupby_sum",
            )
        else:
            numeric_only = self._resolve_numeric_only(numeric_only)

            # If we are grouping on categoricals we want unobserved categories to
            # return zero, rather than the default of NaN which the reindexing in
            # _agg_general() returns. GH #31422
            with com.temp_setattr(self, "observed", True):
                result = self._agg_general(
                    numeric_only=numeric_only,
                    min_count=min_count,
                    alias="add",
                    npfunc=np.sum,
                )

            return self._reindex_output(result, fill_value=0)

    @final
    @doc(_groupby_agg_method_template, fname="prod", no=True, mc=0)
    def prod(
        self, numeric_only: bool | lib.NoDefault = lib.no_default, min_count: int = 0
    ):
        numeric_only = self._resolve_numeric_only(numeric_only)

        return self._agg_general(
            numeric_only=numeric_only, min_count=min_count, alias="prod", npfunc=np.prod
        )

    @final
    @doc(_groupby_agg_method_template, fname="min", no=False, mc=-1)
    def min(self, numeric_only: bool = False, min_count: int = -1):
        return self._agg_general(
            numeric_only=numeric_only, min_count=min_count, alias="min", npfunc=np.min
        )

    @final
    @doc(_groupby_agg_method_template, fname="max", no=False, mc=-1)
    def max(self, numeric_only: bool = False, min_count: int = -1):
        return self._agg_general(
            numeric_only=numeric_only, min_count=min_count, alias="max", npfunc=np.max
        )

    @final
    @doc(_groupby_agg_method_template, fname="first", no=False, mc=-1)
    def first(self, numeric_only: bool = False, min_count: int = -1):
        def first_compat(obj: NDFrameT, axis: int = 0):
            def first(x: Series):
                """Helper function for first item that isn't NA."""
                arr = x.array[notna(x.array)]
                if not len(arr):
                    return np.nan
                return arr[0]

            if isinstance(obj, DataFrame):
                return obj.apply(first, axis=axis)
            elif isinstance(obj, Series):
                return first(obj)
            else:  # pragma: no cover
                raise TypeError(type(obj))

        return self._agg_general(
            numeric_only=numeric_only,
            min_count=min_count,
            alias="first",
            npfunc=first_compat,
        )

    @final
    @doc(_groupby_agg_method_template, fname="last", no=False, mc=-1)
    def last(self, numeric_only: bool = False, min_count: int = -1):
        def last_compat(obj: NDFrameT, axis: int = 0):
            def last(x: Series):
                """Helper function for last item that isn't NA."""
                arr = x.array[notna(x.array)]
                if not len(arr):
                    return np.nan
                return arr[-1]

            if isinstance(obj, DataFrame):
                return obj.apply(last, axis=axis)
            elif isinstance(obj, Series):
                return last(obj)
            else:  # pragma: no cover
                raise TypeError(type(obj))

        return self._agg_general(
            numeric_only=numeric_only,
            min_count=min_count,
            alias="last",
            npfunc=last_compat,
        )

    @final
    @Substitution(name="groupby")
    @Appender(_common_see_also)
    def ohlc(self) -> DataFrame:
        """
        Compute open, high, low and close values of a group, excluding missing values.

        For multiple groupings, the result index will be a MultiIndex

        Returns
        -------
        DataFrame
            Open, high, low and close values within each group.
        """
        if self.obj.ndim == 1:
            # self._iterate_slices() yields only self._selected_obj
            obj = self._selected_obj

            is_numeric = is_numeric_dtype(obj.dtype)
            if not is_numeric:
                raise DataError("No numeric types to aggregate")

            res_values = self.grouper._cython_operation(
                "aggregate", obj._values, "ohlc", axis=0, min_count=-1
            )

            agg_names = ["open", "high", "low", "close"]
            result = self.obj._constructor_expanddim(
                res_values, index=self.grouper.result_index, columns=agg_names
            )
            return self._reindex_output(result)

        return self._apply_to_column_groupbys(
            lambda x: x.ohlc(), self._obj_with_exclusions
        )

    @doc(DataFrame.describe)
    def describe(self, **kwargs):
        with self._group_selection_context():
            result = self.apply(lambda x: x.describe(**kwargs))
            if self.axis == 1:
                return result.T
            return result.unstack()

    @final
    def resample(self, rule, *args, **kwargs):
        """
        Provide resampling when using a TimeGrouper.

        Given a grouper, the function resamples it according to a string
        "string" -> "frequency".

        See the :ref:`frequency aliases <timeseries.offset_aliases>`
        documentation for more details.

        Parameters
        ----------
        rule : str or DateOffset
            The offset string or object representing target grouper conversion.
        *args, **kwargs
            Possible arguments are `how`, `fill_method`, `limit`, `kind` and
            `on`, and other arguments of `TimeGrouper`.

        Returns
        -------
        Grouper
            Return a new grouper with our resampler appended.

        See Also
        --------
        Grouper : Specify a frequency to resample with when
            grouping by a key.
        DatetimeIndex.resample : Frequency conversion and resampling of
            time series.

        Examples
        --------
        >>> idx = pd.date_range('1/1/2000', periods=4, freq='T')
        >>> df = pd.DataFrame(data=4 * [range(2)],
        ...                   index=idx,
        ...                   columns=['a', 'b'])
        >>> df.iloc[2, 0] = 5
        >>> df
                            a  b
        2000-01-01 00:00:00  0  1
        2000-01-01 00:01:00  0  1
        2000-01-01 00:02:00  5  1
        2000-01-01 00:03:00  0  1

        Downsample the DataFrame into 3 minute bins and sum the values of
        the timestamps falling into a bin.

        >>> df.groupby('a').resample('3T').sum()
                                 a  b
        a
        0   2000-01-01 00:00:00  0  2
            2000-01-01 00:03:00  0  1
        5   2000-01-01 00:00:00  5  1

        Upsample the series into 30 second bins.

        >>> df.groupby('a').resample('30S').sum()
                            a  b
        a
        0   2000-01-01 00:00:00  0  1
            2000-01-01 00:00:30  0  0
            2000-01-01 00:01:00  0  1
            2000-01-01 00:01:30  0  0
            2000-01-01 00:02:00  0  0
            2000-01-01 00:02:30  0  0
            2000-01-01 00:03:00  0  1
        5   2000-01-01 00:02:00  5  1

        Resample by month. Values are assigned to the month of the period.

        >>> df.groupby('a').resample('M').sum()
                    a  b
        a
        0   2000-01-31  0  3
        5   2000-01-31  5  1

        Downsample the series into 3 minute bins as above, but close the right
        side of the bin interval.

        >>> df.groupby('a').resample('3T', closed='right').sum()
                                 a  b
        a
        0   1999-12-31 23:57:00  0  1
            2000-01-01 00:00:00  0  2
        5   2000-01-01 00:00:00  5  1

        Downsample the series into 3 minute bins and close the right side of
        the bin interval, but label each bin using the right edge instead of
        the left.

        >>> df.groupby('a').resample('3T', closed='right', label='right').sum()
                                 a  b
        a
        0   2000-01-01 00:00:00  0  1
            2000-01-01 00:03:00  0  2
        5   2000-01-01 00:03:00  5  1
        """
        from pandas.core.resample import get_resampler_for_grouping

        return get_resampler_for_grouping(self, rule, *args, **kwargs)

    @final
    @Substitution(name="groupby")
    @Appender(_common_see_also)
    def rolling(self, *args, **kwargs):
        """
        Return a rolling grouper, providing rolling functionality per group.
        """
        from pandas.core.window import RollingGroupby

        return RollingGroupby(
            self._selected_obj,
            *args,
            _grouper=self.grouper,
            _as_index=self.as_index,
            **kwargs,
        )

    @final
    @Substitution(name="groupby")
    @Appender(_common_see_also)
    def expanding(self, *args, **kwargs):
        """
        Return an expanding grouper, providing expanding
        functionality per group.
        """
        from pandas.core.window import ExpandingGroupby

        return ExpandingGroupby(
            self._selected_obj,
            *args,
            _grouper=self.grouper,
            **kwargs,
        )

    @final
    @Substitution(name="groupby")
    @Appender(_common_see_also)
    def ewm(self, *args, **kwargs):
        """
        Return an ewm grouper, providing ewm functionality per group.
        """
        from pandas.core.window import ExponentialMovingWindowGroupby

        return ExponentialMovingWindowGroupby(
            self._selected_obj,
            *args,
            _grouper=self.grouper,
            **kwargs,
        )

    @final
    def _fill(self, direction: Literal["ffill", "bfill"], limit=None):
        """
        Shared function for `pad` and `backfill` to call Cython method.

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
        pad : Returns Series with minimum number of char in object.
        backfill : Backward fill the missing values in the dataset.
        """
        # Need int value for Cython
        if limit is None:
            limit = -1

        ids, _, _ = self.grouper.group_info
        sorted_labels = np.argsort(ids, kind="mergesort").astype(np.intp, copy=False)
        if direction == "bfill":
            sorted_labels = sorted_labels[::-1]

        col_func = partial(
            libgroupby.group_fillna_indexer,
            labels=ids,
            sorted_labels=sorted_labels,
            direction=direction,
            limit=limit,
            dropna=self.dropna,
        )

        def blk_func(values: ArrayLike) -> ArrayLike:
            mask = isna(values)
            if values.ndim == 1:
                indexer = np.empty(values.shape, dtype=np.intp)
                col_func(out=indexer, mask=mask)
                return algorithms.take_nd(values, indexer)

            else:
                # We broadcast algorithms.take_nd analogous to
                #  np.take_along_axis

                # Note: we only get here with backfill/pad,
                #  so if we have a dtype that cannot hold NAs,
                #  then there will be no -1s in indexer, so we can use
                #  the original dtype (no need to ensure_dtype_can_hold_na)
                if isinstance(values, np.ndarray):
                    out = np.empty(values.shape, dtype=values.dtype)
                else:
                    out = type(values)._empty(values.shape, dtype=values.dtype)

                for i in range(len(values)):
                    # call group_fillna_indexer column-wise
                    indexer = np.empty(values.shape[1], dtype=np.intp)
                    col_func(out=indexer, mask=mask[i])
                    out[i, :] = algorithms.take_nd(values[i], indexer)
                return out

        obj = self._obj_with_exclusions
        if self.axis == 1:
            obj = obj.T
        mgr = obj._mgr
        res_mgr = mgr.apply(blk_func)

        new_obj = obj._constructor(res_mgr)
        if isinstance(new_obj, Series):
            new_obj.name = obj.name

        return self._wrap_transformed_output(new_obj)

    @final
    @Substitution(name="groupby")
    def ffill(self, limit=None):
        """
        Forward fill the values.

        Parameters
        ----------
        limit : int, optional
            Limit of how many values to fill.

        Returns
        -------
        Series or DataFrame
            Object with missing values filled.

        See Also
        --------
        Series.ffill: Returns Series with minimum number of char in object.
        DataFrame.ffill: Object with missing values filled or None if inplace=True.
        Series.fillna: Fill NaN values of a Series.
        DataFrame.fillna: Fill NaN values of a DataFrame.
        """
        return self._fill("ffill", limit=limit)

    def pad(self, limit=None):
        warnings.warn(
            "pad is deprecated and will be removed in a future version. "
            "Use ffill instead.",
            FutureWarning,
            stacklevel=find_stack_level(),
        )
        return self.ffill(limit=limit)

    pad.__doc__ = ffill.__doc__

    @final
    @Substitution(name="groupby")
    def bfill(self, limit=None):
        """
        Backward fill the values.

        Parameters
        ----------
        limit : int, optional
            Limit of how many values to fill.

        Returns
        -------
        Series or DataFrame
            Object with missing values filled.

        See Also
        --------
        Series.bfill :  Backward fill the missing values in the dataset.
        DataFrame.bfill:  Backward fill the missing values in the dataset.
        Series.fillna: Fill NaN values of a Series.
        DataFrame.fillna: Fill NaN values of a DataFrame.
        """
        return self._fill("bfill", limit=limit)

    def backfill(self, limit=None):
        warnings.warn(
            "backfill is deprecated and will be removed in a future version. "
            "Use bfill instead.",
            FutureWarning,
            stacklevel=find_stack_level(),
        )
        return self.bfill(limit=limit)

    backfill.__doc__ = bfill.__doc__

    @final
    @Substitution(name="groupby")
    @Substitution(see_also=_common_see_also)
    def nth(
        self,
        n: PositionalIndexer | tuple,
        dropna: Literal["any", "all", None] = None,
    ) -> NDFrameT:
        """
        Take the nth row from each group if n is an int, otherwise a subset of rows.

        Can be either a call or an index. dropna is not available with index notation.
        Index notation accepts a comma separated list of integers and slices.

        If dropna, will take the nth non-null row, dropna is either
        'all' or 'any'; this is equivalent to calling dropna(how=dropna)
        before the groupby.

        Parameters
        ----------
        n : int, slice or list of ints and slices
            A single nth value for the row or a list of nth values or slices.

            .. versionchanged:: 1.4.0
                Added slice and lists containiing slices.
                Added index notation.

        dropna : {'any', 'all', None}, default None
            Apply the specified dropna operation before counting which row is
            the nth row. Only supported if n is an int.

        Returns
        -------
        Series or DataFrame
            N-th value within each group.
        %(see_also)s
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
        >>> g.nth(slice(None, -1))
             B
        A
        1  NaN
        1  2.0
        2  3.0

        Index notation may also be used

        >>> g.nth[0, 1]
             B
        A
        1  NaN
        1  2.0
        2  3.0
        2  5.0
        >>> g.nth[:-1]
             B
        A
        1  NaN
        1  2.0
        2  3.0

        Specifying `dropna` allows count ignoring ``NaN``

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

        Specifying `as_index=False` in `groupby` keeps the original index.

        >>> df.groupby('A', as_index=False).nth(1)
           A    B
        1  1  2.0
        4  2  5.0
        """
        if not dropna:
            with self._group_selection_context():
                mask = self._make_mask_from_positional_indexer(n)

                ids, _, _ = self.grouper.group_info

                # Drop NA values in grouping
                mask = mask & (ids != -1)

                out = self._mask_selected_obj(mask)
                if not self.as_index:
                    return out

                result_index = self.grouper.result_index
                if self.axis == 0:
                    out.index = result_index[ids[mask]]
                    if not self.observed and isinstance(result_index, CategoricalIndex):
                        out = out.reindex(result_index)

                    out = self._reindex_output(out)
                else:
                    out.columns = result_index[ids[mask]]

                return out.sort_index(axis=self.axis) if self.sort else out

        # dropna is truthy
        if not is_integer(n):
            raise ValueError("dropna option only supported for an integer argument")

        if dropna not in ["any", "all"]:
            # Note: when agg-ing picker doesn't raise this, just returns NaN
            raise ValueError(
                "For a DataFrame or Series groupby.nth, dropna must be "
                "either None, 'any' or 'all', "
                f"(was passed {dropna})."
            )

        # old behaviour, but with all and any support for DataFrames.
        # modified in GH 7559 to have better perf
        n = cast(int, n)
        max_len = n if n >= 0 else -1 - n
        dropped = self.obj.dropna(how=dropna, axis=self.axis)

        # get a new grouper for our dropped obj
        if self.keys is None and self.level is None:

            # we don't have the grouper info available
            # (e.g. we have selected out
            # a column that is not in the current object)
            axis = self.grouper.axis
            grouper = axis[axis.isin(dropped.index)]

        else:

            # create a grouper with the original parameters, but on dropped
            # object
            from pandas.core.groupby.grouper import get_grouper

            grouper, _, _ = get_grouper(
                dropped,
                key=self.keys,
                axis=self.axis,
                level=self.level,
                sort=self.sort,
                mutated=self.mutated,
            )

        grb = dropped.groupby(
            grouper, as_index=self.as_index, sort=self.sort, axis=self.axis
        )
        sizes, result = grb.size(), grb.nth(n)
        mask = (sizes < max_len)._values

        # set the results which don't meet the criteria
        if len(result) and mask.any():
            result.loc[mask] = np.nan

        # reset/reindex to the original groups
        if len(self.obj) == len(dropped) or len(result) == len(
            self.grouper.result_index
        ):
            result.index = self.grouper.result_index
        else:
            result = result.reindex(self.grouper.result_index)

        return result

    @final
    def quantile(self, q=0.5, interpolation: str = "linear"):
        """
        Return group values at the given quantile, a la numpy.percentile.

        Parameters
        ----------
        q : float or array-like, default 0.5 (50% quantile)
            Value(s) between 0 and 1 providing the quantile(s) to compute.
        interpolation : {'linear', 'lower', 'higher', 'midpoint', 'nearest'}
            Method to use when the desired quantile falls between two points.

        Returns
        -------
        Series or DataFrame
            Return type determined by caller of GroupBy object.

        See Also
        --------
        Series.quantile : Similar method for Series.
        DataFrame.quantile : Similar method for DataFrame.
        numpy.percentile : NumPy method to compute qth percentile.

        Examples
        --------
        >>> df = pd.DataFrame([
        ...     ['a', 1], ['a', 2], ['a', 3],
        ...     ['b', 1], ['b', 3], ['b', 5]
        ... ], columns=['key', 'val'])
        >>> df.groupby('key').quantile()
            val
        key
        a    2.0
        b    3.0
        """

        def pre_processor(vals: ArrayLike) -> tuple[np.ndarray, np.dtype | None]:
            if is_object_dtype(vals):
                raise TypeError(
                    "'quantile' cannot be performed against 'object' dtypes!"
                )

            inference: np.dtype | None = None
            if is_integer_dtype(vals.dtype):
                if isinstance(vals, ExtensionArray):
                    out = vals.to_numpy(dtype=float, na_value=np.nan)
                else:
                    out = vals
                inference = np.dtype(np.int64)
            elif is_bool_dtype(vals.dtype) and isinstance(vals, ExtensionArray):
                out = vals.to_numpy(dtype=float, na_value=np.nan)
            elif is_datetime64_dtype(vals.dtype):
                inference = np.dtype("datetime64[ns]")
                out = np.asarray(vals).astype(float)
            elif is_timedelta64_dtype(vals.dtype):
                inference = np.dtype("timedelta64[ns]")
                out = np.asarray(vals).astype(float)
            elif isinstance(vals, ExtensionArray) and is_float_dtype(vals):
                inference = np.dtype(np.float64)
                out = vals.to_numpy(dtype=float, na_value=np.nan)
            else:
                out = np.asarray(vals)

            return out, inference

        def post_processor(vals: np.ndarray, inference: np.dtype | None) -> np.ndarray:
            if inference:
                # Check for edge case
                if not (
                    is_integer_dtype(inference)
                    and interpolation in {"linear", "midpoint"}
                ):
                    vals = vals.astype(inference)

            return vals

        orig_scalar = is_scalar(q)
        if orig_scalar:
            q = [q]

        qs = np.array(q, dtype=np.float64)
        ids, _, ngroups = self.grouper.group_info
        nqs = len(qs)

        func = partial(
            libgroupby.group_quantile, labels=ids, qs=qs, interpolation=interpolation
        )

        # Put '-1' (NaN) labels as the last group so it does not interfere
        # with the calculations. Note: length check avoids failure on empty
        # labels. In that case, the value doesn't matter
        na_label_for_sorting = ids.max() + 1 if len(ids) > 0 else 0
        labels_for_lexsort = np.where(ids == -1, na_label_for_sorting, ids)

        def blk_func(values: ArrayLike) -> ArrayLike:
            mask = isna(values)
            vals, inference = pre_processor(values)

            ncols = 1
            if vals.ndim == 2:
                ncols = vals.shape[0]
                shaped_labels = np.broadcast_to(
                    labels_for_lexsort, (ncols, len(labels_for_lexsort))
                )
            else:
                shaped_labels = labels_for_lexsort

            out = np.empty((ncols, ngroups, nqs), dtype=np.float64)

            # Get an index of values sorted by values and then labels
            order = (vals, shaped_labels)
            sort_arr = np.lexsort(order).astype(np.intp, copy=False)

            if vals.ndim == 1:
                func(out[0], values=vals, mask=mask, sort_indexer=sort_arr)
            else:
                for i in range(ncols):
                    func(out[i], values=vals[i], mask=mask[i], sort_indexer=sort_arr[i])

            if vals.ndim == 1:
                out = out.ravel("K")
            else:
                out = out.reshape(ncols, ngroups * nqs)
            return post_processor(out, inference)

        obj = self._obj_with_exclusions
        is_ser = obj.ndim == 1
        mgr = self._get_data_to_aggregate()

        res_mgr = mgr.grouped_reduce(blk_func, ignore_failures=True)
        if not is_ser and len(res_mgr.items) != len(mgr.items):
            warn_dropping_nuisance_columns_deprecated(type(self), "quantile")

            if len(res_mgr.items) == 0:
                # re-call grouped_reduce to get the desired exception message
                mgr.grouped_reduce(blk_func, ignore_failures=False)
                # grouped_reduce _should_ raise, so this should not be reached
                raise TypeError(  # pragma: no cover
                    "All columns were dropped in grouped_reduce"
                )

        if is_ser:
            res = self._wrap_agged_manager(res_mgr)
        else:
            res = obj._constructor(res_mgr)

        if orig_scalar:
            # Avoid expensive MultiIndex construction
            return self._wrap_aggregated_output(res)
        return self._wrap_aggregated_output(res, qs=qs)

    @final
    @Substitution(name="groupby")
    def ngroup(self, ascending: bool = True):
        """
        Number each group from 0 to the number of groups - 1.

        This is the enumerative complement of cumcount.  Note that the
        numbers given to the groups match the order in which the groups
        would be seen when iterating over the groupby object, not the
        order they are first observed.

        Parameters
        ----------
        ascending : bool, default True
            If False, number in reverse, from number of group - 1 to 0.

        Returns
        -------
        Series
            Unique numbers for each group.

        See Also
        --------
        .cumcount : Number the rows in each group.

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
        """
        with self._group_selection_context():
            index = self._selected_obj.index
            result = self._obj_1d_constructor(
                self.grouper.group_info[0], index, dtype=np.int64
            )
            if not ascending:
                result = self.ngroups - 1 - result
            return result

    @final
    @Substitution(name="groupby")
    def cumcount(self, ascending: bool = True):
        """
        Number each item in each group from 0 to the length of that group - 1.

        Essentially this is equivalent to

        .. code-block:: python

            self.apply(lambda x: pd.Series(np.arange(len(x)), x.index))

        Parameters
        ----------
        ascending : bool, default True
            If False, number in reverse, from length of group - 1 to 0.

        Returns
        -------
        Series
            Sequence number of each element within each group.

        See Also
        --------
        .ngroup : Number the groups themselves.

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
        """
        with self._group_selection_context():
            index = self._selected_obj._get_axis(self.axis)
            cumcounts = self._cumcount_array(ascending=ascending)
            return self._obj_1d_constructor(cumcounts, index)

    @final
    @Substitution(name="groupby")
    @Substitution(see_also=_common_see_also)
    def rank(
        self,
        method: str = "average",
        ascending: bool = True,
        na_option: str = "keep",
        pct: bool = False,
        axis: int = 0,
    ):
        """
        Provide the rank of values within each group.

        Parameters
        ----------
        method : {'average', 'min', 'max', 'first', 'dense'}, default 'average'
            * average: average rank of group.
            * min: lowest rank in group.
            * max: highest rank in group.
            * first: ranks assigned in order they appear in the array.
            * dense: like 'min', but rank always increases by 1 between groups.
        ascending : bool, default True
            False for ranks by high (1) to low (N).
        na_option : {'keep', 'top', 'bottom'}, default 'keep'
            * keep: leave NA values where they are.
            * top: smallest rank if ascending.
            * bottom: smallest rank if descending.
        pct : bool, default False
            Compute percentage rank of data within each group.
        axis : int, default 0
            The axis of the object over which to compute the rank.

        Returns
        -------
        DataFrame with ranking of values within each group
        %(see_also)s
        Examples
        --------
        >>> df = pd.DataFrame(
        ...     {
        ...         "group": ["a", "a", "a", "a", "a", "b", "b", "b", "b", "b"],
        ...         "value": [2, 4, 2, 3, 5, 1, 2, 4, 1, 5],
        ...     }
        ... )
        >>> df
          group  value
        0     a      2
        1     a      4
        2     a      2
        3     a      3
        4     a      5
        5     b      1
        6     b      2
        7     b      4
        8     b      1
        9     b      5
        >>> for method in ['average', 'min', 'max', 'dense', 'first']:
        ...     df[f'{method}_rank'] = df.groupby('group')['value'].rank(method)
        >>> df
          group  value  average_rank  min_rank  max_rank  dense_rank  first_rank
        0     a      2           1.5       1.0       2.0         1.0         1.0
        1     a      4           4.0       4.0       4.0         3.0         4.0
        2     a      2           1.5       1.0       2.0         1.0         2.0
        3     a      3           3.0       3.0       3.0         2.0         3.0
        4     a      5           5.0       5.0       5.0         4.0         5.0
        5     b      1           1.5       1.0       2.0         1.0         1.0
        6     b      2           3.0       3.0       3.0         2.0         3.0
        7     b      4           4.0       4.0       4.0         3.0         4.0
        8     b      1           1.5       1.0       2.0         1.0         2.0
        9     b      5           5.0       5.0       5.0         4.0         5.0
        """
        if na_option not in {"keep", "top", "bottom"}:
            msg = "na_option must be one of 'keep', 'top', or 'bottom'"
            raise ValueError(msg)

        kwargs = {
            "ties_method": method,
            "ascending": ascending,
            "na_option": na_option,
            "pct": pct,
        }
        if axis != 0:
            # DataFrame uses different keyword name
            kwargs["method"] = kwargs.pop("ties_method")
            return self.apply(lambda x: x.rank(axis=axis, numeric_only=False, **kwargs))

        return self._cython_transform(
            "rank",
            numeric_only=False,
            axis=axis,
            **kwargs,
        )

    @final
    @Substitution(name="groupby")
    @Appender(_common_see_also)
    def cumprod(self, axis=0, *args, **kwargs):
        """
        Cumulative product for each group.

        Returns
        -------
        Series or DataFrame
        """
        nv.validate_groupby_func("cumprod", args, kwargs, ["numeric_only", "skipna"])
        if axis != 0:
            return self.apply(lambda x: x.cumprod(axis=axis, **kwargs))

        return self._cython_transform("cumprod", **kwargs)

    @final
    @Substitution(name="groupby")
    @Appender(_common_see_also)
    def cumsum(self, axis=0, *args, **kwargs):
        """
        Cumulative sum for each group.

        Returns
        -------
        Series or DataFrame
        """
        nv.validate_groupby_func("cumsum", args, kwargs, ["numeric_only", "skipna"])
        if axis != 0:
            return self.apply(lambda x: x.cumsum(axis=axis, **kwargs))

        return self._cython_transform("cumsum", **kwargs)

    @final
    @Substitution(name="groupby")
    @Appender(_common_see_also)
    def cummin(self, axis=0, **kwargs):
        """
        Cumulative min for each group.

        Returns
        -------
        Series or DataFrame
        """
        skipna = kwargs.get("skipna", True)
        if axis != 0:
            return self.apply(lambda x: np.minimum.accumulate(x, axis))

        return self._cython_transform("cummin", numeric_only=False, skipna=skipna)

    @final
    @Substitution(name="groupby")
    @Appender(_common_see_also)
    def cummax(self, axis=0, **kwargs):
        """
        Cumulative max for each group.

        Returns
        -------
        Series or DataFrame
        """
        skipna = kwargs.get("skipna", True)
        if axis != 0:
            return self.apply(lambda x: np.maximum.accumulate(x, axis))

        return self._cython_transform("cummax", numeric_only=False, skipna=skipna)

    @final
    def _get_cythonized_result(
        self,
        base_func: Callable,
        cython_dtype: np.dtype,
        numeric_only: bool | lib.NoDefault = lib.no_default,
        needs_counts: bool = False,
        needs_nullable: bool = False,
        needs_mask: bool = False,
        pre_processing=None,
        post_processing=None,
        **kwargs,
    ):
        """
        Get result for Cythonized functions.

        Parameters
        ----------
        base_func : callable, Cythonized function to be called
        cython_dtype : np.dtype
            Type of the array that will be modified by the Cython call.
        numeric_only : bool, default True
            Whether only numeric datatypes should be computed
        needs_counts : bool, default False
            Whether the counts should be a part of the Cython call
        needs_mask : bool, default False
            Whether boolean mask needs to be part of the Cython call
            signature
        needs_nullable : bool, default False
            Whether a bool specifying if the input is nullable is part
            of the Cython call signature
        pre_processing : function, default None
            Function to be applied to `values` prior to passing to Cython.
            Function should return a tuple where the first element is the
            values to be passed to Cython and the second element is an optional
            type which the values should be converted to after being returned
            by the Cython operation. This function is also responsible for
            raising a TypeError if the values have an invalid type. Raises
            if `needs_values` is False.
        post_processing : function, default None
            Function to be applied to result of Cython function. Should accept
            an array of values as the first argument and type inferences as its
            second argument, i.e. the signature should be
            (ndarray, Type). If `needs_nullable=True`, a third argument should be
            `nullable`, to allow for processing specific to nullable values.
        **kwargs : dict
            Extra arguments to be passed back to Cython funcs

        Returns
        -------
        `Series` or `DataFrame`  with filled values
        """
        numeric_only = self._resolve_numeric_only(numeric_only)

        if post_processing and not callable(post_processing):
            raise ValueError("'post_processing' must be a callable!")
        if pre_processing and not callable(pre_processing):
            raise ValueError("'pre_processing' must be a callable!")

        grouper = self.grouper

        ids, _, ngroups = grouper.group_info

        how = base_func.__name__
        base_func = partial(base_func, labels=ids)

        def blk_func(values: ArrayLike) -> ArrayLike:
            values = values.T
            ncols = 1 if values.ndim == 1 else values.shape[1]

            result: ArrayLike
            result = np.zeros(ngroups * ncols, dtype=cython_dtype)
            result = result.reshape((ngroups, ncols))

            func = partial(base_func, out=result)

            inferences = None

            if needs_counts:
                counts = np.zeros(self.ngroups, dtype=np.int64)
                func = partial(func, counts=counts)

            vals = values
            if pre_processing:
                vals, inferences = pre_processing(vals)

            vals = vals.astype(cython_dtype, copy=False)
            if vals.ndim == 1:
                vals = vals.reshape((-1, 1))
            func = partial(func, values=vals)

            if needs_mask:
                mask = isna(values).view(np.uint8)
                if mask.ndim == 1:
                    mask = mask.reshape(-1, 1)
                func = partial(func, mask=mask)

            if needs_nullable:
                is_nullable = isinstance(values, BaseMaskedArray)
                func = partial(func, nullable=is_nullable)

            func(**kwargs)  # Call func to modify indexer values in place

            if values.ndim == 1:
                assert result.shape[1] == 1, result.shape
                result = result[:, 0]

            if post_processing:
                pp_kwargs = {}
                if needs_nullable:
                    pp_kwargs["nullable"] = isinstance(values, BaseMaskedArray)

                result = post_processing(result, inferences, **pp_kwargs)

            return result.T

        obj = self._obj_with_exclusions

        # Operate block-wise instead of column-by-column
        is_ser = obj.ndim == 1
        mgr = self._get_data_to_aggregate()

        if numeric_only:
            mgr = mgr.get_numeric_data()

        res_mgr = mgr.grouped_reduce(blk_func, ignore_failures=True)

        if not is_ser and len(res_mgr.items) != len(mgr.items):
            howstr = how.replace("group_", "")
            warn_dropping_nuisance_columns_deprecated(type(self), howstr)

            if len(res_mgr.items) == 0:
                # We re-call grouped_reduce to get the right exception message
                mgr.grouped_reduce(blk_func, ignore_failures=False)
                # grouped_reduce _should_ raise, so this should not be reached
                raise TypeError(  # pragma: no cover
                    "All columns were dropped in grouped_reduce"
                )

        if is_ser:
            out = self._wrap_agged_manager(res_mgr)
        else:
            out = obj._constructor(res_mgr)

        return self._wrap_aggregated_output(out)

    @final
    @Substitution(name="groupby")
    def shift(self, periods=1, freq=None, axis=0, fill_value=None):
        """
        Shift each group by periods observations.

        If freq is passed, the index will be increased using the periods and the freq.

        Parameters
        ----------
        periods : int, default 1
            Number of periods to shift.
        freq : str, optional
            Frequency string.
        axis : axis to shift, default 0
            Shift direction.
        fill_value : optional
            The scalar value to use for newly introduced missing values.

        Returns
        -------
        Series or DataFrame
            Object shifted within each group.

        See Also
        --------
        Index.shift : Shift values of Index.
        tshift : Shift the time index, using the index’s frequency
            if available.
        """
        if freq is not None or axis != 0:
            return self.apply(lambda x: x.shift(periods, freq, axis, fill_value))

        ids, _, ngroups = self.grouper.group_info
        res_indexer = np.zeros(len(ids), dtype=np.int64)

        libgroupby.group_shift_indexer(res_indexer, ids, ngroups, periods)

        obj = self._obj_with_exclusions

        res = obj._reindex_with_indexers(
            {self.axis: (obj.axes[self.axis], res_indexer)},
            fill_value=fill_value,
            allow_dups=True,
        )
        return res

    @final
    @Substitution(name="groupby")
    @Appender(_common_see_also)
    def pct_change(self, periods=1, fill_method="ffill", limit=None, freq=None, axis=0):
        """
        Calculate pct_change of each value to previous entry in group.

        Returns
        -------
        Series or DataFrame
            Percentage changes within each group.
        """
        # TODO(GH#23918): Remove this conditional for SeriesGroupBy when
        #  GH#23918 is fixed
        if freq is not None or axis != 0:
            return self.apply(
                lambda x: x.pct_change(
                    periods=periods,
                    fill_method=fill_method,
                    limit=limit,
                    freq=freq,
                    axis=axis,
                )
            )
        if fill_method is None:  # GH30463
            fill_method = "ffill"
            limit = 0
        filled = getattr(self, fill_method)(limit=limit)
        fill_grp = filled.groupby(self.grouper.codes, axis=self.axis)
        shifted = fill_grp.shift(periods=periods, freq=freq, axis=self.axis)
        return (filled / shifted) - 1

    @final
    @Substitution(name="groupby")
    @Substitution(see_also=_common_see_also)
    def head(self, n=5):
        """
        Return first n rows of each group.

        Similar to ``.apply(lambda x: x.head(n))``, but it returns a subset of rows
        from the original DataFrame with original index and order preserved
        (``as_index`` flag is ignored).

        Parameters
        ----------
        n : int
            If positive: number of entries to include from start of each group.
            If negative: number of entries to exclude from end of each group.

        Returns
        -------
        Series or DataFrame
            Subset of original Series or DataFrame as determined by n.
        %(see_also)s
        Examples
        --------

        >>> df = pd.DataFrame([[1, 2], [1, 4], [5, 6]],
        ...                   columns=['A', 'B'])
        >>> df.groupby('A').head(1)
           A  B
        0  1  2
        2  5  6
        >>> df.groupby('A').head(-1)
           A  B
        0  1  2
        """
        self._reset_group_selection()
        mask = self._make_mask_from_positional_indexer(slice(None, n))
        return self._mask_selected_obj(mask)

    @final
    @Substitution(name="groupby")
    @Substitution(see_also=_common_see_also)
    def tail(self, n=5):
        """
        Return last n rows of each group.

        Similar to ``.apply(lambda x: x.tail(n))``, but it returns a subset of rows
        from the original DataFrame with original index and order preserved
        (``as_index`` flag is ignored).

        Parameters
        ----------
        n : int
            If positive: number of entries to include from end of each group.
            If negative: number of entries to exclude from start of each group.

        Returns
        -------
        Series or DataFrame
            Subset of original Series or DataFrame as determined by n.
        %(see_also)s
        Examples
        --------

        >>> df = pd.DataFrame([['a', 1], ['a', 2], ['b', 1], ['b', 2]],
        ...                   columns=['A', 'B'])
        >>> df.groupby('A').tail(1)
           A  B
        1  a  2
        3  b  2
        >>> df.groupby('A').tail(-1)
           A  B
        1  a  2
        3  b  2
        """
        self._reset_group_selection()
        if n:
            mask = self._make_mask_from_positional_indexer(slice(-n, None))
        else:
            mask = self._make_mask_from_positional_indexer([])

        return self._mask_selected_obj(mask)

    @final
    def _mask_selected_obj(self, mask: np.ndarray) -> NDFrameT:
        """
        Return _selected_obj with mask applied to the correct axis.

        Parameters
        ----------
        mask : np.ndarray
            Boolean mask to apply.

        Returns
        -------
        Series or DataFrame
            Filtered _selected_obj.
        """
        ids = self.grouper.group_info[0]
        mask = mask & (ids != -1)

        if self.axis == 0:
            return self._selected_obj[mask]
        else:
            return self._selected_obj.iloc[:, mask]

    @final
    def _reindex_output(
        self,
        output: OutputFrameOrSeries,
        fill_value: Scalar = np.NaN,
        qs: npt.NDArray[np.float64] | None = None,
    ) -> OutputFrameOrSeries:
        """
        If we have categorical groupers, then we might want to make sure that
        we have a fully re-indexed output to the levels. This means expanding
        the output space to accommodate all values in the cartesian product of
        our groups, regardless of whether they were observed in the data or
        not. This will expand the output space if there are missing groups.

        The method returns early without modifying the input if the number of
        groupings is less than 2, self.observed == True or none of the groupers
        are categorical.

        Parameters
        ----------
        output : Series or DataFrame
            Object resulting from grouping and applying an operation.
        fill_value : scalar, default np.NaN
            Value to use for unobserved categories if self.observed is False.
        qs : np.ndarray[float64] or None, default None
            quantile values, only relevant for quantile.

        Returns
        -------
        Series or DataFrame
            Object (potentially) re-indexed to include all possible groups.
        """
        groupings = self.grouper.groupings
        if len(groupings) == 1:
            return output

        # if we only care about the observed values
        # we are done
        elif self.observed:
            return output

        # reindexing only applies to a Categorical grouper
        elif not any(
            isinstance(ping.grouping_vector, (Categorical, CategoricalIndex))
            for ping in groupings
        ):
            return output

        levels_list = [ping.group_index for ping in groupings]
        names = self.grouper.names
        if qs is not None:
            # error: Argument 1 to "append" of "list" has incompatible type
            # "ndarray[Any, dtype[floating[_64Bit]]]"; expected "Index"
            levels_list.append(qs)  # type: ignore[arg-type]
            names = names + [None]
        index, _ = MultiIndex.from_product(levels_list, names=names).sortlevel()

        if self.as_index:
            d = {
                self.obj._get_axis_name(self.axis): index,
                "copy": False,
                "fill_value": fill_value,
            }
            return output.reindex(**d)

        # GH 13204
        # Here, the categorical in-axis groupers, which need to be fully
        # expanded, are columns in `output`. An idea is to do:
        # output = output.set_index(self.grouper.names)
        #                .reindex(index).reset_index()
        # but special care has to be taken because of possible not-in-axis
        # groupers.
        # So, we manually select and drop the in-axis grouper columns,
        # reindex `output`, and then reset the in-axis grouper columns.

        # Select in-axis groupers
        in_axis_grps = (
            (i, ping.name) for (i, ping) in enumerate(groupings) if ping.in_axis
        )
        g_nums, g_names = zip(*in_axis_grps)

        output = output.drop(labels=list(g_names), axis=1)

        # Set a temp index and reindex (possibly expanding)
        output = output.set_index(self.grouper.result_index).reindex(
            index, copy=False, fill_value=fill_value
        )

        # Reset in-axis grouper columns
        # (using level numbers `g_nums` because level names may not be unique)
        output = output.reset_index(level=g_nums)

        return output.reset_index(drop=True)

    @final
    def sample(
        self,
        n: int | None = None,
        frac: float | None = None,
        replace: bool = False,
        weights: Sequence | Series | None = None,
        random_state: RandomState | None = None,
    ):
        """
        Return a random sample of items from each group.

        You can use `random_state` for reproducibility.

        .. versionadded:: 1.1.0

        Parameters
        ----------
        n : int, optional
            Number of items to return for each group. Cannot be used with
            `frac` and must be no larger than the smallest group unless
            `replace` is True. Default is one if `frac` is None.
        frac : float, optional
            Fraction of items to return. Cannot be used with `n`.
        replace : bool, default False
            Allow or disallow sampling of the same row more than once.
        weights : list-like, optional
            Default None results in equal probability weighting.
            If passed a list-like then values must have the same length as
            the underlying DataFrame or Series object and will be used as
            sampling probabilities after normalization within each group.
            Values must be non-negative with at least one positive element
            within each group.
        random_state : int, array-like, BitGenerator, np.random.RandomState, np.random.Generator, optional
            If int, array-like, or BitGenerator, seed for random number generator.
            If np.random.RandomState or np.random.Generator, use as given.

            .. versionchanged:: 1.4.0

                np.random.Generator objects now accepted

        Returns
        -------
        Series or DataFrame
            A new object of same type as caller containing items randomly
            sampled within each group from the caller object.

        See Also
        --------
        DataFrame.sample: Generate random samples from a DataFrame object.
        numpy.random.choice: Generate a random sample from a given 1-D numpy
            array.

        Examples
        --------
        >>> df = pd.DataFrame(
        ...     {"a": ["red"] * 2 + ["blue"] * 2 + ["black"] * 2, "b": range(6)}
        ... )
        >>> df
               a  b
        0    red  0
        1    red  1
        2   blue  2
        3   blue  3
        4  black  4
        5  black  5

        Select one row at random for each distinct value in column a. The
        `random_state` argument can be used to guarantee reproducibility:

        >>> df.groupby("a").sample(n=1, random_state=1)
               a  b
        4  black  4
        2   blue  2
        1    red  1

        Set `frac` to sample fixed proportions rather than counts:

        >>> df.groupby("a")["b"].sample(frac=0.5, random_state=2)
        5    5
        2    2
        0    0
        Name: b, dtype: int64

        Control sample probabilities within groups by setting weights:

        >>> df.groupby("a").sample(
        ...     n=1,
        ...     weights=[1, 1, 1, 0, 0, 1],
        ...     random_state=1,
        ... )
               a  b
        5  black  5
        2   blue  2
        0    red  0
        """  # noqa:E501
        size = sample.process_sampling_size(n, frac, replace)
        if weights is not None:
            weights_arr = sample.preprocess_weights(
                self._selected_obj, weights, axis=self.axis
            )

        random_state = com.random_state(random_state)

        group_iterator = self.grouper.get_iterator(self._selected_obj, self.axis)

        sampled_indices = []
        for labels, obj in group_iterator:
            grp_indices = self.indices[labels]
            group_size = len(grp_indices)
            if size is not None:
                sample_size = size
            else:
                assert frac is not None
                sample_size = round(frac * group_size)

            grp_sample = sample.sample(
                group_size,
                size=sample_size,
                replace=replace,
                weights=None if weights is None else weights_arr[grp_indices],
                random_state=random_state,
            )
            sampled_indices.append(grp_indices[grp_sample])

        sampled_indices = np.concatenate(sampled_indices)
        return self._selected_obj.take(sampled_indices, axis=self.axis)


@doc(GroupBy)
def get_groupby(
    obj: NDFrame,
    by: _KeysArgType | None = None,
    axis: int = 0,
    level=None,
    grouper: ops.BaseGrouper | None = None,
    exclusions=None,
    selection=None,
    as_index: bool = True,
    sort: bool = True,
    group_keys: bool = True,
    squeeze: bool = False,
    observed: bool = False,
    mutated: bool = False,
    dropna: bool = True,
) -> GroupBy:

    klass: type[GroupBy]
    if isinstance(obj, Series):
        from pandas.core.groupby.generic import SeriesGroupBy

        klass = SeriesGroupBy
    elif isinstance(obj, DataFrame):
        from pandas.core.groupby.generic import DataFrameGroupBy

        klass = DataFrameGroupBy
    else:  # pragma: no cover
        raise TypeError(f"invalid type: {obj}")

    return klass(
        obj=obj,
        keys=by,
        axis=axis,
        level=level,
        grouper=grouper,
        exclusions=exclusions,
        selection=selection,
        as_index=as_index,
        sort=sort,
        group_keys=group_keys,
        squeeze=squeeze,
        observed=observed,
        mutated=mutated,
        dropna=dropna,
    )


def _insert_quantile_level(idx: Index, qs: npt.NDArray[np.float64]) -> MultiIndex:
    """
    Insert the sequence 'qs' of quantiles as the inner-most level of a MultiIndex.

    The quantile level in the MultiIndex is a repeated copy of 'qs'.

    Parameters
    ----------
    idx : Index
    qs : np.ndarray[float64]

    Returns
    -------
    MultiIndex
    """
    nqs = len(qs)

    if idx._is_multi:
        idx = cast(MultiIndex, idx)
        lev_codes, lev = Index(qs).factorize()
        levels = list(idx.levels) + [lev]
        codes = [np.repeat(x, nqs) for x in idx.codes] + [np.tile(lev_codes, len(idx))]
        mi = MultiIndex(levels=levels, codes=codes, names=idx.names + [None])
    else:
        mi = MultiIndex.from_product([idx, qs])
    return mi


def warn_dropping_nuisance_columns_deprecated(cls, how: str) -> None:
    warnings.warn(
        "Dropping invalid columns in "
        f"{cls.__name__}.{how} is deprecated. "
        "In a future version, a TypeError will be raised. "
        f"Before calling .{how}, select only columns which "
        "should be valid for the function.",
        FutureWarning,
        stacklevel=find_stack_level(),
    )
