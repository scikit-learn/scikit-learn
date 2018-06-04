from __future__ import absolute_import, division, print_function

from collections import Iterator
from functools import wraps, partial
import operator
from operator import getitem
from pprint import pformat
import warnings

from toolz import merge, first, unique, partition_all, remove
import pandas as pd
import numpy as np
from numbers import Number

try:
    from chest import Chest as Cache
except ImportError:
    Cache = dict

from .. import array as da
from .. import core

from ..utils import partial_by_order
from .. import threaded
from ..compatibility import apply, operator_div, bind_method, string_types
from ..context import globalmethod
from ..utils import (random_state_data, pseudorandom, derived_from, funcname,
                     memory_repr, put_lines, M, key_split, OperatorMethodMixin,
                     is_arraylike)
from ..array import Array
from ..base import DaskMethodsMixin, tokenize, dont_optimize, is_dask_collection
from ..delayed import Delayed, to_task_dask

from . import methods
from .accessor import DatetimeAccessor, StringAccessor
from .categorical import CategoricalAccessor, categorize
from .hashing import hash_pandas_object
from .optimize import optimize
from .utils import (meta_nonempty, make_meta, insert_meta_param_description,
                    raise_on_meta_error, clear_known_categories,
                    is_categorical_dtype, has_known_categories, PANDAS_VERSION,
                    index_summary)

no_default = '__no_default__'

if PANDAS_VERSION >= '0.20.0':
    from pandas.util import cache_readonly
    pd.set_option('compute.use_numexpr', False)
else:
    from pandas.util.decorators import cache_readonly
    pd.computation.expressions.set_use_numexpr(False)


def _concat(args):
    if not args:
        return args
    if isinstance(first(core.flatten(args)), np.ndarray):
        return da.core.concatenate3(args)
    if not isinstance(args[0], (pd.DataFrame, pd.Series, pd.Index)):
        try:
            return pd.Series(args)
        except Exception:
            return args
    # We filter out empty partitions here because pandas frequently has
    # inconsistent dtypes in results between empty and non-empty frames.
    # Ideally this would be handled locally for each operation, but in practice
    # this seems easier. TODO: don't do this.
    args2 = [i for i in args if len(i)]
    return args[0] if not args2 else methods.concat(args2, uniform=True)


def _get_return_type(meta):
    if isinstance(meta, _Frame):
        meta = meta._meta

    if isinstance(meta, pd.Series):
        return Series
    elif isinstance(meta, pd.DataFrame):
        return DataFrame
    elif isinstance(meta, pd.Index):
        return Index
    return Scalar


def new_dd_object(dsk, name, meta, divisions):
    """Generic constructor for dask.dataframe objects.

    Decides the appropriate output class based on the type of `meta` provided.
    """
    if isinstance(meta, (pd.Series, pd.DataFrame, pd.Index)):
        return _get_return_type(meta)(dsk, name, meta, divisions)
    elif is_arraylike(meta):
        import dask.array as da
        chunks = (((np.nan,) * (len(divisions) - 1),) +
                  tuple((d,) for d in meta.shape[1:]))
        if len(chunks) > 1:
            dsk = dsk.copy()
            suffix = (0,) * (len(chunks) - 1)
            for i in range(len(chunks[0])):
                dsk[(name, i) + suffix] = dsk.pop((name, i))
        return da.Array(dsk, name=name, chunks=chunks, dtype=meta.dtype)
    else:
        return _get_return_type(meta)(dsk, name, meta, divisions)


def finalize(results):
    return _concat(results)


class Scalar(DaskMethodsMixin, OperatorMethodMixin):
    """ A Dask object to represent a pandas scalar"""
    def __init__(self, dsk, name, meta, divisions=None):
        # divisions is ignored, only present to be compatible with other
        # objects.
        self.dask = dsk
        self._name = name
        meta = make_meta(meta)
        if isinstance(meta, (pd.DataFrame, pd.Series, pd.Index)):
            raise TypeError("Expected meta to specify scalar, got "
                            "{0}".format(type(meta).__name__))
        self._meta = meta

    def __dask_graph__(self):
        return self.dask

    def __dask_keys__(self):
        return [self.key]

    def __dask_tokenize__(self):
        return self._name

    __dask_optimize__ = globalmethod(optimize, key='dataframe_optimize',
                                     falsey=dont_optimize)
    __dask_scheduler__ = staticmethod(threaded.get)

    def __dask_postcompute__(self):
        return first, ()

    def __dask_postpersist__(self):
        return Scalar, (self._name, self._meta, self.divisions)

    @property
    def _meta_nonempty(self):
        return self._meta

    @property
    def dtype(self):
        return self._meta.dtype

    def __dir__(self):
        o = set(dir(type(self)))
        o.update(self.__dict__)
        if not hasattr(self._meta, 'dtype'):
            o.remove('dtype')  # dtype only in `dir` if available
        return list(o)

    @property
    def divisions(self):
        """Dummy divisions to be compat with Series and DataFrame"""
        return [None, None]

    def __repr__(self):
        name = self._name if len(self._name) < 10 else self._name[:7] + '...'
        if hasattr(self._meta, 'dtype'):
            extra = ', dtype=%s' % self._meta.dtype
        else:
            extra = ', type=%s' % type(self._meta).__name__
        return "dd.Scalar<%s%s>" % (name, extra)

    def __array__(self):
        # array interface is required to support pandas instance + Scalar
        # Otherwise, above op results in pd.Series of Scalar (object dtype)
        return np.asarray(self.compute())

    @property
    def _args(self):
        return (self.dask, self._name, self._meta)

    def __getstate__(self):
        return self._args

    def __setstate__(self, state):
        self.dask, self._name, self._meta = state

    @property
    def key(self):
        return (self._name, 0)

    @classmethod
    def _get_unary_operator(cls, op):
        def f(self):
            name = funcname(op) + '-' + tokenize(self)
            dsk = {(name, 0): (op, (self._name, 0))}
            meta = op(self._meta_nonempty)
            return Scalar(merge(dsk, self.dask), name, meta)
        return f

    @classmethod
    def _get_binary_operator(cls, op, inv=False):
        return lambda self, other: _scalar_binary(op, self, other, inv=inv)

    def to_delayed(self, optimize_graph=True):
        """Convert into a ``dask.delayed`` object.

        Parameters
        ----------
        optimize_graph : bool, optional
            If True [default], the graph is optimized before converting into
            ``dask.delayed`` objects.
        """
        from dask.delayed import Delayed
        dsk = self.__dask_graph__()
        if optimize_graph:
            dsk = self.__dask_optimize__(dsk, self.__dask_keys__())
        return Delayed(self.key, dsk)


def _scalar_binary(op, self, other, inv=False):
    name = '{0}-{1}'.format(funcname(op), tokenize(self, other))

    dsk = self.dask
    return_type = _get_return_type(other)

    if isinstance(other, Scalar):
        dsk = merge(dsk, other.dask)
        other_key = (other._name, 0)
    elif is_dask_collection(other):
        return NotImplemented
    else:
        other_key = other

    if inv:
        dsk.update({(name, 0): (op, other_key, (self._name, 0))})
    else:
        dsk.update({(name, 0): (op, (self._name, 0), other_key)})

    other_meta = make_meta(other)
    other_meta_nonempty = meta_nonempty(other_meta)
    if inv:
        meta = op(other_meta_nonempty, self._meta_nonempty)
    else:
        meta = op(self._meta_nonempty, other_meta_nonempty)

    if return_type is not Scalar:
        return return_type(dsk, name, meta,
                           [other.index.min(), other.index.max()])
    else:
        return Scalar(dsk, name, meta)


class _Frame(DaskMethodsMixin, OperatorMethodMixin):
    """ Superclass for DataFrame and Series

    Parameters
    ----------

    dsk: dict
        The dask graph to compute this DataFrame
    name: str
        The key prefix that specifies which keys in the dask comprise this
        particular DataFrame / Series
    meta: pandas.DataFrame, pandas.Series, or pandas.Index
        An empty pandas object with names, dtypes, and indices matching the
        expected output.
    divisions: tuple of index values
        Values along which we partition our blocks on the index
    """
    def __init__(self, dsk, name, meta, divisions):
        self.dask = dsk
        self._name = name
        meta = make_meta(meta)
        if not isinstance(meta, self._partition_type):
            raise TypeError("Expected meta to specify type {0}, got type "
                            "{1}".format(self._partition_type.__name__,
                                         type(meta).__name__))
        self._meta = meta
        self.divisions = tuple(divisions)

    def __dask_graph__(self):
        return self.dask

    def __dask_keys__(self):
        return [(self._name, i) for i in range(self.npartitions)]

    def __dask_tokenize__(self):
        return self._name

    __dask_optimize__ = globalmethod(optimize, key='dataframe_optimize',
                                     falsey=dont_optimize)
    __dask_scheduler__ = staticmethod(threaded.get)

    def __dask_postcompute__(self):
        return finalize, ()

    def __dask_postpersist__(self):
        return type(self), (self._name, self._meta, self.divisions)

    @property
    def _constructor(self):
        return new_dd_object

    @property
    def npartitions(self):
        """Return number of partitions"""
        return len(self.divisions) - 1

    @property
    def size(self):
        """ Size of the series """
        return self.reduction(methods.size, np.sum, token='size', meta=int,
                              split_every=False)

    @property
    def _meta_nonempty(self):
        """ A non-empty version of `_meta` with fake data."""
        return meta_nonempty(self._meta)

    @property
    def _args(self):
        return (self.dask, self._name, self._meta, self.divisions)

    def __getstate__(self):
        return self._args

    def __setstate__(self, state):
        self.dask, self._name, self._meta, self.divisions = state

    def copy(self):
        """ Make a copy of the dataframe

        This is strictly a shallow copy of the underlying computational graph.
        It does not affect the underlying data
        """
        return new_dd_object(self.dask, self._name,
                             self._meta, self.divisions)

    def __array__(self, dtype=None, **kwargs):
        self._computed = self.compute()
        x = np.array(self._computed)
        return x

    def __array_wrap__(self, array, context=None):
        raise NotImplementedError

    def __array_ufunc__(self, numpy_ufunc, method, *inputs, **kwargs):
        out = kwargs.get('out', ())
        for x in inputs + out:
            # ufuncs work with 0-dimensional NumPy ndarrays
            # so we don't want to raise NotImplemented
            if isinstance(x, np.ndarray) and x.shape == ():
                continue
            elif not isinstance(x, (Number, Scalar, _Frame, Array,
                                    pd.DataFrame, pd.Series, pd.Index)):
                return NotImplemented

        if method == '__call__':
            if numpy_ufunc.signature is not None:
                return NotImplemented
            if numpy_ufunc.nout > 1:
                # ufuncs with multiple output values
                # are not yet supported for frames
                return NotImplemented
            else:
                return elemwise(numpy_ufunc, *inputs, **kwargs)
        else:
            # ufunc methods are not yet supported for frames
            return NotImplemented

    @property
    def _elemwise(self):
        return elemwise

    @property
    def _repr_data(self):
        raise NotImplementedError

    @property
    def _repr_divisions(self):
        name = "npartitions={0}".format(self.npartitions)
        if self.known_divisions:
            divisions = pd.Index(self.divisions, name=name)
        else:
            # avoid to be converted to NaN
            divisions = pd.Index([''] * (self.npartitions + 1), name=name)
        return divisions

    def __repr__(self):
        data = self._repr_data.to_string(max_rows=5, show_dimensions=False)
        return """Dask {klass} Structure:
{data}
Dask Name: {name}, {task} tasks""".format(klass=self.__class__.__name__,
                                          data=data, name=key_split(self._name),
                                          task=len(self.dask))

    @property
    def index(self):
        """Return dask Index instance"""
        name = self._name + '-index'
        dsk = {(name, i): (getattr, key, 'index')
               for i, key in enumerate(self.__dask_keys__())}

        return Index(merge(dsk, self.dask), name,
                     self._meta.index, self.divisions)

    def reset_index(self, drop=False):
        """Reset the index to the default index.

        Note that unlike in ``pandas``, the reset ``dask.dataframe`` index will
        not be monotonically increasing from 0. Instead, it will restart at 0
        for each partition (e.g. ``index1 = [0, ..., 10], index2 = [0, ...]``).
        This is due to the inability to statically know the full length of the
        index.

        For DataFrame with multi-level index, returns a new DataFrame with
        labeling information in the columns under the index names, defaulting
        to 'level_0', 'level_1', etc. if any are None. For a standard index,
        the index name will be used (if set), otherwise a default 'index' or
        'level_0' (if 'index' is already taken) will be used.

        Parameters
        ----------
        drop : boolean, default False
            Do not try to insert index into dataframe columns.
        """
        return self.map_partitions(M.reset_index, drop=drop).clear_divisions()

    @property
    def known_divisions(self):
        """Whether divisions are already known"""
        return len(self.divisions) > 0 and self.divisions[0] is not None

    def clear_divisions(self):
        """ Forget division information """
        divisions = (None,) * (self.npartitions + 1)
        return type(self)(self.dask, self._name, self._meta, divisions)

    def get_partition(self, n):
        """Get a dask DataFrame/Series representing the `nth` partition."""
        if 0 <= n < self.npartitions:
            name = 'get-partition-%s-%s' % (str(n), self._name)
            dsk = {(name, 0): (self._name, n)}
            divisions = self.divisions[n:n + 2]
            return new_dd_object(merge(self.dask, dsk), name,
                                 self._meta, divisions)
        else:
            msg = "n must be 0 <= n < {0}".format(self.npartitions)
            raise ValueError(msg)

    @derived_from(pd.DataFrame)
    def drop_duplicates(self, split_every=None, split_out=1, **kwargs):
        # Let pandas error on bad inputs
        self._meta_nonempty.drop_duplicates(**kwargs)
        if 'subset' in kwargs and kwargs['subset'] is not None:
            split_out_setup = split_out_on_cols
            split_out_setup_kwargs = {'cols': kwargs['subset']}
        else:
            split_out_setup = split_out_setup_kwargs = None

        if kwargs.get('keep', True) is False:
            raise NotImplementedError("drop_duplicates with keep=False")

        chunk = M.drop_duplicates
        return aca(self, chunk=chunk, aggregate=chunk, meta=self._meta,
                   token='drop-duplicates', split_every=split_every,
                   split_out=split_out, split_out_setup=split_out_setup,
                   split_out_setup_kwargs=split_out_setup_kwargs, **kwargs)

    def __len__(self):
        return self.reduction(len, np.sum, token='len', meta=int,
                              split_every=False).compute()

    def __bool__(self):
        raise ValueError("The truth value of a {0} is ambiguous. "
                         "Use a.any() or a.all()."
                         .format(self.__class__.__name__))

    __nonzero__ = __bool__  # python 2

    def _scalarfunc(self, cast_type):
        def wrapper():
            raise TypeError("cannot convert the series to "
                            "{0}".format(str(cast_type)))

        return wrapper

    def __float__(self):
        return self._scalarfunc(float)

    def __int__(self):
        return self._scalarfunc(int)

    __long__ = __int__  # python 2

    def __complex__(self):
        return self._scalarfunc(complex)

    @insert_meta_param_description(pad=12)
    def map_partitions(self, func, *args, **kwargs):
        """ Apply Python function on each DataFrame partition.

        Note that the index and divisions are assumed to remain unchanged.

        Parameters
        ----------
        func : function
            Function applied to each partition.
        args, kwargs :
            Arguments and keywords to pass to the function. The partition will
            be the first argument, and these will be passed *after*. Arguments
            and keywords may contain ``Scalar``, ``Delayed`` or regular
            python objects.
        $META

        Examples
        --------
        Given a DataFrame, Series, or Index, such as:

        >>> import dask.dataframe as dd
        >>> df = pd.DataFrame({'x': [1, 2, 3, 4, 5],
        ...                    'y': [1., 2., 3., 4., 5.]})
        >>> ddf = dd.from_pandas(df, npartitions=2)

        One can use ``map_partitions`` to apply a function on each partition.
        Extra arguments and keywords can optionally be provided, and will be
        passed to the function after the partition.

        Here we apply a function with arguments and keywords to a DataFrame,
        resulting in a Series:

        >>> def myadd(df, a, b=1):
        ...     return df.x + df.y + a + b
        >>> res = ddf.map_partitions(myadd, 1, b=2)
        >>> res.dtype
        dtype('float64')

        By default, dask tries to infer the output metadata by running your
        provided function on some fake data. This works well in many cases, but
        can sometimes be expensive, or even fail. To avoid this, you can
        manually specify the output metadata with the ``meta`` keyword. This
        can be specified in many forms, for more information see
        ``dask.dataframe.utils.make_meta``.

        Here we specify the output is a Series with no name, and dtype
        ``float64``:

        >>> res = ddf.map_partitions(myadd, 1, b=2, meta=(None, 'f8'))

        Here we map a function that takes in a DataFrame, and returns a
        DataFrame with a new column:

        >>> res = ddf.map_partitions(lambda df: df.assign(z=df.x * df.y))
        >>> res.dtypes
        x      int64
        y    float64
        z    float64
        dtype: object

        As before, the output metadata can also be specified manually. This
        time we pass in a ``dict``, as the output is a DataFrame:

        >>> res = ddf.map_partitions(lambda df: df.assign(z=df.x * df.y),
        ...                          meta={'x': 'i8', 'y': 'f8', 'z': 'f8'})

        In the case where the metadata doesn't change, you can also pass in
        the object itself directly:

        >>> res = ddf.map_partitions(lambda df: df.head(), meta=df)

        Also note that the index and divisions are assumed to remain unchanged.
        If the function you're mapping changes the index/divisions, you'll need
        to clear them afterwards:

        >>> ddf.map_partitions(func).clear_divisions()  # doctest: +SKIP
        """
        return map_partitions(func, self, *args, **kwargs)

    @insert_meta_param_description(pad=12)
    def map_overlap(self, func, before, after, *args, **kwargs):
        """Apply a function to each partition, sharing rows with adjacent partitions.

        This can be useful for implementing windowing functions such as
        ``df.rolling(...).mean()`` or ``df.diff()``.

        Parameters
        ----------
        func : function
            Function applied to each partition.
        before : int
            The number of rows to prepend to partition ``i`` from the end of
            partition ``i - 1``.
        after : int
            The number of rows to append to partition ``i`` from the beginning
            of partition ``i + 1``.
        args, kwargs :
            Arguments and keywords to pass to the function. The partition will
            be the first argument, and these will be passed *after*.
        $META

        Notes
        -----
        Given positive integers ``before`` and ``after``, and a function
        ``func``, ``map_overlap`` does the following:

        1. Prepend ``before`` rows to each partition ``i`` from the end of
           partition ``i - 1``. The first partition has no rows prepended.

        2. Append ``after`` rows to each partition ``i`` from the beginning of
           partition ``i + 1``. The last partition has no rows appended.

        3. Apply ``func`` to each partition, passing in any extra ``args`` and
           ``kwargs`` if provided.

        4. Trim ``before`` rows from the beginning of all but the first
           partition.

        5. Trim ``after`` rows from the end of all but the last partition.

        Note that the index and divisions are assumed to remain unchanged.

        Examples
        --------
        Given a DataFrame, Series, or Index, such as:

        >>> import dask.dataframe as dd
        >>> df = pd.DataFrame({'x': [1, 2, 4, 7, 11],
        ...                    'y': [1., 2., 3., 4., 5.]})
        >>> ddf = dd.from_pandas(df, npartitions=2)

        A rolling sum with a trailing moving window of size 2 can be computed by
        overlapping 2 rows before each partition, and then mapping calls to
        ``df.rolling(2).sum()``:

        >>> ddf.compute()
            x    y
        0   1  1.0
        1   2  2.0
        2   4  3.0
        3   7  4.0
        4  11  5.0
        >>> ddf.map_overlap(lambda df: df.rolling(2).sum(), 2, 0).compute()
              x    y
        0   NaN  NaN
        1   3.0  3.0
        2   6.0  5.0
        3  11.0  7.0
        4  18.0  9.0

        The pandas ``diff`` method computes a discrete difference shifted by a
        number of periods (can be positive or negative). This can be
        implemented by mapping calls to ``df.diff`` to each partition after
        prepending/appending that many rows, depending on sign:

        >>> def diff(df, periods=1):
        ...     before, after = (periods, 0) if periods > 0 else (0, -periods)
        ...     return df.map_overlap(lambda df, periods=1: df.diff(periods),
        ...                           periods, 0, periods=periods)
        >>> diff(ddf, 1).compute()
             x    y
        0  NaN  NaN
        1  1.0  1.0
        2  2.0  1.0
        3  3.0  1.0
        4  4.0  1.0

        If you have a ``DatetimeIndex``, you can use a ``pd.Timedelta`` for time-
        based windows.

        >>> ts = pd.Series(range(10), index=pd.date_range('2017', periods=10))
        >>> dts = dd.from_pandas(ts, npartitions=2)
        >>> dts.map_overlap(lambda df: df.rolling('2D').sum(),
        ...                 pd.Timedelta('2D'), 0).compute()
        2017-01-01     0.0
        2017-01-02     1.0
        2017-01-03     3.0
        2017-01-04     5.0
        2017-01-05     7.0
        2017-01-06     9.0
        2017-01-07    11.0
        2017-01-08    13.0
        2017-01-09    15.0
        2017-01-10    17.0
        dtype: float64
        """
        from .rolling import map_overlap
        return map_overlap(func, self, before, after, *args, **kwargs)

    @insert_meta_param_description(pad=12)
    def reduction(self, chunk, aggregate=None, combine=None, meta=no_default,
                  token=None, split_every=None, chunk_kwargs=None,
                  aggregate_kwargs=None, combine_kwargs=None, **kwargs):
        """Generic row-wise reductions.

        Parameters
        ----------
        chunk : callable
            Function to operate on each partition. Should return a
            ``pandas.DataFrame``, ``pandas.Series``, or a scalar.
        aggregate : callable, optional
            Function to operate on the concatenated result of ``chunk``. If not
            specified, defaults to ``chunk``. Used to do the final aggregation
            in a tree reduction.

            The input to ``aggregate`` depends on the output of ``chunk``.
            If the output of ``chunk`` is a:

            - scalar: Input is a Series, with one row per partition.
            - Series: Input is a DataFrame, with one row per partition. Columns
              are the rows in the output series.
            - DataFrame: Input is a DataFrame, with one row per partition.
              Columns are the columns in the output dataframes.

            Should return a ``pandas.DataFrame``, ``pandas.Series``, or a
            scalar.
        combine : callable, optional
            Function to operate on intermediate concatenated results of
            ``chunk`` in a tree-reduction. If not provided, defaults to
            ``aggregate``. The input/output requirements should match that of
            ``aggregate`` described above.
        $META
        token : str, optional
            The name to use for the output keys.
        split_every : int, optional
            Group partitions into groups of this size while performing a
            tree-reduction. If set to False, no tree-reduction will be used,
            and all intermediates will be concatenated and passed to
            ``aggregate``. Default is 8.
        chunk_kwargs : dict, optional
            Keyword arguments to pass on to ``chunk`` only.
        aggregate_kwargs : dict, optional
            Keyword arguments to pass on to ``aggregate`` only.
        combine_kwargs : dict, optional
            Keyword arguments to pass on to ``combine`` only.
        kwargs :
            All remaining keywords will be passed to ``chunk``, ``combine``,
            and ``aggregate``.

        Examples
        --------
        >>> import pandas as pd
        >>> import dask.dataframe as dd
        >>> df = pd.DataFrame({'x': range(50), 'y': range(50, 100)})
        >>> ddf = dd.from_pandas(df, npartitions=4)

        Count the number of rows in a DataFrame. To do this, count the number
        of rows in each partition, then sum the results:

        >>> res = ddf.reduction(lambda x: x.count(),
        ...                     aggregate=lambda x: x.sum())
        >>> res.compute()
        x    50
        y    50
        dtype: int64

        Count the number of rows in a Series with elements greater than or
        equal to a value (provided via a keyword).

        >>> def count_greater(x, value=0):
        ...     return (x >= value).sum()
        >>> res = ddf.x.reduction(count_greater, aggregate=lambda x: x.sum(),
        ...                       chunk_kwargs={'value': 25})
        >>> res.compute()
        25

        Aggregate both the sum and count of a Series at the same time:

        >>> def sum_and_count(x):
        ...     return pd.Series({'sum': x.sum(), 'count': x.count()})
        >>> res = ddf.x.reduction(sum_and_count, aggregate=lambda x: x.sum())
        >>> res.compute()
        count      50
        sum      1225
        dtype: int64

        Doing the same, but for a DataFrame. Here ``chunk`` returns a
        DataFrame, meaning the input to ``aggregate`` is a DataFrame with an
        index with non-unique entries for both 'x' and 'y'. We groupby the
        index, and sum each group to get the final result.

        >>> def sum_and_count(x):
        ...     return pd.DataFrame({'sum': x.sum(), 'count': x.count()})
        >>> res = ddf.reduction(sum_and_count,
        ...                     aggregate=lambda x: x.groupby(level=0).sum())
        >>> res.compute()
           count   sum
        x     50  1225
        y     50  3725
        """
        if aggregate is None:
            aggregate = chunk

        if combine is None:
            if combine_kwargs:
                raise ValueError("`combine_kwargs` provided with no `combine`")
            combine = aggregate
            combine_kwargs = aggregate_kwargs

        chunk_kwargs = chunk_kwargs.copy() if chunk_kwargs else {}
        chunk_kwargs['aca_chunk'] = chunk

        combine_kwargs = combine_kwargs.copy() if combine_kwargs else {}
        combine_kwargs['aca_combine'] = combine

        aggregate_kwargs = aggregate_kwargs.copy() if aggregate_kwargs else {}
        aggregate_kwargs['aca_aggregate'] = aggregate

        return aca(self, chunk=_reduction_chunk, aggregate=_reduction_aggregate,
                   combine=_reduction_combine, meta=meta, token=token,
                   split_every=split_every, chunk_kwargs=chunk_kwargs,
                   aggregate_kwargs=aggregate_kwargs,
                   combine_kwargs=combine_kwargs, **kwargs)

    @derived_from(pd.DataFrame)
    def pipe(self, func, *args, **kwargs):
        # Taken from pandas:
        # https://github.com/pydata/pandas/blob/master/pandas/core/generic.py#L2698-L2707
        if isinstance(func, tuple):
            func, target = func
            if target in kwargs:
                raise ValueError('%s is both the pipe target and a keyword '
                                 'argument' % target)
            kwargs[target] = self
            return func(*args, **kwargs)
        else:
            return func(self, *args, **kwargs)

    def random_split(self, frac, random_state=None):
        """ Pseudorandomly split dataframe into different pieces row-wise

        Parameters
        ----------
        frac : list
            List of floats that should sum to one.
        random_state: int or np.random.RandomState
            If int create a new RandomState with this as the seed
        Otherwise draw from the passed RandomState

        Examples
        --------

        50/50 split

        >>> a, b = df.random_split([0.5, 0.5])  # doctest: +SKIP

        80/10/10 split, consistent random_state

        >>> a, b, c = df.random_split([0.8, 0.1, 0.1], random_state=123)  # doctest: +SKIP

        See Also
        --------
        dask.DataFrame.sample
        """
        if not np.allclose(sum(frac), 1):
            raise ValueError("frac should sum to 1")
        state_data = random_state_data(self.npartitions, random_state)
        token = tokenize(self, frac, random_state)
        name = 'split-' + token
        dsk = {(name, i): (pd_split, (self._name, i), frac, state)
               for i, state in enumerate(state_data)}

        out = []
        for i in range(len(frac)):
            name2 = 'split-%d-%s' % (i, token)
            dsk2 = {(name2, j): (getitem, (name, j), i)
                    for j in range(self.npartitions)}
            out.append(type(self)(merge(self.dask, dsk, dsk2), name2,
                                  self._meta, self.divisions))
        return out

    def head(self, n=5, npartitions=1, compute=True):
        """ First n rows of the dataset

        Parameters
        ----------
        n : int, optional
            The number of rows to return. Default is 5.
        npartitions : int, optional
            Elements are only taken from the first ``npartitions``, with a
            default of 1. If there are fewer than ``n`` rows in the first
            ``npartitions`` a warning will be raised and any found rows
            returned. Pass -1 to use all partitions.
        compute : bool, optional
            Whether to compute the result, default is True.
        """
        if npartitions <= -1:
            npartitions = self.npartitions
        if npartitions > self.npartitions:
            msg = "only {} partitions, head received {}"
            raise ValueError(msg.format(self.npartitions, npartitions))

        name = 'head-%d-%d-%s' % (npartitions, n, self._name)

        if npartitions > 1:
            name_p = 'head-partial-%d-%s' % (n, self._name)

            dsk = {}
            for i in range(npartitions):
                dsk[(name_p, i)] = (M.head, (self._name, i), n)

            concat = (_concat, [(name_p, i) for i in range(npartitions)])
            dsk[(name, 0)] = (safe_head, concat, n)
        else:
            dsk = {(name, 0): (safe_head, (self._name, 0), n)}

        result = new_dd_object(merge(self.dask, dsk), name, self._meta,
                               [self.divisions[0], self.divisions[npartitions]])

        if compute:
            result = result.compute()
        return result

    def tail(self, n=5, compute=True):
        """ Last n rows of the dataset

        Caveat, the only checks the last n rows of the last partition.
        """
        name = 'tail-%d-%s' % (n, self._name)
        dsk = {(name, 0): (M.tail, (self._name, self.npartitions - 1), n)}

        result = new_dd_object(merge(self.dask, dsk), name,
                               self._meta, self.divisions[-2:])

        if compute:
            result = result.compute()
        return result

    @property
    def loc(self):
        """ Purely label-location based indexer for selection by label.

        >>> df.loc["b"]  # doctest: +SKIP
        >>> df.loc["b":"d"]  # doctest: +SKIP"""
        from .indexing import _LocIndexer
        return _LocIndexer(self)

    # NOTE: `iloc` is not implemented because of performance concerns.
    # see https://github.com/dask/dask/pull/507

    def repartition(self, divisions=None, npartitions=None, freq=None, force=False):
        """ Repartition dataframe along new divisions

        Parameters
        ----------
        divisions : list, optional
            List of partitions to be used. If specified npartitions will be
            ignored.
        npartitions : int, optional
            Number of partitions of output. Only used if divisions isn't
            specified.
        freq : str, pd.Timedelta
            A period on which to partition timeseries data like ``'7D'`` or
            ``'12h'`` or ``pd.Timedelta(hours=12)``.  Assumes a datetime index.
        force : bool, default False
            Allows the expansion of the existing divisions.
            If False then the new divisions lower and upper bounds must be
            the same as the old divisions.

        Examples
        --------
        >>> df = df.repartition(npartitions=10)  # doctest: +SKIP
        >>> df = df.repartition(divisions=[0, 5, 10, 20])  # doctest: +SKIP
        >>> df = df.repartition(freq='7d')  # doctest: +SKIP
        """
        if npartitions is not None and divisions is not None:
            warnings.warn("When providing both npartitions and divisions to "
                          "repartition only npartitions is used.")

        if npartitions is not None:
            return repartition_npartitions(self, npartitions)
        elif divisions is not None:
            return repartition(self, divisions, force=force)
        elif freq is not None:
            return repartition_freq(self, freq=freq)
        else:
            raise ValueError(
                "Provide either divisions= or npartitions= to repartition")

    @derived_from(pd.DataFrame)
    def fillna(self, value=None, method=None, limit=None, axis=None):
        axis = self._validate_axis(axis)
        if method is None and limit is not None:
            raise NotImplementedError("fillna with set limit and method=None")
        if isinstance(value, _Frame):
            test_value = value._meta_nonempty.values[0]
        else:
            test_value = value
        meta = self._meta_nonempty.fillna(value=test_value, method=method,
                                          limit=limit, axis=axis)

        if axis == 1 or method is None:
            # Control whether or not dask's partition alignment happens.
            # We don't want for a pandas Series.
            # We do want it for a dask Series
            if isinstance(value, pd.Series):
                args = ()
                kwargs = {'value': value}
            else:
                args = (value,)
                kwargs = {}
            return self.map_partitions(M.fillna, *args, method=method,
                                       limit=limit, axis=axis, meta=meta,
                                       **kwargs)

        if method in ('pad', 'ffill'):
            method = 'ffill'
            skip_check = 0
            before, after = 1 if limit is None else limit, 0
        else:
            method = 'bfill'
            skip_check = self.npartitions - 1
            before, after = 0, 1 if limit is None else limit

        if limit is None:
            name = 'fillna-chunk-' + tokenize(self, method)
            dsk = {(name, i): (methods.fillna_check, (self._name, i),
                               method, i != skip_check)
                   for i in range(self.npartitions)}
            parts = new_dd_object(merge(dsk, self.dask), name, meta,
                                  self.divisions)
        else:
            parts = self

        return parts.map_overlap(M.fillna, before, after, method=method,
                                 limit=limit, meta=meta)

    @derived_from(pd.DataFrame)
    def ffill(self, axis=None, limit=None):
        return self.fillna(method='ffill', limit=limit, axis=axis)

    @derived_from(pd.DataFrame)
    def bfill(self, axis=None, limit=None):
        return self.fillna(method='bfill', limit=limit, axis=axis)

    def sample(self, frac, replace=False, random_state=None):
        """ Random sample of items

        Parameters
        ----------
        frac : float, optional
            Fraction of axis items to return.
        replace: boolean, optional
            Sample with or without replacement. Default = False.
        random_state: int or ``np.random.RandomState``
            If int we create a new RandomState with this as the seed
            Otherwise we draw from the passed RandomState

        See Also
        --------
        DataFrame.random_split
        pandas.DataFrame.sample
        """

        if random_state is None:
            random_state = np.random.RandomState()

        name = 'sample-' + tokenize(self, frac, replace, random_state)

        state_data = random_state_data(self.npartitions, random_state)
        dsk = {(name, i): (methods.sample, (self._name, i), state, frac, replace)
               for i, state in enumerate(state_data)}

        return new_dd_object(merge(self.dask, dsk), name,
                             self._meta, self.divisions)

    def to_hdf(self, path_or_buf, key, mode='a', append=False, get=None, **kwargs):
        """ See dd.to_hdf docstring for more information """
        from .io import to_hdf
        return to_hdf(self, path_or_buf, key, mode, append, get=get, **kwargs)

    def to_parquet(self, path, *args, **kwargs):
        """ See dd.to_parquet docstring for more information """
        from .io import to_parquet
        return to_parquet(self, path, *args, **kwargs)

    def to_csv(self, filename, **kwargs):
        """ See dd.to_csv docstring for more information """
        from .io import to_csv
        return to_csv(self, filename, **kwargs)

    def to_delayed(self, optimize_graph=True):
        """Convert into a list of ``dask.delayed`` objects, one per partition.

        Parameters
        ----------
        optimize_graph : bool, optional
            If True [default], the graph is optimized before converting into
            ``dask.delayed`` objects.

        Examples
        --------
        >>> partitions = df.to_delayed()  # doctest: +SKIP

        See Also
        --------
        dask.dataframe.from_delayed
        """
        from dask.delayed import Delayed
        keys = self.__dask_keys__()
        dsk = self.__dask_graph__()
        if optimize_graph:
            dsk = self.__dask_optimize__(dsk, keys)
        return [Delayed(k, dsk) for k in keys]

    @classmethod
    def _get_unary_operator(cls, op):
        return lambda self: elemwise(op, self)

    @classmethod
    def _get_binary_operator(cls, op, inv=False):
        if inv:
            return lambda self, other: elemwise(op, other, self)
        else:
            return lambda self, other: elemwise(op, self, other)

    def rolling(self, window, min_periods=None, freq=None, center=False,
                win_type=None, axis=0):
        """Provides rolling transformations.

        Parameters
        ----------
        window : int, str, offset
           Size of the moving window. This is the number of observations used
           for calculating the statistic. The window size must not be so large
           as to span more than one adjacent partition. If using an offset
           or offset alias like '5D', the data must have a ``DatetimeIndex``

           .. versionchanged:: 0.15.0

              Now accepts offsets and string offset aliases

        min_periods : int, default None
            Minimum number of observations in window required to have a value
            (otherwise result is NA).
        center : boolean, default False
            Set the labels at the center of the window.
        win_type : string, default None
            Provide a window type. The recognized window types are identical
            to pandas.
        axis : int, default 0

        Returns
        -------
        a Rolling object on which to call a method to compute a statistic

        Notes
        -----
        The `freq` argument is not supported.
        """
        from dask.dataframe.rolling import Rolling

        if isinstance(window, int):
            if window < 0:
                raise ValueError('window must be >= 0')

        if min_periods is not None:
            if not isinstance(min_periods, int):
                raise ValueError('min_periods must be an integer')
            if min_periods < 0:
                raise ValueError('min_periods must be >= 0')

        return Rolling(self, window=window, min_periods=min_periods,
                       freq=freq, center=center, win_type=win_type, axis=axis)

    @derived_from(pd.DataFrame)
    def diff(self, periods=1, axis=0):
        axis = self._validate_axis(axis)
        if not isinstance(periods, int):
            raise TypeError("periods must be an integer")

        if axis == 1:
            return self.map_partitions(M.diff, token='diff', periods=periods,
                                       axis=1)

        before, after = (periods, 0) if periods > 0 else (0, -periods)
        return self.map_overlap(M.diff, before, after, token='diff',
                                periods=periods)

    @derived_from(pd.DataFrame)
    def shift(self, periods=1, freq=None, axis=0):
        axis = self._validate_axis(axis)
        if not isinstance(periods, int):
            raise TypeError("periods must be an integer")

        if axis == 1:
            return self.map_partitions(M.shift, token='shift', periods=periods,
                                       freq=freq, axis=1)

        if freq is None:
            before, after = (periods, 0) if periods > 0 else (0, -periods)
            return self.map_overlap(M.shift, before, after, token='shift',
                                    periods=periods)

        # Let pandas error on invalid arguments
        meta = self._meta_nonempty.shift(periods, freq=freq)
        out = self.map_partitions(M.shift, token='shift', periods=periods,
                                  freq=freq, meta=meta)
        return maybe_shift_divisions(out, periods, freq=freq)

    def _reduction_agg(self, name, axis=None, skipna=True,
                       split_every=False, out=None):
        axis = self._validate_axis(axis)

        meta = getattr(self._meta_nonempty, name)(axis=axis, skipna=skipna)
        token = self._token_prefix + name

        method = getattr(M, name)
        if axis == 1:
            result = self.map_partitions(method, meta=meta,
                                         token=token, skipna=skipna, axis=axis)
            return handle_out(out, result)
        else:
            result = self.reduction(method, meta=meta, token=token,
                                    skipna=skipna, axis=axis,
                                    split_every=split_every)
            if isinstance(self, DataFrame):
                result.divisions = (min(self.columns), max(self.columns))
            return handle_out(out, result)

    @derived_from(pd.DataFrame)
    def abs(self):
        _raise_if_object_series(self, "abs")
        meta = self._meta_nonempty.abs()
        return self.map_partitions(M.abs, meta=meta)

    @derived_from(pd.DataFrame)
    def all(self, axis=None, skipna=True, split_every=False, out=None):
        return self._reduction_agg('all', axis=axis, skipna=skipna,
                                   split_every=split_every, out=out)

    @derived_from(pd.DataFrame)
    def any(self, axis=None, skipna=True, split_every=False, out=None):
        return self._reduction_agg('any', axis=axis, skipna=skipna,
                                   split_every=split_every, out=out)

    @derived_from(pd.DataFrame)
    def sum(self, axis=None, skipna=True, split_every=False, dtype=None, out=None):
        return self._reduction_agg('sum', axis=axis, skipna=skipna,
                                   split_every=split_every, out=out)

    @derived_from(pd.DataFrame)
    def prod(self, axis=None, skipna=True, split_every=False, dtype=None, out=None):
        return self._reduction_agg('prod', axis=axis, skipna=skipna,
                                   split_every=split_every, out=out)

    @derived_from(pd.DataFrame)
    def max(self, axis=None, skipna=True, split_every=False, out=None):
        return self._reduction_agg('max', axis=axis, skipna=skipna,
                                   split_every=split_every, out=out)

    @derived_from(pd.DataFrame)
    def min(self, axis=None, skipna=True, split_every=False, out=None):
        return self._reduction_agg('min', axis=axis, skipna=skipna,
                                   split_every=split_every, out=out)

    @derived_from(pd.DataFrame)
    def idxmax(self, axis=None, skipna=True, split_every=False):
        fn = 'idxmax'
        axis = self._validate_axis(axis)
        meta = self._meta_nonempty.idxmax(axis=axis, skipna=skipna)
        if axis == 1:
            return map_partitions(M.idxmax, self, meta=meta,
                                  token=self._token_prefix + fn,
                                  skipna=skipna, axis=axis)
        else:
            scalar = not isinstance(meta, pd.Series)
            result = aca([self], chunk=idxmaxmin_chunk, aggregate=idxmaxmin_agg,
                         combine=idxmaxmin_combine, meta=meta,
                         aggregate_kwargs={'scalar': scalar},
                         token=self._token_prefix + fn, split_every=split_every,
                         skipna=skipna, fn=fn)
            if isinstance(self, DataFrame):
                result.divisions = (min(self.columns), max(self.columns))
            return result

    @derived_from(pd.DataFrame)
    def idxmin(self, axis=None, skipna=True, split_every=False):
        fn = 'idxmin'
        axis = self._validate_axis(axis)
        meta = self._meta_nonempty.idxmax(axis=axis)
        if axis == 1:
            return map_partitions(M.idxmin, self, meta=meta,
                                  token=self._token_prefix + fn,
                                  skipna=skipna, axis=axis)
        else:
            scalar = not isinstance(meta, pd.Series)
            result = aca([self], chunk=idxmaxmin_chunk, aggregate=idxmaxmin_agg,
                         combine=idxmaxmin_combine, meta=meta,
                         aggregate_kwargs={'scalar': scalar},
                         token=self._token_prefix + fn, split_every=split_every,
                         skipna=skipna, fn=fn)
            if isinstance(self, DataFrame):
                result.divisions = (min(self.columns), max(self.columns))
            return result

    @derived_from(pd.DataFrame)
    def count(self, axis=None, split_every=False):
        axis = self._validate_axis(axis)
        token = self._token_prefix + 'count'
        if axis == 1:
            meta = self._meta_nonempty.count(axis=axis)
            return self.map_partitions(M.count, meta=meta, token=token,
                                       axis=axis)
        else:
            meta = self._meta_nonempty.count()
            result = self.reduction(M.count, aggregate=M.sum, meta=meta,
                                    token=token, split_every=split_every)
            if isinstance(self, DataFrame):
                result.divisions = (min(self.columns), max(self.columns))
            return result

    @derived_from(pd.DataFrame)
    def mean(self, axis=None, skipna=True, split_every=False, dtype=None, out=None):
        axis = self._validate_axis(axis)
        _raise_if_object_series(self, "mean")
        meta = self._meta_nonempty.mean(axis=axis, skipna=skipna)
        if axis == 1:
            result = map_partitions(M.mean, self, meta=meta,
                                    token=self._token_prefix + 'mean',
                                    axis=axis, skipna=skipna)
            return handle_out(out, result)
        else:
            num = self._get_numeric_data()
            s = num.sum(skipna=skipna, split_every=split_every)
            n = num.count(split_every=split_every)
            name = self._token_prefix + 'mean-%s' % tokenize(self, axis, skipna)
            result = map_partitions(methods.mean_aggregate, s, n,
                                    token=name, meta=meta)
            if isinstance(self, DataFrame):
                result.divisions = (min(self.columns), max(self.columns))
            return handle_out(out, result)

    @derived_from(pd.DataFrame)
    def var(self, axis=None, skipna=True, ddof=1, split_every=False, dtype=None, out=None):
        axis = self._validate_axis(axis)
        _raise_if_object_series(self, "var")
        meta = self._meta_nonempty.var(axis=axis, skipna=skipna)
        if axis == 1:
            result = map_partitions(M.var, self, meta=meta,
                                    token=self._token_prefix + 'var',
                                    axis=axis, skipna=skipna, ddof=ddof)
            return handle_out(out, result)
        else:
            num = self._get_numeric_data()
            x = 1.0 * num.sum(skipna=skipna, split_every=split_every)
            x2 = 1.0 * (num ** 2).sum(skipna=skipna, split_every=split_every)
            n = num.count(split_every=split_every)
            name = self._token_prefix + 'var'
            result = map_partitions(methods.var_aggregate, x2, x, n,
                                    token=name, meta=meta, ddof=ddof)
            if isinstance(self, DataFrame):
                result.divisions = (min(self.columns), max(self.columns))
            return handle_out(out, result)

    @derived_from(pd.DataFrame)
    def std(self, axis=None, skipna=True, ddof=1, split_every=False, dtype=None, out=None):
        axis = self._validate_axis(axis)
        _raise_if_object_series(self, "std")
        meta = self._meta_nonempty.std(axis=axis, skipna=skipna)
        if axis == 1:
            result = map_partitions(M.std, self, meta=meta,
                                    token=self._token_prefix + 'std',
                                    axis=axis, skipna=skipna, ddof=ddof)
            return handle_out(out, result)
        else:
            v = self.var(skipna=skipna, ddof=ddof, split_every=split_every)
            name = self._token_prefix + 'std'
            result = map_partitions(np.sqrt, v, meta=meta, token=name)
            return handle_out(out, result)

    @derived_from(pd.DataFrame)
    def sem(self, axis=None, skipna=None, ddof=1, split_every=False):
        axis = self._validate_axis(axis)
        _raise_if_object_series(self, "sem")
        meta = self._meta_nonempty.sem(axis=axis, skipna=skipna, ddof=ddof)
        if axis == 1:
            return map_partitions(M.sem, self, meta=meta,
                                  token=self._token_prefix + 'sem',
                                  axis=axis, skipna=skipna, ddof=ddof)
        else:
            num = self._get_numeric_data()
            v = num.var(skipna=skipna, ddof=ddof, split_every=split_every)
            n = num.count(split_every=split_every)
            name = self._token_prefix + 'sem'
            result = map_partitions(np.sqrt, v / n, meta=meta, token=name)
            if isinstance(self, DataFrame):
                result.divisions = (min(self.columns), max(self.columns))
            return result

    def quantile(self, q=0.5, axis=0):
        """ Approximate row-wise and precise column-wise quantiles of DataFrame

        Parameters
        ----------

        q : list/array of floats, default 0.5 (50%)
            Iterable of numbers ranging from 0 to 1 for the desired quantiles
        axis : {0, 1, 'index', 'columns'} (default 0)
            0 or 'index' for row-wise, 1 or 'columns' for column-wise
        """
        axis = self._validate_axis(axis)
        keyname = 'quantiles-concat--' + tokenize(self, q, axis)

        if axis == 1:
            if isinstance(q, list):
                # Not supported, the result will have current index as columns
                raise ValueError("'q' must be scalar when axis=1 is specified")
            return map_partitions(M.quantile, self, q, axis,
                                  token=keyname, meta=(q, 'f8'))
        else:
            _raise_if_object_series(self, "quantile")
            meta = self._meta.quantile(q, axis=axis)
            num = self._get_numeric_data()
            quantiles = tuple(quantile(self[c], q) for c in num.columns)

            dask = {}
            dask = merge(dask, *[_q.dask for _q in quantiles])
            qnames = [(_q._name, 0) for _q in quantiles]

            if isinstance(quantiles[0], Scalar):
                dask[(keyname, 0)] = (pd.Series, qnames, num.columns,
                                      None, meta.name)
                divisions = (min(num.columns), max(num.columns))
                return Series(dask, keyname, meta, divisions)
            else:
                dask[(keyname, 0)] = (methods.concat, qnames, 1)
                return DataFrame(dask, keyname, meta, quantiles[0].divisions)

    @derived_from(pd.DataFrame)
    def describe(self, split_every=False):
        # currently, only numeric describe is supported
        num = self._get_numeric_data()
        if self.ndim == 2 and len(num.columns) == 0:
            raise ValueError("DataFrame contains only non-numeric data.")
        elif self.ndim == 1 and self.dtype == 'object':
            raise ValueError("Cannot compute ``describe`` on object dtype.")

        stats = [num.count(split_every=split_every),
                 num.mean(split_every=split_every),
                 num.std(split_every=split_every),
                 num.min(split_every=split_every),
                 num.quantile([0.25, 0.5, 0.75]),
                 num.max(split_every=split_every)]
        stats_names = [(s._name, 0) for s in stats]

        name = 'describe--' + tokenize(self, split_every)
        dsk = merge(num.dask, *(s.dask for s in stats))
        dsk[(name, 0)] = (methods.describe_aggregate, stats_names)

        return new_dd_object(dsk, name, num._meta, divisions=[None, None])

    def _cum_agg(self, op_name, chunk, aggregate, axis, skipna=True,
                 chunk_kwargs=None, out=None):
        """ Wrapper for cumulative operation """

        axis = self._validate_axis(axis)

        if axis == 1:
            name = '{0}{1}(axis=1)'.format(self._token_prefix, op_name)
            result = self.map_partitions(chunk, token=name, **chunk_kwargs)
            return handle_out(out, result)
        else:
            # cumulate each partitions
            name1 = '{0}{1}-map'.format(self._token_prefix, op_name)
            cumpart = map_partitions(chunk, self, token=name1, meta=self,
                                     **chunk_kwargs)

            name2 = '{0}{1}-take-last'.format(self._token_prefix, op_name)
            cumlast = map_partitions(_take_last, cumpart, skipna,
                                     meta=pd.Series([]), token=name2)

            suffix = tokenize(self)
            name = '{0}{1}-{2}'.format(self._token_prefix, op_name, suffix)
            cname = '{0}{1}-cum-last-{2}'.format(self._token_prefix, op_name,
                                                 suffix)

            # aggregate cumulated partisions and its previous last element
            dask = {}
            dask[(name, 0)] = (cumpart._name, 0)

            for i in range(1, self.npartitions):
                # store each cumulative step to graph to reduce computation
                if i == 1:
                    dask[(cname, i)] = (cumlast._name, i - 1)
                else:
                    # aggregate with previous cumulation results
                    dask[(cname, i)] = (aggregate, (cname, i - 1),
                                        (cumlast._name, i - 1))
                dask[(name, i)] = (aggregate, (cumpart._name, i), (cname, i))
            result = new_dd_object(merge(dask, cumpart.dask, cumlast.dask),
                                   name, chunk(self._meta), self.divisions)
            return handle_out(out, result)

    @derived_from(pd.DataFrame)
    def cumsum(self, axis=None, skipna=True, dtype=None, out=None):
        return self._cum_agg('cumsum',
                             chunk=M.cumsum,
                             aggregate=operator.add,
                             axis=axis, skipna=skipna,
                             chunk_kwargs=dict(axis=axis, skipna=skipna),
                             out=out)

    @derived_from(pd.DataFrame)
    def cumprod(self, axis=None, skipna=True, dtype=None, out=None):
        return self._cum_agg('cumprod',
                             chunk=M.cumprod,
                             aggregate=operator.mul,
                             axis=axis, skipna=skipna,
                             chunk_kwargs=dict(axis=axis, skipna=skipna),
                             out=out)

    @derived_from(pd.DataFrame)
    def cummax(self, axis=None, skipna=True, out=None):
        return self._cum_agg('cummax',
                             chunk=M.cummax,
                             aggregate=methods.cummax_aggregate,
                             axis=axis, skipna=skipna,
                             chunk_kwargs=dict(axis=axis, skipna=skipna),
                             out=out)

    @derived_from(pd.DataFrame)
    def cummin(self, axis=None, skipna=True, out=None):
        return self._cum_agg('cummin',
                             chunk=M.cummin,
                             aggregate=methods.cummin_aggregate,
                             axis=axis, skipna=skipna,
                             chunk_kwargs=dict(axis=axis, skipna=skipna),
                             out=out)

    @derived_from(pd.DataFrame)
    def where(self, cond, other=np.nan):
        # cond and other may be dask instance,
        # passing map_partitions via keyword will not be aligned
        return map_partitions(M.where, self, cond, other)

    @derived_from(pd.DataFrame)
    def mask(self, cond, other=np.nan):
        return map_partitions(M.mask, self, cond, other)

    @derived_from(pd.DataFrame)
    def notnull(self):
        return self.map_partitions(M.notnull)

    @derived_from(pd.DataFrame)
    def isnull(self):
        return self.map_partitions(M.isnull)

    @derived_from(pd.DataFrame)
    def isna(self):
        if hasattr(pd, 'isna'):
            return self.map_partitions(M.isna)
        else:
            raise NotImplementedError("Need more recent version of Pandas "
                                      "to support isna. "
                                      "Please use isnull instead.")

    @derived_from(pd.DataFrame)
    def isin(self, values):
        return elemwise(M.isin, self, list(values))

    @derived_from(pd.DataFrame)
    def astype(self, dtype):
        # XXX: Pandas will segfault for empty dataframes when setting
        # categorical dtypes. This operation isn't allowed currently anyway. We
        # get the metadata with a non-empty frame to throw the error instead of
        # segfaulting.
        if isinstance(self._meta, pd.DataFrame) and is_categorical_dtype(dtype):
            meta = self._meta_nonempty.astype(dtype)
        else:
            meta = self._meta.astype(dtype)
        if hasattr(dtype, 'items'):
            # Pandas < 0.21.0, no `categories` attribute, so unknown
            # Pandas >= 0.21.0, known if `categories` attribute is not None
            set_unknown = [k for k, v in dtype.items()
                           if (is_categorical_dtype(v) and
                               getattr(v, 'categories', None) is None)]
            meta = clear_known_categories(meta, cols=set_unknown)
        elif (is_categorical_dtype(dtype) and
              getattr(dtype, 'categories', None) is None):
            meta = clear_known_categories(meta)
        return self.map_partitions(M.astype, dtype=dtype, meta=meta)

    @derived_from(pd.Series)
    def append(self, other):
        # because DataFrame.append will override the method,
        # wrap by pd.Series.append docstring
        from .multi import concat

        if isinstance(other, (list, dict)):
            msg = "append doesn't support list or dict input"
            raise NotImplementedError(msg)

        return concat([self, other], join='outer', interleave_partitions=False)

    @derived_from(pd.DataFrame)
    def align(self, other, join='outer', axis=None, fill_value=None):
        meta1, meta2 = _emulate(M.align, self, other, join, axis=axis,
                                fill_value=fill_value)
        aligned = self.map_partitions(M.align, other, join=join, axis=axis,
                                      fill_value=fill_value)

        token = tokenize(self, other, join, axis, fill_value)

        name1 = 'align1-' + token
        dsk1 = {(name1, i): (getitem, key, 0)
                for i, key in enumerate(aligned.__dask_keys__())}
        dsk1.update(aligned.dask)
        result1 = new_dd_object(dsk1, name1, meta1, aligned.divisions)

        name2 = 'align2-' + token
        dsk2 = {(name2, i): (getitem, key, 1)
                for i, key in enumerate(aligned.__dask_keys__())}
        dsk2.update(aligned.dask)
        result2 = new_dd_object(dsk2, name2, meta2, aligned.divisions)

        return result1, result2

    @derived_from(pd.DataFrame)
    def combine(self, other, func, fill_value=None, overwrite=True):
        return self.map_partitions(M.combine, other, func,
                                   fill_value=fill_value, overwrite=overwrite)

    @derived_from(pd.DataFrame)
    def combine_first(self, other):
        return self.map_partitions(M.combine_first, other)

    @classmethod
    def _bind_operator_method(cls, name, op):
        """ bind operator method like DataFrame.add to this class """
        raise NotImplementedError

    @derived_from(pd.DataFrame)
    def resample(self, rule, how=None, closed=None, label=None):
        from .tseries.resample import _resample
        return _resample(self, rule, how=how, closed=closed, label=label)

    @derived_from(pd.DataFrame)
    def first(self, offset):
        # Let pandas error on bad args
        self._meta_nonempty.first(offset)

        if not self.known_divisions:
            raise ValueError("`first` is not implemented for unknown divisions")

        offset = pd.tseries.frequencies.to_offset(offset)
        date = self.divisions[0] + offset
        end = self.loc._get_partitions(date)

        include_right = offset.isAnchored() or not hasattr(offset, '_inc')

        if end == self.npartitions - 1:
            divs = self.divisions
        else:
            divs = self.divisions[:end + 1] + (date,)

        name = 'first-' + tokenize(self, offset)
        dsk = {(name, i): (self._name, i) for i in range(end)}
        dsk[(name, end)] = (methods.boundary_slice, (self._name, end),
                            None, date, include_right, True, 'ix')
        return new_dd_object(merge(self.dask, dsk), name, self, divs)

    @derived_from(pd.DataFrame)
    def last(self, offset):
        # Let pandas error on bad args
        self._meta_nonempty.first(offset)

        if not self.known_divisions:
            raise ValueError("`last` is not implemented for unknown divisions")

        offset = pd.tseries.frequencies.to_offset(offset)
        date = self.divisions[-1] - offset
        start = self.loc._get_partitions(date)

        if start == 0:
            divs = self.divisions
        else:
            divs = (date,) + self.divisions[start + 1:]

        name = 'last-' + tokenize(self, offset)
        dsk = {(name, i + 1): (self._name, j + 1)
               for i, j in enumerate(range(start, self.npartitions))}
        dsk[(name, 0)] = (methods.boundary_slice, (self._name, start),
                          date, None, True, False, 'ix')
        return new_dd_object(merge(self.dask, dsk), name, self, divs)

    def nunique_approx(self, split_every=None):
        """Approximate number of unique rows.

        This method uses the HyperLogLog algorithm for cardinality
        estimation to compute the approximate number of unique rows.
        The approximate error is 0.406%.

        Parameters
        ----------
        split_every : int, optional
            Group partitions into groups of this size while performing a
            tree-reduction. If set to False, no tree-reduction will be used.
            Default is 8.

        Returns
        -------
        a float representing the approximate number of elements
        """
        from . import hyperloglog # here to avoid circular import issues

        return aca([self], chunk=hyperloglog.compute_hll_array,
                   combine=hyperloglog.reduce_state,
                   aggregate=hyperloglog.estimate_count,
                   split_every=split_every, b=16, meta=float)

    @property
    def values(self):
        """ Return a dask.array of the values of this dataframe

        Warning: This creates a dask.array without precise shape information.
        Operations that depend on shape information, like slicing or reshaping,
        will not work.
        """
        return self.map_partitions(methods.values)


def _raise_if_object_series(x, funcname):
    """
    Utility function to raise an error if an object column does not support
    a certain operation like `mean`.
    """
    if isinstance(x, Series) and hasattr(x, "dtype") and x.dtype == object:
        raise ValueError("`%s` not supported with object series" % funcname)


class Series(_Frame):
    """ Parallel Pandas Series

    Do not use this class directly.  Instead use functions like
    ``dd.read_csv``, ``dd.read_parquet``, or ``dd.from_pandas``.

    Parameters
    ----------

    dsk: dict
        The dask graph to compute this Series
    _name: str
        The key prefix that specifies which keys in the dask comprise this
        particular Series
    meta: pandas.Series
        An empty ``pandas.Series`` with names, dtypes, and index matching the
        expected output.
    divisions: tuple of index values
        Values along which we partition our blocks on the index

    See Also
    --------
    dask.dataframe.DataFrame
    """

    _partition_type = pd.Series
    _token_prefix = 'series-'

    def __array_wrap__(self, array, context=None):
        if isinstance(context, tuple) and len(context) > 0:
            if isinstance(context[1][0], np.ndarray) and context[1][0].shape == ():
                index = None
            else:
                index = context[1][0].index

        return pd.Series(array, index=index, name=self.name)

    @property
    def name(self):
        return self._meta.name

    @name.setter
    def name(self, name):
        self._meta.name = name
        renamed = _rename_dask(self, name)
        # update myself
        self.dask.update(renamed.dask)
        self._name = renamed._name

    @property
    def ndim(self):
        """ Return dimensionality """
        return 1

    @property
    def dtype(self):
        """ Return data type """
        return self._meta.dtype

    @cache_readonly
    def dt(self):
        """ Namespace of datetime methods """
        return DatetimeAccessor(self)

    @cache_readonly
    def cat(self):
        return CategoricalAccessor(self)

    @cache_readonly
    def str(self):
        """ Namespace for string methods """
        return StringAccessor(self)

    def __dir__(self):
        o = set(dir(type(self)))
        o.update(self.__dict__)
        # Remove the `cat` and `str` accessors if not available. We can't
        # decide this statically for the `dt` accessor, as it works on
        # datetime-like things as well.
        for accessor in ['cat', 'str']:
            if not hasattr(self._meta, accessor):
                o.remove(accessor)
        return list(o)

    @property
    def nbytes(self):
        """ Number of bytes """
        return self.reduction(methods.nbytes, np.sum, token='nbytes',
                              meta=int, split_every=False)

    @property
    def _repr_data(self):
        return _repr_data_series(self._meta, self._repr_divisions)

    def __repr__(self):
        """ have to overwrite footer """
        if self.name is not None:
            footer = "Name: {name}, dtype: {dtype}".format(name=self.name,
                                                           dtype=self.dtype)
        else:
            footer = "dtype: {dtype}".format(dtype=self.dtype)

        return """Dask {klass} Structure:
{data}
{footer}
Dask Name: {name}, {task} tasks""".format(klass=self.__class__.__name__,
                                          data=self.to_string(),
                                          footer=footer,
                                          name=key_split(self._name),
                                          task=len(self.dask))

    def rename(self, index=None, inplace=False, sorted_index=False):
        """Alter Series index labels or name

        Function / dict values must be unique (1-to-1). Labels not contained in
        a dict / Series will be left as-is. Extra labels listed don't throw an
        error.

        Alternatively, change ``Series.name`` with a scalar value.

        Parameters
        ----------
        index : scalar, hashable sequence, dict-like or callable, optional
            If dict-like or callable, the transformation is applied to the
            index. Scalar or hashable sequence-like will alter the
            ``Series.name`` attribute.
        inplace : boolean, default False
            Whether to return a new Series or modify this one inplace.
        sorted_index : bool, default False
            If true, the output ``Series`` will have known divisions inferred
            from the input series and the transformation. Ignored for
            non-callable/dict-like ``index`` or when the input series has
            unknown divisions. Note that this may only be set to ``True`` if
            you know that the transformed index is monotonicly increasing. Dask
            will check that transformed divisions are monotonic, but cannot
            check all the values between divisions, so incorrectly setting this
            can result in bugs.

        Returns
        -------
        renamed : Series

        See Also
        --------
        pandas.Series.rename
        """
        from pandas.api.types import is_scalar, is_list_like, is_dict_like
        if is_scalar(index) or (is_list_like(index) and not is_dict_like(index)):
            res = self if inplace else self.copy()
            res.name = index
        else:
            res = self.map_partitions(M.rename, index)
            if self.known_divisions:
                if sorted_index and (callable(index) or is_dict_like(index)):
                    old = pd.Series(range(self.npartitions + 1),
                                    index=self.divisions)
                    new = old.rename(index).index
                    if not new.is_monotonic_increasing:
                        msg = ("sorted_index=True, but the transformed index "
                               "isn't monotonic_increasing")
                        raise ValueError(msg)
                    res.divisions = tuple(new.tolist())
                else:
                    res = res.clear_divisions()
            if inplace:
                self.dask = res.dask
                self._name = res._name
                self.divisions = res.divisions
                self._meta = res._meta
                res = self
        return res

    @derived_from(pd.Series)
    def round(self, decimals=0):
        return elemwise(M.round, self, decimals)

    @derived_from(pd.DataFrame)
    def to_timestamp(self, freq=None, how='start', axis=0):
        df = elemwise(M.to_timestamp, self, freq, how, axis)
        df.divisions = tuple(pd.Index(self.divisions).to_timestamp())
        return df

    def quantile(self, q=0.5):
        """ Approximate quantiles of Series

        q : list/array of floats, default 0.5 (50%)
            Iterable of numbers ranging from 0 to 1 for the desired quantiles
        """
        return quantile(self, q)

    def _repartition_quantiles(self, npartitions, upsample=1.0):
        """ Approximate quantiles of Series used for repartitioning
        """
        from .partitionquantiles import partition_quantiles
        return partition_quantiles(self, npartitions, upsample=upsample)

    def __getitem__(self, key):
        if isinstance(key, Series) and self.divisions == key.divisions:
            name = 'index-%s' % tokenize(self, key)
            dsk = dict(((name, i), (operator.getitem, (self._name, i),
                                    (key._name, i)))
                       for i in range(self.npartitions))
            return Series(merge(self.dask, key.dask, dsk), name,
                          self._meta, self.divisions)
        raise NotImplementedError()

    @derived_from(pd.DataFrame)
    def _get_numeric_data(self, how='any', subset=None):
        return self

    @derived_from(pd.Series)
    def iteritems(self):
        for i in range(self.npartitions):
            s = self.get_partition(i).compute()
            for item in s.iteritems():
                yield item

    @classmethod
    def _validate_axis(cls, axis=0):
        if axis not in (0, 'index', None):
            raise ValueError('No axis named {0}'.format(axis))
        # convert to numeric axis
        return {None: 0, 'index': 0}.get(axis, axis)

    @derived_from(pd.Series)
    def groupby(self, by=None, **kwargs):
        from dask.dataframe.groupby import SeriesGroupBy
        return SeriesGroupBy(self, by=by, **kwargs)

    @derived_from(pd.Series)
    def count(self, split_every=False):
        return super(Series, self).count(split_every=split_every)

    def unique(self, split_every=None, split_out=1):
        """
        Return Series of unique values in the object. Includes NA values.

        Returns
        -------
        uniques : Series
        """
        return aca(self, chunk=methods.unique, aggregate=methods.unique,
                   meta=self._meta, token='unique', split_every=split_every,
                   series_name=self.name, split_out=split_out)

    @derived_from(pd.Series)
    def nunique(self, split_every=None):
        return self.drop_duplicates(split_every=split_every).count()

    @derived_from(pd.Series)
    def value_counts(self, split_every=None, split_out=1):
        return aca(self, chunk=M.value_counts,
                   aggregate=methods.value_counts_aggregate,
                   combine=methods.value_counts_combine,
                   meta=self._meta.value_counts(), token='value-counts',
                   split_every=split_every, split_out=split_out,
                   split_out_setup=split_out_on_index)

    @derived_from(pd.Series)
    def nlargest(self, n=5, split_every=None):
        return aca(self, chunk=M.nlargest, aggregate=M.nlargest,
                   meta=self._meta, token='series-nlargest',
                   split_every=split_every, n=n)

    @derived_from(pd.Series)
    def nsmallest(self, n=5, split_every=None):
        return aca(self, chunk=M.nsmallest, aggregate=M.nsmallest,
                   meta=self._meta, token='series-nsmallest',
                   split_every=split_every, n=n)

    @derived_from(pd.Series)
    def isin(self, values):
        return elemwise(M.isin, self, list(values))

    @insert_meta_param_description(pad=12)
    @derived_from(pd.Series)
    def map(self, arg, na_action=None, meta=no_default):
        if not (isinstance(arg, (pd.Series, dict)) or callable(arg)):
            raise TypeError("arg must be pandas.Series, dict or callable."
                            " Got {0}".format(type(arg)))
        name = 'map-' + tokenize(self, arg, na_action)
        dsk = {(name, i): (M.map, k, arg, na_action) for i, k in
               enumerate(self.__dask_keys__())}
        dsk.update(self.dask)
        if meta is no_default:
            meta = _emulate(M.map, self, arg, na_action=na_action, udf=True)
        else:
            meta = make_meta(meta)

        return Series(dsk, name, meta, self.divisions)

    @derived_from(pd.Series)
    def dropna(self):
        return self.map_partitions(M.dropna)

    @derived_from(pd.Series)
    def between(self, left, right, inclusive=True):
        return self.map_partitions(M.between, left=left,
                                   right=right, inclusive=inclusive)

    @derived_from(pd.Series)
    def clip(self, lower=None, upper=None, out=None):
        if out is not None:
            raise ValueError("'out' must be None")
        # np.clip may pass out
        return self.map_partitions(M.clip, lower=lower, upper=upper)

    @derived_from(pd.Series)
    def clip_lower(self, threshold):
        return self.map_partitions(M.clip_lower, threshold=threshold)

    @derived_from(pd.Series)
    def clip_upper(self, threshold):
        return self.map_partitions(M.clip_upper, threshold=threshold)

    @derived_from(pd.Series)
    def align(self, other, join='outer', axis=None, fill_value=None):
        return super(Series, self).align(other, join=join, axis=axis,
                                         fill_value=fill_value)

    @derived_from(pd.Series)
    def combine(self, other, func, fill_value=None):
        return self.map_partitions(M.combine, other, func,
                                   fill_value=fill_value)

    @derived_from(pd.Series)
    def squeeze(self):
        return self

    @derived_from(pd.Series)
    def combine_first(self, other):
        return self.map_partitions(M.combine_first, other)

    def to_bag(self, index=False):
        """ Craeate a Dask Bag from a Series """
        from .io import to_bag
        return to_bag(self, index)

    @derived_from(pd.Series)
    def to_frame(self, name=None):
        return self.map_partitions(M.to_frame, name,
                                   meta=self._meta.to_frame(name))

    @derived_from(pd.Series)
    def to_string(self, max_rows=5):
        # option_context doesn't affect
        return self._repr_data.to_string(max_rows=max_rows)

    @classmethod
    def _bind_operator_method(cls, name, op):
        """ bind operator method like DataFrame.add to this class """

        def meth(self, other, level=None, fill_value=None, axis=0):
            if level is not None:
                raise NotImplementedError('level must be None')
            axis = self._validate_axis(axis)
            meta = _emulate(op, self, other, axis=axis, fill_value=fill_value)
            return map_partitions(op, self, other, meta=meta,
                                  axis=axis, fill_value=fill_value)
        meth.__doc__ = op.__doc__
        bind_method(cls, name, meth)

    @classmethod
    def _bind_comparison_method(cls, name, comparison):
        """ bind comparison method like DataFrame.add to this class """

        def meth(self, other, level=None, axis=0):
            if level is not None:
                raise NotImplementedError('level must be None')
            axis = self._validate_axis(axis)
            return elemwise(comparison, self, other, axis=axis)

        meth.__doc__ = comparison.__doc__
        bind_method(cls, name, meth)

    @insert_meta_param_description(pad=12)
    def apply(self, func, convert_dtype=True, meta=no_default, args=(), **kwds):
        """ Parallel version of pandas.Series.apply

        Parameters
        ----------
        func : function
            Function to apply
        convert_dtype : boolean, default True
            Try to find better dtype for elementwise function results.
            If False, leave as dtype=object.
        $META
        args : tuple
            Positional arguments to pass to function in addition to the value.

        Additional keyword arguments will be passed as keywords to the function.

        Returns
        -------
        applied : Series or DataFrame if func returns a Series.

        Examples
        --------
        >>> import dask.dataframe as dd
        >>> s = pd.Series(range(5), name='x')
        >>> ds = dd.from_pandas(s, npartitions=2)

        Apply a function elementwise across the Series, passing in extra
        arguments in ``args`` and ``kwargs``:

        >>> def myadd(x, a, b=1):
        ...     return x + a + b
        >>> res = ds.apply(myadd, args=(2,), b=1.5)

        By default, dask tries to infer the output metadata by running your
        provided function on some fake data. This works well in many cases, but
        can sometimes be expensive, or even fail. To avoid this, you can
        manually specify the output metadata with the ``meta`` keyword. This
        can be specified in many forms, for more information see
        ``dask.dataframe.utils.make_meta``.

        Here we specify the output is a Series with name ``'x'``, and dtype
        ``float64``:

        >>> res = ds.apply(myadd, args=(2,), b=1.5, meta=('x', 'f8'))

        In the case where the metadata doesn't change, you can also pass in
        the object itself directly:

        >>> res = ds.apply(lambda x: x + 1, meta=ds)

        See Also
        --------
        dask.Series.map_partitions
        """
        if meta is no_default:
            msg = ("`meta` is not specified, inferred from partial data. "
                   "Please provide `meta` if the result is unexpected.\n"
                   "  Before: .apply(func)\n"
                   "  After:  .apply(func, meta={'x': 'f8', 'y': 'f8'}) for dataframe result\n"
                   "  or:     .apply(func, meta=('x', 'f8'))            for series result")
            warnings.warn(msg)

            meta = _emulate(M.apply, self._meta_nonempty, func,
                            convert_dtype=convert_dtype,
                            args=args, udf=True, **kwds)

        return map_partitions(M.apply, self, func,
                              convert_dtype, args, meta=meta, **kwds)

    @derived_from(pd.Series)
    def cov(self, other, min_periods=None, split_every=False):
        from .multi import concat
        if not isinstance(other, Series):
            raise TypeError("other must be a dask.dataframe.Series")
        df = concat([self, other], axis=1)
        return cov_corr(df, min_periods, scalar=True, split_every=split_every)

    @derived_from(pd.Series)
    def corr(self, other, method='pearson', min_periods=None,
             split_every=False):
        from .multi import concat
        if not isinstance(other, Series):
            raise TypeError("other must be a dask.dataframe.Series")
        if method != 'pearson':
            raise NotImplementedError("Only Pearson correlation has been "
                                      "implemented")
        df = concat([self, other], axis=1)
        return cov_corr(df, min_periods, corr=True, scalar=True,
                        split_every=split_every)

    @derived_from(pd.Series)
    def autocorr(self, lag=1, split_every=False):
        if not isinstance(lag, int):
            raise TypeError("lag must be an integer")
        return self.corr(self if lag == 0 else self.shift(lag),
                         split_every=split_every)

    @derived_from(pd.Series)
    def memory_usage(self, index=True, deep=False):
        from ..delayed import delayed
        result = self.map_partitions(M.memory_usage, index=index, deep=deep)
        return delayed(sum)(result.to_delayed())


class Index(Series):

    _partition_type = pd.Index
    _token_prefix = 'index-'

    _dt_attributes = {'nanosecond', 'microsecond', 'millisecond', 'dayofyear',
                      'minute', 'hour', 'day', 'dayofweek', 'second', 'week',
                      'weekday', 'weekofyear', 'month', 'quarter', 'year'}

    _cat_attributes = {'known', 'as_known', 'as_unknown', 'add_categories',
                       'categories', 'remove_categories', 'reorder_categories',
                       'as_ordered', 'codes', 'remove_unused_categories',
                       'set_categories', 'as_unordered', 'ordered',
                       'rename_categories'}

    def __getattr__(self, key):
        if is_categorical_dtype(self.dtype) and key in self._cat_attributes:
            return getattr(self.cat, key)
        elif key in self._dt_attributes:
            return getattr(self.dt, key)
        raise AttributeError("'Index' object has no attribute %r" % key)

    def __dir__(self):
        out = super(Index, self).__dir__()
        out.extend(self._dt_attributes)
        if is_categorical_dtype(self.dtype):
            out.extend(self._cat_attributes)
        return out

    @property
    def index(self):
        msg = "'{0}' object has no attribute 'index'"
        raise AttributeError(msg.format(self.__class__.__name__))

    def __array_wrap__(self, array, context=None):
        return pd.Index(array, name=self.name)

    def head(self, n=5, compute=True):
        """ First n items of the Index.

        Caveat, this only checks the first partition.
        """
        name = 'head-%d-%s' % (n, self._name)
        dsk = {(name, 0): (operator.getitem, (self._name, 0), slice(0, n))}

        result = new_dd_object(merge(self.dask, dsk), name,
                               self._meta, self.divisions[:2])

        if compute:
            result = result.compute()
        return result

    @derived_from(pd.Index)
    def max(self, split_every=False):
        return self.reduction(M.max, meta=self._meta_nonempty.max(),
                              token=self._token_prefix + 'max',
                              split_every=split_every)

    @derived_from(pd.Index)
    def min(self, split_every=False):
        return self.reduction(M.min, meta=self._meta_nonempty.min(),
                              token=self._token_prefix + 'min',
                              split_every=split_every)

    def count(self, split_every=False):
        return self.reduction(methods.index_count, np.sum,
                              token='index-count', meta=int,
                              split_every=split_every)

    @derived_from(pd.Index)
    def shift(self, periods=1, freq=None):
        if isinstance(self._meta, pd.PeriodIndex):
            if freq is not None:
                raise ValueError("PeriodIndex doesn't accept `freq` argument")
            meta = self._meta_nonempty.shift(periods)
            out = self.map_partitions(M.shift, periods, meta=meta,
                                      token='shift')
        else:
            # Pandas will raise for other index types that don't implement shift
            meta = self._meta_nonempty.shift(periods, freq=freq)
            out = self.map_partitions(M.shift, periods, token='shift',
                                      meta=meta, freq=freq)
        if freq is None:
            freq = meta.freq
        return maybe_shift_divisions(out, periods, freq=freq)


class DataFrame(_Frame):
    """
    Parallel Pandas DataFrame

    Do not use this class directly.  Instead use functions like
    ``dd.read_csv``, ``dd.read_parquet``, or ``dd.from_pandas``.

    Parameters
    ----------
    dsk: dict
        The dask graph to compute this DataFrame
    name: str
        The key prefix that specifies which keys in the dask comprise this
        particular DataFrame
    meta: pandas.DataFrame
        An empty ``pandas.DataFrame`` with names, dtypes, and index matching
        the expected output.
    divisions: tuple of index values
        Values along which we partition our blocks on the index
    """

    _partition_type = pd.DataFrame
    _token_prefix = 'dataframe-'

    def __array_wrap__(self, array, context=None):
        if isinstance(context, tuple) and len(context) > 0:
            if isinstance(context[1][0], np.ndarray) and context[1][0].shape == ():
                index = None
            else:
                index = context[1][0].index

        return pd.DataFrame(array, index=index, columns=self.columns)

    @property
    def columns(self):
        return self._meta.columns

    @columns.setter
    def columns(self, columns):
        renamed = _rename_dask(self, columns)
        self._meta = renamed._meta
        self._name = renamed._name
        self.dask.update(renamed.dask)

    def __getitem__(self, key):
        name = 'getitem-%s' % tokenize(self, key)
        if np.isscalar(key) or isinstance(key, (tuple, string_types)):

            if isinstance(self._meta.index, (pd.DatetimeIndex, pd.PeriodIndex)):
                if key not in self._meta.columns:
                    return self.loc[key]

            # error is raised from pandas
            meta = self._meta[_extract_meta(key)]
            dsk = dict(((name, i), (operator.getitem, (self._name, i), key))
                       for i in range(self.npartitions))
            return new_dd_object(merge(self.dask, dsk), name,
                                 meta, self.divisions)
        elif isinstance(key, slice):
            return self.loc[key]

        if isinstance(key, list):
            # error is raised from pandas
            meta = self._meta[_extract_meta(key)]

            dsk = dict(((name, i), (operator.getitem, (self._name, i), key))
                       for i in range(self.npartitions))
            return new_dd_object(merge(self.dask, dsk), name,
                                 meta, self.divisions)
        if isinstance(key, Series):
            # do not perform dummy calculation, as columns will not be changed.
            #
            if self.divisions != key.divisions:
                from .multi import _maybe_align_partitions
                self, key = _maybe_align_partitions([self, key])
            dsk = {(name, i): (M._getitem_array, (self._name, i), (key._name, i))
                   for i in range(self.npartitions)}
            return new_dd_object(merge(self.dask, key.dask, dsk), name,
                                 self, self.divisions)
        raise NotImplementedError(key)

    def __setitem__(self, key, value):
        if isinstance(key, (tuple, list)):
            df = self.assign(**{k: value[c]
                                for k, c in zip(key, value.columns)})
        else:
            df = self.assign(**{key: value})

        self.dask = df.dask
        self._name = df._name
        self._meta = df._meta
        self.divisions = df.divisions

    def __delitem__(self, key):
        result = self.drop([key], axis=1)
        self.dask = result.dask
        self._name = result._name
        self._meta = result._meta

    def __setattr__(self, key, value):
        try:
            columns = object.__getattribute__(self, '_meta').columns
        except AttributeError:
            columns = ()

        if key in columns:
            self[key] = value
        else:
            object.__setattr__(self, key, value)

    def __getattr__(self, key):
        if key in self.columns:
            meta = self._meta[key]
            name = 'getitem-%s' % tokenize(self, key)
            dsk = dict(((name, i), (operator.getitem, (self._name, i), key))
                       for i in range(self.npartitions))
            return new_dd_object(merge(self.dask, dsk), name,
                                 meta, self.divisions)
        raise AttributeError("'DataFrame' object has no attribute %r" % key)

    def __dir__(self):
        o = set(dir(type(self)))
        o.update(self.__dict__)
        o.update(c for c in self.columns if
                 (isinstance(c, pd.compat.string_types) and
                  pd.compat.isidentifier(c)))
        return list(o)

    @property
    def ndim(self):
        """ Return dimensionality """
        return 2

    @property
    def dtypes(self):
        """ Return data types """
        return self._meta.dtypes

    @derived_from(pd.DataFrame)
    def get_dtype_counts(self):
        return self._meta.get_dtype_counts()

    @derived_from(pd.DataFrame)
    def get_ftype_counts(self):
        return self._meta.get_ftype_counts()

    @derived_from(pd.DataFrame)
    def select_dtypes(self, include=None, exclude=None):
        cs = self._meta.select_dtypes(include=include, exclude=exclude).columns
        return self[list(cs)]

    def set_index(self, other, drop=True, sorted=False, npartitions=None,
                  divisions=None, **kwargs):
        """Set the DataFrame index (row labels) using an existing column

        This realigns the dataset to be sorted by a new column.  This can have a
        significant impact on performance, because joins, groupbys, lookups, etc.
        are all much faster on that column.  However, this performance increase
        comes with a cost, sorting a parallel dataset requires expensive shuffles.
        Often we ``set_index`` once directly after data ingest and filtering and
        then perform many cheap computations off of the sorted dataset.

        This function operates exactly like ``pandas.set_index`` except with
        different performance costs (it is much more expensive).  Under normal
        operation this function does an initial pass over the index column to
        compute approximate qunatiles to serve as future divisions.  It then passes
        over the data a second time, splitting up each input partition into several
        pieces and sharing those pieces to all of the output partitions now in
        sorted order.

        In some cases we can alleviate those costs, for example if your dataset is
        sorted already then we can avoid making many small pieces or if you know
        good values to split the new index column then we can avoid the initial
        pass over the data.  For example if your new index is a datetime index and
        your data is already sorted by day then this entire operation can be done
        for free.  You can control these options with the following parameters.

        Parameters
        ----------
        df: Dask DataFrame
        index: string or Dask Series
        npartitions: int, None, or 'auto'
            The ideal number of output partitions.   If None use the same as
            the input.  If 'auto' then decide by memory use.
        shuffle: string, optional
            Either ``'disk'`` for single-node operation or ``'tasks'`` for
            distributed operation.  Will be inferred by your current scheduler.
        sorted: bool, optional
            If the index column is already sorted in increasing order.
            Defaults to False
        divisions: list, optional
            Known values on which to separate index values of the partitions.
            See http://dask.pydata.org/en/latest/dataframe-design.html#partitions
            Defaults to computing this with a single pass over the data. Note
            that if ``sorted=True``, specified divisions are assumed to match
            the existing partitions in the data. If this is untrue, you should
            leave divisions empty and call ``repartition`` after ``set_index``.
        compute: bool
            Whether or not to trigger an immediate computation. Defaults to False.

        Examples
        --------
        >>> df2 = df.set_index('x')  # doctest: +SKIP
        >>> df2 = df.set_index(d.x)  # doctest: +SKIP
        >>> df2 = df.set_index(d.timestamp, sorted=True)  # doctest: +SKIP

        A common case is when we have a datetime column that we know to be
        sorted and is cleanly divided by day.  We can set this index for free
        by specifying both that the column is pre-sorted and the particular
        divisions along which is is separated

        >>> import pandas as pd
        >>> divisions = pd.date_range('2000', '2010', freq='1D')
        >>> df2 = df.set_index('timestamp', sorted=True, divisions=divisions)  # doctest: +SKIP
        """
        pre_sorted = sorted
        del sorted

        if divisions is not None:
            check_divisions(divisions)

        if pre_sorted:
            from .shuffle import set_sorted_index
            return set_sorted_index(self, other, drop=drop, divisions=divisions,
                                    **kwargs)
        else:
            from .shuffle import set_index
            return set_index(self, other, drop=drop, npartitions=npartitions,
                             divisions=divisions, **kwargs)

    @derived_from(pd.DataFrame)
    def nlargest(self, n=5, columns=None, split_every=None):
        token = 'dataframe-nlargest'
        return aca(self, chunk=M.nlargest, aggregate=M.nlargest,
                   meta=self._meta, token=token, split_every=split_every,
                   n=n, columns=columns)

    @derived_from(pd.DataFrame)
    def nsmallest(self, n=5, columns=None, split_every=None):
        token = 'dataframe-nsmallest'
        return aca(self, chunk=M.nsmallest, aggregate=M.nsmallest,
                   meta=self._meta, token=token, split_every=split_every,
                   n=n, columns=columns)

    @derived_from(pd.DataFrame)
    def groupby(self, by=None, **kwargs):
        from dask.dataframe.groupby import DataFrameGroupBy
        return DataFrameGroupBy(self, by=by, **kwargs)

    @wraps(categorize)
    def categorize(self, columns=None, index=None, split_every=None, **kwargs):
        return categorize(self, columns=columns, index=index,
                          split_every=split_every, **kwargs)

    @derived_from(pd.DataFrame)
    def assign(self, **kwargs):
        for k, v in kwargs.items():
            if not (isinstance(v, (Series, Scalar, pd.Series)) or
                    callable(v) or pd.api.types.is_scalar(v)):
                raise TypeError("Column assignment doesn't support type "
                                "{0}".format(type(v).__name__))
        pairs = list(sum(kwargs.items(), ()))

        # Figure out columns of the output
        df2 = self._meta.assign(**_extract_meta(kwargs))
        return elemwise(methods.assign, self, *pairs, meta=df2)

    @derived_from(pd.DataFrame)
    def rename(self, index=None, columns=None):
        if index is not None:
            raise ValueError("Cannot rename index.")

        # *args here is index, columns but columns arg is already used
        return self.map_partitions(M.rename, None, columns=columns)

    def query(self, expr, **kwargs):
        """ Filter dataframe with complex expression

        Blocked version of pd.DataFrame.query

        This is like the sequential version except that this will also happen
        in many threads.  This may conflict with ``numexpr`` which will use
        multiple threads itself.  We recommend that you set numexpr to use a
        single thread

            import numexpr
            numexpr.set_nthreads(1)

        See also
        --------
        pandas.DataFrame.query
        """
        return self.map_partitions(M.query, expr, **kwargs)

    @derived_from(pd.DataFrame)
    def eval(self, expr, inplace=None, **kwargs):
        if inplace is None:
            if PANDAS_VERSION >= '0.21.0':
                inplace = False
        if '=' in expr and inplace in (True, None):
            raise NotImplementedError("Inplace eval not supported."
                                      " Please use inplace=False")
        meta = self._meta.eval(expr, inplace=inplace, **kwargs)
        return self.map_partitions(M.eval, expr, meta=meta, inplace=inplace, **kwargs)

    @derived_from(pd.DataFrame)
    def dropna(self, how='any', subset=None):
        return self.map_partitions(M.dropna, how=how, subset=subset)

    @derived_from(pd.DataFrame)
    def clip(self, lower=None, upper=None, out=None):
        if out is not None:
            raise ValueError("'out' must be None")
        return self.map_partitions(M.clip, lower=lower, upper=upper)

    @derived_from(pd.DataFrame)
    def clip_lower(self, threshold):
        return self.map_partitions(M.clip_lower, threshold=threshold)

    @derived_from(pd.DataFrame)
    def clip_upper(self, threshold):
        return self.map_partitions(M.clip_upper, threshold=threshold)

    @derived_from(pd.DataFrame)
    def squeeze(self, axis=None):
        if axis in [None, 1]:
            if len(self.columns) == 1:
                return self[self.columns[0]]
            else:
                return self

        elif axis == 0:
            raise NotImplementedError("{0} does not support "
                                      "squeeze along axis 0".format(type(self)))

        elif axis not in [0, 1, None]:
            raise ValueError('No axis {0} for object type {1}'.format(
                axis, type(self)))

    @derived_from(pd.DataFrame)
    def to_timestamp(self, freq=None, how='start', axis=0):
        df = elemwise(M.to_timestamp, self, freq, how, axis)
        df.divisions = tuple(pd.Index(self.divisions).to_timestamp())
        return df

    def to_bag(self, index=False):
        """Convert to a dask Bag of tuples of each row.

        Parameters
        ----------
        index : bool, optional
            If True, the index is included as the first element of each tuple.
            Default is False.
        """
        from .io import to_bag
        return to_bag(self, index)

    @derived_from(pd.DataFrame)
    def to_string(self, max_rows=5):
        # option_context doesn't affect
        return self._repr_data.to_string(max_rows=max_rows,
                                         show_dimensions=False)

    def _get_numeric_data(self, how='any', subset=None):
        # calculate columns to avoid unnecessary calculation
        numerics = self._meta._get_numeric_data()

        if len(numerics.columns) < len(self.columns):
            name = self._token_prefix + '-get_numeric_data'
            return self.map_partitions(M._get_numeric_data,
                                       meta=numerics, token=name)
        else:
            # use myself if all numerics
            return self

    @classmethod
    def _validate_axis(cls, axis=0):
        if axis not in (0, 1, 'index', 'columns', None):
            raise ValueError('No axis named {0}'.format(axis))
        # convert to numeric axis
        return {None: 0, 'index': 0, 'columns': 1}.get(axis, axis)

    @derived_from(pd.DataFrame)
    def drop(self, labels, axis=0, errors='raise'):
        axis = self._validate_axis(axis)
        if axis == 1:
            return self.map_partitions(M.drop, labels, axis=axis, errors=errors)
        raise NotImplementedError("Drop currently only works for axis=1")

    @derived_from(pd.DataFrame)
    def merge(self, right, how='inner', on=None, left_on=None, right_on=None,
              left_index=False, right_index=False, suffixes=('_x', '_y'),
              indicator=False, npartitions=None, shuffle=None):

        if not isinstance(right, (DataFrame, pd.DataFrame)):
            raise ValueError('right must be DataFrame')

        from .multi import merge
        return merge(self, right, how=how, on=on, left_on=left_on,
                     right_on=right_on, left_index=left_index,
                     right_index=right_index, suffixes=suffixes,
                     npartitions=npartitions, indicator=indicator,
                     shuffle=shuffle)

    @derived_from(pd.DataFrame)
    def join(self, other, on=None, how='left',
             lsuffix='', rsuffix='', npartitions=None, shuffle=None):

        if not isinstance(other, (DataFrame, pd.DataFrame)):
            raise ValueError('other must be DataFrame')

        from .multi import merge
        return merge(self, other, how=how,
                     left_index=on is None, right_index=True,
                     left_on=on, suffixes=[lsuffix, rsuffix],
                     npartitions=npartitions, shuffle=shuffle)

    @derived_from(pd.DataFrame)
    def append(self, other):
        if isinstance(other, Series):
            msg = ('Unable to appending dd.Series to dd.DataFrame.'
                   'Use pd.Series to append as row.')
            raise ValueError(msg)
        elif isinstance(other, pd.Series):
            other = other.to_frame().T
        return super(DataFrame, self).append(other)

    @derived_from(pd.DataFrame)
    def iterrows(self):
        for i in range(self.npartitions):
            df = self.get_partition(i).compute()
            for row in df.iterrows():
                yield row

    @derived_from(pd.DataFrame)
    def itertuples(self):
        for i in range(self.npartitions):
            df = self.get_partition(i).compute()
            for row in df.itertuples():
                yield row

    @classmethod
    def _bind_operator_method(cls, name, op):
        """ bind operator method like DataFrame.add to this class """

        # name must be explicitly passed for div method whose name is truediv

        def meth(self, other, axis='columns', level=None, fill_value=None):
            if level is not None:
                raise NotImplementedError('level must be None')

            axis = self._validate_axis(axis)

            if axis in (1, 'columns'):
                # When axis=1 and other is a series, `other` is transposed
                # and the operator is applied broadcast across rows. This
                # isn't supported with dd.Series.
                if isinstance(other, Series):
                    msg = 'Unable to {0} dd.Series with axis=1'.format(name)
                    raise ValueError(msg)
                elif isinstance(other, pd.Series):
                    # Special case for pd.Series to avoid unwanted partitioning
                    # of other. We pass it in as a kwarg to prevent this.
                    meta = _emulate(op, self, other=other, axis=axis,
                                    fill_value=fill_value)
                    return map_partitions(op, self, other=other, meta=meta,
                                          axis=axis, fill_value=fill_value)

            meta = _emulate(op, self, other, axis=axis, fill_value=fill_value)
            return map_partitions(op, self, other, meta=meta,
                                  axis=axis, fill_value=fill_value)
        meth.__doc__ = op.__doc__
        bind_method(cls, name, meth)

    @classmethod
    def _bind_comparison_method(cls, name, comparison):
        """ bind comparison method like DataFrame.add to this class """

        def meth(self, other, axis='columns', level=None):
            if level is not None:
                raise NotImplementedError('level must be None')
            axis = self._validate_axis(axis)
            return elemwise(comparison, self, other, axis=axis)

        meth.__doc__ = comparison.__doc__
        bind_method(cls, name, meth)

    @insert_meta_param_description(pad=12)
    def apply(self, func, axis=0, broadcast=None, raw=False, reduce=None,
              args=(), meta=no_default, **kwds):
        """ Parallel version of pandas.DataFrame.apply

        This mimics the pandas version except for the following:

        1.  Only ``axis=1`` is supported (and must be specified explicitly).
        2.  The user should provide output metadata via the `meta` keyword.

        Parameters
        ----------
        func : function
            Function to apply to each column/row
        axis : {0 or 'index', 1 or 'columns'}, default 0
            - 0 or 'index': apply function to each column (NOT SUPPORTED)
            - 1 or 'columns': apply function to each row
        $META
        args : tuple
            Positional arguments to pass to function in addition to the array/series

        Additional keyword arguments will be passed as keywords to the function

        Returns
        -------
        applied : Series or DataFrame

        Examples
        --------
        >>> import dask.dataframe as dd
        >>> df = pd.DataFrame({'x': [1, 2, 3, 4, 5],
        ...                    'y': [1., 2., 3., 4., 5.]})
        >>> ddf = dd.from_pandas(df, npartitions=2)

        Apply a function to row-wise passing in extra arguments in ``args`` and
        ``kwargs``:

        >>> def myadd(row, a, b=1):
        ...     return row.sum() + a + b
        >>> res = ddf.apply(myadd, axis=1, args=(2,), b=1.5)

        By default, dask tries to infer the output metadata by running your
        provided function on some fake data. This works well in many cases, but
        can sometimes be expensive, or even fail. To avoid this, you can
        manually specify the output metadata with the ``meta`` keyword. This
        can be specified in many forms, for more information see
        ``dask.dataframe.utils.make_meta``.

        Here we specify the output is a Series with name ``'x'``, and dtype
        ``float64``:

        >>> res = ddf.apply(myadd, axis=1, args=(2,), b=1.5, meta=('x', 'f8'))

        In the case where the metadata doesn't change, you can also pass in
        the object itself directly:

        >>> res = ddf.apply(lambda row: row + 1, axis=1, meta=ddf)

        See Also
        --------
        dask.DataFrame.map_partitions
        """

        axis = self._validate_axis(axis)
        pandas_kwargs = {
            'axis': axis,
            'broadcast': broadcast,
            'raw': raw,
            'reduce': None,
        }

        if PANDAS_VERSION >= '0.23.0':
            kwds.setdefault('result_type', None)

        kwds.update(pandas_kwargs)

        if axis == 0:
            msg = ("dd.DataFrame.apply only supports axis=1\n"
                   "  Try: df.apply(func, axis=1)")
            raise NotImplementedError(msg)

        if meta is no_default:
            msg = ("`meta` is not specified, inferred from partial data. "
                   "Please provide `meta` if the result is unexpected.\n"
                   "  Before: .apply(func)\n"
                   "  After:  .apply(func, meta={'x': 'f8', 'y': 'f8'}) for dataframe result\n"
                   "  or:     .apply(func, meta=('x', 'f8'))            for series result")
            warnings.warn(msg)

            meta = _emulate(M.apply, self._meta_nonempty, func,
                            args=args, udf=True, **kwds)

        return map_partitions(M.apply, self, func, args=args, meta=meta, **kwds)

    @derived_from(pd.DataFrame)
    def applymap(self, func, meta='__no_default__'):
        return elemwise(M.applymap, self, func, meta=meta)

    @derived_from(pd.DataFrame)
    def round(self, decimals=0):
        return elemwise(M.round, self, decimals)

    @derived_from(pd.DataFrame)
    def cov(self, min_periods=None, split_every=False):
        return cov_corr(self, min_periods, split_every=split_every)

    @derived_from(pd.DataFrame)
    def corr(self, method='pearson', min_periods=None, split_every=False):
        if method != 'pearson':
            raise NotImplementedError("Only Pearson correlation has been "
                                      "implemented")
        return cov_corr(self, min_periods, True, split_every=split_every)

    def info(self, buf=None, verbose=False, memory_usage=False):
        """
        Concise summary of a Dask DataFrame.
        """

        if buf is None:
            import sys
            buf = sys.stdout

        lines = [str(type(self))]

        if len(self.columns) == 0:
            lines.append('Index: 0 entries')
            lines.append('Empty %s' % type(self).__name__)
            put_lines(buf, lines)
            return

        # Group and execute the required computations
        computations = {}
        if verbose:
            computations.update({'index': self.index, 'count': self.count()})
        if memory_usage:
            computations.update({'memory_usage': self.map_partitions(M.memory_usage, index=True)})
        computations = dict(zip(computations.keys(), da.compute(*computations.values())))

        if verbose:
            index = computations['index']
            counts = computations['count']
            lines.append(index_summary(index))
            lines.append('Data columns (total {} columns):'.format(len(self.columns)))

            if PANDAS_VERSION >= '0.20.0':
                from pandas.io.formats.printing import pprint_thing
            else:
                from pandas.formats.printing import pprint_thing
            space = max([len(pprint_thing(k)) for k in self.columns]) + 3
            column_template = '{!s:<%d} {} non-null {}' % space
            column_info = [column_template.format(pprint_thing(x[0]), x[1], x[2])
                           for x in zip(self.columns, counts, self.dtypes)]
        else:
            column_info = [index_summary(self.columns, name='Columns')]

        lines.extend(column_info)
        dtype_counts = ['%s(%d)' % k for k in sorted(self.dtypes.value_counts().iteritems(), key=str)]
        lines.append('dtypes: {}'.format(', '.join(dtype_counts)))

        if memory_usage:
            memory_int = computations['memory_usage'].sum()
            lines.append('memory usage: {}\n'.format(memory_repr(memory_int)))

        put_lines(buf, lines)

    @derived_from(pd.DataFrame)
    def memory_usage(self, index=True, deep=False):
        result = self.map_partitions(M.memory_usage, index=index, deep=deep)
        result = result.groupby(result.index).sum()
        return result

    def pivot_table(self, index=None, columns=None,
                    values=None, aggfunc='mean'):
        """
        Create a spreadsheet-style pivot table as a DataFrame. Target ``columns``
        must have category dtype to infer result's ``columns``.
        ``index``, ``columns``, ``values`` and ``aggfunc`` must be all scalar.

        Parameters
        ----------
        values : scalar
            column to aggregate
        index : scalar
            column to be index
        columns : scalar
            column to be columns
        aggfunc : {'mean', 'sum', 'count'}, default 'mean'

        Returns
        -------
        table : DataFrame
        """
        from .reshape import pivot_table
        return pivot_table(self, index=index, columns=columns, values=values,
                           aggfunc=aggfunc)

    def to_records(self, index=False):
        from .io import to_records
        return to_records(self)

    @derived_from(pd.DataFrame)
    def to_html(self, max_rows=5):
        # pd.Series doesn't have html repr
        data = self._repr_data.to_html(max_rows=max_rows,
                                       show_dimensions=False)
        return self._HTML_FMT.format(data=data, name=key_split(self._name),
                                     task=len(self.dask))

    @property
    def _repr_data(self):
        meta = self._meta
        index = self._repr_divisions
        values = {c: _repr_data_series(meta[c], index) for c in meta.columns}
        return pd.DataFrame(values, columns=meta.columns)

    _HTML_FMT = """<div><strong>Dask DataFrame Structure:</strong></div>
{data}
<div>Dask Name: {name}, {task} tasks</div>"""

    def _repr_html_(self):
        data = self._repr_data.to_html(max_rows=5,
                                       show_dimensions=False, notebook=True)
        return self._HTML_FMT.format(data=data, name=key_split(self._name),
                                     task=len(self.dask))

    def _select_columns_or_index(self, columns_or_index):
        """
        Parameters
        ----------
        columns_or_index
            Column or index name, or a list of these

        Returns
        -------
        dd.DataFrame
            Dask DataFrame with columns corresponding to each column or
            index level in columns_or_index.  If included, the column
            corresponding to the index level is named _index
        """

        # Ensure columns_or_index is a list
        columns_or_index = (columns_or_index
                            if isinstance(columns_or_index, list)
                            else [columns_or_index])

        column_names = [n for n in columns_or_index if self._is_column_label_reference(n)]

        selected_df = self[column_names]
        if self._contains_index_name(columns_or_index):
            # Index name was included
            selected_df = selected_df.assign(_index=self.index)

        return selected_df

    def _is_column_label_reference(self, key):
        """
        Test whether a key is a column label reference

        To be considered a column label reference, `key` must match the name of at
        least one column.
        """
        return (not is_dask_collection(key) and
                (np.isscalar(key) or isinstance(key, tuple)) and
                key in self.columns)

    def _is_index_level_reference(self, key):
        """
        Test whether a key is an index level reference

        To be considered an index level reference, `key` must match the index name
        and must NOT match the name of any column.
        """
        return (self.index.name is not None and
                not is_dask_collection(key) and
                (np.isscalar(key) or isinstance(key, tuple)) and
                key == self.index.name and
                key not in self.columns)

    def _contains_index_name(self, columns_or_index):
        """
        Test whether the input contains a reference to the index of the DataFrame
        """
        if isinstance(columns_or_index, list):
            return any(self._is_index_level_reference(n) for n in columns_or_index)
        else:
            return self._is_index_level_reference(columns_or_index)


# bind operators
for op in [operator.abs, operator.add, operator.and_, operator_div,
           operator.eq, operator.gt, operator.ge, operator.inv,
           operator.lt, operator.le, operator.mod, operator.mul,
           operator.ne, operator.neg, operator.or_, operator.pow,
           operator.sub, operator.truediv, operator.floordiv, operator.xor]:
    _Frame._bind_operator(op)
    Scalar._bind_operator(op)

for name in ['add', 'sub', 'mul', 'div',
             'truediv', 'floordiv', 'mod', 'pow',
             'radd', 'rsub', 'rmul', 'rdiv',
             'rtruediv', 'rfloordiv', 'rmod', 'rpow']:
    meth = getattr(pd.DataFrame, name)
    DataFrame._bind_operator_method(name, meth)

    meth = getattr(pd.Series, name)
    Series._bind_operator_method(name, meth)

for name in ['lt', 'gt', 'le', 'ge', 'ne', 'eq']:
    meth = getattr(pd.DataFrame, name)
    DataFrame._bind_comparison_method(name, meth)

    meth = getattr(pd.Series, name)
    Series._bind_comparison_method(name, meth)


def is_broadcastable(dfs, s):
    """
    This Series is broadcastable against another dataframe in the sequence
    """
    return (isinstance(s, Series) and
            s.npartitions == 1 and
            s.known_divisions and
            any(s.divisions == (min(df.columns), max(df.columns))
                for df in dfs if isinstance(df, DataFrame)))


def elemwise(op, *args, **kwargs):
    """ Elementwise operation for Dask dataframes

    Parameters
    ----------
    op: callable
        Function to apply across input dataframes
    *args: DataFrames, Series, Scalars, Arrays,
        The arguments of the operation
    **kwrags: scalars
    meta: pd.DataFrame, pd.Series (optional)
        Valid metadata for the operation.  Will evaluate on a small piece of
        data if not provided.

    Examples
    --------
    >>> elemwise(operator.add, df.x, df.y)  # doctest: +SKIP
    """
    meta = kwargs.pop('meta', no_default)
    out = kwargs.pop('out', None)

    _name = funcname(op) + '-' + tokenize(op, *args, **kwargs)

    args = _maybe_from_pandas(args)

    from .multi import _maybe_align_partitions
    args = _maybe_align_partitions(args)
    dasks = [arg for arg in args if isinstance(arg, (_Frame, Scalar, Array))]
    dfs = [df for df in dasks if isinstance(df, _Frame)]

    # Clean up dask arrays if present
    for i, a in enumerate(dasks):
        if not isinstance(a, Array):
            continue
        # Ensure that they have similar-ish chunk structure
        if not all(len(a.chunks[0]) == df.npartitions for df in dfs):
            msg = ("When combining dask arrays with dataframes they must "
                   "match chunking exactly.  Operation: %s" % funcname(op))
            raise ValueError(msg)
        # Rechunk to have a single chunk along all other axes
        if a.ndim > 1:
            a = a.rechunk({i + 1: d for i, d in enumerate(a.shape[1:])})
            dasks[i] = a

    divisions = dfs[0].divisions
    _is_broadcastable = partial(is_broadcastable, dfs)
    dfs = list(remove(_is_broadcastable, dfs))
    n = len(divisions) - 1

    other = [(i, arg) for i, arg in enumerate(args)
             if not isinstance(arg, (_Frame, Scalar, Array))]

    # adjust the key length of Scalar
    keys = [d.__dask_keys__() * n
            if isinstance(d, Scalar) or _is_broadcastable(d)
            else core.flatten(d.__dask_keys__()) for d in dasks]

    if other:
        dsk = {(_name, i):
               (apply, partial_by_order, list(frs),
                {'function': op, 'other': other})
               for i, frs in enumerate(zip(*keys))}
    else:
        dsk = {(_name, i): (op,) + frs for i, frs in enumerate(zip(*keys))}
    dsk = merge(dsk, *[d.dask for d in dasks])

    if meta is no_default:
        if len(dfs) >= 2 and not all(hasattr(d, 'npartitions') for d in dasks):
            # should not occur in current funcs
            msg = 'elemwise with 2 or more DataFrames and Scalar is not supported'
            raise NotImplementedError(msg)
        # For broadcastable series, use no rows.
        parts = [d._meta if _is_broadcastable(d) or isinstance(d, Array)
                 else d._meta_nonempty for d in dasks]
        with raise_on_meta_error(funcname(op)):
            meta = partial_by_order(*parts, function=op, other=other)

    result = new_dd_object(dsk, _name, meta, divisions)
    return handle_out(out, result)


def handle_out(out, result):
    """ Handle out parameters

    If out is a dask.DataFrame, dask.Series or dask.Scalar then
    this overwrites the contents of it with the result
    """
    if isinstance(out, tuple):
        if len(out) == 1:
            out = out[0]
        elif len(out) > 1:
            raise NotImplementedError("The out parameter is not fully supported")
        else:
            out = None

    if out is not None and type(out) != type(result):
        raise TypeError(
            "Mismatched types between result and out parameter. "
            "out=%s, result=%s" % (str(type(out)), str(type(result))))

    if isinstance(out, DataFrame):
        if len(out.columns) != len(result.columns):
            raise ValueError(
                "Mismatched columns count between result and out parameter. "
                "out=%s, result=%s" % (str(len(out.columns)), str(len(result.columns))))

    if isinstance(out, (Series, DataFrame, Scalar)):
        out._meta = result._meta
        out._name = result._name
        out.dask = result.dask

        if not isinstance(out, Scalar):
            out.divisions = result.divisions
    elif out is not None:
        msg = ("The out parameter is not fully supported."
               " Received type %s, expected %s " % ( type(out).__name__, type(result).__name__))
        raise NotImplementedError(msg)
    else:
        return result


def _maybe_from_pandas(dfs):
    from .io import from_pandas
    dfs = [from_pandas(df, 1) if isinstance(df, (pd.Series, pd.DataFrame))
           else df for df in dfs]
    return dfs


def hash_shard(df, nparts, split_out_setup=None, split_out_setup_kwargs=None):
    if split_out_setup:
        h = split_out_setup(df, **(split_out_setup_kwargs or {}))
    else:
        h = df
    h = hash_pandas_object(h, index=False)
    if isinstance(h, pd.Series):
        h = h._values
    h %= nparts
    return {i: df.iloc[h == i] for i in range(nparts)}


def split_evenly(df, k):
    """ Split dataframe into k roughly equal parts """
    divisions = np.linspace(0, len(df), k + 1).astype(int)
    return {i: df.iloc[divisions[i]: divisions[i + 1]] for i in range(k)}


def split_out_on_index(df):
    h = df.index
    if isinstance(h, pd.MultiIndex):
        h = pd.DataFrame([], index=h).reset_index()
    return h


def split_out_on_cols(df, cols=None):
    return df[cols]


@insert_meta_param_description
def apply_concat_apply(args, chunk=None, aggregate=None, combine=None,
                       meta=no_default, token=None, chunk_kwargs=None,
                       aggregate_kwargs=None, combine_kwargs=None,
                       split_every=None, split_out=None, split_out_setup=None,
                       split_out_setup_kwargs=None, **kwargs):
    """Apply a function to blocks, then concat, then apply again

    Parameters
    ----------
    args :
        Positional arguments for the `chunk` function. All `dask.dataframe`
        objects should be partitioned and indexed equivalently.
    chunk : function [block-per-arg] -> block
        Function to operate on each block of data
    aggregate : function concatenated-block -> block
        Function to operate on the concatenated result of chunk
    combine : function concatenated-block -> block, optional
        Function to operate on intermediate concatenated results of chunk
        in a tree-reduction. If not provided, defaults to aggregate.
    $META
    token : str, optional
        The name to use for the output keys.
    chunk_kwargs : dict, optional
        Keywords for the chunk function only.
    aggregate_kwargs : dict, optional
        Keywords for the aggregate function only.
    combine_kwargs : dict, optional
        Keywords for the combine function only.
    split_every : int, optional
        Group partitions into groups of this size while performing a
        tree-reduction. If set to False, no tree-reduction will be used,
        and all intermediates will be concatenated and passed to ``aggregate``.
        Default is 8.
    split_out : int, optional
        Number of output partitions. Split occurs after first chunk reduction.
    split_out_setup : callable, optional
        If provided, this function is called on each chunk before performing
        the hash-split. It should return a pandas object, where each row
        (excluding the index) is hashed. If not provided, the chunk is hashed
        as is.
    split_out_setup_kwargs : dict, optional
        Keywords for the `split_out_setup` function only.
    kwargs :
        All remaining keywords will be passed to ``chunk``, ``aggregate``, and
        ``combine``.

    Examples
    --------
    >>> def chunk(a_block, b_block):
    ...     pass

    >>> def agg(df):
    ...     pass

    >>> apply_concat_apply([a, b], chunk=chunk, aggregate=agg)  # doctest: +SKIP
    """
    if chunk_kwargs is None:
        chunk_kwargs = dict()
    if aggregate_kwargs is None:
        aggregate_kwargs = dict()
    chunk_kwargs.update(kwargs)
    aggregate_kwargs.update(kwargs)

    if combine is None:
        if combine_kwargs:
            raise ValueError("`combine_kwargs` provided with no `combine`")
        combine = aggregate
        combine_kwargs = aggregate_kwargs
    else:
        if combine_kwargs is None:
            combine_kwargs = dict()
        combine_kwargs.update(kwargs)

    if not isinstance(args, (tuple, list)):
        args = [args]

    npartitions = set(arg.npartitions for arg in args
                      if isinstance(arg, _Frame))
    if len(npartitions) > 1:
        raise ValueError("All arguments must have same number of partitions")
    npartitions = npartitions.pop()

    if split_every is None:
        split_every = 8
    elif split_every is False:
        split_every = npartitions
    elif split_every < 2 or not isinstance(split_every, int):
        raise ValueError("split_every must be an integer >= 2")

    token_key = tokenize(token or (chunk, aggregate), meta, args,
                         chunk_kwargs, aggregate_kwargs, combine_kwargs,
                         split_every, split_out, split_out_setup,
                         split_out_setup_kwargs)

    # Chunk
    a = '{0}-chunk-{1}'.format(token or funcname(chunk), token_key)
    if len(args) == 1 and isinstance(args[0], _Frame) and not chunk_kwargs:
        dsk = {(a, 0, i, 0): (chunk, key)
               for i, key in enumerate(args[0].__dask_keys__())}
    else:
        dsk = {(a, 0, i, 0): (apply, chunk,
                              [(x._name, i) if isinstance(x, _Frame)
                               else x for x in args], chunk_kwargs)
               for i in range(args[0].npartitions)}

    # Split
    if split_out and split_out > 1:
        split_prefix = 'split-%s' % token_key
        shard_prefix = 'shard-%s' % token_key
        for i in range(args[0].npartitions):
            dsk[(split_prefix, i)] = (hash_shard, (a, 0, i, 0), split_out,
                                      split_out_setup, split_out_setup_kwargs)
            for j in range(split_out):
                dsk[(shard_prefix, 0, i, j)] = (getitem, (split_prefix, i), j)
        a = shard_prefix
    else:
        split_out = 1

    # Combine
    b = '{0}-combine-{1}'.format(token or funcname(combine), token_key)
    k = npartitions
    depth = 0
    while k > split_every:
        for part_i, inds in enumerate(partition_all(split_every, range(k))):
            for j in range(split_out):
                conc = (_concat, [(a, depth, i, j) for i in inds])
                if combine_kwargs:
                    dsk[(b, depth + 1, part_i, j)] = (apply, combine, [conc], combine_kwargs)
                else:
                    dsk[(b, depth + 1, part_i, j)] = (combine, conc)
        k = part_i + 1
        a = b
        depth += 1

    # Aggregate
    for j in range(split_out):
        b = '{0}-agg-{1}'.format(token or funcname(aggregate), token_key)
        conc = (_concat, [(a, depth, i, j) for i in range(k)])
        if aggregate_kwargs:
            dsk[(b, j)] = (apply, aggregate, [conc], aggregate_kwargs)
        else:
            dsk[(b, j)] = (aggregate, conc)

    if meta is no_default:
        meta_chunk = _emulate(chunk, *args, udf=True, **chunk_kwargs)
        meta = _emulate(aggregate, _concat([meta_chunk]), udf=True,
                        **aggregate_kwargs)
    meta = make_meta(meta)

    for arg in args:
        if isinstance(arg, _Frame):
            dsk.update(arg.dask)

    divisions = [None] * (split_out + 1)

    return new_dd_object(dsk, b, meta, divisions)


aca = apply_concat_apply


def _extract_meta(x, nonempty=False):
    """
    Extract internal cache data (``_meta``) from dd.DataFrame / dd.Series
    """
    if isinstance(x, (Scalar, _Frame)):
        return x._meta_nonempty if nonempty else x._meta
    elif isinstance(x, list):
        return [_extract_meta(_x, nonempty) for _x in x]
    elif isinstance(x, tuple):
        return tuple([_extract_meta(_x, nonempty) for _x in x])
    elif isinstance(x, dict):
        res = {}
        for k in x:
            res[k] = _extract_meta(x[k], nonempty)
        return res
    else:
        return x


def _emulate(func, *args, **kwargs):
    """
    Apply a function using args / kwargs. If arguments contain dd.DataFrame /
    dd.Series, using internal cache (``_meta``) for calculation
    """
    with raise_on_meta_error(funcname(func), udf=kwargs.pop('udf', False)):
        return func(*_extract_meta(args, True), **_extract_meta(kwargs, True))


@insert_meta_param_description
def map_partitions(func, *args, **kwargs):
    """ Apply Python function on each DataFrame partition.

    Parameters
    ----------
    func : function
        Function applied to each partition.
    args, kwargs :
        Arguments and keywords to pass to the function.  At least one of the
        args should be a Dask.dataframe. Arguments and keywords may contain
        ``Scalar``, ``Delayed`` or regular python objects.
    $META
    """
    meta = kwargs.pop('meta', no_default)

    if meta is not no_default:
        meta = make_meta(meta)

    assert callable(func)
    if 'token' in kwargs:
        name = kwargs.pop('token')
        token = tokenize(meta, *args, **kwargs)
    else:
        name = funcname(func)
        token = tokenize(func, meta, *args, **kwargs)
    name = '{0}-{1}'.format(name, token)

    from .multi import _maybe_align_partitions
    args = _maybe_from_pandas(args)
    args = _maybe_align_partitions(args)

    if meta is no_default:
        meta = _emulate(func, *args, udf=True, **kwargs)

    if all(isinstance(arg, Scalar) for arg in args):
        dask = {(name, 0):
                (apply, func, (tuple, [(arg._name, 0) for arg in args]), kwargs)}
        return Scalar(merge(dask, *[arg.dask for arg in args]), name, meta)
    elif not (isinstance(meta, (pd.Series, pd.DataFrame, pd.Index)) or is_arraylike(meta)):
        # If `meta` is not a pandas object, the concatenated results will be a
        # different type
        meta = _concat([meta])
    meta = make_meta(meta)

    args, args_dasks = _process_lazy_args(args)
    kwargs_task, kwargs_dsk = to_task_dask(kwargs)
    args_dasks.append(kwargs_dsk)

    dfs = [df for df in args if isinstance(df, _Frame)]
    dsk = {}
    for i in range(dfs[0].npartitions):
        values = [(arg._name, i if isinstance(arg, _Frame) else 0)
                  if isinstance(arg, (_Frame, Scalar)) else arg for arg in args]
        dsk[(name, i)] = (apply_and_enforce, func, values, kwargs_task, meta)

    args_dasks.extend([arg.dask for arg in args if isinstance(arg, (_Frame, Scalar))])

    return new_dd_object(merge(dsk, *args_dasks), name, meta, args[0].divisions)


def _process_lazy_args(args):
    # NOTE: we use this function instead of to_task_dask to avoid
    # manipulating _Frame instances that need to be aligned
    dsk = [arg.dask for arg in args if isinstance(arg, Delayed)]
    args = [arg._key if isinstance(arg, Delayed) else arg for arg in args ]

    return args, dsk


def apply_and_enforce(func, args, kwargs, meta):
    """Apply a function, and enforce the output to match meta

    Ensures the output has the same columns, even if empty."""
    df = func(*args, **kwargs)
    if isinstance(df, (pd.DataFrame, pd.Series, pd.Index)):
        if len(df) == 0:
            return meta
        c = meta.columns if isinstance(df, pd.DataFrame) else meta.name
        return _rename(c, df)
    return df


def _rename(columns, df):
    """
    Rename columns of pd.DataFrame or name of pd.Series.
    Not for dd.DataFrame or dd.Series.

    Parameters
    ----------
    columns : tuple, string, pd.DataFrame or pd.Series
        Column names, Series name or pandas instance which has the
        target column names / name.
    df : pd.DataFrame or pd.Series
        target DataFrame / Series to be renamed
    """
    assert not isinstance(df, _Frame)

    if columns is no_default:
        return df

    if isinstance(columns, Iterator):
        columns = list(columns)

    if isinstance(df, pd.DataFrame):
        if isinstance(columns, pd.DataFrame):
            columns = columns.columns
        if not isinstance(columns, pd.Index):
            columns = pd.Index(columns)
        if (len(columns) == len(df.columns) and
                type(columns) is type(df.columns) and
                columns.equals(df.columns)):
            # if target is identical, rename is not necessary
            return df
        # deep=False doesn't doesn't copy any data/indices, so this is cheap
        df = df.copy(deep=False)
        df.columns = columns
        return df
    elif isinstance(df, (pd.Series, pd.Index)):
        if isinstance(columns, (pd.Series, pd.Index)):
            columns = columns.name
        if df.name == columns:
            return df
        return df.rename(columns)
    # map_partition may pass other types
    return df


def _rename_dask(df, names):
    """
    Destructively rename columns of dd.DataFrame or name of dd.Series.
    Not for pd.DataFrame or pd.Series.

    Internaly used to overwrite dd.DataFrame.columns and dd.Series.name
    We can't use map_partition because it applies function then rename

    Parameters
    ----------
    df : dd.DataFrame or dd.Series
        target DataFrame / Series to be renamed
    names : tuple, string
        Column names/Series name
    """

    assert isinstance(df, _Frame)
    metadata = _rename(names, df._meta)
    name = 'rename-{0}'.format(tokenize(df, metadata))

    dsk = {}
    for i in range(df.npartitions):
        dsk[name, i] = (_rename, metadata, (df._name, i))
    return new_dd_object(merge(dsk, df.dask), name, metadata, df.divisions)


def quantile(df, q):
    """Approximate quantiles of Series.

    Parameters
    ----------
    q : list/array of floats
        Iterable of numbers ranging from 0 to 100 for the desired quantiles
    """
    assert isinstance(df, Series)
    from dask.array.percentile import _percentile, merge_percentiles

    # currently, only Series has quantile method
    if isinstance(df, Index):
        meta = pd.Series(df._meta_nonempty).quantile(q)
    else:
        meta = df._meta_nonempty.quantile(q)

    if isinstance(meta, pd.Series):
        # Index.quantile(list-like) must be pd.Series, not pd.Index
        df_name = df.name
        finalize_tsk = lambda tsk: (pd.Series, tsk, q, None, df_name)
        return_type = Series
    else:
        finalize_tsk = lambda tsk: (getitem, tsk, 0)
        return_type = Scalar
        q = [q]

    # pandas uses quantile in [0, 1]
    # numpy / everyone else uses [0, 100]
    qs = np.asarray(q) * 100
    token = tokenize(df, qs)

    if len(qs) == 0:
        name = 'quantiles-' + token
        empty_index = pd.Index([], dtype=float)
        return Series({(name, 0): pd.Series([], name=df.name, index=empty_index)},
                      name, df._meta, [None, None])
    else:
        new_divisions = [np.min(q), np.max(q)]

    df = df.dropna()
    name = 'quantiles-1-' + token
    val_dsk = {(name, i): (_percentile, (getattr, key, 'values'), qs)
               for i, key in enumerate(df.__dask_keys__())}

    name3 = 'quantiles-3-' + token
    merge_dsk = {(name3, 0): finalize_tsk((merge_percentiles, qs,
                                           [qs] * df.npartitions,
                                           sorted(val_dsk)))}
    dsk = merge(df.dask, val_dsk, merge_dsk)
    return return_type(dsk, name3, meta, new_divisions)


def cov_corr(df, min_periods=None, corr=False, scalar=False, split_every=False):
    """DataFrame covariance and pearson correlation.

    Computes pairwise covariance or correlation of columns, excluding NA/null
    values.

    Parameters
    ----------
    df : DataFrame
    min_periods : int, optional
        Minimum number of observations required per pair of columns
        to have a valid result.
    corr : bool, optional
        If True, compute the Pearson correlation. If False [default], compute
        the covariance.
    scalar : bool, optional
        If True, compute covariance between two variables as a scalar. Only
        valid if `df` has 2 columns.  If False [default], compute the entire
        covariance/correlation matrix.
    split_every : int, optional
        Group partitions into groups of this size while performing a
        tree-reduction. If set to False, no tree-reduction will be used.
        Default is False.
    """
    if min_periods is None:
        min_periods = 2
    elif min_periods < 2:
        raise ValueError("min_periods must be >= 2")

    if split_every is False:
        split_every = df.npartitions
    elif split_every < 2 or not isinstance(split_every, int):
        raise ValueError("split_every must be an integer >= 2")

    df = df._get_numeric_data()

    if scalar and len(df.columns) != 2:
        raise ValueError("scalar only valid for 2 column dataframe")

    token = tokenize(df, min_periods, scalar, split_every)

    funcname = 'corr' if corr else 'cov'
    a = '{0}-chunk-{1}'.format(funcname, df._name)
    dsk = {(a, i): (cov_corr_chunk, f, corr)
           for (i, f) in enumerate(df.__dask_keys__())}

    prefix = '{0}-combine-{1}-'.format(funcname, df._name)
    k = df.npartitions
    b = a
    depth = 0
    while k > split_every:
        b = prefix + str(depth)
        for part_i, inds in enumerate(partition_all(split_every, range(k))):
            dsk[(b, part_i)] = (cov_corr_combine, [(a, i) for i in inds], corr)
        k = part_i + 1
        a = b
        depth += 1

    name = '{0}-{1}'.format(funcname, token)
    dsk[(name, 0)] = (cov_corr_agg, [(a, i) for i in range(k)],
                      df.columns, min_periods, corr, scalar)
    dsk.update(df.dask)
    if scalar:
        return Scalar(dsk, name, 'f8')
    meta = make_meta([(c, 'f8') for c in df.columns], index=df.columns)
    return DataFrame(dsk, name, meta, (df.columns[0], df.columns[-1]))


def cov_corr_chunk(df, corr=False):
    """Chunk part of a covariance or correlation computation"""
    mat = df.astype('float64', copy=False).values
    mask = np.isfinite(mat)
    keep = np.bitwise_and(mask[:, None, :], mask[:, :, None])

    x = np.where(keep, mat[:, None, :], np.nan)
    sums = np.nansum(x, 0)
    counts = keep.astype('int').sum(0)
    cov = df.cov().values
    dtype = [('sum', sums.dtype), ('count', counts.dtype), ('cov', cov.dtype)]
    if corr:
        m = np.nansum((x - sums / np.where(counts, counts, np.nan)) ** 2, 0)
        dtype.append(('m', m.dtype))

    out = np.empty(counts.shape, dtype=dtype)
    out['sum'] = sums
    out['count'] = counts
    out['cov'] = cov * (counts - 1)
    if corr:
        out['m'] = m
    return out


def cov_corr_combine(data, corr=False):
    data = np.concatenate(data).reshape((len(data),) + data[0].shape)
    sums = np.nan_to_num(data['sum'])
    counts = data['count']

    cum_sums = np.cumsum(sums, 0)
    cum_counts = np.cumsum(counts, 0)

    s1 = cum_sums[:-1]
    s2 = sums[1:]
    n1 = cum_counts[:-1]
    n2 = counts[1:]
    with np.errstate(invalid='ignore'):
        d = (s2 / n2) - (s1 / n1)
        C = (np.nansum((n1 * n2) / (n1 + n2) * (d * d.transpose((0, 2, 1))), 0) +
             np.nansum(data['cov'], 0))

    out = np.empty(C.shape, dtype=data.dtype)
    out['sum'] = cum_sums[-1]
    out['count'] = cum_counts[-1]
    out['cov'] = C

    if corr:
        nobs = np.where(cum_counts[-1], cum_counts[-1], np.nan)
        mu = cum_sums[-1] / nobs
        counts_na = np.where(counts, counts, np.nan)
        m = np.nansum(data['m'] + counts * (sums / counts_na - mu) ** 2,
                      axis=0)
        out['m'] = m
    return out


def cov_corr_agg(data, cols, min_periods=2, corr=False, scalar=False):
    out = cov_corr_combine(data, corr)
    counts = out['count']
    C = out['cov']
    C[counts < min_periods] = np.nan
    if corr:
        m2 = out['m']
        den = np.sqrt(m2 * m2.T)
    else:
        den = np.where(counts, counts, np.nan) - 1
    with np.errstate(invalid='ignore', divide='ignore'):
        mat = C / den
    if scalar:
        return mat[0, 1]
    return pd.DataFrame(mat, columns=cols, index=cols)


def pd_split(df, p, random_state=None):
    """ Split DataFrame into multiple pieces pseudorandomly

    >>> df = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6],
    ...                    'b': [2, 3, 4, 5, 6, 7]})

    >>> a, b = pd_split(df, [0.5, 0.5], random_state=123)  # roughly 50/50 split
    >>> a
       a  b
    1  2  3
    2  3  4
    5  6  7
    >>> b
       a  b
    0  1  2
    3  4  5
    4  5  6
    """
    p = list(p)
    index = pseudorandom(len(df), p, random_state)
    return [df.iloc[index == i] for i in range(len(p))]


def _take_last(a, skipna=True):
    """
    take last row (Series) of DataFrame / last value of Series
    considering NaN.

    Parameters
    ----------
    a : pd.DataFrame or pd.Series
    skipna : bool, default True
        Whether to exclude NaN

    """
    if skipna is False:
        return a.iloc[-1]
    else:
        # take last valid value excluding NaN, NaN location may be different
        # in each columns
        group_dummy = np.ones(len(a.index))
        last_row = a.groupby(group_dummy).last()
        if isinstance(a, pd.DataFrame):
            return pd.Series(last_row.values[0], index=a.columns)
        else:
            return last_row.values[0]


def check_divisions(divisions):
    if not isinstance(divisions, (list, tuple)):
        raise ValueError('New division must be list or tuple')
    divisions = list(divisions)
    if divisions != sorted(divisions):
        raise ValueError('New division must be sorted')
    if len(divisions[:-1]) != len(list(unique(divisions[:-1]))):
        msg = 'New division must be unique, except for the last element'
        raise ValueError(msg)


def repartition_divisions(a, b, name, out1, out2, force=False):
    """ dask graph to repartition dataframe by new divisions

    Parameters
    ----------
    a : tuple
        old divisions
    b : tuple, list
        new divisions
    name : str
        name of old dataframe
    out1 : str
        name of temporary splits
    out2 : str
        name of new dataframe
    force : bool, default False
        Allows the expansion of the existing divisions.
        If False then the new divisions lower and upper bounds must be
        the same as the old divisions.

    Examples
    --------
    >>> repartition_divisions([1, 3, 7], [1, 4, 6, 7], 'a', 'b', 'c')  # doctest: +SKIP
    {('b', 0): (<function boundary_slice at ...>, ('a', 0), 1, 3, False),
     ('b', 1): (<function boundary_slice at ...>, ('a', 1), 3, 4, False),
     ('b', 2): (<function boundary_slice at ...>, ('a', 1), 4, 6, False),
     ('b', 3): (<function boundary_slice at ...>, ('a', 1), 6, 7, False)
     ('c', 0): (<function concat at ...>,
                (<type 'list'>, [('b', 0), ('b', 1)])),
     ('c', 1): ('b', 2),
     ('c', 2): ('b', 3)}
    """
    check_divisions(b)

    if len(b) < 2:
        # minimum division is 2 elements, like [0, 0]
        raise ValueError('New division must be longer than 2 elements')

    if force:
        if a[0] < b[0]:
            msg = ('left side of the new division must be equal or smaller '
                   'than old division')
            raise ValueError(msg)
        if a[-1] > b[-1]:
            msg = ('right side of the new division must be equal or larger '
                   'than old division')
            raise ValueError(msg)
    else:
        if a[0] != b[0]:
            msg = 'left side of old and new divisions are different'
            raise ValueError(msg)
        if a[-1] != b[-1]:
            msg = 'right side of old and new divisions are different'
            raise ValueError(msg)

    def _is_single_last_div(x):
        """Whether last division only contains single label"""
        return len(x) >= 2 and x[-1] == x[-2]

    c = [a[0]]
    d = dict()
    low = a[0]

    i, j = 1, 1     # indices for old/new divisions
    k = 0           # index for temp divisions

    last_elem = _is_single_last_div(a)

    # process through old division
    # left part of new division can be processed in this loop
    while (i < len(a) and j < len(b)):
        if a[i] < b[j]:
            # tuple is something like:
            # (methods.boundary_slice, ('from_pandas-#', 0), 3, 4, False))
            d[(out1, k)] = (methods.boundary_slice, (name, i - 1), low, a[i], False)
            low = a[i]
            i += 1
        elif a[i] > b[j]:
            d[(out1, k)] = (methods.boundary_slice, (name, i - 1), low, b[j], False)
            low = b[j]
            j += 1
        else:
            d[(out1, k)] = (methods.boundary_slice, (name, i - 1), low, b[j], False)
            low = b[j]
            i += 1
            j += 1
        c.append(low)
        k += 1

    # right part of new division can remain
    if a[-1] < b[-1] or b[-1] == b[-2]:
        for _j in range(j, len(b)):
            # always use right-most of old division
            # because it may contain last element
            m = len(a) - 2
            d[(out1, k)] = (methods.boundary_slice, (name, m), low, b[_j], False)
            low = b[_j]
            c.append(low)
            k += 1
    else:
        # even if new division is processed through,
        # right-most element of old division can remain
        if last_elem and i < len(a):
            d[(out1, k)] = (methods.boundary_slice, (name, i - 1), a[i], a[i], False)
            k += 1
        c.append(a[-1])

    # replace last element of tuple with True
    d[(out1, k - 1)] = d[(out1, k - 1)][:-1] + (True,)

    i, j = 0, 1

    last_elem = _is_single_last_div(c)

    while j < len(b):
        tmp = []
        while c[i] < b[j]:
            tmp.append((out1, i))
            i += 1
        if last_elem and c[i] == b[-1] and (b[-1] != b[-2] or j == len(b) - 1) and i < k:
            # append if last split is not included
            tmp.append((out1, i))
            i += 1
        if len(tmp) == 0:
            # dummy slice to return empty DataFrame or Series,
            # which retain original data attributes (columns / name)
            d[(out2, j - 1)] = (methods.boundary_slice, (name, 0), a[0], a[0], False)
        elif len(tmp) == 1:
            d[(out2, j - 1)] = tmp[0]
        else:
            if not tmp:
                raise ValueError('check for duplicate partitions\nold:\n%s\n\n'
                                 'new:\n%s\n\ncombined:\n%s'
                                 % (pformat(a), pformat(b), pformat(c)))
            d[(out2, j - 1)] = (methods.concat, tmp)
        j += 1
    return d


def repartition_freq(df, freq=None):
    """ Repartition a timeseries dataframe by a new frequency """
    if not isinstance(df.divisions[0], pd.Timestamp):
        raise TypeError("Can only repartition on frequency for timeseries")
    try:
        start = df.divisions[0].ceil(freq)
    except ValueError:
        start = df.divisions[0]
    divisions = pd.DatetimeIndex(start=start,
                                 end=df.divisions[-1],
                                 freq=freq).tolist()
    if not len(divisions):
        divisions = [df.divisions[0], df.divisions[-1]]
    else:
        if divisions[-1] != df.divisions[-1]:
            divisions.append(df.divisions[-1])
        if divisions[0] != df.divisions[0]:
            divisions = [df.divisions[0]] + divisions

    return df.repartition(divisions=divisions)


def repartition_npartitions(df, npartitions):
    """ Repartition dataframe to a smaller number of partitions """
    new_name = 'repartition-%d-%s' % (npartitions, tokenize(df))
    if df.npartitions == npartitions:
        return df
    elif df.npartitions > npartitions:
        npartitions_ratio = df.npartitions / npartitions
        new_partitions_boundaries = [int(new_partition_index * npartitions_ratio)
                                     for new_partition_index in range(npartitions + 1)]
        dsk = {}
        for new_partition_index in range(npartitions):
            value = (methods.concat,
                     [(df._name, old_partition_index) for old_partition_index in
                      range(new_partitions_boundaries[new_partition_index],
                            new_partitions_boundaries[new_partition_index + 1])])
            dsk[new_name, new_partition_index] = value
        divisions = [df.divisions[new_partition_index]
                     for new_partition_index in new_partitions_boundaries]
        return new_dd_object(merge(df.dask, dsk), new_name, df._meta, divisions)
    else:
        original_divisions = divisions = pd.Series(df.divisions)
        if (df.known_divisions and (np.issubdtype(divisions.dtype, np.datetime64) or
                                    np.issubdtype(divisions.dtype, np.number))):
            if np.issubdtype(divisions.dtype, np.datetime64):
                divisions = divisions.values.astype('float64')

            if isinstance(divisions, pd.Series):
                divisions = divisions.values

            n = len(divisions)
            divisions = np.interp(x=np.linspace(0, n, npartitions + 1),
                                  xp=np.linspace(0, n, n),
                                  fp=divisions)
            if np.issubdtype(original_divisions.dtype, np.datetime64):
                divisions = pd.Series(divisions).astype(original_divisions.dtype).tolist()
            elif np.issubdtype(original_divisions.dtype, np.integer):
                divisions = divisions.astype(original_divisions.dtype)

            if isinstance(divisions, np.ndarray):
                divisions = divisions.tolist()

            divisions = list(divisions)
            divisions[0] = df.divisions[0]
            divisions[-1] = df.divisions[-1]

            return df.repartition(divisions=divisions)
        else:
            ratio = npartitions / df.npartitions
            split_name = 'split-%s' % tokenize(df, npartitions)
            dsk = {}
            last = 0
            j = 0
            for i in range(df.npartitions):
                new = last + ratio
                if i == df.npartitions - 1:
                    k = npartitions - j
                else:
                    k = int(new - last)
                dsk[(split_name, i)] = (split_evenly, (df._name, i), k)
                for jj in range(k):
                    dsk[(new_name, j)] = (getitem, (split_name, i), jj)
                    j += 1
                last = new

            divisions = [None] * (npartitions + 1)
            return new_dd_object(merge(df.dask, dsk), new_name, df._meta, divisions)


def repartition(df, divisions=None, force=False):
    """ Repartition dataframe along new divisions

    Dask.DataFrame objects are partitioned along their index.  Often when
    multiple dataframes interact we need to align these partitionings.  The
    ``repartition`` function constructs a new DataFrame object holding the same
    data but partitioned on different values.  It does this by performing a
    sequence of ``loc`` and ``concat`` calls to split and merge the previous
    generation of partitions.

    Parameters
    ----------

    divisions : list
        List of partitions to be used
    force : bool, default False
        Allows the expansion of the existing divisions.
        If False then the new divisions lower and upper bounds must be
        the same as the old divisions.

    Examples
    --------

    >>> df = df.repartition([0, 5, 10, 20])  # doctest: +SKIP

    Also works on Pandas objects

    >>> ddf = dd.repartition(df, [0, 5, 10, 20])  # doctest: +SKIP
    """

    token = tokenize(df, divisions)
    if isinstance(df, _Frame):
        tmp = 'repartition-split-' + token
        out = 'repartition-merge-' + token
        dsk = repartition_divisions(df.divisions, divisions,
                                    df._name, tmp, out, force=force)
        return new_dd_object(merge(df.dask, dsk), out,
                             df._meta, divisions)
    elif isinstance(df, (pd.Series, pd.DataFrame)):
        name = 'repartition-dataframe-' + token
        from .utils import shard_df_on_index
        dfs = shard_df_on_index(df, divisions[1:-1])
        dsk = dict(((name, i), df) for i, df in enumerate(dfs))
        return new_dd_object(dsk, name, df, divisions)
    raise ValueError('Data must be DataFrame or Series')


def _reduction_chunk(x, aca_chunk=None, **kwargs):
    o = aca_chunk(x, **kwargs)
    # Return a dataframe so that the concatenated version is also a dataframe
    return o.to_frame().T if isinstance(o, pd.Series) else o


def _reduction_combine(x, aca_combine=None, **kwargs):
    if isinstance(x, list):
        x = pd.Series(x)
    o = aca_combine(x, **kwargs)
    # Return a dataframe so that the concatenated version is also a dataframe
    return o.to_frame().T if isinstance(o, pd.Series) else o


def _reduction_aggregate(x, aca_aggregate=None, **kwargs):
    if isinstance(x, list):
        x = pd.Series(x)
    return aca_aggregate(x, **kwargs)


def idxmaxmin_chunk(x, fn=None, skipna=True):
    minmax = 'max' if fn == 'idxmax' else 'min'
    if len(x) > 0:
        idx = getattr(x, fn)(skipna=skipna)
        value = getattr(x, minmax)(skipna=skipna)
    else:
        idx = value = pd.Series([], dtype='i8')
    if isinstance(idx, pd.Series):
        return pd.DataFrame({'idx': idx, 'value': value})
    return pd.DataFrame({'idx': [idx], 'value': [value]})


def idxmaxmin_row(x, fn=None, skipna=True):
    minmax = 'max' if fn == 'idxmax' else 'min'
    if len(x) > 0:
        x = x.set_index('idx')
        idx = [getattr(x.value, fn)(skipna=skipna)]
        value = [getattr(x.value, minmax)(skipna=skipna)]
    else:
        idx = value = pd.Series([], dtype='i8')
    return pd.DataFrame({'idx': idx, 'value': value})


def idxmaxmin_combine(x, fn=None, skipna=True):
    if len(x) == 0:
        return x
    return (x.groupby(level=0)
             .apply(idxmaxmin_row, fn=fn, skipna=skipna)
             .reset_index(level=1, drop=True))


def idxmaxmin_agg(x, fn=None, skipna=True, scalar=False):
    res = idxmaxmin_combine(x, fn, skipna=skipna)['idx']
    if len(res) == 0:
        raise ValueError("attempt to get argmax of an empty sequence")
    if scalar:
        return res[0]
    res.name = None
    return res


def safe_head(df, n):
    r = df.head(n=n)
    if len(r) != n:
        msg = ("Insufficient elements for `head`. {0} elements "
               "requested, only {1} elements available. Try passing larger "
               "`npartitions` to `head`.")
        warnings.warn(msg.format(n, len(r)))
    return r


def maybe_shift_divisions(df, periods, freq):
    """Maybe shift divisions by periods of size freq

    Used to shift the divisions for the `shift` method. If freq isn't a fixed
    size (not anchored or relative), then the divisions are shifted
    appropriately. Otherwise the divisions are cleared.

    Parameters
    ----------
    df : dd.DataFrame, dd.Series, or dd.Index
    periods : int
        The number of periods to shift.
    freq : DateOffset, timedelta, or time rule string
        The frequency to shift by.
    """
    if isinstance(freq, str):
        freq = pd.tseries.frequencies.to_offset(freq)
    if (isinstance(freq, pd.DateOffset) and
            (freq.isAnchored() or not hasattr(freq, 'delta'))):
        # Can't infer divisions on relative or anchored offsets, as
        # divisions may now split identical index value.
        # (e.g. index_partitions = [[1, 2, 3], [3, 4, 5]])
        return df.clear_divisions()
    if df.known_divisions:
        divs = pd.Series(range(len(df.divisions)), index=df.divisions)
        divisions = divs.shift(periods, freq=freq).index
        return type(df)(df.dask, df._name, df._meta, divisions)
    return df


def to_delayed(df, optimize_graph=True):
    """Convert into a list of ``dask.delayed`` objects, one per partition.

    Deprecated, please use the equivalent ``df.to_delayed`` method instead.

    Parameters
    ----------
    optimize_graph : bool, optional
        If True [default], the graph is optimized before converting into
        ``dask.delayed`` objects.

    See Also
    --------
    dask.dataframe.from_delayed
    """
    warnings.warn("DeprecationWarning: The `dd.to_delayed` function is "
                  "deprecated, please use the `.to_delayed()` method instead.")
    return df.to_delayed(optimize_graph=optimize_graph)


@wraps(pd.to_datetime)
def to_datetime(arg, **kwargs):
    meta = pd.Series([pd.Timestamp('2000')])
    return map_partitions(pd.to_datetime, arg, meta=meta, **kwargs)


@wraps(pd.to_timedelta)
def to_timedelta(arg, unit='ns', errors='raise'):
    meta = pd.Series([pd.Timedelta(1, unit=unit)])
    return map_partitions(pd.to_timedelta, arg, unit=unit, errors=errors,
                          meta=meta)


if hasattr(pd, 'isna'):
    @wraps(pd.isna)
    def isna(arg):
        return map_partitions(pd.isna, arg)


def _repr_data_series(s, index):
    """A helper for creating the ``_repr_data`` property"""
    npartitions = len(index) - 1
    if is_categorical_dtype(s):
        if has_known_categories(s):
            dtype = 'category[known]'
        else:
            dtype = 'category[unknown]'
    else:
        dtype = str(s.dtype)
    return pd.Series([dtype] + ['...'] * npartitions, index=index, name=s.name)
