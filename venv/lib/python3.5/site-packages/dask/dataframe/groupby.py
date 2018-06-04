from __future__ import absolute_import, division, print_function

import collections
import itertools as it
import operator
import warnings

import numpy as np
import pandas as pd

from .core import (DataFrame, Series, aca, map_partitions, merge,
                   new_dd_object, no_default, split_out_on_index)
from .methods import drop_columns
from .shuffle import shuffle
from .utils import make_meta, insert_meta_param_description, raise_on_meta_error
from ..base import tokenize
from ..utils import derived_from, M, funcname, itemgetter


# #############################################
#
# GroupBy implementation notes
#
# Dask groupby supports reductions, i.e., mean, sum and alike, and apply. The
# former do not shuffle the data and are efficiently implemented as tree
# reductions. The latter is implemented by shuffling the underlying partiitons
# such that all items of a group can be found in the same parititon.
#
# The argument to ``.groupby``, the index, can be a ``str``, ``dd.DataFrame``,
# ``dd.Series``, or a list thereof. In operations on the grouped object, the
# divisions of the the grouped object and the items of index have to align.
# Currently, there is no support to shuffle the index values as part of the
# groupby operation. Therefore, the alignment has to be guaranteed by the
# caller.
#
# To operate on matchings paritions, most groupby operations exploit the
# corresponding support in ``apply_concat_apply``. Specifically, this function
# operates on matching paritiotns of frame-like objects passed as varargs.
#
# After the inital chunk step, the passed index is implicitly passed along to
# subsequent operations as the index of the parittions. Groupby operations on
# the individual parttions can then access the index via the ``levels``
# parameter of the ``groupby`` function. The correct arguments is determined by
# the ``_determine_levels`` function.
#
# To minimize overhead, series in an index that were obtained by getitem on the
# object to group are not passed as series to the various operations, but as
# columnn keys. This transformation is implemented as ``_normalize_index``.
#
# #############################################

def _determine_levels(index):
    """Determine the correct levels argument to groupby.
    """
    if isinstance(index, (tuple, list)) and len(index) > 1:
        return list(range(len(index)))
    else:
        return 0


def _normalize_index(df, index):
    """Replace series with column names in an index wherever possible.
    """
    if not isinstance(df, DataFrame):
        return index

    elif isinstance(index, list):
        return [_normalize_index(df, col) for col in index]

    elif (isinstance(index, Series) and index.name in df.columns and
          index._name == df[index.name]._name):
        return index.name

    elif (isinstance(index, DataFrame) and
          set(index.columns).issubset(df.columns) and
          index._name == df[index.columns]._name):
        return list(index.columns)

    else:
        return index


def _maybe_slice(grouped, columns):
    """
    Slice columns if grouped is pd.DataFrameGroupBy
    """
    if isinstance(grouped, pd.core.groupby.DataFrameGroupBy):
        if columns is not None:
            if isinstance(columns, (tuple, list, set, pd.Index)):
                columns = list(columns)
            return grouped[columns]
    return grouped


def _is_aligned(df, by):
    """Check if `df` and `by` have aligned indices"""
    if isinstance(by, (pd.Series, pd.DataFrame)):
        return df.index.equals(by.index)
    elif isinstance(by, (list, tuple)):
        return all(_is_aligned(df, i) for i in by)
    else:
        return True


def _groupby_raise_unaligned(df, **kwargs):
    """Groupby, but raise if df and `by` key are unaligned.

    Pandas supports grouping by a column that doesn't align with the input
    frame/series/index. However, the reindexing this causes doesn't seem to be
    threadsafe, and can result in incorrect results. Since grouping by an
    unaligned key is generally a bad idea, we just error loudly in dask.

    For more information see pandas GH issue #15244 and Dask GH issue #1876."""
    by = kwargs.get('by', None)
    if by is not None and not _is_aligned(df, by):
        msg = ("Grouping by an unaligned index is unsafe and unsupported.\n"
               "This can be caused by filtering only one of the object or\n"
               "grouping key. For example, the following works in pandas,\n"
               "but not in dask:\n"
               "\n"
               "df[df.foo < 0].groupby(df.bar)\n"
               "\n"
               "This can be avoided by either filtering beforehand, or\n"
               "passing in the name of the column instead:\n"
               "\n"
               "df2 = df[df.foo < 0]\n"
               "df2.groupby(df2.bar)\n"
               "# or\n"
               "df[df.foo < 0].groupby('bar')\n"
               "\n"
               "For more information see dask GH issue #1876.")
        raise ValueError(msg)
    elif by is not None and len(by):
        # since we're coming through apply, `by` will be a tuple.
        # Pandas treats tuples as a single key, and lists as multiple keys
        # We want multiple keys
        kwargs.update(by=list(by))
    return df.groupby(**kwargs)


def _groupby_slice_apply(df, grouper, key, func, *args, **kwargs):
    # No need to use raise if unaligned here - this is only called after
    # shuffling, which makes everything aligned already
    g = df.groupby(grouper)
    if key:
        g = g[key]
    return g.apply(func, *args, **kwargs)


def _groupby_get_group(df, by_key, get_key, columns):
    # SeriesGroupBy may pass df which includes group key
    grouped = _groupby_raise_unaligned(df, by=by_key)

    if get_key in grouped.groups:
        if isinstance(df, pd.DataFrame):
            grouped = grouped[columns]
        return grouped.get_group(get_key)

    else:
        # to create empty DataFrame/Series, which has the same
        # dtype as the original
        if isinstance(df, pd.DataFrame):
            # may be SeriesGroupBy
            df = df[columns]
        return df.iloc[0:0]


###############################################################
# Aggregation
###############################################################

# Implementation detail: use class to make it easier to pass inside spec
class Aggregation(object):
    """A user defined aggregation.

    Parameters
    ----------
    name : str
        the name of the aggregation. It should be unique, since intermediate
        result will be identified by this name.
    chunk : callable
        a function that will be called with the grouped column of each
        partition. It can either return a single series or a tuple of series.
        The index has to be equal to the groups.
    agg : callable
        a function that will be called to aggregate the results of each chunk.
        Again the argument(s) will be grouped series. If ``chunk`` returned a
        tuple, ``agg`` will be called with all of them as individual positional
        arguments.
    finalize : callable
        an optional finalizer that will be called with the results from the
        aggregation.

    Examples
    --------

    ``sum`` can be implemented as::

        custom_sum = dd.Aggregation('custom_sum', lambda s: s.sum(), lambda s0: s0.sum())
        df.groupby('g').agg(custom_sum)

    and ``mean`` can be implemented as::

        custom_mean = dd.Aggregation(
            'custom_mean',
            lambda s: (s.count(), s.sum()),
            lambda count, sum: (count.sum(), sum.sum()),
            lambda count, sum: sum / count,
        )
        df.groupby('g').agg(custom_mean)

    """
    def __init__(self, name, chunk, agg, finalize=None):
        self.chunk = chunk
        self.agg = agg
        self.finalize = finalize
        self.__name__ = name


def _groupby_aggregate(df, aggfunc=None, levels=None):
    return aggfunc(df.groupby(level=levels, sort=False))


def _apply_chunk(df, *index, **kwargs):
    func = kwargs.pop('chunk')
    columns = kwargs.pop('columns')

    g = _groupby_raise_unaligned(df, by=index)

    if isinstance(df, pd.Series) or columns is None:
        return func(g)
    else:
        if isinstance(columns, (tuple, list, set, pd.Index)):
            columns = list(columns)
        return func(g[columns])


def _var_chunk(df, *index):
    if isinstance(df, pd.Series):
        df = df.to_frame()
    g = _groupby_raise_unaligned(df, by=index)
    x = g.sum()
    x2 = g.agg(lambda x: (x**2).sum()).rename(columns=lambda c: c + '-x2')
    n = g.count().rename(columns=lambda c: c + '-count')
    return pd.concat([x, x2, n], axis=1)


def _var_combine(g, levels):
    return g.groupby(level=levels, sort=False).sum()


def _var_agg(g, levels, ddof):
    g = g.groupby(level=levels, sort=False).sum()
    nc = len(g.columns)
    x = g[g.columns[:nc // 3]]
    x2 = g[g.columns[nc // 3:2 * nc // 3]].rename(columns=lambda c: c[:-3])
    n = g[g.columns[-nc // 3:]].rename(columns=lambda c: c[:-6])

    # TODO: replace with _finalize_var?
    result = x2 - x ** 2 / n
    div = (n - ddof)
    div[div < 0] = 0
    result /= div
    result[(n - ddof) == 0] = np.nan
    assert isinstance(result, pd.DataFrame)
    return result


###############################################################
# nunique
###############################################################

def _nunique_df_chunk(df, *index, **kwargs):
    levels = kwargs.pop('levels')
    name = kwargs.pop('name')

    g = _groupby_raise_unaligned(df, by=index)
    grouped = g[[name]].apply(pd.DataFrame.drop_duplicates)

    # we set the index here to force a possibly duplicate index
    # for our reduce step
    if isinstance(levels, list):
        grouped.index = pd.MultiIndex.from_arrays([
            grouped.index.get_level_values(level=level) for level in levels
        ])
    else:
        grouped.index = grouped.index.get_level_values(level=levels)

    return grouped


def _drop_duplicates_rename(df):
    # Avoid duplicate index labels in a groupby().apply() context
    # https://github.com/dask/dask/issues/3039
    # https://github.com/pandas-dev/pandas/pull/18882
    names = [None] * df.index.nlevels
    return df.drop_duplicates().rename_axis(names, copy=False)


def _nunique_df_combine(df, levels):
    result = df.groupby(level=levels,
                        sort=False).apply(_drop_duplicates_rename)

    if isinstance(levels, list):
        result.index = pd.MultiIndex.from_arrays([
            result.index.get_level_values(level=level) for level in levels
        ])
    else:
        result.index = result.index.get_level_values(level=levels)

    return result


def _nunique_df_aggregate(df, levels, name):
    return df.groupby(level=levels, sort=False)[name].nunique()


def _nunique_series_chunk(df, *index, **_ignored_):
    # convert series to data frame, then hand over to dataframe code path
    assert isinstance(df, pd.Series)

    df = df.to_frame()
    kwargs = dict(name=df.columns[0], levels=_determine_levels(index))
    return _nunique_df_chunk(df, *index, **kwargs)


###############################################################
# Aggregate support
#
# Aggregate is implemented as:
#
# 1. group-by-aggregate all partitions into intermediate values
# 2. collect all partitions into a single partition
# 3. group-by-aggregate the result into intermediate values
# 4. transform all intermediate values into the result
#
# In Step 1 and 3 the dataframe is grouped on the same columns.
#
###############################################################
def _make_agg_id(func, column):
    return '{!s}-{!s}-{}'.format(func, column, tokenize(func, column))


def _normalize_spec(spec, non_group_columns):
    """
    Return a list of ``(result_column, func, input_column)`` tuples.

    Spec can be

    - a function
    - a list of functions
    - a dictionary that maps input-columns to functions
    - a dictionary that maps input-columns to a lists of functions
    - a dictionary that maps input-columns to a dictionaries that map
      output-columns to functions.

    The non-group columns are a list of all column names that are not used in
    the groupby operation.

    Usually, the result columns are mutli-level names, returned as tuples.
    If only a single function is supplied or dictionary mapping columns
    to single functions, simple names are returned as strings (see the first
    two examples below).

    Examples
    --------
    >>> _normalize_spec('mean', ['a', 'b', 'c'])
    [('a', 'mean', 'a'), ('b', 'mean', 'b'), ('c', 'mean', 'c')]

    >>> spec = collections.OrderedDict([('a', 'mean'), ('b', 'count')])
    >>> _normalize_spec(spec, ['a', 'b', 'c'])
    [('a', 'mean', 'a'), ('b', 'count', 'b')]

    >>> _normalize_spec(['var', 'mean'], ['a', 'b', 'c'])
    ... # doctest: +NORMALIZE_WHITESPACE
    [(('a', 'var'), 'var', 'a'), (('a', 'mean'), 'mean', 'a'), \
     (('b', 'var'), 'var', 'b'), (('b', 'mean'), 'mean', 'b'), \
     (('c', 'var'), 'var', 'c'), (('c', 'mean'), 'mean', 'c')]

    >>> spec = collections.OrderedDict([('a', 'mean'), ('b', ['sum', 'count'])])
    >>> _normalize_spec(spec, ['a', 'b', 'c'])
    ... # doctest: +NORMALIZE_WHITESPACE
    [(('a', 'mean'), 'mean', 'a'), (('b', 'sum'), 'sum', 'b'), \
      (('b', 'count'), 'count', 'b')]

    >>> spec = collections.OrderedDict()
    >>> spec['a'] = ['mean', 'size']
    >>> spec['b'] = collections.OrderedDict([('e', 'count'), ('f', 'var')])
    >>> _normalize_spec(spec, ['a', 'b', 'c'])
    ... # doctest: +NORMALIZE_WHITESPACE
    [(('a', 'mean'), 'mean', 'a'), (('a', 'size'), 'size', 'a'), \
     (('b', 'e'), 'count', 'b'), (('b', 'f'), 'var', 'b')]
    """
    if not isinstance(spec, dict):
        spec = collections.OrderedDict(zip(non_group_columns, it.repeat(spec)))

    res = []

    if isinstance(spec, dict):
        for input_column, subspec in spec.items():
            if isinstance(subspec, dict):
                res.extend(((input_column, result_column), func, input_column)
                           for result_column, func in subspec.items())

            else:
                if not isinstance(subspec, list):
                    subspec = [subspec]

                res.extend(((input_column, funcname(func)), func, input_column)
                           for func in subspec)

    else:
        raise ValueError("unsupported agg spec of type {}".format(type(spec)))

    compounds = (list, tuple, dict)
    use_flat_columns = not any(isinstance(subspec, compounds)
                               for subspec in spec.values())

    if use_flat_columns:
        res = [(input_col, func, input_col) for (_, func, input_col) in res]

    return res


def _build_agg_args(spec):
    """
    Create transformation functions for a normalized aggregate spec.

    Parameters
    ----------
    spec: a list of (result-column, aggregation-function, input-column) triples.
        To work with all arugment forms understood by pandas use
        ``_normalize_spec`` to normalize the argment before passing it on to
        ``_build_agg_args``.

    Returns
    -------
    chunk_funcs: a list of (intermediate-column, function, keyword) triples
        that are applied on grouped chunks of the initial dataframe.

    agg_funcs: a list of (intermediate-column, functions, keword) triples that
        are applied on the grouped concatination of the preprocessed chunks.

    finalizers: a list of (result-column, function, keyword) triples that are
        applied after the ``agg_funcs``. They are used to create final results
        from intermediate representations.
    """
    known_np_funcs = {np.min: 'min', np.max: 'max'}

    # check that there are no name conflicts for a single input column
    by_name = {}
    for _, func, input_column in spec:
        key = funcname(known_np_funcs.get(func, func)), input_column
        by_name.setdefault(key, []).append((func, input_column))

    for funcs in by_name.values():
        if len(funcs) != 1:
            raise ValueError('conflicting aggregation functions: {}'.format(funcs))

    chunks = {}
    aggs = {}
    finalizers = []

    for (result_column, func, input_column) in spec:
        if not isinstance(func, Aggregation):
            func = funcname(known_np_funcs.get(func, func))

        impls = _build_agg_args_single(result_column, func, input_column)

        # overwrite existing result-columns, generate intermediates only once
        chunks.update((spec[0], spec) for spec in impls['chunk_funcs'])
        aggs.update((spec[0], spec) for spec in impls['aggregate_funcs'])

        finalizers.append(impls['finalizer'])

    chunks = sorted(chunks.values())
    aggs = sorted(aggs.values())

    return chunks, aggs, finalizers


def _build_agg_args_single(result_column, func, input_column):
    simple_impl = {
        'sum': (M.sum, M.sum),
        'min': (M.min, M.min),
        'max': (M.max, M.max),
        'count': (M.count, M.sum),
        'size': (M.size, M.sum),
        'first': (M.first, M.first),
        'last': (M.last, M.last)
    }

    if func in simple_impl.keys():
        return _build_agg_args_simple(result_column, func, input_column,
                                      simple_impl[func])

    elif func == 'var':
        return _build_agg_args_var(result_column, func, input_column)

    elif func == 'std':
        return _build_agg_args_std(result_column, func, input_column)

    elif func == 'mean':
        return _build_agg_args_mean(result_column, func, input_column)

    elif isinstance(func, Aggregation):
        return _build_agg_args_custom(result_column, func, input_column)

    else:
        raise ValueError("unknown aggregate {}".format(func))


def _build_agg_args_simple(result_column, func, input_column, impl_pair):
    intermediate = _make_agg_id(func, input_column)
    chunk_impl, agg_impl = impl_pair

    return dict(
        chunk_funcs=[(intermediate, _apply_func_to_column,
                     dict(column=input_column, func=chunk_impl))],
        aggregate_funcs=[(intermediate, _apply_func_to_column,
                         dict(column=intermediate, func=agg_impl))],
        finalizer=(result_column, itemgetter(intermediate), dict()),
    )


def _build_agg_args_var(result_column, func, input_column):
    int_sum = _make_agg_id('sum', input_column)
    int_sum2 = _make_agg_id('sum2', input_column)
    int_count = _make_agg_id('count', input_column)

    return dict(
        chunk_funcs=[
            (int_sum, _apply_func_to_column,
             dict(column=input_column, func=M.sum)),
            (int_count, _apply_func_to_column,
             dict(column=input_column, func=M.count)),
            (int_sum2, _compute_sum_of_squares,
             dict(column=input_column)),
        ],
        aggregate_funcs=[
            (col, _apply_func_to_column, dict(column=col, func=M.sum))
            for col in (int_sum, int_count, int_sum2)
        ],
        finalizer=(result_column, _finalize_var,
                   dict(sum_column=int_sum, count_column=int_count,
                        sum2_column=int_sum2)),
    )


def _build_agg_args_std(result_column, func, input_column):
    impls = _build_agg_args_var(result_column, func, input_column)

    result_column, _, kwargs = impls['finalizer']
    impls['finalizer'] = (result_column, _finalize_std, kwargs)

    return impls


def _build_agg_args_mean(result_column, func, input_column):
    int_sum = _make_agg_id('sum', input_column)
    int_count = _make_agg_id('count', input_column)

    return dict(
        chunk_funcs=[
            (int_sum, _apply_func_to_column,
             dict(column=input_column, func=M.sum)),
            (int_count, _apply_func_to_column,
             dict(column=input_column, func=M.count)),
        ],
        aggregate_funcs=[
            (col, _apply_func_to_column, dict(column=col, func=M.sum))
            for col in (int_sum, int_count)
        ],
        finalizer=(result_column, _finalize_mean,
                   dict(sum_column=int_sum, count_column=int_count)),
    )


def _build_agg_args_custom(result_column, func, input_column):
    col = _make_agg_id(funcname(func), input_column)

    if func.finalize is None:
        finalizer = (result_column, operator.itemgetter(col), dict())

    else:
        finalizer = (
            result_column, _apply_func_to_columns,
            dict(func=func.finalize, prefix=col)
        )

    return dict(
        chunk_funcs=[
            (col, _apply_func_to_column,
             dict(func=func.chunk, column=input_column))
        ],
        aggregate_funcs=[
            (col, _apply_func_to_columns, dict(func=func.agg, prefix=col))
        ],
        finalizer=finalizer
    )


def _groupby_apply_funcs(df, *index, **kwargs):
    """
    Group a dataframe and apply multiple aggregation functions.

    Parameters
    ----------
    df: pandas.DataFrame
        The dataframe to work on.
    index: list of groupers
        If given, they are added to the keyword arguments as the ``by``
        argument.
    funcs: list of result-colum, function, keywordargument triples
        The list of functions that are applied on the grouped data frame.
        Has to be passed as a keyword argument.
    kwargs:
        All keyword arguments, but ``funcs``, are passed verbatim to the groupby
        operation of the dataframe

    Returns
    -------
    aggregated:
        the aggregated dataframe.
    """
    if len(index):
        # since we're coming through apply, `by` will be a tuple.
        # Pandas treats tuples as a single key, and lists as multiple keys
        # We want multiple keys
        kwargs.update(by=list(index))

    funcs = kwargs.pop('funcs')
    grouped = _groupby_raise_unaligned(df, **kwargs)

    result = collections.OrderedDict()
    for result_column, func, func_kwargs in funcs:
        r = func(grouped, **func_kwargs)

        if isinstance(r, tuple):
            for idx, s in enumerate(r):
                result['{}-{}'.format(result_column, idx)] = s

        else:
            result[result_column] = r

    return pd.DataFrame(result)


def _compute_sum_of_squares(grouped, column):
    base = grouped[column] if column is not None else grouped
    return base.apply(lambda x: (x ** 2).sum())


def _agg_finalize(df, aggregate_funcs, finalize_funcs, level):
    # finish the final aggregation level
    df = _groupby_apply_funcs(df, funcs=aggregate_funcs, level=level)

    # and finalize the result
    result = collections.OrderedDict()
    for result_column, func, kwargs in finalize_funcs:
        result[result_column] = func(df, **kwargs)

    return pd.DataFrame(result)


def _apply_func_to_column(df_like, column, func):
    if column is None:
        return func(df_like)

    return func(df_like[column])


def _apply_func_to_columns(df_like, prefix, func):
    if isinstance(df_like, pd.DataFrame):
        columns = df_like.columns
    else:
        # handle GroupBy objects
        columns = df_like._selected_obj.columns

    columns = sorted(col for col in columns if col.startswith(prefix))

    columns = [df_like[col] for col in columns]
    return func(*columns)


def _finalize_mean(df, sum_column, count_column):
    return df[sum_column] / df[count_column]


def _finalize_var(df, count_column, sum_column, sum2_column, ddof=1):
    n = df[count_column]
    x = df[sum_column]
    x2 = df[sum2_column]

    result = x2 - x ** 2 / n
    div = (n - ddof)
    div[div < 0] = 0
    result /= div
    result[(n - ddof) == 0] = np.nan

    return result


def _finalize_std(df, count_column, sum_column, sum2_column, ddof=1):
    result = _finalize_var(df, count_column, sum_column, sum2_column, ddof)
    return np.sqrt(result)


def _cum_agg_aligned(part, cum_last, index, columns, func, initial):
    align = cum_last.reindex(part.set_index(index).index, fill_value=initial)
    align.index = part.index
    return func(part[columns], align)


def _cum_agg_filled(a, b, func, initial):
    union = a.index.union(b.index)
    return func(a.reindex(union, fill_value=initial),
                b.reindex(union, fill_value=initial), fill_value=initial)


def _cumcount_aggregate(a, b, fill_value=None):
    return a.add(b, fill_value=fill_value) + 1


class _GroupBy(object):
    """ Superclass for DataFrameGroupBy and SeriesGroupBy

    Parameters
    ----------

    obj: DataFrame or Series
        DataFrame or Series to be grouped
    by: str, list or Series
        The key for grouping
    slice: str, list
        The slice keys applied to GroupBy result
    """
    def __init__(self, df, by=None, slice=None):

        assert isinstance(df, (DataFrame, Series))
        self.obj = df
        # grouping key passed via groupby method
        self.index = _normalize_index(df, by)

        if isinstance(self.index, list):
            do_index_partition_align = all(
                item.divisions == df.divisions if isinstance(item, Series) else True
                for item in self.index
            )
        elif isinstance(self.index, Series):
            do_index_partition_align = df.divisions == self.index.divisions
        else:
            do_index_partition_align = True

        if not do_index_partition_align:
            raise NotImplementedError("The grouped object and index of the "
                                      "groupby must have the same divisions.")

        # slicing key applied to _GroupBy instance
        self._slice = slice

        if isinstance(self.index, list):
            index_meta = [item._meta if isinstance(item, Series) else item for item in self.index]

        elif isinstance(self.index, Series):
            index_meta = self.index._meta

        else:
            index_meta = self.index

        self._meta = self.obj._meta.groupby(index_meta)

    @property
    def _meta_nonempty(self):
        """
        Return a pd.DataFrameGroupBy / pd.SeriesGroupBy which contains sample data.
        """
        sample = self.obj._meta_nonempty

        if isinstance(self.index, list):
            index_meta = [item._meta_nonempty if isinstance(item, Series) else item for item in self.index]

        elif isinstance(self.index, Series):
            index_meta = self.index._meta_nonempty

        else:
            index_meta = self.index

        grouped = sample.groupby(index_meta)
        return _maybe_slice(grouped, self._slice)

    def _aca_agg(self, token, func, aggfunc=None, split_every=None,
                 split_out=1):
        if aggfunc is None:
            aggfunc = func

        meta = func(self._meta)
        columns = meta.name if isinstance(meta, pd.Series) else meta.columns

        token = self._token_prefix + token
        levels = _determine_levels(self.index)

        return aca([self.obj, self.index] if not isinstance(self.index, list) else [self.obj] + self.index,
                   chunk=_apply_chunk,
                   chunk_kwargs=dict(chunk=func, columns=columns),
                   aggregate=_groupby_aggregate,
                   meta=meta, token=token, split_every=split_every,
                   aggregate_kwargs=dict(aggfunc=aggfunc, levels=levels),
                   split_out=split_out, split_out_setup=split_out_on_index)

    def _cum_agg(self, token, chunk, aggregate, initial):
        """ Wrapper for cumulative groupby operation """
        meta = chunk(self._meta)
        columns = meta.name if isinstance(meta, pd.Series) else meta.columns
        index = self.index if isinstance(self.index, list) else [self.index]

        name = self._token_prefix + token
        name_part = name + '-map'
        name_last = name + '-take-last'
        name_cum = name + '-cum-last'

        # cumulate each partitions
        cumpart_raw = map_partitions(_apply_chunk, self.obj, *index,
                                     chunk=chunk,
                                     columns=columns,
                                     token=name_part,
                                     meta=meta)

        cumpart_raw_frame = (cumpart_raw.to_frame()
                             if isinstance(meta, pd.Series)
                             else cumpart_raw)

        cumpart_ext = cumpart_raw_frame.assign(
            **{i: self.obj[i]
               if np.isscalar(i) and i in self.obj.columns
               else self.obj.index
               for i in index})

        # Use pd.Grouper objects to specify that we are grouping by columns.
        # Otherwise, pandas will throw an ambiguity warning if the
        # DataFrame's index (self.obj.index) was included in the grouping
        # specification (self.index). See pandas #14432
        index_groupers = [pd.Grouper(key=ind) for ind in index]
        cumlast = map_partitions(_apply_chunk, cumpart_ext, *index_groupers,
                                 columns=0 if columns is None else columns,
                                 chunk=M.last,
                                 meta=meta,
                                 token=name_last)

        # aggregate cumulated partisions and its previous last element
        dask = {}
        dask[(name, 0)] = (cumpart_raw._name, 0)

        for i in range(1, self.obj.npartitions):
            # store each cumulative step to graph to reduce computation
            if i == 1:
                dask[(name_cum, i)] = (cumlast._name, i - 1)
            else:
                # aggregate with previous cumulation results
                dask[(name_cum, i)] = (_cum_agg_filled,
                                       (name_cum, i - 1),
                                       (cumlast._name, i - 1),
                                       aggregate, initial)
            dask[(name, i)] = (_cum_agg_aligned,
                               (cumpart_ext._name, i), (name_cum, i),
                               index, 0 if columns is None else columns,
                               aggregate, initial)
        return new_dd_object(merge(dask, cumpart_ext.dask, cumlast.dask),
                             name, chunk(self._meta), self.obj.divisions)

    @derived_from(pd.core.groupby.GroupBy)
    def cumsum(self, axis=0):
        if axis:
            return self.obj.cumsum(axis=axis)
        else:
            return self._cum_agg('cumsum',
                                 chunk=M.cumsum,
                                 aggregate=M.add,
                                 initial=0)

    @derived_from(pd.core.groupby.GroupBy)
    def cumprod(self, axis=0):
        if axis:
            return self.obj.cumprod(axis=axis)
        else:
            return self._cum_agg('cumprod',
                                 chunk=M.cumprod,
                                 aggregate=M.mul,
                                 initial=1)

    @derived_from(pd.core.groupby.GroupBy)
    def cumcount(self, axis=None):
        return self._cum_agg('cumcount',
                             chunk=M.cumcount,
                             aggregate=_cumcount_aggregate,
                             initial=-1)

    @derived_from(pd.core.groupby.GroupBy)
    def sum(self, split_every=None, split_out=1):
        return self._aca_agg(token='sum', func=M.sum, split_every=split_every,
                             split_out=split_out)

    @derived_from(pd.core.groupby.GroupBy)
    def min(self, split_every=None, split_out=1):
        return self._aca_agg(token='min', func=M.min, split_every=split_every,
                             split_out=split_out)

    @derived_from(pd.core.groupby.GroupBy)
    def max(self, split_every=None, split_out=1):
        return self._aca_agg(token='max', func=M.max, split_every=split_every,
                             split_out=split_out)

    @derived_from(pd.core.groupby.GroupBy)
    def count(self, split_every=None, split_out=1):
        return self._aca_agg(token='count', func=M.count,
                             aggfunc=M.sum, split_every=split_every,
                             split_out=split_out)

    @derived_from(pd.core.groupby.GroupBy)
    def mean(self, split_every=None, split_out=1):
        return (self.sum(split_every=split_every, split_out=split_out) /
                self.count(split_every=split_every, split_out=split_out))

    @derived_from(pd.core.groupby.GroupBy)
    def size(self, split_every=None, split_out=1):
        return self._aca_agg(token='size', func=M.size, aggfunc=M.sum,
                             split_every=split_every, split_out=split_out)

    @derived_from(pd.core.groupby.GroupBy)
    def var(self, ddof=1, split_every=None, split_out=1):
        levels = _determine_levels(self.index)
        result = aca([self.obj, self.index] if not isinstance(self.index, list) else [self.obj] + self.index,
                     chunk=_var_chunk,
                     aggregate=_var_agg, combine=_var_combine,
                     token=self._token_prefix + 'var',
                     aggregate_kwargs={'ddof': ddof, 'levels': levels},
                     combine_kwargs={'levels': levels},
                     split_every=split_every, split_out=split_out,
                     split_out_setup=split_out_on_index)

        if isinstance(self.obj, Series):
            result = result[result.columns[0]]
        if self._slice:
            result = result[self._slice]

        return result

    @derived_from(pd.core.groupby.GroupBy)
    def std(self, ddof=1, split_every=None, split_out=1):
        v = self.var(ddof, split_every=split_every, split_out=split_out)
        result = map_partitions(np.sqrt, v, meta=v)
        return result

    @derived_from(pd.core.groupby.GroupBy)
    def first(self, split_every=None, split_out=1):
        return self._aca_agg(token='first', func=M.first, split_every=split_every,
                             split_out=split_out)

    @derived_from(pd.core.groupby.GroupBy)
    def last(self, split_every=None, split_out=1):
        return self._aca_agg(token='last', func=M.last, split_every=split_every,
                             split_out=split_out)

    @derived_from(pd.core.groupby.GroupBy)
    def get_group(self, key):
        token = self._token_prefix + 'get_group'

        meta = self._meta.obj
        if isinstance(meta, pd.DataFrame) and self._slice is not None:
            meta = meta[self._slice]
        columns = meta.columns if isinstance(meta, pd.DataFrame) else meta.name

        return map_partitions(_groupby_get_group, self.obj, self.index, key,
                              columns, meta=meta, token=token)

    def aggregate(self, arg, split_every, split_out=1):
        if isinstance(self.obj, DataFrame):
            if isinstance(self.index, tuple) or np.isscalar(self.index):
                group_columns = {self.index}

            elif isinstance(self.index, list):
                group_columns = {i for i in self.index
                                 if isinstance(i, tuple) or np.isscalar(i)}

            else:
                group_columns = set()

            if self._slice:
                # pandas doesn't exclude the grouping column in a SeriesGroupBy
                # like df.groupby('a')['a'].agg(...)
                non_group_columns = self._slice
                if not isinstance(non_group_columns, list):
                    non_group_columns = [non_group_columns]
            else:
                # NOTE: this step relies on the index normalization to replace
                #       series with their name in an index.
                non_group_columns = [col for col in self.obj.columns
                                     if col not in group_columns]

            spec = _normalize_spec(arg, non_group_columns)

        elif isinstance(self.obj, Series):
            if isinstance(arg, (list, tuple, dict)):
                # implementation detail: if self.obj is a series, a pseudo column
                # None is used to denote the series itself. This pseudo column is
                # removed from the result columns before passing the spec along.
                spec = _normalize_spec({None: arg}, [])
                spec = [(result_column, func, input_column)
                        for ((_, result_column), func, input_column) in spec]

            else:
                spec = _normalize_spec({None: arg}, [])
                spec = [(self.obj.name, func, input_column)
                        for (_, func, input_column) in spec]

        else:
            raise ValueError("aggregate on unknown object {}".format(self.obj))

        chunk_funcs, aggregate_funcs, finalizers = _build_agg_args(spec)

        if isinstance(self.index, (tuple, list)) and len(self.index) > 1:
            levels = list(range(len(self.index)))
        else:
            levels = 0

        if not isinstance(self.index, list):
            chunk_args = [self.obj, self.index]

        else:
            chunk_args = [self.obj] + self.index

        return aca(chunk_args,
                   chunk=_groupby_apply_funcs,
                   chunk_kwargs=dict(funcs=chunk_funcs),
                   combine=_groupby_apply_funcs,
                   combine_kwargs=dict(funcs=aggregate_funcs, level=levels),
                   aggregate=_agg_finalize,
                   aggregate_kwargs=dict(
                       aggregate_funcs=aggregate_funcs,
                       finalize_funcs=finalizers,
                       level=levels,
                   ),
                   token='aggregate', split_every=split_every,
                   split_out=split_out, split_out_setup=split_out_on_index)

    @insert_meta_param_description(pad=12)
    def apply(self, func, *args, **kwargs):
        """ Parallel version of pandas GroupBy.apply

        This mimics the pandas version except for the following:

        1.  The user should provide output metadata.
        2.  If the grouper does not align with the index then this causes a full
            shuffle.  The order of rows within each group may not be preserved.

        Parameters
        ----------
        func: function
            Function to apply
        args, kwargs : Scalar, Delayed or object
            Arguments and keywords to pass to the function.
        $META

        Returns
        -------
        applied : Series or DataFrame depending on columns keyword
        """
        meta = kwargs.get('meta', no_default)

        if meta is no_default:
            msg = ("`meta` is not specified, inferred from partial data. "
                   "Please provide `meta` if the result is unexpected.\n"
                   "  Before: .apply(func)\n"
                   "  After:  .apply(func, meta={'x': 'f8', 'y': 'f8'}) for dataframe result\n"
                   "  or:     .apply(func, meta=('x', 'f8'))            for series result")
            warnings.warn(msg, stacklevel=2)

            with raise_on_meta_error("groupby.apply({0})".format(funcname(func))):
                meta = self._meta_nonempty.apply(func, *args, **kwargs)

        meta = make_meta(meta)

        # Validate self.index
        if (isinstance(self.index, list) and
                any(isinstance(item, Series) for item in self.index)):
            raise NotImplementedError("groupby-apply with a multiple Series "
                                      "is currently not supported")

        df = self.obj
        should_shuffle = not (df.known_divisions and
                              df._contains_index_name(self.index))

        if should_shuffle:
            if isinstance(self.index, DataFrame):  # add index columns to dataframe
                df2 = df.assign(**{'_index_' + c: self.index[c]
                                   for c in self.index.columns})
                index = self.index
            elif isinstance(self.index, Series):
                df2 = df.assign(_index=self.index)
                index = self.index
            else:
                df2 = df
                index = df._select_columns_or_index(self.index)

            df3 = shuffle(df2, index)  # shuffle dataframe and index
        else:
            df3 = df

        if should_shuffle and isinstance(self.index, DataFrame):
            # extract index from dataframe
            cols = ['_index_' + c for c in self.index.columns]
            index2 = df3[cols]
            if isinstance(meta, pd.DataFrame):
                df4 = df3.map_partitions(drop_columns, cols, meta.columns.dtype)
            else:
                df4 = df3.drop(cols, axis=1)
        elif should_shuffle and isinstance(self.index, Series):
            index2 = df3['_index']
            index2.name = self.index.name
            if isinstance(meta, pd.DataFrame):
                df4 = df3.map_partitions(drop_columns, '_index',
                                         meta.columns.dtype)
            else:
                df4 = df3.drop('_index', axis=1)
        else:
            df4 = df3
            index2 = self.index

        # Perform embarrassingly parallel groupby-apply
        kwargs['meta'] = meta
        df5 = map_partitions(_groupby_slice_apply, df4, index2,
                             self._slice, func, *args, **kwargs)

        return df5


class DataFrameGroupBy(_GroupBy):

    _token_prefix = 'dataframe-groupby-'

    def __getitem__(self, key):
        if isinstance(key, list):
            g = DataFrameGroupBy(self.obj, by=self.index, slice=key)
        else:
            g = SeriesGroupBy(self.obj, by=self.index, slice=key)

        # error is raised from pandas
        g._meta = g._meta[key]
        return g

    def __dir__(self):
        return sorted(set(dir(type(self)) + list(self.__dict__) +
                      list(filter(pd.compat.isidentifier, self.obj.columns))))

    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError as e:
            raise AttributeError(e)

    @derived_from(pd.core.groupby.DataFrameGroupBy)
    def aggregate(self, arg, split_every=None, split_out=1):
        if arg == 'size':
            return self.size()

        return super(DataFrameGroupBy, self).aggregate(arg, split_every=split_every, split_out=split_out)

    @derived_from(pd.core.groupby.DataFrameGroupBy)
    def agg(self, arg, split_every=None, split_out=1):
        return self.aggregate(arg, split_every=split_every, split_out=split_out)


class SeriesGroupBy(_GroupBy):

    _token_prefix = 'series-groupby-'

    def __init__(self, df, by=None, slice=None):
        # for any non series object, raise pandas-compat error message

        if isinstance(df, Series):
            if isinstance(by, Series):
                pass
            elif isinstance(by, list):
                if len(by) == 0:
                    raise ValueError("No group keys passed!")

                non_series_items = [item for item in by
                                    if not isinstance(item, Series)]
                # raise error from pandas, if applicable
                df._meta.groupby(non_series_items)
            else:
                # raise error from pandas, if applicable
                df._meta.groupby(by)

        super(SeriesGroupBy, self).__init__(df, by=by, slice=slice)

    def nunique(self, split_every=None, split_out=1):
        name = self._meta.obj.name
        levels = _determine_levels(self.index)

        if isinstance(self.obj, DataFrame):
            chunk = _nunique_df_chunk

        else:
            chunk = _nunique_series_chunk

        return aca([self.obj, self.index] if not isinstance(self.index, list) else [self.obj] + self.index,
                   chunk=chunk,
                   aggregate=_nunique_df_aggregate,
                   combine=_nunique_df_combine,
                   token='series-groupby-nunique',
                   chunk_kwargs={'levels': levels, 'name': name},
                   aggregate_kwargs={'levels': levels, 'name': name},
                   combine_kwargs={'levels': levels},
                   split_every=split_every, split_out=split_out,
                   split_out_setup=split_out_on_index)

    @derived_from(pd.core.groupby.SeriesGroupBy)
    def aggregate(self, arg, split_every=None, split_out=1):
        result = super(SeriesGroupBy, self).aggregate(arg, split_every=split_every, split_out=split_out)
        if self._slice:
            result = result[self._slice]

        if not isinstance(arg, (list, dict)) and isinstance(result, DataFrame):
            result = result[result.columns[0]]

        return result

    @derived_from(pd.core.groupby.SeriesGroupBy)
    def agg(self, arg, split_every=None, split_out=1):
        return self.aggregate(arg, split_every=split_every, split_out=split_out)
