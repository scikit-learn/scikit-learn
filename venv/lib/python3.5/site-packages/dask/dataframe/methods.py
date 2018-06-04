from __future__ import print_function, absolute_import, division

import warnings

import numpy as np
import pandas as pd
from pandas.api.types import is_categorical_dtype
from toolz import partition

from .utils import PANDAS_VERSION
if PANDAS_VERSION >= '0.20.0':
    from pandas.api.types import union_categoricals
else:
    from pandas.types.concat import union_categoricals

# ---------------------------------
# indexing
# ---------------------------------


def loc(df, iindexer, cindexer=None):
    """
    .loc for known divisions
    """
    if cindexer is None:
        return df.loc[iindexer]
    else:
        return df.loc[iindexer, cindexer]


def try_loc(df, iindexer, cindexer=None):
    """
    .loc for unknown divisions
    """
    try:
        return loc(df, iindexer, cindexer)
    except KeyError:
        return df.head(0).loc[:, cindexer]


def boundary_slice(df, start, stop, right_boundary=True, left_boundary=True,
                   kind='loc'):
    """Index slice start/stop. Can switch include/exclude boundaries.

    >>> df = pd.DataFrame({'x': [10, 20, 30, 40, 50]}, index=[1, 2, 2, 3, 4])
    >>> boundary_slice(df, 2, None)
        x
    2  20
    2  30
    3  40
    4  50
    >>> boundary_slice(df, 1, 3)
        x
    1  10
    2  20
    2  30
    3  40
    >>> boundary_slice(df, 1, 3, right_boundary=False)
        x
    1  10
    2  20
    2  30
    """
    if kind == 'loc' and not df.index.is_monotonic:
        # Pandas treats missing keys differently for label-slicing
        # on monotonic vs. non-monotonic indexes
        # If the index is monotonic, `df.loc[start:stop]` is fine.
        # If it's not, `df.loc[start:stop]` raises when `start` is missing
        if start is not None:
            if left_boundary:
                df = df[df.index >= start]
            else:
                df = df[df.index > start]
        if stop is not None:
            if right_boundary:
                df = df[df.index <= stop]
            else:
                df = df[df.index < stop]
        return df
    else:
        result = getattr(df, kind)[start:stop]
    if not right_boundary:
        right_index = result.index.get_slice_bound(stop, 'left', kind)
        result = result.iloc[:right_index]
    if not left_boundary:
        left_index = result.index.get_slice_bound(start, 'right', kind)
        result = result.iloc[left_index:]
    return result


def index_count(x):
    # Workaround since Index doesn't implement `.count`
    return pd.notnull(x).sum()


def mean_aggregate(s, n):
    try:
        return s / n
    except ZeroDivisionError:
        return np.float64(np.nan)


def var_aggregate(x2, x, n, ddof):
    try:
        result = (x2 / n) - (x / n)**2
        if ddof != 0:
            result = result * n / (n - ddof)
        return result
    except ZeroDivisionError:
        return np.float64(np.nan)


def describe_aggregate(values):
    assert len(values) == 6
    count, mean, std, min, q, max = values
    typ = pd.DataFrame if isinstance(count, pd.Series) else pd.Series
    part1 = typ([count, mean, std, min],
                index=['count', 'mean', 'std', 'min'])
    q.index = ['25%', '50%', '75%']
    part3 = typ([max], index=['max'])
    return pd.concat([part1, q, part3])


def cummin_aggregate(x, y):
    if isinstance(x, (pd.Series, pd.DataFrame)):
        return x.where((x < y) | x.isnull(), y, axis=x.ndim - 1)
    else:       # scalar
        return x if x < y else y


def cummax_aggregate(x, y):
    if isinstance(x, (pd.Series, pd.DataFrame)):
        return x.where((x > y) | x.isnull(), y, axis=x.ndim - 1)
    else:       # scalar
        return x if x > y else y


def assign(df, *pairs):
    kwargs = dict(partition(2, pairs))
    return df.assign(**kwargs)


def unique(x, series_name=None):
    # unique returns np.ndarray, it must be wrapped
    return pd.Series(x.unique(), name=series_name)


def value_counts_combine(x):
    return x.groupby(level=0).sum()


def value_counts_aggregate(x):
    return x.groupby(level=0).sum().sort_values(ascending=False)


def nbytes(x):
    return x.nbytes


def size(x):
    return x.size


def values(df):
    return df.values


def sample(df, state, frac, replace):
    rs = np.random.RandomState(state)
    return df.sample(random_state=rs, frac=frac, replace=replace) if len(df) > 0 else df


def drop_columns(df, columns, dtype):
    df = df.drop(columns, axis=1)
    df.columns = df.columns.astype(dtype)
    return df


def fillna_check(df, method, check=True):
    out = df.fillna(method=method)
    if check and out.isnull().values.all(axis=0).any():
        raise ValueError("All NaN partition encountered in `fillna`. Try "
                         "using ``df.repartition`` to increase the partition "
                         "size, or specify `limit` in `fillna`.")
    return out


# ---------------------------------
# reshape
# ---------------------------------


def pivot_agg(df):
    return df.groupby(level=0).sum()


def pivot_sum(df, index, columns, values):
    return pd.pivot_table(df, index=index, columns=columns,
                          values=values, aggfunc='sum')


def pivot_count(df, index, columns, values):
    # we cannot determine dtype until concatenationg all partitions.
    # make dtype deterministic, always coerce to np.float64
    return pd.pivot_table(df, index=index, columns=columns,
                          values=values, aggfunc='count').astype(np.float64)


# ---------------------------------
# concat
# ---------------------------------

if PANDAS_VERSION < '0.20.0':
    def _get_level_values(x, n):
        return x.get_level_values(n)
else:
    def _get_level_values(x, n):
        return x._get_level_values(n)


def concat(dfs, axis=0, join='outer', uniform=False):
    """Concatenate, handling some edge cases:

    - Unions categoricals between partitions
    - Ignores empty partitions

    Parameters
    ----------
    dfs : list of DataFrame, Series, or Index
    axis : int or str, optional
    join : str, optional
    uniform : bool, optional
        Whether to treat ``dfs[0]`` as representative of ``dfs[1:]``. Set to
        True if all arguments have the same columns and dtypes (but not
        necessarily categories). Default is False.
    """
    if axis == 1:
        return pd.concat(dfs, axis=axis, join=join)

    if len(dfs) == 1:
        return dfs[0]

    # Support concatenating indices along axis 0
    if isinstance(dfs[0], pd.Index):
        if isinstance(dfs[0], pd.CategoricalIndex):
            return pd.CategoricalIndex(union_categoricals(dfs),
                                       name=dfs[0].name)
        elif isinstance(dfs[0], pd.MultiIndex):
            first, rest = dfs[0], dfs[1:]
            if all((isinstance(o, pd.MultiIndex) and o.nlevels >= first.nlevels)
                    for o in rest):
                arrays = [concat([_get_level_values(i, n) for i in dfs])
                          for n in range(first.nlevels)]
                return pd.MultiIndex.from_arrays(arrays, names=first.names)

            to_concat = (first.values, ) + tuple(k._values for k in rest)
            new_tuples = np.concatenate(to_concat)
            try:
                return pd.MultiIndex.from_tuples(new_tuples, names=first.names)
            except Exception:
                return pd.Index(new_tuples)
        return dfs[0].append(dfs[1:])

    # Handle categorical index separately
    dfs0_index = dfs[0].index

    has_categoricalindex = (
        isinstance(dfs0_index, pd.CategoricalIndex) or
        (isinstance(dfs0_index, pd.MultiIndex) and
         any(isinstance(i, pd.CategoricalIndex) for i in dfs0_index.levels)))

    if has_categoricalindex:
        dfs2 = [df.reset_index(drop=True) for df in dfs]
        ind = concat([df.index for df in dfs])
    else:
        dfs2 = dfs
        ind = None

    # Concatenate the partitions together, handling categories as needed
    if (isinstance(dfs2[0], pd.DataFrame) if uniform else
            any(isinstance(df, pd.DataFrame) for df in dfs2)):
        if uniform:
            dfs3 = dfs2
            cat_mask = dfs2[0].dtypes == 'category'
        else:
            # When concatenating mixed dataframes and series on axis 1, Pandas
            # converts series to dataframes with a single column named 0, then
            # concatenates.
            dfs3 = [df if isinstance(df, pd.DataFrame) else
                    df.to_frame().rename(columns={df.name: 0}) for df in dfs2]
            # pandas may raise a RuntimeWarning for comparing ints and strs
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                cat_mask = pd.concat([(df.dtypes == 'category').to_frame().T
                                      for df in dfs3], join=join).any()

        if cat_mask.any():
            not_cat = cat_mask[~cat_mask].index
            out = pd.concat([df[df.columns.intersection(not_cat)]
                             for df in dfs3], join=join)
            temp_ind = out.index
            for col in cat_mask.index.difference(not_cat):
                # Find an example of categoricals in this column
                for df in dfs3:
                    sample = df.get(col)
                    if sample is not None:
                        break
                # Extract partitions, subbing in missing if needed
                parts = []
                for df in dfs3:
                    if col in df.columns:
                        parts.append(df[col])
                    else:
                        codes = np.full(len(df), -1, dtype='i8')
                        data = pd.Categorical.from_codes(codes,
                                                         sample.cat.categories,
                                                         sample.cat.ordered)
                        parts.append(data)
                out[col] = union_categoricals(parts)
                # Pandas resets index type on assignment if frame is empty
                # https://github.com/pandas-dev/pandas/issues/17101
                if not len(temp_ind):
                    out.index = temp_ind
            out = out.reindex(columns=cat_mask.index)
        else:
            # pandas may raise a RuntimeWarning for comparing ints and strs
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                out = pd.concat(dfs3, join=join)
    else:
        if is_categorical_dtype(dfs2[0].dtype):
            if ind is None:
                ind = concat([df.index for df in dfs2])
            return pd.Series(union_categoricals(dfs2), index=ind,
                             name=dfs2[0].name)
        out = pd.concat(dfs2, join=join)
    # Re-add the index if needed
    if ind is not None:
        out.index = ind
    return out


def merge(left, right, how, left_on, right_on,
          left_index, right_index, indicator, suffixes,
          default_left, default_right):

    if not len(left):
        left = default_left

    if not len(right):
        right = default_right

    return pd.merge(left, right, how=how, left_on=left_on, right_on=right_on,
                    left_index=left_index, right_index=right_index,
                    suffixes=suffixes, indicator=indicator)
