from __future__ import absolute_import, division, print_function

import numpy as np
import pandas as pd

from .core import Series, DataFrame, map_partitions, apply_concat_apply
from . import methods
from .utils import is_categorical_dtype, is_scalar, has_known_categories


###############################################################
# Dummies
###############################################################


def get_dummies(data, prefix=None, prefix_sep='_', dummy_na=False,
                columns=None, sparse=False, drop_first=False):
    """
    Convert categorical variable into dummy/indicator variables. Data must
    have category dtype to infer result's ``columns``

    Parameters
    ----------
    data : Series or DataFrame with category dtype
    prefix : string, list of strings, or dict of strings, default None
        String to append DataFrame column names
        Pass a list with length equal to the number of columns
        when calling get_dummies on a DataFrame. Alternativly, `prefix`
        can be a dictionary mapping column names to prefixes.
    prefix_sep : string, default '_'
        If appending prefix, separator/delimiter to use. Or pass a
        list or dictionary as with `prefix.`
    dummy_na : bool, default False
        Add a column to indicate NaNs, if False NaNs are ignored.
    columns : list-like, default None
        Column names in the DataFrame to be encoded.
        If `columns` is None then all the columns with
        `category` dtype will be converted.
    drop_first : bool, default False
        Whether to get k-1 dummies out of k categorical levels by removing the
        first level.
    Returns
    -------
    dummies : DataFrame
    """

    if isinstance(data, (pd.Series, pd.DataFrame)):
        return pd.get_dummies(data, prefix=prefix,
                              prefix_sep=prefix_sep, dummy_na=dummy_na,
                              columns=columns, sparse=sparse,
                              drop_first=drop_first)

    not_cat_msg = ("`get_dummies` with non-categorical dtypes is not "
                   "supported. Please use `df.categorize()` beforehand to "
                   "convert to categorical dtype.")

    unknown_cat_msg = ("`get_dummies` with unknown categories is not "
                       "supported. Please use `column.cat.as_known()` or "
                       "`df.categorize()` beforehand to ensure known "
                       "categories")

    if isinstance(data, Series):
        if not is_categorical_dtype(data):
            raise NotImplementedError(not_cat_msg)
        if not has_known_categories(data):
            raise NotImplementedError(unknown_cat_msg)
    elif isinstance(data, DataFrame):
        if columns is None:
            if (data.dtypes == 'object').any():
                raise NotImplementedError(not_cat_msg)
            columns = data._meta.select_dtypes(include=['category']).columns
        else:
            if not all(is_categorical_dtype(data[c]) for c in columns):
                raise NotImplementedError(not_cat_msg)

        if not all(has_known_categories(data[c]) for c in columns):
            raise NotImplementedError(unknown_cat_msg)

    if sparse:
        raise NotImplementedError('sparse=True is not supported')

    return map_partitions(pd.get_dummies, data, prefix=prefix,
                          prefix_sep=prefix_sep, dummy_na=dummy_na,
                          columns=columns, sparse=sparse,
                          drop_first=drop_first)


###############################################################
# Pivot table
###############################################################


def pivot_table(df, index=None, columns=None,
                values=None, aggfunc='mean'):
    """
    Create a spreadsheet-style pivot table as a DataFrame. Target ``columns``
    must have category dtype to infer result's ``columns``.
    ``index``, ``columns``, ``values`` and ``aggfunc`` must be all scalar.

    Parameters
    ----------
    data : DataFrame
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

    if not is_scalar(index) or index is None:
        raise ValueError("'index' must be the name of an existing column")
    if not is_scalar(columns) or columns is None:
        raise ValueError("'columns' must be the name of an existing column")
    if not is_categorical_dtype(df[columns]):
        raise ValueError("'columns' must be category dtype")
    if not has_known_categories(df[columns]):
        raise ValueError("'columns' must have known categories. Please use "
                         "`df[columns].cat.as_known()` beforehand to ensure "
                         "known categories")
    if not is_scalar(values) or values is None:
        raise ValueError("'values' must be the name of an existing column")
    if not is_scalar(aggfunc) or aggfunc not in ('mean', 'sum', 'count'):
        raise ValueError("aggfunc must be either 'mean', 'sum' or 'count'")

    # _emulate can't work for empty data
    # the result must have CategoricalIndex columns
    new_columns = pd.CategoricalIndex(df[columns].cat.categories, name=columns)
    meta = pd.DataFrame(columns=new_columns, dtype=np.float64)
    meta.index.name = index

    kwargs = {'index': index, 'columns': columns, 'values': values}

    pv_sum = apply_concat_apply([df],
                                chunk=methods.pivot_sum,
                                aggregate=methods.pivot_agg,
                                meta=meta,
                                token='pivot_table_sum',
                                chunk_kwargs=kwargs)
    pv_count = apply_concat_apply([df],
                                  chunk=methods.pivot_count,
                                  aggregate=methods.pivot_agg,
                                  meta=meta,
                                  token='pivot_table_count',
                                  chunk_kwargs=kwargs)

    if aggfunc == 'sum':
        return pv_sum
    elif aggfunc == 'count':
        return pv_count
    elif aggfunc == 'mean':
        return pv_sum / pv_count
    else:
        raise ValueError


###############################################################
# Melt
###############################################################


def melt(frame, id_vars=None, value_vars=None, var_name=None,
         value_name='value', col_level=None):

    from dask.dataframe.core import no_default

    return frame.map_partitions(pd.melt, meta=no_default, id_vars=id_vars,
                                value_vars=value_vars,
                                var_name=var_name, value_name=value_name,
                                col_level=col_level, token='melt')
