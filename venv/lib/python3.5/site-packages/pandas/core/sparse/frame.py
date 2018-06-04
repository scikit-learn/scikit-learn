"""
Data structures for sparse float data. Life is made simpler by dealing only
with float64 data
"""
from __future__ import division
# pylint: disable=E1101,E1103,W0231,E0202

import warnings
from pandas.compat import lmap
from pandas import compat
import numpy as np

from pandas.core.dtypes.missing import isna, notna
from pandas.core.dtypes.cast import maybe_upcast, find_common_type
from pandas.core.dtypes.common import _ensure_platform_int, is_scipy_sparse

from pandas.compat.numpy import function as nv
from pandas.core.index import Index, MultiIndex, _ensure_index
from pandas.core.series import Series
from pandas.core.frame import DataFrame, extract_index, _prep_ndarray
import pandas.core.algorithms as algos
from pandas.core.internals import (BlockManager,
                                   create_block_manager_from_arrays)
import pandas.core.generic as generic
from pandas.core.sparse.series import SparseSeries, SparseArray
from pandas._libs.sparse import BlockIndex, get_blocks
from pandas.util._decorators import Appender
import pandas.core.ops as ops
import pandas.core.common as com

_shared_doc_kwargs = dict(klass='SparseDataFrame')


class SparseDataFrame(DataFrame):
    """
    DataFrame containing sparse floating point data in the form of SparseSeries
    objects

    Parameters
    ----------
    data : same types as can be passed to DataFrame or scipy.sparse.spmatrix
        .. versionchanged :: 0.23.0
           If data is a dict, argument order is maintained for Python 3.6
           and later.

    index : array-like, optional
    column : array-like, optional
    default_kind : {'block', 'integer'}, default 'block'
        Default sparse kind for converting Series to SparseSeries. Will not
        override SparseSeries passed into constructor
    default_fill_value : float
        Default fill_value for converting Series to SparseSeries
        (default: nan). Will not override SparseSeries passed in.
    """
    _subtyp = 'sparse_frame'

    def __init__(self, data=None, index=None, columns=None, default_kind=None,
                 default_fill_value=None, dtype=None, copy=False):

        # pick up the defaults from the Sparse structures
        if isinstance(data, SparseDataFrame):
            if index is None:
                index = data.index
            if columns is None:
                columns = data.columns
            if default_fill_value is None:
                default_fill_value = data.default_fill_value
            if default_kind is None:
                default_kind = data.default_kind
        elif isinstance(data, (SparseSeries, SparseArray)):
            if index is None:
                index = data.index
            if default_fill_value is None:
                default_fill_value = data.fill_value
            if columns is None and hasattr(data, 'name'):
                columns = [data.name]
            if columns is None:
                raise Exception("cannot pass a series w/o a name or columns")
            data = {columns[0]: data}

        if default_fill_value is None:
            default_fill_value = np.nan
        if default_kind is None:
            default_kind = 'block'

        self._default_kind = default_kind
        self._default_fill_value = default_fill_value

        if is_scipy_sparse(data):
            mgr = self._init_spmatrix(data, index, columns, dtype=dtype,
                                      fill_value=default_fill_value)
        elif isinstance(data, dict):
            mgr = self._init_dict(data, index, columns, dtype=dtype)
        elif isinstance(data, (np.ndarray, list)):
            mgr = self._init_matrix(data, index, columns, dtype=dtype)
        elif isinstance(data, SparseDataFrame):
            mgr = self._init_mgr(data._data,
                                 dict(index=index, columns=columns),
                                 dtype=dtype, copy=copy)
        elif isinstance(data, DataFrame):
            mgr = self._init_dict(data, data.index, data.columns, dtype=dtype)
        elif isinstance(data, Series):
            mgr = self._init_dict(data.to_frame(), data.index,
                                  columns=None, dtype=dtype)
        elif isinstance(data, BlockManager):
            mgr = self._init_mgr(data, axes=dict(index=index, columns=columns),
                                 dtype=dtype, copy=copy)
        elif data is None:
            data = DataFrame()

            if index is None:
                index = Index([])
            else:
                index = _ensure_index(index)

            if columns is None:
                columns = Index([])
            else:
                for c in columns:
                    data[c] = SparseArray(np.nan, index=index,
                                          kind=self._default_kind,
                                          fill_value=self._default_fill_value)
            mgr = to_manager(data, columns, index)
            if dtype is not None:
                mgr = mgr.astype(dtype)
        else:
            msg = ('SparseDataFrame called with unknown type "{data_type}" '
                   'for data argument')
            raise TypeError(msg.format(data_type=type(data).__name__))

        generic.NDFrame.__init__(self, mgr)

    @property
    def _constructor(self):
        return SparseDataFrame

    _constructor_sliced = SparseSeries

    def _init_dict(self, data, index, columns, dtype=None):
        # pre-filter out columns if we passed it
        if columns is not None:
            columns = _ensure_index(columns)
            data = {k: v for k, v in compat.iteritems(data) if k in columns}
        else:
            keys = com._dict_keys_to_ordered_list(data)
            columns = Index(keys)

        if index is None:
            index = extract_index(list(data.values()))

        sp_maker = lambda x: SparseArray(x, kind=self._default_kind,
                                         fill_value=self._default_fill_value,
                                         copy=True, dtype=dtype)
        sdict = {}
        for k, v in compat.iteritems(data):
            if isinstance(v, Series):
                # Force alignment, no copy necessary
                if not v.index.equals(index):
                    v = v.reindex(index)

                if not isinstance(v, SparseSeries):
                    v = sp_maker(v.values)
            elif isinstance(v, SparseArray):
                v = v.copy()
            else:
                if isinstance(v, dict):
                    v = [v.get(i, np.nan) for i in index]

                v = sp_maker(v)
            sdict[k] = v

        # TODO: figure out how to handle this case, all nan's?
        # add in any other columns we want to have (completeness)
        nan_arr = np.empty(len(index), dtype='float64')
        nan_arr.fill(np.nan)
        nan_arr = sp_maker(nan_arr)
        sdict.update((c, nan_arr) for c in columns if c not in sdict)

        return to_manager(sdict, columns, index)

    def _init_matrix(self, data, index, columns, dtype=None):
        """ Init self from ndarray or list of lists """
        data = _prep_ndarray(data, copy=False)
        index, columns = self._prep_index(data, index, columns)
        data = {idx: data[:, i] for i, idx in enumerate(columns)}
        return self._init_dict(data, index, columns, dtype)

    def _init_spmatrix(self, data, index, columns, dtype=None,
                       fill_value=None):
        """ Init self from scipy.sparse matrix """
        index, columns = self._prep_index(data, index, columns)
        data = data.tocoo()
        N = len(index)

        # Construct a dict of SparseSeries
        sdict = {}
        values = Series(data.data, index=data.row, copy=False)
        for col, rowvals in values.groupby(data.col):
            # get_blocks expects int32 row indices in sorted order
            rowvals = rowvals.sort_index()
            rows = rowvals.index.values.astype(np.int32)
            blocs, blens = get_blocks(rows)

            sdict[columns[col]] = SparseSeries(
                rowvals.values, index=index,
                fill_value=fill_value,
                sparse_index=BlockIndex(N, blocs, blens))

        # Add any columns that were empty and thus not grouped on above
        sdict.update({column: SparseSeries(index=index,
                                           fill_value=fill_value,
                                           sparse_index=BlockIndex(N, [], []))
                      for column in columns
                      if column not in sdict})

        return self._init_dict(sdict, index, columns, dtype)

    def _prep_index(self, data, index, columns):
        N, K = data.shape
        if index is None:
            index = com._default_index(N)
        if columns is None:
            columns = com._default_index(K)

        if len(columns) != K:
            raise ValueError('Column length mismatch: {columns} vs. {K}'
                             .format(columns=len(columns), K=K))
        if len(index) != N:
            raise ValueError('Index length mismatch: {index} vs. {N}'
                             .format(index=len(index), N=N))
        return index, columns

    def to_coo(self):
        """
        Return the contents of the frame as a sparse SciPy COO matrix.

        .. versionadded:: 0.20.0

        Returns
        -------
        coo_matrix : scipy.sparse.spmatrix
            If the caller is heterogeneous and contains booleans or objects,
            the result will be of dtype=object. See Notes.

        Notes
        -----
        The dtype will be the lowest-common-denominator type (implicit
        upcasting); that is to say if the dtypes (even of numeric types)
        are mixed, the one that accommodates all will be chosen.

        e.g. If the dtypes are float16 and float32, dtype will be upcast to
        float32. By numpy.find_common_type convention, mixing int64 and
        and uint64 will result in a float64 dtype.
        """
        try:
            from scipy.sparse import coo_matrix
        except ImportError:
            raise ImportError('Scipy is not installed')

        dtype = find_common_type(self.dtypes)
        cols, rows, datas = [], [], []
        for col, name in enumerate(self):
            s = self[name]
            row = s.sp_index.to_int_index().indices
            cols.append(np.repeat(col, len(row)))
            rows.append(row)
            datas.append(s.sp_values.astype(dtype, copy=False))

        cols = np.concatenate(cols)
        rows = np.concatenate(rows)
        datas = np.concatenate(datas)
        return coo_matrix((datas, (rows, cols)), shape=self.shape)

    def __array_wrap__(self, result):
        return self._constructor(
            result, index=self.index, columns=self.columns,
            default_kind=self._default_kind,
            default_fill_value=self._default_fill_value).__finalize__(self)

    def __getstate__(self):
        # pickling
        return dict(_typ=self._typ, _subtyp=self._subtyp, _data=self._data,
                    _default_fill_value=self._default_fill_value,
                    _default_kind=self._default_kind)

    def _unpickle_sparse_frame_compat(self, state):
        """ original pickle format """
        series, cols, idx, fv, kind = state

        if not isinstance(cols, Index):  # pragma: no cover
            from pandas.io.pickle import _unpickle_array
            columns = _unpickle_array(cols)
        else:
            columns = cols

        if not isinstance(idx, Index):  # pragma: no cover
            from pandas.io.pickle import _unpickle_array
            index = _unpickle_array(idx)
        else:
            index = idx

        series_dict = DataFrame()
        for col, (sp_index, sp_values) in compat.iteritems(series):
            series_dict[col] = SparseSeries(sp_values, sparse_index=sp_index,
                                            fill_value=fv)

        self._data = to_manager(series_dict, columns, index)
        self._default_fill_value = fv
        self._default_kind = kind

    def to_dense(self):
        """
        Convert to dense DataFrame

        Returns
        -------
        df : DataFrame
        """
        data = {k: v.to_dense() for k, v in compat.iteritems(self)}
        return DataFrame(data, index=self.index, columns=self.columns)

    def _apply_columns(self, func):
        """ get new SparseDataFrame applying func to each columns """

        new_data = {}
        for col, series in compat.iteritems(self):
            new_data[col] = func(series)

        return self._constructor(
            data=new_data, index=self.index, columns=self.columns,
            default_fill_value=self.default_fill_value).__finalize__(self)

    def astype(self, dtype):
        return self._apply_columns(lambda x: x.astype(dtype))

    def copy(self, deep=True):
        """
        Make a copy of this SparseDataFrame
        """
        result = super(SparseDataFrame, self).copy(deep=deep)
        result._default_fill_value = self._default_fill_value
        result._default_kind = self._default_kind
        return result

    @property
    def default_fill_value(self):
        return self._default_fill_value

    @property
    def default_kind(self):
        return self._default_kind

    @property
    def density(self):
        """
        Ratio of non-sparse points to total (dense) data points
        represented in the frame
        """
        tot_nonsparse = sum(ser.sp_index.npoints
                            for _, ser in compat.iteritems(self))
        tot = len(self.index) * len(self.columns)
        return tot_nonsparse / float(tot)

    def fillna(self, value=None, method=None, axis=0, inplace=False,
               limit=None, downcast=None):
        new_self = super(SparseDataFrame,
                         self).fillna(value=value, method=method, axis=axis,
                                      inplace=inplace, limit=limit,
                                      downcast=downcast)
        if not inplace:
            self = new_self

        # set the fill value if we are filling as a scalar with nothing special
        # going on
        if (value is not None and value == value and method is None and
                limit is None):
            self._default_fill_value = value

        if not inplace:
            return self

    # ----------------------------------------------------------------------
    # Support different internal representation of SparseDataFrame

    def _sanitize_column(self, key, value, **kwargs):
        """
        Creates a new SparseArray from the input value.

        Parameters
        ----------
        key : object
        value : scalar, Series, or array-like
        kwargs : dict

        Returns
        -------
        sanitized_column : SparseArray

        """
        sp_maker = lambda x, index=None: SparseArray(
            x, index=index, fill_value=self._default_fill_value,
            kind=self._default_kind)
        if isinstance(value, SparseSeries):
            clean = value.reindex(self.index).as_sparse_array(
                fill_value=self._default_fill_value, kind=self._default_kind)

        elif isinstance(value, SparseArray):
            if len(value) != len(self.index):
                raise AssertionError('Length of values does not match '
                                     'length of index')
            clean = value

        elif hasattr(value, '__iter__'):
            if isinstance(value, Series):
                clean = value.reindex(self.index)
                if not isinstance(value, SparseSeries):
                    clean = sp_maker(clean)
            else:
                if len(value) != len(self.index):
                    raise AssertionError('Length of values does not match '
                                         'length of index')
                clean = sp_maker(value)

        # Scalar
        else:
            clean = sp_maker(value, self.index)

        # always return a SparseArray!
        return clean

    def __getitem__(self, key):
        """
        Retrieve column or slice from DataFrame
        """
        if isinstance(key, slice):
            date_rng = self.index[key]
            return self.reindex(date_rng)
        elif isinstance(key, (np.ndarray, list, Series)):
            return self._getitem_array(key)
        else:
            return self._get_item_cache(key)

    def get_value(self, index, col, takeable=False):
        """
        Quickly retrieve single value at passed column and index

        .. deprecated:: 0.21.0

        Please use .at[] or .iat[] accessors.

        Parameters
        ----------
        index : row label
        col : column label
        takeable : interpret the index/col as indexers, default False

        Returns
        -------
        value : scalar value
        """
        warnings.warn("get_value is deprecated and will be removed "
                      "in a future release. Please use "
                      ".at[] or .iat[] accessors instead", FutureWarning,
                      stacklevel=2)
        return self._get_value(index, col, takeable=takeable)

    def _get_value(self, index, col, takeable=False):
        if takeable is True:
            series = self._iget_item_cache(col)
        else:
            series = self._get_item_cache(col)

        return series._get_value(index, takeable=takeable)
    _get_value.__doc__ = get_value.__doc__

    def set_value(self, index, col, value, takeable=False):
        """
        Put single value at passed column and index

        .. deprecated:: 0.21.0

        Please use .at[] or .iat[] accessors.

        Parameters
        ----------
        index : row label
        col : column label
        value : scalar value
        takeable : interpret the index/col as indexers, default False

        Notes
        -----
        This method *always* returns a new object. It is currently not
        particularly efficient (and potentially very expensive) but is provided
        for API compatibility with DataFrame

        Returns
        -------
        frame : DataFrame
        """
        warnings.warn("set_value is deprecated and will be removed "
                      "in a future release. Please use "
                      ".at[] or .iat[] accessors instead", FutureWarning,
                      stacklevel=2)
        return self._set_value(index, col, value, takeable=takeable)

    def _set_value(self, index, col, value, takeable=False):
        dense = self.to_dense()._set_value(
            index, col, value, takeable=takeable)
        return dense.to_sparse(kind=self._default_kind,
                               fill_value=self._default_fill_value)
    _set_value.__doc__ = set_value.__doc__

    def _slice(self, slobj, axis=0, kind=None):
        if axis == 0:
            new_index = self.index[slobj]
            new_columns = self.columns
        else:
            new_index = self.index
            new_columns = self.columns[slobj]

        return self.reindex(index=new_index, columns=new_columns)

    def xs(self, key, axis=0, copy=False):
        """
        Returns a row (cross-section) from the SparseDataFrame as a Series
        object.

        Parameters
        ----------
        key : some index contained in the index

        Returns
        -------
        xs : Series
        """
        if axis == 1:
            data = self[key]
            return data

        i = self.index.get_loc(key)
        data = self.take([i]).get_values()[0]
        return Series(data, index=self.columns)

    # ----------------------------------------------------------------------
    # Arithmetic-related methods

    def _combine_frame(self, other, func, fill_value=None, level=None):
        this, other = self.align(other, join='outer', level=level, copy=False)
        new_index, new_columns = this.index, this.columns

        if level is not None:
            raise NotImplementedError("'level' argument is not supported")

        if self.empty and other.empty:
            return self._constructor(index=new_index).__finalize__(self)

        new_data = {}
        if fill_value is not None:
            # TODO: be a bit more intelligent here
            for col in new_columns:
                if col in this and col in other:
                    dleft = this[col].to_dense()
                    dright = other[col].to_dense()
                    result = dleft._binop(dright, func, fill_value=fill_value)
                    result = result.to_sparse(fill_value=this[col].fill_value)
                    new_data[col] = result
        else:

            for col in new_columns:
                if col in this and col in other:
                    new_data[col] = func(this[col], other[col])

        # if the fill values are the same use them? or use a valid one
        new_fill_value = None
        other_fill_value = getattr(other, 'default_fill_value', np.nan)
        if self.default_fill_value == other_fill_value:
            new_fill_value = self.default_fill_value
        elif np.isnan(self.default_fill_value) and not np.isnan(
                other_fill_value):
            new_fill_value = other_fill_value
        elif not np.isnan(self.default_fill_value) and np.isnan(
                other_fill_value):
            new_fill_value = self.default_fill_value

        return self._constructor(data=new_data, index=new_index,
                                 columns=new_columns,
                                 default_fill_value=new_fill_value
                                 ).__finalize__(self)

    def _combine_match_index(self, other, func, level=None):
        new_data = {}

        if level is not None:
            raise NotImplementedError("'level' argument is not supported")

        new_index = self.index.union(other.index)
        this = self
        if self.index is not new_index:
            this = self.reindex(new_index)

        if other.index is not new_index:
            other = other.reindex(new_index)

        for col, series in compat.iteritems(this):
            new_data[col] = func(series.values, other.values)

        # fill_value is a function of our operator
        fill_value = None
        if isna(other.fill_value) or isna(self.default_fill_value):
            fill_value = np.nan
        else:
            fill_value = func(np.float64(self.default_fill_value),
                              np.float64(other.fill_value))

        return self._constructor(
            new_data, index=new_index, columns=self.columns,
            default_fill_value=fill_value).__finalize__(self)

    def _combine_match_columns(self, other, func, level=None, try_cast=True):
        # patched version of DataFrame._combine_match_columns to account for
        # NumPy circumventing __rsub__ with float64 types, e.g.: 3.0 - series,
        # where 3.0 is numpy.float64 and series is a SparseSeries. Still
        # possible for this to happen, which is bothersome

        if level is not None:
            raise NotImplementedError("'level' argument is not supported")

        new_data = {}

        union = intersection = self.columns

        if not union.equals(other.index):
            union = other.index.union(self.columns)
            intersection = other.index.intersection(self.columns)

        for col in intersection:
            new_data[col] = func(self[col], float(other[col]))

        return self._constructor(
            new_data, index=self.index, columns=union,
            default_fill_value=self.default_fill_value).__finalize__(self)

    def _combine_const(self, other, func, errors='raise', try_cast=True):
        return self._apply_columns(lambda x: func(x, other))

    def _reindex_index(self, index, method, copy, level, fill_value=np.nan,
                       limit=None, takeable=False):
        if level is not None:
            raise TypeError('Reindex by level not supported for sparse')

        if self.index.equals(index):
            if copy:
                return self.copy()
            else:
                return self

        if len(self.index) == 0:
            return self._constructor(
                index=index, columns=self.columns).__finalize__(self)

        indexer = self.index.get_indexer(index, method, limit=limit)
        indexer = _ensure_platform_int(indexer)
        mask = indexer == -1
        need_mask = mask.any()

        new_series = {}
        for col, series in self.iteritems():
            if mask.all():
                continue

            values = series.values
            # .take returns SparseArray
            new = values.take(indexer)
            if need_mask:
                new = new.values
                # convert integer to float if necessary. need to do a lot
                # more than that, handle boolean etc also
                new, fill_value = maybe_upcast(new, fill_value=fill_value)
                np.putmask(new, mask, fill_value)

            new_series[col] = new

        return self._constructor(
            new_series, index=index, columns=self.columns,
            default_fill_value=self._default_fill_value).__finalize__(self)

    def _reindex_columns(self, columns, method, copy, level, fill_value=None,
                         limit=None, takeable=False):
        if level is not None:
            raise TypeError('Reindex by level not supported for sparse')

        if notna(fill_value):
            raise NotImplementedError("'fill_value' argument is not supported")

        if limit:
            raise NotImplementedError("'limit' argument is not supported")

        if method is not None:
            raise NotImplementedError("'method' argument is not supported")

        # TODO: fill value handling
        sdict = {k: v for k, v in compat.iteritems(self) if k in columns}
        return self._constructor(
            sdict, index=self.index, columns=columns,
            default_fill_value=self._default_fill_value).__finalize__(self)

    def _reindex_with_indexers(self, reindexers, method=None, fill_value=None,
                               limit=None, copy=False, allow_dups=False):

        if method is not None or limit is not None:
            raise NotImplementedError("cannot reindex with a method or limit "
                                      "with sparse")

        if fill_value is None:
            fill_value = np.nan

        reindexers = {self._get_axis_number(a): val
                      for (a, val) in compat.iteritems(reindexers)}

        index, row_indexer = reindexers.get(0, (None, None))
        columns, col_indexer = reindexers.get(1, (None, None))

        if columns is None:
            columns = self.columns

        new_arrays = {}
        for col in columns:
            if col not in self:
                continue
            if row_indexer is not None:
                new_arrays[col] = algos.take_1d(self[col].get_values(),
                                                row_indexer,
                                                fill_value=fill_value)
            else:
                new_arrays[col] = self[col]

        return self._constructor(new_arrays, index=index,
                                 columns=columns).__finalize__(self)

    def _join_compat(self, other, on=None, how='left', lsuffix='', rsuffix='',
                     sort=False):
        if on is not None:
            raise NotImplementedError("'on' keyword parameter is not yet "
                                      "implemented")
        return self._join_index(other, how, lsuffix, rsuffix)

    def _join_index(self, other, how, lsuffix, rsuffix):
        if isinstance(other, Series):
            if other.name is None:
                raise ValueError('Other Series must have a name')

            other = SparseDataFrame(
                {other.name: other},
                default_fill_value=self._default_fill_value)

        join_index = self.index.join(other.index, how=how)

        this = self.reindex(join_index)
        other = other.reindex(join_index)

        this, other = this._maybe_rename_join(other, lsuffix, rsuffix)

        from pandas import concat
        return concat([this, other], axis=1, verify_integrity=True)

    def _maybe_rename_join(self, other, lsuffix, rsuffix):
        to_rename = self.columns.intersection(other.columns)
        if len(to_rename) > 0:
            if not lsuffix and not rsuffix:
                raise ValueError('columns overlap but no suffix specified: '
                                 '{to_rename}'.format(to_rename=to_rename))

            def lrenamer(x):
                if x in to_rename:
                    return '{x}{lsuffix}'.format(x=x, lsuffix=lsuffix)
                return x

            def rrenamer(x):
                if x in to_rename:
                    return '{x}{rsuffix}'.format(x=x, rsuffix=rsuffix)
                return x

            this = self.rename(columns=lrenamer)
            other = other.rename(columns=rrenamer)
        else:
            this = self

        return this, other

    def transpose(self, *args, **kwargs):
        """
        Returns a DataFrame with the rows/columns switched.
        """
        nv.validate_transpose(args, kwargs)
        return self._constructor(
            self.values.T, index=self.columns, columns=self.index,
            default_fill_value=self._default_fill_value,
            default_kind=self._default_kind).__finalize__(self)

    T = property(transpose)

    @Appender(DataFrame.count.__doc__)
    def count(self, axis=0, **kwds):
        if axis is None:
            axis = self._stat_axis_number

        return self.apply(lambda x: x.count(), axis=axis)

    def cumsum(self, axis=0, *args, **kwargs):
        """
        Return SparseDataFrame of cumulative sums over requested axis.

        Parameters
        ----------
        axis : {0, 1}
            0 for row-wise, 1 for column-wise

        Returns
        -------
        y : SparseDataFrame
        """
        nv.validate_cumsum(args, kwargs)

        if axis is None:
            axis = self._stat_axis_number

        return self.apply(lambda x: x.cumsum(), axis=axis)

    @Appender(generic._shared_docs['isna'] % _shared_doc_kwargs)
    def isna(self):
        return self._apply_columns(lambda x: x.isna())
    isnull = isna

    @Appender(generic._shared_docs['notna'] % _shared_doc_kwargs)
    def notna(self):
        return self._apply_columns(lambda x: x.notna())
    notnull = notna

    def apply(self, func, axis=0, broadcast=None, reduce=None,
              result_type=None):
        """
        Analogous to DataFrame.apply, for SparseDataFrame

        Parameters
        ----------
        func : function
            Function to apply to each column
        axis : {0, 1, 'index', 'columns'}
        broadcast : bool, default False
            For aggregation functions, return object of same size with values
            propagated

            .. deprecated:: 0.23.0
               This argument will be removed in a future version, replaced
               by result_type='broadcast'.

        reduce : boolean or None, default None
            Try to apply reduction procedures. If the DataFrame is empty,
            apply will use reduce to determine whether the result should be a
            Series or a DataFrame. If reduce is None (the default), apply's
            return value will be guessed by calling func an empty Series (note:
            while guessing, exceptions raised by func will be ignored). If
            reduce is True a Series will always be returned, and if False a
            DataFrame will always be returned.

            .. deprecated:: 0.23.0
               This argument will be removed in a future version, replaced
               by result_type='reduce'.

        result_type : {'expand', 'reduce', 'broadcast, None}
            These only act when axis=1 {columns}:

            * 'expand' : list-like results will be turned into columns.
            * 'reduce' : return a Series if possible rather than expanding
              list-like results. This is the opposite to 'expand'.
            * 'broadcast' : results will be broadcast to the original shape
              of the frame, the original index & columns will be retained.

            The default behaviour (None) depends on the return value of the
            applied function: list-like results will be returned as a Series
            of those. However if the apply function returns a Series these
            are expanded to columns.

            .. versionadded:: 0.23.0

        Returns
        -------
        applied : Series or SparseDataFrame
        """
        if not len(self.columns):
            return self
        axis = self._get_axis_number(axis)

        if isinstance(func, np.ufunc):
            new_series = {}
            for k, v in compat.iteritems(self):
                applied = func(v)
                applied.fill_value = func(v.fill_value)
                new_series[k] = applied
            return self._constructor(
                new_series, index=self.index, columns=self.columns,
                default_fill_value=self._default_fill_value,
                default_kind=self._default_kind).__finalize__(self)

        from pandas.core.apply import frame_apply
        op = frame_apply(self,
                         func=func,
                         axis=axis,
                         reduce=reduce,
                         broadcast=broadcast,
                         result_type=result_type)
        return op.get_result()

    def applymap(self, func):
        """
        Apply a function to a DataFrame that is intended to operate
        elementwise, i.e. like doing map(func, series) for each series in the
        DataFrame

        Parameters
        ----------
        func : function
            Python function, returns a single value from a single value

        Returns
        -------
        applied : DataFrame
        """
        return self.apply(lambda x: lmap(func, x))


def to_manager(sdf, columns, index):
    """ create and return the block manager from a dataframe of series,
    columns, index
    """

    # from BlockManager perspective
    axes = [_ensure_index(columns), _ensure_index(index)]

    return create_block_manager_from_arrays(
        [sdf[c] for c in columns], columns, axes)


def stack_sparse_frame(frame):
    """
    Only makes sense when fill_value is NaN
    """
    lengths = [s.sp_index.npoints for _, s in compat.iteritems(frame)]
    nobs = sum(lengths)

    # this is pretty fast
    minor_labels = np.repeat(np.arange(len(frame.columns)), lengths)

    inds_to_concat = []
    vals_to_concat = []
    # TODO: Figure out whether this can be reached.
    # I think this currently can't be reached because you can't build a
    # SparseDataFrame with a non-np.NaN fill value (fails earlier).
    for _, series in compat.iteritems(frame):
        if not np.isnan(series.fill_value):
            raise TypeError('This routine assumes NaN fill value')

        int_index = series.sp_index.to_int_index()
        inds_to_concat.append(int_index.indices)
        vals_to_concat.append(series.sp_values)

    major_labels = np.concatenate(inds_to_concat)
    stacked_values = np.concatenate(vals_to_concat)
    index = MultiIndex(levels=[frame.index, frame.columns],
                       labels=[major_labels, minor_labels],
                       verify_integrity=False)

    lp = DataFrame(stacked_values.reshape((nobs, 1)), index=index,
                   columns=['foo'])
    return lp.sort_index(level=0)


def homogenize(series_dict):
    """
    Conform a set of SparseSeries (with NaN fill_value) to a common SparseIndex
    corresponding to the locations where they all have data

    Parameters
    ----------
    series_dict : dict or DataFrame

    Notes
    -----
    Using the dumbest algorithm I could think of. Should put some more thought
    into this

    Returns
    -------
    homogenized : dict of SparseSeries
    """
    index = None

    need_reindex = False

    for _, series in compat.iteritems(series_dict):
        if not np.isnan(series.fill_value):
            raise TypeError('this method is only valid with NaN fill values')

        if index is None:
            index = series.sp_index
        elif not series.sp_index.equals(index):
            need_reindex = True
            index = index.intersect(series.sp_index)

    if need_reindex:
        output = {}
        for name, series in compat.iteritems(series_dict):
            if not series.sp_index.equals(index):
                series = series.sparse_reindex(index)

            output[name] = series
    else:
        output = series_dict

    return output


# use unaccelerated ops for sparse objects
ops.add_flex_arithmetic_methods(SparseDataFrame)
ops.add_special_arithmetic_methods(SparseDataFrame)
