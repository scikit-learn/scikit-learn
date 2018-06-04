"""
Data structures for sparse float data. Life is made simpler by dealing only
with float64 data
"""

# pylint: disable=E1101,E1103,W0231

import numpy as np
import warnings

from pandas.core.dtypes.missing import isna, notna

from pandas.compat.numpy import function as nv
from pandas.core.index import Index, _ensure_index, InvalidIndexError
from pandas.core.series import Series
from pandas.core.internals import SingleBlockManager
from pandas.core import generic
import pandas.core.common as com
import pandas.core.ops as ops
import pandas._libs.index as libindex
from pandas.util._decorators import Appender

from pandas.core.sparse.array import (
    make_sparse, SparseArray,
    _make_index)
from pandas._libs.sparse import BlockIndex, IntIndex
import pandas._libs.sparse as splib

from pandas.core.sparse.scipy_sparse import (
    _sparse_series_to_coo,
    _coo_to_sparse_series)


_shared_doc_kwargs = dict(axes='index', klass='SparseSeries',
                          axes_single_arg="{0, 'index'}",
                          optional_labels='', optional_axis='')


class SparseSeries(Series):
    """Data structure for labeled, sparse floating point data

    Parameters
    ----------
    data : {array-like, Series, SparseSeries, dict}
        .. versionchanged :: 0.23.0
           If data is a dict, argument order is maintained for Python 3.6
           and later.

    kind : {'block', 'integer'}
    fill_value : float
        Code for missing value. Defaults depends on dtype.
        0 for int dtype, False for bool dtype, and NaN for other dtypes
    sparse_index : {BlockIndex, IntIndex}, optional
        Only if you have one. Mainly used internally

    Notes
    -----
    SparseSeries objects are immutable via the typical Python means. If you
    must change values, convert to dense, make your changes, then convert back
    to sparse
    """
    _subtyp = 'sparse_series'

    def __init__(self, data=None, index=None, sparse_index=None, kind='block',
                 fill_value=None, name=None, dtype=None, copy=False,
                 fastpath=False):

        # we are called internally, so short-circuit
        if fastpath:

            # data is an ndarray, index is defined

            if not isinstance(data, SingleBlockManager):
                data = SingleBlockManager(data, index, fastpath=True)
            if copy:
                data = data.copy()

        else:

            if data is None:
                data = []

            if isinstance(data, Series) and name is None:
                name = data.name

            if isinstance(data, SparseArray):
                if index is not None:
                    assert (len(index) == len(data))
                sparse_index = data.sp_index
                if fill_value is None:
                    fill_value = data.fill_value

                data = np.asarray(data)

            elif isinstance(data, SparseSeries):
                if index is None:
                    index = data.index.view()
                if fill_value is None:
                    fill_value = data.fill_value
                # extract the SingleBlockManager
                data = data._data

            elif isinstance(data, (Series, dict)):
                data = Series(data, index=index)
                index = data.index.view()

                res = make_sparse(data, kind=kind, fill_value=fill_value)
                data, sparse_index, fill_value = res

            elif isinstance(data, (tuple, list, np.ndarray)):
                # array-like
                if sparse_index is None:
                    res = make_sparse(data, kind=kind, fill_value=fill_value)
                    data, sparse_index, fill_value = res
                else:
                    assert (len(data) == sparse_index.npoints)

            elif isinstance(data, SingleBlockManager):
                if dtype is not None:
                    data = data.astype(dtype)
                if index is None:
                    index = data.index.view()
                elif not data.index.equals(index) or copy:  # pragma: no cover
                    # GH#19275 SingleBlockManager input should only be called
                    # internally
                    raise AssertionError('Cannot pass both SingleBlockManager '
                                         '`data` argument and a different '
                                         '`index` argument.  `copy` must '
                                         'be False.')

            else:
                length = len(index)

                if data == fill_value or (isna(data) and isna(fill_value)):
                    if kind == 'block':
                        sparse_index = BlockIndex(length, [], [])
                    else:
                        sparse_index = IntIndex(length, [])
                    data = np.array([])

                else:
                    if kind == 'block':
                        locs, lens = ([0], [length]) if length else ([], [])
                        sparse_index = BlockIndex(length, locs, lens)
                    else:
                        sparse_index = IntIndex(length, index)
                    v = data
                    data = np.empty(length)
                    data.fill(v)

            if index is None:
                index = com._default_index(sparse_index.length)
            index = _ensure_index(index)

            # create/copy the manager
            if isinstance(data, SingleBlockManager):

                if copy:
                    data = data.copy()
            else:

                # create a sparse array
                if not isinstance(data, SparseArray):
                    data = SparseArray(data, sparse_index=sparse_index,
                                       fill_value=fill_value, dtype=dtype,
                                       copy=copy)

                data = SingleBlockManager(data, index)

        generic.NDFrame.__init__(self, data)

        self.index = index
        self.name = name

    @property
    def values(self):
        """ return the array """
        return self.block.values

    def __array__(self, result=None):
        """ the array interface, return my values """
        return self.block.values

    def get_values(self):
        """ same as values """
        return self.block.to_dense().view()

    @property
    def block(self):
        return self._data._block

    @property
    def fill_value(self):
        return self.block.fill_value

    @fill_value.setter
    def fill_value(self, v):
        self.block.fill_value = v

    @property
    def sp_index(self):
        return self.block.sp_index

    @property
    def sp_values(self):
        return self.values.sp_values

    @property
    def npoints(self):
        return self.sp_index.npoints

    @classmethod
    def from_array(cls, arr, index=None, name=None, copy=False,
                   fill_value=None, fastpath=False):
        """Construct SparseSeries from array.

        .. deprecated:: 0.23.0
            Use the pd.SparseSeries(..) constructor instead.
        """
        warnings.warn("'from_array' is deprecated and will be removed in a "
                      "future version. Please use the pd.SparseSeries(..) "
                      "constructor instead.", FutureWarning, stacklevel=2)
        return cls(arr, index=index, name=name, copy=copy,
                   fill_value=fill_value, fastpath=fastpath)

    @property
    def _constructor(self):
        return SparseSeries

    @property
    def _constructor_expanddim(self):
        from pandas.core.sparse.api import SparseDataFrame
        return SparseDataFrame

    @property
    def kind(self):
        if isinstance(self.sp_index, BlockIndex):
            return 'block'
        elif isinstance(self.sp_index, IntIndex):
            return 'integer'

    def as_sparse_array(self, kind=None, fill_value=None, copy=False):
        """ return my self as a sparse array, do not copy by default """

        if fill_value is None:
            fill_value = self.fill_value
        if kind is None:
            kind = self.kind
        return SparseArray(self.values, sparse_index=self.sp_index,
                           fill_value=fill_value, kind=kind, copy=copy)

    def __len__(self):
        return len(self.block)

    @property
    def shape(self):
        return self._data.shape

    def __unicode__(self):
        # currently, unicode is same as repr...fixes infinite loop
        series_rep = Series.__unicode__(self)
        rep = '{series}\n{index!r}'.format(series=series_rep,
                                           index=self.sp_index)
        return rep

    def __array_wrap__(self, result, context=None):
        """
        Gets called prior to a ufunc (and after)

        See SparseArray.__array_wrap__ for detail.
        """
        if isinstance(context, tuple) and len(context) == 3:
            ufunc, args, domain = context
            args = [getattr(a, 'fill_value', a) for a in args]
            with np.errstate(all='ignore'):
                fill_value = ufunc(self.fill_value, *args[1:])
        else:
            fill_value = self.fill_value

        return self._constructor(result, index=self.index,
                                 sparse_index=self.sp_index,
                                 fill_value=fill_value,
                                 copy=False).__finalize__(self)

    def __array_finalize__(self, obj):
        """
        Gets called after any ufunc or other array operations, necessary
        to pass on the index.
        """
        self.name = getattr(obj, 'name', None)
        self.fill_value = getattr(obj, 'fill_value', None)

    def _reduce(self, op, name, axis=0, skipna=True, numeric_only=None,
                filter_type=None, **kwds):
        """ perform a reduction operation """
        return op(self.get_values(), skipna=skipna, **kwds)

    def __getstate__(self):
        # pickling
        return dict(_typ=self._typ, _subtyp=self._subtyp, _data=self._data,
                    fill_value=self.fill_value, name=self.name)

    def _unpickle_series_compat(self, state):

        nd_state, own_state = state

        # recreate the ndarray
        data = np.empty(nd_state[1], dtype=nd_state[2])
        np.ndarray.__setstate__(data, nd_state)

        index, fill_value, sp_index = own_state[:3]
        name = None
        if len(own_state) > 3:
            name = own_state[3]

        # create a sparse array
        if not isinstance(data, SparseArray):
            data = SparseArray(data, sparse_index=sp_index,
                               fill_value=fill_value, copy=False)

        # recreate
        data = SingleBlockManager(data, index, fastpath=True)
        generic.NDFrame.__init__(self, data)

        self._set_axis(0, index)
        self.name = name

    def __iter__(self):
        """ forward to the array """
        return iter(self.values)

    def _set_subtyp(self, is_all_dates):
        if is_all_dates:
            object.__setattr__(self, '_subtyp', 'sparse_time_series')
        else:
            object.__setattr__(self, '_subtyp', 'sparse_series')

    def _ixs(self, i, axis=0):
        """
        Return the i-th value or values in the SparseSeries by location

        Parameters
        ----------
        i : int, slice, or sequence of integers

        Returns
        -------
        value : scalar (int) or Series (slice, sequence)
        """
        label = self.index[i]
        if isinstance(label, Index):
            return self.take(i, axis=axis)
        else:
            return self._get_val_at(i)

    def _get_val_at(self, loc):
        """ forward to the array """
        return self.block.values._get_val_at(loc)

    def __getitem__(self, key):
        try:
            return self.index.get_value(self, key)

        except InvalidIndexError:
            pass
        except KeyError:
            if isinstance(key, (int, np.integer)):
                return self._get_val_at(key)
            elif key is Ellipsis:
                return self
            raise Exception('Requested index not in this series!')

        except TypeError:
            # Could not hash item, must be array-like?
            pass

        key = com._values_from_object(key)
        if self.index.nlevels > 1 and isinstance(key, tuple):
            # to handle MultiIndex labels
            key = self.index.get_loc(key)
        return self._constructor(self.values[key],
                                 index=self.index[key]).__finalize__(self)

    def _get_values(self, indexer):
        try:
            return self._constructor(self._data.get_slice(indexer),
                                     fastpath=True).__finalize__(self)
        except Exception:
            return self[indexer]

    def _set_with_engine(self, key, value):
        return self._set_value(key, value)

    def abs(self):
        """
        Return an object with absolute value taken. Only applicable to objects
        that are all numeric

        Returns
        -------
        abs: type of caller
        """
        return self._constructor(np.abs(self.values),
                                 index=self.index).__finalize__(self)

    def get(self, label, default=None):
        """
        Returns value occupying requested label, default to specified
        missing value if not present. Analogous to dict.get

        Parameters
        ----------
        label : object
            Label value looking for
        default : object, optional
            Value to return if label not in index

        Returns
        -------
        y : scalar
        """
        if label in self.index:
            loc = self.index.get_loc(label)
            return self._get_val_at(loc)
        else:
            return default

    def get_value(self, label, takeable=False):
        """
        Retrieve single value at passed index label

        .. deprecated:: 0.21.0

        Please use .at[] or .iat[] accessors.

        Parameters
        ----------
        index : label
        takeable : interpret the index as indexers, default False

        Returns
        -------
        value : scalar value
        """
        warnings.warn("get_value is deprecated and will be removed "
                      "in a future release. Please use "
                      ".at[] or .iat[] accessors instead", FutureWarning,
                      stacklevel=2)

        return self._get_value(label, takeable=takeable)

    def _get_value(self, label, takeable=False):
        loc = label if takeable is True else self.index.get_loc(label)
        return self._get_val_at(loc)
    _get_value.__doc__ = get_value.__doc__

    def set_value(self, label, value, takeable=False):
        """
        Quickly set single value at passed label. If label is not contained, a
        new object is created with the label placed at the end of the result
        index

        .. deprecated:: 0.21.0

        Please use .at[] or .iat[] accessors.

        Parameters
        ----------
        label : object
            Partial indexing with MultiIndex not allowed
        value : object
            Scalar value
        takeable : interpret the index as indexers, default False

        Notes
        -----
        This method *always* returns a new object. It is not particularly
        efficient but is provided for API compatibility with Series

        Returns
        -------
        series : SparseSeries
        """
        warnings.warn("set_value is deprecated and will be removed "
                      "in a future release. Please use "
                      ".at[] or .iat[] accessors instead", FutureWarning,
                      stacklevel=2)
        return self._set_value(label, value, takeable=takeable)

    def _set_value(self, label, value, takeable=False):
        values = self.to_dense()

        # if the label doesn't exist, we will create a new object here
        # and possibly change the index
        new_values = values._set_value(label, value, takeable=takeable)
        if new_values is not None:
            values = new_values
        new_index = values.index
        values = SparseArray(values, fill_value=self.fill_value,
                             kind=self.kind)
        self._data = SingleBlockManager(values, new_index)
        self._index = new_index
    _set_value.__doc__ = set_value.__doc__

    def _set_values(self, key, value):

        # this might be inefficient as we have to recreate the sparse array
        # rather than setting individual elements, but have to convert
        # the passed slice/boolean that's in dense space into a sparse indexer
        # not sure how to do that!
        if isinstance(key, Series):
            key = key.values

        values = self.values.to_dense()
        values[key] = libindex.convert_scalar(values, value)
        values = SparseArray(values, fill_value=self.fill_value,
                             kind=self.kind)
        self._data = SingleBlockManager(values, self.index)

    def to_dense(self, sparse_only=False):
        """
        Convert SparseSeries to a Series.

        Parameters
        ----------
        sparse_only : bool, default False
            .. deprecated:: 0.20.0
                This argument will be removed in a future version.

            If True, return just the non-sparse values, or the dense version
            of `self.values` if False.

        Returns
        -------
        s : Series
        """
        if sparse_only:
            warnings.warn(("The 'sparse_only' parameter has been deprecated "
                           "and will be removed in a future version."),
                          FutureWarning, stacklevel=2)
            int_index = self.sp_index.to_int_index()
            index = self.index.take(int_index.indices)
            return Series(self.sp_values, index=index, name=self.name)
        else:
            return Series(self.values.to_dense(), index=self.index,
                          name=self.name)

    @property
    def density(self):
        r = float(self.sp_index.npoints) / float(self.sp_index.length)
        return r

    def copy(self, deep=True):
        """
        Make a copy of the SparseSeries. Only the actual sparse values need to
        be copied
        """
        new_data = self._data
        if deep:
            new_data = self._data.copy()

        return self._constructor(new_data, sparse_index=self.sp_index,
                                 fill_value=self.fill_value).__finalize__(self)

    @Appender(generic._shared_docs['reindex'] % _shared_doc_kwargs)
    def reindex(self, index=None, method=None, copy=True, limit=None,
                **kwargs):

        return super(SparseSeries, self).reindex(index=index, method=method,
                                                 copy=copy, limit=limit,
                                                 **kwargs)

    def sparse_reindex(self, new_index):
        """
        Conform sparse values to new SparseIndex

        Parameters
        ----------
        new_index : {BlockIndex, IntIndex}

        Returns
        -------
        reindexed : SparseSeries
        """
        if not isinstance(new_index, splib.SparseIndex):
            raise TypeError('new index must be a SparseIndex')

        block = self.block.sparse_reindex(new_index)
        new_data = SingleBlockManager(block, self.index)
        return self._constructor(new_data, index=self.index,
                                 sparse_index=new_index,
                                 fill_value=self.fill_value).__finalize__(self)

    @Appender(generic._shared_docs['take'])
    def take(self, indices, axis=0, convert=None, *args, **kwargs):
        if convert is not None:
            msg = ("The 'convert' parameter is deprecated "
                   "and will be removed in a future version.")
            warnings.warn(msg, FutureWarning, stacklevel=2)
        else:
            convert = True

        nv.validate_take_with_convert(convert, args, kwargs)
        new_values = SparseArray.take(self.values, indices)
        new_index = self.index.take(indices)
        return self._constructor(new_values,
                                 index=new_index).__finalize__(self)

    def cumsum(self, axis=0, *args, **kwargs):
        """
        Cumulative sum of non-NA/null values.

        When performing the cumulative summation, any non-NA/null values will
        be skipped. The resulting SparseSeries will preserve the locations of
        NaN values, but the fill value will be `np.nan` regardless.

        Parameters
        ----------
        axis : {0}

        Returns
        -------
        cumsum : SparseSeries
        """
        nv.validate_cumsum(args, kwargs)
        if axis is not None:
            axis = self._get_axis_number(axis)

        new_array = self.values.cumsum()

        return self._constructor(
            new_array, index=self.index,
            sparse_index=new_array.sp_index).__finalize__(self)

    @Appender(generic._shared_docs['isna'] % _shared_doc_kwargs)
    def isna(self):
        arr = SparseArray(isna(self.values.sp_values),
                          sparse_index=self.values.sp_index,
                          fill_value=isna(self.fill_value))
        return self._constructor(arr, index=self.index).__finalize__(self)
    isnull = isna

    @Appender(generic._shared_docs['notna'] % _shared_doc_kwargs)
    def notna(self):
        arr = SparseArray(notna(self.values.sp_values),
                          sparse_index=self.values.sp_index,
                          fill_value=notna(self.fill_value))
        return self._constructor(arr, index=self.index).__finalize__(self)
    notnull = notna

    def dropna(self, axis=0, inplace=False, **kwargs):
        """
        Analogous to Series.dropna. If fill_value=NaN, returns a dense Series
        """
        # TODO: make more efficient
        axis = self._get_axis_number(axis or 0)
        dense_valid = self.to_dense().dropna()
        if inplace:
            raise NotImplementedError("Cannot perform inplace dropna"
                                      " operations on a SparseSeries")
        if isna(self.fill_value):
            return dense_valid
        else:
            dense_valid = dense_valid[dense_valid != self.fill_value]
            return dense_valid.to_sparse(fill_value=self.fill_value)

    @Appender(generic._shared_docs['shift'] % _shared_doc_kwargs)
    def shift(self, periods, freq=None, axis=0):
        if periods == 0:
            return self.copy()

        # no special handling of fill values yet
        if not isna(self.fill_value):
            shifted = self.to_dense().shift(periods, freq=freq,
                                            axis=axis)
            return shifted.to_sparse(fill_value=self.fill_value,
                                     kind=self.kind)

        if freq is not None:
            return self._constructor(
                self.sp_values, sparse_index=self.sp_index,
                index=self.index.shift(periods, freq),
                fill_value=self.fill_value).__finalize__(self)

        int_index = self.sp_index.to_int_index()
        new_indices = int_index.indices + periods
        start, end = new_indices.searchsorted([0, int_index.length])

        new_indices = new_indices[start:end]
        new_sp_index = _make_index(len(self), new_indices, self.sp_index)

        arr = self.values._simple_new(self.sp_values[start:end].copy(),
                                      new_sp_index, fill_value=np.nan)
        return self._constructor(arr, index=self.index).__finalize__(self)

    def combine_first(self, other):
        """
        Combine Series values, choosing the calling Series's values
        first. Result index will be the union of the two indexes

        Parameters
        ----------
        other : Series

        Returns
        -------
        y : Series
        """
        if isinstance(other, SparseSeries):
            other = other.to_dense()

        dense_combined = self.to_dense().combine_first(other)
        return dense_combined.to_sparse(fill_value=self.fill_value)

    def to_coo(self, row_levels=(0, ), column_levels=(1, ), sort_labels=False):
        """
        Create a scipy.sparse.coo_matrix from a SparseSeries with MultiIndex.

        Use row_levels and column_levels to determine the row and column
        coordinates respectively. row_levels and column_levels are the names
        (labels) or numbers of the levels. {row_levels, column_levels} must be
        a partition of the MultiIndex level names (or numbers).

        Parameters
        ----------
        row_levels : tuple/list
        column_levels : tuple/list
        sort_labels : bool, default False
            Sort the row and column labels before forming the sparse matrix.

        Returns
        -------
        y : scipy.sparse.coo_matrix
        rows : list (row labels)
        columns : list (column labels)

        Examples
        --------
        >>> from numpy import nan
        >>> s = Series([3.0, nan, 1.0, 3.0, nan, nan])
        >>> s.index = MultiIndex.from_tuples([(1, 2, 'a', 0),
                                              (1, 2, 'a', 1),
                                              (1, 1, 'b', 0),
                                              (1, 1, 'b', 1),
                                              (2, 1, 'b', 0),
                                              (2, 1, 'b', 1)],
                                              names=['A', 'B', 'C', 'D'])
        >>> ss = s.to_sparse()
        >>> A, rows, columns = ss.to_coo(row_levels=['A', 'B'],
                                         column_levels=['C', 'D'],
                                         sort_labels=True)
        >>> A
        <3x4 sparse matrix of type '<class 'numpy.float64'>'
                with 3 stored elements in COOrdinate format>
        >>> A.todense()
        matrix([[ 0.,  0.,  1.,  3.],
        [ 3.,  0.,  0.,  0.],
        [ 0.,  0.,  0.,  0.]])
        >>> rows
        [(1, 1), (1, 2), (2, 1)]
        >>> columns
        [('a', 0), ('a', 1), ('b', 0), ('b', 1)]
        """
        A, rows, columns = _sparse_series_to_coo(self, row_levels,
                                                 column_levels,
                                                 sort_labels=sort_labels)
        return A, rows, columns

    @classmethod
    def from_coo(cls, A, dense_index=False):
        """
        Create a SparseSeries from a scipy.sparse.coo_matrix.

        Parameters
        ----------
        A : scipy.sparse.coo_matrix
        dense_index : bool, default False
            If False (default), the SparseSeries index consists of only the
            coords of the non-null entries of the original coo_matrix.
            If True, the SparseSeries index consists of the full sorted
            (row, col) coordinates of the coo_matrix.

        Returns
        -------
        s : SparseSeries

        Examples
        ---------
        >>> from scipy import sparse
        >>> A = sparse.coo_matrix(([3.0, 1.0, 2.0], ([1, 0, 0], [0, 2, 3])),
                               shape=(3, 4))
        >>> A
        <3x4 sparse matrix of type '<class 'numpy.float64'>'
                with 3 stored elements in COOrdinate format>
        >>> A.todense()
        matrix([[ 0.,  0.,  1.,  2.],
                [ 3.,  0.,  0.,  0.],
                [ 0.,  0.,  0.,  0.]])
        >>> ss = SparseSeries.from_coo(A)
        >>> ss
        0  2    1
           3    2
        1  0    3
        dtype: float64
        BlockIndex
        Block locations: array([0], dtype=int32)
        Block lengths: array([3], dtype=int32)
        """
        return _coo_to_sparse_series(A, dense_index=dense_index)


# overwrite series methods with unaccelerated Sparse-specific versions
ops.add_flex_arithmetic_methods(SparseSeries)
ops.add_special_arithmetic_methods(SparseSeries)
