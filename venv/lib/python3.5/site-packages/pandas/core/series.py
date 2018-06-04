"""
Data structure for 1-dimensional cross-sectional and time series data
"""
from __future__ import division

# pylint: disable=E1101,E1103
# pylint: disable=W0703,W0622,W0613,W0201

import types
import warnings
from textwrap import dedent

import numpy as np
import numpy.ma as ma

from pandas.core.accessor import CachedAccessor
from pandas.core.arrays import ExtensionArray
from pandas.core.dtypes.common import (
    is_categorical_dtype,
    is_bool,
    is_integer, is_integer_dtype,
    is_float_dtype,
    is_extension_type,
    is_extension_array_dtype,
    is_datetime64tz_dtype,
    is_timedelta64_dtype,
    is_object_dtype,
    is_list_like,
    is_hashable,
    is_iterator,
    is_dict_like,
    is_scalar,
    _is_unorderable_exception,
    _ensure_platform_int,
    pandas_dtype)
from pandas.core.dtypes.generic import (
    ABCSparseArray, ABCDataFrame, ABCIndexClass)
from pandas.core.dtypes.cast import (
    maybe_upcast, infer_dtype_from_scalar,
    maybe_convert_platform,
    maybe_cast_to_datetime, maybe_castable,
    construct_1d_arraylike_from_scalar,
    construct_1d_object_array_from_listlike)
from pandas.core.dtypes.missing import (
    isna,
    notna,
    remove_na_arraylike,
    na_value_for_dtype)

from pandas.core.index import (Index, MultiIndex, InvalidIndexError,
                               Float64Index, _ensure_index)
from pandas.core.indexing import check_bool_indexer, maybe_convert_indices
from pandas.core import generic, base
from pandas.core.internals import SingleBlockManager
from pandas.core.arrays.categorical import Categorical, CategoricalAccessor
from pandas.core.indexes.accessors import CombinedDatetimelikeProperties
from pandas.core.indexes.datetimes import DatetimeIndex
from pandas.core.indexes.timedeltas import TimedeltaIndex
from pandas.core.indexes.period import PeriodIndex
from pandas import compat
from pandas.io.formats.terminal import get_terminal_size
from pandas.compat import (
    zip, u, OrderedDict, StringIO, range, get_range_parameters, PY36)
from pandas.compat.numpy import function as nv

import pandas.core.ops as ops
import pandas.core.algorithms as algorithms

import pandas.core.common as com
import pandas.core.nanops as nanops
import pandas.io.formats.format as fmt
from pandas.util._decorators import (
    Appender, deprecate, deprecate_kwarg, Substitution)
from pandas.util._validators import validate_bool_kwarg

from pandas._libs import index as libindex, tslib as libts, lib, iNaT
from pandas.core.config import get_option
from pandas.core.strings import StringMethods

import pandas.plotting._core as gfx

__all__ = ['Series']

_shared_doc_kwargs = dict(
    axes='index', klass='Series', axes_single_arg="{0 or 'index'}",
    axis="""
    axis : {0 or 'index'}
        Parameter needed for compatibility with DataFrame.
    """,
    inplace="""inplace : boolean, default False
        If True, performs operation inplace and returns None.""",
    unique='np.ndarray', duplicated='Series',
    optional_by='', optional_mapper='', optional_labels='', optional_axis='',
    versionadded_to_excel='\n    .. versionadded:: 0.20.0\n')


# see gh-16971
def remove_na(arr):
    """Remove null values from array like structure.

    .. deprecated:: 0.21.0
        Use s[s.notnull()] instead.
    """

    warnings.warn("remove_na is deprecated and is a private "
                  "function. Do not use.", FutureWarning, stacklevel=2)
    return remove_na_arraylike(arr)


def _coerce_method(converter):
    """ install the scalar coercion methods """

    def wrapper(self):
        if len(self) == 1:
            return converter(self.iloc[0])
        raise TypeError("cannot convert the series to "
                        "{0}".format(str(converter)))

    return wrapper

# ----------------------------------------------------------------------
# Series class


class Series(base.IndexOpsMixin, generic.NDFrame):
    """
    One-dimensional ndarray with axis labels (including time series).

    Labels need not be unique but must be a hashable type. The object
    supports both integer- and label-based indexing and provides a host of
    methods for performing operations involving the index. Statistical
    methods from ndarray have been overridden to automatically exclude
    missing data (currently represented as NaN).

    Operations between Series (+, -, /, *, **) align values based on their
    associated index values-- they need not be the same length. The result
    index will be the sorted union of the two indexes.

    Parameters
    ----------
    data : array-like, dict, or scalar value
        Contains data stored in Series

        .. versionchanged :: 0.23.0
           If data is a dict, argument order is maintained for Python 3.6
           and later.

    index : array-like or Index (1d)
        Values must be hashable and have the same length as `data`.
        Non-unique index values are allowed. Will default to
        RangeIndex (0, 1, 2, ..., n) if not provided. If both a dict and index
        sequence are used, the index will override the keys found in the
        dict.
    dtype : numpy.dtype or None
        If None, dtype will be inferred
    copy : boolean, default False
        Copy input data
    """
    _metadata = ['name']
    _accessors = set(['dt', 'cat', 'str'])
    _deprecations = generic.NDFrame._deprecations | frozenset(
        ['asobject', 'sortlevel', 'reshape', 'get_value', 'set_value',
         'from_csv', 'valid'])

    def __init__(self, data=None, index=None, dtype=None, name=None,
                 copy=False, fastpath=False):

        # we are called internally, so short-circuit
        if fastpath:

            # data is an ndarray, index is defined
            if not isinstance(data, SingleBlockManager):
                data = SingleBlockManager(data, index, fastpath=True)
            if copy:
                data = data.copy()
            if index is None:
                index = data.index

        else:

            if index is not None:
                index = _ensure_index(index)

            if data is None:
                data = {}
            if dtype is not None:
                dtype = self._validate_dtype(dtype)

            if isinstance(data, MultiIndex):
                raise NotImplementedError("initializing a Series from a "
                                          "MultiIndex is not supported")
            elif isinstance(data, Index):
                if name is None:
                    name = data.name

                if dtype is not None:
                    # astype copies
                    data = data.astype(dtype)
                else:
                    # need to copy to avoid aliasing issues
                    data = data._values.copy()
                copy = False

            elif isinstance(data, np.ndarray):
                pass
            elif isinstance(data, Series):
                if name is None:
                    name = data.name
                if index is None:
                    index = data.index
                else:
                    data = data.reindex(index, copy=copy)
                data = data._data
            elif isinstance(data, dict):
                data, index = self._init_dict(data, index, dtype)
                dtype = None
                copy = False
            elif isinstance(data, SingleBlockManager):
                if index is None:
                    index = data.index
                elif not data.index.equals(index) or copy:
                    # GH#19275 SingleBlockManager input should only be called
                    # internally
                    raise AssertionError('Cannot pass both SingleBlockManager '
                                         '`data` argument and a different '
                                         '`index` argument.  `copy` must '
                                         'be False.')

            elif is_extension_array_dtype(data) and dtype is not None:
                if not data.dtype.is_dtype(dtype):
                    raise ValueError("Cannot specify a dtype '{}' with an "
                                     "extension array of a different "
                                     "dtype ('{}').".format(dtype,
                                                            data.dtype))

            elif (isinstance(data, types.GeneratorType) or
                  (compat.PY3 and isinstance(data, map))):
                data = list(data)
            elif isinstance(data, (set, frozenset)):
                raise TypeError("{0!r} type is unordered"
                                "".format(data.__class__.__name__))
            else:

                # handle sparse passed here (and force conversion)
                if isinstance(data, ABCSparseArray):
                    data = data.to_dense()

            if index is None:
                if not is_list_like(data):
                    data = [data]
                index = com._default_index(len(data))
            elif is_list_like(data):

                # a scalar numpy array is list-like but doesn't
                # have a proper length
                try:
                    if len(index) != len(data):
                        raise ValueError(
                            'Length of passed values is {val}, '
                            'index implies {ind}'
                            .format(val=len(data), ind=len(index)))
                except TypeError:
                    pass

            # create/copy the manager
            if isinstance(data, SingleBlockManager):
                if dtype is not None:
                    data = data.astype(dtype=dtype, errors='ignore',
                                       copy=copy)
                elif copy:
                    data = data.copy()
            else:
                data = _sanitize_array(data, index, dtype, copy,
                                       raise_cast_failure=True)

                data = SingleBlockManager(data, index, fastpath=True)

        generic.NDFrame.__init__(self, data, fastpath=True)

        self.name = name
        self._set_axis(0, index, fastpath=True)

    def _init_dict(self, data, index=None, dtype=None):
        """
        Derive the "_data" and "index" attributes of a new Series from a
        dictionary input.

        Parameters
        ----------
        data : dict or dict-like
            Data used to populate the new Series
        index : Index or index-like, default None
            index for the new Series: if None, use dict keys
        dtype : dtype, default None
            dtype for the new Series: if None, infer from data

        Returns
        -------
        _data : BlockManager for the new Series
        index : index for the new Series
        """
        # Looking for NaN in dict doesn't work ({np.nan : 1}[float('nan')]
        # raises KeyError), so we iterate the entire dict, and align
        if data:
            keys, values = zip(*compat.iteritems(data))
            values = list(values)
        elif index is not None:
            # fastpath for Series(data=None). Just use broadcasting a scalar
            # instead of reindexing.
            values = na_value_for_dtype(dtype)
            keys = index
        else:
            keys, values = [], []

        # Input is now list-like, so rely on "standard" construction:
        s = Series(values, index=keys, dtype=dtype)

        # Now we just make sure the order is respected, if any
        if data and index is not None:
            s = s.reindex(index, copy=False)
        elif not PY36 and not isinstance(data, OrderedDict) and data:
            # Need the `and data` to avoid sorting Series(None, index=[...])
            # since that isn't really dict-like
            try:
                s = s.sort_index()
            except TypeError:
                pass
        return s._data, s.index

    @classmethod
    def from_array(cls, arr, index=None, name=None, dtype=None, copy=False,
                   fastpath=False):
        """Construct Series from array.

        .. deprecated :: 0.23.0
            Use pd.Series(..) constructor instead.

        """
        warnings.warn("'from_array' is deprecated and will be removed in a "
                      "future version. Please use the pd.Series(..) "
                      "constructor instead.", FutureWarning, stacklevel=2)
        if isinstance(arr, ABCSparseArray):
            from pandas.core.sparse.series import SparseSeries
            cls = SparseSeries
        return cls(arr, index=index, name=name, dtype=dtype,
                   copy=copy, fastpath=fastpath)

    @property
    def _constructor(self):
        return Series

    @property
    def _constructor_expanddim(self):
        from pandas.core.frame import DataFrame
        return DataFrame

    # types
    @property
    def _can_hold_na(self):
        return self._data._can_hold_na

    _index = None

    def _set_axis(self, axis, labels, fastpath=False):
        """ override generic, we want to set the _typ here """

        if not fastpath:
            labels = _ensure_index(labels)

        is_all_dates = labels.is_all_dates
        if is_all_dates:
            if not isinstance(labels,
                              (DatetimeIndex, PeriodIndex, TimedeltaIndex)):
                try:
                    labels = DatetimeIndex(labels)
                    # need to set here because we changed the index
                    if fastpath:
                        self._data.set_axis(axis, labels)
                except (libts.OutOfBoundsDatetime, ValueError):
                    # labels may exceeds datetime bounds,
                    # or not be a DatetimeIndex
                    pass

        self._set_subtyp(is_all_dates)

        object.__setattr__(self, '_index', labels)
        if not fastpath:
            self._data.set_axis(axis, labels)

    def _set_subtyp(self, is_all_dates):
        if is_all_dates:
            object.__setattr__(self, '_subtyp', 'time_series')
        else:
            object.__setattr__(self, '_subtyp', 'series')

    def _update_inplace(self, result, **kwargs):
        # we want to call the generic version and not the IndexOpsMixin
        return generic.NDFrame._update_inplace(self, result, **kwargs)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        if value is not None and not is_hashable(value):
            raise TypeError('Series.name must be a hashable type')
        object.__setattr__(self, '_name', value)

    # ndarray compatibility
    @property
    def dtype(self):
        """ return the dtype object of the underlying data """
        return self._data.dtype

    @property
    def dtypes(self):
        """ return the dtype object of the underlying data """
        return self._data.dtype

    @property
    def ftype(self):
        """ return if the data is sparse|dense """
        return self._data.ftype

    @property
    def ftypes(self):
        """ return if the data is sparse|dense """
        return self._data.ftype

    @property
    def values(self):
        """
        Return Series as ndarray or ndarray-like
        depending on the dtype

        Returns
        -------
        arr : numpy.ndarray or ndarray-like

        Examples
        --------
        >>> pd.Series([1, 2, 3]).values
        array([1, 2, 3])

        >>> pd.Series(list('aabc')).values
        array(['a', 'a', 'b', 'c'], dtype=object)

        >>> pd.Series(list('aabc')).astype('category').values
        [a, a, b, c]
        Categories (3, object): [a, b, c]

        Timezone aware datetime data is converted to UTC:

        >>> pd.Series(pd.date_range('20130101', periods=3,
        ...                         tz='US/Eastern')).values
        array(['2013-01-01T05:00:00.000000000',
               '2013-01-02T05:00:00.000000000',
               '2013-01-03T05:00:00.000000000'], dtype='datetime64[ns]')

        """
        return self._data.external_values()

    @property
    def _values(self):
        """ return the internal repr of this data """
        return self._data.internal_values()

    def _formatting_values(self):
        """Return the values that can be formatted (used by SeriesFormatter
        and DataFrameFormatter)
        """
        return self._data.formatting_values()

    def get_values(self):
        """ same as values (but handles sparseness conversions); is a view """
        return self._data.get_values()

    @property
    def asobject(self):
        """Return object Series which contains boxed values.

        .. deprecated :: 0.23.0

           Use ``astype(object)`` instead.

        *this is an internal non-public method*
        """
        warnings.warn("'asobject' is deprecated. Use 'astype(object)'"
                      " instead", FutureWarning, stacklevel=2)
        return self.astype(object).values

    # ops
    def ravel(self, order='C'):
        """
        Return the flattened underlying data as an ndarray

        See also
        --------
        numpy.ndarray.ravel
        """
        return self._values.ravel(order=order)

    def compress(self, condition, *args, **kwargs):
        """
        Return selected slices of an array along given axis as a Series

        See also
        --------
        numpy.ndarray.compress
        """
        nv.validate_compress(args, kwargs)
        return self[condition]

    def nonzero(self):
        """
        Return the *integer* indices of the elements that are non-zero

        This method is equivalent to calling `numpy.nonzero` on the
        series data. For compatibility with NumPy, the return value is
        the same (a tuple with an array of indices for each dimension),
        but it will always be a one-item tuple because series only have
        one dimension.

        Examples
        --------
        >>> s = pd.Series([0, 3, 0, 4])
        >>> s.nonzero()
        (array([1, 3]),)
        >>> s.iloc[s.nonzero()[0]]
        1    3
        3    4
        dtype: int64

        >>> s = pd.Series([0, 3, 0, 4], index=['a', 'b', 'c', 'd'])
        # same return although index of s is different
        >>> s.nonzero()
        (array([1, 3]),)
        >>> s.iloc[s.nonzero()[0]]
        b    3
        d    4
        dtype: int64

        See Also
        --------
        numpy.nonzero
        """
        return self._values.nonzero()

    def put(self, *args, **kwargs):
        """
        Applies the `put` method to its `values` attribute
        if it has one.

        See also
        --------
        numpy.ndarray.put
        """
        self._values.put(*args, **kwargs)

    def __len__(self):
        """
        return the length of the Series
        """
        return len(self._data)

    def view(self, dtype=None):
        """
        Create a new view of the Series.

        This function will return a new Series with a view of the same
        underlying values in memory, optionally reinterpreted with a new data
        type. The new data type must preserve the same size in bytes as to not
        cause index misalignment.

        Parameters
        ----------
        dtype : data type
            Data type object or one of their string representations.

        Returns
        -------
        Series
            A new Series object as a view of the same data in memory.

        See Also
        --------
        numpy.ndarray.view : Equivalent numpy function to create a new view of
            the same data in memory.

        Notes
        -----
        Series are instantiated with ``dtype=float64`` by default. While
        ``numpy.ndarray.view()`` will return a view with the same data type as
        the original array, ``Series.view()`` (without specified dtype)
        will try using ``float64`` and may fail if the original data type size
        in bytes is not the same.

        Examples
        --------
        >>> s = pd.Series([-2, -1, 0, 1, 2], dtype='int8')
        >>> s
        0   -2
        1   -1
        2    0
        3    1
        4    2
        dtype: int8

        The 8 bit signed integer representation of `-1` is `0b11111111`, but
        the same bytes represent 255 if read as an 8 bit unsigned integer:

        >>> us = s.view('uint8')
        >>> us
        0    254
        1    255
        2      0
        3      1
        4      2
        dtype: uint8

        The views share the same underlying values:

        >>> us[0] = 128
        >>> s
        0   -128
        1     -1
        2      0
        3      1
        4      2
        dtype: int8
        """
        return self._constructor(self._values.view(dtype),
                                 index=self.index).__finalize__(self)

    def __array__(self, result=None):
        """
        the array interface, return my values
        """
        return self.get_values()

    def __array_wrap__(self, result, context=None):
        """
        Gets called after a ufunc
        """
        return self._constructor(result, index=self.index,
                                 copy=False).__finalize__(self)

    def __array_prepare__(self, result, context=None):
        """
        Gets called prior to a ufunc
        """

        # nice error message for non-ufunc types
        if context is not None and not isinstance(self._values, np.ndarray):
            obj = context[1][0]
            raise TypeError("{obj} with dtype {dtype} cannot perform "
                            "the numpy op {op}".format(
                                obj=type(obj).__name__,
                                dtype=getattr(obj, 'dtype', None),
                                op=context[0].__name__))
        return result

    # complex
    @property
    def real(self):
        return self.values.real

    @real.setter
    def real(self, v):
        self.values.real = v

    @property
    def imag(self):
        return self.values.imag

    @imag.setter
    def imag(self, v):
        self.values.imag = v

    # coercion
    __float__ = _coerce_method(float)
    __long__ = _coerce_method(int)
    __int__ = _coerce_method(int)

    def _unpickle_series_compat(self, state):
        if isinstance(state, dict):
            self._data = state['_data']
            self.name = state['name']
            self.index = self._data.index

        elif isinstance(state, tuple):

            # < 0.12 series pickle

            nd_state, own_state = state

            # recreate the ndarray
            data = np.empty(nd_state[1], dtype=nd_state[2])
            np.ndarray.__setstate__(data, nd_state)

            # backwards compat
            index, name = own_state[0], None
            if len(own_state) > 1:
                name = own_state[1]

            # recreate
            self._data = SingleBlockManager(data, index, fastpath=True)
            self._index = index
            self.name = name

        else:
            raise Exception("cannot unpickle legacy formats -> [%s]" % state)

    # indexers
    @property
    def axes(self):
        """Return a list of the row axis labels"""
        return [self.index]

    def _ixs(self, i, axis=0):
        """
        Return the i-th value or values in the Series by location

        Parameters
        ----------
        i : int, slice, or sequence of integers

        Returns
        -------
        value : scalar (int) or Series (slice, sequence)
        """
        try:

            # dispatch to the values if we need
            values = self._values
            if isinstance(values, np.ndarray):
                return libindex.get_value_at(values, i)
            else:
                return values[i]
        except IndexError:
            raise
        except Exception:
            if isinstance(i, slice):
                indexer = self.index._convert_slice_indexer(i, kind='iloc')
                return self._get_values(indexer)
            else:
                label = self.index[i]
                if isinstance(label, Index):
                    return self.take(i, axis=axis, convert=True)
                else:
                    return libindex.get_value_at(self, i)

    @property
    def _is_mixed_type(self):
        return False

    def _slice(self, slobj, axis=0, kind=None):
        slobj = self.index._convert_slice_indexer(slobj,
                                                  kind=kind or 'getitem')
        return self._get_values(slobj)

    def __getitem__(self, key):
        key = com._apply_if_callable(key, self)
        try:
            result = self.index.get_value(self, key)

            if not is_scalar(result):
                if is_list_like(result) and not isinstance(result, Series):

                    # we need to box if loc of the key isn't scalar here
                    # otherwise have inline ndarray/lists
                    try:
                        if not is_scalar(self.index.get_loc(key)):
                            result = self._constructor(
                                result, index=[key] * len(result),
                                dtype=self.dtype).__finalize__(self)
                    except KeyError:
                        pass
            return result
        except InvalidIndexError:
            pass
        except (KeyError, ValueError):
            if isinstance(key, tuple) and isinstance(self.index, MultiIndex):
                # kludge
                pass
            elif key is Ellipsis:
                return self
            elif com.is_bool_indexer(key):
                pass
            else:

                # we can try to coerce the indexer (or this will raise)
                new_key = self.index._convert_scalar_indexer(key,
                                                             kind='getitem')
                if type(new_key) != type(key):
                    return self.__getitem__(new_key)
                raise

        except Exception:
            raise

        if is_iterator(key):
            key = list(key)

        if com.is_bool_indexer(key):
            key = check_bool_indexer(self.index, key)

        return self._get_with(key)

    def _get_with(self, key):
        # other: fancy integer or otherwise
        if isinstance(key, slice):
            indexer = self.index._convert_slice_indexer(key, kind='getitem')
            return self._get_values(indexer)
        elif isinstance(key, ABCDataFrame):
            raise TypeError('Indexing a Series with DataFrame is not '
                            'supported, use the appropriate DataFrame column')
        else:
            if isinstance(key, tuple):
                try:
                    return self._get_values_tuple(key)
                except Exception:
                    if len(key) == 1:
                        key = key[0]
                        if isinstance(key, slice):
                            return self._get_values(key)
                    raise

            # pragma: no cover
            if not isinstance(key, (list, np.ndarray, Series, Index)):
                key = list(key)

            if isinstance(key, Index):
                key_type = key.inferred_type
            else:
                key_type = lib.infer_dtype(key)

            if key_type == 'integer':
                if self.index.is_integer() or self.index.is_floating():
                    return self.loc[key]
                else:
                    return self._get_values(key)
            elif key_type == 'boolean':
                return self._get_values(key)
            else:
                try:
                    # handle the dup indexing case (GH 4246)
                    if isinstance(key, (list, tuple)):
                        return self.loc[key]

                    return self.reindex(key)
                except Exception:
                    # [slice(0, 5, None)] will break if you convert to ndarray,
                    # e.g. as requested by np.median
                    # hack
                    if isinstance(key[0], slice):
                        return self._get_values(key)
                    raise

    def _get_values_tuple(self, key):
        # mpl hackaround
        if com._any_none(*key):
            return self._get_values(key)

        if not isinstance(self.index, MultiIndex):
            raise ValueError('Can only tuple-index with a MultiIndex')

        # If key is contained, would have returned by now
        indexer, new_index = self.index.get_loc_level(key)
        return self._constructor(self._values[indexer],
                                 index=new_index).__finalize__(self)

    def _get_values(self, indexer):
        try:
            return self._constructor(self._data.get_slice(indexer),
                                     fastpath=True).__finalize__(self)
        except Exception:
            return self._values[indexer]

    def __setitem__(self, key, value):
        key = com._apply_if_callable(key, self)

        def setitem(key, value):
            try:
                self._set_with_engine(key, value)
                return
            except com.SettingWithCopyError:
                raise
            except (KeyError, ValueError):
                values = self._values
                if (is_integer(key) and
                        not self.index.inferred_type == 'integer'):

                    values[key] = value
                    return
                elif key is Ellipsis:
                    self[:] = value
                    return
                elif com.is_bool_indexer(key):
                    pass
                elif is_timedelta64_dtype(self.dtype):
                    # reassign a null value to iNaT
                    if isna(value):
                        value = iNaT

                        try:
                            self.index._engine.set_value(self._values, key,
                                                         value)
                            return
                        except TypeError:
                            pass

                self.loc[key] = value
                return

            except TypeError as e:
                if (isinstance(key, tuple) and
                        not isinstance(self.index, MultiIndex)):
                    raise ValueError("Can only tuple-index with a MultiIndex")

                # python 3 type errors should be raised
                if _is_unorderable_exception(e):
                    raise IndexError(key)

            if com.is_bool_indexer(key):
                key = check_bool_indexer(self.index, key)
                try:
                    self._where(~key, value, inplace=True)
                    return
                except InvalidIndexError:
                    pass

            self._set_with(key, value)

        # do the setitem
        cacher_needs_updating = self._check_is_chained_assignment_possible()
        setitem(key, value)
        if cacher_needs_updating:
            self._maybe_update_cacher()

    def _set_with_engine(self, key, value):
        values = self._values
        try:
            self.index._engine.set_value(values, key, value)
            return
        except KeyError:
            values[self.index.get_loc(key)] = value
            return

    def _set_with(self, key, value):
        # other: fancy integer or otherwise
        if isinstance(key, slice):
            indexer = self.index._convert_slice_indexer(key, kind='getitem')
            return self._set_values(indexer, value)
        else:
            if isinstance(key, tuple):
                try:
                    self._set_values(key, value)
                except Exception:
                    pass

            if not isinstance(key, (list, Series, np.ndarray, Series)):
                try:
                    key = list(key)
                except Exception:
                    key = [key]

            if isinstance(key, Index):
                key_type = key.inferred_type
            else:
                key_type = lib.infer_dtype(key)

            if key_type == 'integer':
                if self.index.inferred_type == 'integer':
                    self._set_labels(key, value)
                else:
                    return self._set_values(key, value)
            elif key_type == 'boolean':
                self._set_values(key.astype(np.bool_), value)
            else:
                self._set_labels(key, value)

    def _set_labels(self, key, value):
        if isinstance(key, Index):
            key = key.values
        else:
            key = com._asarray_tuplesafe(key)
        indexer = self.index.get_indexer(key)
        mask = indexer == -1
        if mask.any():
            raise ValueError('%s not contained in the index' % str(key[mask]))
        self._set_values(indexer, value)

    def _set_values(self, key, value):
        if isinstance(key, Series):
            key = key._values
        self._data = self._data.setitem(indexer=key, value=value)
        self._maybe_update_cacher()

    @deprecate_kwarg(old_arg_name='reps', new_arg_name='repeats')
    def repeat(self, repeats, *args, **kwargs):
        """
        Repeat elements of an Series. Refer to `numpy.ndarray.repeat`
        for more information about the `repeats` argument.

        See also
        --------
        numpy.ndarray.repeat
        """
        nv.validate_repeat(args, kwargs)
        new_index = self.index.repeat(repeats)
        new_values = self._values.repeat(repeats)
        return self._constructor(new_values,
                                 index=new_index).__finalize__(self)

    def get_value(self, label, takeable=False):
        """Quickly retrieve single value at passed index label

        .. deprecated:: 0.21.0
            Please use .at[] or .iat[] accessors.

        Parameters
        ----------
        label : object
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
        if takeable is True:
            return com._maybe_box_datetimelike(self._values[label])
        return self.index.get_value(self._values, label)
    _get_value.__doc__ = get_value.__doc__

    def set_value(self, label, value, takeable=False):
        """Quickly set single value at passed label. If label is not contained,
        a new object is created with the label placed at the end of the result
        index.

        .. deprecated:: 0.21.0
            Please use .at[] or .iat[] accessors.

        Parameters
        ----------
        label : object
            Partial indexing with MultiIndex not allowed
        value : object
            Scalar value
        takeable : interpret the index as indexers, default False

        Returns
        -------
        series : Series
            If label is contained, will be reference to calling Series,
            otherwise a new object
        """
        warnings.warn("set_value is deprecated and will be removed "
                      "in a future release. Please use "
                      ".at[] or .iat[] accessors instead", FutureWarning,
                      stacklevel=2)
        return self._set_value(label, value, takeable=takeable)

    def _set_value(self, label, value, takeable=False):
        try:
            if takeable:
                self._values[label] = value
            else:
                self.index._engine.set_value(self._values, label, value)
        except KeyError:

            # set using a non-recursive method
            self.loc[label] = value

        return self
    _set_value.__doc__ = set_value.__doc__

    def reset_index(self, level=None, drop=False, name=None, inplace=False):
        """
        Generate a new DataFrame or Series with the index reset.

        This is useful when the index needs to be treated as a column, or
        when the index is meaningless and needs to be reset to the default
        before another operation.

        Parameters
        ----------
        level : int, str, tuple, or list, default optional
            For a Series with a MultiIndex, only remove the specified levels
            from the index. Removes all levels by default.
        drop : bool, default False
            Just reset the index, without inserting it as a column in
            the new DataFrame.
        name : object, optional
            The name to use for the column containing the original Series
            values. Uses ``self.name`` by default. This argument is ignored
            when `drop` is True.
        inplace : bool, default False
            Modify the Series in place (do not create a new object).

        Returns
        -------
        Series or DataFrame
            When `drop` is False (the default), a DataFrame is returned.
            The newly created columns will come first in the DataFrame,
            followed by the original Series values.
            When `drop` is True, a `Series` is returned.
            In either case, if ``inplace=True``, no value is returned.

        See Also
        --------
        DataFrame.reset_index: Analogous function for DataFrame.

        Examples
        --------

        >>> s = pd.Series([1, 2, 3, 4], name='foo',
        ...               index=pd.Index(['a', 'b', 'c', 'd'], name='idx'))

        Generate a DataFrame with default index.

        >>> s.reset_index()
          idx  foo
        0   a    1
        1   b    2
        2   c    3
        3   d    4

        To specify the name of the new column use `name`.

        >>> s.reset_index(name='values')
          idx  values
        0   a       1
        1   b       2
        2   c       3
        3   d       4

        To generate a new Series with the default set `drop` to True.

        >>> s.reset_index(drop=True)
        0    1
        1    2
        2    3
        3    4
        Name: foo, dtype: int64

        To update the Series in place, without generating a new one
        set `inplace` to True. Note that it also requires ``drop=True``.

        >>> s.reset_index(inplace=True, drop=True)
        >>> s
        0    1
        1    2
        2    3
        3    4
        Name: foo, dtype: int64

        The `level` parameter is interesting for Series with a multi-level
        index.

        >>> arrays = [np.array(['bar', 'bar', 'baz', 'baz']),
        ...           np.array(['one', 'two', 'one', 'two'])]
        >>> s2 = pd.Series(
        ...     range(4), name='foo',
        ...     index=pd.MultiIndex.from_arrays(arrays,
        ...                                     names=['a', 'b']))

        To remove a specific level from the Index, use `level`.

        >>> s2.reset_index(level='a')
               a  foo
        b
        one  bar    0
        two  bar    1
        one  baz    2
        two  baz    3

        If `level` is not set, all levels are removed from the Index.

        >>> s2.reset_index()
             a    b  foo
        0  bar  one    0
        1  bar  two    1
        2  baz  one    2
        3  baz  two    3
        """
        inplace = validate_bool_kwarg(inplace, 'inplace')
        if drop:
            new_index = com._default_index(len(self))
            if level is not None and isinstance(self.index, MultiIndex):
                if not isinstance(level, (tuple, list)):
                    level = [level]
                level = [self.index._get_level_number(lev) for lev in level]
                if len(level) < len(self.index.levels):
                    new_index = self.index.droplevel(level)

            if inplace:
                self.index = new_index
                # set name if it was passed, otherwise, keep the previous name
                self.name = name or self.name
            else:
                return self._constructor(self._values.copy(),
                                         index=new_index).__finalize__(self)
        elif inplace:
            raise TypeError('Cannot reset_index inplace on a Series '
                            'to create a DataFrame')
        else:
            df = self.to_frame(name)
            return df.reset_index(level=level, drop=drop)

    def __unicode__(self):
        """
        Return a string representation for a particular DataFrame

        Invoked by unicode(df) in py2 only. Yields a Unicode String in both
        py2/py3.
        """
        buf = StringIO(u(""))
        width, height = get_terminal_size()
        max_rows = (height if get_option("display.max_rows") == 0 else
                    get_option("display.max_rows"))
        show_dimensions = get_option("display.show_dimensions")

        self.to_string(buf=buf, name=self.name, dtype=self.dtype,
                       max_rows=max_rows, length=show_dimensions)
        result = buf.getvalue()

        return result

    def to_string(self, buf=None, na_rep='NaN', float_format=None, header=True,
                  index=True, length=False, dtype=False, name=False,
                  max_rows=None):
        """
        Render a string representation of the Series

        Parameters
        ----------
        buf : StringIO-like, optional
            buffer to write to
        na_rep : string, optional
            string representation of NAN to use, default 'NaN'
        float_format : one-parameter function, optional
            formatter function to apply to columns' elements if they are floats
            default None
        header: boolean, default True
            Add the Series header (index name)
        index : bool, optional
            Add index (row) labels, default True
        length : boolean, default False
            Add the Series length
        dtype : boolean, default False
            Add the Series dtype
        name : boolean, default False
            Add the Series name if not None
        max_rows : int, optional
            Maximum number of rows to show before truncating. If None, show
            all.

        Returns
        -------
        formatted : string (if not buffer passed)
        """

        formatter = fmt.SeriesFormatter(self, name=name, length=length,
                                        header=header, index=index,
                                        dtype=dtype, na_rep=na_rep,
                                        float_format=float_format,
                                        max_rows=max_rows)
        result = formatter.to_string()

        # catch contract violations
        if not isinstance(result, compat.text_type):
            raise AssertionError("result must be of type unicode, type"
                                 " of result is {0!r}"
                                 "".format(result.__class__.__name__))

        if buf is None:
            return result
        else:
            try:
                buf.write(result)
            except AttributeError:
                with open(buf, 'w') as f:
                    f.write(result)

    def iteritems(self):
        """
        Lazily iterate over (index, value) tuples
        """
        return zip(iter(self.index), iter(self))

    items = iteritems

    # ----------------------------------------------------------------------
    # Misc public methods

    def keys(self):
        """Alias for index"""
        return self.index

    def to_dict(self, into=dict):
        """
        Convert Series to {label -> value} dict or dict-like object.

        Parameters
        ----------
        into : class, default dict
            The collections.Mapping subclass to use as the return
            object. Can be the actual class or an empty
            instance of the mapping type you want.  If you want a
            collections.defaultdict, you must pass it initialized.

            .. versionadded:: 0.21.0

        Returns
        -------
        value_dict : collections.Mapping

        Examples
        --------
        >>> s = pd.Series([1, 2, 3, 4])
        >>> s.to_dict()
        {0: 1, 1: 2, 2: 3, 3: 4}
        >>> from collections import OrderedDict, defaultdict
        >>> s.to_dict(OrderedDict)
        OrderedDict([(0, 1), (1, 2), (2, 3), (3, 4)])
        >>> dd = defaultdict(list)
        >>> s.to_dict(dd)
        defaultdict(<type 'list'>, {0: 1, 1: 2, 2: 3, 3: 4})
        """
        # GH16122
        into_c = com.standardize_mapping(into)
        return into_c(compat.iteritems(self))

    def to_frame(self, name=None):
        """
        Convert Series to DataFrame

        Parameters
        ----------
        name : object, default None
            The passed name should substitute for the series name (if it has
            one).

        Returns
        -------
        data_frame : DataFrame
        """
        if name is None:
            df = self._constructor_expanddim(self)
        else:
            df = self._constructor_expanddim({name: self})

        return df

    def to_sparse(self, kind='block', fill_value=None):
        """
        Convert Series to SparseSeries

        Parameters
        ----------
        kind : {'block', 'integer'}
        fill_value : float, defaults to NaN (missing)

        Returns
        -------
        sp : SparseSeries
        """
        from pandas.core.sparse.series import SparseSeries
        return SparseSeries(self, kind=kind,
                            fill_value=fill_value).__finalize__(self)

    def _set_name(self, name, inplace=False):
        """
        Set the Series name.

        Parameters
        ----------
        name : str
        inplace : bool
            whether to modify `self` directly or return a copy
        """
        inplace = validate_bool_kwarg(inplace, 'inplace')
        ser = self if inplace else self.copy()
        ser.name = name
        return ser

    # ----------------------------------------------------------------------
    # Statistics, overridden ndarray methods

    # TODO: integrate bottleneck

    def count(self, level=None):
        """
        Return number of non-NA/null observations in the Series

        Parameters
        ----------
        level : int or level name, default None
            If the axis is a MultiIndex (hierarchical), count along a
            particular level, collapsing into a smaller Series

        Returns
        -------
        nobs : int or Series (if level specified)
        """
        if level is None:
            return notna(com._values_from_object(self)).sum()

        if isinstance(level, compat.string_types):
            level = self.index._get_level_number(level)

        lev = self.index.levels[level]
        lab = np.array(self.index.labels[level], subok=False, copy=True)

        mask = lab == -1
        if mask.any():
            lab[mask] = cnt = len(lev)
            lev = lev.insert(cnt, lev._na_value)

        obs = lab[notna(self.values)]
        out = np.bincount(obs, minlength=len(lev) or None)
        return self._constructor(out, index=lev,
                                 dtype='int64').__finalize__(self)

    def mode(self):
        """Return the mode(s) of the dataset.

        Always returns Series even if only one value is returned.

        Returns
        -------
        modes : Series (sorted)
        """
        # TODO: Add option for bins like value_counts()
        return algorithms.mode(self)

    def unique(self):
        """
        Return unique values of Series object.

        Uniques are returned in order of appearance. Hash table-based unique,
        therefore does NOT sort.

        Returns
        -------
        ndarray or Categorical
            The unique values returned as a NumPy array. In case of categorical
            data type, returned as a Categorical.

        See Also
        --------
        pandas.unique : top-level unique method for any 1-d array-like object.
        Index.unique : return Index with unique values from an Index object.

        Examples
        --------
        >>> pd.Series([2, 1, 3, 3], name='A').unique()
        array([2, 1, 3])

        >>> pd.Series([pd.Timestamp('2016-01-01') for _ in range(3)]).unique()
        array(['2016-01-01T00:00:00.000000000'], dtype='datetime64[ns]')

        >>> pd.Series([pd.Timestamp('2016-01-01', tz='US/Eastern')
        ...            for _ in range(3)]).unique()
        array([Timestamp('2016-01-01 00:00:00-0500', tz='US/Eastern')],
              dtype=object)

        An unordered Categorical will return categories in the order of
        appearance.

        >>> pd.Series(pd.Categorical(list('baabc'))).unique()
        [b, a, c]
        Categories (3, object): [b, a, c]

        An ordered Categorical preserves the category ordering.

        >>> pd.Series(pd.Categorical(list('baabc'), categories=list('abc'),
        ...                          ordered=True)).unique()
        [b, a, c]
        Categories (3, object): [a < b < c]
        """
        result = super(Series, self).unique()

        if is_datetime64tz_dtype(self.dtype):
            # we are special casing datetime64tz_dtype
            # to return an object array of tz-aware Timestamps

            # TODO: it must return DatetimeArray with tz in pandas 2.0
            result = result.astype(object).values

        return result

    def drop_duplicates(self, keep='first', inplace=False):
        """
        Return Series with duplicate values removed.

        Parameters
        ----------
        keep : {'first', 'last', ``False``}, default 'first'
            - 'first' : Drop duplicates except for the first occurrence.
            - 'last' : Drop duplicates except for the last occurrence.
            - ``False`` : Drop all duplicates.
        inplace : boolean, default ``False``
            If ``True``, performs operation inplace and returns None.

        Returns
        -------
        deduplicated : Series

        See Also
        --------
        Index.drop_duplicates : equivalent method on Index
        DataFrame.drop_duplicates : equivalent method on DataFrame
        Series.duplicated : related method on Series, indicating duplicate
            Series values.

        Examples
        --------
        Generate an Series with duplicated entries.

        >>> s = pd.Series(['lama', 'cow', 'lama', 'beetle', 'lama', 'hippo'],
        ...               name='animal')
        >>> s
        0      lama
        1       cow
        2      lama
        3    beetle
        4      lama
        5     hippo
        Name: animal, dtype: object

        With the 'keep' parameter, the selection behaviour of duplicated values
        can be changed. The value 'first' keeps the first occurrence for each
        set of duplicated entries. The default value of keep is 'first'.

        >>> s.drop_duplicates()
        0      lama
        1       cow
        3    beetle
        5     hippo
        Name: animal, dtype: object

        The value 'last' for parameter 'keep' keeps the last occurrence for
        each set of duplicated entries.

        >>> s.drop_duplicates(keep='last')
        1       cow
        3    beetle
        4      lama
        5     hippo
        Name: animal, dtype: object

        The value ``False`` for parameter 'keep' discards all sets of
        duplicated entries. Setting the value of 'inplace' to ``True`` performs
        the operation inplace and returns ``None``.

        >>> s.drop_duplicates(keep=False, inplace=True)
        >>> s
        1       cow
        3    beetle
        5     hippo
        Name: animal, dtype: object
        """
        return super(Series, self).drop_duplicates(keep=keep, inplace=inplace)

    def duplicated(self, keep='first'):
        """
        Indicate duplicate Series values.

        Duplicated values are indicated as ``True`` values in the resulting
        Series. Either all duplicates, all except the first or all except the
        last occurrence of duplicates can be indicated.

        Parameters
        ----------
        keep : {'first', 'last', False}, default 'first'
            - 'first' : Mark duplicates as ``True`` except for the first
              occurrence.
            - 'last' : Mark duplicates as ``True`` except for the last
              occurrence.
            - ``False`` : Mark all duplicates as ``True``.

        Examples
        --------
        By default, for each set of duplicated values, the first occurrence is
        set on False and all others on True:

        >>> animals = pd.Series(['lama', 'cow', 'lama', 'beetle', 'lama'])
        >>> animals.duplicated()
        0    False
        1    False
        2     True
        3    False
        4     True
        dtype: bool

        which is equivalent to

        >>> animals.duplicated(keep='first')
        0    False
        1    False
        2     True
        3    False
        4     True
        dtype: bool

        By using 'last', the last occurrence of each set of duplicated values
        is set on False and all others on True:

        >>> animals.duplicated(keep='last')
        0     True
        1    False
        2     True
        3    False
        4    False
        dtype: bool

        By setting keep on ``False``, all duplicates are True:

        >>> animals.duplicated(keep=False)
        0     True
        1    False
        2     True
        3    False
        4     True
        dtype: bool

        Returns
        -------
        pandas.core.series.Series

        See Also
        --------
        pandas.Index.duplicated : Equivalent method on pandas.Index
        pandas.DataFrame.duplicated : Equivalent method on pandas.DataFrame
        pandas.Series.drop_duplicates : Remove duplicate values from Series
        """
        return super(Series, self).duplicated(keep=keep)

    def idxmin(self, axis=None, skipna=True, *args, **kwargs):
        """
        Return the row label of the minimum value.

        If multiple values equal the minimum, the first row label with that
        value is returned.

        Parameters
        ----------
        skipna : boolean, default True
            Exclude NA/null values. If the entire Series is NA, the result
            will be NA.
        axis : int, default 0
            For compatibility with DataFrame.idxmin. Redundant for application
            on Series.
        *args, **kwargs
            Additional keywors have no effect but might be accepted
            for compatibility with NumPy.

        Returns
        -------
        idxmin : Index of minimum of values.

        Raises
        ------
        ValueError
            If the Series is empty.

        Notes
        -----
        This method is the Series version of ``ndarray.argmin``. This method
        returns the label of the minimum, while ``ndarray.argmin`` returns
        the position. To get the position, use ``series.values.argmin()``.

        See Also
        --------
        numpy.argmin : Return indices of the minimum values
            along the given axis.
        DataFrame.idxmin : Return index of first occurrence of minimum
            over requested axis.
        Series.idxmax : Return index *label* of the first occurrence
            of maximum of values.

        Examples
        --------
        >>> s = pd.Series(data=[1, None, 4, 1],
        ...               index=['A' ,'B' ,'C' ,'D'])
        >>> s
        A    1.0
        B    NaN
        C    4.0
        D    1.0
        dtype: float64

        >>> s.idxmin()
        'A'

        If `skipna` is False and there is an NA value in the data,
        the function returns ``nan``.

        >>> s.idxmin(skipna=False)
        nan
        """
        skipna = nv.validate_argmin_with_skipna(skipna, args, kwargs)
        i = nanops.nanargmin(com._values_from_object(self), skipna=skipna)
        if i == -1:
            return np.nan
        return self.index[i]

    def idxmax(self, axis=0, skipna=True, *args, **kwargs):
        """
        Return the row label of the maximum value.

        If multiple values equal the maximum, the first row label with that
        value is returned.

        Parameters
        ----------
        skipna : boolean, default True
            Exclude NA/null values. If the entire Series is NA, the result
            will be NA.
        axis : int, default 0
            For compatibility with DataFrame.idxmax. Redundant for application
            on Series.
        *args, **kwargs
            Additional keywors have no effect but might be accepted
            for compatibility with NumPy.

        Returns
        -------
        idxmax : Index of maximum of values.

        Raises
        ------
        ValueError
            If the Series is empty.

        Notes
        -----
        This method is the Series version of ``ndarray.argmax``. This method
        returns the label of the maximum, while ``ndarray.argmax`` returns
        the position. To get the position, use ``series.values.argmax()``.

        See Also
        --------
        numpy.argmax : Return indices of the maximum values
            along the given axis.
        DataFrame.idxmax : Return index of first occurrence of maximum
            over requested axis.
        Series.idxmin : Return index *label* of the first occurrence
            of minimum of values.

        Examples
        --------
        >>> s = pd.Series(data=[1, None, 4, 3, 4],
        ...               index=['A', 'B', 'C', 'D', 'E'])
        >>> s
        A    1.0
        B    NaN
        C    4.0
        D    3.0
        E    4.0
        dtype: float64

        >>> s.idxmax()
        'C'

        If `skipna` is False and there is an NA value in the data,
        the function returns ``nan``.

        >>> s.idxmax(skipna=False)
        nan
        """
        skipna = nv.validate_argmax_with_skipna(skipna, args, kwargs)
        i = nanops.nanargmax(com._values_from_object(self), skipna=skipna)
        if i == -1:
            return np.nan
        return self.index[i]

    # ndarray compat
    argmin = deprecate(
        'argmin', idxmin, '0.21.0',
        msg=dedent("""\
        'argmin' is deprecated, use 'idxmin' instead. The behavior of 'argmin'
        will be corrected to return the positional minimum in the future.
        Use 'series.values.argmin' to get the position of the minimum now.""")
    )
    argmax = deprecate(
        'argmax', idxmax, '0.21.0',
        msg=dedent("""\
        'argmax' is deprecated, use 'idxmax' instead. The behavior of 'argmax'
        will be corrected to return the positional maximum in the future.
        Use 'series.values.argmax' to get the position of the maximum now.""")
    )

    def round(self, decimals=0, *args, **kwargs):
        """
        Round each value in a Series to the given number of decimals.

        Parameters
        ----------
        decimals : int
            Number of decimal places to round to (default: 0).
            If decimals is negative, it specifies the number of
            positions to the left of the decimal point.

        Returns
        -------
        Series object

        See Also
        --------
        numpy.around
        DataFrame.round

        """
        nv.validate_round(args, kwargs)
        result = com._values_from_object(self).round(decimals)
        result = self._constructor(result, index=self.index).__finalize__(self)

        return result

    def quantile(self, q=0.5, interpolation='linear'):
        """
        Return value at the given quantile, a la numpy.percentile.

        Parameters
        ----------
        q : float or array-like, default 0.5 (50% quantile)
            0 <= q <= 1, the quantile(s) to compute
        interpolation : {'linear', 'lower', 'higher', 'midpoint', 'nearest'}
            .. versionadded:: 0.18.0

            This optional parameter specifies the interpolation method to use,
            when the desired quantile lies between two data points `i` and `j`:

                * linear: `i + (j - i) * fraction`, where `fraction` is the
                  fractional part of the index surrounded by `i` and `j`.
                * lower: `i`.
                * higher: `j`.
                * nearest: `i` or `j` whichever is nearest.
                * midpoint: (`i` + `j`) / 2.

        Returns
        -------
        quantile : float or Series
            if ``q`` is an array, a Series will be returned where the
            index is ``q`` and the values are the quantiles.

        Examples
        --------
        >>> s = Series([1, 2, 3, 4])
        >>> s.quantile(.5)
        2.5
        >>> s.quantile([.25, .5, .75])
        0.25    1.75
        0.50    2.50
        0.75    3.25
        dtype: float64

        See Also
        --------
        pandas.core.window.Rolling.quantile
        """

        self._check_percentile(q)

        result = self._data.quantile(qs=q, interpolation=interpolation)

        if is_list_like(q):
            return self._constructor(result,
                                     index=Float64Index(q),
                                     name=self.name)
        else:
            # scalar
            return result

    def corr(self, other, method='pearson', min_periods=None):
        """
        Compute correlation with `other` Series, excluding missing values

        Parameters
        ----------
        other : Series
        method : {'pearson', 'kendall', 'spearman'}
            * pearson : standard correlation coefficient
            * kendall : Kendall Tau correlation coefficient
            * spearman : Spearman rank correlation
        min_periods : int, optional
            Minimum number of observations needed to have a valid result


        Returns
        -------
        correlation : float
        """
        this, other = self.align(other, join='inner', copy=False)
        if len(this) == 0:
            return np.nan
        return nanops.nancorr(this.values, other.values, method=method,
                              min_periods=min_periods)

    def cov(self, other, min_periods=None):
        """
        Compute covariance with Series, excluding missing values

        Parameters
        ----------
        other : Series
        min_periods : int, optional
            Minimum number of observations needed to have a valid result

        Returns
        -------
        covariance : float

        Normalized by N-1 (unbiased estimator).
        """
        this, other = self.align(other, join='inner', copy=False)
        if len(this) == 0:
            return np.nan
        return nanops.nancov(this.values, other.values,
                             min_periods=min_periods)

    def diff(self, periods=1):
        """
        First discrete difference of element.

        Calculates the difference of a Series element compared with another
        element in the Series (default is element in previous row).

        Parameters
        ----------
        periods : int, default 1
            Periods to shift for calculating difference, accepts negative
            values.

        Returns
        -------
        diffed : Series

        See Also
        --------
        Series.pct_change: Percent change over given number of periods.
        Series.shift: Shift index by desired number of periods with an
            optional time freq.
        DataFrame.diff: First discrete difference of object

        Examples
        --------
        Difference with previous row

        >>> s = pd.Series([1, 1, 2, 3, 5, 8])
        >>> s.diff()
        0    NaN
        1    0.0
        2    1.0
        3    1.0
        4    2.0
        5    3.0
        dtype: float64

        Difference with 3rd previous row

        >>> s.diff(periods=3)
        0    NaN
        1    NaN
        2    NaN
        3    2.0
        4    4.0
        5    6.0
        dtype: float64

        Difference with following row

        >>> s.diff(periods=-1)
        0    0.0
        1   -1.0
        2   -1.0
        3   -2.0
        4   -3.0
        5    NaN
        dtype: float64
        """
        result = algorithms.diff(com._values_from_object(self), periods)
        return self._constructor(result, index=self.index).__finalize__(self)

    def autocorr(self, lag=1):
        """
        Lag-N autocorrelation

        Parameters
        ----------
        lag : int, default 1
            Number of lags to apply before performing autocorrelation.

        Returns
        -------
        autocorr : float
        """
        return self.corr(self.shift(lag))

    def dot(self, other):
        """
        Matrix multiplication with DataFrame or inner-product with Series
        objects. Can also be called using `self @ other` in Python >= 3.5.

        Parameters
        ----------
        other : Series or DataFrame

        Returns
        -------
        dot_product : scalar or Series
        """
        from pandas.core.frame import DataFrame
        if isinstance(other, (Series, DataFrame)):
            common = self.index.union(other.index)
            if (len(common) > len(self.index) or
                    len(common) > len(other.index)):
                raise ValueError('matrices are not aligned')

            left = self.reindex(index=common, copy=False)
            right = other.reindex(index=common, copy=False)
            lvals = left.values
            rvals = right.values
        else:
            left = self
            lvals = self.values
            rvals = np.asarray(other)
            if lvals.shape[0] != rvals.shape[0]:
                raise Exception('Dot product shape mismatch, %s vs %s' %
                                (lvals.shape, rvals.shape))

        if isinstance(other, DataFrame):
            return self._constructor(np.dot(lvals, rvals),
                                     index=other.columns).__finalize__(self)
        elif isinstance(other, Series):
            return np.dot(lvals, rvals)
        elif isinstance(rvals, np.ndarray):
            return np.dot(lvals, rvals)
        else:  # pragma: no cover
            raise TypeError('unsupported type: %s' % type(other))

    def __matmul__(self, other):
        """ Matrix multiplication using binary `@` operator in Python>=3.5 """
        return self.dot(other)

    def __rmatmul__(self, other):
        """ Matrix multiplication using binary `@` operator in Python>=3.5 """
        return self.dot(other)

    @Substitution(klass='Series')
    @Appender(base._shared_docs['searchsorted'])
    @deprecate_kwarg(old_arg_name='v', new_arg_name='value')
    def searchsorted(self, value, side='left', sorter=None):
        if sorter is not None:
            sorter = _ensure_platform_int(sorter)
        return self._values.searchsorted(Series(value)._values,
                                         side=side, sorter=sorter)

    # -------------------------------------------------------------------
    # Combination

    def append(self, to_append, ignore_index=False, verify_integrity=False):
        """
        Concatenate two or more Series.

        Parameters
        ----------
        to_append : Series or list/tuple of Series
        ignore_index : boolean, default False
            If True, do not use the index labels.

            .. versionadded:: 0.19.0

        verify_integrity : boolean, default False
            If True, raise Exception on creating index with duplicates

        Notes
        -----
        Iteratively appending to a Series can be more computationally intensive
        than a single concatenate. A better solution is to append values to a
        list and then concatenate the list with the original Series all at
        once.

        See also
        --------
        pandas.concat : General function to concatenate DataFrame, Series
            or Panel objects

        Returns
        -------
        appended : Series

        Examples
        --------
        >>> s1 = pd.Series([1, 2, 3])
        >>> s2 = pd.Series([4, 5, 6])
        >>> s3 = pd.Series([4, 5, 6], index=[3,4,5])
        >>> s1.append(s2)
        0    1
        1    2
        2    3
        0    4
        1    5
        2    6
        dtype: int64

        >>> s1.append(s3)
        0    1
        1    2
        2    3
        3    4
        4    5
        5    6
        dtype: int64

        With `ignore_index` set to True:

        >>> s1.append(s2, ignore_index=True)
        0    1
        1    2
        2    3
        3    4
        4    5
        5    6
        dtype: int64

        With `verify_integrity` set to True:

        >>> s1.append(s2, verify_integrity=True)
        Traceback (most recent call last):
        ...
        ValueError: Indexes have overlapping values: [0, 1, 2]


        """
        from pandas.core.reshape.concat import concat

        if isinstance(to_append, (list, tuple)):
            to_concat = [self] + to_append
        else:
            to_concat = [self, to_append]
        return concat(to_concat, ignore_index=ignore_index,
                      verify_integrity=verify_integrity)

    def _binop(self, other, func, level=None, fill_value=None):
        """
        Perform generic binary operation with optional fill value

        Parameters
        ----------
        other : Series
        func : binary operator
        fill_value : float or object
            Value to substitute for NA/null values. If both Series are NA in a
            location, the result will be NA regardless of the passed fill value
        level : int or level name, default None
            Broadcast across a level, matching Index values on the
            passed MultiIndex level

        Returns
        -------
        combined : Series
        """
        if not isinstance(other, Series):
            raise AssertionError('Other operand must be Series')

        new_index = self.index
        this = self

        if not self.index.equals(other.index):
            this, other = self.align(other, level=level, join='outer',
                                     copy=False)
            new_index = this.index

        this_vals, other_vals = ops.fill_binop(this.values, other.values,
                                               fill_value)

        with np.errstate(all='ignore'):
            result = func(this_vals, other_vals)
        name = ops.get_op_result_name(self, other)
        result = self._constructor(result, index=new_index, name=name)
        result = result.__finalize__(self)
        if name is None:
            # When name is None, __finalize__ overwrites current name
            result.name = None
        return result

    def combine(self, other, func, fill_value=np.nan):
        """
        Perform elementwise binary operation on two Series using given function
        with optional fill value when an index is missing from one Series or
        the other

        Parameters
        ----------
        other : Series or scalar value
        func : function
            Function that takes two scalars as inputs and return a scalar
        fill_value : scalar value

        Returns
        -------
        result : Series

        Examples
        --------
        >>> s1 = Series([1, 2])
        >>> s2 = Series([0, 3])
        >>> s1.combine(s2, lambda x1, x2: x1 if x1 < x2 else x2)
        0    0
        1    2
        dtype: int64

        See Also
        --------
        Series.combine_first : Combine Series values, choosing the calling
            Series's values first
        """
        if isinstance(other, Series):
            new_index = self.index.union(other.index)
            new_name = ops.get_op_result_name(self, other)
            new_values = np.empty(len(new_index), dtype=self.dtype)
            for i, idx in enumerate(new_index):
                lv = self.get(idx, fill_value)
                rv = other.get(idx, fill_value)
                with np.errstate(all='ignore'):
                    new_values[i] = func(lv, rv)
        else:
            new_index = self.index
            with np.errstate(all='ignore'):
                new_values = func(self._values, other)
            new_name = self.name
        return self._constructor(new_values, index=new_index, name=new_name)

    def combine_first(self, other):
        """
        Combine Series values, choosing the calling Series's values
        first. Result index will be the union of the two indexes

        Parameters
        ----------
        other : Series

        Returns
        -------
        combined : Series

        Examples
        --------
        >>> s1 = pd.Series([1, np.nan])
        >>> s2 = pd.Series([3, 4])
        >>> s1.combine_first(s2)
        0    1.0
        1    4.0
        dtype: float64

        See Also
        --------
        Series.combine : Perform elementwise operation on two Series
            using a given function
        """
        new_index = self.index.union(other.index)
        this = self.reindex(new_index, copy=False)
        other = other.reindex(new_index, copy=False)
        # TODO: do we need name?
        name = ops.get_op_result_name(self, other)  # noqa
        rs_vals = com._where_compat(isna(this), other._values, this._values)
        return self._constructor(rs_vals, index=new_index).__finalize__(self)

    def update(self, other):
        """
        Modify Series in place using non-NA values from passed
        Series. Aligns on index

        Parameters
        ----------
        other : Series

        Examples
        --------
        >>> s = pd.Series([1, 2, 3])
        >>> s.update(pd.Series([4, 5, 6]))
        >>> s
        0    4
        1    5
        2    6
        dtype: int64

        >>> s = pd.Series(['a', 'b', 'c'])
        >>> s.update(pd.Series(['d', 'e'], index=[0, 2]))
        >>> s
        0    d
        1    b
        2    e
        dtype: object

        >>> s = pd.Series([1, 2, 3])
        >>> s.update(pd.Series([4, 5, 6, 7, 8]))
        >>> s
        0    4
        1    5
        2    6
        dtype: int64

        If ``other`` contains NaNs the corresponding values are not updated
        in the original Series.

        >>> s = pd.Series([1, 2, 3])
        >>> s.update(pd.Series([4, np.nan, 6]))
        >>> s
        0    4
        1    2
        2    6
        dtype: int64

        """
        other = other.reindex_like(self)
        mask = notna(other)

        self._data = self._data.putmask(mask=mask, new=other, inplace=True)
        self._maybe_update_cacher()

    # ----------------------------------------------------------------------
    # Reindexing, sorting

    def sort_values(self, axis=0, ascending=True, inplace=False,
                    kind='quicksort', na_position='last'):
        """
        Sort by the values.

        Sort a Series in ascending or descending order by some
        criterion.

        Parameters
        ----------
        axis : {0 or 'index'}, default 0
            Axis to direct sorting. The value 'index' is accepted for
            compatibility with DataFrame.sort_values.
        ascending : bool, default True
            If True, sort values in ascending order, otherwise descending.
        inplace : bool, default False
            If True, perform operation in-place.
        kind : {'quicksort', 'mergesort' or 'heapsort'}, default 'quicksort'
            Choice of sorting algorithm. See also :func:`numpy.sort` for more
            information. 'mergesort' is the only stable  algorithm.
        na_position : {'first' or 'last'}, default 'last'
            Argument 'first' puts NaNs at the beginning, 'last' puts NaNs at
            the end.

        Returns
        -------
        Series
            Series ordered by values.

        See Also
        --------
        Series.sort_index : Sort by the Series indices.
        DataFrame.sort_values : Sort DataFrame by the values along either axis.
        DataFrame.sort_index : Sort DataFrame by indices.

        Examples
        --------
        >>> s = pd.Series([np.nan, 1, 3, 10, 5])
        >>> s
        0     NaN
        1     1.0
        2     3.0
        3     10.0
        4     5.0
        dtype: float64

        Sort values ascending order (default behaviour)

        >>> s.sort_values(ascending=True)
        1     1.0
        2     3.0
        4     5.0
        3    10.0
        0     NaN
        dtype: float64

        Sort values descending order

        >>> s.sort_values(ascending=False)
        3    10.0
        4     5.0
        2     3.0
        1     1.0
        0     NaN
        dtype: float64

        Sort values inplace

        >>> s.sort_values(ascending=False, inplace=True)
        >>> s
        3    10.0
        4     5.0
        2     3.0
        1     1.0
        0     NaN
        dtype: float64

        Sort values putting NAs first

        >>> s.sort_values(na_position='first')
        0     NaN
        1     1.0
        2     3.0
        4     5.0
        3    10.0
        dtype: float64

        Sort a series of strings

        >>> s = pd.Series(['z', 'b', 'd', 'a', 'c'])
        >>> s
        0    z
        1    b
        2    d
        3    a
        4    c
        dtype: object

        >>> s.sort_values()
        3    a
        1    b
        4    c
        2    d
        0    z
        dtype: object
        """
        inplace = validate_bool_kwarg(inplace, 'inplace')
        axis = self._get_axis_number(axis)

        # GH 5856/5853
        if inplace and self._is_cached:
            raise ValueError("This Series is a view of some other array, to "
                             "sort in-place you must create a copy")

        def _try_kind_sort(arr):
            # easier to ask forgiveness than permission
            try:
                # if kind==mergesort, it can fail for object dtype
                return arr.argsort(kind=kind)
            except TypeError:
                # stable sort not available for object dtype
                # uses the argsort default quicksort
                return arr.argsort(kind='quicksort')

        arr = self._values
        sortedIdx = np.empty(len(self), dtype=np.int32)

        bad = isna(arr)

        good = ~bad
        idx = com._default_index(len(self))

        argsorted = _try_kind_sort(arr[good])

        if is_list_like(ascending):
            if len(ascending) != 1:
                raise ValueError('Length of ascending (%d) must be 1 '
                                 'for Series' % (len(ascending)))
            ascending = ascending[0]

        if not is_bool(ascending):
            raise ValueError('ascending must be boolean')

        if not ascending:
            argsorted = argsorted[::-1]

        if na_position == 'last':
            n = good.sum()
            sortedIdx[:n] = idx[good][argsorted]
            sortedIdx[n:] = idx[bad]
        elif na_position == 'first':
            n = bad.sum()
            sortedIdx[n:] = idx[good][argsorted]
            sortedIdx[:n] = idx[bad]
        else:
            raise ValueError('invalid na_position: {!r}'.format(na_position))

        result = self._constructor(arr[sortedIdx], index=self.index[sortedIdx])

        if inplace:
            self._update_inplace(result)
        else:
            return result.__finalize__(self)

    def sort_index(self, axis=0, level=None, ascending=True, inplace=False,
                   kind='quicksort', na_position='last', sort_remaining=True):
        """
        Sort Series by index labels.

        Returns a new Series sorted by label if `inplace` argument is
        ``False``, otherwise updates the original series and returns None.

        Parameters
        ----------
        axis : int, default 0
            Axis to direct sorting. This can only be 0 for Series.
        level : int, optional
            If not None, sort on values in specified index level(s).
        ascending : bool, default true
            Sort ascending vs. descending.
        inplace : bool, default False
            If True, perform operation in-place.
        kind : {'quicksort', 'mergesort', 'heapsort'}, default 'quicksort'
            Choice of sorting algorithm. See also :func:`numpy.sort` for more
            information.  'mergesort' is the only stable algorithm. For
            DataFrames, this option is only applied when sorting on a single
            column or label.
        na_position : {'first', 'last'}, default 'last'
            If 'first' puts NaNs at the beginning, 'last' puts NaNs at the end.
            Not implemented for MultiIndex.
        sort_remaining : bool, default True
            If true and sorting by level and index is multilevel, sort by other
            levels too (in order) after sorting by specified level.

        Returns
        -------
        pandas.Series
            The original Series sorted by the labels

        See Also
        --------
        DataFrame.sort_index: Sort DataFrame by the index
        DataFrame.sort_values: Sort DataFrame by the value
        Series.sort_values : Sort Series by the value

        Examples
        --------
        >>> s = pd.Series(['a', 'b', 'c', 'd'], index=[3, 2, 1, 4])
        >>> s.sort_index()
        1    c
        2    b
        3    a
        4    d
        dtype: object

        Sort Descending

        >>> s.sort_index(ascending=False)
        4    d
        3    a
        2    b
        1    c
        dtype: object

        Sort Inplace

        >>> s.sort_index(inplace=True)
        >>> s
        1    c
        2    b
        3    a
        4    d
        dtype: object

        By default NaNs are put at the end, but use `na_position` to place
        them at the beginning

        >>> s = pd.Series(['a', 'b', 'c', 'd'], index=[3, 2, 1, np.nan])
        >>> s.sort_index(na_position='first')
        NaN     d
         1.0    c
         2.0    b
         3.0    a
        dtype: object

        Specify index level to sort

        >>> arrays = [np.array(['qux', 'qux', 'foo', 'foo',
        ...                     'baz', 'baz', 'bar', 'bar']),
        ...           np.array(['two', 'one', 'two', 'one',
        ...                     'two', 'one', 'two', 'one'])]
        >>> s = pd.Series([1, 2, 3, 4, 5, 6, 7, 8], index=arrays)
        >>> s.sort_index(level=1)
        bar  one    8
        baz  one    6
        foo  one    4
        qux  one    2
        bar  two    7
        baz  two    5
        foo  two    3
        qux  two    1
        dtype: int64

        Does not sort by remaining levels when sorting by levels

        >>> s.sort_index(level=1, sort_remaining=False)
        qux  one    2
        foo  one    4
        baz  one    6
        bar  one    8
        qux  two    1
        foo  two    3
        baz  two    5
        bar  two    7
        dtype: int64
        """
        # TODO: this can be combined with DataFrame.sort_index impl as
        # almost identical
        inplace = validate_bool_kwarg(inplace, 'inplace')
        axis = self._get_axis_number(axis)
        index = self.index

        if level:
            new_index, indexer = index.sortlevel(level, ascending=ascending,
                                                 sort_remaining=sort_remaining)
        elif isinstance(index, MultiIndex):
            from pandas.core.sorting import lexsort_indexer
            labels = index._sort_levels_monotonic()
            indexer = lexsort_indexer(labels._get_labels_for_sorting(),
                                      orders=ascending,
                                      na_position=na_position)
        else:
            from pandas.core.sorting import nargsort

            # Check monotonic-ness before sort an index
            # GH11080
            if ((ascending and index.is_monotonic_increasing) or
                    (not ascending and index.is_monotonic_decreasing)):
                if inplace:
                    return
                else:
                    return self.copy()

            indexer = nargsort(index, kind=kind, ascending=ascending,
                               na_position=na_position)

        indexer = _ensure_platform_int(indexer)
        new_index = index.take(indexer)
        new_index = new_index._sort_levels_monotonic()

        new_values = self._values.take(indexer)
        result = self._constructor(new_values, index=new_index)

        if inplace:
            self._update_inplace(result)
        else:
            return result.__finalize__(self)

    def argsort(self, axis=0, kind='quicksort', order=None):
        """
        Overrides ndarray.argsort. Argsorts the value, omitting NA/null values,
        and places the result in the same locations as the non-NA values

        Parameters
        ----------
        axis : int (can only be zero)
        kind : {'mergesort', 'quicksort', 'heapsort'}, default 'quicksort'
            Choice of sorting algorithm. See np.sort for more
            information. 'mergesort' is the only stable algorithm
        order : ignored

        Returns
        -------
        argsorted : Series, with -1 indicated where nan values are present

        See also
        --------
        numpy.ndarray.argsort
        """
        values = self._values
        mask = isna(values)

        if mask.any():
            result = Series(-1, index=self.index, name=self.name,
                            dtype='int64')
            notmask = ~mask
            result[notmask] = np.argsort(values[notmask], kind=kind)
            return self._constructor(result,
                                     index=self.index).__finalize__(self)
        else:
            return self._constructor(
                np.argsort(values, kind=kind), index=self.index,
                dtype='int64').__finalize__(self)

    def nlargest(self, n=5, keep='first'):
        """
        Return the largest `n` elements.

        Parameters
        ----------
        n : int
            Return this many descending sorted values
        keep : {'first', 'last'}, default 'first'
            Where there are duplicate values:
            - ``first`` : take the first occurrence.
            - ``last`` : take the last occurrence.

        Returns
        -------
        top_n : Series
            The n largest values in the Series, in sorted order

        Notes
        -----
        Faster than ``.sort_values(ascending=False).head(n)`` for small `n`
        relative to the size of the ``Series`` object.

        See Also
        --------
        Series.nsmallest

        Examples
        --------
        >>> import pandas as pd
        >>> import numpy as np
        >>> s = pd.Series(np.random.randn(10**6))
        >>> s.nlargest(10)  # only sorts up to the N requested
        219921    4.644710
        82124     4.608745
        421689    4.564644
        425277    4.447014
        718691    4.414137
        43154     4.403520
        283187    4.313922
        595519    4.273635
        503969    4.250236
        121637    4.240952
        dtype: float64
        """
        return algorithms.SelectNSeries(self, n=n, keep=keep).nlargest()

    def nsmallest(self, n=5, keep='first'):
        """
        Return the smallest `n` elements.

        Parameters
        ----------
        n : int
            Return this many ascending sorted values
        keep : {'first', 'last'}, default 'first'
            Where there are duplicate values:
            - ``first`` : take the first occurrence.
            - ``last`` : take the last occurrence.

        Returns
        -------
        bottom_n : Series
            The n smallest values in the Series, in sorted order

        Notes
        -----
        Faster than ``.sort_values().head(n)`` for small `n` relative to
        the size of the ``Series`` object.

        See Also
        --------
        Series.nlargest

        Examples
        --------
        >>> import pandas as pd
        >>> import numpy as np
        >>> s = pd.Series(np.random.randn(10**6))
        >>> s.nsmallest(10)  # only sorts up to the N requested
        288532   -4.954580
        732345   -4.835960
        64803    -4.812550
        446457   -4.609998
        501225   -4.483945
        669476   -4.472935
        973615   -4.401699
        621279   -4.355126
        773916   -4.347355
        359919   -4.331927
        dtype: float64
        """
        return algorithms.SelectNSeries(self, n=n, keep=keep).nsmallest()

    def sortlevel(self, level=0, ascending=True, sort_remaining=True):
        """Sort Series with MultiIndex by chosen level. Data will be
        lexicographically sorted by the chosen level followed by the other
        levels (in order),

        .. deprecated:: 0.20.0
            Use :meth:`Series.sort_index`

        Parameters
        ----------
        level : int or level name, default None
        ascending : bool, default True

        Returns
        -------
        sorted : Series

        See Also
        --------
        Series.sort_index(level=...)

        """
        warnings.warn("sortlevel is deprecated, use sort_index(level=...)",
                      FutureWarning, stacklevel=2)
        return self.sort_index(level=level, ascending=ascending,
                               sort_remaining=sort_remaining)

    def swaplevel(self, i=-2, j=-1, copy=True):
        """
        Swap levels i and j in a MultiIndex

        Parameters
        ----------
        i, j : int, string (can be mixed)
            Level of index to be swapped. Can pass level name as string.

        Returns
        -------
        swapped : Series

        .. versionchanged:: 0.18.1

           The indexes ``i`` and ``j`` are now optional, and default to
           the two innermost levels of the index.

        """
        new_index = self.index.swaplevel(i, j)
        return self._constructor(self._values, index=new_index,
                                 copy=copy).__finalize__(self)

    def reorder_levels(self, order):
        """
        Rearrange index levels using input order. May not drop or duplicate
        levels

        Parameters
        ----------
        order : list of int representing new level order.
               (reference level by number or key)
        axis : where to reorder levels

        Returns
        -------
        type of caller (new object)
        """
        if not isinstance(self.index, MultiIndex):  # pragma: no cover
            raise Exception('Can only reorder levels on a hierarchical axis.')

        result = self.copy()
        result.index = result.index.reorder_levels(order)
        return result

    def unstack(self, level=-1, fill_value=None):
        """
        Unstack, a.k.a. pivot, Series with MultiIndex to produce DataFrame.
        The level involved will automatically get sorted.

        Parameters
        ----------
        level : int, string, or list of these, default last level
            Level(s) to unstack, can pass level name
        fill_value : replace NaN with this value if the unstack produces
            missing values

            .. versionadded:: 0.18.0

        Examples
        --------
        >>> s = pd.Series([1, 2, 3, 4],
        ...     index=pd.MultiIndex.from_product([['one', 'two'], ['a', 'b']]))
        >>> s
        one  a    1
             b    2
        two  a    3
             b    4
        dtype: int64

        >>> s.unstack(level=-1)
             a  b
        one  1  2
        two  3  4

        >>> s.unstack(level=0)
           one  two
        a    1    3
        b    2    4

        Returns
        -------
        unstacked : DataFrame
        """
        from pandas.core.reshape.reshape import unstack
        return unstack(self, level, fill_value)

    # ----------------------------------------------------------------------
    # function application

    def map(self, arg, na_action=None):
        """
        Map values of Series using input correspondence (a dict, Series, or
        function).

        Parameters
        ----------
        arg : function, dict, or Series
            Mapping correspondence.
        na_action : {None, 'ignore'}
            If 'ignore', propagate NA values, without passing them to the
            mapping correspondence.

        Returns
        -------
        y : Series
            Same index as caller.

        Examples
        --------

        Map inputs to outputs (both of type `Series`):

        >>> x = pd.Series([1,2,3], index=['one', 'two', 'three'])
        >>> x
        one      1
        two      2
        three    3
        dtype: int64

        >>> y = pd.Series(['foo', 'bar', 'baz'], index=[1,2,3])
        >>> y
        1    foo
        2    bar
        3    baz

        >>> x.map(y)
        one   foo
        two   bar
        three baz

        If `arg` is a dictionary, return a new Series with values converted
        according to the dictionary's mapping:

        >>> z = {1: 'A', 2: 'B', 3: 'C'}

        >>> x.map(z)
        one   A
        two   B
        three C

        Use na_action to control whether NA values are affected by the mapping
        function.

        >>> s = pd.Series([1, 2, 3, np.nan])

        >>> s2 = s.map('this is a string {}'.format, na_action=None)
        0    this is a string 1.0
        1    this is a string 2.0
        2    this is a string 3.0
        3    this is a string nan
        dtype: object

        >>> s3 = s.map('this is a string {}'.format, na_action='ignore')
        0    this is a string 1.0
        1    this is a string 2.0
        2    this is a string 3.0
        3                     NaN
        dtype: object

        See Also
        --------
        Series.apply : For applying more complex functions on a Series.
        DataFrame.apply : Apply a function row-/column-wise.
        DataFrame.applymap : Apply a function elementwise on a whole DataFrame.

        Notes
        -----
        When `arg` is a dictionary, values in Series that are not in the
        dictionary (as keys) are converted to ``NaN``. However, if the
        dictionary is a ``dict`` subclass that defines ``__missing__`` (i.e.
        provides a method for default values), then this default is used
        rather than ``NaN``:

        >>> from collections import Counter
        >>> counter = Counter()
        >>> counter['bar'] += 1
        >>> y.map(counter)
        1    0
        2    1
        3    0
        dtype: int64
        """
        new_values = super(Series, self)._map_values(
            arg, na_action=na_action)
        return self._constructor(new_values,
                                 index=self.index).__finalize__(self)

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
        return self

    _agg_doc = dedent("""
    Examples
    --------

    >>> s = Series(np.random.randn(10))

    >>> s.agg('min')
    -1.3018049988556679

    >>> s.agg(['min', 'max'])
    min   -1.301805
    max    1.127688
    dtype: float64

    See also
    --------
    pandas.Series.apply
    pandas.Series.transform

    """)

    @Appender(_agg_doc)
    @Appender(generic._shared_docs['aggregate'] % dict(
        versionadded='.. versionadded:: 0.20.0',
        **_shared_doc_kwargs))
    def aggregate(self, func, axis=0, *args, **kwargs):
        axis = self._get_axis_number(axis)
        result, how = self._aggregate(func, *args, **kwargs)
        if result is None:

            # we can be called from an inner function which
            # passes this meta-data
            kwargs.pop('_axis', None)
            kwargs.pop('_level', None)

            # try a regular apply, this evaluates lambdas
            # row-by-row; however if the lambda is expected a Series
            # expression, e.g.: lambda x: x-x.quantile(0.25)
            # this will fail, so we can try a vectorized evaluation

            # we cannot FIRST try the vectorized evaluation, because
            # then .agg and .apply would have different semantics if the
            # operation is actually defined on the Series, e.g. str
            try:
                result = self.apply(func, *args, **kwargs)
            except (ValueError, AttributeError, TypeError):
                result = func(self, *args, **kwargs)

        return result

    agg = aggregate

    def apply(self, func, convert_dtype=True, args=(), **kwds):
        """
        Invoke function on values of Series. Can be ufunc (a NumPy function
        that applies to the entire Series) or a Python function that only works
        on single values

        Parameters
        ----------
        func : function
        convert_dtype : boolean, default True
            Try to find better dtype for elementwise function results. If
            False, leave as dtype=object
        args : tuple
            Positional arguments to pass to function in addition to the value
        Additional keyword arguments will be passed as keywords to the function

        Returns
        -------
        y : Series or DataFrame if func returns a Series

        See also
        --------
        Series.map: For element-wise operations
        Series.agg: only perform aggregating type operations
        Series.transform: only perform transformating type operations

        Examples
        --------

        Create a series with typical summer temperatures for each city.

        >>> import pandas as pd
        >>> import numpy as np
        >>> series = pd.Series([20, 21, 12], index=['London',
        ... 'New York','Helsinki'])
        >>> series
        London      20
        New York    21
        Helsinki    12
        dtype: int64

        Square the values by defining a function and passing it as an
        argument to ``apply()``.

        >>> def square(x):
        ...     return x**2
        >>> series.apply(square)
        London      400
        New York    441
        Helsinki    144
        dtype: int64

        Square the values by passing an anonymous function as an
        argument to ``apply()``.

        >>> series.apply(lambda x: x**2)
        London      400
        New York    441
        Helsinki    144
        dtype: int64

        Define a custom function that needs additional positional
        arguments and pass these additional arguments using the
        ``args`` keyword.

        >>> def subtract_custom_value(x, custom_value):
        ...     return x-custom_value

        >>> series.apply(subtract_custom_value, args=(5,))
        London      15
        New York    16
        Helsinki     7
        dtype: int64

        Define a custom function that takes keyword arguments
        and pass these arguments to ``apply``.

        >>> def add_custom_values(x, **kwargs):
        ...     for month in kwargs:
        ...         x+=kwargs[month]
        ...     return x

        >>> series.apply(add_custom_values, june=30, july=20, august=25)
        London      95
        New York    96
        Helsinki    87
        dtype: int64

        Use a function from the Numpy library.

        >>> series.apply(np.log)
        London      2.995732
        New York    3.044522
        Helsinki    2.484907
        dtype: float64


        """
        if len(self) == 0:
            return self._constructor(dtype=self.dtype,
                                     index=self.index).__finalize__(self)

        # dispatch to agg
        if isinstance(func, (list, dict)):
            return self.aggregate(func, *args, **kwds)

        # if we are a string, try to dispatch
        if isinstance(func, compat.string_types):
            return self._try_aggregate_string_function(func, *args, **kwds)

        # handle ufuncs and lambdas
        if kwds or args and not isinstance(func, np.ufunc):
            f = lambda x: func(x, *args, **kwds)
        else:
            f = func

        with np.errstate(all='ignore'):
            if isinstance(f, np.ufunc):
                return f(self)

            # row-wise access
            if is_extension_type(self.dtype):
                mapped = self._values.map(f)
            else:
                values = self.astype(object).values
                mapped = lib.map_infer(values, f, convert=convert_dtype)

        if len(mapped) and isinstance(mapped[0], Series):
            from pandas.core.frame import DataFrame
            return DataFrame(mapped.tolist(), index=self.index)
        else:
            return self._constructor(mapped,
                                     index=self.index).__finalize__(self)

    def _reduce(self, op, name, axis=0, skipna=True, numeric_only=None,
                filter_type=None, **kwds):
        """
        perform a reduction operation

        if we have an ndarray as a value, then simply perform the operation,
        otherwise delegate to the object

        """
        delegate = self._values
        if isinstance(delegate, np.ndarray):
            # Validate that 'axis' is consistent with Series's single axis.
            self._get_axis_number(axis)
            if numeric_only:
                raise NotImplementedError('Series.{0} does not implement '
                                          'numeric_only.'.format(name))
            with np.errstate(all='ignore'):
                return op(delegate, skipna=skipna, **kwds)

        return delegate._reduce(op=op, name=name, axis=axis, skipna=skipna,
                                numeric_only=numeric_only,
                                filter_type=filter_type, **kwds)

    def _reindex_indexer(self, new_index, indexer, copy):
        if indexer is None:
            if copy:
                return self.copy()
            return self

        new_values = algorithms.take_1d(self._values, indexer,
                                        allow_fill=True, fill_value=None)
        return self._constructor(new_values, index=new_index)

    def _needs_reindex_multi(self, axes, method, level):
        """ check if we do need a multi reindex; this is for compat with
        higher dims
        """
        return False

    @Appender(generic._shared_docs['align'] % _shared_doc_kwargs)
    def align(self, other, join='outer', axis=None, level=None, copy=True,
              fill_value=None, method=None, limit=None, fill_axis=0,
              broadcast_axis=None):
        return super(Series, self).align(other, join=join, axis=axis,
                                         level=level, copy=copy,
                                         fill_value=fill_value, method=method,
                                         limit=limit, fill_axis=fill_axis,
                                         broadcast_axis=broadcast_axis)

    def rename(self, index=None, **kwargs):
        """Alter Series index labels or name

        Function / dict values must be unique (1-to-1). Labels not contained in
        a dict / Series will be left as-is. Extra labels listed don't throw an
        error.

        Alternatively, change ``Series.name`` with a scalar value.

        See the :ref:`user guide <basics.rename>` for more.

        Parameters
        ----------
        index : scalar, hashable sequence, dict-like or function, optional
            dict-like or functions are transformations to apply to
            the index.
            Scalar or hashable sequence-like will alter the ``Series.name``
            attribute.
        copy : boolean, default True
            Also copy underlying data
        inplace : boolean, default False
            Whether to return a new %(klass)s. If True then value of copy is
            ignored.
        level : int or level name, default None
            In case of a MultiIndex, only rename labels in the specified
            level.

        Returns
        -------
        renamed : Series (new object)

        See Also
        --------
        pandas.Series.rename_axis

        Examples
        --------

        >>> s = pd.Series([1, 2, 3])
        >>> s
        0    1
        1    2
        2    3
        dtype: int64
        >>> s.rename("my_name") # scalar, changes Series.name
        0    1
        1    2
        2    3
        Name: my_name, dtype: int64
        >>> s.rename(lambda x: x ** 2)  # function, changes labels
        0    1
        1    2
        4    3
        dtype: int64
        >>> s.rename({1: 3, 2: 5})  # mapping, changes labels
        0    1
        3    2
        5    3
        dtype: int64

        """
        kwargs['inplace'] = validate_bool_kwarg(kwargs.get('inplace', False),
                                                'inplace')

        non_mapping = is_scalar(index) or (is_list_like(index) and
                                           not is_dict_like(index))
        if non_mapping:
            return self._set_name(index, inplace=kwargs.get('inplace'))
        return super(Series, self).rename(index=index, **kwargs)

    @Appender(generic._shared_docs['reindex'] % _shared_doc_kwargs)
    def reindex(self, index=None, **kwargs):
        return super(Series, self).reindex(index=index, **kwargs)

    def drop(self, labels=None, axis=0, index=None, columns=None,
             level=None, inplace=False, errors='raise'):
        """
        Return Series with specified index labels removed.

        Remove elements of a Series based on specifying the index labels.
        When using a multi-index, labels on different levels can be removed
        by specifying the level.

        Parameters
        ----------
        labels : single label or list-like
            Index labels to drop.
        axis : 0, default 0
            Redundant for application on Series.
        index, columns : None
            Redundant for application on Series, but index can be used instead
            of labels.

            .. versionadded:: 0.21.0
        level : int or level name, optional
            For MultiIndex, level for which the labels will be removed.
        inplace : bool, default False
            If True, do operation inplace and return None.
        errors : {'ignore', 'raise'}, default 'raise'
            If 'ignore', suppress error and only existing labels are dropped.

        Returns
        -------
        dropped : pandas.Series

        See Also
        --------
        Series.reindex : Return only specified index labels of Series.
        Series.dropna : Return series without null values.
        Series.drop_duplicates : Return Series with duplicate values removed.
        DataFrame.drop : Drop specified labels from rows or columns.

        Raises
        ------
        KeyError
            If none of the labels are found in the index.

        Examples
        --------
        >>> s = pd.Series(data=np.arange(3), index=['A','B','C'])
        >>> s
        A  0
        B  1
        C  2
        dtype: int64

        Drop labels B en C

        >>> s.drop(labels=['B','C'])
        A  0
        dtype: int64

        Drop 2nd level label in MultiIndex Series

        >>> midx = pd.MultiIndex(levels=[['lama', 'cow', 'falcon'],
        ...                              ['speed', 'weight', 'length']],
        ...                      labels=[[0, 0, 0, 1, 1, 1, 2, 2, 2],
        ...                              [0, 1, 2, 0, 1, 2, 0, 1, 2]])
        >>> s = pd.Series([45, 200, 1.2, 30, 250, 1.5, 320, 1, 0.3],
        ...               index=midx)
        >>> s
        lama    speed      45.0
                weight    200.0
                length      1.2
        cow     speed      30.0
                weight    250.0
                length      1.5
        falcon  speed     320.0
                weight      1.0
                length      0.3
        dtype: float64

        >>> s.drop(labels='weight', level=1)
        lama    speed      45.0
                length      1.2
        cow     speed      30.0
                length      1.5
        falcon  speed     320.0
                length      0.3
        dtype: float64
        """
        return super(Series, self).drop(labels=labels, axis=axis, index=index,
                                        columns=columns, level=level,
                                        inplace=inplace, errors=errors)

    @Substitution(**_shared_doc_kwargs)
    @Appender(generic.NDFrame.fillna.__doc__)
    def fillna(self, value=None, method=None, axis=None, inplace=False,
               limit=None, downcast=None, **kwargs):
        return super(Series, self).fillna(value=value, method=method,
                                          axis=axis, inplace=inplace,
                                          limit=limit, downcast=downcast,
                                          **kwargs)

    @Appender(generic._shared_docs['replace'] % _shared_doc_kwargs)
    def replace(self, to_replace=None, value=None, inplace=False, limit=None,
                regex=False, method='pad'):
        return super(Series, self).replace(to_replace=to_replace, value=value,
                                           inplace=inplace, limit=limit,
                                           regex=regex, method=method)

    @Appender(generic._shared_docs['shift'] % _shared_doc_kwargs)
    def shift(self, periods=1, freq=None, axis=0):
        return super(Series, self).shift(periods=periods, freq=freq, axis=axis)

    def reindex_axis(self, labels, axis=0, **kwargs):
        """Conform Series to new index with optional filling logic.

        .. deprecated:: 0.21.0
            Use ``Series.reindex`` instead.
        """
        # for compatibility with higher dims
        if axis != 0:
            raise ValueError("cannot reindex series on non-zero axis!")
        msg = ("'.reindex_axis' is deprecated and will be removed in a future "
               "version. Use '.reindex' instead.")
        warnings.warn(msg, FutureWarning, stacklevel=2)

        return self.reindex(index=labels, **kwargs)

    def memory_usage(self, index=True, deep=False):
        """
        Return the memory usage of the Series.

        The memory usage can optionally include the contribution of
        the index and of elements of `object` dtype.

        Parameters
        ----------
        index : bool, default True
            Specifies whether to include the memory usage of the Series index.
        deep : bool, default False
            If True, introspect the data deeply by interrogating
            `object` dtypes for system-level memory consumption, and include
            it in the returned value.

        Returns
        -------
        int
            Bytes of memory consumed.

        See Also
        --------
        numpy.ndarray.nbytes : Total bytes consumed by the elements of the
            array.
        DataFrame.memory_usage : Bytes consumed by a DataFrame.

        Examples
        --------

        >>> s = pd.Series(range(3))
        >>> s.memory_usage()
        104

        Not including the index gives the size of the rest of the data, which
        is necessarily smaller:

        >>> s.memory_usage(index=False)
        24

        The memory footprint of `object` values is ignored by default:

        >>> s = pd.Series(["a", "b"])
        >>> s.values
        array(['a', 'b'], dtype=object)
        >>> s.memory_usage()
        96
        >>> s.memory_usage(deep=True)
        212
        """
        v = super(Series, self).memory_usage(deep=deep)
        if index:
            v += self.index.memory_usage(deep=deep)
        return v

    @Appender(generic._shared_docs['_take'])
    def _take(self, indices, axis=0, is_copy=False):

        indices = _ensure_platform_int(indices)
        new_index = self.index.take(indices)

        if is_categorical_dtype(self):
            # https://github.com/pandas-dev/pandas/issues/20664
            # TODO: remove when the default Categorical.take behavior changes
            indices = maybe_convert_indices(indices, len(self._get_axis(axis)))
            kwargs = {'allow_fill': False}
        else:
            kwargs = {}
        new_values = self._values.take(indices, **kwargs)

        result = (self._constructor(new_values, index=new_index,
                                    fastpath=True).__finalize__(self))

        # Maybe set copy if we didn't actually change the index.
        if is_copy:
            if not result._get_axis(axis).equals(self._get_axis(axis)):
                result._set_is_copy(self)

        return result

    def isin(self, values):
        """
        Check whether `values` are contained in Series.

        Return a boolean Series showing whether each element in the Series
        matches an element in the passed sequence of `values` exactly.

        Parameters
        ----------
        values : set or list-like
            The sequence of values to test. Passing in a single string will
            raise a ``TypeError``. Instead, turn a single string into a
            list of one element.

            .. versionadded:: 0.18.1

              Support for values as a set.

        Returns
        -------
        isin : Series (bool dtype)

        Raises
        ------
        TypeError
          * If `values` is a string

        See Also
        --------
        pandas.DataFrame.isin : equivalent method on DataFrame

        Examples
        --------

        >>> s = pd.Series(['lama', 'cow', 'lama', 'beetle', 'lama',
        ...                'hippo'], name='animal')
        >>> s.isin(['cow', 'lama'])
        0     True
        1     True
        2     True
        3    False
        4     True
        5    False
        Name: animal, dtype: bool

        Passing a single string as ``s.isin('lama')`` will raise an error. Use
        a list of one element instead:

        >>> s.isin(['lama'])
        0     True
        1    False
        2     True
        3    False
        4     True
        5    False
        Name: animal, dtype: bool
        """
        result = algorithms.isin(self, values)
        return self._constructor(result, index=self.index).__finalize__(self)

    def between(self, left, right, inclusive=True):
        """
        Return boolean Series equivalent to left <= series <= right.

        This function returns a boolean vector containing `True` wherever the
        corresponding Series element is between the boundary values `left` and
        `right`. NA values are treated as `False`.

        Parameters
        ----------
        left : scalar
            Left boundary.
        right : scalar
            Right boundary.
        inclusive : bool, default True
            Include boundaries.

        Returns
        -------
        Series
            Each element will be a boolean.

        Notes
        -----
        This function is equivalent to ``(left <= ser) & (ser <= right)``

        See Also
        --------
        pandas.Series.gt : Greater than of series and other
        pandas.Series.lt : Less than of series and other

        Examples
        --------
        >>> s = pd.Series([2, 0, 4, 8, np.nan])

        Boundary values are included by default:

        >>> s.between(1, 4)
        0     True
        1    False
        2     True
        3    False
        4    False
        dtype: bool

        With `inclusive` set to ``False`` boundary values are excluded:

        >>> s.between(1, 4, inclusive=False)
        0     True
        1    False
        2    False
        3    False
        4    False
        dtype: bool

        `left` and `right` can be any scalar value:

        >>> s = pd.Series(['Alice', 'Bob', 'Carol', 'Eve'])
        >>> s.between('Anna', 'Daniel')
        0    False
        1     True
        2     True
        3    False
        dtype: bool
        """
        if inclusive:
            lmask = self >= left
            rmask = self <= right
        else:
            lmask = self > left
            rmask = self < right

        return lmask & rmask

    @classmethod
    def from_csv(cls, path, sep=',', parse_dates=True, header=None,
                 index_col=0, encoding=None, infer_datetime_format=False):
        """Read CSV file.

        .. deprecated:: 0.21.0
            Use :func:`pandas.read_csv` instead.

        It is preferable to use the more powerful :func:`pandas.read_csv`
        for most general purposes, but ``from_csv`` makes for an easy
        roundtrip to and from a file (the exact counterpart of
        ``to_csv``), especially with a time Series.

        This method only differs from :func:`pandas.read_csv` in some defaults:

        - `index_col` is ``0`` instead of ``None`` (take first column as index
          by default)
        - `header` is ``None`` instead of ``0`` (the first row is not used as
          the column names)
        - `parse_dates` is ``True`` instead of ``False`` (try parsing the index
          as datetime by default)

        With :func:`pandas.read_csv`, the option ``squeeze=True`` can be used
        to return a Series like ``from_csv``.

        Parameters
        ----------
        path : string file path or file handle / StringIO
        sep : string, default ','
            Field delimiter
        parse_dates : boolean, default True
            Parse dates. Different default from read_table
        header : int, default None
            Row to use as header (skip prior rows)
        index_col : int or sequence, default 0
            Column to use for index. If a sequence is given, a MultiIndex
            is used. Different default from read_table
        encoding : string, optional
            a string representing the encoding to use if the contents are
            non-ascii, for python versions prior to 3
        infer_datetime_format: boolean, default False
            If True and `parse_dates` is True for a column, try to infer the
            datetime format based on the first datetime string. If the format
            can be inferred, there often will be a large parsing speed-up.

        See also
        --------
        pandas.read_csv

        Returns
        -------
        y : Series
        """

        # We're calling `DataFrame.from_csv` in the implementation,
        # which will propagate a warning regarding `from_csv` deprecation.
        from pandas.core.frame import DataFrame
        df = DataFrame.from_csv(path, header=header, index_col=index_col,
                                sep=sep, parse_dates=parse_dates,
                                encoding=encoding,
                                infer_datetime_format=infer_datetime_format)
        result = df.iloc[:, 0]
        if header is None:
            result.index.name = result.name = None

        return result

    def to_csv(self, path=None, index=True, sep=",", na_rep='',
               float_format=None, header=False, index_label=None,
               mode='w', encoding=None, compression=None, date_format=None,
               decimal='.'):
        """
        Write Series to a comma-separated values (csv) file

        Parameters
        ----------
        path : string or file handle, default None
            File path or object, if None is provided the result is returned as
            a string.
        na_rep : string, default ''
            Missing data representation
        float_format : string, default None
            Format string for floating point numbers
        header : boolean, default False
            Write out series name
        index : boolean, default True
            Write row names (index)
        index_label : string or sequence, default None
            Column label for index column(s) if desired. If None is given, and
            `header` and `index` are True, then the index names are used. A
            sequence should be given if the DataFrame uses MultiIndex.
        mode : Python write mode, default 'w'
        sep : character, default ","
            Field delimiter for the output file.
        encoding : string, optional
            a string representing the encoding to use if the contents are
            non-ascii, for python versions prior to 3
        compression : string, optional
            A string representing the compression to use in the output file.
            Allowed values are 'gzip', 'bz2', 'zip', 'xz'. This input is only
            used when the first argument is a filename.
        date_format: string, default None
            Format string for datetime objects.
        decimal: string, default '.'
            Character recognized as decimal separator. E.g. use ',' for
            European data
        """
        from pandas.core.frame import DataFrame
        df = DataFrame(self)
        # result is only a string if no path provided, otherwise None
        result = df.to_csv(path, index=index, sep=sep, na_rep=na_rep,
                           float_format=float_format, header=header,
                           index_label=index_label, mode=mode,
                           encoding=encoding, compression=compression,
                           date_format=date_format, decimal=decimal)
        if path is None:
            return result

    @Appender(generic._shared_docs['to_excel'] % _shared_doc_kwargs)
    def to_excel(self, excel_writer, sheet_name='Sheet1', na_rep='',
                 float_format=None, columns=None, header=True, index=True,
                 index_label=None, startrow=0, startcol=0, engine=None,
                 merge_cells=True, encoding=None, inf_rep='inf', verbose=True):
        df = self.to_frame()
        df.to_excel(excel_writer=excel_writer, sheet_name=sheet_name,
                    na_rep=na_rep, float_format=float_format, columns=columns,
                    header=header, index=index, index_label=index_label,
                    startrow=startrow, startcol=startcol, engine=engine,
                    merge_cells=merge_cells, encoding=encoding,
                    inf_rep=inf_rep, verbose=verbose)

    @Appender(generic._shared_docs['isna'] % _shared_doc_kwargs)
    def isna(self):
        return super(Series, self).isna()

    @Appender(generic._shared_docs['isna'] % _shared_doc_kwargs)
    def isnull(self):
        return super(Series, self).isnull()

    @Appender(generic._shared_docs['notna'] % _shared_doc_kwargs)
    def notna(self):
        return super(Series, self).notna()

    @Appender(generic._shared_docs['notna'] % _shared_doc_kwargs)
    def notnull(self):
        return super(Series, self).notnull()

    def dropna(self, axis=0, inplace=False, **kwargs):
        """
        Return a new Series with missing values removed.

        See the :ref:`User Guide <missing_data>` for more on which values are
        considered missing, and how to work with missing data.

        Parameters
        ----------
        axis : {0 or 'index'}, default 0
            There is only one axis to drop values from.
        inplace : bool, default False
            If True, do operation inplace and return None.
        **kwargs
            Not in use.

        Returns
        -------
        Series
            Series with NA entries dropped from it.

        See Also
        --------
        Series.isna: Indicate missing values.
        Series.notna : Indicate existing (non-missing) values.
        Series.fillna : Replace missing values.
        DataFrame.dropna : Drop rows or columns which contain NA values.
        Index.dropna : Drop missing indices.

        Examples
        --------
        >>> ser = pd.Series([1., 2., np.nan])
        >>> ser
        0    1.0
        1    2.0
        2    NaN
        dtype: float64

        Drop NA values from a Series.

        >>> ser.dropna()
        0    1.0
        1    2.0
        dtype: float64

        Keep the Series with valid entries in the same variable.

        >>> ser.dropna(inplace=True)
        >>> ser
        0    1.0
        1    2.0
        dtype: float64

        Empty strings are not considered NA values. ``None`` is considered an
        NA value.

        >>> ser = pd.Series([np.NaN, 2, pd.NaT, '', None, 'I stay'])
        >>> ser
        0       NaN
        1         2
        2       NaT
        3
        4      None
        5    I stay
        dtype: object
        >>> ser.dropna()
        1         2
        3
        5    I stay
        dtype: object
        """
        inplace = validate_bool_kwarg(inplace, 'inplace')
        kwargs.pop('how', None)
        if kwargs:
            raise TypeError('dropna() got an unexpected keyword '
                            'argument "{0}"'.format(list(kwargs.keys())[0]))

        axis = self._get_axis_number(axis or 0)

        if self._can_hold_na:
            result = remove_na_arraylike(self)
            if inplace:
                self._update_inplace(result)
            else:
                return result
        else:
            if inplace:
                # do nothing
                pass
            else:
                return self.copy()

    def valid(self, inplace=False, **kwargs):
        """Return Series without null values.

        .. deprecated:: 0.23.0
            Use :meth:`Series.dropna` instead.
        """
        warnings.warn("Method .valid will be removed in a future version. "
                      "Use .dropna instead.", FutureWarning, stacklevel=2)
        return self.dropna(inplace=inplace, **kwargs)

    # ----------------------------------------------------------------------
    # Time series-oriented methods

    def to_timestamp(self, freq=None, how='start', copy=True):
        """
        Cast to datetimeindex of timestamps, at *beginning* of period

        Parameters
        ----------
        freq : string, default frequency of PeriodIndex
            Desired frequency
        how : {'s', 'e', 'start', 'end'}
            Convention for converting period to timestamp; start of period
            vs. end

        Returns
        -------
        ts : Series with DatetimeIndex
        """
        new_values = self._values
        if copy:
            new_values = new_values.copy()

        new_index = self.index.to_timestamp(freq=freq, how=how)
        return self._constructor(new_values,
                                 index=new_index).__finalize__(self)

    def to_period(self, freq=None, copy=True):
        """
        Convert Series from DatetimeIndex to PeriodIndex with desired
        frequency (inferred from index if not passed)

        Parameters
        ----------
        freq : string, default

        Returns
        -------
        ts : Series with PeriodIndex
        """
        new_values = self._values
        if copy:
            new_values = new_values.copy()

        new_index = self.index.to_period(freq=freq)
        return self._constructor(new_values,
                                 index=new_index).__finalize__(self)

    # ----------------------------------------------------------------------
    # Accessor Methods
    # ----------------------------------------------------------------------
    str = CachedAccessor("str", StringMethods)
    dt = CachedAccessor("dt", CombinedDatetimelikeProperties)
    cat = CachedAccessor("cat", CategoricalAccessor)
    plot = CachedAccessor("plot", gfx.SeriesPlotMethods)

    # ----------------------------------------------------------------------
    # Add plotting methods to Series
    hist = gfx.hist_series


Series._setup_axes(['index'], info_axis=0, stat_axis=0, aliases={'rows': 0},
                   docs={'index': 'The index (axis labels) of the Series.'})
Series._add_numeric_operations()
Series._add_series_only_operations()
Series._add_series_or_dataframe_operations()

# Add arithmetic!
ops.add_flex_arithmetic_methods(Series)
ops.add_special_arithmetic_methods(Series)


# -----------------------------------------------------------------------------
# Supplementary functions


def _sanitize_index(data, index, copy=False):
    """ sanitize an index type to return an ndarray of the underlying, pass
    thru a non-Index
    """

    if index is None:
        return data

    if len(data) != len(index):
        raise ValueError('Length of values does not match length of ' 'index')

    if isinstance(data, ABCIndexClass) and not copy:
        pass
    elif isinstance(data, (PeriodIndex, DatetimeIndex)):
        data = data._values
        if copy:
            data = data.copy()

    elif isinstance(data, np.ndarray):

        # coerce datetimelike types
        if data.dtype.kind in ['M', 'm']:
            data = _sanitize_array(data, index, copy=copy)

    return data


def _sanitize_array(data, index, dtype=None, copy=False,
                    raise_cast_failure=False):
    """ sanitize input data to an ndarray, copy if specified, coerce to the
    dtype if specified
    """

    if dtype is not None:
        dtype = pandas_dtype(dtype)

    if isinstance(data, ma.MaskedArray):
        mask = ma.getmaskarray(data)
        if mask.any():
            data, fill_value = maybe_upcast(data, copy=True)
            data[mask] = fill_value
        else:
            data = data.copy()

    def _try_cast(arr, take_fast_path):

        # perf shortcut as this is the most common case
        if take_fast_path:
            if maybe_castable(arr) and not copy and dtype is None:
                return arr

        try:
            subarr = maybe_cast_to_datetime(arr, dtype)
            # Take care in creating object arrays (but iterators are not
            # supported):
            if is_object_dtype(dtype) and (is_list_like(subarr) and
                                           not (is_iterator(subarr) or
                                           isinstance(subarr, np.ndarray))):
                subarr = construct_1d_object_array_from_listlike(subarr)
            elif not is_extension_type(subarr):
                subarr = np.array(subarr, dtype=dtype, copy=copy)
        except (ValueError, TypeError):
            if is_categorical_dtype(dtype):
                # We *do* allow casting to categorical, since we know
                # that Categorical is the only array type for 'category'.
                subarr = Categorical(arr, dtype.categories,
                                     ordered=dtype.ordered)
            elif is_extension_array_dtype(dtype):
                # We don't allow casting to third party dtypes, since we don't
                # know what array belongs to which type.
                msg = ("Cannot cast data to extension dtype '{}'. "
                       "Pass the extension array directly.".format(dtype))
                raise ValueError(msg)

            elif dtype is not None and raise_cast_failure:
                raise
            else:
                subarr = np.array(arr, dtype=object, copy=copy)
        return subarr

    # GH #846
    if isinstance(data, (np.ndarray, Index, Series)):

        if dtype is not None:
            subarr = np.array(data, copy=False)

            # possibility of nan -> garbage
            if is_float_dtype(data.dtype) and is_integer_dtype(dtype):
                if not isna(data).any():
                    subarr = _try_cast(data, True)
                elif copy:
                    subarr = data.copy()
            else:
                subarr = _try_cast(data, True)
        elif isinstance(data, Index):
            # don't coerce Index types
            # e.g. indexes can have different conversions (so don't fast path
            # them)
            # GH 6140
            subarr = _sanitize_index(data, index, copy=copy)
        else:

            # we will try to copy be-definition here
            subarr = _try_cast(data, True)

    elif isinstance(data, ExtensionArray):
        subarr = data

        if dtype is not None and not data.dtype.is_dtype(dtype):
            msg = ("Cannot coerce extension array to dtype '{typ}'. "
                   "Do the coercion before passing to the constructor "
                   "instead.".format(typ=dtype))
            raise ValueError(msg)

        if copy:
            subarr = data.copy()
        return subarr

    elif isinstance(data, (list, tuple)) and len(data) > 0:
        if dtype is not None:
            try:
                subarr = _try_cast(data, False)
            except Exception:
                if raise_cast_failure:  # pragma: no cover
                    raise
                subarr = np.array(data, dtype=object, copy=copy)
                subarr = lib.maybe_convert_objects(subarr)

        else:
            subarr = maybe_convert_platform(data)

        subarr = maybe_cast_to_datetime(subarr, dtype)

    elif isinstance(data, range):
        # GH 16804
        start, stop, step = get_range_parameters(data)
        arr = np.arange(start, stop, step, dtype='int64')
        subarr = _try_cast(arr, False)
    else:
        subarr = _try_cast(data, False)

    # scalar like, GH
    if getattr(subarr, 'ndim', 0) == 0:
        if isinstance(data, list):  # pragma: no cover
            subarr = np.array(data, dtype=object)
        elif index is not None:
            value = data

            # figure out the dtype from the value (upcast if necessary)
            if dtype is None:
                dtype, value = infer_dtype_from_scalar(value)
            else:
                # need to possibly convert the value here
                value = maybe_cast_to_datetime(value, dtype)

            subarr = construct_1d_arraylike_from_scalar(
                value, len(index), dtype)

        else:
            return subarr.item()

    # the result that we want
    elif subarr.ndim == 1:
        if index is not None:

            # a 1-element ndarray
            if len(subarr) != len(index) and len(subarr) == 1:
                subarr = construct_1d_arraylike_from_scalar(
                    subarr[0], len(index), subarr.dtype)

    elif subarr.ndim > 1:
        if isinstance(data, np.ndarray):
            raise Exception('Data must be 1-dimensional')
        else:
            subarr = com._asarray_tuplesafe(data, dtype=dtype)

    # This is to prevent mixed-type Series getting all casted to
    # NumPy string type, e.g. NaN --> '-1#IND'.
    if issubclass(subarr.dtype.type, compat.string_types):
        # GH 16605
        # If not empty convert the data to dtype
        # GH 19853: If data is a scalar, subarr has already the result
        if not is_scalar(data):
            if not np.all(isna(data)):
                data = np.array(data, dtype=dtype, copy=False)
            subarr = np.array(data, dtype=object, copy=copy)

    return subarr
