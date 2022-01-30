from __future__ import annotations

from functools import wraps
from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
    Sequence,
    TypeVar,
    cast,
    overload,
)

import numpy as np

from pandas._libs import lib
from pandas._libs.arrays import NDArrayBacked
from pandas._typing import (
    ArrayLike,
    Dtype,
    F,
    PositionalIndexer2D,
    PositionalIndexerTuple,
    ScalarIndexer,
    SequenceIndexer,
    Shape,
    TakeIndexer,
    npt,
    type_t,
)
from pandas.errors import AbstractMethodError
from pandas.util._decorators import doc
from pandas.util._validators import (
    validate_bool_kwarg,
    validate_fillna_kwargs,
    validate_insert_loc,
)

from pandas.core.dtypes.common import (
    is_dtype_equal,
    pandas_dtype,
)
from pandas.core.dtypes.dtypes import (
    DatetimeTZDtype,
    ExtensionDtype,
    PeriodDtype,
)
from pandas.core.dtypes.missing import array_equivalent

from pandas.core import missing
from pandas.core.algorithms import (
    take,
    unique,
    value_counts,
)
from pandas.core.array_algos.quantile import quantile_with_mask
from pandas.core.array_algos.transforms import shift
from pandas.core.arrays.base import ExtensionArray
from pandas.core.construction import extract_array
from pandas.core.indexers import check_array_indexer
from pandas.core.sorting import nargminmax

NDArrayBackedExtensionArrayT = TypeVar(
    "NDArrayBackedExtensionArrayT", bound="NDArrayBackedExtensionArray"
)

if TYPE_CHECKING:

    from pandas._typing import (
        NumpySorter,
        NumpyValueArrayLike,
    )


def ravel_compat(meth: F) -> F:
    """
    Decorator to ravel a 2D array before passing it to a cython operation,
    then reshape the result to our own shape.
    """

    @wraps(meth)
    def method(self, *args, **kwargs):
        if self.ndim == 1:
            return meth(self, *args, **kwargs)

        flags = self._ndarray.flags
        flat = self.ravel("K")
        result = meth(flat, *args, **kwargs)
        order = "F" if flags.f_contiguous else "C"
        return result.reshape(self.shape, order=order)

    return cast(F, method)


class NDArrayBackedExtensionArray(NDArrayBacked, ExtensionArray):
    """
    ExtensionArray that is backed by a single NumPy ndarray.
    """

    _ndarray: np.ndarray

    def _box_func(self, x):
        """
        Wrap numpy type in our dtype.type if necessary.
        """
        return x

    def _validate_scalar(self, value):
        # used by NDArrayBackedExtensionIndex.insert
        raise AbstractMethodError(self)

    # ------------------------------------------------------------------------

    def view(self, dtype: Dtype | None = None) -> ArrayLike:
        # We handle datetime64, datetime64tz, timedelta64, and period
        #  dtypes here. Everything else we pass through to the underlying
        #  ndarray.
        if dtype is None or dtype is self.dtype:
            return self._from_backing_data(self._ndarray)

        if isinstance(dtype, type):
            # we sometimes pass non-dtype objects, e.g np.ndarray;
            #  pass those through to the underlying ndarray
            return self._ndarray.view(dtype)

        dtype = pandas_dtype(dtype)
        arr = self._ndarray

        if isinstance(dtype, (PeriodDtype, DatetimeTZDtype)):
            cls = dtype.construct_array_type()
            return cls(arr.view("i8"), dtype=dtype)
        elif dtype == "M8[ns]":
            from pandas.core.arrays import DatetimeArray

            return DatetimeArray(arr.view("i8"), dtype=dtype)
        elif dtype == "m8[ns]":
            from pandas.core.arrays import TimedeltaArray

            return TimedeltaArray(arr.view("i8"), dtype=dtype)

        # error: Argument "dtype" to "view" of "_ArrayOrScalarCommon" has incompatible
        # type "Union[ExtensionDtype, dtype[Any]]"; expected "Union[dtype[Any], None,
        # type, _SupportsDType, str, Union[Tuple[Any, int], Tuple[Any, Union[int,
        # Sequence[int]]], List[Any], _DTypeDict, Tuple[Any, Any]]]"
        return arr.view(dtype=dtype)  # type: ignore[arg-type]

    def take(
        self: NDArrayBackedExtensionArrayT,
        indices: TakeIndexer,
        *,
        allow_fill: bool = False,
        fill_value: Any = None,
        axis: int = 0,
    ) -> NDArrayBackedExtensionArrayT:
        if allow_fill:
            fill_value = self._validate_scalar(fill_value)

        new_data = take(
            self._ndarray,
            indices,
            allow_fill=allow_fill,
            fill_value=fill_value,
            axis=axis,
        )
        return self._from_backing_data(new_data)

    # ------------------------------------------------------------------------

    def equals(self, other) -> bool:
        if type(self) is not type(other):
            return False
        if not is_dtype_equal(self.dtype, other.dtype):
            return False
        return bool(array_equivalent(self._ndarray, other._ndarray))

    def _values_for_argsort(self) -> np.ndarray:
        return self._ndarray

    # Signature of "argmin" incompatible with supertype "ExtensionArray"
    def argmin(self, axis: int = 0, skipna: bool = True):  # type:ignore[override]
        # override base class by adding axis keyword
        validate_bool_kwarg(skipna, "skipna")
        if not skipna and self.isna().any():
            raise NotImplementedError
        return nargminmax(self, "argmin", axis=axis)

    # Signature of "argmax" incompatible with supertype "ExtensionArray"
    def argmax(self, axis: int = 0, skipna: bool = True):  # type:ignore[override]
        # override base class by adding axis keyword
        validate_bool_kwarg(skipna, "skipna")
        if not skipna and self.isna().any():
            raise NotImplementedError
        return nargminmax(self, "argmax", axis=axis)

    def unique(self: NDArrayBackedExtensionArrayT) -> NDArrayBackedExtensionArrayT:
        new_data = unique(self._ndarray)
        return self._from_backing_data(new_data)

    @classmethod
    @doc(ExtensionArray._concat_same_type)
    def _concat_same_type(
        cls: type[NDArrayBackedExtensionArrayT],
        to_concat: Sequence[NDArrayBackedExtensionArrayT],
        axis: int = 0,
    ) -> NDArrayBackedExtensionArrayT:
        dtypes = {str(x.dtype) for x in to_concat}
        if len(dtypes) != 1:
            raise ValueError("to_concat must have the same dtype (tz)", dtypes)

        new_values = [x._ndarray for x in to_concat]
        new_values = np.concatenate(new_values, axis=axis)
        # error: Argument 1 to "_from_backing_data" of "NDArrayBackedExtensionArray" has
        # incompatible type "List[ndarray]"; expected "ndarray"
        return to_concat[0]._from_backing_data(new_values)  # type: ignore[arg-type]

    @doc(ExtensionArray.searchsorted)
    def searchsorted(
        self,
        value: NumpyValueArrayLike | ExtensionArray,
        side: Literal["left", "right"] = "left",
        sorter: NumpySorter = None,
    ) -> npt.NDArray[np.intp] | np.intp:
        npvalue = self._validate_searchsorted_value(value)
        return self._ndarray.searchsorted(npvalue, side=side, sorter=sorter)

    def _validate_searchsorted_value(
        self, value: NumpyValueArrayLike | ExtensionArray
    ) -> NumpyValueArrayLike:
        if isinstance(value, ExtensionArray):
            return value.to_numpy()
        else:
            return value

    @doc(ExtensionArray.shift)
    def shift(self, periods=1, fill_value=None, axis=0):

        fill_value = self._validate_shift_value(fill_value)
        new_values = shift(self._ndarray, periods, axis, fill_value)

        return self._from_backing_data(new_values)

    def _validate_shift_value(self, fill_value):
        # TODO(2.0): after deprecation in datetimelikearraymixin is enforced,
        #  we can remove this and use validate_fill_value directly
        return self._validate_scalar(fill_value)

    def __setitem__(self, key, value):
        key = check_array_indexer(self, key)
        value = self._validate_setitem_value(value)
        self._ndarray[key] = value

    def _validate_setitem_value(self, value):
        return value

    @overload
    def __getitem__(self, key: ScalarIndexer) -> Any:
        ...

    @overload
    def __getitem__(
        self: NDArrayBackedExtensionArrayT,
        key: SequenceIndexer | PositionalIndexerTuple,
    ) -> NDArrayBackedExtensionArrayT:
        ...

    def __getitem__(
        self: NDArrayBackedExtensionArrayT,
        key: PositionalIndexer2D,
    ) -> NDArrayBackedExtensionArrayT | Any:
        if lib.is_integer(key):
            # fast-path
            result = self._ndarray[key]
            if self.ndim == 1:
                return self._box_func(result)
            return self._from_backing_data(result)

        # error: Incompatible types in assignment (expression has type "ExtensionArray",
        # variable has type "Union[int, slice, ndarray]")
        key = extract_array(key, extract_numpy=True)  # type: ignore[assignment]
        key = check_array_indexer(self, key)
        result = self._ndarray[key]
        if lib.is_scalar(result):
            return self._box_func(result)

        result = self._from_backing_data(result)
        return result

    def _fill_mask_inplace(
        self, method: str, limit, mask: npt.NDArray[np.bool_]
    ) -> None:
        # (for now) when self.ndim == 2, we assume axis=0
        func = missing.get_fill_func(method, ndim=self.ndim)
        func(self._ndarray.T, limit=limit, mask=mask.T)
        return

    @doc(ExtensionArray.fillna)
    def fillna(
        self: NDArrayBackedExtensionArrayT, value=None, method=None, limit=None
    ) -> NDArrayBackedExtensionArrayT:
        value, method = validate_fillna_kwargs(
            value, method, validate_scalar_dict_value=False
        )

        mask = self.isna()
        # error: Argument 2 to "check_value_size" has incompatible type
        # "ExtensionArray"; expected "ndarray"
        value = missing.check_value_size(
            value, mask, len(self)  # type: ignore[arg-type]
        )

        if mask.any():
            if method is not None:
                # TODO: check value is None
                # (for now) when self.ndim == 2, we assume axis=0
                func = missing.get_fill_func(method, ndim=self.ndim)
                new_values, _ = func(self._ndarray.T.copy(), limit=limit, mask=mask.T)
                new_values = new_values.T

                # TODO: PandasArray didn't used to copy, need tests for this
                new_values = self._from_backing_data(new_values)
            else:
                # fill with value
                new_values = self.copy()
                new_values[mask] = value
        else:
            # We validate the fill_value even if there is nothing to fill
            if value is not None:
                self._validate_setitem_value(value)

            new_values = self.copy()
        return new_values

    # ------------------------------------------------------------------------
    # Reductions

    def _wrap_reduction_result(self, axis: int | None, result):
        if axis is None or self.ndim == 1:
            return self._box_func(result)
        return self._from_backing_data(result)

    # ------------------------------------------------------------------------
    # __array_function__ methods

    def _putmask(self, mask: npt.NDArray[np.bool_], value) -> None:
        """
        Analogue to np.putmask(self, mask, value)

        Parameters
        ----------
        mask : np.ndarray[bool]
        value : scalar or listlike

        Raises
        ------
        TypeError
            If value cannot be cast to self.dtype.
        """
        value = self._validate_setitem_value(value)

        np.putmask(self._ndarray, mask, value)

    def _where(
        self: NDArrayBackedExtensionArrayT, mask: np.ndarray, value
    ) -> NDArrayBackedExtensionArrayT:
        """
        Analogue to np.where(mask, self, value)

        Parameters
        ----------
        mask : np.ndarray[bool]
        value : scalar or listlike

        Raises
        ------
        TypeError
            If value cannot be cast to self.dtype.
        """
        value = self._validate_setitem_value(value)

        res_values = np.where(mask, self._ndarray, value)
        return self._from_backing_data(res_values)

    # ------------------------------------------------------------------------
    # Index compat methods

    def insert(
        self: NDArrayBackedExtensionArrayT, loc: int, item
    ) -> NDArrayBackedExtensionArrayT:
        """
        Make new ExtensionArray inserting new item at location. Follows
        Python list.append semantics for negative values.

        Parameters
        ----------
        loc : int
        item : object

        Returns
        -------
        type(self)
        """
        loc = validate_insert_loc(loc, len(self))

        code = self._validate_scalar(item)

        new_vals = np.concatenate(
            (
                self._ndarray[:loc],
                np.asarray([code], dtype=self._ndarray.dtype),
                self._ndarray[loc:],
            )
        )
        return self._from_backing_data(new_vals)

    # ------------------------------------------------------------------------
    # Additional array methods
    #  These are not part of the EA API, but we implement them because
    #  pandas assumes they're there.

    def value_counts(self, dropna: bool = True):
        """
        Return a Series containing counts of unique values.

        Parameters
        ----------
        dropna : bool, default True
            Don't include counts of NA values.

        Returns
        -------
        Series
        """
        if self.ndim != 1:
            raise NotImplementedError

        from pandas import (
            Index,
            Series,
        )

        if dropna:
            # error: Unsupported operand type for ~ ("ExtensionArray")
            values = self[~self.isna()]._ndarray  # type: ignore[operator]
        else:
            values = self._ndarray

        result = value_counts(values, sort=False, dropna=dropna)

        index_arr = self._from_backing_data(np.asarray(result.index._data))
        index = Index(index_arr, name=result.index.name)
        return Series(result._values, index=index, name=result.name)

    def _quantile(
        self: NDArrayBackedExtensionArrayT,
        qs: npt.NDArray[np.float64],
        interpolation: str,
    ) -> NDArrayBackedExtensionArrayT:
        # TODO: disable for Categorical if not ordered?

        # asarray needed for Sparse, see GH#24600
        mask = np.asarray(self.isna())
        mask = np.atleast_2d(mask)

        arr = np.atleast_2d(self._ndarray)
        # TODO: something NDArrayBacked-specific instead of _values_for_factorize[1]?
        fill_value = self._values_for_factorize()[1]

        res_values = quantile_with_mask(arr, mask, fill_value, qs, interpolation)

        result = type(self)._from_factorized(res_values, self)
        if self.ndim == 1:
            assert result.shape == (1, len(qs)), result.shape
            result = result[0]

        return result

    # ------------------------------------------------------------------------
    # numpy-like methods

    @classmethod
    def _empty(
        cls: type_t[NDArrayBackedExtensionArrayT], shape: Shape, dtype: ExtensionDtype
    ) -> NDArrayBackedExtensionArrayT:
        """
        Analogous to np.empty(shape, dtype=dtype)

        Parameters
        ----------
        shape : tuple[int]
        dtype : ExtensionDtype
        """
        # The base implementation uses a naive approach to find the dtype
        #  for the backing ndarray
        arr = cls._from_sequence([], dtype=dtype)
        backing = np.empty(shape, dtype=arr._ndarray.dtype)
        return arr._from_backing_data(backing)
