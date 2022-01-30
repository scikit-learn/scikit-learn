from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Any,
    Sequence,
    TypeVar,
    overload,
)
import warnings

import numpy as np

from pandas._libs import (
    iNaT,
    lib,
    missing as libmissing,
)
from pandas._typing import (
    ArrayLike,
    AstypeArg,
    NpDtype,
    PositionalIndexer,
    Scalar,
    ScalarIndexer,
    SequenceIndexer,
    Shape,
    npt,
    type_t,
)
from pandas.errors import AbstractMethodError
from pandas.util._decorators import (
    cache_readonly,
    doc,
)
from pandas.util._validators import validate_fillna_kwargs

from pandas.core.dtypes.base import ExtensionDtype
from pandas.core.dtypes.common import (
    is_bool,
    is_bool_dtype,
    is_dtype_equal,
    is_float,
    is_float_dtype,
    is_integer_dtype,
    is_list_like,
    is_numeric_dtype,
    is_object_dtype,
    is_scalar,
    is_string_dtype,
    pandas_dtype,
)
from pandas.core.dtypes.inference import is_array_like
from pandas.core.dtypes.missing import (
    array_equivalent,
    isna,
    notna,
)

from pandas.core import (
    arraylike,
    missing,
    nanops,
    ops,
)
from pandas.core.algorithms import (
    factorize_array,
    isin,
    take,
)
from pandas.core.array_algos import masked_reductions
from pandas.core.array_algos.quantile import quantile_with_mask
from pandas.core.arraylike import OpsMixin
from pandas.core.arrays import ExtensionArray
from pandas.core.indexers import check_array_indexer
from pandas.core.ops import invalid_comparison

if TYPE_CHECKING:
    from pandas import Series
    from pandas.core.arrays import BooleanArray

from pandas.compat.numpy import function as nv

BaseMaskedArrayT = TypeVar("BaseMaskedArrayT", bound="BaseMaskedArray")


class BaseMaskedDtype(ExtensionDtype):
    """
    Base class for dtypes for BasedMaskedArray subclasses.
    """

    name: str
    base = None
    type: type

    na_value = libmissing.NA

    @cache_readonly
    def numpy_dtype(self) -> np.dtype:
        """Return an instance of our numpy dtype"""
        return np.dtype(self.type)

    @cache_readonly
    def kind(self) -> str:
        return self.numpy_dtype.kind

    @cache_readonly
    def itemsize(self) -> int:
        """Return the number of bytes in this dtype"""
        return self.numpy_dtype.itemsize

    @classmethod
    def construct_array_type(cls) -> type_t[BaseMaskedArray]:
        """
        Return the array type associated with this dtype.

        Returns
        -------
        type
        """
        raise NotImplementedError


class BaseMaskedArray(OpsMixin, ExtensionArray):
    """
    Base class for masked arrays (which use _data and _mask to store the data).

    numpy based
    """

    # The value used to fill '_data' to avoid upcasting
    _internal_fill_value: Scalar
    # our underlying data and mask are each ndarrays
    _data: np.ndarray
    _mask: np.ndarray

    # Fill values used for any/all
    _truthy_value = Scalar  # bool(_truthy_value) = True
    _falsey_value = Scalar  # bool(_falsey_value) = False

    def __init__(self, values: np.ndarray, mask: np.ndarray, copy: bool = False):
        # values is supposed to already be validated in the subclass
        if not (isinstance(mask, np.ndarray) and mask.dtype == np.bool_):
            raise TypeError(
                "mask should be boolean numpy array. Use "
                "the 'pd.array' function instead"
            )
        if values.shape != mask.shape:
            raise ValueError("values.shape must match mask.shape")

        if copy:
            values = values.copy()
            mask = mask.copy()

        self._data = values
        self._mask = mask

    @property
    def dtype(self) -> BaseMaskedDtype:
        raise AbstractMethodError(self)

    @overload
    def __getitem__(self, item: ScalarIndexer) -> Any:
        ...

    @overload
    def __getitem__(self: BaseMaskedArrayT, item: SequenceIndexer) -> BaseMaskedArrayT:
        ...

    def __getitem__(
        self: BaseMaskedArrayT, item: PositionalIndexer
    ) -> BaseMaskedArrayT | Any:
        item = check_array_indexer(self, item)

        newmask = self._mask[item]
        if is_bool(newmask):
            # This is a scalar indexing
            if newmask:
                return self.dtype.na_value
            return self._data[item]

        return type(self)(self._data[item], newmask)

    @doc(ExtensionArray.fillna)
    def fillna(
        self: BaseMaskedArrayT, value=None, method=None, limit=None
    ) -> BaseMaskedArrayT:
        value, method = validate_fillna_kwargs(value, method)

        mask = self._mask

        if is_array_like(value):
            if len(value) != len(self):
                raise ValueError(
                    f"Length of 'value' does not match. Got ({len(value)}) "
                    f" expected {len(self)}"
                )
            value = value[mask]

        if mask.any():
            if method is not None:
                func = missing.get_fill_func(method, ndim=self.ndim)
                new_values, new_mask = func(
                    self._data.copy().T,
                    limit=limit,
                    mask=mask.copy().T,
                )
                return type(self)(new_values.T, new_mask.view(np.bool_).T)
            else:
                # fill with value
                new_values = self.copy()
                new_values[mask] = value
        else:
            new_values = self.copy()
        return new_values

    def _coerce_to_array(self, values) -> tuple[np.ndarray, np.ndarray]:
        raise AbstractMethodError(self)

    def __setitem__(self, key, value) -> None:
        _is_scalar = is_scalar(value)
        if _is_scalar:
            value = [value]
        value, mask = self._coerce_to_array(value)

        if _is_scalar:
            value = value[0]
            mask = mask[0]

        key = check_array_indexer(self, key)
        self._data[key] = value
        self._mask[key] = mask

    def __iter__(self):
        if self.ndim == 1:
            for i in range(len(self)):
                if self._mask[i]:
                    yield self.dtype.na_value
                else:
                    yield self._data[i]
        else:
            for i in range(len(self)):
                yield self[i]

    def __len__(self) -> int:
        return len(self._data)

    @property
    def shape(self) -> Shape:
        return self._data.shape

    @property
    def ndim(self) -> int:
        return self._data.ndim

    def swapaxes(self: BaseMaskedArrayT, axis1, axis2) -> BaseMaskedArrayT:
        data = self._data.swapaxes(axis1, axis2)
        mask = self._mask.swapaxes(axis1, axis2)
        return type(self)(data, mask)

    def delete(self: BaseMaskedArrayT, loc, axis: int = 0) -> BaseMaskedArrayT:
        data = np.delete(self._data, loc, axis=axis)
        mask = np.delete(self._mask, loc, axis=axis)
        return type(self)(data, mask)

    def reshape(self: BaseMaskedArrayT, *args, **kwargs) -> BaseMaskedArrayT:
        data = self._data.reshape(*args, **kwargs)
        mask = self._mask.reshape(*args, **kwargs)
        return type(self)(data, mask)

    def ravel(self: BaseMaskedArrayT, *args, **kwargs) -> BaseMaskedArrayT:
        # TODO: need to make sure we have the same order for data/mask
        data = self._data.ravel(*args, **kwargs)
        mask = self._mask.ravel(*args, **kwargs)
        return type(self)(data, mask)

    @property
    def T(self: BaseMaskedArrayT) -> BaseMaskedArrayT:
        return type(self)(self._data.T, self._mask.T)

    def __invert__(self: BaseMaskedArrayT) -> BaseMaskedArrayT:
        return type(self)(~self._data, self._mask.copy())

    def to_numpy(
        self,
        dtype: npt.DTypeLike | None = None,
        copy: bool = False,
        na_value: Scalar = lib.no_default,
    ) -> np.ndarray:
        """
        Convert to a NumPy Array.

        By default converts to an object-dtype NumPy array. Specify the `dtype` and
        `na_value` keywords to customize the conversion.

        Parameters
        ----------
        dtype : dtype, default object
            The numpy dtype to convert to.
        copy : bool, default False
            Whether to ensure that the returned value is a not a view on
            the array. Note that ``copy=False`` does not *ensure* that
            ``to_numpy()`` is no-copy. Rather, ``copy=True`` ensure that
            a copy is made, even if not strictly necessary. This is typically
            only possible when no missing values are present and `dtype`
            is the equivalent numpy dtype.
        na_value : scalar, optional
             Scalar missing value indicator to use in numpy array. Defaults
             to the native missing value indicator of this array (pd.NA).

        Returns
        -------
        numpy.ndarray

        Examples
        --------
        An object-dtype is the default result

        >>> a = pd.array([True, False, pd.NA], dtype="boolean")
        >>> a.to_numpy()
        array([True, False, <NA>], dtype=object)

        When no missing values are present, an equivalent dtype can be used.

        >>> pd.array([True, False], dtype="boolean").to_numpy(dtype="bool")
        array([ True, False])
        >>> pd.array([1, 2], dtype="Int64").to_numpy("int64")
        array([1, 2])

        However, requesting such dtype will raise a ValueError if
        missing values are present and the default missing value :attr:`NA`
        is used.

        >>> a = pd.array([True, False, pd.NA], dtype="boolean")
        >>> a
        <BooleanArray>
        [True, False, <NA>]
        Length: 3, dtype: boolean

        >>> a.to_numpy(dtype="bool")
        Traceback (most recent call last):
        ...
        ValueError: cannot convert to bool numpy array in presence of missing values

        Specify a valid `na_value` instead

        >>> a.to_numpy(dtype="bool", na_value=False)
        array([ True, False, False])
        """
        if na_value is lib.no_default:
            na_value = libmissing.NA
        if dtype is None:
            dtype = object
        if self._hasna:
            if (
                not is_object_dtype(dtype)
                and not is_string_dtype(dtype)
                and na_value is libmissing.NA
            ):
                raise ValueError(
                    f"cannot convert to '{dtype}'-dtype NumPy array "
                    "with missing values. Specify an appropriate 'na_value' "
                    "for this dtype."
                )
            # don't pass copy to astype -> always need a copy since we are mutating
            data = self._data.astype(dtype)
            data[self._mask] = na_value
        else:
            data = self._data.astype(dtype, copy=copy)
        return data

    @overload
    def astype(self, dtype: npt.DTypeLike, copy: bool = ...) -> np.ndarray:
        ...

    @overload
    def astype(self, dtype: ExtensionDtype, copy: bool = ...) -> ExtensionArray:
        ...

    @overload
    def astype(self, dtype: AstypeArg, copy: bool = ...) -> ArrayLike:
        ...

    def astype(self, dtype: AstypeArg, copy: bool = True) -> ArrayLike:
        dtype = pandas_dtype(dtype)

        if is_dtype_equal(dtype, self.dtype):
            if copy:
                return self.copy()
            return self

        # if we are astyping to another nullable masked dtype, we can fastpath
        if isinstance(dtype, BaseMaskedDtype):
            # TODO deal with NaNs for FloatingArray case
            data = self._data.astype(dtype.numpy_dtype, copy=copy)
            # mask is copied depending on whether the data was copied, and
            # not directly depending on the `copy` keyword
            mask = self._mask if data is self._data else self._mask.copy()
            cls = dtype.construct_array_type()
            return cls(data, mask, copy=False)

        if isinstance(dtype, ExtensionDtype):
            eacls = dtype.construct_array_type()
            return eacls._from_sequence(self, dtype=dtype, copy=copy)

        raise NotImplementedError("subclass must implement astype to np.dtype")

    __array_priority__ = 1000  # higher than ndarray so ops dispatch to us

    def __array__(self, dtype: NpDtype | None = None) -> np.ndarray:
        """
        the array interface, return my values
        We return an object array here to preserve our scalar values
        """
        return self.to_numpy(dtype=dtype)

    _HANDLED_TYPES: tuple[type, ...]

    def __array_ufunc__(self, ufunc: np.ufunc, method: str, *inputs, **kwargs):
        # For MaskedArray inputs, we apply the ufunc to ._data
        # and mask the result.

        out = kwargs.get("out", ())

        for x in inputs + out:
            if not isinstance(x, self._HANDLED_TYPES + (BaseMaskedArray,)):
                return NotImplemented

        # for binary ops, use our custom dunder methods
        result = ops.maybe_dispatch_ufunc_to_dunder_op(
            self, ufunc, method, *inputs, **kwargs
        )
        if result is not NotImplemented:
            return result

        if "out" in kwargs:
            # e.g. test_ufunc_with_out
            return arraylike.dispatch_ufunc_with_out(
                self, ufunc, method, *inputs, **kwargs
            )

        if method == "reduce":
            result = arraylike.dispatch_reduction_ufunc(
                self, ufunc, method, *inputs, **kwargs
            )
            if result is not NotImplemented:
                return result

        mask = np.zeros(len(self), dtype=bool)
        inputs2 = []
        for x in inputs:
            if isinstance(x, BaseMaskedArray):
                mask |= x._mask
                inputs2.append(x._data)
            else:
                inputs2.append(x)

        def reconstruct(x):
            # we don't worry about scalar `x` here, since we
            # raise for reduce up above.
            from pandas.core.arrays import (
                BooleanArray,
                FloatingArray,
                IntegerArray,
            )

            if is_bool_dtype(x.dtype):
                m = mask.copy()
                return BooleanArray(x, m)
            elif is_integer_dtype(x.dtype):
                m = mask.copy()
                return IntegerArray(x, m)
            elif is_float_dtype(x.dtype):
                m = mask.copy()
                if x.dtype == np.float16:
                    # reached in e.g. np.sqrt on BooleanArray
                    # we don't support float16
                    x = x.astype(np.float32)
                return FloatingArray(x, m)
            else:
                x[mask] = np.nan
            return x

        result = getattr(ufunc, method)(*inputs2, **kwargs)
        if ufunc.nout > 1:
            # e.g. np.divmod
            return tuple(reconstruct(x) for x in result)
        elif method == "reduce":
            # e.g. np.add.reduce; test_ufunc_reduce_raises
            if self._mask.any():
                return self._na_value
            return result
        else:
            return reconstruct(result)

    def __arrow_array__(self, type=None):
        """
        Convert myself into a pyarrow Array.
        """
        import pyarrow as pa

        return pa.array(self._data, mask=self._mask, type=type)

    @property
    def _hasna(self) -> bool:
        # Note: this is expensive right now! The hope is that we can
        # make this faster by having an optional mask, but not have to change
        # source code using it..

        # error: Incompatible return value type (got "bool_", expected "bool")
        return self._mask.any()  # type: ignore[return-value]

    def _cmp_method(self, other, op) -> BooleanArray:
        from pandas.core.arrays import BooleanArray

        mask = None

        if isinstance(other, BaseMaskedArray):
            other, mask = other._data, other._mask

        elif is_list_like(other):
            other = np.asarray(other)
            if other.ndim > 1:
                raise NotImplementedError("can only perform ops with 1-d structures")
            if len(self) != len(other):
                raise ValueError("Lengths must match to compare")

        if other is libmissing.NA:
            # numpy does not handle pd.NA well as "other" scalar (it returns
            # a scalar False instead of an array)
            # This may be fixed by NA.__array_ufunc__. Revisit this check
            # once that's implemented.
            result = np.zeros(self._data.shape, dtype="bool")
            mask = np.ones(self._data.shape, dtype="bool")
        else:
            with warnings.catch_warnings():
                # numpy may show a FutureWarning:
                #     elementwise comparison failed; returning scalar instead,
                #     but in the future will perform elementwise comparison
                # before returning NotImplemented. We fall back to the correct
                # behavior today, so that should be fine to ignore.
                warnings.filterwarnings("ignore", "elementwise", FutureWarning)
                with np.errstate(all="ignore"):
                    method = getattr(self._data, f"__{op.__name__}__")
                    result = method(other)

                if result is NotImplemented:
                    result = invalid_comparison(self._data, other, op)

        # nans propagate
        if mask is None:
            mask = self._mask.copy()
        else:
            mask = self._mask | mask

        return BooleanArray(result, mask, copy=False)

    def _maybe_mask_result(self, result, mask, other, op_name: str):
        """
        Parameters
        ----------
        result : array-like
        mask : array-like bool
        other : scalar or array-like
        op_name : str
        """
        # if we have a float operand we are by-definition
        # a float result
        # or our op is a divide
        if (
            (is_float_dtype(other) or is_float(other))
            or (op_name in ["rtruediv", "truediv"])
            or (is_float_dtype(self.dtype) and is_numeric_dtype(result.dtype))
        ):
            from pandas.core.arrays import FloatingArray

            return FloatingArray(result, mask, copy=False)

        elif is_bool_dtype(result):
            from pandas.core.arrays import BooleanArray

            return BooleanArray(result, mask, copy=False)

        elif result.dtype == "timedelta64[ns]":
            # e.g. test_numeric_arr_mul_tdscalar_numexpr_path
            from pandas.core.arrays import TimedeltaArray

            result[mask] = iNaT
            return TimedeltaArray._simple_new(result)

        elif is_integer_dtype(result):
            from pandas.core.arrays import IntegerArray

            return IntegerArray(result, mask, copy=False)

        else:
            result[mask] = np.nan
            return result

    def isna(self) -> np.ndarray:
        return self._mask.copy()

    @property
    def _na_value(self):
        return self.dtype.na_value

    @property
    def nbytes(self) -> int:
        return self._data.nbytes + self._mask.nbytes

    @classmethod
    def _concat_same_type(
        cls: type[BaseMaskedArrayT],
        to_concat: Sequence[BaseMaskedArrayT],
        axis: int = 0,
    ) -> BaseMaskedArrayT:
        data = np.concatenate([x._data for x in to_concat], axis=axis)
        mask = np.concatenate([x._mask for x in to_concat], axis=axis)
        return cls(data, mask)

    def take(
        self: BaseMaskedArrayT,
        indexer,
        *,
        allow_fill: bool = False,
        fill_value: Scalar | None = None,
        axis: int = 0,
    ) -> BaseMaskedArrayT:
        # we always fill with 1 internally
        # to avoid upcasting
        data_fill_value = self._internal_fill_value if isna(fill_value) else fill_value
        result = take(
            self._data,
            indexer,
            fill_value=data_fill_value,
            allow_fill=allow_fill,
            axis=axis,
        )

        mask = take(
            self._mask, indexer, fill_value=True, allow_fill=allow_fill, axis=axis
        )

        # if we are filling
        # we only fill where the indexer is null
        # not existing missing values
        # TODO(jreback) what if we have a non-na float as a fill value?
        if allow_fill and notna(fill_value):
            fill_mask = np.asarray(indexer) == -1
            result[fill_mask] = fill_value
            mask = mask ^ fill_mask

        return type(self)(result, mask, copy=False)

    # error: Return type "BooleanArray" of "isin" incompatible with return type
    # "ndarray" in supertype "ExtensionArray"
    def isin(self, values) -> BooleanArray:  # type: ignore[override]

        from pandas.core.arrays import BooleanArray

        # algorithms.isin will eventually convert values to an ndarray, so no extra
        # cost to doing it here first
        values_arr = np.asarray(values)
        result = isin(self._data, values_arr)

        if self._hasna:
            values_have_NA = is_object_dtype(values_arr.dtype) and any(
                val is self.dtype.na_value for val in values_arr
            )

            # For now, NA does not propagate so set result according to presence of NA,
            # see https://github.com/pandas-dev/pandas/pull/38379 for some discussion
            result[self._mask] = values_have_NA

        mask = np.zeros(self._data.shape, dtype=bool)
        return BooleanArray(result, mask, copy=False)

    def copy(self: BaseMaskedArrayT) -> BaseMaskedArrayT:
        data, mask = self._data, self._mask
        data = data.copy()
        mask = mask.copy()
        return type(self)(data, mask, copy=False)

    @doc(ExtensionArray.factorize)
    def factorize(self, na_sentinel: int = -1) -> tuple[np.ndarray, ExtensionArray]:
        arr = self._data
        mask = self._mask

        codes, uniques = factorize_array(arr, na_sentinel=na_sentinel, mask=mask)

        # the hashtables don't handle all different types of bits
        uniques = uniques.astype(self.dtype.numpy_dtype, copy=False)
        uniques_ea = type(self)(uniques, np.zeros(len(uniques), dtype=bool))
        return codes, uniques_ea

    def value_counts(self, dropna: bool = True) -> Series:
        """
        Returns a Series containing counts of each unique value.

        Parameters
        ----------
        dropna : bool, default True
            Don't include counts of missing values.

        Returns
        -------
        counts : Series

        See Also
        --------
        Series.value_counts
        """
        from pandas import (
            Index,
            Series,
        )
        from pandas.arrays import IntegerArray

        # compute counts on the data with no nans
        data = self._data[~self._mask]
        value_counts = Index(data).value_counts()

        index = value_counts.index

        # if we want nans, count the mask
        if dropna:
            counts = value_counts._values
        else:
            counts = np.empty(len(value_counts) + 1, dtype="int64")
            counts[:-1] = value_counts
            counts[-1] = self._mask.sum()

            index = index.insert(len(index), self.dtype.na_value)

        index = index.astype(self.dtype)

        mask = np.zeros(len(counts), dtype="bool")
        counts = IntegerArray(counts, mask)

        return Series(counts, index=index)

    @doc(ExtensionArray.equals)
    def equals(self, other) -> bool:
        if type(self) != type(other):
            return False
        if other.dtype != self.dtype:
            return False

        # GH#44382 if e.g. self[1] is np.nan and other[1] is pd.NA, we are NOT
        #  equal.
        if not np.array_equal(self._mask, other._mask):
            return False

        left = self._data[~self._mask]
        right = other._data[~other._mask]
        return array_equivalent(left, right, dtype_equal=True)

    def _quantile(
        self: BaseMaskedArrayT, qs: npt.NDArray[np.float64], interpolation: str
    ) -> BaseMaskedArrayT:
        """
        Dispatch to quantile_with_mask, needed because we do not have
        _from_factorized.

        Notes
        -----
        We assume that all impacted cases are 1D-only.
        """
        mask = np.atleast_2d(np.asarray(self.isna()))
        npvalues: np.ndarray = np.atleast_2d(np.asarray(self))

        res = quantile_with_mask(
            npvalues,
            mask=mask,
            fill_value=self.dtype.na_value,
            qs=qs,
            interpolation=interpolation,
        )
        assert res.ndim == 2
        assert res.shape[0] == 1
        res = res[0]
        try:
            out = type(self)._from_sequence(res, dtype=self.dtype)
        except TypeError:
            # GH#42626: not able to safely cast Int64
            # for floating point output
            out = np.asarray(res, dtype=np.float64)
        return out

    def _reduce(self, name: str, *, skipna: bool = True, **kwargs):
        if name in {"any", "all", "min", "max", "sum", "prod"}:
            return getattr(self, name)(skipna=skipna, **kwargs)

        data = self._data
        mask = self._mask

        if name in {"mean"}:
            op = getattr(masked_reductions, name)
            result = op(data, mask, skipna=skipna, **kwargs)
            return result

        # coerce to a nan-aware float if needed
        # (we explicitly use NaN within reductions)
        if self._hasna:
            data = self.to_numpy("float64", na_value=np.nan)

        # median, var, std, skew, kurt, idxmin, idxmax
        op = getattr(nanops, "nan" + name)
        result = op(data, axis=0, skipna=skipna, mask=mask, **kwargs)

        if np.isnan(result):
            return libmissing.NA

        return result

    def _wrap_reduction_result(self, name: str, result, skipna, **kwargs):
        if isinstance(result, np.ndarray):
            axis = kwargs["axis"]
            if skipna:
                # we only retain mask for all-NA rows/columns
                mask = self._mask.all(axis=axis)
            else:
                mask = self._mask.any(axis=axis)

            return self._maybe_mask_result(result, mask, other=None, op_name=name)
        return result

    def sum(self, *, skipna=True, min_count=0, axis: int | None = 0, **kwargs):
        nv.validate_sum((), kwargs)

        # TODO: do this in validate_sum?
        if "out" in kwargs:
            # np.sum; test_floating_array_numpy_sum
            if kwargs["out"] is not None:
                raise NotImplementedError
            kwargs.pop("out")

        result = masked_reductions.sum(
            self._data,
            self._mask,
            skipna=skipna,
            min_count=min_count,
            axis=axis,
        )
        return self._wrap_reduction_result(
            "sum", result, skipna=skipna, axis=axis, **kwargs
        )

    def prod(self, *, skipna=True, min_count=0, axis: int | None = 0, **kwargs):
        nv.validate_prod((), kwargs)
        result = masked_reductions.prod(
            self._data,
            self._mask,
            skipna=skipna,
            min_count=min_count,
            axis=axis,
        )
        return self._wrap_reduction_result(
            "prod", result, skipna=skipna, axis=axis, **kwargs
        )

    def min(self, *, skipna=True, axis: int | None = 0, **kwargs):
        nv.validate_min((), kwargs)
        return masked_reductions.min(
            self._data,
            self._mask,
            skipna=skipna,
            axis=axis,
        )

    def max(self, *, skipna=True, axis: int | None = 0, **kwargs):
        nv.validate_max((), kwargs)
        return masked_reductions.max(
            self._data,
            self._mask,
            skipna=skipna,
            axis=axis,
        )

    def any(self, *, skipna: bool = True, **kwargs):
        """
        Return whether any element is truthy.

        Returns False unless there is at least one element that is truthy.
        By default, NAs are skipped. If ``skipna=False`` is specified and
        missing values are present, similar :ref:`Kleene logic <boolean.kleene>`
        is used as for logical operations.

        .. versionchanged:: 1.4.0

        Parameters
        ----------
        skipna : bool, default True
            Exclude NA values. If the entire array is NA and `skipna` is
            True, then the result will be False, as for an empty array.
            If `skipna` is False, the result will still be True if there is
            at least one element that is truthy, otherwise NA will be returned
            if there are NA's present.
        **kwargs : any, default None
            Additional keywords have no effect but might be accepted for
            compatibility with NumPy.

        Returns
        -------
        bool or :attr:`pandas.NA`

        See Also
        --------
        numpy.any : Numpy version of this method.
        BaseMaskedArray.all : Return whether all elements are truthy.

        Examples
        --------
        The result indicates whether any element is truthy (and by default
        skips NAs):

        >>> pd.array([True, False, True]).any()
        True
        >>> pd.array([True, False, pd.NA]).any()
        True
        >>> pd.array([False, False, pd.NA]).any()
        False
        >>> pd.array([], dtype="boolean").any()
        False
        >>> pd.array([pd.NA], dtype="boolean").any()
        False
        >>> pd.array([pd.NA], dtype="Float64").any()
        False

        With ``skipna=False``, the result can be NA if this is logically
        required (whether ``pd.NA`` is True or False influences the result):

        >>> pd.array([True, False, pd.NA]).any(skipna=False)
        True
        >>> pd.array([1, 0, pd.NA]).any(skipna=False)
        True
        >>> pd.array([False, False, pd.NA]).any(skipna=False)
        <NA>
        >>> pd.array([0, 0, pd.NA]).any(skipna=False)
        <NA>
        """
        kwargs.pop("axis", None)
        nv.validate_any((), kwargs)

        values = self._data.copy()
        # Argument 3 to "putmask" has incompatible type "object"; expected
        # "Union[_SupportsArray[dtype[Any]], _NestedSequence[
        # _SupportsArray[dtype[Any]]], bool, int, float, complex, str, bytes, _Nested
        # Sequence[Union[bool, int, float, complex, str, bytes]]]"  [arg-type]
        np.putmask(values, self._mask, self._falsey_value)  # type: ignore[arg-type]
        result = values.any()
        if skipna:
            return result
        else:
            if result or len(self) == 0 or not self._mask.any():
                return result
            else:
                return self.dtype.na_value

    def all(self, *, skipna: bool = True, **kwargs):
        """
        Return whether all elements are truthy.

        Returns True unless there is at least one element that is falsey.
        By default, NAs are skipped. If ``skipna=False`` is specified and
        missing values are present, similar :ref:`Kleene logic <boolean.kleene>`
        is used as for logical operations.

        .. versionchanged:: 1.4.0

        Parameters
        ----------
        skipna : bool, default True
            Exclude NA values. If the entire array is NA and `skipna` is
            True, then the result will be True, as for an empty array.
            If `skipna` is False, the result will still be False if there is
            at least one element that is falsey, otherwise NA will be returned
            if there are NA's present.
        **kwargs : any, default None
            Additional keywords have no effect but might be accepted for
            compatibility with NumPy.

        Returns
        -------
        bool or :attr:`pandas.NA`

        See Also
        --------
        numpy.all : Numpy version of this method.
        BooleanArray.any : Return whether any element is truthy.

        Examples
        --------
        The result indicates whether all elements are truthy (and by default
        skips NAs):

        >>> pd.array([True, True, pd.NA]).all()
        True
        >>> pd.array([1, 1, pd.NA]).all()
        True
        >>> pd.array([True, False, pd.NA]).all()
        False
        >>> pd.array([], dtype="boolean").all()
        True
        >>> pd.array([pd.NA], dtype="boolean").all()
        True
        >>> pd.array([pd.NA], dtype="Float64").all()
        True

        With ``skipna=False``, the result can be NA if this is logically
        required (whether ``pd.NA`` is True or False influences the result):

        >>> pd.array([True, True, pd.NA]).all(skipna=False)
        <NA>
        >>> pd.array([1, 1, pd.NA]).all(skipna=False)
        <NA>
        >>> pd.array([True, False, pd.NA]).all(skipna=False)
        False
        >>> pd.array([1, 0, pd.NA]).all(skipna=False)
        False
        """
        kwargs.pop("axis", None)
        nv.validate_all((), kwargs)

        values = self._data.copy()
        # Argument 3 to "putmask" has incompatible type "object"; expected
        # "Union[_SupportsArray[dtype[Any]], _NestedSequence[
        # _SupportsArray[dtype[Any]]], bool, int, float, complex, str, bytes, _Neste
        # dSequence[Union[bool, int, float, complex, str, bytes]]]"  [arg-type]
        np.putmask(values, self._mask, self._truthy_value)  # type: ignore[arg-type]
        result = values.all()

        if skipna:
            return result
        else:
            if not result or len(self) == 0 or not self._mask.any():
                return result
            else:
                return self.dtype.na_value
