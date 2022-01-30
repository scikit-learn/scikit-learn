from __future__ import annotations

from typing import (
    Callable,
    Hashable,
)
import warnings

import numpy as np

from pandas._libs import (
    index as libindex,
    lib,
)
from pandas._typing import (
    Dtype,
    DtypeObj,
    npt,
)
from pandas.util._decorators import (
    cache_readonly,
    doc,
)
from pandas.util._exceptions import find_stack_level

from pandas.core.dtypes.cast import astype_nansafe
from pandas.core.dtypes.common import (
    is_dtype_equal,
    is_extension_array_dtype,
    is_float,
    is_float_dtype,
    is_integer_dtype,
    is_numeric_dtype,
    is_scalar,
    is_signed_integer_dtype,
    is_unsigned_integer_dtype,
    needs_i8_conversion,
    pandas_dtype,
)
from pandas.core.dtypes.generic import ABCSeries

from pandas.core.indexes.base import (
    Index,
    maybe_extract_name,
)


class NumericIndex(Index):
    """
    Immutable sequence used for indexing and alignment. The basic object
    storing axis labels for all pandas objects. NumericIndex is a special case
    of `Index` with purely numpy int/uint/float labels.

    .. versionadded:: 1.4.0

    Parameters
    ----------
    data : array-like (1-dimensional)
    dtype : NumPy dtype (default: None)
    copy : bool
        Make a copy of input ndarray.
    name : object
        Name to be stored in the index.

    Attributes
    ----------
    None

    Methods
    ----------
    None

    See Also
    --------
    Index : The base pandas Index type.
    Int64Index : Index of purely int64 labels (deprecated).
    UInt64Index : Index of purely uint64 labels (deprecated).
    Float64Index : Index of  purely float64 labels (deprecated).

    Notes
    -----
    An NumericIndex instance can **only** contain numpy int64/32/16/8, uint64/32/16/8 or
    float64/32/16 dtype. In particular, ``NumericIndex`` *can not* hold Pandas numeric
    dtypes (:class:`Int64Dtype`, :class:`Int32Dtype` etc.).
    """

    _typ = "numericindex"
    _values: np.ndarray
    _default_dtype: np.dtype | None = None
    _dtype_validation_metadata: tuple[Callable[..., bool], str] = (
        is_numeric_dtype,
        "numeric type",
    )
    _is_numeric_dtype = True
    _can_hold_strings = False
    _is_backward_compat_public_numeric_index: bool = True

    # error: Signature of "_can_hold_na" incompatible with supertype "Index"
    @cache_readonly
    def _can_hold_na(self) -> bool:  # type: ignore[override]
        if is_float_dtype(self.dtype):
            return True
        else:
            return False

    _engine_types: dict[np.dtype, type[libindex.IndexEngine]] = {
        np.dtype(np.int8): libindex.Int8Engine,
        np.dtype(np.int16): libindex.Int16Engine,
        np.dtype(np.int32): libindex.Int32Engine,
        np.dtype(np.int64): libindex.Int64Engine,
        np.dtype(np.uint8): libindex.UInt8Engine,
        np.dtype(np.uint16): libindex.UInt16Engine,
        np.dtype(np.uint32): libindex.UInt32Engine,
        np.dtype(np.uint64): libindex.UInt64Engine,
        np.dtype(np.float32): libindex.Float32Engine,
        np.dtype(np.float64): libindex.Float64Engine,
    }

    @property
    def _engine_type(self):
        # error: Invalid index type "Union[dtype[Any], ExtensionDtype]" for
        # "Dict[dtype[Any], Type[IndexEngine]]"; expected type "dtype[Any]"
        return self._engine_types[self.dtype]  # type: ignore[index]

    @cache_readonly
    def inferred_type(self) -> str:
        return {
            "i": "integer",
            "u": "integer",
            "f": "floating",
        }[self.dtype.kind]

    def __new__(cls, data=None, dtype: Dtype | None = None, copy=False, name=None):
        name = maybe_extract_name(name, data, cls)

        subarr = cls._ensure_array(data, dtype, copy)
        return cls._simple_new(subarr, name=name)

    @classmethod
    def _ensure_array(cls, data, dtype, copy: bool):
        """
        Ensure we have a valid array to pass to _simple_new.
        """
        cls._validate_dtype(dtype)

        if not isinstance(data, (np.ndarray, Index)):
            # Coerce to ndarray if not already ndarray or Index
            if is_scalar(data):
                raise cls._scalar_data_error(data)

            # other iterable of some kind
            if not isinstance(data, (ABCSeries, list, tuple)):
                data = list(data)

            orig = data
            data = np.asarray(data, dtype=dtype)
            if dtype is None and data.dtype.kind == "f":
                if cls is UInt64Index and (data >= 0).all():
                    # https://github.com/numpy/numpy/issues/19146
                    data = np.asarray(orig, dtype=np.uint64)

        if issubclass(data.dtype.type, str):
            cls._string_data_error(data)

        dtype = cls._ensure_dtype(dtype)

        if copy or not is_dtype_equal(data.dtype, dtype):
            # TODO: the try/except below is because it's difficult to predict the error
            # and/or error message from different combinations of data and dtype.
            # Efforts to avoid this try/except welcome.
            # See https://github.com/pandas-dev/pandas/pull/41153#discussion_r676206222
            try:
                subarr = np.array(data, dtype=dtype, copy=copy)
                cls._validate_dtype(subarr.dtype)
            except (TypeError, ValueError):
                raise ValueError(f"data is not compatible with {cls.__name__}")
            cls._assert_safe_casting(data, subarr)
        else:
            subarr = data

        if subarr.ndim > 1:
            # GH#13601, GH#20285, GH#27125
            raise ValueError("Index data must be 1-dimensional")

        subarr = np.asarray(subarr)
        return subarr

    @classmethod
    def _validate_dtype(cls, dtype: Dtype | None) -> None:
        if dtype is None:
            return

        validation_func, expected = cls._dtype_validation_metadata
        if not validation_func(dtype):
            raise ValueError(
                f"Incorrect `dtype` passed: expected {expected}, received {dtype}"
            )

    @classmethod
    def _ensure_dtype(cls, dtype: Dtype | None) -> np.dtype | None:
        """
        Ensure int64 dtype for Int64Index etc. but allow int32 etc. for NumericIndex.

        Assumes dtype has already been validated.
        """
        if dtype is None:
            return cls._default_dtype

        dtype = pandas_dtype(dtype)
        assert isinstance(dtype, np.dtype)

        if cls._is_backward_compat_public_numeric_index:
            # dtype for NumericIndex
            return dtype
        else:
            # dtype for Int64Index, UInt64Index etc. Needed for backwards compat.
            return cls._default_dtype

    def __contains__(self, key) -> bool:
        """
        Check if key is a float and has a decimal. If it has, return False.
        """
        if not is_integer_dtype(self.dtype):
            return super().__contains__(key)

        hash(key)
        try:
            if is_float(key) and int(key) != key:
                # otherwise the `key in self._engine` check casts e.g. 1.1 -> 1
                return False
            return key in self._engine
        except (OverflowError, TypeError, ValueError):
            return False

    @doc(Index.astype)
    def astype(self, dtype, copy: bool = True):
        dtype = pandas_dtype(dtype)
        if is_float_dtype(self.dtype):
            if needs_i8_conversion(dtype):
                raise TypeError(
                    f"Cannot convert Float64Index to dtype {dtype}; integer "
                    "values are required for conversion"
                )
            elif is_integer_dtype(dtype) and not is_extension_array_dtype(dtype):
                # TODO(ExtensionIndex); this can change once we have an EA Index type
                # GH 13149
                arr = astype_nansafe(self._values, dtype=dtype)
                if isinstance(self, Float64Index):
                    return Int64Index(arr, name=self.name)
                else:
                    return NumericIndex(arr, name=self.name, dtype=dtype)
        elif self._is_backward_compat_public_numeric_index:
            # this block is needed so e.g. NumericIndex[int8].astype("int32") returns
            # NumericIndex[int32] and not Int64Index with dtype int64.
            # When Int64Index etc. are removed from the code base, removed this also.
            if not is_extension_array_dtype(dtype) and is_numeric_dtype(dtype):
                return self._constructor(self, dtype=dtype, copy=copy)

        return super().astype(dtype, copy=copy)

    # ----------------------------------------------------------------
    # Indexing Methods

    # error: Decorated property not supported
    @cache_readonly  # type: ignore[misc]
    @doc(Index._should_fallback_to_positional)
    def _should_fallback_to_positional(self) -> bool:
        return False

    @doc(Index._convert_slice_indexer)
    def _convert_slice_indexer(self, key: slice, kind: str):
        if is_float_dtype(self.dtype):
            assert kind in ["loc", "getitem"]

            # We always treat __getitem__ slicing as label-based
            # translate to locations
            return self.slice_indexer(key.start, key.stop, key.step)

        return super()._convert_slice_indexer(key, kind=kind)

    @doc(Index._maybe_cast_slice_bound)
    def _maybe_cast_slice_bound(self, label, side: str, kind=lib.no_default):
        assert kind in ["loc", "getitem", None, lib.no_default]
        self._deprecated_arg(kind, "kind", "_maybe_cast_slice_bound")

        # we will try to coerce to integers
        return self._maybe_cast_indexer(label)

    # ----------------------------------------------------------------

    @doc(Index._shallow_copy)
    def _shallow_copy(self, values, name: Hashable = lib.no_default):
        if not self._can_hold_na and values.dtype.kind == "f":
            name = self._name if name is lib.no_default else name
            # Ensure we are not returning an Int64Index with float data:
            return Float64Index._simple_new(values, name=name)
        return super()._shallow_copy(values=values, name=name)

    def _convert_tolerance(self, tolerance, target):
        tolerance = super()._convert_tolerance(tolerance, target)

        if not np.issubdtype(tolerance.dtype, np.number):
            if tolerance.ndim > 0:
                raise ValueError(
                    f"tolerance argument for {type(self).__name__} must contain "
                    "numeric elements if it is list type"
                )
            else:
                raise ValueError(
                    f"tolerance argument for {type(self).__name__} must be numeric "
                    f"if it is a scalar: {repr(tolerance)}"
                )
        return tolerance

    def _is_comparable_dtype(self, dtype: DtypeObj) -> bool:
        # If we ever have BoolIndex or ComplexIndex, this may need to be tightened
        return is_numeric_dtype(dtype)

    @classmethod
    def _assert_safe_casting(cls, data: np.ndarray, subarr: np.ndarray) -> None:
        """
        Ensure incoming data can be represented with matching signed-ness.

        Needed if the process of casting data from some accepted dtype to the internal
        dtype(s) bears the risk of truncation (e.g. float to int).
        """
        if is_integer_dtype(subarr.dtype):
            if not np.array_equal(data, subarr):
                raise TypeError("Unsafe NumPy casting, you must explicitly cast")

    @property
    def _is_all_dates(self) -> bool:
        """
        Checks that all the labels are datetime objects.
        """
        return False

    def _format_native_types(
        self, *, na_rep="", float_format=None, decimal=".", quoting=None, **kwargs
    ):
        from pandas.io.formats.format import FloatArrayFormatter

        if is_float_dtype(self.dtype):
            formatter = FloatArrayFormatter(
                self._values,
                na_rep=na_rep,
                float_format=float_format,
                decimal=decimal,
                quoting=quoting,
                fixed_width=False,
            )
            return formatter.get_result_as_array()

        return super()._format_native_types(
            na_rep=na_rep,
            float_format=float_format,
            decimal=decimal,
            quoting=quoting,
            **kwargs,
        )


_num_index_shared_docs = {}


_num_index_shared_docs[
    "class_descr"
] = """
    Immutable sequence used for indexing and alignment. The basic object
    storing axis labels for all pandas objects. %(klass)s is a special case
    of `Index` with purely %(ltype)s labels. %(extra)s.

    .. deprecated:: 1.4.0
        In pandas v2.0 %(klass)s will be removed and :class:`NumericIndex` used instead.
        %(klass)s will remain fully functional for the duration of pandas 1.x.

    Parameters
    ----------
    data : array-like (1-dimensional)
    dtype : NumPy dtype (default: %(dtype)s)
    copy : bool
        Make a copy of input ndarray.
    name : object
        Name to be stored in the index.

    Attributes
    ----------
    None

    Methods
    ----------
    None

    See Also
    --------
    Index : The base pandas Index type.
    NumericIndex : Index of numpy int/uint/float data.

    Notes
    -----
    An Index instance can **only** contain hashable objects.
"""


class IntegerIndex(NumericIndex):
    """
    This is an abstract class for Int64Index, UInt64Index.
    """

    _is_backward_compat_public_numeric_index: bool = False

    @property
    def asi8(self) -> npt.NDArray[np.int64]:
        # do not cache or you'll create a memory leak
        warnings.warn(
            "Index.asi8 is deprecated and will be removed in a future version.",
            FutureWarning,
            stacklevel=find_stack_level(),
        )
        return self._values.view(self._default_dtype)

    def _validate_fill_value(self, value):
        # e.g. np.array([1.0]) we want np.array([1], dtype=self.dtype)
        #  see TestSetitemFloatNDarrayIntoIntegerSeries
        super()._validate_fill_value(value)
        if hasattr(value, "dtype") and is_float_dtype(value.dtype):
            converted = value.astype(self.dtype)
            if (converted == value).all():
                # See also: can_hold_element
                return converted
            raise TypeError
        return value


class Int64Index(IntegerIndex):
    _index_descr_args = {
        "klass": "Int64Index",
        "ltype": "integer",
        "dtype": "int64",
        "extra": "",
    }
    __doc__ = _num_index_shared_docs["class_descr"] % _index_descr_args

    _typ = "int64index"
    _engine_type = libindex.Int64Engine
    _default_dtype = np.dtype(np.int64)
    _dtype_validation_metadata = (is_signed_integer_dtype, "signed integer")


class UInt64Index(IntegerIndex):
    _index_descr_args = {
        "klass": "UInt64Index",
        "ltype": "unsigned integer",
        "dtype": "uint64",
        "extra": "",
    }
    __doc__ = _num_index_shared_docs["class_descr"] % _index_descr_args

    _typ = "uint64index"
    _engine_type = libindex.UInt64Engine
    _default_dtype = np.dtype(np.uint64)
    _dtype_validation_metadata = (is_unsigned_integer_dtype, "unsigned integer")

    def _validate_fill_value(self, value):
        # e.g. np.array([1]) we want np.array([1], dtype=np.uint64)
        #  see test_where_uin64
        super()._validate_fill_value(value)
        if hasattr(value, "dtype") and is_signed_integer_dtype(value.dtype):
            if (value >= 0).all():
                return value.astype(self.dtype)
            raise TypeError
        return value


class Float64Index(NumericIndex):
    _index_descr_args = {
        "klass": "Float64Index",
        "dtype": "float64",
        "ltype": "float",
        "extra": "",
    }
    __doc__ = _num_index_shared_docs["class_descr"] % _index_descr_args

    _typ = "float64index"
    _engine_type = libindex.Float64Engine
    _default_dtype = np.dtype(np.float64)
    _dtype_validation_metadata = (is_float_dtype, "float")
    _is_backward_compat_public_numeric_index: bool = False
