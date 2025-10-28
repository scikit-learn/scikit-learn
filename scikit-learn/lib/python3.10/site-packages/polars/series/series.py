from __future__ import annotations

import contextlib
import math
import os
import sys
from collections.abc import Iterable, Sequence
from contextlib import nullcontext
from datetime import date, datetime, time, timedelta
from decimal import Decimal as PyDecimal
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    ClassVar,
    Literal,
    NoReturn,
    Union,
    overload,
)

import polars._reexport as pl
from polars import functions as F
from polars._dependencies import (
    _ALTAIR_AVAILABLE,
    _PYARROW_AVAILABLE,
    _check_for_numpy,
    _check_for_pandas,
    _check_for_pyarrow,
    _check_for_torch,
    altair,
    import_optional,
    torch,
)
from polars._dependencies import numpy as np
from polars._dependencies import pandas as pd
from polars._dependencies import pyarrow as pa
from polars._utils.construction import (
    arrow_to_pyseries,
    dataframe_to_pyseries,
    iterable_to_pyseries,
    numpy_to_pyseries,
    pandas_to_pyseries,
    sequence_to_pyseries,
    series_to_pyseries,
)
from polars._utils.convert import (
    date_to_int,
    datetime_to_int,
    time_to_int,
    timedelta_to_int,
)
from polars._utils.deprecation import (
    deprecate_renamed_parameter,
    deprecated,
    issue_deprecation_warning,
)
from polars._utils.getitem import get_series_item_by_key
from polars._utils.unstable import unstable
from polars._utils.various import (
    BUILDING_SPHINX_DOCS,
    _is_generator,
    no_default,
    parse_version,
    qualified_type_name,
    require_same_type,
    scale_bytes,
    sphinx_accessor,
    warn_null_comparison,
)
from polars._utils.wrap import wrap_df, wrap_s
from polars.datatypes import (
    Array,
    Boolean,
    Categorical,
    Date,
    Datetime,
    Decimal,
    Duration,
    Enum,
    Float32,
    Float64,
    Int32,
    Int64,
    List,
    Null,
    Object,
    String,
    Time,
    UInt16,
    UInt32,
    UInt64,
    Unknown,
    is_polars_dtype,
    maybe_cast,
    numpy_char_code_to_dtype,
    parse_into_dtype,
    supported_numpy_char_code,
)
from polars.datatypes._utils import dtype_to_init_repr
from polars.exceptions import ComputeError, ModuleUpgradeRequiredError, ShapeError
from polars.interchange.protocol import CompatLevel
from polars.series.array import ArrayNameSpace
from polars.series.binary import BinaryNameSpace
from polars.series.categorical import CatNameSpace
from polars.series.datetime import DateTimeNameSpace
from polars.series.list import ListNameSpace
from polars.series.plotting import SeriesPlot
from polars.series.string import StringNameSpace
from polars.series.struct import StructNameSpace
from polars.series.utils import expr_dispatch, get_ffi_func

with contextlib.suppress(ImportError):  # Module not available when building docs
    from polars._plr import PyDataFrame, PySeries

if TYPE_CHECKING:
    with contextlib.suppress(ImportError):  # Module not available when building docs
        import polars._plr as plr

    from collections.abc import Collection, Generator, Mapping

    import jax
    import numpy.typing as npt

    from polars import DataFrame, DataType, Expr
    from polars._typing import (
        ArrowArrayExportable,
        ArrowStreamExportable,
        BufferInfo,
        ClosedInterval,
        ComparisonOperator,
        FillNullStrategy,
        InterpolationMethod,
        IntoExpr,
        IntoExprColumn,
        MultiIndexSelector,
        NonNestedLiteral,
        NullBehavior,
        NumericLiteral,
        PolarsDataType,
        PythonLiteral,
        QuantileMethod,
        RankMethod,
        RoundMode,
        SearchSortedSide,
        SeriesBuffers,
        SingleIndexSelector,
        SizeUnit,
        TemporalLiteral,
    )
    from polars._utils.various import NoDefault

    if sys.version_info >= (3, 11):
        from typing import Self
    else:
        from typing_extensions import Self

    if sys.version_info >= (3, 13):
        from warnings import deprecated
    else:
        from typing_extensions import deprecated  # noqa: TC004

elif BUILDING_SPHINX_DOCS:
    # note: we assign this way to work around an autocomplete issue in ipython/jedi
    # (ref: https://github.com/davidhalter/jedi/issues/2057)
    current_module = sys.modules[__name__]
    current_module.property = sphinx_accessor

ArrayLike = Union[
    Sequence[Any],
    "Series",
    "pa.Array",
    "pa.ChunkedArray",
    "np.ndarray[Any, Any]",
    "pd.Series[Any]",
    "pd.DatetimeIndex",
    "ArrowArrayExportable",
    "ArrowStreamExportable",
]


@expr_dispatch
class Series:
    """
    A Series represents a single column in a Polars DataFrame.

    Parameters
    ----------
    name : str, default None
        Name of the Series. Will be used as a column name when used in a DataFrame.
        When not specified, name is set to an empty string.
    values : ArrayLike, default None
        One-dimensional data in various forms. Supported are: Sequence, Series,
        pyarrow Array, and numpy ndarray.
    dtype : DataType, default None
        Data type of the resulting Series. If set to `None` (default), the data type is
        inferred from the `values` input. The strategy for data type inference depends
        on the `strict` parameter:

        - If `strict` is set to True (default), the inferred data type is equal to the
          first non-null value, or `Null` if all values are null.
        - If `strict` is set to False, the inferred data type is the supertype of the
          values, or :class:`Object` if no supertype can be found. **WARNING**: A full
          pass over the values is required to determine the supertype.
        - If no values were passed, the resulting data type is :class:`Null`.

    strict : bool, default True
        Throw an error if any value does not exactly match the given or inferred data
        type. If set to `False`, values that do not match the data type are cast to
        that data type or, if casting is not possible, set to null instead.
    nan_to_null : bool, default False
        In case a numpy array is used to create this Series, indicate how to deal
        with np.nan values. (This parameter is a no-op on non-numpy data).

    Examples
    --------
    Constructing a Series by specifying name and values positionally:

    >>> s = pl.Series("a", [1, 2, 3])
    >>> s
    shape: (3,)
    Series: 'a' [i64]
    [
            1
            2
            3
    ]

    Notice that the dtype is automatically inferred as a polars Int64:

    >>> s.dtype
    Int64

    Constructing a Series with a specific dtype:

    >>> s2 = pl.Series("a", [1, 2, 3], dtype=pl.Float32)
    >>> s2
    shape: (3,)
    Series: 'a' [f32]
    [
        1.0
        2.0
        3.0
    ]

    It is possible to construct a Series with values as the first positional argument.
    This syntax considered an anti-pattern, but it can be useful in certain
    scenarios. You must specify any other arguments through keywords.

    >>> s3 = pl.Series([1, 2, 3])
    >>> s3
    shape: (3,)
    Series: '' [i64]
    [
            1
            2
            3
    ]
    """

    # NOTE: This `= None` is needed to generate the docs with sphinx_accessor.
    _s: PySeries = None  # type: ignore[assignment]
    _accessors: ClassVar[set[str]] = {
        "arr",
        "bin",
        "cat",
        "dt",
        "list",
        "plot",
        "str",
        "struct",
    }

    def __init__(
        self,
        name: str | ArrayLike | None = None,
        values: ArrayLike | None = None,
        dtype: PolarsDataType | None = None,
        *,
        strict: bool = True,
        nan_to_null: bool = False,
    ) -> None:
        # If 'Unknown' treat as None to trigger type inference
        if dtype == Unknown:
            dtype = None
        elif dtype is not None and not is_polars_dtype(dtype):
            dtype = parse_into_dtype(dtype)

        # Handle case where values are passed as the first argument
        original_name: str | None = None
        if name is None:
            name = ""
        elif isinstance(name, str):
            original_name = name
        else:
            if values is None:
                values = name
                name = ""
            else:
                msg = "Series name must be a string"
                raise TypeError(msg)

        if isinstance(values, Sequence):
            self._s = sequence_to_pyseries(
                name,
                values,
                dtype=dtype,
                strict=strict,
                nan_to_null=nan_to_null,
            )

        elif values is None:
            self._s = sequence_to_pyseries(name, [], dtype=dtype)

        elif _check_for_numpy(values) and isinstance(values, np.ndarray):
            self._s = numpy_to_pyseries(
                name, values, strict=strict, nan_to_null=nan_to_null
            )
            if values.dtype.type in [np.datetime64, np.timedelta64]:
                # cast to appropriate dtype, handling NaT values
                input_dtype = _resolve_temporal_dtype(None, values.dtype)
                dtype = _resolve_temporal_dtype(dtype, values.dtype)
                if dtype is not None:
                    self._s = (
                        # `values.dtype` has already been validated in
                        # `numpy_to_pyseries`, so `input_dtype` can't be `None`
                        self.cast(input_dtype, strict=False)  # type: ignore[arg-type]
                        .cast(dtype)
                        .scatter(np.argwhere(np.isnat(values)).flatten(), None)
                        ._s
                    )
                    return

            if dtype is not None:
                self._s = self.cast(dtype, strict=strict)._s

        elif _check_for_torch(values) and isinstance(values, torch.Tensor):
            self._s = numpy_to_pyseries(
                name, values.numpy(force=False), strict=strict, nan_to_null=nan_to_null
            )
            if dtype is not None:
                self._s = self.cast(dtype, strict=strict)._s

        elif _check_for_pyarrow(values) and isinstance(
            values, (pa.Array, pa.ChunkedArray)
        ):
            self._s = arrow_to_pyseries(name, values, dtype=dtype, strict=strict)

        elif _check_for_pandas(values) and isinstance(
            values, (pd.Series, pd.Index, pd.DatetimeIndex)
        ):
            self._s = pandas_to_pyseries(name, values, dtype=dtype, strict=strict)

        elif not hasattr(values, "__arrow_c_stream__") and _is_generator(values):
            self._s = iterable_to_pyseries(name, values, dtype=dtype, strict=strict)

        elif isinstance(values, Series):
            self._s = series_to_pyseries(
                original_name, values, dtype=dtype, strict=strict
            )

        elif isinstance(values, pl.DataFrame):
            self._s = dataframe_to_pyseries(
                original_name, values, dtype=dtype, strict=strict
            )

        elif hasattr(values, "__arrow_c_array__"):
            self._s = PySeries.from_arrow_c_array(values)

        elif hasattr(values, "__arrow_c_stream__"):
            self._s = PySeries.from_arrow_c_stream(values)

        else:
            msg = (
                f"Series constructor called with unsupported type {type(values).__name__!r}"
                " for the `values` parameter"
            )
            raise TypeError(msg)

    @classmethod
    def _from_pyseries(cls, pyseries: PySeries) -> Self:
        series = cls.__new__(cls)
        series._s = pyseries
        return series

    @classmethod
    @deprecated(
        "`_import_from_c` is deprecated; use `_import_arrow_from_c` instead. If "
        "you are using an extension, please compile it with the latest 'pyo3-polars'"
    )
    def _import_from_c(cls, name: str, pointers: list[tuple[int, int]]) -> Self:
        # `_import_from_c` was deprecated in 1.3
        return cls._from_pyseries(PySeries._import_arrow_from_c(name, pointers))

    @classmethod
    def _import_arrow_from_c(cls, name: str, pointers: list[tuple[int, int]]) -> Self:
        """
        Construct a Series from Arrows C interface.

        Parameters
        ----------
        name
            The name that should be given to the `Series`.
        pointers
            A list with tuples containing two entries:
             - The raw pointer to a C ArrowArray struct
             - The raw pointer to a C ArrowSchema struct

        Warning
        -------
        This will read the `array` pointer without moving it. The host process should
        garbage collect the heap pointer, but not its contents.
        """
        return cls._from_pyseries(PySeries._import_arrow_from_c(name, pointers))

    @classmethod
    def _import(cls, pointer: int) -> Self:
        return cls._from_pyseries(PySeries._import(pointer))

    def _export_arrow_to_c(self, out_ptr: int, out_schema_ptr: int) -> None:
        """
        Export to a C ArrowArray and C ArrowSchema struct, given their pointers.

        Parameters
        ----------
        out_ptr: int
            The raw pointer to a C ArrowArray struct.
        out_schema_ptr: int (optional)
            The raw pointer to a C ArrowSchema struct.

        Notes
        -----
        The series should only contain a single chunk. If you want to export all chunks,
        first call `Series.get_chunks` to give you a list of chunks.

        Warning
        -------
        Safety
        This function will write to the pointers given in `out_ptr` and `out_schema_ptr`
        and thus is highly unsafe.

        Leaking
        If you don't pass the ArrowArray struct to a consumer,
        array memory will leak. This is a low-level function intended for
        expert users.
        """
        self._s._export_arrow_to_c(out_ptr, out_schema_ptr)

    def _get_buffer_info(self) -> BufferInfo:
        """
        Return pointer, offset, and length information about the underlying buffer.

        Returns
        -------
        tuple of ints
            Tuple of the form (pointer, offset, length)

        Raises
        ------
        TypeError
            If the `Series` data type is not physical.
        ComputeError
            If the `Series` contains multiple chunks.

        Notes
        -----
        This method is mainly intended for use with the dataframe interchange protocol.
        """
        return self._s._get_buffer_info()

    def _get_buffers(self) -> SeriesBuffers:
        """
        Return the underlying values, validity, and offsets buffers as Series.

        The values buffer always exists.
        The validity buffer may not exist if the column contains no null values.
        The offsets buffer only exists for Series of data type `String` and `List`.

        Returns
        -------
        dict
            Dictionary with `"values"`, `"validity"`, and `"offsets"` keys mapping
            to the corresponding buffer or `None` if the buffer doesn't exist.

        Warnings
        --------
        The underlying buffers for `String` Series cannot be represented in this
        format. Instead, the buffers are converted to a values and offsets buffer.

        Notes
        -----
        This method is mainly intended for use with the dataframe interchange protocol.
        """
        buffers = self._s._get_buffers()
        keys = ("values", "validity", "offsets")
        return {  # type: ignore[return-value]
            k: self._from_pyseries(b) if b is not None else b
            for k, b in zip(keys, buffers)
        }

    @classmethod
    def _from_buffer(
        cls, dtype: PolarsDataType, buffer_info: BufferInfo, owner: Any
    ) -> Self:
        """
        Construct a Series from information about its underlying buffer.

        Parameters
        ----------
        dtype
            The data type of the buffer.
            Must be a physical type (integer, float, or boolean).
        buffer_info
            Tuple containing buffer information in the form `(pointer, offset, length)`.
        owner
            The object owning the buffer.

        Returns
        -------
        Series

        Raises
        ------
        TypeError
            When the given `dtype` is not supported.

        Notes
        -----
        This method is mainly intended for use with the dataframe interchange protocol.
        """
        return cls._from_pyseries(PySeries._from_buffer(dtype, buffer_info, owner))

    @classmethod
    def _from_buffers(
        cls,
        dtype: PolarsDataType,
        data: Series | Sequence[Series],
        validity: Series | None = None,
    ) -> Self:
        """
        Construct a Series from information about its underlying buffers.

        Parameters
        ----------
        dtype
            The data type of the resulting Series.
        data
            Buffers describing the data. For most data types, this is a single Series of
            the physical data type of `dtype`. Some data types require multiple buffers:

            - `String`: A data buffer of type `UInt8` and an offsets buffer
              of type `Int64`. Note that this does not match how the data
              is represented internally and data copy is required to construct
              the Series.
        validity
            Validity buffer. If specified, must be a Series of data type `Boolean`.

        Returns
        -------
        Series

        Raises
        ------
        TypeError
            When the given `dtype` is not supported or the other inputs do not match
            the requirements for constructing a Series of the given `dtype`.

        Warnings
        --------
        Constructing a `String` Series requires specifying a values and offsets buffer,
        which does not match the actual underlying buffers. The values and offsets
        buffer are converted into the actual buffers, which copies data.

        Notes
        -----
        This method is mainly intended for use with the dataframe interchange protocol.
        """
        if isinstance(data, Series):
            data_lst = [data._s]
        else:
            data_lst = [s._s for s in data]
        validity_series: plr.PySeries | None = None
        if validity is not None:
            validity_series = validity._s
        return cls._from_pyseries(
            PySeries._from_buffers(dtype, data_lst, validity_series)
        )

    @staticmethod
    def _newest_compat_level() -> int:
        """
        Get the newest supported compat level.

        This is for pyo3-polars.
        """
        return CompatLevel._newest()._version

    @property
    def dtype(self) -> DataType:
        """
        Get the data type of this Series.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.dtype
        Int64
        """
        return self._s.dtype()

    @property
    def flags(self) -> dict[str, bool]:
        """
        Get flags that are set on the Series.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.flags
        {'SORTED_ASC': False, 'SORTED_DESC': False}
        """
        out = {
            "SORTED_ASC": self._s.is_sorted_ascending_flag(),
            "SORTED_DESC": self._s.is_sorted_descending_flag(),
        }
        if self.dtype == List:
            out["FAST_EXPLODE"] = self._s.can_fast_explode_flag()
        return out

    @property
    def name(self) -> str:
        """
        Get the name of this Series.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.name
        'a'
        """
        return self._s.name()

    @property
    def shape(self) -> tuple[int]:
        """
        Shape of this Series.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.shape
        (3,)
        """
        return (self._s.len(),)

    def __bool__(self) -> NoReturn:
        msg = (
            "the truth value of a Series is ambiguous"
            "\n\n"
            "Here are some things you might want to try:\n"
            "- instead of `if s`, use `if not s.is_empty()`\n"
            "- instead of `s1 and s2`, use `s1 & s2`\n"
            "- instead of `s1 or s2`, use `s1 | s2`\n"
            "- instead of `s in [y, z]`, use `s.is_in([y, z])`\n"
        )
        raise TypeError(msg)

    def __getstate__(self) -> bytes:
        return self._s.__getstate__()

    def __setstate__(self, state: bytes) -> None:
        self._s = Series()._s  # Initialize with a dummy
        self._s.__setstate__(state)

    def __str__(self) -> str:
        s_repr: str = self._s.as_str()
        return s_repr.replace("Series", f"{self.__class__.__name__}", 1)

    def __repr__(self) -> str:
        return self.__str__()

    def __len__(self) -> int:
        return self.len()

    @overload
    def __and__(self, other: Expr) -> Expr: ...

    @overload
    def __and__(self, other: Any) -> Series: ...

    def __and__(self, other: Any) -> Expr | Series:
        if isinstance(other, pl.Expr):
            return F.lit(self) & other
        if not isinstance(other, Series):
            other = Series([other])
        return self._from_pyseries(self._s.bitand(other._s))

    @overload
    def __rand__(self, other: Expr) -> Expr: ...

    @overload
    def __rand__(self, other: Any) -> Series: ...

    def __rand__(self, other: Any) -> Expr | Series:
        if isinstance(other, pl.Expr):
            return other & F.lit(self)
        if not isinstance(other, Series):
            other = Series([other])
        return other & self

    @overload
    def __or__(self, other: Expr) -> Expr: ...

    @overload
    def __or__(self, other: Any) -> Series: ...

    def __or__(self, other: Any) -> Expr | Series:
        if isinstance(other, pl.Expr):
            return F.lit(self) | other
        if not isinstance(other, Series):
            other = Series([other])
        return self._from_pyseries(self._s.bitor(other._s))

    @overload
    def __ror__(self, other: Expr) -> Expr: ...

    @overload
    def __ror__(self, other: Any) -> Series: ...

    def __ror__(self, other: Any) -> Expr | Series:
        if isinstance(other, pl.Expr):
            return other | F.lit(self)
        if not isinstance(other, Series):
            other = Series([other])
        return other | self

    @overload
    def __xor__(self, other: Expr) -> Expr: ...

    @overload
    def __xor__(self, other: Any) -> Series: ...

    def __xor__(self, other: Any) -> Expr | Series:
        if isinstance(other, pl.Expr):
            return F.lit(self) ^ other
        if not isinstance(other, Series):
            other = Series([other])
        return self._from_pyseries(self._s.bitxor(other._s))

    @overload
    def __rxor__(self, other: Expr) -> Expr: ...

    @overload
    def __rxor__(self, other: Any) -> Series: ...

    def __rxor__(self, other: Any) -> Expr | Series:
        if isinstance(other, pl.Expr):
            return other ^ F.lit(self)
        if not isinstance(other, Series):
            other = Series([other])
        return other ^ self

    def _comp(self, other: Any, op: ComparisonOperator) -> Series:
        # special edge-case; boolean broadcast series (eq/neq) is its own result
        if self.dtype == Boolean and isinstance(other, bool) and op in ("eq", "neq"):
            if (other is True and op == "eq") or (other is False and op == "neq"):
                return self.clone()
            elif (other is False and op == "eq") or (other is True and op == "neq"):
                return ~self

        elif isinstance(other, float) and self.dtype.is_integer():
            # require upcast when comparing int series to float value
            self = self.cast(Float64)
            f = get_ffi_func(op + "_<>", Float64, self._s)
            assert f is not None
            return self._from_pyseries(f(other))

        elif isinstance(other, datetime):
            if self.dtype == Date:
                # require upcast when comparing date series to datetime
                self = self.cast(Datetime("us"))
                time_unit = "us"
            elif self.dtype == Datetime:
                # Use local time zone info
                time_zone = self.dtype.time_zone  # type: ignore[attr-defined]
                if str(other.tzinfo) != str(time_zone):
                    msg = f"datetime time zone {other.tzinfo!r} does not match Series timezone {time_zone!r}"
                    raise TypeError(msg)
                time_unit = self.dtype.time_unit  # type: ignore[attr-defined]
            else:
                msg = f"cannot compare datetime.datetime to Series of type {self.dtype}"
                raise ValueError(msg)
            ts = datetime_to_int(other, time_unit)  # type: ignore[arg-type]
            f = get_ffi_func(op + "_<>", Int64, self._s)
            assert f is not None
            return self._from_pyseries(f(ts))

        elif isinstance(other, time) and self.dtype == Time:
            d = time_to_int(other)
            f = get_ffi_func(op + "_<>", Int64, self._s)
            assert f is not None
            return self._from_pyseries(f(d))

        elif isinstance(other, timedelta) and self.dtype == Duration:
            time_unit = self.dtype.time_unit  # type: ignore[attr-defined]
            td = timedelta_to_int(other, time_unit)  # type: ignore[arg-type]
            f = get_ffi_func(op + "_<>", Int64, self._s)
            assert f is not None
            return self._from_pyseries(f(td))

        elif self.dtype in [Categorical, Enum] and not isinstance(other, Series):
            other = Series([other])

        elif isinstance(other, date) and self.dtype == Date:
            d = date_to_int(other)
            f = get_ffi_func(op + "_<>", Int32, self._s)
            assert f is not None
            return self._from_pyseries(f(d))

        if isinstance(other, Sequence) and not isinstance(other, str):
            if self.dtype in (List, Array):
                other = [other]
            other = Series("", other)
            if other.dtype == Null:
                other.cast(self.dtype)

        if isinstance(other, Series):
            return self._from_pyseries(getattr(self._s, op)(other._s))
        try:
            f = get_ffi_func(op + "_<>", self.dtype, self._s)
        except NotImplementedError:
            f = None
        if f is None:
            msg = f"Series of type {self.dtype} does not have {op} operator"
            raise NotImplementedError(msg)
        if other is not None:
            other = maybe_cast(other, self.dtype)

        return self._from_pyseries(f(other))

    @overload  # type: ignore[override]
    def __eq__(self, other: Expr) -> Expr: ...  # type: ignore[overload-overlap]

    @overload
    def __eq__(self, other: object) -> Series: ...

    def __eq__(self, other: object) -> Series | Expr:
        warn_null_comparison(other)
        if isinstance(other, pl.Expr):
            return F.lit(self).__eq__(other)
        return self._comp(other, "eq")

    @overload  # type: ignore[override]
    def __ne__(self, other: Expr) -> Expr: ...  # type: ignore[overload-overlap]

    @overload
    def __ne__(self, other: object) -> Series: ...

    def __ne__(self, other: object) -> Series | Expr:
        warn_null_comparison(other)
        if isinstance(other, pl.Expr):
            return F.lit(self).__ne__(other)
        return self._comp(other, "neq")

    @overload
    def __gt__(self, other: Expr) -> Expr: ...

    @overload
    def __gt__(self, other: Any) -> Series: ...

    def __gt__(self, other: Any) -> Series | Expr:
        warn_null_comparison(other)
        if isinstance(other, pl.Expr):
            return F.lit(self).__gt__(other)
        return self._comp(other, "gt")

    @overload
    def __lt__(self, other: Expr) -> Expr: ...

    @overload
    def __lt__(self, other: Any) -> Series: ...

    def __lt__(self, other: Any) -> Series | Expr:
        warn_null_comparison(other)
        if isinstance(other, pl.Expr):
            return F.lit(self).__lt__(other)
        return self._comp(other, "lt")

    @overload
    def __ge__(self, other: Expr) -> Expr: ...

    @overload
    def __ge__(self, other: Any) -> Series: ...

    def __ge__(self, other: Any) -> Series | Expr:
        warn_null_comparison(other)
        if isinstance(other, pl.Expr):
            return F.lit(self).__ge__(other)
        return self._comp(other, "gt_eq")

    @overload
    def __le__(self, other: Expr) -> Expr: ...

    @overload
    def __le__(self, other: Any) -> Series: ...

    def __le__(self, other: Any) -> Series | Expr:
        warn_null_comparison(other)
        if isinstance(other, pl.Expr):
            return F.lit(self).__le__(other)
        return self._comp(other, "lt_eq")

    @overload
    def le(self, other: Expr) -> Expr: ...

    @overload
    def le(self, other: Any) -> Series: ...

    def le(self, other: Any) -> Series | Expr:
        """Method equivalent of operator expression `series <= other`."""
        return self.__le__(other)

    @overload
    def lt(self, other: Expr) -> Expr: ...

    @overload
    def lt(self, other: Any) -> Series: ...

    def lt(self, other: Any) -> Series | Expr:
        """Method equivalent of operator expression `series < other`."""
        return self.__lt__(other)

    @overload
    def eq(self, other: Expr) -> Expr: ...

    @overload
    def eq(self, other: Any) -> Series: ...

    def eq(self, other: Any) -> Series | Expr:
        """Method equivalent of operator expression `series == other`."""
        return self.__eq__(other)

    @overload
    def eq_missing(self, other: Expr) -> Expr: ...

    @overload
    def eq_missing(self, other: Any) -> Series: ...

    def eq_missing(self, other: Any) -> Series | Expr:
        """
        Method equivalent of equality operator `series == other` where `None == None`.

        This differs from the standard `eq` where null values are propagated.

        Parameters
        ----------
        other
            A literal or expression value to compare with.

        See Also
        --------
        ne_missing
        eq

        Examples
        --------
        >>> s1 = pl.Series("a", [333, 200, None])
        >>> s2 = pl.Series("a", [100, 200, None])
        >>> s1.eq(s2)
        shape: (3,)
        Series: 'a' [bool]
        [
            false
            true
            null
        ]
        >>> s1.eq_missing(s2)
        shape: (3,)
        Series: 'a' [bool]
        [
            false
            true
            true
        ]
        """
        if isinstance(other, pl.Expr):
            return F.lit(self).eq_missing(other)
        return self.to_frame().select(F.col(self.name).eq_missing(other)).to_series()

    @overload
    def ne(self, other: Expr) -> Expr: ...

    @overload
    def ne(self, other: Any) -> Series: ...

    def ne(self, other: Any) -> Series | Expr:
        """Method equivalent of operator expression `series != other`."""
        return self.__ne__(other)

    @overload
    def ne_missing(self, other: Expr) -> Expr: ...

    @overload
    def ne_missing(self, other: Any) -> Series: ...

    def ne_missing(self, other: Any) -> Series | Expr:
        """
        Method equivalent of equality operator `series != other` where `None == None`.

        This differs from the standard `ne` where null values are propagated.

        Parameters
        ----------
        other
            A literal or expression value to compare with.

        See Also
        --------
        eq_missing
        ne

        Examples
        --------
        >>> s1 = pl.Series("a", [333, 200, None])
        >>> s2 = pl.Series("a", [100, 200, None])
        >>> s1.ne(s2)
        shape: (3,)
        Series: 'a' [bool]
        [
            true
            false
            null
        ]
        >>> s1.ne_missing(s2)
        shape: (3,)
        Series: 'a' [bool]
        [
            true
            false
            false
        ]
        """
        if isinstance(other, pl.Expr):
            return F.lit(self).ne_missing(other)
        return self.to_frame().select(F.col(self.name).ne_missing(other)).to_series()

    @overload
    def ge(self, other: Expr) -> Expr: ...

    @overload
    def ge(self, other: Any) -> Series: ...

    def ge(self, other: Any) -> Series | Expr:
        """Method equivalent of operator expression `series >= other`."""
        return self.__ge__(other)

    @overload
    def gt(self, other: Expr) -> Expr: ...

    @overload
    def gt(self, other: Any) -> Series: ...

    def gt(self, other: Any) -> Series | Expr:
        """Method equivalent of operator expression `series > other`."""
        return self.__gt__(other)

    def _arithmetic(self, other: Any, op_s: str, op_ffi: str) -> Self:
        if isinstance(other, pl.Expr):
            # expand pl.lit, pl.datetime, pl.duration Exprs to compatible Series
            other = self.to_frame().select_seq(other).to_series()
        elif other is None:
            other = pl.Series("", [None])

        if isinstance(other, Series):
            return self._from_pyseries(getattr(self._s, op_s)(other._s))
        elif _check_for_numpy(other) and isinstance(other, np.ndarray):
            return self._from_pyseries(getattr(self._s, op_s)(Series(other)._s))
        elif (
            isinstance(other, (float, date, datetime, timedelta, str))
            and not self.dtype.is_float()
        ):
            _s = sequence_to_pyseries(self.name, [other])
            if "rhs" in op_ffi:
                return self._from_pyseries(getattr(_s, op_s)(self._s))
            else:
                return self._from_pyseries(getattr(self._s, op_s)(_s))

        if self.dtype.is_decimal() and isinstance(other, (PyDecimal, int)):
            if isinstance(other, int):
                pyseries = sequence_to_pyseries(self.name, [other])
                _s = self._from_pyseries(pyseries).cast(Decimal(scale=0))._s
            else:
                _s = sequence_to_pyseries(self.name, [other], dtype=Decimal)

            if "rhs" in op_ffi:
                return self._from_pyseries(getattr(_s, op_s)(self._s))
            else:
                return self._from_pyseries(getattr(self._s, op_s)(_s))
        else:
            other = maybe_cast(other, self.dtype)
            f = get_ffi_func(op_ffi, self.dtype, self._s)
        if f is None:
            msg = (
                f"cannot do arithmetic with Series of dtype: {self.dtype!r} and argument"
                f" of type: {type(other).__name__!r}"
            )
            raise TypeError(msg)
        return self._from_pyseries(f(other))

    @overload
    def __add__(self, other: DataFrame) -> DataFrame: ...

    @overload
    def __add__(self, other: Expr) -> Expr: ...

    @overload
    def __add__(self, other: Any) -> Self: ...

    def __add__(self, other: Any) -> Series | DataFrame | Expr:
        if isinstance(other, str):
            other = Series("", [other])
        elif isinstance(other, pl.DataFrame):
            return other + self
        elif isinstance(other, pl.Expr):
            return F.lit(self) + other
        if self.dtype.is_decimal() and isinstance(other, (float, int)):
            return self.to_frame().select(F.col(self.name) + other).to_series()
        return self._arithmetic(other, "add", "add_<>")

    @overload
    def __sub__(self, other: Expr) -> Expr: ...

    @overload
    def __sub__(self, other: Any) -> Self: ...

    def __sub__(self, other: Any) -> Series | Expr:
        if isinstance(other, pl.Expr):
            return F.lit(self) - other
        if self.dtype.is_decimal() and isinstance(other, (float, int)):
            return self.to_frame().select(F.col(self.name) - other).to_series()
        return self._arithmetic(other, "sub", "sub_<>")

    def _recursive_cast_to_dtype(self, leaf_dtype: PolarsDataType) -> Series:
        """
        Convert leaf dtype the to given primitive datatype.

        This is equivalent to logic in DataType::cast_leaf() in Rust.
        """

        def convert_to_primitive(dtype: PolarsDataType) -> PolarsDataType:
            if isinstance(dtype, Array):
                return Array(convert_to_primitive(dtype.inner), shape=dtype.shape)
            if isinstance(dtype, List):
                return List(convert_to_primitive(dtype.inner))
            return leaf_dtype

        return self.cast(convert_to_primitive(self.dtype))

    @overload
    def __truediv__(self, other: Expr) -> Expr: ...

    @overload
    def __truediv__(self, other: Any) -> Series: ...

    def __truediv__(self, other: Any) -> Series | Expr:
        if isinstance(other, pl.Expr):
            return F.lit(self) / other
        if self.dtype.is_temporal() and not isinstance(self.dtype, Duration):
            msg = "first cast to integer before dividing datelike dtypes"
            raise TypeError(msg)
        if isinstance(other, (int, float)) and (
            self.dtype.is_decimal() or isinstance(self.dtype, Duration)
        ):
            return self.to_frame().select(F.col(self.name) / other).to_series()

        self = (
            self
            if (
                self.dtype.is_float()
                or self.dtype.is_decimal()
                or isinstance(self.dtype, (List, Array, Duration))
                or (
                    isinstance(other, Series) and isinstance(other.dtype, (List, Array))
                )
            )
            else self._recursive_cast_to_dtype(Float64())
        )

        return self._arithmetic(other, "div", "div_<>")

    @overload
    def __floordiv__(self, other: Expr) -> Expr: ...

    @overload
    def __floordiv__(self, other: Any) -> Series: ...

    def __floordiv__(self, other: Any) -> Series | Expr:
        if isinstance(other, pl.Expr):
            return F.lit(self) // other
        if self.dtype.is_temporal():
            msg = "first cast to integer before dividing datelike dtypes"
            raise TypeError(msg)
        if self.dtype.is_decimal() and isinstance(other, (float, int)):
            return self.to_frame().select(F.col(self.name) // other).to_series()

        if not isinstance(other, pl.Expr):
            other = F.lit(other)
        return self.to_frame().select_seq(F.col(self.name) // other).to_series()

    def __invert__(self) -> Series:
        return self.not_()

    @overload
    def __mul__(self, other: Expr) -> Expr: ...

    @overload
    def __mul__(self, other: DataFrame) -> DataFrame: ...

    @overload
    def __mul__(self, other: Any) -> Series: ...

    def __mul__(self, other: Any) -> Series | DataFrame | Expr:
        if isinstance(other, pl.Expr):
            return F.lit(self) * other
        if self.dtype.is_temporal() and not isinstance(self.dtype, Duration):
            msg = "first cast to integer before multiplying datelike dtypes"
            raise TypeError(msg)
        if isinstance(other, (int, float)) and (
            self.dtype.is_decimal() or isinstance(self.dtype, Duration)
        ):
            return self.to_frame().select(F.col(self.name) * other).to_series()
        elif isinstance(other, pl.DataFrame):
            return other * self
        else:
            return self._arithmetic(other, "mul", "mul_<>")

    @overload
    def __mod__(self, other: Expr) -> Expr: ...

    @overload
    def __mod__(self, other: Any) -> Series: ...

    def __mod__(self, other: Any) -> Series | Expr:
        if isinstance(other, pl.Expr):
            return F.lit(self).__mod__(other)
        if self.dtype.is_temporal():
            msg = "first cast to integer before applying modulo on datelike dtypes"
            raise TypeError(msg)
        if self.dtype.is_decimal() and isinstance(other, (float, int)):
            return self.to_frame().select(F.col(self.name) % other).to_series()
        return self._arithmetic(other, "rem", "rem_<>")

    def __rmod__(self, other: Any) -> Series:
        if self.dtype.is_temporal():
            msg = "first cast to integer before applying modulo on datelike dtypes"
            raise TypeError(msg)
        return self._arithmetic(other, "rem", "rem_<>_rhs")

    def __radd__(self, other: Any) -> Series:
        if isinstance(other, str) or (
            isinstance(other, (int, float)) and self.dtype.is_decimal()
        ):
            return self.to_frame().select(other + F.col(self.name)).to_series()
        return self._arithmetic(other, "add", "add_<>_rhs")

    def __rsub__(self, other: Any) -> Series:
        if isinstance(other, (int, float)) and self.dtype.is_decimal():
            return self.to_frame().select(other - F.col(self.name)).to_series()
        return self._arithmetic(other, "sub", "sub_<>_rhs")

    def __rtruediv__(self, other: Any) -> Series:
        if self.dtype.is_temporal():
            msg = "first cast to integer before dividing datelike dtypes"
            raise TypeError(msg)
        if self.dtype.is_float():
            self.__rfloordiv__(other)
        if isinstance(other, (int, float)) and self.dtype.is_decimal():
            return self.to_frame().select(other / F.col(self.name)).to_series()

        if isinstance(other, int):
            other = float(other)
        return self.cast(Float64).__rfloordiv__(other)

    def __rfloordiv__(self, other: Any) -> Series:
        if self.dtype.is_temporal():
            msg = "first cast to integer before dividing datelike dtypes"
            raise TypeError(msg)
        return self._arithmetic(other, "div", "div_<>_rhs")

    def __rmul__(self, other: Any) -> Series:
        if self.dtype.is_temporal() and not isinstance(self.dtype, Duration):
            msg = "first cast to integer before multiplying datelike dtypes"
            raise TypeError(msg)
        if isinstance(other, (int, float)) and (
            self.dtype.is_decimal() or isinstance(self.dtype, Duration)
        ):
            return self.to_frame().select(other * F.col(self.name)).to_series()
        return self._arithmetic(other, "mul", "mul_<>")

    def __pow__(self, exponent: int | float | Series) -> Series:
        return self.pow(exponent)

    def __rpow__(self, other: Any) -> Series:
        return (
            self.to_frame()
            .select_seq((other ** F.col(self.name)).alias(self.name))
            .to_series()
        )

    def __matmul__(self, other: Any) -> float | Series | None:
        if isinstance(other, Sequence) or (
            _check_for_numpy(other) and isinstance(other, np.ndarray)
        ):
            other = Series(other)
        # elif isinstance(other, pl.DataFrame):
        #     return other.__rmatmul__(self)  # type: ignore[return-value]
        return self.dot(other)

    def __rmatmul__(self, other: Any) -> float | Series | None:
        if isinstance(other, Sequence) or (
            _check_for_numpy(other) and isinstance(other, np.ndarray)
        ):
            other = Series(other)
        return other.dot(self)

    def __neg__(self) -> Series:
        return self.to_frame().select_seq(-F.col(self.name)).to_series()

    def __pos__(self) -> Series:
        return self

    def __abs__(self) -> Series:
        return self.abs()

    def __copy__(self) -> Self:
        return self.clone()

    def __deepcopy__(self, memo: None = None) -> Self:
        return self.clone()

    def __contains__(self, item: Any) -> bool:
        if item is None:
            return self.has_nulls()
        return self.implode().list.contains(item).item()

    def __iter__(self) -> Generator[Any]:
        if self.dtype in (List, Array):
            # TODO: either make a change and return py-native list data here, or find
            #  a faster way to return nested/List series; sequential 'get_index' calls
            #  make this path a lot slower (~10x) than it needs to be.
            get_index = self._s.get_index
            for idx in range(self.len()):
                yield get_index(idx)
        else:
            buffer_size = 25_000
            for offset in range(0, self.len(), buffer_size):
                yield from self.slice(offset, buffer_size).to_list()

    @overload
    def __getitem__(self, key: SingleIndexSelector) -> Any: ...

    @overload
    def __getitem__(self, key: MultiIndexSelector) -> Series: ...

    def __getitem__(
        self, key: SingleIndexSelector | MultiIndexSelector
    ) -> Any | Series:
        """
        Get part of the Series as a new Series or scalar.

        Parameters
        ----------
        key
            Row(s) to select.

        Returns
        -------
        Series or scalar, depending on `key`.

        Examples
        --------
        >>> s = pl.Series("a", [1, 4, 2])
        >>> s[0]
        1
        >>> s[0:2]
        shape: (2,)
        Series: 'a' [i64]
        [
            1
            4
        ]
        """
        return get_series_item_by_key(self, key)

    def __setitem__(
        self,
        key: int | Series | np.ndarray[Any, Any] | Sequence[object] | tuple[object],
        value: Any,
    ) -> None:
        # do the single idx as first branch as those are likely in a tight loop
        if isinstance(key, int) and not isinstance(key, bool):
            self.scatter(key, value)
            return None
        elif isinstance(value, Sequence) and not isinstance(value, str):
            if self.dtype.is_numeric() or self.dtype.is_temporal():
                self.scatter(key, value)  # type: ignore[arg-type]
                return None
            msg = (
                f"cannot set Series of dtype: {self.dtype!r} with list/tuple as value;"
                " use a scalar value"
            )
            raise TypeError(msg)
        if isinstance(key, Series):
            if key.dtype == Boolean:
                self._s = self.set(key, value)._s
            elif key.dtype == UInt64:
                self._s = self.scatter(key.cast(UInt32), value)._s
            elif key.dtype == UInt32:
                self._s = self.scatter(key, value)._s

        # TODO: implement for these types without casting to series
        elif _check_for_numpy(key) and isinstance(key, np.ndarray):
            if key.dtype == np.bool_:
                # boolean numpy mask
                self._s = self.scatter(np.argwhere(key)[:, 0], value)._s
            else:
                s = self._from_pyseries(
                    PySeries.new_u32("", np.array(key, np.uint32), _strict=True)
                )
                self.__setitem__(s, value)
        elif isinstance(key, (list, tuple)):
            s = self._from_pyseries(sequence_to_pyseries("", key, dtype=UInt32))
            self.__setitem__(s, value)
        else:
            msg = f'cannot use "{key!r}" for indexing'
            raise TypeError(msg)

    def __array__(
        self, dtype: npt.DTypeLike | None = None, copy: bool | None = None
    ) -> np.ndarray[Any, Any]:
        """
        Return a NumPy ndarray with the given data type.

        This method ensures a Polars Series can be treated as a NumPy ndarray.
        It enables `np.asarray` and NumPy universal functions.

        See the NumPy documentation for more information:
        https://numpy.org/doc/stable/user/basics.interoperability.html#the-array-method

        See Also
        --------
        __array_ufunc__
        """
        # Cast String types to fixed-length string to support string ufuncs
        # TODO: Use variable-length strings instead when NumPy 2.0.0 comes out:
        # https://numpy.org/devdocs/reference/routines.dtypes.html#numpy.dtypes.StringDType
        if dtype is None and not self.has_nulls() and self.dtype == String:
            dtype = np.dtype("U")

        if copy is None:
            writable, allow_copy = False, True
        elif copy is True:
            writable, allow_copy = True, True
        elif copy is False:
            writable, allow_copy = False, False
        else:
            msg = f"invalid input for `copy`: {copy!r}"
            raise TypeError(msg)

        arr = self.to_numpy(writable=writable, allow_copy=allow_copy)

        if dtype is not None and dtype != arr.dtype:
            if copy is False:
                # TODO: Only raise when data must be copied
                msg = f"copy not allowed: cast from {arr.dtype} to {dtype} prohibited"
                raise RuntimeError(msg)

            arr = arr.__array__(dtype)

        return arr

    def __array_ufunc__(
        self, ufunc: np.ufunc, method: str, *inputs: Any, **kwargs: Any
    ) -> Series:
        """Numpy universal functions."""
        if self._s.n_chunks() > 1:
            self._s.rechunk(in_place=True)

        s = self._s

        if method == "__call__":
            if ufunc.nout != 1:
                msg = "only ufuncs that return one 1D array are supported"
                raise NotImplementedError(msg)

            args: list[int | float | np.ndarray[Any, Any]] = []
            for arg in inputs:
                if isinstance(arg, (int, float, np.ndarray)):
                    args.append(arg)
                elif isinstance(arg, Series):
                    phys_arg = arg.to_physical()
                    if phys_arg._s.n_chunks() > 1:
                        phys_arg._s.rechunk(in_place=True)
                    args.append(phys_arg._s.to_numpy_view())  # type: ignore[arg-type]
                else:
                    msg = f"unsupported type {qualified_type_name(arg)!r} for {arg!r}"
                    raise TypeError(msg)

            # Get minimum dtype needed to be able to cast all input arguments to the
            # same dtype.
            dtype_char_minimum: str = np.result_type(*args).char

            # Get all possible output dtypes for ufunc.
            # Input dtypes and output dtypes seem to always match for ufunc.types,
            # so pick all the different output dtypes.
            dtypes_ufunc = [
                input_output_type[-1]
                for input_output_type in ufunc.types
                if supported_numpy_char_code(input_output_type[-1])
            ]

            # Get the first ufunc dtype from all possible ufunc dtypes for which
            # the input arguments can be safely cast to that ufunc dtype.
            for dtype_ufunc in dtypes_ufunc:
                if np.can_cast(dtype_char_minimum, dtype_ufunc):
                    dtype_char_minimum = dtype_ufunc
                    break

            # Override minimum dtype if requested.
            dtype_char = (
                np.dtype(kwargs.pop("dtype")).char
                if "dtype" in kwargs
                else dtype_char_minimum
            )

            # Only generalized ufuncs have a signature set:
            is_generalized_ufunc = bool(ufunc.signature)

            if is_generalized_ufunc:
                # Generalized ufuncs will operate on the whole array, so
                # missing data can corrupt the results.
                if self.has_nulls():
                    msg = "can't pass a Series with missing data to a generalized ufunc, as it might give unexpected results. See https://docs.pola.rs/user-guide/expressions/missing-data/ for suggestions on how to remove or fill in missing data."
                    raise ComputeError(msg)
                # If the input and output are the same size, e.g. "(n)->(n)" we
                # can allocate ourselves and save a copy. If they're different,
                # we let the ufunc do the allocation, since only it knows the
                # output size.
                assert ufunc.signature is not None  # pacify MyPy
                ufunc_input, ufunc_output = ufunc.signature.split("->")
                if ufunc_output == "()":
                    # If the result a scalar, just let the function do its
                    # thing, no need for any song and dance involving
                    # allocation:
                    return ufunc(*args, dtype=dtype_char, **kwargs)
                else:
                    allocate_output = ufunc_input == ufunc_output
            else:
                allocate_output = True

            f = get_ffi_func("apply_ufunc_<>", numpy_char_code_to_dtype(dtype_char), s)

            if f is None:
                msg = (
                    "could not find "
                    f"`apply_ufunc_{numpy_char_code_to_dtype(dtype_char)}`"
                )
                raise NotImplementedError(msg)

            series = f(
                lambda out: ufunc(*args, out=out, dtype=dtype_char, **kwargs),
                allocate_output,
            )

            result = self._from_pyseries(series)
            if is_generalized_ufunc:
                # In this case we've disallowed passing in missing data, so no
                # further processing is needed.
                return result

            # We're using a regular ufunc, that operates value by value. That
            # means we allowed missing data in the input, so filter it out:
            validity_mask = self.is_not_null()
            for arg in inputs:
                if isinstance(arg, Series):
                    validity_mask &= arg.is_not_null()
            return (
                result.to_frame()
                .select(F.when(validity_mask).then(F.col(self.name)))
                .to_series(0)
            )
        else:
            msg = (
                "only `__call__` is implemented for numpy ufuncs on a Series, got "
                f"`{method!r}`"
            )
            raise NotImplementedError(msg)

    def __arrow_c_stream__(self, requested_schema: object | None = None) -> object:
        """
        Export a Series via the Arrow PyCapsule Interface.

        https://arrow.apache.org/docs/dev/format/CDataInterface/PyCapsuleInterface.html
        """
        return self._s.__arrow_c_stream__(requested_schema)

    def _repr_html_(self) -> str:
        """Format output data in HTML for display in Jupyter Notebooks."""
        return self.to_frame()._repr_html_(_from_series=True)

    def item(self, index: int | None = None) -> Any:
        """
        Return the Series as a scalar, or return the element at the given index.

        If no index is provided, this is equivalent to `s[0]`, with a check
        that the shape is (1,). With an index, this is equivalent to `s[index]`.

        Examples
        --------
        >>> s1 = pl.Series("a", [1])
        >>> s1.item()
        1
        >>> s2 = pl.Series("a", [9, 8, 7])
        >>> s2.cum_sum().item(-1)
        24
        """
        if index is None:
            if len(self) != 1:
                msg = (
                    "can only call '.item()' if the Series is of length 1,"
                    f" or an explicit index is provided (Series is of length {len(self)})"
                )
                raise ValueError(msg)
            return self._s.get_index(0)

        return self._s.get_index_signed(index)

    def estimated_size(self, unit: SizeUnit = "b") -> int | float:
        """
        Return an estimation of the total (heap) allocated size of the Series.

        Estimated size is given in the specified unit (bytes by default).

        This estimation is the sum of the size of its buffers, validity, including
        nested arrays. Multiple arrays may share buffers and bitmaps. Therefore, the
        size of 2 arrays is not the sum of the sizes computed from this function. In
        particular, [`StructArray`]'s size is an upper bound.

        When an array is sliced, its allocated size remains constant because the buffer
        unchanged. However, this function will yield a smaller number. This is because
        this function returns the visible size of the buffer, not its total capacity.

        FFI buffers are included in this estimation.

        Notes
        -----
        For data with Object dtype, the estimated size only reports the pointer
        size, which is a huge underestimation.

        Parameters
        ----------
        unit : {'b', 'kb', 'mb', 'gb', 'tb'}
            Scale the returned size to the given unit.

        Examples
        --------
        >>> s = pl.Series("values", list(range(1_000_000)), dtype=pl.UInt32)
        >>> s.estimated_size()
        4000000
        >>> s.estimated_size("mb")
        3.814697265625
        """
        sz = self._s.estimated_size()
        return scale_bytes(sz, unit)

    def sqrt(self) -> Series:
        """
        Compute the square root of the elements.

        Syntactic sugar for

        >>> pl.Series([1, 2]) ** 0.5
        shape: (2,)
        Series: '' [f64]
        [
            1.0
            1.414214
        ]

        Examples
        --------
        >>> s = pl.Series([1, 2, 3])
        >>> s.sqrt()
        shape: (3,)
        Series: '' [f64]
        [
            1.0
            1.414214
            1.732051
        ]
        """

    def cbrt(self) -> Series:
        """
        Compute the cube root of the elements.

        Optimization for

        >>> pl.Series([1, 2]) ** (1.0 / 3)
        shape: (2,)
        Series: '' [f64]
        [
            1.0
            1.259921
        ]

        Examples
        --------
        >>> s = pl.Series([1, 2, 3])
        >>> s.cbrt()
        shape: (3,)
        Series: '' [f64]
        [
            1.0
            1.259921
            1.44225
        ]
        """

    @overload
    def any(self, *, ignore_nulls: Literal[True] = ...) -> bool: ...

    @overload
    def any(self, *, ignore_nulls: bool) -> bool | None: ...

    def any(self, *, ignore_nulls: bool = True) -> bool | None:
        """
        Return whether any of the values in the column are `True`.

        Only works on columns of data type :class:`Boolean`.

        Parameters
        ----------
        ignore_nulls
            * If set to `True` (default), null values are ignored. If there
              are no non-null values, the output is `False`.
            * If set to `False`, `Kleene logic`_ is used to deal with nulls:
              if the column contains any null values and no `True` values,
              the output is `None`.

            .. _Kleene logic: https://en.wikipedia.org/wiki/Three-valued_logic

        Returns
        -------
        bool or None

        Examples
        --------
        >>> pl.Series([True, False]).any()
        True
        >>> pl.Series([False, False]).any()
        False
        >>> pl.Series([None, False]).any()
        False

        Enable Kleene logic by setting `ignore_nulls=False`.

        >>> pl.Series([None, False]).any(ignore_nulls=False)  # Returns None
        """
        return self._s.any(ignore_nulls=ignore_nulls)

    @overload
    def all(self, *, ignore_nulls: Literal[True] = ...) -> bool: ...

    @overload
    def all(self, *, ignore_nulls: bool) -> bool | None: ...

    def all(self, *, ignore_nulls: bool = True) -> bool | None:
        """
        Return whether all values in the column are `True`.

        Only works on columns of data type :class:`Boolean`.

        Parameters
        ----------
        ignore_nulls
            * If set to `True` (default), null values are ignored. If there
              are no non-null values, the output is `True`.
            * If set to `False`, `Kleene logic`_ is used to deal with nulls:
              if the column contains any null values and no `False` values,
              the output is `None`.

            .. _Kleene logic: https://en.wikipedia.org/wiki/Three-valued_logic

        Returns
        -------
        bool or None

        Examples
        --------
        >>> pl.Series([True, True]).all()
        True
        >>> pl.Series([False, True]).all()
        False
        >>> pl.Series([None, True]).all()
        True

        Enable Kleene logic by setting `ignore_nulls=False`.

        >>> pl.Series([None, True]).all(ignore_nulls=False)  # Returns None
        """
        return self._s.all(ignore_nulls=ignore_nulls)

    def log(self, base: float | Series = math.e) -> Series:
        """
        Compute the logarithm to a given base.

        Examples
        --------
        >>> s = pl.Series([1, 2, 3])
        >>> s.log()
        shape: (3,)
        Series: '' [f64]
        [
            0.0
            0.693147
            1.098612
        ]
        """

    def log1p(self) -> Series:
        """
        Compute the natural logarithm of the input array plus one, element-wise.

        Examples
        --------
        >>> s = pl.Series([1, 2, 3])
        >>> s.log1p()
        shape: (3,)
        Series: '' [f64]
        [
            0.693147
            1.098612
            1.386294
        ]
        """

    def log10(self) -> Series:
        """
        Compute the base 10 logarithm of the input array, element-wise.

        Examples
        --------
        >>> s = pl.Series([10, 100, 1000])
        >>> s.log10()
        shape: (3,)
        Series: '' [f64]
        [
            1.0
            2.0
            3.0
        ]
        """

    def exp(self) -> Series:
        """
        Compute the exponential, element-wise.

        Examples
        --------
        >>> s = pl.Series([1, 2, 3])
        >>> s.exp()
        shape: (3,)
        Series: '' [f64]
        [
            2.718282
            7.389056
            20.085537
        ]
        """

    def drop_nulls(self) -> Series:
        """
        Drop all null values.

        The original order of the remaining elements is preserved.

        See Also
        --------
        drop_nans

        Notes
        -----
        A null value is not the same as a NaN value.
        To drop NaN values, use :func:`drop_nans`.

        Examples
        --------
        >>> s = pl.Series([1.0, None, 3.0, float("nan")])
        >>> s.drop_nulls()
        shape: (3,)
        Series: '' [f64]
        [
                1.0
                3.0
                NaN
        ]
        """

    def drop_nans(self) -> Series:
        """
        Drop all floating point NaN values.

        The original order of the remaining elements is preserved.

        See Also
        --------
        drop_nulls

        Notes
        -----
        A NaN value is not the same as a null value.
        To drop null values, use :func:`drop_nulls`.

        Examples
        --------
        >>> s = pl.Series([1.0, None, 3.0, float("nan")])
        >>> s.drop_nans()
        shape: (3,)
        Series: '' [f64]
        [
                1.0
                null
                3.0
        ]
        """

    def to_frame(self, name: str | None = None) -> DataFrame:
        """
        Cast this Series to a DataFrame.

        Parameters
        ----------
        name
            optionally name/rename the Series column in the new DataFrame.

        Examples
        --------
        >>> s = pl.Series("a", [123, 456])
        >>> df = s.to_frame()
        >>> df
        shape: (2, 1)
        ┌─────┐
        │ a   │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 123 │
        │ 456 │
        └─────┘

        >>> df = s.to_frame("xyz")
        >>> df
        shape: (2, 1)
        ┌─────┐
        │ xyz │
        │ --- │
        │ i64 │
        ╞═════╡
        │ 123 │
        │ 456 │
        └─────┘
        """
        if isinstance(name, str):
            return wrap_df(PyDataFrame([self.rename(name)._s]))
        return wrap_df(PyDataFrame([self._s]))

    def describe(
        self,
        percentiles: Sequence[float] | float | None = (0.25, 0.50, 0.75),
        interpolation: QuantileMethod = "nearest",
    ) -> DataFrame:
        """
        Quick summary statistics of a Series.

        Series with mixed datatypes will return summary statistics for the datatype of
        the first value.

        Parameters
        ----------
        percentiles
            One or more percentiles to include in the summary statistics (if the
            Series has a numeric dtype). All values must be in the range `[0, 1]`.
        interpolation : {'nearest', 'higher', 'lower', 'midpoint', 'linear', 'equiprobable'}
            Interpolation method used when calculating percentiles.

        Notes
        -----
        The median is included by default as the 50% percentile.

        Returns
        -------
        DataFrame
            Mapping with summary statistics of a Series.

        Examples
        --------
        >>> s = pl.Series([1, 2, 3, 4, 5])
        >>> s.describe()
        shape: (9, 2)
        ┌────────────┬──────────┐
        │ statistic  ┆ value    │
        │ ---        ┆ ---      │
        │ str        ┆ f64      │
        ╞════════════╪══════════╡
        │ count      ┆ 5.0      │
        │ null_count ┆ 0.0      │
        │ mean       ┆ 3.0      │
        │ std        ┆ 1.581139 │
        │ min        ┆ 1.0      │
        │ 25%        ┆ 2.0      │
        │ 50%        ┆ 3.0      │
        │ 75%        ┆ 4.0      │
        │ max        ┆ 5.0      │
        └────────────┴──────────┘

        Non-numeric data types may not have all statistics available.

        >>> s = pl.Series(["aa", "aa", None, "bb", "cc"])
        >>> s.describe()
        shape: (4, 2)
        ┌────────────┬───────┐
        │ statistic  ┆ value │
        │ ---        ┆ ---   │
        │ str        ┆ str   │
        ╞════════════╪═══════╡
        │ count      ┆ 4     │
        │ null_count ┆ 1     │
        │ min        ┆ aa    │
        │ max        ┆ cc    │
        └────────────┴───────┘
        """  # noqa: W505
        stats = self.to_frame().describe(
            percentiles=percentiles,
            interpolation=interpolation,
        )
        stats.columns = ["statistic", "value"]
        return stats.filter(F.col("value").is_not_null())

    def sum(self) -> int | float:
        """
        Reduce this Series to the sum value.

        Notes
        -----
        * Dtypes in {Int8, UInt8, Int16, UInt16} are cast to
          Int64 before summing to prevent overflow issues.
        * If there are no non-null values, then the output is `0`.
          If you would prefer empty sums to return `None`, you can
          use `s.sum() if s.count() else None` instead
          of `s.sum()`.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.sum()
        6
        """
        return self._s.sum()

    def mean(self) -> PythonLiteral | None:
        """
        Reduce this Series to the mean value.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.mean()
        2.0
        """
        return self._s.mean()

    def product(self) -> int | float:
        """
        Reduce this Series to the product value.

        Notes
        -----
        If there are no non-null values, then the output is `1`.
        If you would prefer empty products to return `None`, you can
        use `s.product() if s.count() else None` instead
        of `s.product()`.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.product()
        6
        """
        return self._s.product()

    def pow(self, exponent: int | float | Series) -> Series:
        """
        Raise to the power of the given exponent.

        If the exponent is float, the result follows the dtype of exponent.
        Otherwise, it follows dtype of base.

        Parameters
        ----------
        exponent
            The exponent. Accepts Series input.

        Examples
        --------
        Raising integers to positive integers results in integers:

        >>> s = pl.Series("foo", [1, 2, 3, 4])
        >>> s.pow(3)
        shape: (4,)
        Series: 'foo' [i64]
        [
            1
            8
            27
            64
        ]

        In order to raise integers to negative integers, you can cast either the
        base or the exponent to float:

        >>> s.pow(-3.0)
        shape: (4,)
        Series: 'foo' [f64]
        [
                1.0
                0.125
                0.037037
                0.015625
        ]
        """
        if _check_for_numpy(exponent) and isinstance(exponent, np.ndarray):
            exponent = Series(exponent)
        return self.to_frame().select_seq(F.col(self.name).pow(exponent)).to_series()

    def min(self) -> PythonLiteral | None:
        """
        Get the minimal value in this Series.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.min()
        1
        """
        return self._s.min()

    def max(self) -> PythonLiteral | None:
        """
        Get the maximum value in this Series.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.max()
        3
        """
        return self._s.max()

    def nan_max(self) -> int | float | date | datetime | timedelta | str:
        """
        Get maximum value, but propagate/poison encountered NaN values.

        This differs from numpy's `nanmax` as numpy defaults to propagating NaN values,
        whereas polars defaults to ignoring them.

        Examples
        --------
        >>> s = pl.Series("a", [1, 3, 4])
        >>> s.nan_max()
        4

        >>> s = pl.Series("a", [1.0, float("nan"), 4.0])
        >>> s.nan_max()
        nan
        """
        return self.to_frame().select_seq(F.col(self.name).nan_max()).item()

    def nan_min(self) -> int | float | date | datetime | timedelta | str:
        """
        Get minimum value, but propagate/poison encountered NaN values.

        This differs from numpy's `nanmax` as numpy defaults to propagating NaN values,
        whereas polars defaults to ignoring them.

        Examples
        --------
        >>> s = pl.Series("a", [1, 3, 4])
        >>> s.nan_min()
        1

        >>> s = pl.Series("a", [1.0, float("nan"), 4.0])
        >>> s.nan_min()
        nan
        """
        return self.to_frame().select_seq(F.col(self.name).nan_min()).item()

    def std(self, ddof: int = 1) -> float | timedelta | None:
        """
        Get the standard deviation of this Series.

        Parameters
        ----------
        ddof
            “Delta Degrees of Freedom”: the divisor used in the calculation is N - ddof,
            where N represents the number of elements.
            By default ddof is 1.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.std()
        1.0
        """
        return self._s.std(ddof)

    def var(self, ddof: int = 1) -> float | timedelta | None:
        """
        Get variance of this Series.

        Parameters
        ----------
        ddof
            “Delta Degrees of Freedom”: the divisor used in the calculation is N - ddof,
            where N represents the number of elements.
            By default ddof is 1.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.var()
        1.0
        """
        return self._s.var(ddof)

    def median(self) -> PythonLiteral | None:
        """
        Get the median of this Series.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.median()
        2.0
        """
        return self._s.median()

    def quantile(
        self, quantile: float, interpolation: QuantileMethod = "nearest"
    ) -> float | None:
        """
        Get the quantile value of this Series.

        Parameters
        ----------
        quantile
            Quantile between 0.0 and 1.0.
        interpolation : {'nearest', 'higher', 'lower', 'midpoint', 'linear', 'equiprobable'}
            Interpolation method.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.quantile(0.5)
        2.0
        """  # noqa: W505
        return self._s.quantile(quantile, interpolation)

    def to_dummies(
        self,
        *,
        separator: str = "_",
        drop_first: bool = False,
        drop_nulls: bool = False,
    ) -> DataFrame:
        """
        Get dummy/indicator variables.

        Parameters
        ----------
        separator
            Separator/delimiter used when generating column names.
        drop_first
            Remove the first category from the variable being encoded.
        drop_nulls
            If there are `None` values in the series, a `null` column is not generated

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.to_dummies()
        shape: (3, 3)
        ┌─────┬─────┬─────┐
        │ a_1 ┆ a_2 ┆ a_3 │
        │ --- ┆ --- ┆ --- │
        │ u8  ┆ u8  ┆ u8  │
        ╞═════╪═════╪═════╡
        │ 1   ┆ 0   ┆ 0   │
        │ 0   ┆ 1   ┆ 0   │
        │ 0   ┆ 0   ┆ 1   │
        └─────┴─────┴─────┘

        >>> s.to_dummies(drop_first=True)
        shape: (3, 2)
        ┌─────┬─────┐
        │ a_2 ┆ a_3 │
        │ --- ┆ --- │
        │ u8  ┆ u8  │
        ╞═════╪═════╡
        │ 0   ┆ 0   │
        │ 1   ┆ 0   │
        │ 0   ┆ 1   │
        └─────┴─────┘
        """
        return wrap_df(self._s.to_dummies(separator, drop_first, drop_nulls))

    @unstable()
    def cut(
        self,
        breaks: Sequence[float],
        *,
        labels: Sequence[str] | None = None,
        left_closed: bool = False,
        include_breaks: bool = False,
    ) -> Series:
        """
        Bin continuous values into discrete categories.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Parameters
        ----------
        breaks
            List of unique cut points.
        labels
            Names of the categories. The number of labels must be equal to the number
            of cut points plus one.
        left_closed
            Set the intervals to be left-closed instead of right-closed.
        include_breaks
            Include a column with the right endpoint of the bin each observation falls
            in. This will change the data type of the output from a
            :class:`Categorical` to a :class:`Struct`.

        Returns
        -------
        Series
            Series of data type :class:`Categorical` if `include_breaks` is set to
            `False` (default), otherwise a Series of data type :class:`Struct`.

        See Also
        --------
        qcut

        Examples
        --------
        Divide the column into three categories.

        >>> s = pl.Series("foo", [-2, -1, 0, 1, 2])
        >>> s.cut([-1, 1], labels=["a", "b", "c"])
        shape: (5,)
        Series: 'foo' [cat]
        [
                "a"
                "a"
                "b"
                "b"
                "c"
        ]

        Create a DataFrame with the breakpoint and category for each value.

        >>> cut = s.cut([-1, 1], include_breaks=True).alias("cut")
        >>> s.to_frame().with_columns(cut).unnest("cut")
        shape: (5, 3)
        ┌─────┬────────────┬────────────┐
        │ foo ┆ breakpoint ┆ category   │
        │ --- ┆ ---        ┆ ---        │
        │ i64 ┆ f64        ┆ cat        │
        ╞═════╪════════════╪════════════╡
        │ -2  ┆ -1.0       ┆ (-inf, -1] │
        │ -1  ┆ -1.0       ┆ (-inf, -1] │
        │ 0   ┆ 1.0        ┆ (-1, 1]    │
        │ 1   ┆ 1.0        ┆ (-1, 1]    │
        │ 2   ┆ inf        ┆ (1, inf]   │
        └─────┴────────────┴────────────┘
        """

    @unstable()
    def qcut(
        self,
        quantiles: Sequence[float] | int,
        *,
        labels: Sequence[str] | None = None,
        left_closed: bool = False,
        allow_duplicates: bool = False,
        include_breaks: bool = False,
    ) -> Series:
        """
        Bin continuous values into discrete categories based on their quantiles.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Parameters
        ----------
        quantiles
            Either a list of quantile probabilities between 0 and 1 or a positive
            integer determining the number of bins with uniform probability.
        labels
            Names of the categories. The number of labels must be equal to the number
            of cut points plus one.
        left_closed
            Set the intervals to be left-closed instead of right-closed.
        allow_duplicates
            If set to `True`, duplicates in the resulting quantiles are dropped,
            rather than raising a `DuplicateError`. This can happen even with unique
            probabilities, depending on the data.
        include_breaks
            Include a column with the right endpoint of the bin each observation falls
            in. This will change the data type of the output from a
            :class:`Categorical` to a :class:`Struct`.

        Returns
        -------
        Series
            Series of data type :class:`Categorical` if `include_breaks` is set to
            `False` (default), otherwise a Series of data type :class:`Struct`.

        See Also
        --------
        cut

        Examples
        --------
        Divide a column into three categories according to pre-defined quantile
        probabilities.

        >>> s = pl.Series("foo", [-2, -1, 0, 1, 2])
        >>> s.qcut([0.25, 0.75], labels=["a", "b", "c"])
        shape: (5,)
        Series: 'foo' [cat]
        [
                "a"
                "a"
                "b"
                "b"
                "c"
        ]

        Divide a column into two categories using uniform quantile probabilities.

        >>> s.qcut(2, labels=["low", "high"], left_closed=True)
        shape: (5,)
        Series: 'foo' [cat]
        [
                "low"
                "low"
                "high"
                "high"
                "high"
        ]

        Create a DataFrame with the breakpoint and category for each value.

        >>> cut = s.qcut([0.25, 0.75], include_breaks=True).alias("cut")
        >>> s.to_frame().with_columns(cut).unnest("cut")
        shape: (5, 3)
        ┌─────┬────────────┬────────────┐
        │ foo ┆ breakpoint ┆ category   │
        │ --- ┆ ---        ┆ ---        │
        │ i64 ┆ f64        ┆ cat        │
        ╞═════╪════════════╪════════════╡
        │ -2  ┆ -1.0       ┆ (-inf, -1] │
        │ -1  ┆ -1.0       ┆ (-inf, -1] │
        │ 0   ┆ 1.0        ┆ (-1, 1]    │
        │ 1   ┆ 1.0        ┆ (-1, 1]    │
        │ 2   ┆ inf        ┆ (1, inf]   │
        └─────┴────────────┴────────────┘
        """

    def rle(self) -> Series:
        """
        Compress the Series data using run-length encoding.

        Run-length encoding (RLE) encodes data by storing each *run* of identical values
        as a single value and its length.

        Returns
        -------
        Series
            Series of data type `Struct` with fields `len` of data type `UInt32`
            and `value` of the original data type.

        Examples
        --------
        >>> s = pl.Series("s", [1, 1, 2, 1, None, 1, 3, 3])
        >>> s.rle().struct.unnest()
        shape: (6, 2)
        ┌─────┬───────┐
        │ len ┆ value │
        │ --- ┆ ---   │
        │ u32 ┆ i64   │
        ╞═════╪═══════╡
        │ 2   ┆ 1     │
        │ 1   ┆ 2     │
        │ 1   ┆ 1     │
        │ 1   ┆ null  │
        │ 1   ┆ 1     │
        │ 2   ┆ 3     │
        └─────┴───────┘
        """

    def rle_id(self) -> Series:
        """
        Get a distinct integer ID for each run of identical values.

        The ID starts at 0 and increases by one each time the value of the column
        changes.

        Returns
        -------
        Series
            Series of data type `UInt32`.

        See Also
        --------
        rle

        Notes
        -----
        This functionality is especially useful for defining a new group for every time
        a column's value changes, rather than for every distinct value of that column.

        Examples
        --------
        >>> s = pl.Series("s", [1, 1, 2, 1, None, 1, 3, 3])
        >>> s.rle_id()
        shape: (8,)
        Series: 's' [u32]
        [
            0
            0
            1
            2
            3
            4
            5
            5
        ]
        """

    @unstable()
    def hist(
        self,
        bins: list[float] | None = None,
        *,
        bin_count: int | None = None,
        include_category: bool = True,
        include_breakpoint: bool = True,
    ) -> DataFrame:
        """
        Bin values into buckets and count their occurrences.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Parameters
        ----------
        bins
            Bin edges. If None given, we determine the edges based on the data.
        bin_count
            If `bins` is not provided, `bin_count` uniform bins are created that fully
            encompass the data.
        include_breakpoint
            Include a column that indicates the upper breakpoint.
        include_category
            Include a column that shows the intervals as categories.

        Returns
        -------
        DataFrame

        Examples
        --------
        >>> a = pl.Series("a", [1, 3, 8, 8, 2, 1, 3])
        >>> a.hist(bin_count=4)
        shape: (4, 3)
        ┌────────────┬─────────────┬───────┐
        │ breakpoint ┆ category    ┆ count │
        │ ---        ┆ ---         ┆ ---   │
        │ f64        ┆ cat         ┆ u32   │
        ╞════════════╪═════════════╪═══════╡
        │ 2.75       ┆ [1.0, 2.75] ┆ 3     │
        │ 4.5        ┆ (2.75, 4.5] ┆ 2     │
        │ 6.25       ┆ (4.5, 6.25] ┆ 0     │
        │ 8.0        ┆ (6.25, 8.0] ┆ 2     │
        └────────────┴─────────────┴───────┘
        """
        out = (
            self.to_frame()
            .select_seq(
                F.col(self.name).hist(
                    bins=bins,
                    bin_count=bin_count,
                    include_category=include_category,
                    include_breakpoint=include_breakpoint,
                )
            )
            .to_series()
        )
        if not include_breakpoint and not include_category:
            return out.to_frame()
        else:
            return out.struct.unnest()

    def value_counts(
        self,
        *,
        sort: bool = False,
        parallel: bool = False,
        name: str | None = None,
        normalize: bool = False,
    ) -> DataFrame:
        """
        Count the occurrences of unique values.

        Parameters
        ----------
        sort
            Sort the output by count, in descending order.
            If set to `False` (default), the order is non-deterministic.
        parallel
            Execute the computation in parallel.

            .. note::
                This option should likely *not* be enabled in a `group_by` context,
                as the computation will already be parallelized per group.
        name
            Give the resulting count column a specific name; if `normalize` is
            True this defaults to "proportion", otherwise defaults to "count".
        normalize
            If True, the count is returned as the relative frequency of unique
            values normalized to 1.0.

        Returns
        -------
        DataFrame
            Columns map the unique values to their count (or proportion).

        Examples
        --------
        >>> s = pl.Series("color", ["red", "blue", "red", "green", "blue", "blue"])
        >>> s.value_counts()  # doctest: +IGNORE_RESULT
        shape: (3, 2)
        ┌───────┬───────┐
        │ color ┆ count │
        │ ---   ┆ ---   │
        │ str   ┆ u32   │
        ╞═══════╪═══════╡
        │ red   ┆ 2     │
        │ green ┆ 1     │
        │ blue  ┆ 3     │
        └───────┴───────┘

        Sort the output by count and customize the count column name.

        >>> s.value_counts(sort=True, name="n")
        shape: (3, 2)
        ┌───────┬─────┐
        │ color ┆ n   │
        │ ---   ┆ --- │
        │ str   ┆ u32 │
        ╞═══════╪═════╡
        │ blue  ┆ 3   │
        │ red   ┆ 2   │
        │ green ┆ 1   │
        └───────┴─────┘

        Return the count as a relative frequency, normalized to 1.0:

        >>> s.value_counts(sort=True, normalize=True, name="fraction")
        shape: (3, 2)
        ┌───────┬──────────┐
        │ color ┆ fraction │
        │ ---   ┆ ---      │
        │ str   ┆ f64      │
        ╞═══════╪══════════╡
        │ blue  ┆ 0.5      │
        │ red   ┆ 0.333333 │
        │ green ┆ 0.166667 │
        └───────┴──────────┘
        """
        name = name or ("proportion" if normalize else "count")
        return pl.DataFrame._from_pydf(
            self._s.value_counts(
                sort=sort, parallel=parallel, name=name, normalize=normalize
            )
        )

    def unique_counts(self) -> Series:
        """
        Return a count of the unique values in the order of appearance.

        Examples
        --------
        >>> s = pl.Series("id", ["a", "b", "b", "c", "c", "c"])
        >>> s.unique_counts()
        shape: (3,)
        Series: 'id' [u32]
        [
            1
            2
            3
        ]
        """

    def entropy(self, base: float = math.e, *, normalize: bool = True) -> float | None:
        """
        Computes the entropy.

        Uses the formula `-sum(pk * log(pk))` where `pk` are discrete probabilities.

        Parameters
        ----------
        base
            Given base, defaults to `e`
        normalize
            Normalize pk if it doesn't sum to 1.

        Examples
        --------
        >>> a = pl.Series([0.99, 0.005, 0.005])
        >>> a.entropy(normalize=True)
        0.06293300616044681
        >>> b = pl.Series([0.65, 0.10, 0.25])
        >>> b.entropy(normalize=True)
        0.8568409950394724
        """
        return (
            self.to_frame()
            .select_seq(F.col(self.name).entropy(base, normalize=normalize))
            .to_series()
            .item()
        )

    @unstable()
    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def cumulative_eval(
        self, expr: Expr, *, min_samples: int = 1, parallel: bool = False
    ) -> Series:
        """
        Run an expression over a sliding window that increases `1` slot every iteration.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        expr
            Expression to evaluate
        min_samples
            Number of valid values there should be in the window before the expression
            is evaluated. valid values = `length - null_count`
        parallel
            Run in parallel. Don't do this in a group by or another operation that
            already has much parallelization.

        Warnings
        --------
        This can be really slow as it can have `O(n^2)` complexity. Don't use this
        for operations that visit all elements.

        Examples
        --------
        >>> s = pl.Series("values", [1, 2, 3, 4, 5])
        >>> s.cumulative_eval(pl.element().first() - pl.element().last() ** 2)
        shape: (5,)
        Series: 'values' [i64]
        [
            0
            -3
            -8
            -15
            -24
        ]
        """

    def alias(self, name: str) -> Series:
        """
        Rename the series.

        Parameters
        ----------
        name
            The new name.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.alias("b")
        shape: (3,)
        Series: 'b' [i64]
        [
                1
                2
                3
        ]
        """
        s = self.clone()
        s._s.rename(name)
        return s

    def rename(self, name: str) -> Series:
        """
        Rename this Series.

        Alias for :func:`Series.alias`.

        Parameters
        ----------
        name
            New name.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.rename("b")
        shape: (3,)
        Series: 'b' [i64]
        [
                1
                2
                3
        ]
        """
        return self.alias(name)

    def chunk_lengths(self) -> list[int]:
        """
        Get the length of each individual chunk.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s2 = pl.Series("a", [4, 5, 6])

        Concatenate Series with rechunk = True

        >>> pl.concat([s, s2], rechunk=True).chunk_lengths()
        [6]

        Concatenate Series with rechunk = False

        >>> pl.concat([s, s2], rechunk=False).chunk_lengths()
        [3, 3]
        """
        return self._s.chunk_lengths()

    def n_chunks(self) -> int:
        """
        Get the number of chunks that this Series contains.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.n_chunks()
        1
        >>> s2 = pl.Series("a", [4, 5, 6])

        Concatenate Series with rechunk = True

        >>> pl.concat([s, s2], rechunk=True).n_chunks()
        1

        Concatenate Series with rechunk = False

        >>> pl.concat([s, s2], rechunk=False).n_chunks()
        2
        """
        return self._s.n_chunks()

    def cum_max(self, *, reverse: bool = False) -> Series:
        """
        Get an array with the cumulative max computed at every element.

        Parameters
        ----------
        reverse
            reverse the operation.

        Examples
        --------
        >>> s = pl.Series("s", [3, 5, 1])
        >>> s.cum_max()
        shape: (3,)
        Series: 's' [i64]
        [
            3
            5
            5
        ]
        """

    def cum_min(self, *, reverse: bool = False) -> Series:
        """
        Get an array with the cumulative min computed at every element.

        Parameters
        ----------
        reverse
            reverse the operation.

        Examples
        --------
        >>> s = pl.Series("s", [1, 2, 3])
        >>> s.cum_min()
        shape: (3,)
        Series: 's' [i64]
        [
            1
            1
            1
        ]
        """

    def cum_prod(self, *, reverse: bool = False) -> Series:
        """
        Get an array with the cumulative product computed at every element.

        Parameters
        ----------
        reverse
            reverse the operation.

        Notes
        -----
        Dtypes in {Int8, UInt8, Int16, UInt16} are cast to
        Int64 before summing to prevent overflow issues.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.cum_prod()
        shape: (3,)
        Series: 'a' [i64]
        [
            1
            2
            6
        ]
        """

    def cum_sum(self, *, reverse: bool = False) -> Series:
        """
        Get an array with the cumulative sum computed at every element.

        Parameters
        ----------
        reverse
            reverse the operation.

        Notes
        -----
        Dtypes in {Int8, UInt8, Int16, UInt16} are cast to
        Int64 before summing to prevent overflow issues.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.cum_sum()
        shape: (3,)
        Series: 'a' [i64]
        [
            1
            3
            6
        ]
        """

    def cum_count(self, *, reverse: bool = False) -> Self:
        """
        Return the cumulative count of the non-null values in the column.

        Parameters
        ----------
        reverse
            Reverse the operation.

        Examples
        --------
        >>> s = pl.Series(["x", "k", None, "d"])
        >>> s.cum_count()
        shape: (4,)
        Series: '' [u32]
        [
                1
                2
                2
                3
        ]
        """

    def slice(self, offset: int, length: int | None = None) -> Series:
        """
        Get a slice of this Series.

        Parameters
        ----------
        offset
            Start index. Negative indexing is supported.
        length
            Length of the slice. If set to `None`, all rows starting at the offset
            will be selected.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3, 4])
        >>> s.slice(1, 2)
        shape: (2,)
        Series: 'a' [i64]
        [
                2
                3
        ]
        """
        return self._from_pyseries(self._s.slice(offset=offset, length=length))

    def append(self, other: Series) -> Self:
        """
        Append a Series to this one.

        The resulting series will consist of multiple chunks.

        Parameters
        ----------
        other
            Series to append.

        Warnings
        --------
        This method modifies the series in-place. The series is returned for
        convenience only.

        See Also
        --------
        extend

        Examples
        --------
        >>> a = pl.Series("a", [1, 2, 3])
        >>> b = pl.Series("b", [4, 5])
        >>> a.append(b)
        shape: (5,)
        Series: 'a' [i64]
        [
            1
            2
            3
            4
            5
        ]

        The resulting series will consist of multiple chunks.

        >>> a.n_chunks()
        2
        """
        require_same_type(self, other)
        self._s.append(other._s)
        return self

    def extend(self, other: Series) -> Self:
        """
        Extend the memory backed by this Series with the values from another.

        Different from `append`, which adds the chunks from `other` to the chunks of
        this series, `extend` appends the data from `other` to the underlying memory
        locations and thus may cause a reallocation (which is expensive).

        If this does `not` cause a reallocation, the resulting data structure will not
        have any extra chunks and thus will yield faster queries.

        Prefer `extend` over `append` when you want to do a query after a single
        append. For instance, during online operations where you add `n` rows
        and rerun a query.

        Prefer `append` over `extend` when you want to append many times
        before doing a query. For instance, when you read in multiple files and want
        to store them in a single `Series`. In the latter case, finish the sequence
        of `append` operations with a `rechunk`.

        Parameters
        ----------
        other
            Series to extend the series with.

        Warnings
        --------
        This method modifies the series in-place. The series is returned for
        convenience only.

        See Also
        --------
        append

        Examples
        --------
        >>> a = pl.Series("a", [1, 2, 3])
        >>> b = pl.Series("b", [4, 5])
        >>> a.extend(b)
        shape: (5,)
        Series: 'a' [i64]
        [
            1
            2
            3
            4
            5
        ]

        The resulting series will consist of a single chunk.

        >>> a.n_chunks()
        1
        """
        require_same_type(self, other)
        self._s.extend(other._s)
        return self

    def filter(self, predicate: Series | Iterable[bool]) -> Self:
        """
        Filter elements by a boolean mask.

        The original order of the remaining elements is preserved.

        Elements where the filter does not evaluate to True are discarded, including
        nulls.

        Parameters
        ----------
        predicate
            Boolean mask.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> mask = pl.Series("", [True, False, True])
        >>> s.filter(mask)
        shape: (2,)
        Series: 'a' [i64]
        [
                1
                3
        ]
        """
        if not isinstance(predicate, Series):
            predicate = Series("", predicate)
        return self._from_pyseries(self._s.filter(predicate._s))

    def head(self, n: int = 10) -> Series:
        """
        Get the first `n` elements.

        Parameters
        ----------
        n
            Number of elements to return. If a negative value is passed, return all
            elements except the last `abs(n)`.

        See Also
        --------
        tail, slice

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3, 4, 5])
        >>> s.head(3)
        shape: (3,)
        Series: 'a' [i64]
        [
                1
                2
                3
        ]

        Pass a negative value to get all rows `except` the last `abs(n)`.

        >>> s.head(-3)
        shape: (2,)
        Series: 'a' [i64]
        [
                1
                2
        ]
        """
        if n < 0:
            n = max(0, self.len() + n)
        return self._from_pyseries(self._s.head(n))

    def tail(self, n: int = 10) -> Series:
        """
        Get the last `n` elements.

        Parameters
        ----------
        n
            Number of elements to return. If a negative value is passed, return all
            elements except the first `abs(n)`.

        See Also
        --------
        head, slice

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3, 4, 5])
        >>> s.tail(3)
        shape: (3,)
        Series: 'a' [i64]
        [
                3
                4
                5
        ]

        Pass a negative value to get all rows `except` the first `abs(n)`.

        >>> s.tail(-3)
        shape: (2,)
        Series: 'a' [i64]
        [
                4
                5
        ]
        """
        if n < 0:
            n = max(0, self.len() + n)
        return self._from_pyseries(self._s.tail(n))

    def limit(self, n: int = 10) -> Series:
        """
        Get the first `n` elements.

        Alias for :func:`Series.head`.

        Parameters
        ----------
        n
            Number of elements to return. If a negative value is passed, return all
            elements except the last `abs(n)`.

        See Also
        --------
        head

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3, 4, 5])
        >>> s.limit(3)
        shape: (3,)
        Series: 'a' [i64]
        [
            1
            2
            3
        ]

        Pass a negative value to get all rows `except` the last `abs(n)`.

        >>> s.limit(-3)
        shape: (2,)
        Series: 'a' [i64]
        [
                1
                2
        ]
        """
        return self.head(n)

    def gather_every(self, n: int, offset: int = 0) -> Series:
        """
        Take every nth value in the Series and return as new Series.

        Parameters
        ----------
        n
            Gather every *n*-th row.
        offset
            Start the row index at this offset.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3, 4])
        >>> s.gather_every(2)
        shape: (2,)
        Series: 'a' [i64]
        [
            1
            3
        ]
        >>> s.gather_every(2, offset=1)
        shape: (2,)
        Series: 'a' [i64]
        [
            2
            4
        ]
        """

    def sort(
        self,
        *,
        descending: bool = False,
        nulls_last: bool = False,
        multithreaded: bool = True,
        in_place: bool = False,
    ) -> Self:
        """
        Sort this Series.

        Parameters
        ----------
        descending
            Sort in descending order.
        nulls_last
            Place null values last instead of first.
        multithreaded
            Sort using multiple threads.
        in_place
            Sort in-place.

        Examples
        --------
        >>> s = pl.Series("a", [1, 3, 4, 2])
        >>> s.sort()
        shape: (4,)
        Series: 'a' [i64]
        [
                1
                2
                3
                4
        ]
        >>> s.sort(descending=True)
        shape: (4,)
        Series: 'a' [i64]
        [
                4
                3
                2
                1
        ]
        """
        if in_place:
            self._s = self._s.sort(descending, nulls_last, multithreaded)
            return self
        else:
            return self._from_pyseries(
                self._s.sort(descending, nulls_last, multithreaded)
            )

    def top_k(self, k: int = 5) -> Series:
        r"""
        Return the `k` largest elements.

        Non-null elements are always preferred over null elements. The output is
        not guaranteed to be in any particular order, call :func:`sort` after
        this function if you wish the output to be sorted.

        This has time complexity:

        .. math:: O(n)

        Parameters
        ----------
        k
            Number of elements to return.

        See Also
        --------
        top_k_by
        bottom_k
        bottom_k_by

        Examples
        --------
        >>> s = pl.Series("a", [2, 5, 1, 4, 3])
        >>> s.top_k(3)
        shape: (3,)
        Series: 'a' [i64]
        [
            5
            4
            3
        ]
        """

    def top_k_by(
        self,
        by: IntoExpr | Iterable[IntoExpr],
        k: int = 5,
        *,
        reverse: bool | Sequence[bool] = False,
    ) -> Series:
        r"""
        Return the `k` largest elements of the `by` column.

        Non-null elements are always preferred over null elements, regardless of
        the value of `reverse`. The output is not guaranteed to be in any
        particular order, call :func:`sort` after this function if you wish the
        output to be sorted.

        This has time complexity:

        .. math:: O(n \log{n})

        Parameters
        ----------
        by
            Column used to determine the largest elements.
            Accepts expression input. Strings are parsed as column names.
        k
            Number of elements to return.
        reverse
            Consider the `k` smallest elements of the `by` column (instead of the `k`
            largest). This can be specified per column by passing a sequence of
            booleans.

        See Also
        --------
        top_k
        bottom_k
        bottom_k_by

        Examples
        --------
        >>> s = pl.Series("a", [2, 5, 1, 4, 3])
        >>> s.top_k_by("a", 3)
        shape: (3,)
        Series: 'a' [i64]
        [
            5
            4
            3
        ]
        """

    def bottom_k(self, k: int = 5) -> Series:
        r"""
        Return the `k` smallest elements.

        Non-null elements are always preferred over null elements. The output is
        not guaranteed to be in any particular order, call :func:`sort` after
        this function if you wish the output to be sorted.

        This has time complexity:

        .. math:: O(n)

        Parameters
        ----------
        k
            Number of elements to return.

        See Also
        --------
        top_k
        top_k_by
        bottom_k_by

        Examples
        --------
        >>> s = pl.Series("a", [2, 5, 1, 4, 3])
        >>> s.bottom_k(3)
        shape: (3,)
        Series: 'a' [i64]
        [
            1
            2
            3
        ]
        """

    def bottom_k_by(
        self,
        by: IntoExpr | Iterable[IntoExpr],
        k: int = 5,
        *,
        reverse: bool | Sequence[bool] = False,
    ) -> Series:
        r"""
        Return the `k` smallest elements of the `by` column.

        Non-null elements are always preferred over null elements, regardless of
        the value of `reverse`. The output is not guaranteed to be in any
        particular order, call :func:`sort` after this function if you wish the
        output to be sorted.

        This has time complexity:

        .. math:: O(n \log{n})

        Parameters
        ----------
        by
            Column used to determine the smallest elements.
            Accepts expression input. Strings are parsed as column names.
        k
            Number of elements to return.
        reverse
            Consider the `k` largest elements of the `by` column( (instead of the `k`
            smallest). This can be specified per column by passing a sequence of
            booleans.

        See Also
        --------
        top_k
        top_k_by
        bottom_k

        Examples
        --------
        >>> s = pl.Series("a", [2, 5, 1, 4, 3])
        >>> s.bottom_k_by("a", 3)
        shape: (3,)
        Series: 'a' [i64]
        [
            1
            2
            3
        ]
        """

    def arg_sort(self, *, descending: bool = False, nulls_last: bool = False) -> Series:
        """
        Get the index values that would sort this Series.

        Parameters
        ----------
        descending
            Sort in descending order.
        nulls_last
            Place null values last instead of first.

        See Also
        --------
        Series.gather: Take values by index.
        Series.rank : Get the rank of each row.

        Examples
        --------
        >>> s = pl.Series("a", [5, 3, 4, 1, 2])
        >>> s.arg_sort()
        shape: (5,)
        Series: 'a' [u32]
        [
            3
            4
            1
            2
            0
        ]
        """

    def arg_unique(self) -> Series:
        """
        Get unique index as Series.

        Returns
        -------
        Series

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 2, 3])
        >>> s.arg_unique()
        shape: (3,)
        Series: 'a' [u32]
        [
                0
                1
                3
        ]
        """

    def arg_min(self) -> int | None:
        """
        Get the index of the minimal value.

        Returns
        -------
        int

        Examples
        --------
        >>> s = pl.Series("a", [3, 2, 1])
        >>> s.arg_min()
        2
        """
        return self._s.arg_min()

    def arg_max(self) -> int | None:
        """
        Get the index of the maximal value.

        Returns
        -------
        int

        Examples
        --------
        >>> s = pl.Series("a", [3, 2, 1])
        >>> s.arg_max()
        0
        """
        return self._s.arg_max()

    @overload
    def search_sorted(
        self,
        element: NonNestedLiteral | None,
        side: SearchSortedSide = ...,
        *,
        descending: bool = ...,
    ) -> int: ...

    @overload
    def search_sorted(
        self,
        element: list[NonNestedLiteral | None] | np.ndarray[Any, Any] | Expr | Series,
        side: SearchSortedSide = ...,
        *,
        descending: bool = ...,
    ) -> Series: ...

    def search_sorted(
        self,
        element: IntoExpr | np.ndarray[Any, Any] | None,
        side: SearchSortedSide = "any",
        *,
        descending: bool = False,
    ) -> int | Series:
        """
        Find indices where elements should be inserted to maintain order.

        .. math:: a[i-1] < v <= a[i]

        Parameters
        ----------
        element
            Expression or scalar value.
        side : {'any', 'left', 'right'}
            If 'any', the index of the first suitable location found is given.
            If 'left', the index of the leftmost suitable location found is given.
            If 'right', return the rightmost suitable location found is given.
        descending
            Boolean indicating whether the values are descending or not (they
            are required to be sorted either way).

        Examples
        --------
        >>> s = pl.Series("set", [1, 2, 3, 4, 4, 5, 6, 7])
        >>> s.search_sorted(4)
        3
        >>> s.search_sorted(4, "left")
        3
        >>> s.search_sorted(4, "right")
        5
        >>> s.search_sorted([1, 4, 5])
        shape: (3,)
        Series: 'set' [u32]
        [
                0
                3
                5
        ]
        >>> s.search_sorted([1, 4, 5], "left")
        shape: (3,)
        Series: 'set' [u32]
        [
                0
                3
                5
        ]
        >>> s.search_sorted([1, 4, 5], "right")
        shape: (3,)
        Series: 'set' [u32]
        [
                1
                5
                6
        ]
        """
        df = F.select(F.lit(self).search_sorted(element, side, descending=descending))
        if isinstance(element, (list, Series, pl.Expr)):
            return df.to_series()
        elif _check_for_numpy(element) and isinstance(element, np.ndarray):
            return df.to_series()
        else:
            return df.item()

    def unique(self, *, maintain_order: bool = False) -> Series:
        """
        Get unique elements in series.

        Parameters
        ----------
        maintain_order
            Maintain order of data. This requires more work.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 2, 3])
        >>> s.unique().sort()
        shape: (3,)
        Series: 'a' [i64]
        [
            1
            2
            3
        ]
        """

    def gather(
        self, indices: int | list[int] | Expr | Series | np.ndarray[Any, Any]
    ) -> Series:
        """
        Take values by index.

        Parameters
        ----------
        indices
            Index location used for selection.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3, 4])
        >>> s.gather([1, 3])
        shape: (2,)
        Series: 'a' [i64]
        [
                2
                4
        ]
        """

    def null_count(self) -> int:
        """
        Count the null values in this Series.

        Examples
        --------
        >>> s = pl.Series([1, None, None])
        >>> s.null_count()
        2
        """
        return self._s.null_count()

    def has_nulls(self) -> bool:
        """
        Check whether the Series contains one or more null values.

        Examples
        --------
        >>> s = pl.Series([1, 2, None])
        >>> s.has_nulls()
        True
        >>> s[:2].has_nulls()
        False
        """
        return self._s.has_nulls()

    @deprecated(
        "`has_validity` is deprecated; use `has_nulls` "
        "instead to check for the presence of null values."
    )
    def has_validity(self) -> bool:
        """
        Check whether the Series contains one or more null values.

        .. deprecated:: 0.20.30
            Use the :meth:`has_nulls` method instead.
        """
        return self._s.has_nulls()

    def is_empty(self) -> bool:
        """
        Check if the Series is empty.

        Examples
        --------
        >>> s = pl.Series("a", [], dtype=pl.Float32)
        >>> s.is_empty()
        True
        """
        return self.len() == 0

    def is_sorted(self, *, descending: bool = False, nulls_last: bool = False) -> bool:
        """
        Check if the Series is sorted.

        Parameters
        ----------
        descending
            Check if the Series is sorted in descending order
        nulls_last
            Set nulls at the end of the Series in sorted check.

        Examples
        --------
        >>> s = pl.Series([1, 3, 2])
        >>> s.is_sorted()
        False

        >>> s = pl.Series([3, 2, 1])
        >>> s.is_sorted(descending=True)
        True
        """
        return self._s.is_sorted(descending, nulls_last)

    def not_(self) -> Series:
        """
        Negate a boolean Series.

        Returns
        -------
        Series
            Series of data type :class:`Boolean`.

        Examples
        --------
        >>> s = pl.Series("a", [True, False, False])
        >>> s.not_()
        shape: (3,)
        Series: 'a' [bool]
        [
            false
            true
            true
        ]
        """
        return self._from_pyseries(self._s.not_())

    def is_null(self) -> Series:
        """
        Returns a boolean Series indicating which values are null.

        Returns
        -------
        Series
            Series of data type :class:`Boolean`.

        Examples
        --------
        >>> s = pl.Series("a", [1.0, 2.0, 3.0, None])
        >>> s.is_null()
        shape: (4,)
        Series: 'a' [bool]
        [
            false
            false
            false
            true
        ]
        """

    def is_not_null(self) -> Series:
        """
        Returns a boolean Series indicating which values are not null.

        Returns
        -------
        Series
            Series of data type :class:`Boolean`.

        Examples
        --------
        >>> s = pl.Series("a", [1.0, 2.0, 3.0, None])
        >>> s.is_not_null()
        shape: (4,)
        Series: 'a' [bool]
        [
            true
            true
            true
            false
        ]
        """

    def is_finite(self) -> Series:
        """
        Returns a boolean Series indicating which values are finite.

        Returns
        -------
        Series
            Series of data type :class:`Boolean`.

        Examples
        --------
        >>> import numpy as np
        >>> s = pl.Series("a", [1.0, 2.0, np.inf])
        >>> s.is_finite()
        shape: (3,)
        Series: 'a' [bool]
        [
                true
                true
                false
        ]
        """

    def is_infinite(self) -> Series:
        """
        Returns a boolean Series indicating which values are infinite.

        Returns
        -------
        Series
            Series of data type :class:`Boolean`.

        Examples
        --------
        >>> import numpy as np
        >>> s = pl.Series("a", [1.0, 2.0, np.inf])
        >>> s.is_infinite()
        shape: (3,)
        Series: 'a' [bool]
        [
                false
                false
                true
        ]
        """

    def is_nan(self) -> Series:
        """
        Returns a boolean Series indicating which values are NaN.

        Returns
        -------
        Series
            Series of data type :class:`Boolean`.

        Examples
        --------
        >>> import numpy as np
        >>> s = pl.Series("a", [1.0, 2.0, 3.0, np.nan])
        >>> s.is_nan()
        shape: (4,)
        Series: 'a' [bool]
        [
                false
                false
                false
                true
        ]
        """

    def is_not_nan(self) -> Series:
        """
        Returns a boolean Series indicating which values are not NaN.

        Returns
        -------
        Series
            Series of data type :class:`Boolean`.

        Examples
        --------
        >>> import numpy as np
        >>> s = pl.Series("a", [1.0, 2.0, 3.0, np.nan])
        >>> s.is_not_nan()
        shape: (4,)
        Series: 'a' [bool]
        [
                true
                true
                true
                false
        ]
        """

    def is_in(
        self,
        other: Series | Collection[Any],
        *,
        nulls_equal: bool = False,
    ) -> Series:
        """
        Check if elements of this Series are in the other Series.

        Parameters
        ----------
        other
            A Series or collection to search in.
        nulls_equal : bool, default False
            If True, treat null as a distinct value. Null values will not propagate.

        Returns
        -------
        Series
            Series of data type :class:`Boolean`.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s2 = pl.Series("b", [2, 4, None])
        >>> s2.is_in(s)
        shape: (3,)
        Series: 'b' [bool]
        [
                true
                false
                null
        ]
        >>> # when nulls_equal=True, None is treated as a distinct value
        >>> s2.is_in(s, nulls_equal=True)
        shape: (3,)
        Series: 'b' [bool]
        [
                true
                false
                false
        ]

        >>> # check if some values are a member of sublists
        >>> sets = pl.Series("sets", [[1, 2, 3], [1, 2], [9, 10]])
        >>> optional_members = pl.Series("optional_members", [1, 2, 3])
        >>> print(sets)
        shape: (3,)
        Series: 'sets' [list[i64]]
        [
            [1, 2, 3]
            [1, 2]
            [9, 10]
        ]
        >>> print(optional_members)
        shape: (3,)
        Series: 'optional_members' [i64]
        [
            1
            2
            3
        ]
        >>> optional_members.is_in(sets)
        shape: (3,)
        Series: 'optional_members' [bool]
        [
            true
            true
            false
        ]
        """

    def arg_true(self) -> Series:
        """
        Get index values where Boolean Series evaluate True.

        Returns
        -------
        Series
            Series of data type :class:`UInt32`.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> (s == 2).arg_true()
        shape: (1,)
        Series: 'a' [u32]
        [
                1
        ]
        """
        return F.arg_where(self, eager=True)

    def is_unique(self) -> Series:
        """
        Get mask of all unique values.

        Returns
        -------
        Series
            Series of data type :class:`Boolean`.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 2, 3])
        >>> s.is_unique()
        shape: (4,)
        Series: 'a' [bool]
        [
                true
                false
                false
                true
        ]
        """

    def is_first_distinct(self) -> Series:
        """
        Return a boolean mask indicating the first occurrence of each distinct value.

        Returns
        -------
        Series
            Series of data type :class:`Boolean`.

        Examples
        --------
        >>> s = pl.Series([1, 1, 2, 3, 2])
        >>> s.is_first_distinct()
        shape: (5,)
        Series: '' [bool]
        [
                true
                false
                true
                true
                false
        ]
        """

    def is_last_distinct(self) -> Series:
        """
        Return a boolean mask indicating the last occurrence of each distinct value.

        Returns
        -------
        Series
            Series of data type :class:`Boolean`.

        Examples
        --------
        >>> s = pl.Series([1, 1, 2, 3, 2])
        >>> s.is_last_distinct()
        shape: (5,)
        Series: '' [bool]
        [
                false
                true
                false
                true
                true
        ]
        """

    def is_duplicated(self) -> Series:
        """
        Get mask of all duplicated values.

        Returns
        -------
        Series
            Series of data type :class:`Boolean`.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 2, 3])
        >>> s.is_duplicated()
        shape: (4,)
        Series: 'a' [bool]
        [
                false
                true
                true
                false
        ]
        """

    def explode(self) -> Series:
        """
        Explode a list Series.

        This means that every item is expanded to a new row.

        Returns
        -------
        Series
            Series with the data type of the list elements.

        See Also
        --------
        Series.list.explode : Explode a list column.

        Examples
        --------
        >>> s = pl.Series("a", [[1, 2, 3], [4, 5, 6]])
        >>> s
        shape: (2,)
        Series: 'a' [list[i64]]
        [
                [1, 2, 3]
                [4, 5, 6]
        ]
        >>> s.explode()
        shape: (6,)
        Series: 'a' [i64]
        [
                1
                2
                3
                4
                5
                6
        ]
        """

    @deprecate_renamed_parameter("strict", "check_dtypes", version="0.20.31")
    def equals(
        self,
        other: Series,
        *,
        check_dtypes: bool = False,
        check_names: bool = False,
        null_equal: bool = True,
    ) -> bool:
        """
        Check whether the Series is equal to another Series.

        .. versionchanged:: 0.20.31
            The `strict` parameter was renamed `check_dtypes`.

        Parameters
        ----------
        other
            Series to compare with.
        check_dtypes
            Require data types to match.
        check_names
            Require names to match.
        null_equal
            Consider null values as equal.

        See Also
        --------
        polars.testing.assert_series_equal

        Examples
        --------
        >>> s1 = pl.Series("a", [1, 2, 3])
        >>> s2 = pl.Series("b", [4, 5, 6])
        >>> s1.equals(s1)
        True
        >>> s1.equals(s2)
        False
        """
        require_same_type(self, other)
        return self._s.equals(
            other._s,
            check_dtypes=check_dtypes,
            check_names=check_names,
            null_equal=null_equal,
        )

    def cast(
        self,
        dtype: type[int | float | str | bool] | PolarsDataType,
        *,
        strict: bool = True,
        wrap_numerical: bool = False,
    ) -> Self:
        r"""
        Cast between data types.

        Parameters
        ----------
        dtype
            DataType to cast to.
        strict
            If True invalid casts generate exceptions instead of `null`\s.
        wrap_numerical
            If True numeric casts wrap overflowing values instead of
            marking the cast as invalid.

        Examples
        --------
        >>> s = pl.Series("a", [True, False, True])
        >>> s
        shape: (3,)
        Series: 'a' [bool]
        [
            true
            false
            true
        ]

        >>> s.cast(pl.UInt32)
        shape: (3,)
        Series: 'a' [u32]
        [
            1
            0
            1
        ]
        """
        # Do not dispatch cast as it is expensive and used in other functions.
        dtype = parse_into_dtype(dtype)
        return self._from_pyseries(self._s.cast(dtype, strict, wrap_numerical))

    def to_physical(self) -> Series:
        """
        Cast to physical representation of the logical dtype.

        - :func:`polars.datatypes.Date` -> :func:`polars.datatypes.Int32`
        - :func:`polars.datatypes.Datetime` -> :func:`polars.datatypes.Int64`
        - :func:`polars.datatypes.Time` -> :func:`polars.datatypes.Int64`
        - :func:`polars.datatypes.Duration` -> :func:`polars.datatypes.Int64`
        - :func:`polars.datatypes.Categorical` -> :func:`polars.datatypes.UInt32`
        - `List(inner)` -> `List(physical of inner)`
        - `Array(inner)` -> `Array(physical of inner)`
        - `Struct(fields)` -> `Struct(physical of fields)`
        - Other data types will be left unchanged.

        Warnings
        --------
        The physical representations are an implementation detail
        and not guaranteed to be stable.

        Examples
        --------
        Replicating the pandas
        `pd.Series.factorize
        <https://pandas.pydata.org/docs/reference/api/pandas.Series.factorize.html>`_
        method.

        >>> s = pl.Series("values", ["a", None, "x", "a"])
        >>> s.cast(pl.Categorical).to_physical()
        shape: (4,)
        Series: 'values' [u32]
        [
            0
            null
            1
            0
        ]
        """

    def to_list(self) -> list[Any]:
        """
        Convert this Series to a Python list.

        This operation copies data.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.to_list()
        [1, 2, 3]
        >>> type(s.to_list())
        <class 'list'>
        """
        return self._s.to_list()

    def rechunk(self, *, in_place: bool = False) -> Self:
        """
        Create a single chunk of memory for this Series.

        Parameters
        ----------
        in_place
            In place or not.

        Examples
        --------
        >>> s1 = pl.Series("a", [1, 2, 3])
        >>> s1.n_chunks()
        1
        >>> s2 = pl.Series("a", [4, 5, 6])
        >>> s = pl.concat([s1, s2], rechunk=False)
        >>> s.n_chunks()
        2
        >>> s.rechunk(in_place=True)
        shape: (6,)
        Series: 'a' [i64]
        [
                1
                2
                3
                4
                5
                6
        ]
        >>> s.n_chunks()
        1
        """
        opt_s = self._s.rechunk(in_place)
        if in_place:
            return self
        else:
            assert opt_s is not None
            return self._from_pyseries(opt_s)

    def reverse(self) -> Series:
        """
        Return Series in reverse order.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3], dtype=pl.Int8)
        >>> s.reverse()
        shape: (3,)
        Series: 'a' [i8]
        [
            3
            2
            1
        ]
        """

    def is_between(
        self,
        lower_bound: IntoExpr,
        upper_bound: IntoExpr,
        closed: ClosedInterval = "both",
    ) -> Series:
        """
        Get a boolean mask of the values that are between the given lower/upper bounds.

        Parameters
        ----------
        lower_bound
            Lower bound value. Accepts expression input. Non-expression inputs
            (including strings) are parsed as literals.
        upper_bound
            Upper bound value. Accepts expression input. Non-expression inputs
            (including strings) are parsed as literals.
        closed : {'both', 'left', 'right', 'none'}
            Define which sides of the interval are closed (inclusive).

        Notes
        -----
        If the value of the `lower_bound` is greater than that of the `upper_bound`
        then the result will be False, as no value can satisfy the condition.

        Examples
        --------
        >>> s = pl.Series("num", [1, 2, 3, 4, 5])
        >>> s.is_between(2, 4)
        shape: (5,)
        Series: 'num' [bool]
        [
            false
            true
            true
            true
            false
        ]

        Use the `closed` argument to include or exclude the values at the bounds:

        >>> s.is_between(2, 4, closed="left")
        shape: (5,)
        Series: 'num' [bool]
        [
            false
            true
            true
            false
            false
        ]

        You can also use strings as well as numeric/temporal values:

        >>> s = pl.Series("s", ["a", "b", "c", "d", "e"])
        >>> s.is_between("b", "d", closed="both")
        shape: (5,)
        Series: 's' [bool]
        [
            false
            true
            true
            true
            false
        ]
        """
        if closed == "none":
            out = (self > lower_bound) & (self < upper_bound)
        elif closed == "both":
            out = (self >= lower_bound) & (self <= upper_bound)
        elif closed == "right":
            out = (self > lower_bound) & (self <= upper_bound)
        elif closed == "left":
            out = (self >= lower_bound) & (self < upper_bound)

        if isinstance(out, pl.Expr):
            out = F.select(out).to_series()

        return out

    def is_close(
        self,
        other: IntoExpr,
        *,
        abs_tol: float = 0.0,
        rel_tol: float = 1e-09,
        nans_equal: bool = False,
    ) -> Series:
        r"""
        Get a boolean mask of the values being close to the other values.

        Two values `a` and `b` are considered close if the following condition holds:

        .. math::
            |a-b| \le max \{ \text{rel_tol} \cdot max \{ |a|, |b| \}, \text{abs_tol} \}

        Parameters
        ----------
        other
            A literal or expression value to compare with.
        abs_tol
            Absolute tolerance. This is the maximum allowed absolute difference between
            two values. Must be non-negative.
        rel_tol
            Relative tolerance. This is the maximum allowed difference between two
            values, relative to the larger absolute value. Must be non-negative.
        nans_equal
            Whether NaN values should be considered equal.

        Returns
        -------
        Series
            Series of data type :class:`Boolean`.

        Notes
        -----
            The implementation of this method is symmetric and mirrors the behavior of
            :meth:`math.isclose`. Specifically note that this behavior is different to
            :meth:`numpy.isclose`.

        Examples
        --------
        >>> s = pl.Series("s", [1.0, 1.2, 1.4, 1.45, 1.6])
        >>> s.is_close(1.4, abs_tol=0.1)
        shape: (5,)
        Series: 's' [bool]
        [
            false
            false
            true
            true
            false
        ]
        """
        return F.select(
            F.lit(self).is_close(
                other, abs_tol=abs_tol, rel_tol=rel_tol, nans_equal=nans_equal
            )
        ).to_series()

    def to_numpy(
        self,
        *,
        writable: bool = False,
        allow_copy: bool = True,
        use_pyarrow: bool | None = None,
        zero_copy_only: bool | None = None,
    ) -> np.ndarray[Any, Any]:
        """
        Convert this Series to a NumPy ndarray.

        This operation copies data only when necessary. The conversion is zero copy when
        all of the following hold:

        - The data type is an integer, float, `Datetime`, `Duration`, or `Array`.
        - The Series contains no null values.
        - The Series consists of a single chunk.
        - The `writable` parameter is set to `False` (default).

        Parameters
        ----------
        writable
            Ensure the resulting array is writable. This will force a copy of the data
            if the array was created without copy as the underlying Arrow data is
            immutable.
        allow_copy
            Allow memory to be copied to perform the conversion. If set to `False`,
            causes conversions that are not zero-copy to fail.

        use_pyarrow
            First convert to PyArrow, then call `pyarrow.Array.to_numpy
            <https://arrow.apache.org/docs/python/generated/pyarrow.Array.html#pyarrow.Array.to_numpy>`_
            to convert to NumPy. If set to `False`, Polars' own conversion logic is
            used.

            .. deprecated:: 0.20.28
                Polars now uses its native engine by default for conversion to NumPy.
                To use PyArrow's engine, call `.to_arrow().to_numpy()` instead.

        zero_copy_only
            Raise an exception if the conversion to a NumPy would require copying
            the underlying data. Data copy occurs, for example, when the Series contains
            nulls or non-numeric types.

            .. deprecated:: 0.20.10
                Use the `allow_copy` parameter instead, which is the inverse of this
                one.

        Examples
        --------
        Numeric data without nulls can be converted without copying data.
        The resulting array will not be writable.

        >>> s = pl.Series([1, 2, 3], dtype=pl.Int8)
        >>> arr = s.to_numpy()
        >>> arr
        array([1, 2, 3], dtype=int8)
        >>> arr.flags.writeable
        False

        Set `writable=True` to force data copy to make the array writable.

        >>> s.to_numpy(writable=True).flags.writeable
        True

        Integer Series containing nulls will be cast to a float type with `nan`
        representing a null value. This requires data to be copied.

        >>> s = pl.Series([1, 2, None], dtype=pl.UInt16)
        >>> s.to_numpy()
        array([ 1.,  2., nan], dtype=float32)

        Set `allow_copy=False` to raise an error if data would be copied.

        >>> s.to_numpy(allow_copy=False)  # doctest: +SKIP
        Traceback (most recent call last):
        ...
        RuntimeError: copy not allowed: cannot convert to a NumPy array without copying data

        Series of data type `Array` and `Struct` will result in an array with more than
        one dimension.

        >>> s = pl.Series([[1, 2, 3], [4, 5, 6]], dtype=pl.Array(pl.Int64, 3))
        >>> s.to_numpy()
        array([[1, 2, 3],
               [4, 5, 6]])
        """  # noqa: W505
        if zero_copy_only is not None:
            issue_deprecation_warning(
                "the `zero_copy_only` parameter for `Series.to_numpy` is deprecated."
                " Use the `allow_copy` parameter instead, which is the inverse of `zero_copy_only`.",
                version="0.20.10",
            )
            allow_copy = not zero_copy_only

        if use_pyarrow is not None:
            issue_deprecation_warning(
                "the `use_pyarrow` parameter for `Series.to_numpy` is deprecated."
                " Polars now uses its native engine for conversion to NumPy by default."
                " To use PyArrow's engine, call `.to_arrow().to_numpy()` instead.",
                version="0.20.28",
            )
        else:
            use_pyarrow = False

        if (
            use_pyarrow
            and _PYARROW_AVAILABLE
            and self.dtype not in (Date, Datetime, Duration, Array, Object)
        ):
            if not allow_copy and self.n_chunks() > 1 and not self.is_empty():
                msg = "cannot return a zero-copy array"
                raise ValueError(msg)

            return self.to_arrow().to_numpy(
                zero_copy_only=not allow_copy, writable=writable
            )

        return self._s.to_numpy(writable=writable, allow_copy=allow_copy)

    @unstable()
    def to_jax(self, device: jax.Device | str | None = None) -> jax.Array:
        """
        Convert this Series to a Jax Array.

        .. versionadded:: 0.20.27

        .. warning::
            This functionality is currently considered **unstable**. It may be
            changed at any point without it being considered a breaking change.

        Parameters
        ----------
        device
            Specify the jax `Device` on which the array will be created; can provide
            a string (such as "cpu", "gpu", or "tpu") in which case the device is
            retrieved as `jax.devices(string)[0]`. For more specific control you
            can supply the instantiated `Device` directly. If None, arrays are
            created on the default device.

        Examples
        --------
        >>> s = pl.Series("x", [10.5, 0.0, -10.0, 5.5])
        >>> s.to_jax()
        Array([ 10.5,   0. , -10. ,   5.5], dtype=float32)
        """
        jx = import_optional(
            "jax",
            install_message="Please see `https://jax.readthedocs.io/en/latest/installation.html` "
            "for specific installation recommendations for the Jax package",
        )
        if isinstance(device, str):
            device = jx.devices(device)[0]
        if (
            jx.config.jax_enable_x64
            or bool(int(os.environ.get("JAX_ENABLE_X64", "0")))
            or self.dtype not in {Float64, Int64, UInt64}
        ):
            srs = self
        else:
            single_precision = {Float64: Float32, Int64: Int32, UInt64: UInt32}
            srs = self.cast(single_precision[self.dtype])  # type: ignore[index]

        with nullcontext() if device is None else jx.default_device(device):
            return jx.numpy.asarray(
                # note: jax arrays are immutable, so can avoid a copy (vs torch)
                a=srs.to_numpy(writable=False),
                order="K",
            )

    @unstable()
    def to_torch(self) -> torch.Tensor:
        """
        Convert this Series to a PyTorch Tensor.

        .. versionadded:: 0.20.23

        .. warning::
            This functionality is currently considered **unstable**. It may be
            changed at any point without it being considered a breaking change.

        Notes
        -----
        PyTorch tensors do not support UInt16, UInt32, or UInt64; these dtypes
        will be automatically cast to Int32, Int64, and Int64, respectively.

        Examples
        --------
        >>> s = pl.Series("x", [1, 0, 1, 2, 0], dtype=pl.UInt8)
        >>> s.to_torch()
        tensor([1, 0, 1, 2, 0], dtype=torch.uint8)
        >>> s = pl.Series("x", [5.5, -10.0, 2.5], dtype=pl.Float32)
        >>> s.to_torch()
        tensor([  5.5000, -10.0000,   2.5000])
        """
        torch = import_optional("torch")

        # PyTorch tensors do not support uint16/32/64
        if self.dtype in (UInt32, UInt64):
            srs = self.cast(Int64)
        elif self.dtype == UInt16:
            srs = self.cast(Int32)
        else:
            srs = self

        # we have to build the tensor from a writable array or PyTorch will complain
        # about it (writing to a readonly array results in undefined behavior)
        numpy_array = srs.to_numpy(writable=True)
        try:
            tensor = torch.from_numpy(numpy_array)
        except TypeError:
            if self.dtype == List:
                msg = "cannot convert List dtype to Tensor (use Array dtype instead)"
                raise TypeError(msg) from None
            raise
        # note: named tensors are currently experimental
        # tensor.rename(self.name)
        return tensor

    @deprecate_renamed_parameter("future", "compat_level", version="1.1")
    def to_arrow(self, *, compat_level: CompatLevel | None = None) -> pa.Array:
        """
        Return the underlying Arrow array.

        If the Series contains only a single chunk this operation is zero copy.

        .. versionchanged:: 1.24
            The `future` parameter was renamed `compat_level`.

        Parameters
        ----------
        compat_level
            Use a specific compatibility level
            when exporting Polars' internal data structures.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s = s.to_arrow()
        >>> s
        <pyarrow.lib.Int64Array object at ...>
        [
          1,
          2,
          3
        ]
        """
        compat_level_py: int | bool
        if compat_level is None:
            compat_level_py = False
        elif isinstance(compat_level, CompatLevel):
            compat_level_py = compat_level._version
        else:
            msg = f"`compat_level` has invalid type: {qualified_type_name(compat_level)!r}"
            raise TypeError(msg)
        return self._s.to_arrow(compat_level_py)

    def to_pandas(
        self, *, use_pyarrow_extension_array: bool = False, **kwargs: Any
    ) -> pd.Series[Any]:
        """
        Convert this Series to a pandas Series.

        This operation copies data if `use_pyarrow_extension_array` is not enabled.

        Parameters
        ----------
        use_pyarrow_extension_array
            Use a PyArrow-backed extension array instead of a NumPy array for the pandas
            Series. This allows zero copy operations and preservation of null values.
            Subsequent operations on the resulting pandas Series may trigger conversion
            to NumPy if those operations are not supported by PyArrow compute functions.
        **kwargs
            Additional keyword arguments to be passed to
            :meth:`pyarrow.Array.to_pandas`.

        Returns
        -------
        :class:`pandas.Series`

        Notes
        -----
        This operation requires that both :mod:`pandas` and :mod:`pyarrow` are
        installed.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.to_pandas()
        0    1
        1    2
        2    3
        Name: a, dtype: int64

        Null values are converted to `NaN`.

        >>> s = pl.Series("b", [1, 2, None])
        >>> s.to_pandas()
        0    1.0
        1    2.0
        2    NaN
        Name: b, dtype: float64

        Pass `use_pyarrow_extension_array=True` to get a pandas Series backed by a
        PyArrow extension array. This will preserve null values.

        >>> s.to_pandas(use_pyarrow_extension_array=True)
        0       1
        1       2
        2    <NA>
        Name: b, dtype: int64[pyarrow]
        """
        if self.dtype == Object:
            # Can't convert via PyArrow, so do it via NumPy
            return pd.Series(self.to_numpy(), dtype=object, name=self.name)

        if use_pyarrow_extension_array:
            if parse_version(pd.__version__) < (1, 5):
                msg = f'pandas>=1.5.0 is required for `to_pandas("use_pyarrow_extension_array=True")`, found Pandas {pd.__version__}'
                raise ModuleUpgradeRequiredError(msg)
            if not _PYARROW_AVAILABLE or parse_version(pa.__version__) < (8, 0):
                raise ModuleUpgradeRequiredError(
                    f'pyarrow>=8.0.0 is required for `to_pandas("use_pyarrow_extension_array=True")`'
                    f", found pyarrow {pa.__version__!r}"
                    if _PYARROW_AVAILABLE
                    else ""
                )

        pa_arr = self.to_arrow()
        # pandas does not support unsigned dictionary indices
        if pa.types.is_dictionary(pa_arr.type):
            pa_arr = pa_arr.cast(pa.dictionary(pa.int64(), pa.large_string()))

        if use_pyarrow_extension_array:
            pd_series = pa_arr.to_pandas(
                self_destruct=True,
                split_blocks=True,
                types_mapper=lambda pa_dtype: pd.ArrowDtype(pa_dtype),
                **kwargs,
            )
        else:
            date_as_object = kwargs.pop("date_as_object", False)
            pd_series = pa_arr.to_pandas(date_as_object=date_as_object, **kwargs)

        pd_series.name = self.name
        return pd_series

    def to_init_repr(self, n: int = 1000) -> str:
        """
        Convert Series to instantiable string representation.

        Parameters
        ----------
        n
            Only use first n elements.

        See Also
        --------
        polars.Series.to_init_repr
        polars.from_repr

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, None, 4], dtype=pl.Int16)
        >>> print(s.to_init_repr())
        pl.Series('a', [1, 2, None, 4], dtype=pl.Int16)
        >>> s_from_str_repr = eval(s.to_init_repr())
        >>> s_from_str_repr
        shape: (4,)
        Series: 'a' [i16]
        [
            1
            2
            null
            4
        ]
        """
        values = self.head(n).to_list()
        dtype_init_repr = dtype_to_init_repr(self.dtype)
        return f"pl.Series({self.name!r}, {values}, dtype={dtype_init_repr})"

    def count(self) -> int:
        """
        Return the number of non-null elements in the column.

        See Also
        --------
        len

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, None])
        >>> s.count()
        2
        """
        return self.len() - self.null_count()

    def len(self) -> int:
        """
        Return the number of elements in the Series.

        Null values count towards the total.

        See Also
        --------
        count

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, None])
        >>> s.len()
        3
        """
        return self._s.len()

    def set(self, filter: Series, value: int | float | str | bool | None) -> Series:
        """
        Set masked values.

        Parameters
        ----------
        filter
            Boolean mask.
        value
            Value with which to replace the masked values.

        Notes
        -----
        Use of this function is frequently an anti-pattern, as it can
        block optimisation (predicate pushdown, etc). Consider using
        `pl.when(predicate).then(value).otherwise(self)` instead.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.set(s == 2, 10)
        shape: (3,)
        Series: 'a' [i64]
        [
                1
                10
                3
        ]

        It is better to implement this as follows:

        >>> s.to_frame().select(
        ...     pl.when(pl.col("a") == 2).then(10).otherwise(pl.col("a"))
        ... )
        shape: (3, 1)
        ┌─────────┐
        │ literal │
        │ ---     │
        │ i64     │
        ╞═════════╡
        │ 1       │
        │ 10      │
        │ 3       │
        └─────────┘
        """
        f = get_ffi_func("set_with_mask_<>", self.dtype, self._s)
        if f is None:
            msg = f"Series of type {self.dtype} can not be set"
            raise NotImplementedError(msg)
        return self._from_pyseries(f(filter._s, value))

    def scatter(
        self,
        indices: Series | Iterable[int] | int | np.ndarray[Any, Any],
        values: Series | Iterable[PythonLiteral] | PythonLiteral | None,
    ) -> Series:
        """
        Set values at the index locations.

        Parameters
        ----------
        indices
            Integers representing the index locations.
        values
            Replacement values.

        Notes
        -----
        Use of this function is frequently an anti-pattern, as it can
        block optimization (predicate pushdown, etc). Consider using
        `pl.when(predicate).then(value).otherwise(self)` instead.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.scatter(1, 10)
        shape: (3,)
        Series: 'a' [i64]
        [
                1
                10
                3
        ]

        It is better to implement this as follows:

        >>> s.to_frame().with_row_index().select(
        ...     pl.when(pl.col("index") == 1).then(10).otherwise(pl.col("a"))
        ... )
        shape: (3, 1)
        ┌─────────┐
        │ literal │
        │ ---     │
        │ i64     │
        ╞═════════╡
        │ 1       │
        │ 10      │
        │ 3       │
        └─────────┘
        """
        if not isinstance(indices, Iterable):
            index: Any = indices  # Workaround for older NumPy versions
            indices = [index]
        indices = Series(values=indices)
        if indices.is_empty():
            return self

        if not isinstance(values, Series):
            if not isinstance(values, Iterable) or isinstance(values, str):
                values = [values]
            values = Series(values=values)

        self._s.scatter(indices._s, values._s)
        return self

    def index_of(self, element: IntoExpr) -> int | None:
        """
        Get the index of the first occurrence of a value, or ``None`` if it's not found.

        Parameters
        ----------
        element
            Value to find.

        Examples
        --------
        >>> s = pl.Series("a", [1, None, 17])
        >>> s.index_of(17)
        2
        >>> s.index_of(None)  # search for a null
        1
        >>> s.index_of(55) is None
        True
        """
        return F.select(F.lit(self).index_of(element)).item()

    def clear(self, n: int = 0) -> Series:
        """
        Create an empty copy of the current Series, with zero to 'n' elements.

        The copy has an identical name/dtype, but no data.

        Parameters
        ----------
        n
            Number of (empty) elements to return in the cleared frame.

        See Also
        --------
        clone : Cheap deepcopy/clone.

        Examples
        --------
        >>> s = pl.Series("a", [None, True, False])
        >>> s.clear()
        shape: (0,)
        Series: 'a' [bool]
        [
        ]

        >>> s.clear(n=2)
        shape: (2,)
        Series: 'a' [bool]
        [
            null
            null
        ]
        """
        if n < 0:
            msg = f"`n` should be greater than or equal to 0, got {n}"
            raise ValueError(msg)
        # faster path
        if n == 0:
            return self._from_pyseries(self._s.clear())
        s = (
            self.__class__(name=self.name, values=[], dtype=self.dtype)
            if len(self) > 0
            else self.clone()
        )
        return s.extend_constant(None, n=n) if n > 0 else s

    def clone(self) -> Self:
        """
        Create a copy of this Series.

        This is a cheap operation that does not copy data.

        See Also
        --------
        clear : Create an empty copy of the current Series, with identical
            schema but no data.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.clone()
        shape: (3,)
        Series: 'a' [i64]
        [
                1
                2
                3
        ]
        """
        return self._from_pyseries(self._s.clone())

    def fill_nan(self, value: int | float | Expr | None) -> Series:
        """
        Fill floating point NaN value with a fill value.

        Parameters
        ----------
        value
            Value used to fill NaN values.

        See Also
        --------
        fill_null

        Notes
        -----
        A NaN value is not the same as a null value.
        To fill null values, use :func:`fill_null`.

        Examples
        --------
        >>> s = pl.Series("a", [1.0, 2.0, 3.0, float("nan")])
        >>> s.fill_nan(0)
        shape: (4,)
        Series: 'a' [f64]
        [
                1.0
                2.0
                3.0
                0.0
        ]
        """

    def fill_null(
        self,
        value: Any | Expr | None = None,
        strategy: FillNullStrategy | None = None,
        limit: int | None = None,
    ) -> Series:
        """
        Fill null values using the specified value or strategy.

        Parameters
        ----------
        value
            Value used to fill null values.
        strategy : {None, 'forward', 'backward', 'min', 'max', 'mean', 'zero', 'one'}
            Strategy used to fill null values.
        limit
            Number of consecutive null values to fill when using the 'forward' or
            'backward' strategy.

        See Also
        --------
        backward_fill
        fill_nan
        forward_fill

        Notes
        -----
        A null value is not the same as a NaN value.
        To fill NaN values, use :func:`fill_nan`.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3, None])
        >>> s.fill_null(strategy="forward")
        shape: (4,)
        Series: 'a' [i64]
        [
            1
            2
            3
            3
        ]
        >>> s.fill_null(strategy="min")
        shape: (4,)
        Series: 'a' [i64]
        [
            1
            2
            3
            1
        ]
        >>> s = pl.Series("b", ["x", None, "z"])
        >>> s.fill_null(pl.lit(""))
        shape: (3,)
        Series: 'b' [str]
        [
            "x"
            ""
            "z"
        ]
        """

    def backward_fill(self, limit: int | None = None) -> Series:
        """
        Fill missing values with the next non-null value.

        This is an alias of `.fill_null(strategy="backward")`.

        Parameters
        ----------
        limit
            The number of consecutive null values to backward fill.

        See Also
        --------
        fill_null
        forward_fill
        shift
        """
        return self.fill_null(strategy="backward", limit=limit)

    def forward_fill(self, limit: int | None = None) -> Series:
        """
        Fill missing values with the last non-null value.

        This is an alias of `.fill_null(strategy="forward")`.

        Parameters
        ----------
        limit
            The number of consecutive null values to forward fill.

        See Also
        --------
        backward_fill
        fill_null
        shift
        """
        return self.fill_null(strategy="forward", limit=limit)

    def floor(self) -> Series:
        """
        Rounds down to the nearest integer value.

        Only works on floating point Series.

        Examples
        --------
        >>> s = pl.Series("a", [1.12345, 2.56789, 3.901234])
        >>> s.floor()
        shape: (3,)
        Series: 'a' [f64]
        [
                1.0
                2.0
                3.0
        ]
        """

    def ceil(self) -> Series:
        """
        Rounds up to the nearest integer value.

        Only works on floating point Series.

        Examples
        --------
        >>> s = pl.Series("a", [1.12345, 2.56789, 3.901234])
        >>> s.ceil()
        shape: (3,)
        Series: 'a' [f64]
        [
                2.0
                3.0
                4.0
        ]
        """

    def round(self, decimals: int = 0, mode: RoundMode = "half_to_even") -> Series:
        """
        Round underlying floating point data by `decimals` digits.

        The default rounding mode is "half to even" (also known as "bankers' rounding").

        Parameters
        ----------
        decimals
            Number of decimals to round by.
        mode : {'half_to_even', 'half_away_from_zero'}
            Rounding mode.

        Examples
        --------
        >>> s = pl.Series("a", [1.12345, 2.56789, 3.901234])
        >>> s.round(2)
        shape: (3,)
        Series: 'a' [f64]
        [
                1.12
                2.57
                3.9
        ]

        >>> s = pl.Series([-3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5])
        >>> s.round(mode="half_to_even")
        shape: (8,)
        Series: '' [f64]
        [
            -4.0
            -2.0
            -2.0
            -0.0
            0.0
            2.0
            2.0
            4.0
        ]
        """

    def round_sig_figs(self, digits: int) -> Series:
        """
        Round to a number of significant figures.

        Parameters
        ----------
        digits
            Number of significant figures to round to.

        Examples
        --------
        >>> s = pl.Series([0.01234, 3.333, 3450.0])
        >>> s.round_sig_figs(2)
        shape: (3,)
        Series: '' [f64]
        [
                0.012
                3.3
                3500.0
        ]
        """

    def dot(self, other: Series | ArrayLike) -> int | float | None:
        """
        Compute the dot/inner product between two Series.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s2 = pl.Series("b", [4.0, 5.0, 6.0])
        >>> s.dot(s2)
        32.0

        Parameters
        ----------
        other
            Series (or array) to compute dot product with.
        """
        if not isinstance(other, Series):
            other = Series(other)
        if len(self) != len(other):
            n, m = len(self), len(other)
            msg = f"Series length mismatch: expected {n!r}, found {m!r}"
            raise ShapeError(msg)
        return self._s.dot(other._s)

    def mode(self) -> Series:
        """
        Compute the most occurring value(s).

        Can return multiple Values.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 2, 3])
        >>> s.mode()
        shape: (1,)
        Series: 'a' [i64]
        [
                2
        ]
        """

    def sign(self) -> Series:
        """
        Compute the element-wise sign function on numeric types.

        The returned value is computed as follows:

        * -1 if x < 0.
        *  1 if x > 0.
        *  x otherwise (typically 0, but could be NaN if the input is).

        Null values are preserved as-is, and the dtype of the input is preserved.

        Examples
        --------
        >>> s = pl.Series("a", [-9.0, -0.0, 0.0, 4.0, float("nan"), None])
        >>> s.sign()
        shape: (6,)
        Series: 'a' [f64]
        [
                -1.0
                -0.0
                0.0
                1.0
                NaN
                null
        ]
        """

    def sin(self) -> Series:
        """
        Compute the element-wise value for the sine.

        Examples
        --------
        >>> import math
        >>> s = pl.Series("a", [0.0, math.pi / 2.0, math.pi])
        >>> s.sin()
        shape: (3,)
        Series: 'a' [f64]
        [
            0.0
            1.0
            1.2246e-16
        ]
        """

    def cos(self) -> Series:
        """
        Compute the element-wise value for the cosine.

        Examples
        --------
        >>> import math
        >>> s = pl.Series("a", [0.0, math.pi / 2.0, math.pi])
        >>> s.cos()
        shape: (3,)
        Series: 'a' [f64]
        [
            1.0
            6.1232e-17
            -1.0
        ]
        """

    def tan(self) -> Series:
        """
        Compute the element-wise value for the tangent.

        Examples
        --------
        >>> import math
        >>> s = pl.Series("a", [0.0, math.pi / 2.0, math.pi])
        >>> s.tan()
        shape: (3,)
        Series: 'a' [f64]
        [
            0.0
            1.6331e16
            -1.2246e-16
        ]
        """

    def cot(self) -> Series:
        """
        Compute the element-wise value for the cotangent.

        Examples
        --------
        >>> import math
        >>> s = pl.Series("a", [0.0, math.pi / 2.0, math.pi])
        >>> s.cot()
        shape: (3,)
        Series: 'a' [f64]
        [
            inf
            6.1232e-17
            -8.1656e15
        ]
        """

    def arcsin(self) -> Series:
        """
        Compute the element-wise value for the inverse sine.

        Examples
        --------
        >>> s = pl.Series("a", [1.0, 0.0, -1.0])
        >>> s.arcsin()
        shape: (3,)
        Series: 'a' [f64]
        [
            1.570796
            0.0
            -1.570796
        ]
        """

    def arccos(self) -> Series:
        """
        Compute the element-wise value for the inverse cosine.

        Examples
        --------
        >>> s = pl.Series("a", [1.0, 0.0, -1.0])
        >>> s.arccos()
        shape: (3,)
        Series: 'a' [f64]
        [
            0.0
            1.570796
            3.141593
        ]
        """

    def arctan(self) -> Series:
        """
        Compute the element-wise value for the inverse tangent.

        Examples
        --------
        >>> s = pl.Series("a", [1.0, 0.0, -1.0])
        >>> s.arctan()
        shape: (3,)
        Series: 'a' [f64]
        [
            0.785398
            0.0
            -0.785398
        ]
        """

    def arcsinh(self) -> Series:
        """
        Compute the element-wise value for the inverse hyperbolic sine.

        Examples
        --------
        >>> s = pl.Series("a", [1.0, 0.0, -1.0])
        >>> s.arcsinh()
        shape: (3,)
        Series: 'a' [f64]
        [
            0.881374
            0.0
            -0.881374
        ]
        """

    def arccosh(self) -> Series:
        """
        Compute the element-wise value for the inverse hyperbolic cosine.

        Examples
        --------
        >>> s = pl.Series("a", [5.0, 1.0, 0.0, -1.0])
        >>> s.arccosh()
        shape: (4,)
        Series: 'a' [f64]
        [
            2.292432
            0.0
            NaN
            NaN
        ]
        """

    def arctanh(self) -> Series:
        """
        Compute the element-wise value for the inverse hyperbolic tangent.

        Examples
        --------
        >>> s = pl.Series("a", [2.0, 1.0, 0.5, 0.0, -0.5, -1.0, -1.1])
        >>> s.arctanh()
        shape: (7,)
        Series: 'a' [f64]
        [
            NaN
            inf
            0.549306
            0.0
            -0.549306
            -inf
            NaN
        ]
        """

    def sinh(self) -> Series:
        """
        Compute the element-wise value for the hyperbolic sine.

        Examples
        --------
        >>> s = pl.Series("a", [1.0, 0.0, -1.0])
        >>> s.sinh()
        shape: (3,)
        Series: 'a' [f64]
        [
            1.175201
            0.0
            -1.175201
        ]
        """

    def cosh(self) -> Series:
        """
        Compute the element-wise value for the hyperbolic cosine.

        Examples
        --------
        >>> s = pl.Series("a", [1.0, 0.0, -1.0])
        >>> s.cosh()
        shape: (3,)
        Series: 'a' [f64]
        [
            1.543081
            1.0
            1.543081
        ]
        """

    def tanh(self) -> Series:
        """
        Compute the element-wise value for the hyperbolic tangent.

        Examples
        --------
        >>> s = pl.Series("a", [1.0, 0.0, -1.0])
        >>> s.tanh()
        shape: (3,)
        Series: 'a' [f64]
        [
            0.761594
            0.0
            -0.761594
        ]
        """

    def map_elements(
        self,
        function: Callable[[Any], Any],
        return_dtype: PolarsDataType | None = None,
        *,
        skip_nulls: bool = True,
    ) -> Self:
        """
        Map a custom/user-defined function (UDF) over elements in this Series.

        .. warning::
            This method is much slower than the native expressions API.
            Only use it if you cannot implement your logic otherwise.

            Suppose that the function is: `x ↦ sqrt(x)`:

            - For mapping elements of a series, consider: `s.sqrt()`.
            - For mapping inner elements of lists, consider:
              `s.list.eval(pl.element().sqrt())`.
            - For mapping elements of struct fields, consider:
              `s.struct.field("field_name").sqrt()`.

        If the function returns a different datatype, the return_dtype arg should
        be set, otherwise the method will fail.

        Implementing logic using a Python function is almost always *significantly*
        slower and more memory intensive than implementing the same logic using
        the native expression API because:

        - The native expression engine runs in Rust; UDFs run in Python.
        - Use of Python UDFs forces the DataFrame to be materialized in memory.
        - Polars-native expressions can be parallelised (UDFs typically cannot).
        - Polars-native expressions can be logically optimised (UDFs cannot).

        Wherever possible you should strongly prefer the native expression API
        to achieve the best performance.

        Parameters
        ----------
        function
            Custom function or lambda.
        return_dtype
            Output datatype.
            If not set, the dtype will be inferred based on the first non-null value
            that is returned by the function.
        skip_nulls
            Nulls will be skipped and not passed to the python function.
            This is faster because python can be skipped and because we call
            more specialized functions.

        Warnings
        --------
        If `return_dtype` is not provided, this may lead to unexpected results.
        We allow this, but it is considered a bug in the user's query.

        Notes
        -----
        * If your function is expensive and you don't want it to be called more than
          once for a given input, consider applying an `@lru_cache` decorator to it.
          If your data is suitable you may achieve *significant* speedups.

        * A UDF passed to `map_elements` must be pure, meaning that it cannot modify
          or depend on state other than its arguments.


        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.map_elements(lambda x: x + 10, return_dtype=pl.Int64)  # doctest: +SKIP
        shape: (3,)
        Series: 'a' [i64]
        [
                11
                12
                13
        ]

        Returns
        -------
        Series
        """
        from polars._utils.udfs import warn_on_inefficient_map

        if return_dtype is None:
            pl_return_dtype = None
        else:
            pl_return_dtype = parse_into_dtype(return_dtype)

        warn_on_inefficient_map(function, columns=[self.name], map_target="series")
        return self._from_pyseries(
            self._s.map_elements(
                function, return_dtype=pl_return_dtype, skip_nulls=skip_nulls
            )
        )

    def shift(self, n: int = 1, *, fill_value: IntoExpr | None = None) -> Series:
        """
        Shift values by the given number of indices.

        Parameters
        ----------
        n
            Number of indices to shift forward. If a negative value is passed, values
            are shifted in the opposite direction instead.
        fill_value
            Fill the resulting null values with this value. Accepts scalar expression
            input. Non-expression inputs are parsed as literals.

        Notes
        -----
        This method is similar to the `LAG` operation in SQL when the value for `n`
        is positive. With a negative value for `n`, it is similar to `LEAD`.

        Examples
        --------
        By default, values are shifted forward by one index.

        >>> s = pl.Series([1, 2, 3, 4])
        >>> s.shift()
        shape: (4,)
        Series: '' [i64]
        [
                null
                1
                2
                3
        ]

        Pass a negative value to shift in the opposite direction instead.

        >>> s.shift(-2)
        shape: (4,)
        Series: '' [i64]
        [
                3
                4
                null
                null
        ]

        Specify `fill_value` to fill the resulting null values.

        >>> s.shift(-2, fill_value=100)
        shape: (4,)
        Series: '' [i64]
        [
                3
                4
                100
                100
        ]
        """

    def zip_with(self, mask: Series, other: Series) -> Self:
        """
        Take values from self or other based on the given mask.

        Where mask evaluates true, take values from self. Where mask evaluates false,
        take values from other.

        Parameters
        ----------
        mask
            Boolean Series.
        other
            Series of same type.

        Returns
        -------
        Series

        Examples
        --------
        >>> s1 = pl.Series([1, 2, 3, 4, 5])
        >>> s2 = pl.Series([5, 4, 3, 2, 1])
        >>> s1.zip_with(s1 < s2, s2)
        shape: (5,)
        Series: '' [i64]
        [
                1
                2
                3
                2
                1
        ]
        >>> mask = pl.Series([True, False, True, False, True])
        >>> s1.zip_with(mask, s2)
        shape: (5,)
        Series: '' [i64]
        [
                1
                4
                3
                2
                5
        ]
        """
        require_same_type(self, other)
        return self._from_pyseries(self._s.zip_with(mask._s, other._s))

    @unstable()
    def rolling_min_by(
        self,
        by: IntoExpr,
        window_size: timedelta | str,
        *,
        min_samples: int = 1,
        closed: ClosedInterval = "right",
    ) -> Self:
        """
        Compute a rolling min based on another series.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Given a `by` column `<t_0, t_1, ..., t_n>`, then `closed="right"`
        (the default) means the windows will be:

            - (t_0 - window_size, t_0]
            - (t_1 - window_size, t_1]
            - ...
            - (t_n - window_size, t_n]

        Parameters
        ----------
        by
            Should be ``DateTime``, ``Date``, ``UInt64``, ``UInt32``, ``Int64``,
            or ``Int32`` data type (note that the integral ones require using `'i'`
            in `window size`).
        window_size
            The length of the window. Can be a dynamic temporal
            size indicated by a timedelta or the following string language:

            - 1ns   (1 nanosecond)
            - 1us   (1 microsecond)
            - 1ms   (1 millisecond)
            - 1s    (1 second)
            - 1m    (1 minute)
            - 1h    (1 hour)
            - 1d    (1 calendar day)
            - 1w    (1 calendar week)
            - 1mo   (1 calendar month)
            - 1q    (1 calendar quarter)
            - 1y    (1 calendar year)
            - 1i    (1 index count)

            By "calendar day", we mean the corresponding time on the next day
            (which may not be 24 hours, due to daylight savings). Similarly for
            "calendar week", "calendar month", "calendar quarter", and
            "calendar year".
        min_samples
            The number of values in the window that should be non-null before computing
            a result.
        closed : {'left', 'right', 'both', 'none'}
            Define which sides of the temporal interval are closed (inclusive),
            defaults to `'right'`.

        Notes
        -----
        If you want to compute multiple aggregation statistics over the same dynamic
        window, consider using `rolling` - this method can cache the window size
        computation.

        Examples
        --------
        Create a series with a row index value

        >>> from datetime import timedelta, datetime
        >>> start = datetime(2001, 1, 1)
        >>> stop = datetime(2001, 1, 2)
        >>> s = pl.Series("index", range(25))
        >>> s
        shape: (25,)
        Series: 'index' [i64]
        [
            0
            1
            2
            3
            4
            …
            20
            21
            22
            23
            24
        ]

        Create another series to apply the window mask:

        >>> d = pl.Series("date", pl.datetime_range(start, stop, "1h", eager=True))
        >>> d
        shape: (25,)
        Series: 'date' [datetime[μs]]
        [
            2001-01-01 00:00:00
            2001-01-01 01:00:00
            2001-01-01 02:00:00
            2001-01-01 03:00:00
            2001-01-01 04:00:00
            …
            2001-01-01 20:00:00
            2001-01-01 21:00:00
            2001-01-01 22:00:00
            2001-01-01 23:00:00
            2001-01-02 00:00:00
        ]

        Compute the rolling min with the temporal windows
        from the second series closed on the right:

        >>> s.rolling_min_by(d, "3h")
        shape: (25,)
        Series: 'index' [i64]
        [
            0
            0
            0
            1
            2
            …
            18
            19
            20
            21
            22
        ]
        """

    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def rolling_min(
        self,
        window_size: int,
        weights: list[float] | None = None,
        *,
        min_samples: int | None = None,
        center: bool = False,
    ) -> Series:
        """
        Apply a rolling min (moving min) over the values in this array.

        A window of length `window_size` will traverse the array. The values that fill
        this window will (optionally) be multiplied with the weights given by the
        `weight` vector. The resulting values will be aggregated to their min.

        The window at a given row will include the row itself and the `window_size - 1`
        elements before it.

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        window_size
            The length of the window in number of elements.
        weights
            An optional slice with the same length as the window that will be multiplied
            elementwise with the values in the window.
        min_samples
            The number of values in the window that should be non-null before computing
            a result. If set to `None` (default), it will be set equal to `window_size`.
        center
            Set the labels at the center of the window.

        Examples
        --------
        >>> s = pl.Series("a", [100, 200, 300, 400, 500])
        >>> s.rolling_min(window_size=3)
        shape: (5,)
        Series: 'a' [i64]
        [
            null
            null
            100
            200
            300
        ]
        """

    @unstable()
    def rolling_max_by(
        self,
        by: IntoExpr,
        window_size: timedelta | str,
        *,
        min_samples: int = 1,
        closed: ClosedInterval = "right",
    ) -> Self:
        """
        Compute a rolling max based on another series.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Given a `by` column `<t_0, t_1, ..., t_n>`, then `closed="right"`
        (the default) means the windows will be:

            - (t_0 - window_size, t_0]
            - (t_1 - window_size, t_1]
            - ...
            - (t_n - window_size, t_n]

        Parameters
        ----------
        by
            Should be ``DateTime``, ``Date``, ``UInt64``, ``UInt32``, ``Int64``,
            or ``Int32`` data type (note that the integral ones require using `'i'`
            in `window size`).
        window_size
            The length of the window. Can be a dynamic temporal
            size indicated by a timedelta or the following string language:

            - 1ns   (1 nanosecond)
            - 1us   (1 microsecond)
            - 1ms   (1 millisecond)
            - 1s    (1 second)
            - 1m    (1 minute)
            - 1h    (1 hour)
            - 1d    (1 calendar day)
            - 1w    (1 calendar week)
            - 1mo   (1 calendar month)
            - 1q    (1 calendar quarter)
            - 1y    (1 calendar year)
            - 1i    (1 index count)

            By "calendar day", we mean the corresponding time on the next day
            (which may not be 24 hours, due to daylight savings). Similarly for
            "calendar week", "calendar month", "calendar quarter", and
            "calendar year".
        min_samples
            The number of values in the window that should be non-null before computing
            a result.
        closed : {'left', 'right', 'both', 'none'}
            Define which sides of the temporal interval are closed (inclusive),
            defaults to `'right'`.

        Notes
        -----
        If you want to compute multiple aggregation statistics over the same dynamic
        window, consider using `rolling` - this method can cache the window size
        computation.

        Examples
        --------
        Create a series with a row index value

        >>> from datetime import timedelta, datetime
        >>> start = datetime(2001, 1, 1)
        >>> stop = datetime(2001, 1, 2)
        >>> s = pl.Series("index", range(25))
        >>> s
        shape: (25,)
        Series: 'index' [i64]
        [
            0
            1
            2
            3
            4
            …
            20
            21
            22
            23
            24
        ]

        Create another series to apply the window mask:

        >>> d = pl.Series("date", pl.datetime_range(start, stop, "1h", eager=True))
        >>> d
        shape: (25,)
        Series: 'date' [datetime[μs]]
        [
            2001-01-01 00:00:00
            2001-01-01 01:00:00
            2001-01-01 02:00:00
            2001-01-01 03:00:00
            2001-01-01 04:00:00
            …
            2001-01-01 20:00:00
            2001-01-01 21:00:00
            2001-01-01 22:00:00
            2001-01-01 23:00:00
            2001-01-02 00:00:00
        ]

        Compute the rolling max with the temporal windows
        from the second series closed on the right:

        >>> s.rolling_max_by(d, "3h")
        shape: (25,)
        Series: 'index' [i64]
        [
            0
            1
            2
            3
            4
            …
            20
            21
            22
            23
            24
        ]
        """

    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def rolling_max(
        self,
        window_size: int,
        weights: list[float] | None = None,
        *,
        min_samples: int | None = None,
        center: bool = False,
    ) -> Series:
        """
        Apply a rolling max (moving max) over the values in this array.

        A window of length `window_size` will traverse the array. The values that fill
        this window will (optionally) be multiplied with the weights given by the
        `weight` vector. The resulting values will be aggregated to their max.

        The window at a given row will include the row itself and the `window_size - 1`
        elements before it.

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        window_size
            The length of the window in number of elements.
        weights
            An optional slice with the same length as the window that will be multiplied
            elementwise with the values in the window.
        min_samples
            The number of values in the window that should be non-null before computing
            a result. If set to `None` (default), it will be set equal to `window_size`.
        center
            Set the labels at the center of the window.

        Examples
        --------
        >>> s = pl.Series("a", [100, 200, 300, 400, 500])
        >>> s.rolling_max(window_size=2)
        shape: (5,)
        Series: 'a' [i64]
        [
            null
            200
            300
            400
            500
        ]
        """

    @unstable()
    def rolling_mean_by(
        self,
        by: IntoExpr,
        window_size: timedelta | str,
        *,
        min_samples: int = 1,
        closed: ClosedInterval = "right",
    ) -> Self:
        """
        Compute a rolling mean based on another series.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Given a `by` column `<t_0, t_1, ..., t_n>`, then `closed="right"`
        (the default) means the windows will be:

            - (t_0 - window_size, t_0]
            - (t_1 - window_size, t_1]
            - ...
            - (t_n - window_size, t_n]

        Parameters
        ----------
        by
            Should be ``DateTime``, ``Date``, ``UInt64``, ``UInt32``, ``Int64``,
            or ``Int32`` data type (note that the integral ones require using `'i'`
            in `window size`).
        window_size
            The length of the window. Can be a dynamic temporal
            size indicated by a timedelta or the following string language:

            - 1ns   (1 nanosecond)
            - 1us   (1 microsecond)
            - 1ms   (1 millisecond)
            - 1s    (1 second)
            - 1m    (1 minute)
            - 1h    (1 hour)
            - 1d    (1 calendar day)
            - 1w    (1 calendar week)
            - 1mo   (1 calendar month)
            - 1q    (1 calendar quarter)
            - 1y    (1 calendar year)
            - 1i    (1 index count)

            By "calendar day", we mean the corresponding time on the next day
            (which may not be 24 hours, due to daylight savings). Similarly for
            "calendar week", "calendar month", "calendar quarter", and
            "calendar year".
        min_samples
            The number of values in the window that should be non-null before computing
            a result.
        closed : {'left', 'right', 'both', 'none'}
            Define which sides of the temporal interval are closed (inclusive),
            defaults to `'right'`.

        Notes
        -----
        If you want to compute multiple aggregation statistics over the same dynamic
        window, consider using `rolling` - this method can cache the window size
        computation.

        Examples
        --------
        Create a series with a row index value

        >>> from datetime import timedelta, datetime
        >>> start = datetime(2001, 1, 1)
        >>> stop = datetime(2001, 1, 2)
        >>> s = pl.Series("index", range(25))
        >>> s
        shape: (25,)
        Series: 'index' [i64]
        [
            0
            1
            2
            3
            4
            …
            20
            21
            22
            23
            24
        ]

        Create another series to apply the window mask:

        >>> d = pl.Series("date", pl.datetime_range(start, stop, "1h", eager=True))
        >>> d
        shape: (25,)
        Series: 'date' [datetime[μs]]
        [
            2001-01-01 00:00:00
            2001-01-01 01:00:00
            2001-01-01 02:00:00
            2001-01-01 03:00:00
            2001-01-01 04:00:00
            …
            2001-01-01 20:00:00
            2001-01-01 21:00:00
            2001-01-01 22:00:00
            2001-01-01 23:00:00
            2001-01-02 00:00:00
        ]

        Compute the rolling mean with the temporal windows
        from the second series closed on the right:

        >>> s.rolling_mean_by(d, "3h")
        shape: (25,)
        Series: 'index' [f64]
        [
            0.0
            0.5
            1.0
            2.0
            3.0
            …
            19.0
            20.0
            21.0
            22.0
            23.0
        ]
        """

    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def rolling_mean(
        self,
        window_size: int,
        weights: list[float] | None = None,
        *,
        min_samples: int | None = None,
        center: bool = False,
    ) -> Series:
        """
        Apply a rolling mean (moving mean) over the values in this array.

        A window of length `window_size` will traverse the array. The values that fill
        this window will (optionally) be multiplied with the weights given by the
        `weight` vector. The resulting values will be aggregated to their mean.

        The window at a given row will include the row itself and the `window_size - 1`
        elements before it.

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        window_size
            The length of the window in number of elements.
        weights
            An optional slice with the same length as the window that will be multiplied
            elementwise with the values in the window.
        min_samples
            The number of values in the window that should be non-null before computing
            a result. If set to `None` (default), it will be set equal to `window_size`.
        center
            Set the labels at the center of the window.

        Examples
        --------
        >>> s = pl.Series("a", [100, 200, 300, 400, 500])
        >>> s.rolling_mean(window_size=2)
        shape: (5,)
        Series: 'a' [f64]
        [
            null
            150.0
            250.0
            350.0
            450.0
        ]
        """

    @unstable()
    def rolling_sum_by(
        self,
        by: IntoExpr,
        window_size: timedelta | str,
        *,
        min_samples: int = 1,
        closed: ClosedInterval = "right",
    ) -> Self:
        """
        Compute a rolling sum based on another series.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Given a `by` column `<t_0, t_1, ..., t_n>`, then `closed="right"`
        (the default) means the windows will be:

            - (t_0 - window_size, t_0]
            - (t_1 - window_size, t_1]
            - ...
            - (t_n - window_size, t_n]

        Parameters
        ----------
        window_size
            The length of the window. Can be a dynamic temporal
            size indicated by a timedelta or the following string language:

            - 1ns   (1 nanosecond)
            - 1us   (1 microsecond)
            - 1ms   (1 millisecond)
            - 1s    (1 second)
            - 1m    (1 minute)
            - 1h    (1 hour)
            - 1d    (1 calendar day)
            - 1w    (1 calendar week)
            - 1mo   (1 calendar month)
            - 1q    (1 calendar quarter)
            - 1y    (1 calendar year)
            - 1i    (1 index count)

            By "calendar day", we mean the corresponding time on the next day
            (which may not be 24 hours, due to daylight savings). Similarly for
            "calendar week", "calendar month", "calendar quarter", and
            "calendar year".
        min_samples
            The number of values in the window that should be non-null before computing
            a result.
        by
            Should be ``DateTime``, ``Date``, ``UInt64``, ``UInt32``, ``Int64``,
            or ``Int32`` data type (note that the integral ones require using `'i'`
            in `window size`).
        closed : {'left', 'right', 'both', 'none'}
            Define which sides of the temporal interval are closed (inclusive),
            defaults to `'right'`.

        Notes
        -----
        If you want to compute multiple aggregation statistics over the same dynamic
        window, consider using `rolling` - this method can cache the window size
        computation.

        Examples
        --------
        Create a series with a row index value

        >>> from datetime import timedelta, datetime
        >>> start = datetime(2001, 1, 1)
        >>> stop = datetime(2001, 1, 2)
        >>> s = pl.Series("index", range(25))
        >>> s
        shape: (25,)
        Series: 'index' [i64]
        [
            0
            1
            2
            3
            4
            …
            20
            21
            22
            23
            24
        ]

        Create another series to apply the window mask:

        >>> d = pl.Series("date", pl.datetime_range(start, stop, "1h", eager=True))
        >>> d
        shape: (25,)
        Series: 'date' [datetime[μs]]
        [
            2001-01-01 00:00:00
            2001-01-01 01:00:00
            2001-01-01 02:00:00
            2001-01-01 03:00:00
            2001-01-01 04:00:00
            …
            2001-01-01 20:00:00
            2001-01-01 21:00:00
            2001-01-01 22:00:00
            2001-01-01 23:00:00
            2001-01-02 00:00:00
        ]

        Compute the rolling mean with the temporal windows
        from the second series closed on the right:

        >>> s.rolling_sum_by(d, "3h")
        shape: (25,)
        Series: 'index' [i64]
        [
            0
            1
            3
            6
            9
            …
            57
            60
            63
            66
            69
        ]
        """

    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def rolling_sum(
        self,
        window_size: int,
        weights: list[float] | None = None,
        *,
        min_samples: int | None = None,
        center: bool = False,
    ) -> Series:
        """
        Apply a rolling sum (moving sum) over the values in this array.

        A window of length `window_size` will traverse the array. The values that fill
        this window will (optionally) be multiplied with the weights given by the
        `weight` vector. The resulting values will be aggregated to their sum.

        The window at a given row will include the row itself and the `window_size - 1`
        elements before it.

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        window_size
            The length of the window in number of elements.
        weights
            An optional slice with the same length as the window that will be multiplied
            elementwise with the values in the window.
        min_samples
            The number of values in the window that should be non-null before computing
            a result. If set to `None` (default), it will be set equal to `window_size`.
        center
            Set the labels at the center of the window.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3, 4, 5])
        >>> s.rolling_sum(window_size=2)
        shape: (5,)
        Series: 'a' [i64]
        [
                null
                3
                5
                7
                9
        ]
        """

    @unstable()
    def rolling_std_by(
        self,
        by: IntoExpr,
        window_size: timedelta | str,
        *,
        min_samples: int = 1,
        closed: ClosedInterval = "right",
        ddof: int = 1,
    ) -> Self:
        """
        Compute a rolling standard deviation based on another series.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Given a `by` column `<t_0, t_1, ..., t_n>`, then `closed="right"`
        (the default) means the windows will be:

            - (t_0 - window_size, t_0]
            - (t_1 - window_size, t_1]
            - ...
            - (t_n - window_size, t_n]

        Parameters
        ----------
        by
            Should be ``DateTime``, ``Date``, ``UInt64``, ``UInt32``, ``Int64``,
            or ``Int32`` data type (note that the integral ones require using `'i'`
            in `window size`).
        window_size
            The length of the window. Can be a dynamic temporal
            size indicated by a timedelta or the following string language:

            - 1ns   (1 nanosecond)
            - 1us   (1 microsecond)
            - 1ms   (1 millisecond)
            - 1s    (1 second)
            - 1m    (1 minute)
            - 1h    (1 hour)
            - 1d    (1 calendar day)
            - 1w    (1 calendar week)
            - 1mo   (1 calendar month)
            - 1q    (1 calendar quarter)
            - 1y    (1 calendar year)
            - 1i    (1 index count)

            By "calendar day", we mean the corresponding time on the next day
            (which may not be 24 hours, due to daylight savings). Similarly for
            "calendar week", "calendar month", "calendar quarter", and
            "calendar year".
        min_samples
            The number of values in the window that should be non-null before computing
            a result.
        closed : {'left', 'right', 'both', 'none'}
            Define which sides of the temporal interval are closed (inclusive),
            defaults to `'right'`.
        ddof
            "Delta Degrees of Freedom": The divisor for a length N window is N - ddof

        Notes
        -----
        If you want to compute multiple aggregation statistics over the same dynamic
        window, consider using `rolling` - this method can cache the window size
        computation.

        Examples
        --------
        Create a series with a row index value

        >>> from datetime import timedelta, datetime
        >>> start = datetime(2001, 1, 1)
        >>> stop = datetime(2001, 1, 2)
        >>> s = pl.Series("index", range(25))
        >>> s
        shape: (25,)
        Series: 'index' [i64]
        [
            0
            1
            2
            3
            4
            …
            20
            21
            22
            23
            24
        ]

        Create another series to apply the window mask:

        >>> d = pl.Series("date", pl.datetime_range(start, stop, "1h", eager=True))
        >>> d
        shape: (25,)
        Series: 'date' [datetime[μs]]
        [
                2001-01-01 00:00:00
                2001-01-01 01:00:00
                2001-01-01 02:00:00
                2001-01-01 03:00:00
                2001-01-01 04:00:00
                …
                2001-01-01 20:00:00
                2001-01-01 21:00:00
                2001-01-01 22:00:00
                2001-01-01 23:00:00
                2001-01-02 00:00:00
        ]

        Compute the rolling std with the temporal windows
        from the second series closed on the right:

        >>> s.rolling_std_by(d, "3h")
        shape: (25,)
        Series: 'index' [f64]
        [
            null
            0.707107
            1.0
            1.0
            1.0
            …
            1.0
            1.0
            1.0
            1.0
            1.0
        ]
        """

    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def rolling_std(
        self,
        window_size: int,
        weights: list[float] | None = None,
        *,
        min_samples: int | None = None,
        center: bool = False,
        ddof: int = 1,
    ) -> Series:
        """
        Compute a rolling std dev.

        A window of length `window_size` will traverse the array. The values that fill
        this window will (optionally) be multiplied with the weights given by the
        `weight` vector. The resulting values will be aggregated to their std dev.

        The window at a given row will include the row itself and the `window_size - 1`
        elements before it.

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        window_size
            The length of the window in number of elements.
        weights
            An optional slice with the same length as the window that will be multiplied
            elementwise with the values in the window.
        min_samples
            The number of values in the window that should be non-null before computing
            a result. If set to `None` (default), it will be set equal to `window_size`.
        center
            Set the labels at the center of the window.
        ddof
            "Delta Degrees of Freedom": The divisor for a length N window is N - ddof

        Examples
        --------
        >>> s = pl.Series("a", [1.0, 2.0, 3.0, 4.0, 6.0, 8.0])
        >>> s.rolling_std(window_size=3)
        shape: (6,)
        Series: 'a' [f64]
        [
                null
                null
                1.0
                1.0
                1.527525
                2.0
        ]
        """

    @unstable()
    def rolling_var_by(
        self,
        by: IntoExpr,
        window_size: timedelta | str,
        *,
        min_samples: int = 1,
        closed: ClosedInterval = "right",
        ddof: int = 1,
    ) -> Self:
        """
        Compute a rolling variance based on another series.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Given a `by` column `<t_0, t_1, ..., t_n>`, then `closed="right"`
        (the default) means the windows will be:

            - (t_0 - window_size, t_0]
            - (t_1 - window_size, t_1]
            - ...
            - (t_n - window_size, t_n]

        Parameters
        ----------
        by
            Should be ``DateTime``, ``Date``, ``UInt64``, ``UInt32``, ``Int64``,
            or ``Int32`` data type (note that the integral ones require using `'i'`
            in `window size`).
        window_size
            The length of the window. Can be a dynamic temporal
            size indicated by a timedelta or the following string language:

            - 1ns   (1 nanosecond)
            - 1us   (1 microsecond)
            - 1ms   (1 millisecond)
            - 1s    (1 second)
            - 1m    (1 minute)
            - 1h    (1 hour)
            - 1d    (1 calendar day)
            - 1w    (1 calendar week)
            - 1mo   (1 calendar month)
            - 1q    (1 calendar quarter)
            - 1y    (1 calendar year)
            - 1i    (1 index count)

            By "calendar day", we mean the corresponding time on the next day
            (which may not be 24 hours, due to daylight savings). Similarly for
            "calendar week", "calendar month", "calendar quarter", and
            "calendar year".
        min_samples
            The number of values in the window that should be non-null before computing
            a result.
        closed : {'left', 'right', 'both', 'none'}
            Define which sides of the temporal interval are closed (inclusive),
            defaults to `'right'`.
        ddof
            "Delta Degrees of Freedom": The divisor for a length N window is N - ddof

        Notes
        -----
        If you want to compute multiple aggregation statistics over the same dynamic
        window, consider using `rolling` - this method can cache the window size
        computation.

        Examples
        --------
        Create a series with a row index value

        >>> from datetime import timedelta, datetime
        >>> start = datetime(2001, 1, 1)
        >>> stop = datetime(2001, 1, 2)
        >>> s = pl.Series("index", range(25))
        >>> s
        shape: (25,)
        Series: 'index' [i64]
        [
            0
            1
            2
            3
            4
            …
            20
            21
            22
            23
            24
        ]

        Create another series to apply the window mask:

        >>> d = pl.Series("date", pl.datetime_range(start, stop, "1h", eager=True))
        >>> d
        shape: (25,)
        Series: 'date' [datetime[μs]]
        [
            2001-01-01 00:00:00
            2001-01-01 01:00:00
            2001-01-01 02:00:00
            2001-01-01 03:00:00
            2001-01-01 04:00:00
            …
            2001-01-01 20:00:00
            2001-01-01 21:00:00
            2001-01-01 22:00:00
            2001-01-01 23:00:00
            2001-01-02 00:00:00
        ]

        Compute the rolling std with the temporal windows
        from the second series closed on the right:

        >>> s.rolling_std_by(d, "3h")
        shape: (25,)
        Series: 'index' [f64]
        [
            null
            0.707107
            1.0
            1.0
            1.0
            …
            1.0
            1.0
            1.0
            1.0
            1.0
        ]
        """

    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def rolling_var(
        self,
        window_size: int,
        weights: list[float] | None = None,
        *,
        min_samples: int | None = None,
        center: bool = False,
        ddof: int = 1,
    ) -> Series:
        """
        Compute a rolling variance.

        A window of length `window_size` will traverse the array. The values that fill
        this window will (optionally) be multiplied with the weights given by the
        `weight` vector. The resulting values will be aggregated to their variance.

        The window at a given row will include the row itself and the `window_size - 1`
        elements before it.

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        window_size
            The length of the window in number of elements.
        weights
            An optional slice with the same length as the window that will be multiplied
            elementwise with the values in the window.
        min_samples
            The number of values in the window that should be non-null before computing
            a result. If set to `None` (default), it will be set equal to `window_size`.
        center
            Set the labels at the center of the window.
        ddof
            "Delta Degrees of Freedom": The divisor for a length N window is N - ddof

        Examples
        --------
        >>> s = pl.Series("a", [1.0, 2.0, 3.0, 4.0, 6.0, 8.0])
        >>> s.rolling_var(window_size=3)
        shape: (6,)
        Series: 'a' [f64]
        [
                null
                null
                1.0
                1.0
                2.333333
                4.0
        ]
        """

    @unstable()
    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def rolling_map(
        self,
        function: Callable[[Series], Any],
        window_size: int,
        weights: list[float] | None = None,
        *,
        min_samples: int | None = None,
        center: bool = False,
    ) -> Series:
        """
        Compute a custom rolling window function.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        function
            Custom aggregation function.
        window_size
            The length of the window in number of elements.
        weights
            An optional slice with the same length as the window that will be multiplied
            elementwise with the values in the window.
        min_samples
            The number of values in the window that should be non-null before computing
            a result. If set to `None` (default), it will be set equal to `window_size`.
        center
            Set the labels at the center of the window.

        Warnings
        --------
        Computing custom functions is extremely slow. Use specialized rolling
        functions such as :func:`Series.rolling_sum` if at all possible.

        Examples
        --------
        >>> from numpy import nansum
        >>> s = pl.Series([11.0, 2.0, 9.0, float("nan"), 8.0])
        >>> s.rolling_map(nansum, window_size=3)
        shape: (5,)
        Series: '' [f64]
        [
                null
                null
                22.0
                11.0
                17.0
        ]
        """

    @unstable()
    def rolling_median_by(
        self,
        by: IntoExpr,
        window_size: timedelta | str,
        *,
        min_samples: int = 1,
        closed: ClosedInterval = "right",
    ) -> Self:
        """
        Compute a rolling median based on another series.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Given a `by` column `<t_0, t_1, ..., t_n>`, then `closed="right"`
        (the default) means the windows will be:

            - (t_0 - window_size, t_0]
            - (t_1 - window_size, t_1]
            - ...
            - (t_n - window_size, t_n]

        Parameters
        ----------
        by
            Should be ``DateTime``, ``Date``, ``UInt64``, ``UInt32``, ``Int64``,
            or ``Int32`` data type (note that the integral ones require using `'i'`
            in `window size`).
        window_size
            The length of the window. Can be a dynamic temporal
            size indicated by a timedelta or the following string language:

            - 1ns   (1 nanosecond)
            - 1us   (1 microsecond)
            - 1ms   (1 millisecond)
            - 1s    (1 second)
            - 1m    (1 minute)
            - 1h    (1 hour)
            - 1d    (1 calendar day)
            - 1w    (1 calendar week)
            - 1mo   (1 calendar month)
            - 1q    (1 calendar quarter)
            - 1y    (1 calendar year)
            - 1i    (1 index count)

            By "calendar day", we mean the corresponding time on the next day
            (which may not be 24 hours, due to daylight savings). Similarly for
            "calendar week", "calendar month", "calendar quarter", and
            "calendar year".
        min_samples
            The number of values in the window that should be non-null before computing
            a result.
        closed : {'left', 'right', 'both', 'none'}
            Define which sides of the temporal interval are closed (inclusive),
            defaults to `'right'`.

        Notes
        -----
        If you want to compute multiple aggregation statistics over the same dynamic
        window, consider using `rolling` - this method can cache the window size
        computation.

        Examples
        --------
        Create a series with a row index value

        >>> from datetime import timedelta, datetime
        >>> start = datetime(2001, 1, 1)
        >>> stop = datetime(2001, 1, 2)
        >>> s = pl.Series("index", range(25))
        >>> s
        shape: (25,)
        Series: 'index' [i64]
        [
            0
            1
            2
            3
            4
            …
            20
            21
            22
            23
            24
        ]

        Create another series to apply the window mask:

        >>> d = pl.Series("date", pl.datetime_range(start, stop, "1h", eager=True))
        >>> d
        shape: (25,)
        Series: 'date' [datetime[μs]]
        [
            2001-01-01 00:00:00
            2001-01-01 01:00:00
            2001-01-01 02:00:00
            2001-01-01 03:00:00
            2001-01-01 04:00:00
            …
            2001-01-01 20:00:00
            2001-01-01 21:00:00
            2001-01-01 22:00:00
            2001-01-01 23:00:00
            2001-01-02 00:00:00
        ]

        Compute the rolling median with the temporal windows
        from the second series closed on the right:

        >>> s.rolling_median_by(d, "3h")
        shape: (25,)
        Series: 'index' [f64]
        [
            0.0
            0.5
            1.0
            2.0
            3.0
            …
            19.0
            20.0
            21.0
            22.0
            23.0
        ]
        """

    @unstable()
    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def rolling_median(
        self,
        window_size: int,
        weights: list[float] | None = None,
        *,
        min_samples: int | None = None,
        center: bool = False,
    ) -> Series:
        """
        Compute a rolling median.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        The window at a given row will include the row itself and the `window_size - 1`
        elements before it.

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        window_size
            The length of the window in number of elements.
        weights
            An optional slice with the same length as the window that will be multiplied
            elementwise with the values in the window.
        min_samples
            The number of values in the window that should be non-null before computing
            a result. If set to `None` (default), it will be set equal to `window_size`.
        center
            Set the labels at the center of the window.

        Examples
        --------
        >>> s = pl.Series("a", [1.0, 2.0, 3.0, 4.0, 6.0, 8.0])
        >>> s.rolling_median(window_size=3)
        shape: (6,)
        Series: 'a' [f64]
        [
                null
                null
                2.0
                3.0
                4.0
                6.0
        ]
        """

    @unstable()
    def rolling_quantile_by(
        self,
        by: IntoExpr,
        window_size: timedelta | str,
        *,
        quantile: float,
        interpolation: QuantileMethod = "nearest",
        min_samples: int = 1,
        closed: ClosedInterval = "right",
    ) -> Self:
        """
        Compute a rolling quantile based on another series.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Given a `by` column `<t_0, t_1, ..., t_n>`, then `closed="right"`
        (the default) means the windows will be:

            - (t_0 - window_size, t_0]
            - (t_1 - window_size, t_1]
            - ...
            - (t_n - window_size, t_n]

        Parameters
        ----------
        by
            Should be ``DateTime``, ``Date``, ``UInt64``, ``UInt32``, ``Int64``,
            or ``Int32`` data type (note that the integral ones require using `'i'`
            in `window size`).
        quantile
            Quantile between 0.0 and 1.0.
        interpolation : {'nearest', 'higher', 'lower', 'midpoint', 'linear', 'equiprobable'}
            Interpolation method.
        window_size
            The length of the window. Can be a dynamic
            temporal size indicated by a timedelta or the following string language:

            - 1ns   (1 nanosecond)
            - 1us   (1 microsecond)
            - 1ms   (1 millisecond)
            - 1s    (1 second)
            - 1m    (1 minute)
            - 1h    (1 hour)
            - 1d    (1 calendar day)
            - 1w    (1 calendar week)
            - 1mo   (1 calendar month)
            - 1q    (1 calendar quarter)
            - 1y    (1 calendar year)
            - 1i    (1 index count)

            By "calendar day", we mean the corresponding time on the next day
            (which may not be 24 hours, due to daylight savings). Similarly for
            "calendar week", "calendar month", "calendar quarter", and
            "calendar year".
        min_samples
            The number of values in the window that should be non-null before computing
            a result.
        closed : {'left', 'right', 'both', 'none'}
            Define which sides of the temporal interval are closed (inclusive),
            defaults to `'right'`.

        Notes
        -----
        If you want to compute multiple aggregation statistics over the same dynamic
        window, consider using `rolling` - this method can cache the window size
        computation.

        Examples
        --------
        Create a series with a row index value

        >>> from datetime import timedelta, datetime
        >>> start = datetime(2001, 1, 1)
        >>> stop = datetime(2001, 1, 2)
        >>> s = pl.Series("index", range(25))
        >>> s
        shape: (25,)
        Series: 'index' [i64]
        [
            0
            1
            2
            3
            4
            …
            20
            21
            22
            23
            24
        ]

        Create another series to apply the window mask:

        >>> d = pl.Series("date", pl.datetime_range(start, stop, "1h", eager=True))
        >>> d
        shape: (25,)
        Series: 'date' [datetime[μs]]
        [
            2001-01-01 00:00:00
            2001-01-01 01:00:00
            2001-01-01 02:00:00
            2001-01-01 03:00:00
            2001-01-01 04:00:00
            …
            2001-01-01 20:00:00
            2001-01-01 21:00:00
            2001-01-01 22:00:00
            2001-01-01 23:00:00
            2001-01-02 00:00:00
        ]

        Compute the rolling quantile with the temporal windows from the second series closed on the right:

        >>> s.rolling_quantile_by(d, "3h", quantile=0.5)
        shape: (25,)
        Series: 'index' [f64]
        [
            0.0
            1.0
            1.0
            2.0
            3.0
            …
            19.0
            20.0
            21.0
            22.0
            23.0
        ]
        """  # noqa: W505

    @unstable()
    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def rolling_quantile(
        self,
        quantile: float,
        interpolation: QuantileMethod = "nearest",
        window_size: int = 2,
        weights: list[float] | None = None,
        *,
        min_samples: int | None = None,
        center: bool = False,
    ) -> Series:
        """
        Compute a rolling quantile.

        The window at a given row will include the row itself and the `window_size - 1`
        elements before it.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        quantile
            Quantile between 0.0 and 1.0.
        interpolation : {'nearest', 'higher', 'lower', 'midpoint', 'linear', 'equiprobable'}
            Interpolation method.
        window_size
            The length of the window in number of elements.
        weights
            An optional slice with the same length as the window that will be multiplied
            elementwise with the values in the window.
        min_samples
            The number of values in the window that should be non-null before computing
            a result. If set to `None` (default), it will be set equal to `window_size`.
        center
            Set the labels at the center of the window.

        Examples
        --------
        >>> s = pl.Series("a", [1.0, 2.0, 3.0, 4.0, 6.0, 8.0])
        >>> s.rolling_quantile(quantile=0.33, window_size=3)
        shape: (6,)
        Series: 'a' [f64]
        [
                null
                null
                2.0
                3.0
                4.0
                6.0
        ]
        >>> s.rolling_quantile(quantile=0.33, interpolation="linear", window_size=3)
        shape: (6,)
        Series: 'a' [f64]
        [
                null
                null
                1.66
                2.66
                3.66
                5.32
        ]
        """  # noqa: W505

    @unstable()
    def rolling_rank_by(
        self,
        by: IntoExpr,
        window_size: timedelta | str,
        method: RankMethod = "average",
        *,
        seed: int | None = None,
        min_samples: int = 1,
        closed: ClosedInterval = "right",
    ) -> Series:
        """
        Compute a rolling rank based on another column.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        Given a `by` column `<t_0, t_1, ..., t_n>`, then `closed="right"`
        (the default) means the windows will be:

            - (t_0 - window_size, t_0]
            - (t_1 - window_size, t_1]
            - ...
            - (t_n - window_size, t_n]

        Parameters
        ----------
        by
            Should be ``DateTime``, ``Date``, ``UInt64``, ``UInt32``, ``Int64``,
            or ``Int32`` data type (note that the integral ones require using `'i'`
            in `window size`).
        window_size
            The length of the window. Can be a dynamic
            temporal size indicated by a timedelta or the following string language:

            - 1ns   (1 nanosecond)
            - 1us   (1 microsecond)
            - 1ms   (1 millisecond)
            - 1s    (1 second)
            - 1m    (1 minute)
            - 1h    (1 hour)
            - 1d    (1 calendar day)
            - 1w    (1 calendar week)
            - 1mo   (1 calendar month)
            - 1q    (1 calendar quarter)
            - 1y    (1 calendar year)
            - 1i    (1 index count)

            By "calendar day", we mean the corresponding time on the next day
            (which may not be 24 hours, due to daylight savings). Similarly for
            "calendar week", "calendar month", "calendar quarter", and
            "calendar year".
        method : {'average', 'min', 'max', 'dense', 'random'}
            The method used to assign ranks to tied elements.
            The following methods are available (default is 'average'):

            - 'average' : The average of the ranks that would have been assigned to
              all the tied values is assigned to each value.
            - 'min' : The minimum of the ranks that would have been assigned to all
              the tied values is assigned to each value. (This is also referred to
              as "competition" ranking.)
            - 'max' : The maximum of the ranks that would have been assigned to all
              the tied values is assigned to each value.
            - 'dense' : Like 'min', but the rank of the next highest element is
              assigned the rank immediately after those assigned to the tied
              elements.
            - 'random' : Choose a random rank for each value in a tie.
        seed
            Random seed used when `method='random'`. If set to None (default), a
            random seed is generated for each rolling rank operation.
        min_samples
            The number of values in the window that should be non-null before computing
            a result.
        closed : {'left', 'right', 'both', 'none'}
            Define which sides of the temporal interval are closed (inclusive),
            defaults to `'right'`.

        Returns
        -------
        Series
            A Series of data :class:`.Float64` if `method` is `"average"` or,
            the index size (see :func:`.get_index_type()`) otherwise.
        """

    @unstable()
    def rolling_rank(
        self,
        window_size: int,
        method: RankMethod = "average",
        *,
        seed: int | None = None,
        min_samples: int | None = None,
        center: bool = False,
    ) -> Series:
        """
        Compute a rolling rank.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        A window of length `window_size` will traverse the array. The values
        that fill this window will be ranked according to the `method`
        parameter. The resulting values will be the rank of the value that is
        at the end of the sliding window.

        Parameters
        ----------
        window_size
            Integer size of the rolling window.
        method : {'average', 'min', 'max', 'dense', 'random'}
            The method used to assign ranks to tied elements.
            The following methods are available (default is 'average'):

            - 'average' : The average of the ranks that would have been assigned to
              all the tied values is assigned to each value.
            - 'min' : The minimum of the ranks that would have been assigned to all
              the tied values is assigned to each value. (This is also referred to
              as "competition" ranking.)
            - 'max' : The maximum of the ranks that would have been assigned to all
              the tied values is assigned to each value.
            - 'dense' : Like 'min', but the rank of the next highest element is
              assigned the rank immediately after those assigned to the tied
              elements.
            - 'random' : Choose a random rank for each value in a tie.
        seed
            Random seed used when `method='random'`. If set to None (default), a
            random seed is generated for each rolling rank operation.
        min_samples
            The number of values in the window that should be non-null before computing
            a result. If set to `None` (default), it will be set equal to `window_size`.
        center
            Set the labels at the center of the window.

        Returns
        -------
        Series
            A Series of data :class:`.Float64` if `method` is `"average"` or,
            the index size (see :func:`.get_index_type()`) otherwise.

        Examples
        --------
        >>> pl.Series([1, 4, 4, 1, 9]).rolling_rank(3, method="average")
        shape: (5,)
        Series: '' [f64]
        [
            null
            null
            2.5
            1.0
            3.0
        ]
        """

    @unstable()
    def rolling_skew(
        self,
        window_size: int,
        *,
        bias: bool = True,
        min_samples: int | None = None,
        center: bool = False,
    ) -> Series:
        """
        Compute a rolling skew.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        The window at a given row includes the row itself and the
        `window_size - 1` elements before it.

        Parameters
        ----------
        window_size
            Integer size of the rolling window.
        bias
            If False, the calculations are corrected for statistical bias.
        min_samples
            The number of values in the window that should be non-null before computing
            a result. If set to `None` (default), it will be set equal to `window_size`.
        center
            Set the labels at the center of the window.

        See Also
        --------
        Series.skew

        Examples
        --------
        >>> pl.Series([1, 4, 2, 9]).rolling_skew(3)
        shape: (4,)
        Series: '' [f64]
        [
            null
            null
            0.381802
            0.47033
        ]

        Note how the values match

        >>> pl.Series([1, 4, 2]).skew(), pl.Series([4, 2, 9]).skew()
        (0.38180177416060584, 0.47033046033698594)
        """

    @unstable()
    def rolling_kurtosis(
        self,
        window_size: int,
        *,
        fisher: bool = True,
        bias: bool = True,
        min_samples: int | None = None,
        center: bool = False,
    ) -> Series:
        """
        Compute a rolling kurtosis.

        .. warning::
            This functionality is considered **unstable**. It may be changed
            at any point without it being considered a breaking change.

        The window at a given row will include the row itself, and the `window_size - 1`
        elements before it.

        Parameters
        ----------
        window_size
            Integer size of the rolling window.
        fisher : bool, optional
            If True, Fisher's definition is used (normal ==> 0.0). If False,
            Pearson's definition is used (normal ==> 3.0).
        bias : bool, optional
            If False, the calculations are corrected for statistical bias.
        min_samples
            The number of values in the window that should be non-null before computing
            a result. If set to `None` (default), it will be set equal to `window_size`.
        center
            Set the labels at the center of the window.

        See Also
        --------
        Series.kurtosis

        Examples
        --------
        >>> pl.Series([1, 4, 2, 9]).rolling_kurtosis(3)
        shape: (4,)
        Series: '' [f64]
        [
            null
            null
            -1.5
            -1.5
        ]
        """

    def sample(
        self,
        n: int | None = None,
        *,
        fraction: float | None = None,
        with_replacement: bool = False,
        shuffle: bool = False,
        seed: int | None = None,
    ) -> Series:
        """
        Sample from this Series.

        Parameters
        ----------
        n
            Number of items to return. Cannot be used with `fraction`. Defaults to 1 if
            `fraction` is None.
        fraction
            Fraction of items to return. Cannot be used with `n`.
        with_replacement
            Allow values to be sampled more than once.
        shuffle
            Shuffle the order of sampled data points.
        seed
            Seed for the random number generator. If set to None (default), a
            random seed is generated for each sample operation.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3, 4, 5])
        >>> s.sample(2, seed=0)  # doctest: +IGNORE_RESULT
        shape: (2,)
        Series: 'a' [i64]
        [
            1
            5
        ]
        """

    def peak_max(self) -> Self:
        """
        Get a boolean mask of the local maximum peaks.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3, 4, 5])
        >>> s.peak_max()
        shape: (5,)
        Series: 'a' [bool]
        [
                false
                false
                false
                false
                true
        ]
        """

    def peak_min(self) -> Self:
        """
        Get a boolean mask of the local minimum peaks.

        Examples
        --------
        >>> s = pl.Series("a", [4, 1, 3, 2, 5])
        >>> s.peak_min()
        shape: (5,)
        Series: 'a' [bool]
        [
            false
            true
            false
            true
            false
        ]
        """

    def n_unique(self) -> int:
        """
        Count the number of unique values in this Series.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 2, 3])
        >>> s.n_unique()
        3
        """
        return self._s.n_unique()

    def shrink_to_fit(self, *, in_place: bool = False) -> Series:
        """
        Shrink Series memory usage.

        Shrinks the underlying array capacity to exactly fit the actual data.
        (Note that this function does not change the Series data type).
        """
        if in_place:
            self._s.shrink_to_fit()
            return self
        else:
            series = self.clone()
            series._s.shrink_to_fit()
            return series

    def hash(
        self,
        seed: int = 0,
        seed_1: int | None = None,
        seed_2: int | None = None,
        seed_3: int | None = None,
    ) -> Series:
        """
        Hash the Series.

        The hash value is of type `UInt64`.

        Parameters
        ----------
        seed
            Random seed parameter. Defaults to 0.
        seed_1
            Random seed parameter. Defaults to `seed` if not set.
        seed_2
            Random seed parameter. Defaults to `seed` if not set.
        seed_3
            Random seed parameter. Defaults to `seed` if not set.

        Notes
        -----
        This implementation of `hash` does not guarantee stable results
        across different Polars versions. Its stability is only guaranteed within a
        single version.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.hash(seed=42)  # doctest: +IGNORE_RESULT
        shape: (3,)
        Series: 'a' [u64]
        [
            10734580197236529959
            3022416320763508302
            13756996518000038261
        ]
        """

    def reinterpret(self, *, signed: bool = True) -> Series:
        """
        Reinterpret the underlying bits as a signed/unsigned integer.

        This operation is only allowed for 64bit integers. For lower bits integers,
        you can safely use that cast operation.

        Parameters
        ----------
        signed
            If True, reinterpret as `pl.Int64`. Otherwise, reinterpret as `pl.UInt64`.

        Examples
        --------
        >>> s = pl.Series("a", [-(2**60), -2, 3])
        >>> s
        shape: (3,)
        Series: 'a' [i64]
        [
                -1152921504606846976
                -2
                3
        ]
        >>> s.reinterpret(signed=False)
        shape: (3,)
        Series: 'a' [u64]
        [
                17293822569102704640
                18446744073709551614
                3
        ]
        """

    def interpolate(self, method: InterpolationMethod = "linear") -> Series:
        """
        Interpolate intermediate values.

        Nulls at the beginning and end of the series remain null.

        Parameters
        ----------
        method : {'linear', 'nearest'}
            Interpolation method.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, None, None, 5])
        >>> s.interpolate()
        shape: (5,)
        Series: 'a' [f64]
        [
            1.0
            2.0
            3.0
            4.0
            5.0
        ]
        """

    def interpolate_by(self, by: IntoExpr) -> Series:
        """
        Interpolate intermediate values with x-coordinate based on another column.

        Nulls at the beginning and end of the series remain null.

        Parameters
        ----------
        by
            Column to interpolate values based on.

        Examples
        --------
        Fill null values using linear interpolation.

        >>> s = pl.Series([1, None, None, 3])
        >>> by = pl.Series([1, 2, 7, 8])
        >>> s.interpolate_by(by)
        shape: (4,)
        Series: '' [f64]
        [
            1.0
            1.285714
            2.714286
            3.0
        ]
        """

    def abs(self) -> Series:
        """
        Compute absolute values.

        Same as `abs(series)`.

        Examples
        --------
        >>> s = pl.Series([1, -2, -3])
        >>> s.abs()
        shape: (3,)
        Series: '' [i64]
        [
            1
            2
            3
        ]
        """

    def rank(
        self,
        method: RankMethod = "average",
        *,
        descending: bool = False,
        seed: int | None = None,
    ) -> Series:
        """
        Assign ranks to data, dealing with ties appropriately.

        Parameters
        ----------
        method : {'average', 'min', 'max', 'dense', 'ordinal', 'random'}
            The method used to assign ranks to tied elements.
            The following methods are available (default is 'average'):

            - 'average' : The average of the ranks that would have been assigned to
              all the tied values is assigned to each value.
            - 'min' : The minimum of the ranks that would have been assigned to all
              the tied values is assigned to each value. (This is also referred to
              as "competition" ranking.)
            - 'max' : The maximum of the ranks that would have been assigned to all
              the tied values is assigned to each value.
            - 'dense' : Like 'min', but the rank of the next highest element is
              assigned the rank immediately after those assigned to the tied
              elements.
            - 'ordinal' : All values are given a distinct rank, corresponding to
              the order that the values occur in the Series.
            - 'random' : Like 'ordinal', but the rank for ties is not dependent
              on the order that the values occur in the Series.
        descending
            Rank in descending order.
        seed
            If `method="random"`, use this as seed.

        Examples
        --------
        The 'average' method:

        >>> s = pl.Series("a", [3, 6, 1, 1, 6])
        >>> s.rank()
        shape: (5,)
        Series: 'a' [f64]
        [
            3.0
            4.5
            1.5
            1.5
            4.5
        ]

        The 'ordinal' method:

        >>> s = pl.Series("a", [3, 6, 1, 1, 6])
        >>> s.rank("ordinal")
        shape: (5,)
        Series: 'a' [u32]
        [
            3
            4
            1
            2
            5
        ]
        """

    def diff(self, n: int = 1, null_behavior: NullBehavior = "ignore") -> Series:
        """
        Calculate the first discrete difference between shifted items.

        Parameters
        ----------
        n
            Number of slots to shift.
        null_behavior : {'ignore', 'drop'}
            How to handle null values.

        Examples
        --------
        >>> s = pl.Series("s", values=[20, 10, 30, 25, 35], dtype=pl.Int8)
        >>> s.diff()
        shape: (5,)
        Series: 's' [i8]
        [
            null
            -10
            20
            -5
            10
        ]

        >>> s.diff(n=2)
        shape: (5,)
        Series: 's' [i8]
        [
            null
            null
            10
            15
            5
        ]

        >>> s.diff(n=2, null_behavior="drop")
        shape: (3,)
        Series: 's' [i8]
        [
            10
            15
            5
        ]
        """

    def pct_change(self, n: int | IntoExprColumn = 1) -> Series:
        """
        Computes percentage change between values.

        Percentage change (as fraction) between current element and most-recent
        non-null element at least `n` period(s) before the current element.

        Computes the change from the previous row by default.

        Parameters
        ----------
        n
            periods to shift for forming percent change.

        Notes
        -----
        Null values are preserved. If you're coming from pandas, this matches
        their ``fill_method=None`` behaviour.

        Examples
        --------
        >>> pl.Series(range(10)).pct_change()
        shape: (10,)
        Series: '' [f64]
        [
            null
            inf
            1.0
            0.5
            0.333333
            0.25
            0.2
            0.166667
            0.142857
            0.125
        ]

        >>> pl.Series([1, 2, 4, 8, 16, 32, 64, 128, 256, 512]).pct_change(2)
        shape: (10,)
        Series: '' [f64]
        [
            null
            null
            3.0
            3.0
            3.0
            3.0
            3.0
            3.0
            3.0
            3.0
        ]
        """

    def skew(self, *, bias: bool = True) -> float | None:
        r"""
        Compute the sample skewness of a data set.

        For normally distributed data, the skewness should be about zero. For
        unimodal continuous distributions, a skewness value greater than zero means
        that there is more weight in the right tail of the distribution. The
        function `skewtest` can be used to determine if the skewness value
        is close enough to zero, statistically speaking.


        See scipy.stats for more information.

        Parameters
        ----------
        bias : bool, optional
            If False, the calculations are corrected for statistical bias.

        Notes
        -----
        The sample skewness is computed as the Fisher-Pearson coefficient
        of skewness, i.e.

        .. math:: g_1=\frac{m_3}{m_2^{3/2}}

        where

        .. math:: m_i=\frac{1}{N}\sum_{n=1}^N(x[n]-\bar{x})^i

        is the biased sample :math:`i\texttt{th}` central moment, and
        :math:`\bar{x}` is
        the sample mean. If `bias` is False, the calculations are
        corrected for bias and the value computed is the adjusted
        Fisher-Pearson standardized moment coefficient, i.e.

        .. math::
            G_1 = \frac{k_3}{k_2^{3/2}} = \frac{\sqrt{N(N-1)}}{N-2}\frac{m_3}{m_2^{3/2}}

        Examples
        --------
        >>> s = pl.Series([1, 2, 2, 4, 5])
        >>> s.skew()
        0.34776706224699483
        """
        return self._s.skew(bias)

    def kurtosis(self, *, fisher: bool = True, bias: bool = True) -> float | None:
        """
        Compute the kurtosis (Fisher or Pearson) of a dataset.

        Kurtosis is the fourth central moment divided by the square of the
        variance. If Fisher's definition is used, then 3.0 is subtracted from
        the result to give 0.0 for a normal distribution.
        If bias is False then the kurtosis is calculated using k statistics to
        eliminate bias coming from biased moment estimators

        See scipy.stats for more information

        Parameters
        ----------
        fisher : bool, optional
            If True, Fisher's definition is used (normal ==> 0.0). If False,
            Pearson's definition is used (normal ==> 3.0).
        bias : bool, optional
            If False, the calculations are corrected for statistical bias.

        Examples
        --------
        >>> s = pl.Series("grades", [66, 79, 54, 97, 96, 70, 69, 85, 93, 75])
        >>> s.kurtosis()
        -1.0522623626787952
        >>> s.kurtosis(fisher=False)
        1.9477376373212048
        >>> s.kurtosis(fisher=False, bias=False)
        2.1040361802642717
        """
        return self._s.kurtosis(fisher, bias)

    def clip(
        self,
        lower_bound: NumericLiteral | TemporalLiteral | IntoExprColumn | None = None,
        upper_bound: NumericLiteral | TemporalLiteral | IntoExprColumn | None = None,
    ) -> Series:
        """
        Set values outside the given boundaries to the boundary value.

        Parameters
        ----------
        lower_bound
            Lower bound. Accepts expression input.
            Non-expression inputs are parsed as literals.
            If set to `None` (default), no lower bound is applied.
        upper_bound
            Upper bound. Accepts expression input.
            Non-expression inputs are parsed as literals.
            If set to `None` (default), no upper bound is applied.

        See Also
        --------
        when

        Notes
        -----
        This method only works for numeric and temporal columns. To clip other data
        types, consider writing a `when-then-otherwise` expression. See :func:`when`.

        Examples
        --------
        Specifying both a lower and upper bound:

        >>> s = pl.Series([-50, 5, 50, None])
        >>> s.clip(1, 10)
        shape: (4,)
        Series: '' [i64]
        [
                1
                5
                10
                null
        ]

        Specifying only a single bound:

        >>> s.clip(upper_bound=10)
        shape: (4,)
        Series: '' [i64]
        [
                -50
                5
                10
                null
        ]
        """

    def lower_bound(self) -> Self:
        """
        Return the lower bound of this Series' dtype as a unit Series.

        See Also
        --------
        upper_bound : return the upper bound of the given Series' dtype.

        Examples
        --------
        >>> s = pl.Series("s", [-1, 0, 1], dtype=pl.Int32)
        >>> s.lower_bound()
        shape: (1,)
        Series: 's' [i32]
        [
            -2147483648
        ]

        >>> s = pl.Series("s", [1.0, 2.5, 3.0], dtype=pl.Float32)
        >>> s.lower_bound()
        shape: (1,)
        Series: 's' [f32]
        [
            -inf
        ]
        """

    def upper_bound(self) -> Self:
        """
        Return the upper bound of this Series' dtype as a unit Series.

        See Also
        --------
        lower_bound : return the lower bound of the given Series' dtype.

        Examples
        --------
        >>> s = pl.Series("s", [-1, 0, 1], dtype=pl.Int8)
        >>> s.upper_bound()
        shape: (1,)
        Series: 's' [i8]
        [
            127
        ]

        >>> s = pl.Series("s", [1.0, 2.5, 3.0], dtype=pl.Float64)
        >>> s.upper_bound()
        shape: (1,)
        Series: 's' [f64]
        [
            inf
        ]
        """

    def replace(
        self,
        old: IntoExpr | Sequence[Any] | Mapping[Any, Any],
        new: IntoExpr | Sequence[Any] | NoDefault = no_default,
        *,
        default: IntoExpr | NoDefault = no_default,
        return_dtype: PolarsDataType | None = None,
    ) -> Self:
        """
        Replace values by different values of the same data type.

        Parameters
        ----------
        old
            Value or sequence of values to replace.
            Also accepts a mapping of values to their replacement as syntactic sugar for
            `replace(old=Series(mapping.keys()), new=Series(mapping.values()))`.
        new
            Value or sequence of values to replace by.
            Length must match the length of `old` or have length 1.

        default
            Set values that were not replaced to this value.
            Defaults to keeping the original value.
            Accepts expression input. Non-expression inputs are parsed as literals.

            .. deprecated:: 0.20.31
                Use :meth:`replace_all` instead to set a default while replacing values.

        return_dtype
            The data type of the resulting expression. If set to `None` (default),
            the data type is determined automatically based on the other inputs.

            .. deprecated:: 0.20.31
                Use :meth:`replace_all` instead to set a return data type while
                replacing values.


        See Also
        --------
        replace_strict
        str.replace

        Notes
        -----
        The global string cache must be enabled when replacing categorical values.

        Examples
        --------
        Replace a single value by another value. Values that were not replaced remain
        unchanged.

        >>> s = pl.Series([1, 2, 2, 3])
        >>> s.replace(2, 100)
        shape: (4,)
        Series: '' [i64]
        [
                1
                100
                100
                3
        ]

        Replace multiple values by passing sequences to the `old` and `new` parameters.

        >>> s.replace([2, 3], [100, 200])
        shape: (4,)
        Series: '' [i64]
        [
                1
                100
                100
                200
        ]

        Passing a mapping with replacements is also supported as syntactic sugar.

        >>> mapping = {2: 100, 3: 200}
        >>> s.replace(mapping)
        shape: (4,)
        Series: '' [i64]
        [
                1
                100
                100
                200
        ]

        The original data type is preserved when replacing by values of a different
        data type. Use :meth:`replace_strict` to replace and change the return data
        type.

        >>> s = pl.Series(["x", "y", "z"])
        >>> mapping = {"x": 1, "y": 2, "z": 3}
        >>> s.replace(mapping)
        shape: (3,)
        Series: '' [str]
        [
                "1"
                "2"
                "3"
        ]
        """

    def replace_strict(
        self,
        old: IntoExpr | Sequence[Any] | Mapping[Any, Any],
        new: IntoExpr | Sequence[Any] | NoDefault = no_default,
        *,
        default: IntoExpr | NoDefault = no_default,
        return_dtype: PolarsDataType | None = None,
    ) -> Self:
        """
        Replace all values by different values.

        Parameters
        ----------
        old
            Value or sequence of values to replace.
            Also accepts a mapping of values to their replacement as syntactic sugar for
            `replace_strict(old=Series(mapping.keys()), new=Series(mapping.values()))`.
        new
            Value or sequence of values to replace by.
            Length must match the length of `old` or have length 1.
        default
            Set values that were not replaced to this value. If no default is specified,
            (default), an error is raised if any values were not replaced.
            Accepts expression input. Non-expression inputs are parsed as literals.
        return_dtype
            The data type of the resulting Series. If set to `None` (default),
            the data type is determined automatically based on the other inputs.

        Raises
        ------
        InvalidOperationError
            If any non-null values in the original column were not replaced, and no
            `default` was specified.

        See Also
        --------
        replace
        str.replace

        Notes
        -----
        The global string cache must be enabled when replacing categorical values.

        Examples
        --------
        Replace values by passing sequences to the `old` and `new` parameters.

        >>> s = pl.Series([1, 2, 2, 3])
        >>> s.replace_strict([1, 2, 3], [100, 200, 300])
        shape: (4,)
        Series: '' [i64]
        [
                100
                200
                200
                300
        ]

        Passing a mapping with replacements is also supported as syntactic sugar.

        >>> mapping = {1: 100, 2: 200, 3: 300}
        >>> s.replace_strict(mapping)
        shape: (4,)
        Series: '' [i64]
        [
                100
                200
                200
                300
        ]

        By default, an error is raised if any non-null values were not replaced.
        Specify a default to set all values that were not matched.

        >>> mapping = {2: 200, 3: 300}
        >>> s.replace_strict(mapping)  # doctest: +SKIP
        Traceback (most recent call last):
        ...
        polars.exceptions.InvalidOperationError: incomplete mapping specified for `replace_strict`
        >>> s.replace_strict(mapping, default=-1)
        shape: (4,)
        Series: '' [i64]
        [
                -1
                200
                200
                300
        ]

        The default can be another Series.

        >>> default = pl.Series([2.5, 5.0, 7.5, 10.0])
        >>> s.replace_strict(2, 200, default=default)
        shape: (4,)
        Series: '' [f64]
        [
                2.5
                200.0
                200.0
                10.0
        ]

        Replacing by values of a different data type sets the return type based on
        a combination of the `new` data type and the `default` data type.

        >>> s = pl.Series(["x", "y", "z"])
        >>> mapping = {"x": 1, "y": 2, "z": 3}
        >>> s.replace_strict(mapping)
        shape: (3,)
        Series: '' [i64]
        [
                1
                2
                3
        ]
        >>> s.replace_strict(mapping, default="x")
        shape: (3,)
        Series: '' [str]
        [
                "1"
                "2"
                "3"
        ]

        Set the `return_dtype` parameter to control the resulting data type directly.

        >>> s.replace_strict(mapping, return_dtype=pl.UInt8)
        shape: (3,)
        Series: '' [u8]
        [
                1
                2
                3
        ]
        """  # noqa: W505

    def reshape(self, dimensions: tuple[int, ...]) -> Series:
        """
        Reshape this Series to a flat Series or an Array Series.

        Parameters
        ----------
        dimensions
            Tuple of the dimension sizes. If a -1 is used in any of the dimensions, that
            dimension is inferred.

        Returns
        -------
        Series
            If a single dimension is given, results in a Series of the original
            data type.
            If a multiple dimensions are given, results in a Series of data type
            :class:`Array` with shape `dimensions`.

        See Also
        --------
        Series.list.explode : Explode a list column.

        Examples
        --------
        >>> s = pl.Series("foo", [1, 2, 3, 4, 5, 6, 7, 8, 9])
        >>> square = s.reshape((3, 3))
        >>> square
        shape: (3,)
        Series: 'foo' [array[i64, 3]]
        [
                [1, 2, 3]
                [4, 5, 6]
                [7, 8, 9]
        ]
        >>> square.reshape((9,))
        shape: (9,)
        Series: 'foo' [i64]
        [
                1
                2
                3
                4
                5
                6
                7
                8
                9
        ]
        """
        return self._from_pyseries(self._s.reshape(dimensions))

    def shuffle(self, seed: int | None = None) -> Series:
        """
        Shuffle the contents of this Series.

        Parameters
        ----------
        seed
            Seed for the random number generator. If set to None (default), a
            random seed is generated each time the shuffle is called.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.shuffle(seed=1)
        shape: (3,)
        Series: 'a' [i64]
        [
                2
                3
                1
        ]
        """

    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def ewm_mean(
        self,
        *,
        com: float | None = None,
        span: float | None = None,
        half_life: float | None = None,
        alpha: float | None = None,
        adjust: bool = True,
        min_samples: int = 1,
        ignore_nulls: bool = False,
    ) -> Series:
        r"""
        Compute exponentially-weighted moving average.

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        com
            Specify decay in terms of center of mass, :math:`\gamma`, with

                .. math::
                    \alpha = \frac{1}{1 + \gamma} \; \forall \; \gamma \geq 0
        span
            Specify decay in terms of span, :math:`\theta`, with

                .. math::
                    \alpha = \frac{2}{\theta + 1} \; \forall \; \theta \geq 1
        half_life
            Specify decay in terms of half-life, :math:`\tau`, with

                .. math::
                    \alpha = 1 - \exp \left\{ \frac{ -\ln(2) }{ \tau } \right\} \;
                    \forall \; \tau > 0
        alpha
            Specify smoothing factor alpha directly, :math:`0 < \alpha \leq 1`.
        adjust
            Divide by decaying adjustment factor in beginning periods to account for
            imbalance in relative weightings

                - When `adjust=True` (the default) the EW function is calculated
                  using weights :math:`w_i = (1 - \alpha)^i`
                - When `adjust=False` the EW function is calculated
                  recursively by

                  .. math::
                    y_0 &= x_0 \\
                    y_t &= (1 - \alpha)y_{t - 1} + \alpha x_t
        min_samples
            Minimum number of observations in window required to have a value
            (otherwise result is null).
        ignore_nulls
            Ignore missing values when calculating weights.

                - When `ignore_nulls=False` (default), weights are based on absolute
                  positions.
                  For example, the weights of :math:`x_0` and :math:`x_2` used in
                  calculating the final weighted average of
                  [:math:`x_0`, None, :math:`x_2`] are
                  :math:`(1-\alpha)^2` and :math:`1` if `adjust=True`, and
                  :math:`(1-\alpha)^2` and :math:`\alpha` if `adjust=False`.

                - When `ignore_nulls=True`, weights are based
                  on relative positions. For example, the weights of
                  :math:`x_0` and :math:`x_2` used in calculating the final weighted
                  average of [:math:`x_0`, None, :math:`x_2`] are
                  :math:`1-\alpha` and :math:`1` if `adjust=True`,
                  and :math:`1-\alpha` and :math:`\alpha` if `adjust=False`.

        Examples
        --------
        >>> s = pl.Series([1, 2, 3])
        >>> s.ewm_mean(com=1, ignore_nulls=False)
        shape: (3,)
        Series: '' [f64]
        [
                1.0
                1.666667
                2.428571
        ]
        """

    def ewm_mean_by(
        self,
        by: IntoExpr,
        *,
        half_life: str | timedelta,
    ) -> Series:
        r"""
        Compute time-based exponentially weighted moving average.

        Given observations :math:`x_0, x_1, \ldots, x_{n-1}` at times
        :math:`t_0, t_1, \ldots, t_{n-1}`, the EWMA is calculated as

            .. math::

                y_0 &= x_0

                \alpha_i &= 1 - \exp \left\{ \frac{ -\ln(2)(t_i-t_{i-1}) }
                    { \tau } \right\}

                y_i &= \alpha_i x_i + (1 - \alpha_i) y_{i-1}; \quad i > 0

        where :math:`\tau` is the `half_life`.

        Parameters
        ----------
        by
            Times to calculate average by. Should be ``DateTime``, ``Date``, ``UInt64``,
            ``UInt32``, ``Int64``, or ``Int32`` data type.
        half_life
            Unit over which observation decays to half its value.

            Can be created either from a timedelta, or
            by using the following string language:

            - 1ns   (1 nanosecond)
            - 1us   (1 microsecond)
            - 1ms   (1 millisecond)
            - 1s    (1 second)
            - 1m    (1 minute)
            - 1h    (1 hour)
            - 1d    (1 day)
            - 1w    (1 week)
            - 1i    (1 index count)

            Or combine them:
            "3d12h4m25s" # 3 days, 12 hours, 4 minutes, and 25 seconds

            Note that `half_life` is treated as a constant duration - calendar
            durations such as months (or even days in the time-zone-aware case)
            are not supported, please express your duration in an approximately
            equivalent number of hours (e.g. '370h' instead of '1mo').

        Returns
        -------
        Expr
            Float32 if input is Float32, otherwise Float64.

        Examples
        --------
        >>> from datetime import date, timedelta
        >>> df = pl.DataFrame(
        ...     {
        ...         "values": [0, 1, 2, None, 4],
        ...         "times": [
        ...             date(2020, 1, 1),
        ...             date(2020, 1, 3),
        ...             date(2020, 1, 10),
        ...             date(2020, 1, 15),
        ...             date(2020, 1, 17),
        ...         ],
        ...     }
        ... ).sort("times")
        >>> df["values"].ewm_mean_by(df["times"], half_life="4d")
        shape: (5,)
        Series: 'values' [f64]
        [
                0.0
                0.292893
                1.492474
                null
                3.254508
        ]
        """

    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def ewm_std(
        self,
        *,
        com: float | None = None,
        span: float | None = None,
        half_life: float | None = None,
        alpha: float | None = None,
        adjust: bool = True,
        bias: bool = False,
        min_samples: int = 1,
        ignore_nulls: bool = False,
    ) -> Series:
        r"""
        Compute exponentially-weighted moving standard deviation.

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        com
            Specify decay in terms of center of mass, :math:`\gamma`, with

                .. math::
                    \alpha = \frac{1}{1 + \gamma} \; \forall \; \gamma \geq 0
        span
            Specify decay in terms of span, :math:`\theta`, with

                .. math::
                    \alpha = \frac{2}{\theta + 1} \; \forall \; \theta \geq 1
        half_life
            Specify decay in terms of half-life, :math:`\lambda`, with

                .. math::
                    \alpha = 1 - \exp \left\{ \frac{ -\ln(2) }{ \lambda } \right\} \;
                    \forall \; \lambda > 0
        alpha
            Specify smoothing factor alpha directly, :math:`0 < \alpha \leq 1`.
        adjust
            Divide by decaying adjustment factor in beginning periods to account for
            imbalance in relative weightings

                - When `adjust=True` (the default) the EW function is calculated
                  using weights :math:`w_i = (1 - \alpha)^i`
                - When `adjust=False` the EW function is calculated
                  recursively by

                  .. math::
                    y_0 &= x_0 \\
                    y_t &= (1 - \alpha)y_{t - 1} + \alpha x_t
        bias
            When `bias=False`, apply a correction to make the estimate statistically
            unbiased.
        min_samples
            Minimum number of observations in window required to have a value
            (otherwise result is null).
        ignore_nulls
            Ignore missing values when calculating weights.

                - When `ignore_nulls=False` (default), weights are based on absolute
                  positions.
                  For example, the weights of :math:`x_0` and :math:`x_2` used in
                  calculating the final weighted average of
                  [:math:`x_0`, None, :math:`x_2`] are
                  :math:`(1-\alpha)^2` and :math:`1` if `adjust=True`, and
                  :math:`(1-\alpha)^2` and :math:`\alpha` if `adjust=False`.

                - When `ignore_nulls=True`, weights are based
                  on relative positions. For example, the weights of
                  :math:`x_0` and :math:`x_2` used in calculating the final weighted
                  average of [:math:`x_0`, None, :math:`x_2`] are
                  :math:`1-\alpha` and :math:`1` if `adjust=True`,
                  and :math:`1-\alpha` and :math:`\alpha` if `adjust=False`.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.ewm_std(com=1, ignore_nulls=False)
        shape: (3,)
        Series: 'a' [f64]
        [
            0.0
            0.707107
            0.963624
        ]
        """

    @deprecate_renamed_parameter("min_periods", "min_samples", version="1.21.0")
    def ewm_var(
        self,
        *,
        com: float | None = None,
        span: float | None = None,
        half_life: float | None = None,
        alpha: float | None = None,
        adjust: bool = True,
        bias: bool = False,
        min_samples: int = 1,
        ignore_nulls: bool = False,
    ) -> Series:
        r"""
        Compute exponentially-weighted moving variance.

        .. versionchanged:: 1.21.0
            The `min_periods` parameter was renamed `min_samples`.

        Parameters
        ----------
        com
            Specify decay in terms of center of mass, :math:`\gamma`, with

                .. math::
                    \alpha = \frac{1}{1 + \gamma} \; \forall \; \gamma \geq 0
        span
            Specify decay in terms of span, :math:`\theta`, with

                .. math::
                    \alpha = \frac{2}{\theta + 1} \; \forall \; \theta \geq 1
        half_life
            Specify decay in terms of half-life, :math:`\lambda`, with

                .. math::
                    \alpha = 1 - \exp \left\{ \frac{ -\ln(2) }{ \lambda } \right\} \;
                    \forall \; \lambda > 0
        alpha
            Specify smoothing factor alpha directly, :math:`0 < \alpha \leq 1`.
        adjust
            Divide by decaying adjustment factor in beginning periods to account for
            imbalance in relative weightings

                - When `adjust=True` (the default) the EW function is calculated
                  using weights :math:`w_i = (1 - \alpha)^i`
                - When `adjust=False` the EW function is calculated
                  recursively by

                  .. math::
                    y_0 &= x_0 \\
                    y_t &= (1 - \alpha)y_{t - 1} + \alpha x_t
        bias
            When `bias=False`, apply a correction to make the estimate statistically
            unbiased.
        min_samples
            Minimum number of observations in window required to have a value
            (otherwise result is null).
        ignore_nulls
            Ignore missing values when calculating weights.

                - When `ignore_nulls=False` (default), weights are based on absolute
                  positions.
                  For example, the weights of :math:`x_0` and :math:`x_2` used in
                  calculating the final weighted average of
                  [:math:`x_0`, None, :math:`x_2`] are
                  :math:`(1-\alpha)^2` and :math:`1` if `adjust=True`, and
                  :math:`(1-\alpha)^2` and :math:`\alpha` if `adjust=False`.

                - When `ignore_nulls=True`, weights are based
                  on relative positions. For example, the weights of
                  :math:`x_0` and :math:`x_2` used in calculating the final weighted
                  average of [:math:`x_0`, None, :math:`x_2`] are
                  :math:`1-\alpha` and :math:`1` if `adjust=True`,
                  and :math:`1-\alpha` and :math:`\alpha` if `adjust=False`.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.ewm_var(com=1, ignore_nulls=False)
        shape: (3,)
        Series: 'a' [f64]
        [
            0.0
            0.5
            0.928571
        ]
        """

    def extend_constant(self, value: IntoExpr, n: int | IntoExprColumn) -> Series:
        """
        Extremely fast method for extending the Series with 'n' copies of a value.

        Parameters
        ----------
        value
            A constant literal value or a unit expression with which to extend the
            expression result Series; can pass None to extend with nulls.
        n
            The number of additional values that will be added.

        Examples
        --------
        >>> s = pl.Series([1, 2, 3])
        >>> s.extend_constant(99, n=2)
        shape: (5,)
        Series: '' [i64]
        [
                1
                2
                3
                99
                99
        ]
        """

    def set_sorted(self, *, descending: bool = False) -> Self:
        """
        Flags the Series as 'sorted'.

        Enables downstream code to user fast paths for sorted arrays.

        Parameters
        ----------
        descending
            If the `Series` order is descending.

        Warnings
        --------
        This can lead to incorrect results if this `Series` is not sorted!!
        Use with care!

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.set_sorted().max()
        3
        """
        return self._from_pyseries(self._s.set_sorted_flag(descending))

    def new_from_index(self, index: int, length: int) -> Self:
        """
        Create a new Series filled with values from the given index.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3, 4, 5])
        >>> s.new_from_index(1, 3)
        shape: (3,)
        Series: 'a' [i64]
        [
            2
            2
            2
        ]
        """
        return self._from_pyseries(self._s.new_from_index(index, length))

    def shrink_dtype(self) -> Series:
        """
        Shrink numeric columns to the minimal required datatype.

        Shrink to the dtype needed to fit the extrema of this [`Series`].
        This can be used to reduce memory pressure.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3, 4, 5, 6])
        >>> s
        shape: (6,)
        Series: 'a' [i64]
        [
                1
                2
                3
                4
                5
                6
        ]
        >>> s.shrink_dtype()
        shape: (6,)
        Series: 'a' [i8]
        [
                1
                2
                3
                4
                5
                6
        ]
        """
        return wrap_s(self._s.shrink_dtype())

    def get_chunks(self) -> list[Series]:
        """
        Get the chunks of this Series as a list of Series.

        Examples
        --------
        >>> s1 = pl.Series("a", [1, 2, 3])
        >>> s2 = pl.Series("a", [4, 5, 6])
        >>> s = pl.concat([s1, s2], rechunk=False)
        >>> s.get_chunks()
        [shape: (3,)
        Series: 'a' [i64]
        [
                1
                2
                3
        ], shape: (3,)
        Series: 'a' [i64]
        [
                4
                5
                6
        ]]
        """
        return self._s.get_chunks()

    def implode(self) -> Self:
        """
        Aggregate values into a list.

        The returned list itself is a scalar value of `list` dtype.

        Examples
        --------
        >>> s = pl.Series("a", [1, 2, 3])
        >>> s.implode()
        shape: (1,)
        Series: 'a' [list[i64]]
        [
            [1, 2, 3]
        ]
        """

    def bitwise_count_ones(self) -> Self:
        """Evaluate the number of set bits."""

    def bitwise_count_zeros(self) -> Self:
        """Evaluate the number of unset Self."""

    def bitwise_leading_ones(self) -> Self:
        """Evaluate the number most-significant set bits before seeing an unset bit."""

    def bitwise_leading_zeros(self) -> Self:
        """Evaluate the number most-significant unset bits before seeing a set bit."""

    def bitwise_trailing_ones(self) -> Self:
        """Evaluate the number least-significant set bits before seeing an unset bit."""

    def bitwise_trailing_zeros(self) -> Self:
        """Evaluate the number least-significant unset bits before seeing a set bit."""

    def bitwise_and(self) -> PythonLiteral | None:
        """Perform an aggregation of bitwise ANDs."""
        return self._s.bitwise_and()

    def bitwise_or(self) -> PythonLiteral | None:
        """Perform an aggregation of bitwise ORs."""
        return self._s.bitwise_or()

    def bitwise_xor(self) -> PythonLiteral | None:
        """Perform an aggregation of bitwise XORs."""
        return self._s.bitwise_xor()

    def first(self) -> PythonLiteral | None:
        """
        Get the first element of the Series.

        Returns `None` if the Series is empty.
        """
        return self._s.first()

    def last(self) -> PythonLiteral | None:
        """
        Get the last element of the Series.

        Returns `None` if the Series is empty.
        """
        return self._s.last()

    def approx_n_unique(self) -> PythonLiteral | None:
        """
        Approximate count of unique values.

        This is done using the HyperLogLog++ algorithm for cardinality estimation.
        """
        return self._s.approx_n_unique()

    def _row_encode(
        self,
        *,
        unordered: bool = False,
        descending: bool | None = None,
        nulls_last: bool | None = None,
    ) -> Series:
        """Encode to the row encoding."""
        return (
            self.to_frame()
            .select_seq(
                F.col(self.name)._row_encode(
                    unordered=unordered, descending=descending, nulls_last=nulls_last
                )
            )
            .to_series()
        )

    def _row_decode(
        self,
        names: Sequence[str],
        dtypes: Sequence[PolarsDataType],
        *,
        unordered: bool = False,
        descending: Sequence[bool] | None = None,
        nulls_last: Sequence[bool] | None = None,
    ) -> Series:
        """Decode from the row encoding."""
        return (
            self.to_frame()
            .select_seq(
                F.col(self.name)._row_decode(
                    names,
                    dtypes,
                    unordered=unordered,
                    descending=descending,
                    nulls_last=nulls_last,
                )
            )
            .to_series()
        )

    def repeat_by(self, by: int | IntoExprColumn) -> Self:
        """
        Repeat the elements in this Series as specified in the given expression.

        The repeated elements are expanded into a List.

        Parameters
        ----------
        by
            Numeric column that determines how often the values will be repeated.
            The column will be coerced to UInt32. Give this dtype to make the coercion
            a no-op.

        Returns
        -------
        Expr
            Expression of data type List, where the inner data type is equal to the
            original data type.
        """

    # Keep the `list` and `str` properties below at the end of the definition of Series,
    # as to not confuse mypy with the type annotation `str` and `list`

    @property
    def bin(self) -> BinaryNameSpace:
        """Create an object namespace of all binary related methods."""
        return BinaryNameSpace(self)

    @property
    def cat(self) -> CatNameSpace:
        """Create an object namespace of all categorical related methods."""
        return CatNameSpace(self)

    @property
    def dt(self) -> DateTimeNameSpace:
        """Create an object namespace of all datetime related methods."""
        return DateTimeNameSpace(self)

    @property
    def list(self) -> ListNameSpace:
        """Create an object namespace of all list related methods."""
        return ListNameSpace(self)

    @property
    def arr(self) -> ArrayNameSpace:
        """Create an object namespace of all array related methods."""
        return ArrayNameSpace(self)

    @property
    def str(self) -> StringNameSpace:
        """Create an object namespace of all string related methods."""
        return StringNameSpace(self)

    @property
    def struct(self) -> StructNameSpace:
        """Create an object namespace of all struct related methods."""
        return StructNameSpace(self)

    @property
    @unstable()
    def plot(self) -> SeriesPlot:
        """
        Create a plot namespace.

        .. warning::
            This functionality is currently considered **unstable**. It may be
            changed at any point without it being considered a breaking change.

        .. versionchanged:: 1.6.0
            In prior versions of Polars, HvPlot was the plotting backend. If you would
            like to restore the previous plotting functionality, all you need to do
            is add `import hvplot.polars` at the top of your script and replace
            `df.plot` with `df.hvplot`.

        Polars does not implement plotting logic itself, but instead defers to
        Altair:

        - `s.plot.hist(**kwargs)`
          is shorthand for
          `alt.Chart(s.to_frame()).mark_bar(tooltip=True).encode(x=alt.X(f'{s.name}:Q', bin=True), y='count()', **kwargs).interactive()`
        - `s.plot.kde(**kwargs)`
          is shorthand for
          `alt.Chart(s.to_frame()).transform_density(s.name, as_=[s.name, 'density']).mark_area(tooltip=True).encode(x=s.name, y='density:Q', **kwargs).interactive()`
        - for any other attribute `attr`, `s.plot.attr(**kwargs)`
          is shorthand for
          `alt.Chart(s.to_frame().with_row_index()).mark_attr(tooltip=True).encode(x='index', y=s.name, **kwargs).interactive()`

        For configuration, we suggest reading
        `Chart Configuration <https://altair-viz.github.io/altair-tutorial/notebooks/08-Configuration.html>`_.
        For example, you can:

        - Change the width/height/title with ``.properties(width=500, height=350, title="My amazing plot")``.
        - Change the x-axis label rotation with ``.configure_axisX(labelAngle=30)``.
        - Change the opacity of the points in your scatter plot with ``.configure_point(opacity=.5)``.

        Examples
        --------
        Histogram:

        >>> s = pl.Series([1, 4, 4, 6, 2, 4, 3, 5, 5, 7, 1])
        >>> s.plot.hist()  # doctest: +SKIP

        KDE plot:

        >>> s.plot.kde()  # doctest: +SKIP

        Line plot:

        >>> s.plot.line()  # doctest: +SKIP
        """  # noqa: W505
        if not _ALTAIR_AVAILABLE or parse_version(altair.__version__) < (5, 4, 0):
            msg = "altair>=5.4.0 is required for `.plot`"
            raise ModuleUpgradeRequiredError(msg)
        return SeriesPlot(self)


def _resolve_temporal_dtype(
    dtype: PolarsDataType | None,
    ndtype: np.dtype[np.datetime64] | np.dtype[np.timedelta64],
) -> PolarsDataType | None:
    """Given polars/numpy temporal dtypes, resolve to an explicit unit."""
    PolarsType = Duration if ndtype.type == np.timedelta64 else Datetime
    if dtype is None or (dtype == Datetime and not getattr(dtype, "time_unit", None)):
        time_unit = getattr(dtype, "time_unit", None) or np.datetime_data(ndtype)[0]
        # explicit formulation is verbose, but keeps mypy happy
        # (and avoids unsupported timeunits such as "s")
        if time_unit == "ns":
            dtype = PolarsType("ns")
        elif time_unit == "us":
            dtype = PolarsType("us")
        elif time_unit == "ms":
            dtype = PolarsType("ms")
        elif time_unit == "D" and ndtype.type == np.datetime64:
            dtype = Date
    return dtype
