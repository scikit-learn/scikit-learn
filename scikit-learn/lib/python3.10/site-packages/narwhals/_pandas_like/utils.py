from __future__ import annotations

import functools
import operator
import re
from typing import TYPE_CHECKING, Any, Callable, Literal, TypeVar, cast

import pandas as pd

from narwhals._compliant import EagerSeriesNamespace
from narwhals._constants import (
    MS_PER_SECOND,
    NS_PER_MICROSECOND,
    NS_PER_MILLISECOND,
    NS_PER_SECOND,
    SECONDS_PER_DAY,
    US_PER_SECOND,
)
from narwhals._exceptions import issue_warning
from narwhals._utils import (
    Implementation,
    Version,
    _DeferredIterable,
    check_columns_exist,
    isinstance_or_issubclass,
    requires,
)
from narwhals.exceptions import ShapeError

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Mapping
    from types import ModuleType

    from pandas._typing import Dtype as PandasDtype
    from pandas.core.dtypes.dtypes import BaseMaskedDtype
    from typing_extensions import TypeAlias, TypeIs

    from narwhals._duration import IntervalUnit
    from narwhals._pandas_like.expr import PandasLikeExpr
    from narwhals._pandas_like.series import PandasLikeSeries
    from narwhals._pandas_like.typing import (
        NativeDataFrameT,
        NativeNDFrameT,
        NativeSeriesT,
    )
    from narwhals.dtypes import DType
    from narwhals.typing import DTypeBackend, IntoDType, TimeUnit, _1DArray

    ExprT = TypeVar("ExprT", bound=PandasLikeExpr)
    UnitCurrent: TypeAlias = TimeUnit
    UnitTarget: TypeAlias = TimeUnit
    BinOpBroadcast: TypeAlias = Callable[[Any, int], Any]
    IntoRhs: TypeAlias = int


PANDAS_LIKE_IMPLEMENTATION = {
    Implementation.PANDAS,
    Implementation.CUDF,
    Implementation.MODIN,
}
PD_DATETIME_RGX = r"""^
    datetime64\[
        (?P<time_unit>s|ms|us|ns)                 # Match time unit: s, ms, us, or ns
        (?:,                                      # Begin non-capturing group for optional timezone
            \s*                                   # Optional whitespace after comma
            (?P<time_zone>                        # Start named group for timezone
                [a-zA-Z\/]+                       # Match timezone name, e.g., UTC, America/New_York
                (?:[+-]\d{2}:\d{2})?              # Optional offset in format +HH:MM or -HH:MM
                |                                 # OR
                pytz\.FixedOffset\(\d+\)          # Match pytz.FixedOffset with integer offset in parentheses
            )                                     # End time_zone group
        )?                                        # End optional timezone group
    \]                                            # Closing bracket for datetime64
$"""
PATTERN_PD_DATETIME = re.compile(PD_DATETIME_RGX, re.VERBOSE)
PA_DATETIME_RGX = r"""^
    timestamp\[
        (?P<time_unit>s|ms|us|ns)                 # Match time unit: s, ms, us, or ns
        (?:,                                      # Begin non-capturing group for optional timezone
            \s?tz=                                # Match "tz=" prefix
            (?P<time_zone>                        # Start named group for timezone
                [a-zA-Z\/]*                       # Match timezone name (e.g., UTC, America/New_York)
                (?:                               # Begin optional non-capturing group for offset
                    [+-]\d{2}:\d{2}               # Match offset in format +HH:MM or -HH:MM
                )?                                # End optional offset group
            )                                     # End time_zone group
        )?                                        # End optional timezone group
    \]                                            # Closing bracket for timestamp
    \[pyarrow\]                                   # Literal string "[pyarrow]"
$"""
PATTERN_PA_DATETIME = re.compile(PA_DATETIME_RGX, re.VERBOSE)
PD_DURATION_RGX = r"""^
    timedelta64\[
        (?P<time_unit>s|ms|us|ns)                 # Match time unit: s, ms, us, or ns
    \]                                            # Closing bracket for timedelta64
$"""

PATTERN_PD_DURATION = re.compile(PD_DURATION_RGX, re.VERBOSE)
PA_DURATION_RGX = r"""^
    duration\[
        (?P<time_unit>s|ms|us|ns)                 # Match time unit: s, ms, us, or ns
    \]                                            # Closing bracket for duration
    \[pyarrow\]                                   # Literal string "[pyarrow]"
$"""
PATTERN_PA_DURATION = re.compile(PA_DURATION_RGX, re.VERBOSE)

NativeIntervalUnit: TypeAlias = Literal[
    "year",
    "quarter",
    "month",
    "week",
    "day",
    "hour",
    "minute",
    "second",
    "millisecond",
    "microsecond",
    "nanosecond",
]
ALIAS_DICT = {"d": "D", "m": "min"}
UNITS_DICT: Mapping[IntervalUnit, NativeIntervalUnit] = {
    "y": "year",
    "q": "quarter",
    "mo": "month",
    "d": "day",
    "h": "hour",
    "m": "minute",
    "s": "second",
    "ms": "millisecond",
    "us": "microsecond",
    "ns": "nanosecond",
}

PANDAS_VERSION = Implementation.PANDAS._backend_version()
"""Static backend version for `pandas`.

Always available if we reached here, due to a module-level import.
"""


def is_pandas_or_modin(implementation: Implementation) -> bool:
    return implementation in {Implementation.PANDAS, Implementation.MODIN}


def align_and_extract_native(
    lhs: PandasLikeSeries, rhs: PandasLikeSeries | object
) -> tuple[pd.Series[Any] | object, pd.Series[Any] | object]:
    """Validate RHS of binary operation.

    If the comparison isn't supported, return `NotImplemented` so that the
    "right-hand-side" operation (e.g. `__radd__`) can be tried.
    """
    from narwhals._pandas_like.series import PandasLikeSeries

    lhs_index = lhs.native.index

    if lhs._broadcast and isinstance(rhs, PandasLikeSeries) and not rhs._broadcast:
        return lhs.native.iloc[0], rhs.native

    if isinstance(rhs, PandasLikeSeries):
        if rhs._broadcast:
            return (lhs.native, rhs.native.iloc[0])
        if rhs.native.index is not lhs_index:
            return (
                lhs.native,
                set_index(rhs.native, lhs_index, implementation=rhs._implementation),
            )
        return (lhs.native, rhs.native)

    if isinstance(rhs, list):
        msg = "Expected Series or scalar, got list."
        raise TypeError(msg)
    # `rhs` must be scalar, so just leave it as-is
    return lhs.native, rhs


def set_index(
    obj: NativeNDFrameT, index: Any, *, implementation: Implementation
) -> NativeNDFrameT:
    """Wrapper around pandas' set_axis to set object index.

    We can set `copy` / `inplace` based on implementation/version.
    """
    if isinstance(index, implementation.to_native_namespace().Index) and (
        expected_len := len(index)
    ) != (actual_len := len(obj)):
        msg = f"Expected object of length {expected_len}, got length: {actual_len}"
        raise ShapeError(msg)
    if implementation is Implementation.CUDF:
        obj = obj.copy(deep=False)
        obj.index = index
        return obj
    if implementation is Implementation.PANDAS and (
        (1, 5) <= implementation._backend_version() < (3,)
    ):  # pragma: no cover
        return obj.set_axis(index, axis=0, copy=False)
    return obj.set_axis(index, axis=0)  # pragma: no cover


def rename(
    obj: NativeNDFrameT, *args: Any, implementation: Implementation, **kwargs: Any
) -> NativeNDFrameT:
    """Wrapper around pandas' rename so that we can set `copy` based on implementation/version."""
    if implementation is Implementation.PANDAS and (
        implementation._backend_version() >= (3,)
    ):  # pragma: no cover
        result = obj.rename(*args, **kwargs, inplace=False)
    else:
        result = obj.rename(*args, **kwargs, copy=False, inplace=False)
    return cast("NativeNDFrameT", result)  # type: ignore[redundant-cast]


@functools.lru_cache(maxsize=16)
def non_object_native_to_narwhals_dtype(native_dtype: Any, version: Version) -> DType:  # noqa: C901, PLR0912
    dtype = str(native_dtype)

    dtypes = version.dtypes
    if dtype in {"int64", "Int64", "Int64[pyarrow]", "int64[pyarrow]"}:
        return dtypes.Int64()
    if dtype in {"int32", "Int32", "Int32[pyarrow]", "int32[pyarrow]"}:
        return dtypes.Int32()
    if dtype in {"int16", "Int16", "Int16[pyarrow]", "int16[pyarrow]"}:
        return dtypes.Int16()
    if dtype in {"int8", "Int8", "Int8[pyarrow]", "int8[pyarrow]"}:
        return dtypes.Int8()
    if dtype in {"uint64", "UInt64", "UInt64[pyarrow]", "uint64[pyarrow]"}:
        return dtypes.UInt64()
    if dtype in {"uint32", "UInt32", "UInt32[pyarrow]", "uint32[pyarrow]"}:
        return dtypes.UInt32()
    if dtype in {"uint16", "UInt16", "UInt16[pyarrow]", "uint16[pyarrow]"}:
        return dtypes.UInt16()
    if dtype in {"uint8", "UInt8", "UInt8[pyarrow]", "uint8[pyarrow]"}:
        return dtypes.UInt8()
    if dtype in {
        "float64",
        "Float64",
        "Float64[pyarrow]",
        "float64[pyarrow]",
        "double[pyarrow]",
    }:
        return dtypes.Float64()
    if dtype in {
        "float32",
        "Float32",
        "Float32[pyarrow]",
        "float32[pyarrow]",
        "float[pyarrow]",
    }:
        return dtypes.Float32()
    if dtype in {
        # "there is no problem which can't be solved by adding an extra string type" pandas
        "string",
        "string[python]",
        "string[pyarrow]",
        "string[pyarrow_numpy]",
        "large_string[pyarrow]",
        "str",
    }:
        return dtypes.String()
    if dtype in {"bool", "boolean", "boolean[pyarrow]", "bool[pyarrow]"}:
        return dtypes.Boolean()
    if dtype.startswith("dictionary<"):
        return dtypes.Categorical()
    if dtype == "category":
        return native_categorical_to_narwhals_dtype(native_dtype, version)
    if (match_ := PATTERN_PD_DATETIME.match(dtype)) or (
        match_ := PATTERN_PA_DATETIME.match(dtype)
    ):
        dt_time_unit: TimeUnit = match_.group("time_unit")  # type: ignore[assignment]
        dt_time_zone: str | None = match_.group("time_zone")
        return dtypes.Datetime(dt_time_unit, dt_time_zone)
    if (match_ := PATTERN_PD_DURATION.match(dtype)) or (
        match_ := PATTERN_PA_DURATION.match(dtype)
    ):
        du_time_unit: TimeUnit = match_.group("time_unit")  # type: ignore[assignment]
        return dtypes.Duration(du_time_unit)
    if dtype == "date32[day][pyarrow]":
        return dtypes.Date()
    if dtype.startswith("decimal") and dtype.endswith("[pyarrow]"):
        return dtypes.Decimal()
    if dtype.startswith("time") and dtype.endswith("[pyarrow]"):
        return dtypes.Time()
    if dtype.startswith("binary") and dtype.endswith("[pyarrow]"):
        return dtypes.Binary()
    return dtypes.Unknown()  # pragma: no cover


def object_native_to_narwhals_dtype(
    series: PandasLikeSeries | None, version: Version, implementation: Implementation
) -> DType:
    dtypes = version.dtypes
    if implementation is Implementation.CUDF:
        # Per conversations with their maintainers, they don't support arbitrary
        # objects, so we can just return String.
        return dtypes.String()

    infer = pd.api.types.infer_dtype
    # Arbitrary limit of 100 elements to use to sniff dtype.
    inferred_dtype = "empty" if series is None else infer(series.head(100), skipna=True)
    if inferred_dtype == "string":
        return dtypes.String()
    if inferred_dtype == "empty" and version is not Version.V1:
        # Default to String for empty Series.
        return dtypes.String()
    if inferred_dtype == "empty":
        # But preserve returning Object in V1.
        return dtypes.Object()
    return dtypes.Object()


def native_categorical_to_narwhals_dtype(
    native_dtype: pd.CategoricalDtype,
    version: Version,
    implementation: Literal[Implementation.CUDF] | None = None,
) -> DType:
    dtypes = version.dtypes
    if version is Version.V1:
        return dtypes.Categorical()
    if native_dtype.ordered:
        into_iter = (
            _cudf_categorical_to_list(native_dtype)
            if implementation is Implementation.CUDF
            else native_dtype.categories.to_list
        )
        return dtypes.Enum(_DeferredIterable(into_iter))
    return dtypes.Categorical()


def _cudf_categorical_to_list(
    native_dtype: Any,
) -> Callable[[], list[Any]]:  # pragma: no cover
    # NOTE: https://docs.rapids.ai/api/cudf/stable/user_guide/api_docs/api/cudf.core.dtypes.categoricaldtype/#cudf.core.dtypes.CategoricalDtype
    def fn() -> list[Any]:
        return native_dtype.categories.to_arrow().to_pylist()

    return fn


def native_to_narwhals_dtype(
    native_dtype: Any,
    version: Version,
    implementation: Implementation,
    *,
    allow_object: bool = False,
) -> DType:
    str_dtype = str(native_dtype)

    if str_dtype.startswith(("large_list", "list", "struct", "fixed_size_list")):
        from narwhals._arrow.utils import (
            native_to_narwhals_dtype as arrow_native_to_narwhals_dtype,
        )

        if hasattr(native_dtype, "to_arrow"):  # pragma: no cover
            # cudf, cudf.pandas
            return arrow_native_to_narwhals_dtype(native_dtype.to_arrow(), version)
        return arrow_native_to_narwhals_dtype(native_dtype.pyarrow_dtype, version)
    if str_dtype == "category" and implementation.is_cudf():
        # https://github.com/rapidsai/cudf/issues/18536
        # https://github.com/rapidsai/cudf/issues/14027
        return native_categorical_to_narwhals_dtype(
            native_dtype, version, Implementation.CUDF
        )
    if str_dtype != "object":
        return non_object_native_to_narwhals_dtype(native_dtype, version)
    if implementation is Implementation.DASK:
        # Per conversations with their maintainers, they don't support arbitrary
        # objects, so we can just return String.
        return version.dtypes.String()
    if allow_object:
        return object_native_to_narwhals_dtype(None, version, implementation)
    msg = (
        "Unreachable code, object dtype should be handled separately"  # pragma: no cover
    )
    raise AssertionError(msg)


if Implementation.PANDAS._backend_version() >= (1, 2):

    def is_dtype_numpy_nullable(dtype: Any) -> TypeIs[BaseMaskedDtype]:
        """Return `True` if `dtype` is `"numpy_nullable"`."""
        # NOTE: We need a sentinel as the positive case is `BaseMaskedDtype.base = None`
        # See https://github.com/narwhals-dev/narwhals/pull/2740#discussion_r2171667055
        sentinel = object()
        return (
            isinstance(dtype, pd.api.extensions.ExtensionDtype)
            and getattr(dtype, "base", sentinel) is None
        )
else:  # pragma: no cover

    def is_dtype_numpy_nullable(dtype: Any) -> TypeIs[BaseMaskedDtype]:
        # NOTE: `base` attribute was added between 1.1-1.2
        # Checking by isinstance requires using an import path that is no longer valid
        # `1.1`: https://github.com/pandas-dev/pandas/blob/b5958ee1999e9aead1938c0bba2b674378807b3d/pandas/core/arrays/masked.py#L37
        # `1.2`: https://github.com/pandas-dev/pandas/blob/7c48ff4409c622c582c56a5702373f726de08e96/pandas/core/arrays/masked.py#L41
        # `1.5`: https://github.com/pandas-dev/pandas/blob/35b0d1dcadf9d60722c055ee37442dc76a29e64c/pandas/core/dtypes/dtypes.py#L1609
        if isinstance(dtype, pd.api.extensions.ExtensionDtype):
            from pandas.core.arrays.masked import (  # type: ignore[attr-defined]
                BaseMaskedDtype as OldBaseMaskedDtype,  # pyright: ignore[reportAttributeAccessIssue]
            )

            return isinstance(dtype, OldBaseMaskedDtype)
        return False


def get_dtype_backend(dtype: Any, implementation: Implementation) -> DTypeBackend:
    """Get dtype backend for pandas type.

    Matches pandas' `dtype_backend` argument in `convert_dtypes`.
    """
    if implementation is Implementation.CUDF:
        return None
    if is_dtype_pyarrow(dtype):
        return "pyarrow"
    return "numpy_nullable" if is_dtype_numpy_nullable(dtype) else None


# NOTE: Use this to avoid annotating inline
def iter_dtype_backends(
    dtypes: Iterable[Any], implementation: Implementation
) -> Iterator[DTypeBackend]:
    """Yield a `DTypeBackend` per-dtype.

    Matches pandas' `dtype_backend` argument in `convert_dtypes`.
    """
    return (get_dtype_backend(dtype, implementation) for dtype in dtypes)


@functools.lru_cache(maxsize=16)
def is_dtype_pyarrow(dtype: Any) -> TypeIs[pd.ArrowDtype]:
    return hasattr(pd, "ArrowDtype") and isinstance(dtype, pd.ArrowDtype)


dtypes = Version.MAIN.dtypes
NW_TO_PD_DTYPES_INVARIANT: Mapping[type[DType], str] = {
    # TODO(Unassigned): is there no pyarrow-backed categorical?
    # or at least, convert_dtypes(dtype_backend='pyarrow') doesn't
    # convert to it?
    dtypes.Categorical: "category",
    dtypes.Object: "object",
}
NW_TO_PD_DTYPES_BACKEND: Mapping[type[DType], Mapping[DTypeBackend, str | type[Any]]] = {
    dtypes.Float64: {
        "pyarrow": "Float64[pyarrow]",
        "numpy_nullable": "Float64",
        None: "float64",
    },
    dtypes.Float32: {
        "pyarrow": "Float32[pyarrow]",
        "numpy_nullable": "Float32",
        None: "float32",
    },
    dtypes.Int64: {"pyarrow": "Int64[pyarrow]", "numpy_nullable": "Int64", None: "int64"},
    dtypes.Int32: {"pyarrow": "Int32[pyarrow]", "numpy_nullable": "Int32", None: "int32"},
    dtypes.Int16: {"pyarrow": "Int16[pyarrow]", "numpy_nullable": "Int16", None: "int16"},
    dtypes.Int8: {"pyarrow": "Int8[pyarrow]", "numpy_nullable": "Int8", None: "int8"},
    dtypes.UInt64: {
        "pyarrow": "UInt64[pyarrow]",
        "numpy_nullable": "UInt64",
        None: "uint64",
    },
    dtypes.UInt32: {
        "pyarrow": "UInt32[pyarrow]",
        "numpy_nullable": "UInt32",
        None: "uint32",
    },
    dtypes.UInt16: {
        "pyarrow": "UInt16[pyarrow]",
        "numpy_nullable": "UInt16",
        None: "uint16",
    },
    dtypes.UInt8: {"pyarrow": "UInt8[pyarrow]", "numpy_nullable": "UInt8", None: "uint8"},
    dtypes.Boolean: {
        "pyarrow": "boolean[pyarrow]",
        "numpy_nullable": "boolean",
        None: "bool",
    },
}
UNSUPPORTED_DTYPES = (dtypes.Decimal,)


def narwhals_to_native_dtype(  # noqa: C901, PLR0912
    dtype: IntoDType,
    dtype_backend: DTypeBackend,
    implementation: Implementation,
    version: Version,
) -> str | PandasDtype:
    if dtype_backend not in {None, "pyarrow", "numpy_nullable"}:
        msg = f"Expected one of {{None, 'pyarrow', 'numpy_nullable'}}, got: '{dtype_backend}'"
        raise ValueError(msg)
    dtypes = version.dtypes
    base_type = dtype.base_type()
    if pd_type := NW_TO_PD_DTYPES_INVARIANT.get(base_type):
        return pd_type
    if into_pd_type := NW_TO_PD_DTYPES_BACKEND.get(base_type):
        return into_pd_type[dtype_backend]
    if issubclass(base_type, dtypes.String):
        if dtype_backend == "pyarrow":
            import pyarrow as pa  # ignore-banned-import

            # Note: this is different from `string[pyarrow]`, even though the repr
            # looks the same.
            # >>> pd.DataFrame({'a':['foo']}, dtype='string[pyarrow]')['a'].str.len()
            # 0    3
            # Name: a, dtype: Int64
            # >>> pd.DataFrame({'a':['foo']}, dtype=pd.ArrowDtype(pa.string()))['a'].str.len()
            # 0    3
            # Name: a, dtype: int32[pyarrow]
            #
            # `ArrowDType(pa.string())` is what `.convert_dtypes(dtype_backend='pyarrow')` converts to,
            # so we use that here.
            return pd.ArrowDtype(pa.string())
        if dtype_backend == "numpy_nullable":
            return "string"
        return str
    if isinstance_or_issubclass(dtype, dtypes.Datetime):
        if is_pandas_or_modin(implementation) and PANDAS_VERSION < (
            2,
        ):  # pragma: no cover
            if isinstance(dtype, dtypes.Datetime) and dtype.time_unit != "ns":
                found = requires._unparse_version(PANDAS_VERSION)
                available = f"available in 'pandas>=2.0', found version {found!r}."
                changelog_url = "https://pandas.pydata.org/docs/dev/whatsnew/v2.0.0.html#construction-with-datetime64-or-timedelta64-dtype-with-unsupported-resolution"
                msg = (
                    f"`nw.Datetime(time_unit={dtype.time_unit!r})` is only {available}\n"
                    "Narwhals has fallen back to using `time_unit='ns'` to avoid an error.\n\n"
                    "Hint: to avoid this warning, consider either:\n"
                    f"- Upgrading pandas: {changelog_url}\n"
                    f"- Using a bare `nw.Datetime`, if this precision is not important"
                )
                issue_warning(msg, UserWarning)
            dt_time_unit = "ns"
        else:
            dt_time_unit = dtype.time_unit

        if dtype_backend == "pyarrow":
            tz_part = f", tz={tz}" if (tz := dtype.time_zone) else ""
            return f"timestamp[{dt_time_unit}{tz_part}][pyarrow]"
        tz_part = f", {tz}" if (tz := dtype.time_zone) else ""
        return f"datetime64[{dt_time_unit}{tz_part}]"
    if isinstance_or_issubclass(dtype, dtypes.Duration):
        if is_pandas_or_modin(implementation) and PANDAS_VERSION < (
            2,
        ):  # pragma: no cover
            du_time_unit = "ns"
        else:
            du_time_unit = dtype.time_unit
        return (
            f"duration[{du_time_unit}][pyarrow]"
            if dtype_backend == "pyarrow"
            else f"timedelta64[{du_time_unit}]"
        )
    if isinstance_or_issubclass(dtype, dtypes.Date):
        try:
            import pyarrow as pa  # ignore-banned-import
        except ModuleNotFoundError as exc:
            # BUG: Never re-raised?
            msg = "'pyarrow>=13.0.0' is required for `Date` dtype."
            raise ModuleNotFoundError(msg) from exc
        return "date32[pyarrow]"
    if isinstance_or_issubclass(dtype, dtypes.Enum):
        if version is Version.V1:
            msg = "Converting to Enum is not supported in narwhals.stable.v1"
            raise NotImplementedError(msg)
        if isinstance(dtype, dtypes.Enum):
            ns = implementation.to_native_namespace()
            return ns.CategoricalDtype(dtype.categories, ordered=True)
        msg = "Can not cast / initialize Enum without categories present"
        raise ValueError(msg)
    if issubclass(
        base_type, (dtypes.Struct, dtypes.Array, dtypes.List, dtypes.Time, dtypes.Binary)
    ):
        return narwhals_to_native_arrow_dtype(dtype, implementation, version)
    if issubclass(base_type, UNSUPPORTED_DTYPES):
        msg = f"Converting to {base_type.__name__} dtype is not supported for {implementation}."
        raise NotImplementedError(msg)
    msg = f"Unknown dtype: {dtype}"  # pragma: no cover
    raise AssertionError(msg)


def narwhals_to_native_arrow_dtype(
    dtype: IntoDType, implementation: Implementation, version: Version
) -> pd.ArrowDtype:
    if is_pandas_or_modin(implementation) and PANDAS_VERSION >= (2, 2):
        try:
            import pyarrow as pa  # ignore-banned-import  # noqa: F401
        except ImportError as exc:  # pragma: no cover
            msg = f"Unable to convert to {dtype} to to the following exception: {exc.msg}"
            raise ImportError(msg) from exc
        from narwhals._arrow.utils import narwhals_to_native_dtype as _to_arrow_dtype

        return pd.ArrowDtype(_to_arrow_dtype(dtype, version))
    msg = (  # pragma: no cover
        f"Converting to {dtype} dtype is not supported for implementation "
        f"{implementation} and version {version}."
    )
    raise NotImplementedError(msg)


def int_dtype_mapper(dtype: Any) -> str:
    if "pyarrow" in str(dtype):
        return "Int64[pyarrow]"
    if str(dtype).lower() != str(dtype):  # pragma: no cover
        return "Int64"
    return "int64"


_TIMESTAMP_DATETIME_OP_FACTOR: Mapping[
    tuple[UnitCurrent, UnitTarget], tuple[BinOpBroadcast, IntoRhs]
] = {
    ("ns", "us"): (operator.floordiv, 1_000),
    ("ns", "ms"): (operator.floordiv, 1_000_000),
    ("us", "ns"): (operator.mul, NS_PER_MICROSECOND),
    ("us", "ms"): (operator.floordiv, 1_000),
    ("ms", "ns"): (operator.mul, NS_PER_MILLISECOND),
    ("ms", "us"): (operator.mul, 1_000),
    ("s", "ns"): (operator.mul, NS_PER_SECOND),
    ("s", "us"): (operator.mul, US_PER_SECOND),
    ("s", "ms"): (operator.mul, MS_PER_SECOND),
}


def calculate_timestamp_datetime(
    s: NativeSeriesT, current: TimeUnit, time_unit: TimeUnit
) -> NativeSeriesT:
    if current == time_unit:
        return s
    if item := _TIMESTAMP_DATETIME_OP_FACTOR.get((current, time_unit)):
        fn, factor = item
        return fn(s, factor)
    msg = (  # pragma: no cover
        f"unexpected time unit {current}, please report an issue at "
        "https://github.com/narwhals-dev/narwhals"
    )
    raise AssertionError(msg)


_TIMESTAMP_DATE_FACTOR: Mapping[TimeUnit, int] = {
    "ns": NS_PER_SECOND,
    "us": US_PER_SECOND,
    "ms": MS_PER_SECOND,
    "s": 1,
}


def calculate_timestamp_date(s: NativeSeriesT, time_unit: TimeUnit) -> NativeSeriesT:
    return s * SECONDS_PER_DAY * _TIMESTAMP_DATE_FACTOR[time_unit]


def select_columns_by_name(
    df: NativeDataFrameT,
    column_names: list[str] | _1DArray,  # NOTE: Cannot be a tuple!
    implementation: Implementation,
) -> NativeDataFrameT | Any:
    """Select columns by name.

    Prefer this over `df.loc[:, column_names]` as it's
    generally more performant.
    """
    if len(column_names) == df.shape[1] and (df.columns == column_names).all():
        return df
    if (df.columns.dtype.kind == "b") or (
        implementation is Implementation.PANDAS
        and implementation._backend_version() < (1, 5)
    ):
        # See https://github.com/narwhals-dev/narwhals/issues/1349#issuecomment-2470118122
        # for why we need this
        if error := check_columns_exist(column_names, available=df.columns.tolist()):
            raise error
        return df.loc[:, column_names]
    try:
        return df[column_names]
    except KeyError as e:
        if error := check_columns_exist(column_names, available=df.columns.tolist()):
            raise error from e
        raise


def is_non_nullable_boolean(s: PandasLikeSeries) -> bool:
    # cuDF booleans are nullable but the native dtype is still 'bool'.
    return (
        s._implementation
        in {Implementation.PANDAS, Implementation.MODIN, Implementation.DASK}
        and s.native.dtype == "bool"
    )


def import_array_module(implementation: Implementation, /) -> ModuleType:
    """Returns numpy or cupy module depending on the given implementation."""
    if implementation in {Implementation.PANDAS, Implementation.MODIN}:
        import numpy as np

        return np
    if implementation is Implementation.CUDF:
        import cupy as cp  # ignore-banned-import  # cuDF dependency.

        return cp
    msg = f"Expected pandas/modin/cudf, got: {implementation}"  # pragma: no cover
    raise AssertionError(msg)


class PandasLikeSeriesNamespace(EagerSeriesNamespace["PandasLikeSeries", Any]): ...
