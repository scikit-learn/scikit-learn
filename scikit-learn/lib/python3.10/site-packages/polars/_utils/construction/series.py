from __future__ import annotations

import contextlib
from collections.abc import Generator, Iterator, Mapping
from datetime import date, datetime, time, timedelta
from enum import Enum as PyEnum
from itertools import islice
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
)

import polars._reexport as pl
import polars._utils.construction as plc
from polars._dependencies import (
    _PYARROW_AVAILABLE,
    _check_for_numpy,
    dataclasses,
)
from polars._dependencies import numpy as np
from polars._dependencies import pandas as pd
from polars._dependencies import pyarrow as pa
from polars._utils.construction.dataframe import _sequence_of_dict_to_pydf
from polars._utils.construction.utils import (
    get_first_non_none,
    is_namedtuple,
    is_pydantic_model,
    is_simple_numpy_backed_pandas_series,
    is_sqlalchemy_row,
)
from polars._utils.various import (
    range_to_series,
)
from polars._utils.wrap import wrap_s
from polars.datatypes import (
    Array,
    Boolean,
    Categorical,
    Date,
    Datetime,
    Decimal,
    Duration,
    Enum,
    List,
    Null,
    Object,
    String,
    Struct,
    Time,
    Unknown,
    dtype_to_py_type,
    is_polars_dtype,
    numpy_char_code_to_dtype,
    parse_into_dtype,
    try_parse_into_dtype,
)
from polars.datatypes.constructor import (
    numpy_type_to_constructor,
    numpy_values_and_dtype,
    polars_type_to_constructor,
    py_type_to_constructor,
)

with contextlib.suppress(ImportError):  # Module not available when building docs
    from polars._plr import PySeries

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

    from polars import DataFrame, Series
    from polars._dependencies import pandas as pd
    from polars._typing import PolarsDataType


def sequence_to_pyseries(
    name: str,
    values: Sequence[Any],
    dtype: PolarsDataType | None = None,
    *,
    strict: bool = True,
    nan_to_null: bool = False,
) -> PySeries:
    """Construct a PySeries from a sequence."""
    python_dtype: type | None = None

    if isinstance(values, range):
        return range_to_series(name, values, dtype=dtype)._s

    # empty sequence
    if len(values) == 0 and dtype is None:
        # if dtype for empty sequence could be guessed
        # (e.g comparisons between self and other), default to Null
        dtype = Null

    # lists defer to subsequent handling; identify nested type
    elif dtype in (List, Array):
        python_dtype = list

    # infer temporal type handling
    py_temporal_types = {date, datetime, timedelta, time}
    pl_temporal_types = {Date, Datetime, Duration, Time}

    value = get_first_non_none(values)
    if value is not None:
        if (
            dataclasses.is_dataclass(value)
            or is_pydantic_model(value)
            or is_namedtuple(value.__class__)
            or is_sqlalchemy_row(value)
        ) and dtype != Object:
            return pl.DataFrame(values).to_struct(name)._s
        elif (
            not isinstance(value, dict) and isinstance(value, Mapping)
        ) and dtype != Object:
            return _sequence_of_dict_to_pydf(
                value,
                data=values,
                strict=strict,
                schema_overrides=None,
                infer_schema_length=None,
                schema=None,
            ).to_struct(name, [])
        elif isinstance(value, range) and dtype is None:
            values = [range_to_series("", v) for v in values]
        else:
            # for temporal dtypes:
            # * if the values are integer, we take the physical branch.
            # * if the values are python types, take the temporal branch.
            # * if the values are ISO-8601 strings, init then convert via strptime.
            # * if the values are floats/other dtypes, this is an error.
            if dtype in py_temporal_types and isinstance(value, int):
                dtype = parse_into_dtype(dtype)  # construct from integer
            elif (
                dtype in pl_temporal_types or type(dtype) in pl_temporal_types
            ) and not isinstance(value, int):
                python_dtype = dtype_to_py_type(dtype)  # type: ignore[arg-type]

    # if values are enums, infer and load the appropriate dtype/values
    if issubclass(type(value), PyEnum):
        if dtype is None and python_dtype is None:
            with contextlib.suppress(TypeError):
                dtype = Enum(type(value))
        if not isinstance(value, (str, int)):
            values = [v.value for v in values]

    # physical branch
    # flat data
    if (
        dtype is not None
        and is_polars_dtype(dtype)
        and not dtype.is_nested()
        and dtype != Unknown
        and (python_dtype is None)
    ):
        constructor = polars_type_to_constructor(dtype)
        pyseries = _construct_series_with_fallbacks(
            constructor, name, values, dtype, strict=strict
        )
        if dtype in (
            Date,
            Datetime,
            Duration,
            Time,
            Boolean,
            Categorical,
            Enum,
        ) or isinstance(dtype, (Categorical, Decimal)):
            if pyseries.dtype() != dtype:
                pyseries = pyseries.cast(dtype, strict=strict, wrap_numerical=False)

        # Uninstanced Decimal is a bit special and has various inference paths
        if dtype == Decimal:
            if pyseries.dtype() == String:
                pyseries = pyseries.str_to_decimal_infer(inference_length=0)
            elif pyseries.dtype().is_float():
                # Go through string so we infer an appropriate scale.
                pyseries = pyseries.cast(
                    String, strict=strict, wrap_numerical=False
                ).str_to_decimal_infer(inference_length=0)
            elif pyseries.dtype().is_integer() or pyseries.dtype() == Null:
                pyseries = pyseries.cast(
                    Decimal(scale=0), strict=strict, wrap_numerical=False
                )
            elif not isinstance(pyseries.dtype(), Decimal):
                msg = f"can't convert {pyseries.dtype()} to Decimal"
                raise TypeError(msg)

        return pyseries

    elif dtype == Struct:
        # This is very bad. Goes via rows? And needs to do outer nullability separate.
        # It also has two data passes.
        # TODO: eventually go into struct builder
        struct_schema = dtype.to_schema() if isinstance(dtype, Struct) else None
        empty = {}  # type: ignore[var-annotated]

        data = []
        invalid = []
        for i, v in enumerate(values):
            if v is None:
                invalid.append(i)
                data.append(empty)
            else:
                data.append(v)

        return plc.sequence_to_pydf(
            data=data,
            schema=struct_schema,
            orient="row",
        ).to_struct(name, invalid)

    if python_dtype is None:
        if value is None:
            constructor = polars_type_to_constructor(Null)
            return constructor(name, values, strict)

        # generic default dtype
        python_dtype = type(value)

    # temporal branch
    if issubclass(python_dtype, tuple(py_temporal_types)):
        if dtype is None:
            dtype = parse_into_dtype(python_dtype)  # construct from integer
        elif dtype in py_temporal_types:
            dtype = parse_into_dtype(dtype)

        values_dtype = None if value is None else try_parse_into_dtype(type(value))
        if values_dtype is not None and values_dtype.is_float():
            msg = f"'float' object cannot be interpreted as a {python_dtype.__name__!r}"
            raise TypeError(
                # we do not accept float values as temporal; if this is
                # required, the caller should explicitly cast to int first.
                msg
            )

        # We use the AnyValue builder to create the datetime array
        # We store the values internally as UTC and set the timezone
        py_series = PySeries.new_from_any_values(name, values, strict)

        time_unit = getattr(dtype, "time_unit", None)
        time_zone = getattr(dtype, "time_zone", None)

        if dtype.is_temporal() and values_dtype == String and dtype != Duration:
            s = wrap_s(py_series).str.strptime(dtype, strict=strict)  # type: ignore[arg-type]
        elif time_unit is not None and values_dtype != Date:
            s = wrap_s(py_series).dt.cast_time_unit(time_unit)
        else:
            s = wrap_s(py_series)

        if (values_dtype == Date) & (dtype == Datetime):
            s = s.cast(Datetime(time_unit or "us"))

        if dtype == Datetime and time_zone is not None:
            return s.dt.convert_time_zone(time_zone)._s
        return s._s

    elif (
        _check_for_numpy(value)
        and isinstance(value, np.ndarray)
        and len(value.shape) == 1
    ):
        n_elems = len(value)
        if all(len(v) == n_elems for v in values):
            # can take (much) faster path if all lists are the same length
            return numpy_to_pyseries(
                name,
                np.vstack(values),
                strict=strict,
                nan_to_null=nan_to_null,
            )
        else:
            return PySeries.new_series_list(
                name,
                [
                    numpy_to_pyseries("", v, strict=strict, nan_to_null=nan_to_null)
                    for v in values
                ],
                strict,
            )

    elif python_dtype in (list, tuple):
        if dtype is None:
            return PySeries.new_from_any_values(name, values, strict=strict)
        elif dtype == Object:
            return PySeries.new_object(name, values, strict)
        else:
            if (inner_dtype := getattr(dtype, "inner", None)) is not None:
                pyseries_list = [
                    None
                    if value is None
                    else sequence_to_pyseries(
                        "",
                        value,
                        inner_dtype,
                        strict=strict,
                        nan_to_null=nan_to_null,
                    )
                    for value in values
                ]
                pyseries = PySeries.new_series_list(name, pyseries_list, strict)
            else:
                pyseries = PySeries.new_from_any_values_and_dtype(
                    name, values, dtype, strict=strict
                )
            if dtype != pyseries.dtype():
                pyseries = pyseries.cast(dtype, strict=False, wrap_numerical=False)
            return pyseries

    elif python_dtype == pl.Series:
        return PySeries.new_series_list(
            name, [v._s if v is not None else None for v in values], strict
        )

    elif python_dtype == PySeries:
        return PySeries.new_series_list(name, values, strict)
    else:
        constructor = py_type_to_constructor(python_dtype)
        if constructor == PySeries.new_object:
            try:
                srs = PySeries.new_from_any_values(name, values, strict)
                if _check_for_numpy(python_dtype, check_type=False) and isinstance(
                    np.bool_(True), np.generic
                ):
                    dtype = numpy_char_code_to_dtype(np.dtype(python_dtype).char)
                    return srs.cast(dtype, strict=strict, wrap_numerical=False)
                else:
                    return srs

            except RuntimeError:
                return PySeries.new_from_any_values(name, values, strict=strict)

        return _construct_series_with_fallbacks(
            constructor, name, values, dtype, strict=strict
        )


def _construct_series_with_fallbacks(
    constructor: Callable[[str, Sequence[Any], bool], PySeries],
    name: str,
    values: Sequence[Any],
    dtype: PolarsDataType | None,
    *,
    strict: bool,
) -> PySeries:
    """Construct Series, with fallbacks for basic type mismatch (eg: bool/int)."""
    try:
        return constructor(name, values, strict)
    except (TypeError, OverflowError) as e:
        # # This retry with i64 is related to https://github.com/pola-rs/polars/issues/17231
        # # Essentially, when given a [0, u64::MAX] then it would Overflow.
        if (
            isinstance(e, OverflowError)
            and dtype is None
            and constructor == PySeries.new_opt_i64
        ):
            return _construct_series_with_fallbacks(
                PySeries.new_opt_u64, name, values, dtype, strict=strict
            )
        elif dtype is None:
            return PySeries.new_from_any_values(name, values, strict=strict)
        else:
            return PySeries.new_from_any_values_and_dtype(
                name, values, dtype, strict=strict
            )


def iterable_to_pyseries(
    name: str,
    values: Iterable[Any],
    dtype: PolarsDataType | None = None,
    *,
    chunk_size: int = 1_000_000,
    strict: bool = True,
) -> PySeries:
    """Construct a PySeries from an iterable/generator."""
    if not isinstance(values, (Generator, Iterator)):
        values = iter(values)

    def to_series_chunk(values: list[Any], dtype: PolarsDataType | None) -> Series:
        return pl.Series(
            name=name,
            values=values,
            dtype=dtype,
            strict=strict,
        )

    n_chunks = 0
    series: Series = None  # type: ignore[assignment]
    while True:
        slice_values = list(islice(values, chunk_size))
        if not slice_values:
            break
        schunk = to_series_chunk(slice_values, dtype)
        if series is None:
            series = schunk
            dtype = series.dtype
        else:
            series.append(schunk)
            n_chunks += 1

    if series is None:
        series = to_series_chunk([], dtype)
    if n_chunks > 0:
        series.rechunk(in_place=True)

    return series._s


def pandas_to_pyseries(
    name: str,
    values: pd.Series[Any] | pd.Index[Any] | pd.DatetimeIndex,
    dtype: PolarsDataType | None = None,
    *,
    strict: bool = True,
    nan_to_null: bool = True,
) -> PySeries:
    """Construct a PySeries from a pandas Series or DatetimeIndex."""
    if not name and values.name is not None:
        name = str(values.name)
    if is_simple_numpy_backed_pandas_series(values):
        return pl.Series(
            name, values.to_numpy(), dtype=dtype, nan_to_null=nan_to_null, strict=strict
        )._s
    if not _PYARROW_AVAILABLE:
        msg = (
            "pyarrow is required for converting a pandas series to Polars, "
            "unless it is a simple numpy-backed one "
            "(e.g. 'int64', 'bool', 'float32' - not 'Int64')"
        )
        raise ImportError(msg)
    return arrow_to_pyseries(
        name,
        plc.pandas_series_to_arrow(values, nan_to_null=nan_to_null),
        dtype=dtype,
        strict=strict,
    )


def arrow_to_pyseries(
    name: str,
    values: pa.Array,
    dtype: PolarsDataType | None = None,
    *,
    strict: bool = True,
    rechunk: bool = True,
) -> PySeries:
    """Construct a PySeries from an Arrow array."""
    array = plc.coerce_arrow(values)

    # special handling of empty categorical arrays
    if (
        len(array) == 0
        and isinstance(array.type, pa.DictionaryType)
        and array.type.value_type
        in (
            pa.utf8(),
            pa.large_utf8(),
        )
    ):
        pys = pl.Series(name, [], dtype=Categorical)._s

    elif not hasattr(array, "num_chunks"):
        pys = PySeries.from_arrow(name, array)
    else:
        if array.num_chunks > 1:
            # somehow going through ffi with a structarray
            # returns the first chunk every time
            if isinstance(array.type, pa.StructType):
                pys = PySeries.from_arrow(name, array.combine_chunks())
            else:
                it = array.iterchunks()
                pys = PySeries.from_arrow(name, next(it))
                for a in it:
                    pys.append(PySeries.from_arrow(name, a))
        elif array.num_chunks == 0:
            pys = PySeries.from_arrow(name, pa.nulls(0, type=array.type))
        else:
            pys = PySeries.from_arrow(name, array.chunks[0])

        if rechunk:
            pys.rechunk(in_place=True)

    return (
        pys.cast(dtype, strict=strict, wrap_numerical=False)
        if dtype is not None
        else pys
    )


def numpy_to_pyseries(
    name: str,
    values: np.ndarray[Any, Any],
    *,
    strict: bool = True,
    nan_to_null: bool = False,
) -> PySeries:
    """Construct a PySeries from a numpy array."""
    values = np.ascontiguousarray(values)

    if values.ndim == 1:
        values, dtype = numpy_values_and_dtype(values)
        constructor = numpy_type_to_constructor(values, dtype)
        return constructor(
            name, values, nan_to_null if dtype in (np.float32, np.float64) else strict
        )
    else:
        original_shape = values.shape
        values_1d = values.reshape(-1)

        from polars.series.utils import _with_no_check_length

        py_s = _with_no_check_length(
            lambda: numpy_to_pyseries(
                name,
                values_1d,
                strict=strict,
                nan_to_null=nan_to_null,
            )
        )
        return wrap_s(py_s).reshape(original_shape)._s


def series_to_pyseries(
    name: str | None,
    values: Series,
    *,
    dtype: PolarsDataType | None = None,
    strict: bool = True,
) -> PySeries:
    """Construct a new PySeries from a Polars Series."""
    s = values.clone()
    if dtype is not None and dtype != s.dtype:
        s = s.cast(dtype, strict=strict)
    if name is not None:
        s = s.alias(name)
    return s._s


def dataframe_to_pyseries(
    name: str | None,
    values: DataFrame,
    *,
    dtype: PolarsDataType | None = None,
    strict: bool = True,
) -> PySeries:
    """Construct a new PySeries from a Polars DataFrame."""
    if values.width > 1:
        name = name or ""
        s = values.to_struct(name)
    elif values.width == 1:
        s = values.to_series()
        if name is not None:
            s = s.alias(name)
    else:
        msg = "cannot initialize Series from DataFrame without any columns"
        raise TypeError(msg)

    if dtype is not None and dtype != s.dtype:
        s = s.cast(dtype, strict=strict)

    return s._s
