from __future__ import annotations

import contextlib
from collections.abc import Generator, Mapping, Sequence
from datetime import date, datetime, time, timedelta
from functools import singledispatch
from itertools import islice, zip_longest
from operator import itemgetter
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
)

import polars._reexport as pl
import polars._utils.construction as plc
from polars import functions as F
from polars._dependencies import (
    _NUMPY_AVAILABLE,
    _PYARROW_AVAILABLE,
    _check_for_numpy,
    _check_for_pandas,
    dataclasses,
)
from polars._dependencies import numpy as np
from polars._dependencies import pandas as pd
from polars._dependencies import pyarrow as pa
from polars._utils.construction.utils import (
    contains_nested,
    get_first_non_none,
    is_namedtuple,
    is_pydantic_model,
    is_simple_numpy_backed_pandas_series,
    is_sqlalchemy_row,
    nt_unpack,
    try_get_type_hints,
)
from polars._utils.various import (
    _is_generator,
    arrlen,
    issue_warning,
    parse_version,
)
from polars.datatypes import (
    N_INFER_DEFAULT,
    Categorical,
    Duration,
    Enum,
    String,
    Struct,
    Unknown,
    is_polars_dtype,
    parse_into_dtype,
    try_parse_into_dtype,
)
from polars.exceptions import DataOrientationWarning, ShapeError
from polars.meta import thread_pool_size

with contextlib.suppress(ImportError):  # Module not available when building docs
    from polars._plr import PyDataFrame

if TYPE_CHECKING:
    from collections.abc import Iterable, MutableMapping

    from polars import DataFrame, Series
    from polars._plr import PySeries
    from polars._typing import (
        Orientation,
        PolarsDataType,
        SchemaDefinition,
        SchemaDict,
    )

_MIN_NUMPY_SIZE_FOR_MULTITHREADING = 1000


def dict_to_pydf(
    data: Mapping[str, Sequence[object] | Mapping[str, Sequence[object]] | Series],
    schema: SchemaDefinition | None = None,
    *,
    schema_overrides: SchemaDict | None = None,
    strict: bool = True,
    nan_to_null: bool = False,
    allow_multithreaded: bool = True,
) -> PyDataFrame:
    """Construct a PyDataFrame from a dictionary of sequences."""
    if isinstance(schema, Mapping) and data:
        if not all((col in schema) for col in data):
            msg = "the given column-schema names do not match the data dictionary"
            raise ValueError(msg)
        data = {col: data[col] for col in schema}

    column_names, schema_overrides = _unpack_schema(
        schema, lookup_names=data.keys(), schema_overrides=schema_overrides
    )
    if not column_names:
        column_names = list(data)

    if data and _NUMPY_AVAILABLE:
        # if there are 3 or more numpy arrays of sufficient size, we multi-thread:
        count_numpy = sum(
            int(
                allow_multithreaded
                and _check_for_numpy(val)
                and isinstance(val, np.ndarray)
                and len(val) > _MIN_NUMPY_SIZE_FOR_MULTITHREADING
                # integers and non-nan floats are zero-copy
                and nan_to_null
                and val.dtype in (np.float32, np.float64)
            )
            for val in data.values()
        )
        if count_numpy >= 3:
            # yes, multi-threading was easier in python here; we cannot have multiple
            # threads running python and release the gil in pyo3 (it will deadlock).

            # (note: 'dummy' is threaded)
            # We catch FileNotFoundError: see 16675
            try:
                import multiprocessing.dummy

                pool_size = thread_pool_size()
                with multiprocessing.dummy.Pool(pool_size) as pool:
                    data = dict(
                        zip(
                            column_names,
                            pool.map(
                                lambda t: (
                                    pl.Series(t[0], t[1], nan_to_null=nan_to_null)
                                    if isinstance(t[1], np.ndarray)
                                    else t[1]
                                ),
                                list(data.items()),
                            ),
                        )
                    )
            except FileNotFoundError:
                return dict_to_pydf(
                    data=data,
                    schema=schema,
                    schema_overrides=schema_overrides,
                    strict=strict,
                    nan_to_null=nan_to_null,
                    allow_multithreaded=False,
                )

    if not data and schema_overrides:
        data_series = [
            pl.Series(
                name,
                [],
                dtype=schema_overrides.get(name),
                strict=strict,
                nan_to_null=nan_to_null,
            )._s
            for name in column_names
        ]
    else:
        data_series = [
            s._s
            for s in _expand_dict_values(
                data,
                schema_overrides=schema_overrides,
                strict=strict,
                nan_to_null=nan_to_null,
            ).values()
        ]

    data_series = _handle_columns_arg(data_series, columns=column_names, from_dict=True)
    pydf = PyDataFrame(data_series)

    if schema_overrides and pydf.dtypes() != list(schema_overrides.values()):
        pydf = _post_apply_columns(
            pydf, column_names, schema_overrides=schema_overrides, strict=strict
        )
    return pydf


def _unpack_schema(
    schema: SchemaDefinition | None,
    *,
    schema_overrides: SchemaDict | None = None,
    n_expected: int | None = None,
    lookup_names: Iterable[str] | None = None,
) -> tuple[list[str], SchemaDict]:
    """
    Unpack column names and create dtype lookup.

    Works for any (name, dtype) pairs or schema dict input,
    overriding any inferred dtypes with explicit dtypes if supplied.
    """

    def _normalize_dtype(dtype: Any) -> PolarsDataType:
        """Parse non-Polars data types as Polars data types."""
        if is_polars_dtype(dtype, include_unknown=True):
            return dtype
        else:
            return parse_into_dtype(dtype)

    def _parse_schema_overrides(
        schema_overrides: SchemaDict | None = None,
    ) -> dict[str, PolarsDataType]:
        """Parse schema overrides as a dictionary of name to Polars data type."""
        if schema_overrides is None:
            return {}

        return {
            name: _normalize_dtype(dtype) for name, dtype in schema_overrides.items()
        }

    schema_overrides = _parse_schema_overrides(schema_overrides)

    # fast path for empty schema
    if not schema:
        columns = (
            [f"column_{i}" for i in range(n_expected)] if n_expected is not None else []
        )
        return columns, schema_overrides

    # determine column names from schema
    if isinstance(schema, Mapping):
        column_names: list[str] = list(schema)
        schema = list(schema.items())
    else:
        column_names = []
        for i, col in enumerate(schema):
            if isinstance(col, str):
                unnamed = not col and col not in schema_overrides
                col = f"column_{i}" if unnamed else col
            else:
                col = col[0]
            column_names.append(col)

    if n_expected is not None and len(column_names) != n_expected:
        msg = "data does not match the number of columns"
        raise ShapeError(msg)

    # determine column dtypes from schema and lookup_names
    lookup: dict[str, str] | None = (
        {
            col: name
            for col, name in zip_longest(column_names, lookup_names)
            if name is not None
        }
        if lookup_names
        else None
    )

    column_dtypes: dict[str, PolarsDataType] = {}
    for col in schema:
        if isinstance(col, str):
            continue

        name, dtype = col
        if dtype is None:
            continue
        else:
            dtype = _normalize_dtype(dtype)
        name = lookup.get(name, name) if lookup else name
        column_dtypes[name] = dtype  # type: ignore[assignment]

    # apply schema overrides
    if schema_overrides:
        column_dtypes.update(schema_overrides)

    return column_names, column_dtypes


def _handle_columns_arg(
    data: list[PySeries],
    columns: Sequence[str] | None = None,
    *,
    from_dict: bool = False,
) -> list[PySeries]:
    """Rename data according to columns argument."""
    if columns is None:
        return data
    elif not data:
        return [pl.Series(name=c)._s for c in columns]
    elif len(data) != len(columns):
        msg = f"dimensions of columns arg ({len(columns)}) must match data dimensions ({len(data)})"
        raise ValueError(msg)

    if from_dict:
        series_map = {s.name(): s for s in data}
        if all((col in series_map) for col in columns):
            return [series_map[col] for col in columns]

    for i, c in enumerate(columns):
        if c != data[i].name():
            data[i] = data[i].clone()
            data[i].rename(c)

    return data


def _post_apply_columns(
    pydf: PyDataFrame,
    columns: SchemaDefinition | None,
    structs: dict[str, Struct] | None = None,
    schema_overrides: SchemaDict | None = None,
    *,
    strict: bool = True,
) -> PyDataFrame:
    """Apply 'columns' param *after* PyDataFrame creation (if no alternative)."""
    pydf_columns, pydf_dtypes = pydf.columns(), pydf.dtypes()
    columns, dtypes = _unpack_schema(
        (columns or pydf_columns), schema_overrides=schema_overrides
    )
    column_subset: list[str] = []
    if columns != pydf_columns:
        if len(columns) < len(pydf_columns) and columns == pydf_columns[: len(columns)]:
            column_subset = columns
        else:
            pydf.set_column_names(columns)

    column_casts = []
    for i, col in enumerate(columns):
        dtype = dtypes.get(col)
        pydf_dtype = pydf_dtypes[i]
        if dtype == Categorical != pydf_dtype:
            column_casts.append(F.col(col).cast(Categorical, strict=strict)._pyexpr)
        elif dtype == Enum != pydf_dtype:
            column_casts.append(F.col(col).cast(dtype, strict=strict)._pyexpr)
        elif structs and (struct := structs.get(col)) and struct != pydf_dtype:
            column_casts.append(F.col(col).cast(struct, strict=strict)._pyexpr)
        elif dtype is not None and dtype != Unknown and dtype != pydf_dtype:
            if dtype.is_temporal() and dtype != Duration and pydf_dtype == String:
                temporal_cast = F.col(col).str.strptime(dtype, strict=strict)._pyexpr  # type: ignore[arg-type]
                column_casts.append(temporal_cast)
            else:
                column_casts.append(F.col(col).cast(dtype, strict=strict)._pyexpr)

    if column_casts or column_subset:
        pyldf = pydf.lazy()
        if column_casts:
            pyldf = pyldf.with_columns(column_casts)
        if column_subset:
            pyldf = pyldf.select([F.col(col)._pyexpr for col in column_subset])
        pydf = pyldf.collect(engine="in-memory", lambda_post_opt=None)

    return pydf


def _expand_dict_values(
    data: Mapping[str, Sequence[object] | Mapping[str, Sequence[object]] | Series],
    *,
    schema_overrides: SchemaDict | None = None,
    strict: bool = True,
    order: Sequence[str] | None = None,
    nan_to_null: bool = False,
) -> dict[str, Series]:
    """Expand any scalar values in dict data (propagate literal as array)."""
    updated_data = {}
    if data:
        if any(isinstance(val, pl.Expr) for val in data.values()):
            msg = (
                "passing Expr objects to the DataFrame constructor is not supported"
                "\n\nHint: Try evaluating the expression first using `select`,"
                " or if you meant to create an Object column containing expressions,"
                " pass a list of Expr objects instead."
            )
            raise TypeError(msg)

        dtypes = schema_overrides or {}
        data = _expand_dict_data(data, dtypes, strict=strict)
        array_len = max((arrlen(val) or 0) for val in data.values())
        if array_len > 0:
            for name, val in data.items():
                dtype = dtypes.get(name)
                if isinstance(val, dict) and dtype != Struct:
                    vdf = pl.DataFrame(val, strict=strict)
                    if (
                        vdf.height == 1
                        and array_len > 1
                        and all(not d.is_nested() for d in vdf.schema.values())
                    ):
                        s_vals = {
                            nm: vdf[nm].extend_constant(v, n=(array_len - 1))
                            for nm, v in val.items()
                        }
                        st = pl.DataFrame(s_vals).to_struct(name)
                    else:
                        st = vdf.to_struct(name)
                    updated_data[name] = st

                elif isinstance(val, pl.Series):
                    s = val.rename(name) if name != val.name else val
                    if dtype and dtype != s.dtype:
                        s = s.cast(dtype, strict=strict)
                    updated_data[name] = s

                elif arrlen(val) is not None or _is_generator(val):
                    updated_data[name] = pl.Series(
                        name=name,
                        values=val,
                        dtype=dtype,
                        strict=strict,
                        nan_to_null=nan_to_null,
                    )
                elif val is None or isinstance(  # type: ignore[redundant-expr]
                    val, (int, float, str, bool, date, datetime, time, timedelta)
                ):
                    updated_data[name] = F.repeat(
                        val, array_len, dtype=dtype, eager=True
                    ).alias(name)
                else:
                    updated_data[name] = pl.Series(
                        name=name, values=[val] * array_len, dtype=dtype, strict=strict
                    )

        elif all((arrlen(val) == 0) for val in data.values()):
            for name, val in data.items():
                updated_data[name] = pl.Series(
                    name, values=val, dtype=dtypes.get(name), strict=strict
                )

        elif all((arrlen(val) is None) for val in data.values()):
            for name, val in data.items():
                updated_data[name] = pl.Series(
                    name,
                    values=(val if _is_generator(val) else [val]),
                    dtype=dtypes.get(name),
                    strict=strict,
                )
    if order and list(updated_data) != order:
        return {col: updated_data.pop(col) for col in order}
    return updated_data


def _expand_dict_data(
    data: Mapping[str, Sequence[object] | Mapping[str, Sequence[object]] | Series],
    dtypes: SchemaDict,
    *,
    strict: bool = True,
) -> Mapping[str, Sequence[object] | Mapping[str, Sequence[object]] | Series]:
    """
    Expand any unsized generators/iterators.

    (Note that `range` is sized, and will take a fast-path on Series init).
    """
    expanded_data = {}
    for name, val in data.items():
        expanded_data[name] = (
            pl.Series(name, val, dtypes.get(name), strict=strict)
            if _is_generator(val)
            else val
        )
    return expanded_data


def sequence_to_pydf(
    data: Sequence[Any],
    schema: SchemaDefinition | None = None,
    *,
    schema_overrides: SchemaDict | None = None,
    strict: bool = True,
    orient: Orientation | None = None,
    infer_schema_length: int | None = N_INFER_DEFAULT,
    nan_to_null: bool = False,
) -> PyDataFrame:
    """Construct a PyDataFrame from a sequence."""
    if not data:
        return dict_to_pydf({}, schema=schema, schema_overrides=schema_overrides)

    return _sequence_to_pydf_dispatcher(
        get_first_non_none(data),
        data=data,
        schema=schema,
        schema_overrides=schema_overrides,
        strict=strict,
        orient=orient,
        infer_schema_length=infer_schema_length,
        nan_to_null=nan_to_null,
    )


@singledispatch
def _sequence_to_pydf_dispatcher(
    first_element: Any,
    data: Sequence[Any],
    schema: SchemaDefinition | None,
    *,
    schema_overrides: SchemaDict | None,
    strict: bool = True,
    orient: Orientation | None,
    infer_schema_length: int | None,
    nan_to_null: bool = False,
) -> PyDataFrame:
    # note: ONLY python-native data should participate in singledispatch registration
    # via top-level decorators, otherwise we have to import the associated module.
    # third-party libraries (such as numpy/pandas) should be identified inline (below)
    # and THEN registered for dispatch (here) so as not to break lazy-loading behaviour.

    common_params = {
        "data": data,
        "schema": schema,
        "schema_overrides": schema_overrides,
        "strict": strict,
        "orient": orient,
        "infer_schema_length": infer_schema_length,
        "nan_to_null": nan_to_null,
    }
    to_pydf: Callable[..., PyDataFrame]
    register_with_singledispatch = True

    if isinstance(first_element, Generator):
        to_pydf = _sequence_of_sequence_to_pydf
        data = [list(row) for row in data]
        first_element = data[0]
        register_with_singledispatch = False

    elif isinstance(first_element, pl.Series):
        to_pydf = _sequence_of_series_to_pydf

    elif _check_for_numpy(first_element) and isinstance(first_element, np.ndarray):
        to_pydf = _sequence_of_numpy_to_pydf

    elif _check_for_pandas(first_element) and isinstance(
        first_element, (pd.Series, pd.Index, pd.DatetimeIndex)
    ):
        to_pydf = _sequence_of_pandas_to_pydf

    elif dataclasses.is_dataclass(first_element):
        to_pydf = _sequence_of_dataclasses_to_pydf

    elif is_pydantic_model(first_element):
        to_pydf = _sequence_of_pydantic_models_to_pydf

    elif is_sqlalchemy_row(first_element):
        to_pydf = _sequence_of_tuple_to_pydf

    elif isinstance(first_element, Sequence) and not isinstance(first_element, str):
        to_pydf = _sequence_of_sequence_to_pydf
    else:
        to_pydf = _sequence_of_elements_to_pydf

    if register_with_singledispatch:
        _sequence_to_pydf_dispatcher.register(type(first_element), to_pydf)

    common_params["first_element"] = first_element
    return to_pydf(**common_params)


@_sequence_to_pydf_dispatcher.register(list)
def _sequence_of_sequence_to_pydf(
    first_element: Sequence[Any] | np.ndarray[Any, Any],
    data: Sequence[Any],
    schema: SchemaDefinition | None,
    *,
    schema_overrides: SchemaDict | None,
    strict: bool,
    orient: Orientation | None,
    infer_schema_length: int | None,
    nan_to_null: bool = False,
) -> PyDataFrame:
    if orient is None:
        if schema is None:
            orient = "col"
        else:
            # Try to infer orientation from schema length and data dimensions
            is_row_oriented = (len(schema) == len(first_element)) and (
                len(schema) != len(data)
            )
            orient = "row" if is_row_oriented else "col"

            if is_row_oriented:
                issue_warning(
                    "Row orientation inferred during DataFrame construction."
                    ' Explicitly specify the orientation by passing `orient="row"` to silence this warning.',
                    DataOrientationWarning,
                )

    if orient == "row":
        column_names, schema_overrides = _unpack_schema(
            schema, schema_overrides=schema_overrides, n_expected=len(first_element)
        )
        local_schema_override = (
            _include_unknowns(schema_overrides, column_names)
            if schema_overrides
            else {}
        )

        unpack_nested = False
        for col, tp in local_schema_override.items():
            if tp in (Categorical, Enum):
                local_schema_override[col] = String
            elif not unpack_nested and (tp.base_type() in (Unknown, Struct)):
                unpack_nested = contains_nested(
                    getattr(first_element, col, None).__class__, is_namedtuple
                )

        if unpack_nested:
            dicts = [nt_unpack(d) for d in data]
            pydf = PyDataFrame.from_dicts(
                dicts,
                schema=None,
                schema_overrides=None,
                strict=strict,
                infer_schema_length=infer_schema_length,
            )
        else:
            pydf = PyDataFrame.from_rows(
                data,
                schema=local_schema_override or None,
                infer_schema_length=infer_schema_length,
            )
        if column_names or schema_overrides:
            pydf = _post_apply_columns(
                pydf, column_names, schema_overrides=schema_overrides, strict=strict
            )
        return pydf

    elif orient == "col":
        column_names, schema_overrides = _unpack_schema(
            schema, schema_overrides=schema_overrides, n_expected=len(data)
        )
        data_series: list[PySeries] = [
            pl.Series(
                column_names[i],
                element,
                dtype=schema_overrides.get(column_names[i]),
                strict=strict,
                nan_to_null=nan_to_null,
            )._s
            for i, element in enumerate(data)
        ]
        return PyDataFrame(data_series)

    else:
        msg = f"`orient` must be one of {{'col', 'row', None}}, got {orient!r}"
        raise ValueError(msg)


def _sequence_of_series_to_pydf(
    first_element: Series,
    data: Sequence[Any],
    schema: SchemaDefinition | None,
    *,
    schema_overrides: SchemaDict | None,
    strict: bool,
    **kwargs: Any,
) -> PyDataFrame:
    series_names = [s.name for s in data]
    column_names, schema_overrides = _unpack_schema(
        schema or series_names,
        schema_overrides=schema_overrides,
        n_expected=len(data),
    )
    data_series: list[PySeries] = []
    for i, s in enumerate(data):
        if not s.name:
            s = s.alias(column_names[i])
        new_dtype = schema_overrides.get(column_names[i])
        if new_dtype and new_dtype != s.dtype:
            s = s.cast(new_dtype, strict=strict, wrap_numerical=False)
        data_series.append(s._s)

    data_series = _handle_columns_arg(data_series, columns=column_names)
    return PyDataFrame(data_series)


@_sequence_to_pydf_dispatcher.register(tuple)
def _sequence_of_tuple_to_pydf(
    first_element: tuple[Any, ...],
    data: Sequence[Any],
    schema: SchemaDefinition | None,
    *,
    schema_overrides: SchemaDict | None,
    strict: bool,
    orient: Orientation | None,
    infer_schema_length: int | None,
    nan_to_null: bool = False,
) -> PyDataFrame:
    # infer additional meta information if namedtuple
    if is_namedtuple(first_element.__class__) or is_sqlalchemy_row(first_element):
        if schema is None:
            schema = first_element._fields  # type: ignore[attr-defined]
            annotations = getattr(first_element, "__annotations__", None)
            if annotations and len(annotations) == len(schema):
                schema = [
                    (name, try_parse_into_dtype(tp))
                    for name, tp in first_element.__annotations__.items()
                ]
        if orient is None:
            orient = "row"

    # ...then defer to generic sequence processing
    return _sequence_of_sequence_to_pydf(
        first_element,
        data=data,
        schema=schema,
        schema_overrides=schema_overrides,
        strict=strict,
        orient=orient,
        infer_schema_length=infer_schema_length,
        nan_to_null=nan_to_null,
    )


@_sequence_to_pydf_dispatcher.register(Mapping)
@_sequence_to_pydf_dispatcher.register(dict)
def _sequence_of_dict_to_pydf(
    first_element: dict[str, Any],
    data: Sequence[Any],
    schema: SchemaDefinition | None,
    *,
    schema_overrides: SchemaDict | None,
    strict: bool,
    infer_schema_length: int | None,
    **kwargs: Any,
) -> PyDataFrame:
    column_names, schema_overrides = _unpack_schema(
        schema, schema_overrides=schema_overrides
    )
    dicts_schema = (
        _include_unknowns(schema_overrides, column_names or list(schema_overrides))
        if column_names
        else None
    )

    pydf = PyDataFrame.from_dicts(
        data,
        dicts_schema,
        schema_overrides,
        strict=strict,
        infer_schema_length=infer_schema_length,
    )
    return pydf


@_sequence_to_pydf_dispatcher.register(str)
def _sequence_of_elements_to_pydf(
    first_element: Any,
    data: Sequence[Any],
    schema: SchemaDefinition | None,
    schema_overrides: SchemaDict | None,
    *,
    strict: bool,
    **kwargs: Any,
) -> PyDataFrame:
    column_names, schema_overrides = _unpack_schema(
        schema, schema_overrides=schema_overrides, n_expected=1
    )
    data_series: list[PySeries] = [
        pl.Series(
            column_names[0],
            data,
            schema_overrides.get(column_names[0]),
            strict=strict,
        )._s
    ]
    data_series = _handle_columns_arg(data_series, columns=column_names)
    return PyDataFrame(data_series)


def _sequence_of_numpy_to_pydf(
    first_element: np.ndarray[Any, Any],
    **kwargs: Any,
) -> PyDataFrame:
    if first_element.ndim == 1:
        return _sequence_of_sequence_to_pydf(first_element, **kwargs)
    else:
        return _sequence_of_elements_to_pydf(first_element, **kwargs)


def _sequence_of_pandas_to_pydf(
    first_element: pd.Series[Any] | pd.Index[Any] | pd.DatetimeIndex,
    data: Sequence[Any],
    schema: SchemaDefinition | None,
    schema_overrides: SchemaDict | None,
    *,
    strict: bool,
    **kwargs: Any,
) -> PyDataFrame:
    if schema is None:
        column_names: list[str] = []
    else:
        column_names, schema_overrides = _unpack_schema(
            schema, schema_overrides=schema_overrides, n_expected=1
        )

    schema_overrides = schema_overrides or {}
    data_series: list[PySeries] = []
    for i, s in enumerate(data):
        name = column_names[i] if column_names else s.name
        pyseries = plc.pandas_to_pyseries(name=name, values=s)
        dtype = schema_overrides.get(name)
        if dtype is not None and dtype != pyseries.dtype():
            pyseries = pyseries.cast(dtype, strict=strict, wrap_numerical=False)
        data_series.append(pyseries)

    return PyDataFrame(data_series)


def _sequence_of_dataclasses_to_pydf(
    first_element: Any,
    data: Sequence[Any],
    schema: SchemaDefinition | None,
    schema_overrides: SchemaDict | None,
    infer_schema_length: int | None,
    *,
    strict: bool = True,
    **kwargs: Any,
) -> PyDataFrame:
    """Initialize DataFrame from Python dataclasses."""
    from dataclasses import asdict, astuple

    (
        unpack_nested,
        column_names,
        schema_overrides,
        overrides,
    ) = _establish_dataclass_or_model_schema(
        first_element, schema, schema_overrides, model_fields=None
    )
    if unpack_nested:
        dicts = [asdict(md) for md in data]
        pydf = PyDataFrame.from_dicts(
            dicts,
            schema=None,
            schema_overrides=None,
            strict=strict,
            infer_schema_length=infer_schema_length,
        )
    else:
        rows = [astuple(dc) for dc in data]
        pydf = PyDataFrame.from_rows(
            rows,  # type: ignore[arg-type]
            schema=overrides or None,
            infer_schema_length=infer_schema_length,
        )

    if overrides:
        structs = {c: tp for c, tp in overrides.items() if isinstance(tp, Struct)}
        pydf = _post_apply_columns(
            pydf, column_names, structs, schema_overrides, strict=strict
        )

    return pydf


def _sequence_of_pydantic_models_to_pydf(
    first_element: Any,
    data: Sequence[Any],
    schema: SchemaDefinition | None,
    schema_overrides: SchemaDict | None,
    infer_schema_length: int | None,
    *,
    strict: bool,
    **kwargs: Any,
) -> PyDataFrame:
    """Initialise DataFrame from pydantic model objects."""
    import pydantic  # note: must already be available in the env here

    old_pydantic = parse_version(pydantic.__version__) < (2, 0)
    model_fields = list(
        first_element.__fields__
        if old_pydantic
        else first_element.__class__.model_fields
    )
    (
        unpack_nested,
        column_names,
        schema_overrides,
        overrides,
    ) = _establish_dataclass_or_model_schema(
        first_element, schema, schema_overrides, model_fields
    )
    if unpack_nested:
        # note: this is an *extremely* slow path, due to the requirement to
        # use pydantic's 'dict()' method to properly unpack nested models
        dicts = (
            [md.dict() for md in data]
            if old_pydantic
            else [md.model_dump(mode="python") for md in data]
        )
        pydf = PyDataFrame.from_dicts(
            dicts,
            schema=None,
            schema_overrides=None,
            strict=strict,
            infer_schema_length=infer_schema_length,
        )

    elif len(model_fields) > 50:
        # 'from_rows' is the faster codepath for models with a lot of fields...
        get_values = itemgetter(*model_fields)
        rows = [get_values(md.__dict__) for md in data]
        pydf = PyDataFrame.from_rows(
            rows, schema=overrides, infer_schema_length=infer_schema_length
        )
    else:
        # ...and 'from_dicts' is faster otherwise
        dicts = [md.__dict__ for md in data]
        pydf = PyDataFrame.from_dicts(
            dicts,
            schema=overrides,
            schema_overrides=None,
            strict=strict,
            infer_schema_length=infer_schema_length,
        )

    if overrides:
        structs = {c: tp for c, tp in overrides.items() if isinstance(tp, Struct)}
        pydf = _post_apply_columns(
            pydf, column_names, structs, schema_overrides, strict=strict
        )

    return pydf


def _establish_dataclass_or_model_schema(
    first_element: Any,
    schema: SchemaDefinition | None,
    schema_overrides: SchemaDict | None,
    model_fields: list[str] | None,
) -> tuple[bool, list[str], SchemaDict, SchemaDict]:
    """Shared utility code for establishing dataclasses/pydantic model cols/schema."""
    from dataclasses import asdict

    unpack_nested = False
    if schema:
        column_names, schema_overrides = _unpack_schema(
            schema, schema_overrides=schema_overrides
        )
        overrides = {col: schema_overrides.get(col, Unknown) for col in column_names}
    else:
        column_names = []
        overrides = {
            col: (try_parse_into_dtype(tp) or Unknown)
            for col, tp in try_get_type_hints(first_element.__class__).items()
            if ((col in model_fields) if model_fields else (col != "__slots__"))
        }
        if schema_overrides:
            overrides.update(schema_overrides)
        elif not model_fields:
            dc_fields = set(asdict(first_element))
            schema_overrides = overrides = {
                nm: tp for nm, tp in overrides.items() if nm in dc_fields
            }
        else:
            schema_overrides = overrides

    for col, tp in overrides.items():
        if tp in (Categorical, Enum):
            overrides[col] = String
        elif not unpack_nested and (tp.base_type() in (Unknown, Struct)):
            unpack_nested = contains_nested(
                getattr(first_element, col, None),
                is_pydantic_model if model_fields else dataclasses.is_dataclass,  # type: ignore[arg-type]
            )

    if model_fields and len(model_fields) == len(overrides):
        overrides = dict(zip(model_fields, overrides.values()))

    return unpack_nested, column_names, schema_overrides, overrides


def _include_unknowns(
    schema: SchemaDict, cols: Sequence[str]
) -> MutableMapping[str, PolarsDataType]:
    """Complete partial schema dict by including Unknown type."""
    return {
        col: (schema.get(col, Unknown) or Unknown)  # type: ignore[truthy-bool]
        for col in cols
    }


def iterable_to_pydf(
    data: Iterable[Any],
    schema: SchemaDefinition | None = None,
    *,
    schema_overrides: SchemaDict | None = None,
    strict: bool = True,
    orient: Orientation | None = None,
    chunk_size: int | None = None,
    infer_schema_length: int | None = N_INFER_DEFAULT,
    rechunk: bool = True,
) -> PyDataFrame:
    """Construct a PyDataFrame from an iterable/generator."""
    original_schema = schema
    column_names: list[str] = []
    dtypes_by_idx: dict[int, PolarsDataType] = {}
    if schema is not None:
        column_names, schema_overrides = _unpack_schema(
            schema, schema_overrides=schema_overrides
        )
    elif schema_overrides:
        _, schema_overrides = _unpack_schema(schema, schema_overrides=schema_overrides)

    if not isinstance(data, Generator):
        data = iter(data)

    if orient == "col":
        if column_names and schema_overrides:
            dtypes_by_idx = {
                idx: schema_overrides.get(col, Unknown)
                for idx, col in enumerate(column_names)
            }

        return pl.DataFrame(
            {
                (column_names[idx] if column_names else f"column_{idx}"): pl.Series(
                    coldata,
                    dtype=dtypes_by_idx.get(idx),
                    strict=strict,
                )
                for idx, coldata in enumerate(data)
            },
        )._df

    def to_frame_chunk(values: list[Any], schema: SchemaDefinition | None) -> DataFrame:
        return pl.DataFrame(
            data=values,
            schema=schema,
            strict=strict,
            orient="row",
            infer_schema_length=infer_schema_length,
            schema_overrides=schema_overrides,
        )

    n_chunks = 0
    n_chunk_elems = 1_000_000

    if chunk_size:
        adaptive_chunk_size = chunk_size
    elif column_names:
        adaptive_chunk_size = n_chunk_elems // len(column_names)
    else:
        adaptive_chunk_size = None

    df: DataFrame = None  # type: ignore[assignment]
    chunk_size = (
        None
        if infer_schema_length is None
        else max(infer_schema_length, adaptive_chunk_size or 1000)
    )
    while True:
        values = list(islice(data, chunk_size))
        if not values:
            break
        frame_chunk = to_frame_chunk(values, original_schema)
        if df is None:
            df = frame_chunk
            if not original_schema:
                original_schema = list(df.schema.items())
            if chunk_size != adaptive_chunk_size:
                if (n_columns := df.width) > 0:
                    chunk_size = adaptive_chunk_size = n_chunk_elems // n_columns
        else:
            df.vstack(frame_chunk, in_place=True)
            n_chunks += 1

    if df is None:
        df = to_frame_chunk([], original_schema)

    if n_chunks > 0 and rechunk:
        df = df.rechunk()

    return df._df


def _check_pandas_columns(data: pd.DataFrame, *, include_index: bool) -> None:
    """Check pandas dataframe columns can be converted to polars."""
    stringified_cols: set[str] = {str(col) for col in data.columns}
    stringified_index: set[str] = (
        {str(idx) for idx in data.index.names} if include_index else set()
    )

    non_unique_cols: bool = len(stringified_cols) < len(data.columns)
    non_unique_indices: bool = (
        (len(stringified_index) < len(data.index.names)) if include_index else False
    )
    if non_unique_cols or non_unique_indices:
        msg = (
            "Pandas dataframe contains non-unique indices and/or column names. "
            "Polars dataframes require unique string names for columns."
        )
        raise ValueError(msg)

    overlapping_cols_and_indices: set[str] = stringified_cols & stringified_index
    if len(overlapping_cols_and_indices) > 0:
        msg = "Pandas indices and column names must not overlap."
        raise ValueError(msg)


def pandas_to_pydf(
    data: pd.DataFrame,
    schema: SchemaDefinition | None = None,
    *,
    schema_overrides: SchemaDict | None = None,
    strict: bool = True,
    rechunk: bool = True,
    nan_to_null: bool = True,
    include_index: bool = False,
) -> PyDataFrame:
    """Construct a PyDataFrame from a pandas DataFrame."""
    _check_pandas_columns(data, include_index=include_index)

    convert_index = include_index and not _pandas_has_default_index(data)
    if not convert_index and all(
        is_simple_numpy_backed_pandas_series(data[col]) for col in data.columns
    ):
        # Convert via NumPy directly, no PyArrow needed.
        return pl.DataFrame(
            {str(col): data[col].to_numpy() for col in data.columns},
            schema=schema,
            strict=strict,
            schema_overrides=schema_overrides,
            nan_to_null=nan_to_null,
        )._df

    if not _PYARROW_AVAILABLE:
        msg = (
            "pyarrow is required for converting a pandas dataframe to Polars, "
            "unless each of its columns is a simple numpy-backed one "
            "(e.g. 'int64', 'bool', 'float32' - not 'Int64')"
        )
        raise ImportError(msg)
    arrow_dict = {}
    length = data.shape[0]

    if convert_index:
        for idxcol in data.index.names:
            arrow_dict[str(idxcol)] = plc.pandas_series_to_arrow(
                # get_level_values accepts `int | str`
                # but `index.names` returns `Hashable`
                data.index.get_level_values(idxcol),  # type: ignore[arg-type, unused-ignore]
                nan_to_null=nan_to_null,
                length=length,
            )

    for col_idx, col_data in data.items():
        arrow_dict[str(col_idx)] = plc.pandas_series_to_arrow(
            col_data, nan_to_null=nan_to_null, length=length
        )

    arrow_table = pa.table(arrow_dict)
    return arrow_to_pydf(
        arrow_table,
        schema=schema,
        schema_overrides=schema_overrides,
        strict=strict,
        rechunk=rechunk,
    )


def _pandas_has_default_index(df: pd.DataFrame) -> bool:
    """Identify if the pandas frame only has a default (or equivalent) index."""
    from pandas.core.indexes.range import RangeIndex

    index_cols = df.index.names

    if len(index_cols) > 1 or index_cols not in ([None], [""]):
        # not default: more than one index, or index is named
        return False
    elif df.index.equals(RangeIndex(start=0, stop=len(df), step=1)):
        # is default: simple range index
        return True
    else:
        # finally, is the index _equivalent_ to a default unnamed
        # integer index with frame data that was previously sorted
        return (
            str(df.index.dtype).startswith("int")
            and (df.index.sort_values() == np.arange(len(df))).all()
        )


def arrow_to_pydf(
    data: pa.Table | pa.RecordBatch,
    schema: SchemaDefinition | None = None,
    *,
    schema_overrides: SchemaDict | None = None,
    strict: bool = True,
    rechunk: bool = True,
) -> PyDataFrame:
    """Construct a PyDataFrame from an Arrow Table or RecordBatch."""
    column_names, schema_overrides = _unpack_schema(
        (schema or data.schema.names), schema_overrides=schema_overrides
    )
    try:
        if column_names != data.schema.names:
            data = data.rename_columns(column_names)
    except pa.ArrowInvalid as e:
        msg = "dimensions of columns arg must match data dimensions"
        raise ValueError(msg) from e

    batches: list[pa.RecordBatch]
    if isinstance(data, pa.RecordBatch):
        batches = [data]
    else:
        batches = data.to_batches()

    # supply the arrow schema so the metadata is intact
    pydf = PyDataFrame.from_arrow_record_batches(batches, data.schema)

    if rechunk:
        pydf = pydf.rechunk()

    if schema_overrides is not None:
        pydf = _post_apply_columns(
            pydf,
            column_names,
            schema_overrides=schema_overrides,
            strict=strict,
        )

    return pydf


def numpy_to_pydf(
    data: np.ndarray[Any, Any],
    schema: SchemaDefinition | None = None,
    *,
    schema_overrides: SchemaDict | None = None,
    orient: Orientation | None = None,
    strict: bool = True,
    nan_to_null: bool = False,
) -> PyDataFrame:
    """Construct a PyDataFrame from a NumPy ndarray (including structured ndarrays)."""
    shape = data.shape
    two_d = len(shape) == 2

    if data.dtype.names is not None:
        structured_array, orient = True, "col"
        record_names = list(data.dtype.names)
        n_columns = len(record_names)
        for nm in record_names:
            shape = data[nm].shape
        if not schema:
            schema = record_names
    else:
        # Unpack columns
        structured_array, record_names = False, []
        if shape == (0,):
            n_columns = 0

        elif len(shape) == 1:
            n_columns = 1

        elif len(shape) == 2:
            if orient is None and schema is None:
                # default convention; first axis is rows, second axis is columns
                n_columns = shape[1]
                orient = "row"

            elif orient is None and schema is not None:
                # infer orientation from 'schema' param; if square array
                # we check the flags to establish row/column major order
                n_schema_cols = len(schema)
                if n_schema_cols == shape[0] and n_schema_cols != shape[1]:
                    orient = "col"
                    n_columns = shape[0]
                elif data.flags["F_CONTIGUOUS"] and shape[0] == shape[1]:
                    orient = "col"
                    n_columns = n_schema_cols
                else:
                    orient = "row"
                    n_columns = shape[1]

            elif orient == "row":
                n_columns = shape[1]
            elif orient == "col":
                n_columns = shape[0]
            else:
                msg = f"`orient` must be one of {{'col', 'row', None}}, got {orient!r}"
                raise ValueError(msg)
        else:
            if shape == ():
                msg = "cannot create DataFrame from zero-dimensional array"
            else:
                msg = f"cannot create DataFrame from array with more than two dimensions; shape = {shape}"
            raise ValueError(msg)

    if schema is not None and len(schema) != n_columns:
        if (n_schema_cols := len(schema)) != 1:
            msg = f"dimensions of `schema` ({n_schema_cols}) must match data dimensions ({n_columns})"
            raise ValueError(msg)
        n_columns = n_schema_cols

    column_names, schema_overrides = _unpack_schema(
        schema, schema_overrides=schema_overrides, n_expected=n_columns
    )

    # Convert data to series
    if structured_array:
        data_series = [
            pl.Series(
                name=series_name,
                values=data[record_name],
                dtype=schema_overrides.get(record_name),
                strict=strict,
                nan_to_null=nan_to_null,
            )._s
            for series_name, record_name in zip(column_names, record_names)
        ]
    elif shape == (0,) and n_columns == 0:
        data_series = []

    elif len(shape) == 1:
        data_series = [
            pl.Series(
                name=column_names[0],
                values=data,
                dtype=schema_overrides.get(column_names[0]),
                strict=strict,
                nan_to_null=nan_to_null,
            )._s
        ]
    else:
        if orient == "row":
            data_series = [
                pl.Series(
                    name=column_names[i],
                    values=(
                        data
                        if two_d and n_columns == 1 and shape[1] > 1
                        else data[:, i]
                    ),
                    dtype=schema_overrides.get(column_names[i]),
                    strict=strict,
                    nan_to_null=nan_to_null,
                )._s
                for i in range(n_columns)
            ]
        else:
            data_series = [
                pl.Series(
                    name=column_names[i],
                    values=(
                        data if two_d and n_columns == 1 and shape[1] > 1 else data[i]
                    ),
                    dtype=schema_overrides.get(column_names[i]),
                    strict=strict,
                    nan_to_null=nan_to_null,
                )._s
                for i in range(n_columns)
            ]

    data_series = _handle_columns_arg(data_series, columns=column_names)
    return PyDataFrame(data_series)


def series_to_pydf(
    data: Series,
    schema: SchemaDefinition | None = None,
    schema_overrides: SchemaDict | None = None,
    *,
    strict: bool = True,
) -> PyDataFrame:
    """Construct a PyDataFrame from a Polars Series."""
    if schema is None and schema_overrides is None:
        return PyDataFrame([data._s])

    data_series = [data._s]
    series_name = [s.name() for s in data_series]
    column_names, schema_overrides = _unpack_schema(
        schema or series_name, schema_overrides=schema_overrides, n_expected=1
    )
    if schema_overrides:
        new_dtype = next(iter(schema_overrides.values()))
        if new_dtype != data.dtype:
            data_series[0] = data_series[0].cast(
                new_dtype, strict=strict, wrap_numerical=False
            )

    data_series = _handle_columns_arg(data_series, columns=column_names)
    return PyDataFrame(data_series)


def dataframe_to_pydf(
    data: DataFrame,
    schema: SchemaDefinition | None = None,
    *,
    schema_overrides: SchemaDict | None = None,
    strict: bool = True,
) -> PyDataFrame:
    """Construct a PyDataFrame from an existing Polars DataFrame."""
    if schema is None and schema_overrides is None:
        return data._df.clone()

    data_series = {c.name: c._s for c in data}
    column_names, schema_overrides = _unpack_schema(
        schema or data.columns, schema_overrides=schema_overrides
    )
    if schema_overrides:
        existing_schema = data.schema
        for name, new_dtype in schema_overrides.items():
            if new_dtype != existing_schema[name]:
                data_series[name] = data_series[name].cast(
                    new_dtype, strict=strict, wrap_numerical=False
                )

    series_cols = _handle_columns_arg(list(data_series.values()), columns=column_names)
    return PyDataFrame(series_cols)
