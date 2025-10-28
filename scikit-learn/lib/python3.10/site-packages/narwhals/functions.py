from __future__ import annotations

import platform
import sys
from collections.abc import Iterable, Mapping, Sequence
from functools import partial
from typing import TYPE_CHECKING, Any, Callable

from narwhals._expression_parsing import (
    ExprKind,
    ExprMetadata,
    apply_n_ary_operation,
    combine_metadata,
    is_scalar_like,
)
from narwhals._utils import (
    Implementation,
    Version,
    deprecate_native_namespace,
    flatten,
    is_compliant_expr,
    is_eager_allowed,
    is_sequence_but_not_str,
    normalize_path,
    supports_arrow_c_stream,
    validate_laziness,
)
from narwhals.dependencies import (
    is_narwhals_series,
    is_numpy_array,
    is_numpy_array_2d,
    is_pyarrow_table,
)
from narwhals.exceptions import InvalidOperationError
from narwhals.expr import Expr
from narwhals.series import Series
from narwhals.translate import from_native, to_native

if TYPE_CHECKING:
    from types import ModuleType

    from typing_extensions import TypeAlias, TypeIs

    from narwhals._compliant import CompliantExpr, CompliantNamespace
    from narwhals._native import NativeDataFrame, NativeLazyFrame, NativeSeries
    from narwhals._translate import IntoArrowTable
    from narwhals._typing import Backend, EagerAllowed, IntoBackend
    from narwhals.dataframe import DataFrame, LazyFrame
    from narwhals.typing import (
        ConcatMethod,
        FileSource,
        FrameT,
        IntoDType,
        IntoExpr,
        IntoSchema,
        NonNestedLiteral,
        _1DArray,
        _2DArray,
    )

    _IntoSchema: TypeAlias = "IntoSchema | Sequence[str] | None"


def concat(items: Iterable[FrameT], *, how: ConcatMethod = "vertical") -> FrameT:
    """Concatenate multiple DataFrames, LazyFrames into a single entity.

    Arguments:
        items: DataFrames, LazyFrames to concatenate.
        how: concatenating strategy

            - vertical: Concatenate vertically. Column names must match.
            - horizontal: Concatenate horizontally. If lengths don't match, then
                missing rows are filled with null values. This is only supported
                when all inputs are (eager) DataFrames.
            - diagonal: Finds a union between the column schemas and fills missing column
                values with null.

    Raises:
        TypeError: The items to concatenate should either all be eager, or all lazy

    Examples:
        Let's take an example of vertical concatenation:

        >>> import pandas as pd
        >>> import polars as pl
        >>> import pyarrow as pa
        >>> import narwhals as nw

        Let's look at one case a for vertical concatenation (pandas backed):

        >>> df_pd_1 = nw.from_native(pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]}))
        >>> df_pd_2 = nw.from_native(pd.DataFrame({"a": [5, 2], "b": [1, 4]}))
        >>> nw.concat([df_pd_1, df_pd_2], how="vertical")
        ┌──────────────────┐
        |Narwhals DataFrame|
        |------------------|
        |        a  b      |
        |     0  1  4      |
        |     1  2  5      |
        |     2  3  6      |
        |     0  5  1      |
        |     1  2  4      |
        └──────────────────┘

        Let's look at one case a for horizontal concatenation (polars backed):

        >>> df_pl_1 = nw.from_native(pl.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]}))
        >>> df_pl_2 = nw.from_native(pl.DataFrame({"c": [5, 2], "d": [1, 4]}))
        >>> nw.concat([df_pl_1, df_pl_2], how="horizontal")
        ┌───────────────────────────┐
        |    Narwhals DataFrame     |
        |---------------------------|
        |shape: (3, 4)              |
        |┌─────┬─────┬──────┬──────┐|
        |│ a   ┆ b   ┆ c    ┆ d    │|
        |│ --- ┆ --- ┆ ---  ┆ ---  │|
        |│ i64 ┆ i64 ┆ i64  ┆ i64  │|
        |╞═════╪═════╪══════╪══════╡|
        |│ 1   ┆ 4   ┆ 5    ┆ 1    │|
        |│ 2   ┆ 5   ┆ 2    ┆ 4    │|
        |│ 3   ┆ 6   ┆ null ┆ null │|
        |└─────┴─────┴──────┴──────┘|
        └───────────────────────────┘

        Let's look at one case a for diagonal concatenation (pyarrow backed):

        >>> df_pa_1 = nw.from_native(pa.table({"a": [1, 2], "b": [3.5, 4.5]}))
        >>> df_pa_2 = nw.from_native(pa.table({"a": [3, 4], "z": ["x", "y"]}))
        >>> nw.concat([df_pa_1, df_pa_2], how="diagonal")
        ┌──────────────────────────┐
        |    Narwhals DataFrame    |
        |--------------------------|
        |pyarrow.Table             |
        |a: int64                  |
        |b: double                 |
        |z: string                 |
        |----                      |
        |a: [[1,2],[3,4]]          |
        |b: [[3.5,4.5],[null,null]]|
        |z: [[null,null],["x","y"]]|
        └──────────────────────────┘
    """
    from narwhals.dependencies import is_narwhals_lazyframe

    if not items:
        msg = "No items to concatenate."
        raise ValueError(msg)
    items = tuple(items)
    validate_laziness(items)
    if how not in {"horizontal", "vertical", "diagonal"}:  # pragma: no cover
        msg = "Only vertical, horizontal and diagonal concatenations are supported."
        raise NotImplementedError(msg)
    first_item = items[0]
    if is_narwhals_lazyframe(first_item) and how == "horizontal":
        msg = (
            "Horizontal concatenation is not supported for LazyFrames.\n\n"
            "Hint: you may want to use `join` instead."
        )
        raise InvalidOperationError(msg)
    plx = first_item.__narwhals_namespace__()
    return first_item._with_compliant(
        plx.concat([df._compliant_frame for df in items], how=how)
    )


def new_series(
    name: str,
    values: Any,
    dtype: IntoDType | None = None,
    *,
    backend: IntoBackend[EagerAllowed],
) -> Series[Any]:
    """Instantiate Narwhals Series from iterable (e.g. list or array).

    Arguments:
        name: Name of resulting Series.
        values: Values of make Series from.
        dtype: (Narwhals) dtype. If not provided, the native library
            may auto-infer it from `values`.
        backend: specifies which eager backend instantiate to.

            `backend` can be specified in various ways

            - As `Implementation.<BACKEND>` with `BACKEND` being `PANDAS`, `PYARROW`,
                `POLARS`, `MODIN` or `CUDF`.
            - As a string: `"pandas"`, `"pyarrow"`, `"polars"`, `"modin"` or `"cudf"`.
            - Directly as a module `pandas`, `pyarrow`, `polars`, `modin` or `cudf`.

    Examples:
        >>> import pandas as pd
        >>> import narwhals as nw
        >>>
        >>> values = [4, 1, 2, 3]
        >>> nw.new_series(name="a", values=values, dtype=nw.Int32, backend=pd)
        ┌─────────────────────┐
        |   Narwhals Series   |
        |---------------------|
        |0    4               |
        |1    1               |
        |2    2               |
        |3    3               |
        |Name: a, dtype: int32|
        └─────────────────────┘
    """
    return _new_series_impl(name, values, dtype, backend=backend)


def _new_series_impl(
    name: str,
    values: Any,
    dtype: IntoDType | None = None,
    *,
    backend: IntoBackend[EagerAllowed],
) -> Series[Any]:
    implementation = Implementation.from_backend(backend)
    if is_eager_allowed(implementation):
        ns = Version.MAIN.namespace.from_backend(implementation).compliant
        series = ns._series.from_iterable(values, name=name, context=ns, dtype=dtype)
        return series.to_narwhals()
    if implementation is Implementation.UNKNOWN:  # pragma: no cover
        _native_namespace = implementation.to_native_namespace()
        try:
            native_series: NativeSeries = _native_namespace.new_series(
                name, values, dtype
            )
            return from_native(native_series, series_only=True).alias(name)
        except AttributeError as e:
            msg = "Unknown namespace is expected to implement `new_series` constructor."
            raise AttributeError(msg) from e
    msg = (
        f"{implementation} support in Narwhals is lazy-only, but `new_series` is an eager-only function.\n\n"
        "Hint: you may want to use an eager backend and then call `.lazy`, e.g.:\n\n"
        f"    nw.new_series('a', [1,2,3], backend='pyarrow').to_frame().lazy('{implementation}')"
    )
    raise ValueError(msg)


@deprecate_native_namespace(warn_version="1.26.0")
def from_dict(
    data: Mapping[str, Any],
    schema: IntoSchema | None = None,
    *,
    backend: IntoBackend[EagerAllowed] | None = None,
    native_namespace: ModuleType | None = None,  # noqa: ARG001
) -> DataFrame[Any]:
    """Instantiate DataFrame from dictionary.

    Indexes (if present, for pandas-like backends) are aligned following
    the [left-hand-rule](../concepts/pandas_index.md/).

    Notes:
        For pandas-like dataframes, conversion to schema is applied after dataframe
        creation.

    Arguments:
        data: Dictionary to create DataFrame from.
        schema: The DataFrame schema as Schema or dict of {name: type}. If not
            specified, the schema will be inferred by the native library.
        backend: specifies which eager backend instantiate to. Only
            necessary if inputs are not Narwhals Series.

            `backend` can be specified in various ways

            - As `Implementation.<BACKEND>` with `BACKEND` being `PANDAS`, `PYARROW`,
                `POLARS`, `MODIN` or `CUDF`.
            - As a string: `"pandas"`, `"pyarrow"`, `"polars"`, `"modin"` or `"cudf"`.
            - Directly as a module `pandas`, `pyarrow`, `polars`, `modin` or `cudf`.
        native_namespace: deprecated, same as `backend`.

    Examples:
        >>> import pandas as pd
        >>> import narwhals as nw
        >>> data = {"c": [5, 2], "d": [1, 4]}
        >>> nw.from_dict(data, backend="pandas")
        ┌──────────────────┐
        |Narwhals DataFrame|
        |------------------|
        |        c  d      |
        |     0  5  1      |
        |     1  2  4      |
        └──────────────────┘
    """
    if backend is None:
        data, backend = _from_dict_no_backend(data)
    implementation = Implementation.from_backend(backend)
    if is_eager_allowed(implementation):
        ns = Version.MAIN.namespace.from_backend(implementation).compliant
        return ns._dataframe.from_dict(data, schema=schema, context=ns).to_narwhals()
    if implementation is Implementation.UNKNOWN:  # pragma: no cover
        _native_namespace = implementation.to_native_namespace()
        try:
            # implementation is UNKNOWN, Narwhals extension using this feature should
            # implement `from_dict` function in the top-level namespace.
            native_frame: NativeDataFrame = _native_namespace.from_dict(
                data, schema=schema
            )
        except AttributeError as e:
            msg = "Unknown namespace is expected to implement `from_dict` function."
            raise AttributeError(msg) from e
        return from_native(native_frame, eager_only=True)
    msg = (
        f"{implementation} support in Narwhals is lazy-only, but `from_dict` is an eager-only function.\n\n"
        "Hint: you may want to use an eager backend and then call `.lazy`, e.g.:\n\n"
        f"    nw.from_dict({{'a': [1, 2]}}, backend='pyarrow').lazy('{implementation}')"
    )
    raise ValueError(msg)


def _from_dict_no_backend(
    data: Mapping[str, Series[Any] | Any], /
) -> tuple[dict[str, Series[Any] | Any], ModuleType]:
    for val in data.values():
        if is_narwhals_series(val):
            native_namespace = val.__native_namespace__()
            break
    else:
        msg = "Calling `from_dict` without `backend` is only supported if all input values are already Narwhals Series"
        raise TypeError(msg)
    data = {key: to_native(value, pass_through=True) for key, value in data.items()}
    return data, native_namespace


def from_dicts(
    data: Sequence[Mapping[str, Any]],
    schema: IntoSchema | None = None,
    *,
    backend: IntoBackend[EagerAllowed],
) -> DataFrame[Any]:
    """Instantiate DataFrame from a sequence of dictionaries representing rows.

    Notes:
        For pandas-like dataframes, conversion to schema is applied after dataframe
        creation.

    Arguments:
        data: Sequence with dictionaries mapping column name to value.
        schema: The DataFrame schema as Schema or dict of {name: type}. If not
            specified, the schema will be inferred by the native library.
        backend: Specifies which eager backend instantiate to.

            `backend` can be specified in various ways

            - As `Implementation.<BACKEND>` with `BACKEND` being `PANDAS`, `PYARROW`,
                `POLARS`, `MODIN` or `CUDF`.
            - As a string: `"pandas"`, `"pyarrow"`, `"polars"`, `"modin"` or `"cudf"`.
            - Directly as a module `pandas`, `pyarrow`, `polars`, `modin` or `cudf`.

    Tip:
        If you expect non-uniform keys in `data`, consider passing `schema` for
        more consistent results, as **inference varies between backends**:

        - pandas uses all rows
        - polars uses the first 100 rows
        - pyarrow uses only the first row

    Examples:
        >>> import polars as pl
        >>> import narwhals as nw
        >>> data = [
        ...     {"item": "apple", "weight": 80, "price": 0.60},
        ...     {"item": "egg", "weight": 55, "price": 0.40},
        ... ]
        >>> nw.DataFrame.from_dicts(data, backend="polars")
        ┌──────────────────────────┐
        |    Narwhals DataFrame    |
        |--------------------------|
        |shape: (2, 3)             |
        |┌───────┬────────┬───────┐|
        |│ item  ┆ weight ┆ price │|
        |│ ---   ┆ ---    ┆ ---   │|
        |│ str   ┆ i64    ┆ f64   │|
        |╞═══════╪════════╪═══════╡|
        |│ apple ┆ 80     ┆ 0.6   │|
        |│ egg   ┆ 55     ┆ 0.4   │|
        |└───────┴────────┴───────┘|
        └──────────────────────────┘
    """
    return Version.MAIN.dataframe.from_dicts(data, schema, backend=backend)


def from_numpy(
    data: _2DArray,
    schema: IntoSchema | Sequence[str] | None = None,
    *,
    backend: IntoBackend[EagerAllowed],
) -> DataFrame[Any]:
    """Construct a DataFrame from a NumPy ndarray.

    Notes:
        Only row orientation is currently supported.

        For pandas-like dataframes, conversion to schema is applied after dataframe
        creation.

    Arguments:
        data: Two-dimensional data represented as a NumPy ndarray.
        schema: The DataFrame schema as Schema, dict of {name: type}, or a sequence of str.
        backend: specifies which eager backend instantiate to.

            `backend` can be specified in various ways

            - As `Implementation.<BACKEND>` with `BACKEND` being `PANDAS`, `PYARROW`,
                `POLARS`, `MODIN` or `CUDF`.
            - As a string: `"pandas"`, `"pyarrow"`, `"polars"`, `"modin"` or `"cudf"`.
            - Directly as a module `pandas`, `pyarrow`, `polars`, `modin` or `cudf`.

    Examples:
        >>> import numpy as np
        >>> import pyarrow as pa
        >>> import narwhals as nw
        >>>
        >>> arr = np.array([[5, 2, 1], [1, 4, 3]])
        >>> schema = {"c": nw.Int16(), "d": nw.Float32(), "e": nw.Int8()}
        >>> nw.from_numpy(arr, schema=schema, backend="pyarrow")
        ┌──────────────────┐
        |Narwhals DataFrame|
        |------------------|
        |  pyarrow.Table   |
        |  c: int16        |
        |  d: float        |
        |  e: int8         |
        |  ----            |
        |  c: [[5,1]]      |
        |  d: [[2,4]]      |
        |  e: [[1,3]]      |
        └──────────────────┘
    """
    if not is_numpy_array_2d(data):
        msg = "`from_numpy` only accepts 2D numpy arrays"
        raise ValueError(msg)
    if not _is_into_schema(schema):
        msg = (
            "`schema` is expected to be one of the following types: "
            "IntoSchema | Sequence[str]. "
            f"Got {type(schema)}."
        )
        raise TypeError(msg)
    implementation = Implementation.from_backend(backend)
    if is_eager_allowed(implementation):
        ns = Version.MAIN.namespace.from_backend(implementation).compliant
        return ns.from_numpy(data, schema).to_narwhals()
    if implementation is Implementation.UNKNOWN:  # pragma: no cover
        _native_namespace = implementation.to_native_namespace()
        try:
            # implementation is UNKNOWN, Narwhals extension using this feature should
            # implement `from_numpy` function in the top-level namespace.
            native_frame: NativeDataFrame = _native_namespace.from_numpy(
                data, schema=schema
            )
        except AttributeError as e:
            msg = "Unknown namespace is expected to implement `from_numpy` function."
            raise AttributeError(msg) from e
        return from_native(native_frame, eager_only=True)
    msg = (
        f"{implementation} support in Narwhals is lazy-only, but `from_numpy` is an eager-only function.\n\n"
        "Hint: you may want to use an eager backend and then call `.lazy`, e.g.:\n\n"
        f"    nw.from_numpy(arr, backend='pyarrow').lazy('{implementation}')"
    )
    raise ValueError(msg)


def _is_into_schema(obj: Any) -> TypeIs[_IntoSchema]:
    from narwhals.schema import Schema

    return (
        obj is None or isinstance(obj, (Mapping, Schema)) or is_sequence_but_not_str(obj)
    )


def from_arrow(
    native_frame: IntoArrowTable, *, backend: IntoBackend[EagerAllowed]
) -> DataFrame[Any]:  # pragma: no cover
    """Construct a DataFrame from an object which supports the PyCapsule Interface.

    Arguments:
        native_frame: Object which implements `__arrow_c_stream__`.
        backend: specifies which eager backend instantiate to.

            `backend` can be specified in various ways

            - As `Implementation.<BACKEND>` with `BACKEND` being `PANDAS`, `PYARROW`,
                `POLARS`, `MODIN` or `CUDF`.
            - As a string: `"pandas"`, `"pyarrow"`, `"polars"`, `"modin"` or `"cudf"`.
            - Directly as a module `pandas`, `pyarrow`, `polars`, `modin` or `cudf`.

    Examples:
        >>> import pandas as pd
        >>> import polars as pl
        >>> import narwhals as nw
        >>>
        >>> df_native = pd.DataFrame({"a": [1, 2], "b": [4.2, 5.1]})
        >>> nw.from_arrow(df_native, backend="polars")
        ┌──────────────────┐
        |Narwhals DataFrame|
        |------------------|
        |  shape: (2, 2)   |
        |  ┌─────┬─────┐   |
        |  │ a   ┆ b   │   |
        |  │ --- ┆ --- │   |
        |  │ i64 ┆ f64 │   |
        |  ╞═════╪═════╡   |
        |  │ 1   ┆ 4.2 │   |
        |  │ 2   ┆ 5.1 │   |
        |  └─────┴─────┘   |
        └──────────────────┘
    """
    if not (supports_arrow_c_stream(native_frame) or is_pyarrow_table(native_frame)):
        msg = f"Given object of type {type(native_frame)} does not support PyCapsule interface"
        raise TypeError(msg)
    implementation = Implementation.from_backend(backend)
    if is_eager_allowed(implementation):
        ns = Version.MAIN.namespace.from_backend(implementation).compliant
        return ns._dataframe.from_arrow(native_frame, context=ns).to_narwhals()
    if implementation is Implementation.UNKNOWN:  # pragma: no cover
        _native_namespace = implementation.to_native_namespace()
        try:
            # implementation is UNKNOWN, Narwhals extension using this feature should
            # implement PyCapsule support
            native: NativeDataFrame = _native_namespace.DataFrame(native_frame)
        except AttributeError as e:
            msg = "Unknown namespace is expected to implement `DataFrame` class which accepts object which supports PyCapsule Interface."
            raise AttributeError(msg) from e
        return from_native(native, eager_only=True)
    msg = (
        f"{implementation} support in Narwhals is lazy-only, but `from_arrow` is an eager-only function.\n\n"
        "Hint: you may want to use an eager backend and then call `.lazy`, e.g.:\n\n"
        f"    nw.from_arrow(df, backend='pyarrow').lazy('{implementation}')"
    )
    raise ValueError(msg)


def _get_sys_info() -> dict[str, str]:
    """System information.

    Returns system and Python version information

    Copied from sklearn

    Returns:
        Dictionary with system info.
    """
    python = sys.version.replace("\n", " ")

    blob = (
        ("python", python),
        ("executable", sys.executable),
        ("machine", platform.platform()),
    )

    return dict(blob)


def _get_deps_info() -> dict[str, str]:
    """Overview of the installed version of main dependencies.

    This function does not import the modules to collect the version numbers
    but instead relies on standard Python package metadata.

    Returns version information on relevant Python libraries

    This function and show_versions were copied from sklearn and adapted

    Returns:
        Mapping from dependency to version.
    """
    from importlib.metadata import distributions

    extra_names = ("narwhals", "numpy")
    member_names = Implementation._member_names_
    exclude = {"PYSPARK_CONNECT", "UNKNOWN"}
    target_names = tuple(
        name.lower() for name in (*extra_names, *member_names) if name not in exclude
    )
    result = dict.fromkeys(target_names, "")  # Initialize with empty strings

    for dist in distributions():
        dist_name, dist_version = dist.name.lower(), dist.version

        if dist_name in result:  # exact match
            result[dist_name] = dist_version
        else:  # prefix match
            for target in target_names:
                if not result[target] and dist_name.startswith(target):
                    result[target] = dist_version
                    break

    return result


def show_versions() -> None:
    """Print useful debugging information.

    Examples:
        >>> from narwhals import show_versions
        >>> show_versions()  # doctest: +SKIP
    """
    sys_info = _get_sys_info()
    deps_info = _get_deps_info()

    print("\nSystem:")  # noqa: T201
    for k, stat in sys_info.items():
        print(f"{k:>10}: {stat}")  # noqa: T201

    print("\nPython dependencies:")  # noqa: T201
    for k, stat in deps_info.items():
        print(f"{k:>13}: {stat}")  # noqa: T201


def read_csv(
    source: FileSource, *, backend: IntoBackend[EagerAllowed], **kwargs: Any
) -> DataFrame[Any]:
    """Read a CSV file into a DataFrame.

    Arguments:
        source: Path to a file.
        backend: The eager backend for DataFrame creation.
            `backend` can be specified in various ways

            - As `Implementation.<BACKEND>` with `BACKEND` being `PANDAS`, `PYARROW`,
                `POLARS`, `MODIN` or `CUDF`.
            - As a string: `"pandas"`, `"pyarrow"`, `"polars"`, `"modin"` or `"cudf"`.
            - Directly as a module `pandas`, `pyarrow`, `polars`, `modin` or `cudf`.
        kwargs: Extra keyword arguments which are passed to the native CSV reader.
            For example, you could use
            `nw.read_csv('file.csv', backend='pandas', engine='pyarrow')`.

    Examples:
        >>> import narwhals as nw
        >>> nw.read_csv("file.csv", backend="pandas")  # doctest:+SKIP
        ┌──────────────────┐
        |Narwhals DataFrame|
        |------------------|
        |        a   b     |
        |     0  1   4     |
        |     1  2   5     |
        └──────────────────┘
    """
    impl = Implementation.from_backend(backend)
    native_namespace = impl.to_native_namespace()
    native_frame: NativeDataFrame
    if impl in {
        Implementation.POLARS,
        Implementation.PANDAS,
        Implementation.MODIN,
        Implementation.CUDF,
    }:
        native_frame = native_namespace.read_csv(normalize_path(source), **kwargs)
    elif impl is Implementation.PYARROW:
        from pyarrow import csv  # ignore-banned-import

        native_frame = csv.read_csv(source, **kwargs)
    elif impl in {
        Implementation.PYSPARK,
        Implementation.DASK,
        Implementation.DUCKDB,
        Implementation.IBIS,
        Implementation.SQLFRAME,
        Implementation.PYSPARK_CONNECT,
    }:
        msg = (
            f"Expected eager backend, found {impl}.\n\n"
            f"Hint: use nw.scan_csv(source={source}, backend={backend})"
        )
        raise ValueError(msg)
    else:  # pragma: no cover
        try:
            # implementation is UNKNOWN, Narwhals extension using this feature should
            # implement `read_csv` function in the top-level namespace.
            native_frame = native_namespace.read_csv(source=source, **kwargs)
        except AttributeError as e:
            msg = "Unknown namespace is expected to implement `read_csv` function."
            raise AttributeError(msg) from e
    return from_native(native_frame, eager_only=True)


def scan_csv(
    source: FileSource, *, backend: IntoBackend[Backend], **kwargs: Any
) -> LazyFrame[Any]:
    """Lazily read from a CSV file.

    For the libraries that do not support lazy dataframes, the function reads
    a csv file eagerly and then converts the resulting dataframe to a lazyframe.

    Arguments:
        source: Path to a file.
        backend: The eager backend for DataFrame creation.
            `backend` can be specified in various ways

            - As `Implementation.<BACKEND>` with `BACKEND` being `PANDAS`, `PYARROW`,
                `POLARS`, `MODIN` or `CUDF`.
            - As a string: `"pandas"`, `"pyarrow"`, `"polars"`, `"modin"` or `"cudf"`.
            - Directly as a module `pandas`, `pyarrow`, `polars`, `modin` or `cudf`.
        kwargs: Extra keyword arguments which are passed to the native CSV reader.
            For example, you could use
            `nw.scan_csv('file.csv', backend=pd, engine='pyarrow')`.

    Examples:
        >>> import duckdb
        >>> import narwhals as nw
        >>>
        >>> nw.scan_csv("file.csv", backend="duckdb").to_native()  # doctest:+SKIP
        ┌─────────┬───────┐
        │    a    │   b   │
        │ varchar │ int32 │
        ├─────────┼───────┤
        │ x       │     1 │
        │ y       │     2 │
        │ z       │     3 │
        └─────────┴───────┘
    """
    implementation = Implementation.from_backend(backend)
    native_namespace = implementation.to_native_namespace()
    native_frame: NativeDataFrame | NativeLazyFrame
    source = normalize_path(source)
    if implementation is Implementation.POLARS:
        native_frame = native_namespace.scan_csv(source, **kwargs)
    elif implementation in {
        Implementation.PANDAS,
        Implementation.MODIN,
        Implementation.CUDF,
        Implementation.DASK,
        Implementation.DUCKDB,
        Implementation.IBIS,
    }:
        native_frame = native_namespace.read_csv(source, **kwargs)
    elif implementation is Implementation.PYARROW:
        from pyarrow import csv  # ignore-banned-import

        native_frame = csv.read_csv(source, **kwargs)
    elif implementation.is_spark_like():
        if (session := kwargs.pop("session", None)) is None:
            msg = "Spark like backends require a session object to be passed in `kwargs`."
            raise ValueError(msg)
        csv_reader = session.read.format("csv")
        native_frame = (
            csv_reader.load(source)
            if (
                implementation is Implementation.SQLFRAME
                and implementation._backend_version() < (3, 27, 0)
            )
            else csv_reader.options(**kwargs).load(source)
        )
    else:  # pragma: no cover
        try:
            # implementation is UNKNOWN, Narwhals extension using this feature should
            # implement `scan_csv` function in the top-level namespace.
            native_frame = native_namespace.scan_csv(source=source, **kwargs)
        except AttributeError as e:
            msg = "Unknown namespace is expected to implement `scan_csv` function."
            raise AttributeError(msg) from e
    return from_native(native_frame).lazy()


def read_parquet(
    source: FileSource, *, backend: IntoBackend[EagerAllowed], **kwargs: Any
) -> DataFrame[Any]:
    """Read into a DataFrame from a parquet file.

    Arguments:
        source: Path to a file.
        backend: The eager backend for DataFrame creation.
            `backend` can be specified in various ways

            - As `Implementation.<BACKEND>` with `BACKEND` being `PANDAS`, `PYARROW`,
                `POLARS`, `MODIN` or `CUDF`.
            - As a string: `"pandas"`, `"pyarrow"`, `"polars"`, `"modin"` or `"cudf"`.
            - Directly as a module `pandas`, `pyarrow`, `polars`, `modin` or `cudf`.
        kwargs: Extra keyword arguments which are passed to the native parquet reader.
            For example, you could use
            `nw.read_parquet('file.parquet', backend=pd, engine='pyarrow')`.

    Examples:
        >>> import pyarrow as pa
        >>> import narwhals as nw
        >>>
        >>> nw.read_parquet("file.parquet", backend="pyarrow")  # doctest:+SKIP
        ┌──────────────────┐
        |Narwhals DataFrame|
        |------------------|
        |pyarrow.Table     |
        |a: int64          |
        |c: double         |
        |----              |
        |a: [[1,2]]        |
        |c: [[0.2,0.1]]    |
        └──────────────────┘
    """
    impl = Implementation.from_backend(backend)
    native_namespace = impl.to_native_namespace()
    native_frame: NativeDataFrame
    if impl in {
        Implementation.POLARS,
        Implementation.PANDAS,
        Implementation.MODIN,
        Implementation.CUDF,
    }:
        source = normalize_path(source)
        native_frame = native_namespace.read_parquet(source, **kwargs)
    elif impl is Implementation.PYARROW:
        import pyarrow.parquet as pq  # ignore-banned-import

        native_frame = pq.read_table(source, **kwargs)  # type: ignore[arg-type]
    elif impl in {
        Implementation.PYSPARK,
        Implementation.DASK,
        Implementation.DUCKDB,
        Implementation.IBIS,
        Implementation.SQLFRAME,
        Implementation.PYSPARK_CONNECT,
    }:
        msg = (
            f"Expected eager backend, found {impl}.\n\n"
            f"Hint: use nw.scan_parquet(source={source}, backend={backend})"
        )
        raise ValueError(msg)
    else:  # pragma: no cover
        try:
            # implementation is UNKNOWN, Narwhals extension using this feature should
            # implement `read_parquet` function in the top-level namespace.
            native_frame = native_namespace.read_parquet(source=source, **kwargs)
        except AttributeError as e:
            msg = "Unknown namespace is expected to implement `read_parquet` function."
            raise AttributeError(msg) from e
    return from_native(native_frame, eager_only=True)


def scan_parquet(
    source: FileSource, *, backend: IntoBackend[Backend], **kwargs: Any
) -> LazyFrame[Any]:
    """Lazily read from a parquet file.

    For the libraries that do not support lazy dataframes, the function reads
    a parquet file eagerly and then converts the resulting dataframe to a lazyframe.

    Note:
        Spark like backends require a session object to be passed in `kwargs`.

        For instance:

        ```py
        import narwhals as nw
        from sqlframe.duckdb import DuckDBSession

        nw.scan_parquet(source, backend="sqlframe", session=DuckDBSession())
        ```

    Arguments:
        source: Path to a file.
        backend: The eager backend for DataFrame creation.
            `backend` can be specified in various ways

            - As `Implementation.<BACKEND>` with `BACKEND` being `PANDAS`, `PYARROW`,
                `POLARS`, `MODIN`, `CUDF`, `PYSPARK` or `SQLFRAME`.
            - As a string: `"pandas"`, `"pyarrow"`, `"polars"`, `"modin"`, `"cudf"`,
                `"pyspark"` or `"sqlframe"`.
            - Directly as a module `pandas`, `pyarrow`, `polars`, `modin`, `cudf`,
                `pyspark.sql` or `sqlframe`.
        kwargs: Extra keyword arguments which are passed to the native parquet reader.
            For example, you could use
            `nw.scan_parquet('file.parquet', backend=pd, engine='pyarrow')`.

    Examples:
        >>> import dask.dataframe as dd
        >>> from sqlframe.duckdb import DuckDBSession
        >>> import narwhals as nw
        >>>
        >>> nw.scan_parquet("file.parquet", backend="dask").collect()  # doctest:+SKIP
        ┌──────────────────┐
        |Narwhals DataFrame|
        |------------------|
        |        a   b     |
        |     0  1   4     |
        |     1  2   5     |
        └──────────────────┘
        >>> nw.scan_parquet(
        ...     "file.parquet", backend="sqlframe", session=DuckDBSession()
        ... ).collect()  # doctest:+SKIP
        ┌──────────────────┐
        |Narwhals DataFrame|
        |------------------|
        |  pyarrow.Table   |
        |  a: int64        |
        |  b: int64        |
        |  ----            |
        |  a: [[1,2]]      |
        |  b: [[4,5]]      |
        └──────────────────┘
    """
    implementation = Implementation.from_backend(backend)
    native_namespace = implementation.to_native_namespace()
    native_frame: NativeDataFrame | NativeLazyFrame
    source = normalize_path(source)
    if implementation is Implementation.POLARS:
        native_frame = native_namespace.scan_parquet(source, **kwargs)
    elif implementation in {
        Implementation.PANDAS,
        Implementation.MODIN,
        Implementation.CUDF,
        Implementation.DASK,
        Implementation.DUCKDB,
        Implementation.IBIS,
    }:
        native_frame = native_namespace.read_parquet(source, **kwargs)
    elif implementation is Implementation.PYARROW:
        import pyarrow.parquet as pq  # ignore-banned-import

        native_frame = pq.read_table(source, **kwargs)
    elif implementation.is_spark_like():
        if (session := kwargs.pop("session", None)) is None:
            msg = "Spark like backends require a session object to be passed in `kwargs`."
            raise ValueError(msg)
        pq_reader = session.read.format("parquet")
        native_frame = (
            pq_reader.load(source)
            if (
                implementation is Implementation.SQLFRAME
                and implementation._backend_version() < (3, 27, 0)
            )
            else pq_reader.options(**kwargs).load(source)
        )

    else:  # pragma: no cover
        try:
            # implementation is UNKNOWN, Narwhals extension using this feature should
            # implement `scan_parquet` function in the top-level namespace.
            native_frame = native_namespace.scan_parquet(source=source, **kwargs)
        except AttributeError as e:
            msg = "Unknown namespace is expected to implement `scan_parquet` function."
            raise AttributeError(msg) from e
    return from_native(native_frame).lazy()


def col(*names: str | Iterable[str]) -> Expr:
    """Creates an expression that references one or more columns by their name(s).

    Arguments:
        names: Name(s) of the columns to use.

    Examples:
        >>> import polars as pl
        >>> import narwhals as nw
        >>>
        >>> df_native = pl.DataFrame({"a": [1, 2], "b": [3, 4], "c": ["x", "z"]})
        >>> nw.from_native(df_native).select(nw.col("a", "b") * nw.col("b"))
        ┌──────────────────┐
        |Narwhals DataFrame|
        |------------------|
        |  shape: (2, 2)   |
        |  ┌─────┬─────┐   |
        |  │ a   ┆ b   │   |
        |  │ --- ┆ --- │   |
        |  │ i64 ┆ i64 │   |
        |  ╞═════╪═════╡   |
        |  │ 3   ┆ 9   │   |
        |  │ 8   ┆ 16  │   |
        |  └─────┴─────┘   |
        └──────────────────┘
    """
    flat_names = flatten(names)

    def func(plx: Any) -> Any:
        return plx.col(*flat_names)

    return Expr(
        func,
        ExprMetadata.selector_single()
        if len(flat_names) == 1
        else ExprMetadata.selector_multi_named(),
    )


def exclude(*names: str | Iterable[str]) -> Expr:
    """Creates an expression that excludes columns by their name(s).

    Arguments:
        names: Name(s) of the columns to exclude.

    Examples:
        >>> import polars as pl
        >>> import narwhals as nw
        >>>
        >>> df_native = pl.DataFrame({"a": [1, 2], "b": [3, 4], "c": ["x", "z"]})
        >>> nw.from_native(df_native).select(nw.exclude("c", "a"))
        ┌──────────────────┐
        |Narwhals DataFrame|
        |------------------|
        |  shape: (2, 1)   |
        |  ┌─────┐         |
        |  │ b   │         |
        |  │ --- │         |
        |  │ i64 │         |
        |  ╞═════╡         |
        |  │ 3   │         |
        |  │ 4   │         |
        |  └─────┘         |
        └──────────────────┘
    """
    exclude_names = frozenset(flatten(names))

    def func(plx: Any) -> Any:
        return plx.exclude(exclude_names)

    return Expr(func, ExprMetadata.selector_multi_unnamed())


def nth(*indices: int | Sequence[int]) -> Expr:
    """Creates an expression that references one or more columns by their index(es).

    Notes:
        `nth` is not supported for Polars version<1.0.0. Please use
        [`narwhals.col`][] instead.

    Arguments:
        indices: One or more indices representing the columns to retrieve.

    Examples:
        >>> import pyarrow as pa
        >>> import narwhals as nw
        >>>
        >>> df_native = pa.table({"a": [1, 2], "b": [3, 4], "c": [0.123, 3.14]})
        >>> nw.from_native(df_native).select(nw.nth(0, 2) * 2)
        ┌──────────────────┐
        |Narwhals DataFrame|
        |------------------|
        |pyarrow.Table     |
        |a: int64          |
        |c: double         |
        |----              |
        |a: [[2,4]]        |
        |c: [[0.246,6.28]] |
        └──────────────────┘
    """
    flat_indices = flatten(indices)

    def func(plx: Any) -> Any:
        return plx.nth(*flat_indices)

    return Expr(
        func,
        ExprMetadata.selector_single()
        if len(flat_indices) == 1
        else ExprMetadata.selector_multi_unnamed(),
    )


# Add underscore so it doesn't conflict with builtin `all`
def all_() -> Expr:
    """Instantiate an expression representing all columns.

    Examples:
        >>> import pandas as pd
        >>> import narwhals as nw
        >>>
        >>> df_native = pd.DataFrame({"a": [1, 2], "b": [3.14, 0.123]})
        >>> nw.from_native(df_native).select(nw.all() * 2)
        ┌──────────────────┐
        |Narwhals DataFrame|
        |------------------|
        |      a      b    |
        |   0  2  6.280    |
        |   1  4  0.246    |
        └──────────────────┘
    """
    return Expr(lambda plx: plx.all(), ExprMetadata.selector_multi_unnamed())


# Add underscore so it doesn't conflict with builtin `len`
def len_() -> Expr:
    """Return the number of rows.

    Examples:
        >>> import polars as pl
        >>> import narwhals as nw
        >>>
        >>> df_native = pl.DataFrame({"a": [1, 2], "b": [5, None]})
        >>> nw.from_native(df_native).select(nw.len())
        ┌──────────────────┐
        |Narwhals DataFrame|
        |------------------|
        |  shape: (1, 1)   |
        |  ┌─────┐         |
        |  │ len │         |
        |  │ --- │         |
        |  │ u32 │         |
        |  ╞═════╡         |
        |  │ 2   │         |
        |  └─────┘         |
        └──────────────────┘
    """

    def func(plx: Any) -> Any:
        return plx.len()

    return Expr(func, ExprMetadata.aggregation())


def sum(*columns: str) -> Expr:
    """Sum all values.

    Note:
        Syntactic sugar for ``nw.col(columns).sum()``

    Arguments:
        columns: Name(s) of the columns to use in the aggregation function

    Examples:
        >>> import pandas as pd
        >>> import narwhals as nw
        >>>
        >>> df_native = pd.DataFrame({"a": [1, 2], "b": [-1.4, 6.2]})
        >>> nw.from_native(df_native).select(nw.sum("a", "b"))
        ┌──────────────────┐
        |Narwhals DataFrame|
        |------------------|
        |       a    b     |
        |    0  3  4.8     |
        └──────────────────┘
    """
    return col(*columns).sum()


def mean(*columns: str) -> Expr:
    """Get the mean value.

    Note:
        Syntactic sugar for ``nw.col(columns).mean()``

    Arguments:
        columns: Name(s) of the columns to use in the aggregation function

    Examples:
        >>> import pyarrow as pa
        >>> import narwhals as nw
        >>>
        >>> df_native = pa.table({"a": [1, 8, 3], "b": [3.14, 6.28, 42.1]})
        >>> nw.from_native(df_native).select(nw.mean("a", "b"))
        ┌─────────────────────────┐
        |   Narwhals DataFrame    |
        |-------------------------|
        |pyarrow.Table            |
        |a: double                |
        |b: double                |
        |----                     |
        |a: [[4]]                 |
        |b: [[17.173333333333336]]|
        └─────────────────────────┘
    """
    return col(*columns).mean()


def median(*columns: str) -> Expr:
    """Get the median value.

    Notes:
        - Syntactic sugar for ``nw.col(columns).median()``
        - Results might slightly differ across backends due to differences in the
            underlying algorithms used to compute the median.

    Arguments:
        columns: Name(s) of the columns to use in the aggregation function

    Examples:
        >>> import polars as pl
        >>> import narwhals as nw
        >>>
        >>> df_native = pl.DataFrame({"a": [4, 5, 2]})
        >>> nw.from_native(df_native).select(nw.median("a"))
        ┌──────────────────┐
        |Narwhals DataFrame|
        |------------------|
        |  shape: (1, 1)   |
        |  ┌─────┐         |
        |  │ a   │         |
        |  │ --- │         |
        |  │ f64 │         |
        |  ╞═════╡         |
        |  │ 4.0 │         |
        |  └─────┘         |
        └──────────────────┘
    """
    return col(*columns).median()


def min(*columns: str) -> Expr:
    """Return the minimum value.

    Note:
       Syntactic sugar for ``nw.col(columns).min()``.

    Arguments:
        columns: Name(s) of the columns to use in the aggregation function.

    Examples:
        >>> import pyarrow as pa
        >>> import narwhals as nw
        >>>
        >>> df_native = pa.table({"a": [1, 2], "b": [5, 10]})
        >>> nw.from_native(df_native).select(nw.min("a", "b"))
        ┌──────────────────┐
        |Narwhals DataFrame|
        |------------------|
        |  pyarrow.Table   |
        |  a: int64        |
        |  b: int64        |
        |  ----            |
        |  a: [[1]]        |
        |  b: [[5]]        |
        └──────────────────┘
    """
    return col(*columns).min()


def max(*columns: str) -> Expr:
    """Return the maximum value.

    Note:
       Syntactic sugar for ``nw.col(columns).max()``.

    Arguments:
        columns: Name(s) of the columns to use in the aggregation function.

    Examples:
        >>> import pandas as pd
        >>> import narwhals as nw
        >>>
        >>> df_native = pd.DataFrame({"a": [1, 2], "b": [5, 10]})
        >>> nw.from_native(df_native).select(nw.max("a", "b"))
        ┌──────────────────┐
        |Narwhals DataFrame|
        |------------------|
        |        a   b     |
        |     0  2  10     |
        └──────────────────┘
    """
    return col(*columns).max()


def _expr_with_n_ary_op(
    func_name: str,
    operation_factory: Callable[
        [CompliantNamespace[Any, Any]], Callable[..., CompliantExpr[Any, Any]]
    ],
    *exprs: IntoExpr,
) -> Expr:
    if not exprs:
        msg = f"At least one expression must be passed to `{func_name}`"
        raise ValueError(msg)
    return Expr(
        lambda plx: apply_n_ary_operation(
            plx, operation_factory(plx), *exprs, str_as_lit=False
        ),
        ExprMetadata.from_horizontal_op(*exprs),
    )


def sum_horizontal(*exprs: IntoExpr | Iterable[IntoExpr]) -> Expr:
    """Sum all values horizontally across columns.

    Warning:
        Unlike Polars, we support horizontal sum over numeric columns only.

    Arguments:
        exprs: Name(s) of the columns to use in the aggregation function. Accepts
            expression input.

    Examples:
        >>> import polars as pl
        >>> import narwhals as nw
        >>>
        >>> df_native = pl.DataFrame({"a": [1, 2, 3], "b": [5, 10, None]})
        >>> nw.from_native(df_native).with_columns(sum=nw.sum_horizontal("a", "b"))
        ┌────────────────────┐
        | Narwhals DataFrame |
        |--------------------|
        |shape: (3, 3)       |
        |┌─────┬──────┬─────┐|
        |│ a   ┆ b    ┆ sum │|
        |│ --- ┆ ---  ┆ --- │|
        |│ i64 ┆ i64  ┆ i64 │|
        |╞═════╪══════╪═════╡|
        |│ 1   ┆ 5    ┆ 6   │|
        |│ 2   ┆ 10   ┆ 12  │|
        |│ 3   ┆ null ┆ 3   │|
        |└─────┴──────┴─────┘|
        └────────────────────┘
    """
    return _expr_with_n_ary_op(
        "sum_horizontal", lambda plx: plx.sum_horizontal, *flatten(exprs)
    )


def min_horizontal(*exprs: IntoExpr | Iterable[IntoExpr]) -> Expr:
    """Get the minimum value horizontally across columns.

    Notes:
        We support `min_horizontal` over numeric columns only.

    Arguments:
        exprs: Name(s) of the columns to use in the aggregation function. Accepts
            expression input.

    Examples:
        >>> import pyarrow as pa
        >>> import narwhals as nw
        >>>
        >>> df_native = pa.table({"a": [1, 8, 3], "b": [4, 5, None]})
        >>> nw.from_native(df_native).with_columns(h_min=nw.min_horizontal("a", "b"))
        ┌──────────────────┐
        |Narwhals DataFrame|
        |------------------|
        | pyarrow.Table    |
        | a: int64         |
        | b: int64         |
        | h_min: int64     |
        | ----             |
        | a: [[1,8,3]]     |
        | b: [[4,5,null]]  |
        | h_min: [[1,5,3]] |
        └──────────────────┘
    """
    return _expr_with_n_ary_op(
        "min_horizontal", lambda plx: plx.min_horizontal, *flatten(exprs)
    )


def max_horizontal(*exprs: IntoExpr | Iterable[IntoExpr]) -> Expr:
    """Get the maximum value horizontally across columns.

    Notes:
        We support `max_horizontal` over numeric columns only.

    Arguments:
        exprs: Name(s) of the columns to use in the aggregation function. Accepts
            expression input.

    Examples:
        >>> import polars as pl
        >>> import narwhals as nw
        >>>
        >>> df_native = pl.DataFrame({"a": [1, 8, 3], "b": [4, 5, None]})
        >>> nw.from_native(df_native).with_columns(h_max=nw.max_horizontal("a", "b"))
        ┌──────────────────────┐
        |  Narwhals DataFrame  |
        |----------------------|
        |shape: (3, 3)         |
        |┌─────┬──────┬───────┐|
        |│ a   ┆ b    ┆ h_max │|
        |│ --- ┆ ---  ┆ ---   │|
        |│ i64 ┆ i64  ┆ i64   │|
        |╞═════╪══════╪═══════╡|
        |│ 1   ┆ 4    ┆ 4     │|
        |│ 8   ┆ 5    ┆ 8     │|
        |│ 3   ┆ null ┆ 3     │|
        |└─────┴──────┴───────┘|
        └──────────────────────┘
    """
    return _expr_with_n_ary_op(
        "max_horizontal", lambda plx: plx.max_horizontal, *flatten(exprs)
    )


class When:
    def __init__(self, *predicates: IntoExpr | Iterable[IntoExpr]) -> None:
        self._predicate = all_horizontal(*flatten(predicates), ignore_nulls=False)

    def then(self, value: IntoExpr | NonNestedLiteral | _1DArray) -> Then:
        kind = ExprKind.from_into_expr(value, str_as_lit=False)
        if self._predicate._metadata.is_scalar_like and not kind.is_scalar_like:
            msg = (
                "If you pass a scalar-like predicate to `nw.when`, then "
                "the `then` value must also be scalar-like."
            )
            raise InvalidOperationError(msg)

        return Then(
            lambda plx: apply_n_ary_operation(
                plx,
                lambda *args: plx.when(args[0]).then(args[1]),
                self._predicate,
                value,
                str_as_lit=False,
            ),
            combine_metadata(
                self._predicate,
                value,
                str_as_lit=False,
                allow_multi_output=False,
                to_single_output=False,
            ),
        )


class Then(Expr):
    def otherwise(self, value: IntoExpr | NonNestedLiteral | _1DArray) -> Expr:
        kind = ExprKind.from_into_expr(value, str_as_lit=False)
        if self._metadata.is_scalar_like and not is_scalar_like(kind):
            msg = (
                "If you pass a scalar-like predicate to `nw.when`, then "
                "the `otherwise` value must also be scalar-like."
            )
            raise InvalidOperationError(msg)

        def func(plx: CompliantNamespace[Any, Any]) -> CompliantExpr[Any, Any]:
            compliant_expr = self._to_compliant_expr(plx)
            compliant_value = plx.parse_into_expr(value, str_as_lit=False)
            if (
                not self._metadata.is_scalar_like
                and is_scalar_like(kind)
                and is_compliant_expr(compliant_value)
            ):
                compliant_value = compliant_value.broadcast(kind)
            return compliant_expr.otherwise(compliant_value)  # type: ignore[attr-defined, no-any-return]

        return Expr(
            func,
            combine_metadata(
                self,
                value,
                str_as_lit=False,
                allow_multi_output=False,
                to_single_output=False,
            ),
        )


def when(*predicates: IntoExpr | Iterable[IntoExpr]) -> When:
    """Start a `when-then-otherwise` expression.

    Expression similar to an `if-else` statement in Python. Always initiated by a
    `pl.when(<condition>).then(<value if condition>)`, and optionally followed by a
    `.otherwise(<value if condition is false>)` can be appended at the end. If not
    appended, and the condition is not `True`, `None` will be returned.

    Info:
        Chaining multiple `.when(<condition>).then(<value>)` statements is currently
        not supported.
        See [Narwhals#668](https://github.com/narwhals-dev/narwhals/issues/668).

    Arguments:
        predicates: Condition(s) that must be met in order to apply the subsequent
            statement. Accepts one or more boolean expressions, which are implicitly
            combined with `&`. String input is parsed as a column name.

    Returns:
        A "when" object, which `.then` can be called on.

    Examples:
        >>> import pandas as pd
        >>> import narwhals as nw
        >>>
        >>> data = {"a": [1, 2, 3], "b": [5, 10, 15]}
        >>> df_native = pd.DataFrame(data)
        >>> nw.from_native(df_native).with_columns(
        ...     nw.when(nw.col("a") < 3).then(5).otherwise(6).alias("a_when")
        ... )
        ┌──────────────────┐
        |Narwhals DataFrame|
        |------------------|
        |    a   b  a_when |
        | 0  1   5       5 |
        | 1  2  10       5 |
        | 2  3  15       6 |
        └──────────────────┘
    """
    return When(*predicates)


def all_horizontal(*exprs: IntoExpr | Iterable[IntoExpr], ignore_nulls: bool) -> Expr:
    r"""Compute the bitwise AND horizontally across columns.

    Arguments:
        exprs: Name(s) of the columns to use in the aggregation function. Accepts
            expression input.
        ignore_nulls: Whether to ignore nulls:

            - If `True`, null values are ignored. If there are no elements, the result
              is `True`.
            - If `False`, Kleene logic is followed. Note that this is not allowed for
              pandas with classical NumPy dtypes when null values are present.

    Examples:
        >>> import pyarrow as pa
        >>> import narwhals as nw
        >>>
        >>> data = {
        ...     "a": [False, False, True, True, False, None],
        ...     "b": [False, True, True, None, None, None],
        ... }
        >>> df_native = pa.table(data)
        >>> nw.from_native(df_native).select(
        ...     "a", "b", all=nw.all_horizontal("a", "b", ignore_nulls=False)
        ... )
        ┌─────────────────────────────────────────┐
        |           Narwhals DataFrame            |
        |-----------------------------------------|
        |pyarrow.Table                            |
        |a: bool                                  |
        |b: bool                                  |
        |all: bool                                |
        |----                                     |
        |a: [[false,false,true,true,false,null]]  |
        |b: [[false,true,true,null,null,null]]    |
        |all: [[false,false,true,null,false,null]]|
        └─────────────────────────────────────────┘
    """
    return _expr_with_n_ary_op(
        "all_horizontal",
        lambda plx: partial(plx.all_horizontal, ignore_nulls=ignore_nulls),
        *flatten(exprs),
    )


def lit(value: NonNestedLiteral, dtype: IntoDType | None = None) -> Expr:
    """Return an expression representing a literal value.

    Arguments:
        value: The value to use as literal.
        dtype: The data type of the literal value. If not provided, the data type will
            be inferred by the native library.

    Examples:
        >>> import pandas as pd
        >>> import narwhals as nw
        >>>
        >>> df_native = pd.DataFrame({"a": [1, 2]})
        >>> nw.from_native(df_native).with_columns(nw.lit(3))
        ┌──────────────────┐
        |Narwhals DataFrame|
        |------------------|
        |     a  literal   |
        |  0  1        3   |
        |  1  2        3   |
        └──────────────────┘
    """
    if is_numpy_array(value):
        msg = (
            "numpy arrays are not supported as literal values. "
            "Consider using `with_columns` to create a new column from the array."
        )
        raise ValueError(msg)

    if isinstance(value, (list, tuple)):
        msg = f"Nested datatypes are not supported yet. Got {value}"
        raise NotImplementedError(msg)

    return Expr(lambda plx: plx.lit(value, dtype), ExprMetadata.literal())


def any_horizontal(*exprs: IntoExpr | Iterable[IntoExpr], ignore_nulls: bool) -> Expr:
    r"""Compute the bitwise OR horizontally across columns.

    Arguments:
        exprs: Name(s) of the columns to use in the aggregation function. Accepts
            expression input.
        ignore_nulls: Whether to ignore nulls:

            - If `True`, null values are ignored. If there are no elements, the result
              is `False`.
            - If `False`, Kleene logic is followed. Note that this is not allowed for
              pandas with classical NumPy dtypes when null values are present.

    Examples:
        >>> import polars as pl
        >>> import narwhals as nw
        >>>
        >>> data = {
        ...     "a": [False, False, True, True, False, None],
        ...     "b": [False, True, True, None, None, None],
        ... }
        >>> df_native = pl.DataFrame(data)
        >>> nw.from_native(df_native).select(
        ...     "a", "b", any=nw.any_horizontal("a", "b", ignore_nulls=False)
        ... )
        ┌─────────────────────────┐
        |   Narwhals DataFrame    |
        |-------------------------|
        |shape: (6, 3)            |
        |┌───────┬───────┬───────┐|
        |│ a     ┆ b     ┆ any   │|
        |│ ---   ┆ ---   ┆ ---   │|
        |│ bool  ┆ bool  ┆ bool  │|
        |╞═══════╪═══════╪═══════╡|
        |│ false ┆ false ┆ false │|
        |│ false ┆ true  ┆ true  │|
        |│ true  ┆ true  ┆ true  │|
        |│ true  ┆ null  ┆ true  │|
        |│ false ┆ null  ┆ null  │|
        |│ null  ┆ null  ┆ null  │|
        |└───────┴───────┴───────┘|
        └─────────────────────────┘
    """
    return _expr_with_n_ary_op(
        "any_horizontal",
        lambda plx: partial(plx.any_horizontal, ignore_nulls=ignore_nulls),
        *flatten(exprs),
    )


def mean_horizontal(*exprs: IntoExpr | Iterable[IntoExpr]) -> Expr:
    """Compute the mean of all values horizontally across columns.

    Arguments:
        exprs: Name(s) of the columns to use in the aggregation function. Accepts
            expression input.

    Examples:
        >>> import pyarrow as pa
        >>> import narwhals as nw
        >>>
        >>> data = {"a": [1, 8, 3], "b": [4, 5, None], "c": ["x", "y", "z"]}
        >>> df_native = pa.table(data)

        We define a dataframe-agnostic function that computes the horizontal mean of "a"
        and "b" columns:

        >>> nw.from_native(df_native).select(nw.mean_horizontal("a", "b"))
        ┌──────────────────┐
        |Narwhals DataFrame|
        |------------------|
        | pyarrow.Table    |
        | a: double        |
        | ----             |
        | a: [[2.5,6.5,3]] |
        └──────────────────┘
    """
    return _expr_with_n_ary_op(
        "mean_horizontal", lambda plx: plx.mean_horizontal, *flatten(exprs)
    )


def concat_str(
    exprs: IntoExpr | Iterable[IntoExpr],
    *more_exprs: IntoExpr,
    separator: str = "",
    ignore_nulls: bool = False,
) -> Expr:
    r"""Horizontally concatenate columns into a single string column.

    Arguments:
        exprs: Columns to concatenate into a single string column. Accepts expression
            input. Strings are parsed as column names, other non-expression inputs are
            parsed as literals. Non-`String` columns are cast to `String`.
        *more_exprs: Additional columns to concatenate into a single string column,
            specified as positional arguments.
        separator: String that will be used to separate the values of each column.
        ignore_nulls: Ignore null values (default is `False`).
            If set to `False`, null values will be propagated and if the row contains any
            null values, the output is null.

    Examples:
        >>> import pandas as pd
        >>> import narwhals as nw
        >>>
        >>> data = {
        ...     "a": [1, 2, 3],
        ...     "b": ["dogs", "cats", None],
        ...     "c": ["play", "swim", "walk"],
        ... }
        >>> df_native = pd.DataFrame(data)
        >>> (
        ...     nw.from_native(df_native).select(
        ...         nw.concat_str(
        ...             [nw.col("a") * 2, nw.col("b"), nw.col("c")], separator=" "
        ...         ).alias("full_sentence")
        ...     )
        ... )
        ┌──────────────────┐
        |Narwhals DataFrame|
        |------------------|
        |   full_sentence  |
        | 0   2 dogs play  |
        | 1   4 cats swim  |
        | 2          None  |
        └──────────────────┘
    """
    flat_exprs = flatten([*flatten([exprs]), *more_exprs])
    return _expr_with_n_ary_op(
        "concat_str",
        lambda plx: lambda *args: plx.concat_str(
            *args, separator=separator, ignore_nulls=ignore_nulls
        ),
        *flat_exprs,
    )


def coalesce(
    exprs: IntoExpr | Iterable[IntoExpr], *more_exprs: IntoExpr | NonNestedLiteral
) -> Expr:
    """Folds the columns from left to right, keeping the first non-null value.

    Arguments:
        exprs: Columns to coalesce, must be a str, nw.Expr, or nw.Series
            where strings are parsed as column names and both nw.Expr/nw.Series
            are passed through as-is. Scalar values must be wrapped in `nw.lit`.

        *more_exprs: Additional columns to coalesce, specified as positional arguments.

    Raises:
        TypeError: If any of the inputs are not a str, nw.Expr, or nw.Series.

    Examples:
        >>> import polars as pl
        >>> import narwhals as nw
        >>> data = [
        ...     (1, 5, None),
        ...     (None, 6, None),
        ...     (None, None, 9),
        ...     (4, 8, 10),
        ...     (None, None, None),
        ... ]
        >>> df = pl.DataFrame(data, schema=["a", "b", "c"], orient="row")
        >>> nw.from_native(df).select(nw.coalesce("a", "b", "c", nw.lit(-1)))
        ┌──────────────────┐
        |Narwhals DataFrame|
        |------------------|
        |  shape: (5, 1)   |
        |  ┌─────┐         |
        |  │ a   │         |
        |  │ --- │         |
        |  │ i64 │         |
        |  ╞═════╡         |
        |  │ 1   │         |
        |  │ 6   │         |
        |  │ 9   │         |
        |  │ 4   │         |
        |  │ -1  │         |
        |  └─────┘         |
        └──────────────────┘
    """
    flat_exprs = flatten([*flatten([exprs]), *more_exprs])

    non_exprs = [expr for expr in flat_exprs if not isinstance(expr, (str, Expr, Series))]
    if non_exprs:
        msg = (
            f"All arguments to `coalesce` must be of type {str!r}, {Expr!r}, or {Series!r}."
            "\nGot the following invalid arguments (type, value):"
            f"\n    {', '.join(repr((type(e), e)) for e in non_exprs)}"
        )
        raise TypeError(msg)

    return Expr(
        lambda plx: apply_n_ary_operation(
            plx, lambda *args: plx.coalesce(*args), *flat_exprs, str_as_lit=False
        ),
        ExprMetadata.from_horizontal_op(*flat_exprs),
    )


def format(f_string: str, *args: IntoExpr) -> Expr:
    """Format expressions as a string.

    Arguments:
        f_string: A string that with placeholders.
        args: Expression(s) that fill the placeholders.

    Examples:
        >>> import duckdb
        >>> import narwhals as nw
        >>> rel = duckdb.sql("select * from values ('a', 1), ('b', 2), ('c', 3) df(a, b)")
        >>> df = nw.from_native(rel)
        >>> df.with_columns(formatted=nw.format("foo_{}_bar_{}", nw.col("a"), "b"))
        ┌─────────────────────────────────┐
        |       Narwhals LazyFrame        |
        |---------------------------------|
        |┌─────────┬───────┬─────────────┐|
        |│    a    │   b   │  formatted  │|
        |│ varchar │ int32 │   varchar   │|
        |├─────────┼───────┼─────────────┤|
        |│ a       │     1 │ foo_a_bar_1 │|
        |│ b       │     2 │ foo_b_bar_2 │|
        |│ c       │     3 │ foo_c_bar_3 │|
        |└─────────┴───────┴─────────────┘|
        └─────────────────────────────────┘
    """
    if (n_placeholders := f_string.count("{}")) != len(args):
        msg = f"number of placeholders should equal the number of arguments. Expected {n_placeholders} arguments, got {len(args)}."
        raise ValueError(msg)

    exprs = []
    it = iter(args)
    for i, s in enumerate(f_string.split("{}")):
        if i > 0:
            exprs.append(next(it))
        if len(s) > 0:
            exprs.append(lit(s))
    return concat_str(exprs, separator="")
