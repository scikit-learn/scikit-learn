from __future__ import annotations

import ast
from _ast import GtE, Lt, LtE
from ast import (
    Attribute,
    BinOp,
    BitAnd,
    BitOr,
    Call,
    Compare,
    Constant,
    Eq,
    Gt,
    Invert,
    List,
    Name,
    UnaryOp,
)
from functools import partial, singledispatch
from typing import TYPE_CHECKING, Any, Callable

import polars._reexport as pl
from polars._utils.convert import to_py_date, to_py_datetime
from polars.dependencies import pyiceberg

if TYPE_CHECKING:
    from datetime import date, datetime

    from pyiceberg.table import Table

    from polars import DataFrame, LazyFrame, Series

__all__ = ["scan_iceberg"]

_temporal_conversions: dict[str, Callable[..., datetime | date]] = {
    "to_py_date": to_py_date,
    "to_py_datetime": to_py_datetime,
}


def scan_iceberg(
    source: str | Table,
    *,
    snapshot_id: int | None = None,
    storage_options: dict[str, Any] | None = None,
) -> LazyFrame:
    """
    Lazily read from an Apache Iceberg table.

    Parameters
    ----------
    source
        A PyIceberg table, or a direct path to the metadata.

        Note: For Local filesystem, absolute and relative paths are supported but
        for the supported object storages - GCS, Azure and S3 full URI must be provided.
    snapshot_id
        The snapshot ID to scan from.
    storage_options
        Extra options for the storage backends supported by `pyiceberg`.
        For cloud storages, this may include configurations for authentication etc.

        More info is available `here <https://py.iceberg.apache.org/configuration/>`__.

    Returns
    -------
    LazyFrame

    Examples
    --------
    Creates a scan for an Iceberg table from local filesystem, or object store.

    >>> table_path = "file:/path/to/iceberg-table/metadata.json"
    >>> pl.scan_iceberg(table_path).collect()  # doctest: +SKIP

    Creates a scan for an Iceberg table from S3.
    See a list of supported storage options for S3 `here
    <https://py.iceberg.apache.org/configuration/#fileio>`__.

    >>> table_path = "s3://bucket/path/to/iceberg-table/metadata.json"
    >>> storage_options = {
    ...     "s3.region": "eu-central-1",
    ...     "s3.access-key-id": "THE_AWS_ACCESS_KEY_ID",
    ...     "s3.secret-access-key": "THE_AWS_SECRET_ACCESS_KEY",
    ... }
    >>> pl.scan_iceberg(
    ...     table_path, storage_options=storage_options
    ... ).collect()  # doctest: +SKIP

    Creates a scan for an Iceberg table from Azure.
    Supported options for Azure are available `here
    <https://py.iceberg.apache.org/configuration/#azure-data-lake>`__.

    Following type of table paths are supported:

    * az://<container>/<path>/metadata.json
    * adl://<container>/<path>/metadata.json
    * abfs[s]://<container>/<path>/metadata.json

    >>> table_path = "az://container/path/to/iceberg-table/metadata.json"
    >>> storage_options = {
    ...     "adlfs.account-name": "AZURE_STORAGE_ACCOUNT_NAME",
    ...     "adlfs.account-key": "AZURE_STORAGE_ACCOUNT_KEY",
    ... }
    >>> pl.scan_iceberg(
    ...     table_path, storage_options=storage_options
    ... ).collect()  # doctest: +SKIP

    Creates a scan for an Iceberg table from Google Cloud Storage.
    Supported options for GCS are available `here
    <https://py.iceberg.apache.org/configuration/#google-cloud-storage>`__.

    >>> table_path = "s3://bucket/path/to/iceberg-table/metadata.json"
    >>> storage_options = {
    ...     "gcs.project-id": "my-gcp-project",
    ...     "gcs.oauth.token": "ya29.dr.AfM...",
    ... }
    >>> pl.scan_iceberg(
    ...     table_path, storage_options=storage_options
    ... ).collect()  # doctest: +SKIP

    Creates a scan for an Iceberg table with additional options.
    In the below example, `without_files` option is used which loads the table without
    file tracking information.

    >>> table_path = "/path/to/iceberg-table/metadata.json"
    >>> storage_options = {"py-io-impl": "pyiceberg.io.fsspec.FsspecFileIO"}
    >>> pl.scan_iceberg(
    ...     table_path, storage_options=storage_options
    ... ).collect()  # doctest: +SKIP

    Creates a scan for an Iceberg table using a specific snapshot ID.

    >>> table_path = "/path/to/iceberg-table/metadata.json"
    >>> snapshot_id = 7051579356916758811
    >>> pl.scan_iceberg(table_path, snapshot_id=snapshot_id).collect()  # doctest: +SKIP
    """
    from pyiceberg.io.pyarrow import schema_to_pyarrow
    from pyiceberg.table import StaticTable

    if isinstance(source, str):
        source = StaticTable.from_metadata(
            metadata_location=source, properties=storage_options or {}
        )

    if snapshot_id is not None:
        snapshot = source.snapshot_by_id(snapshot_id)
        if snapshot is None:
            msg = f"Snapshot ID not found: {snapshot_id}"
            raise ValueError(msg)

    func = partial(_scan_pyarrow_dataset_impl, source, snapshot_id=snapshot_id)
    arrow_schema = schema_to_pyarrow(source.schema())
    return pl.LazyFrame._scan_python_function(arrow_schema, func, pyarrow=True)


def _scan_pyarrow_dataset_impl(
    tbl: Table,
    with_columns: list[str] | None = None,
    predicate: str = "",
    n_rows: int | None = None,
    snapshot_id: int | None = None,
    **kwargs: Any,
) -> DataFrame | Series:
    """
    Take the projected columns and materialize an arrow table.

    Parameters
    ----------
    tbl
        pyarrow dataset
    with_columns
        Columns that are projected
    predicate
        pyarrow expression that can be evaluated with eval
    n_rows:
        Materialize only n rows from the arrow dataset.
    snapshot_id:
        The snapshot ID to scan from.
    batch_size
        The maximum row count for scanned pyarrow record batches.
    kwargs:
        For backward compatibility

    Returns
    -------
    DataFrame
    """
    from polars import from_arrow

    scan = tbl.scan(limit=n_rows, snapshot_id=snapshot_id)

    if with_columns is not None:
        scan = scan.select(*with_columns)

    if predicate is not None:
        try:
            expr_ast = _to_ast(predicate)
            pyiceberg_expr = _convert_predicate(expr_ast)
        except ValueError as e:
            msg = f"Could not convert predicate to PyIceberg: {predicate}"
            raise ValueError(msg) from e

        scan = scan.filter(pyiceberg_expr)

    return from_arrow(scan.to_arrow())


def _to_ast(expr: str) -> ast.expr:
    """
    Converts a Python string to an AST.

    This will take the Python Arrow expression (as a string), and it will
    be converted into a Python AST that can be traversed to convert it to a PyIceberg
    expression.

    The reason to convert it to an AST is because the PyArrow expression
    itself doesn't have any methods/properties to traverse the expression.
    We need this to convert it into a PyIceberg expression.

    Parameters
    ----------
    expr
        The string expression

    Returns
    -------
    The AST representing the Arrow expression
    """
    return ast.parse(expr, mode="eval").body


@singledispatch
def _convert_predicate(a: Any) -> Any:
    """Walks the AST to convert the PyArrow expression to a PyIceberg expression."""
    msg = f"Unexpected symbol: {a}"
    raise ValueError(msg)


@_convert_predicate.register(Constant)
def _(a: Constant) -> Any:
    return a.value


@_convert_predicate.register(Name)
def _(a: Name) -> Any:
    return a.id


@_convert_predicate.register(UnaryOp)
def _(a: UnaryOp) -> Any:
    if isinstance(a.op, Invert):
        return pyiceberg.expressions.Not(_convert_predicate(a.operand))
    else:
        msg = f"Unexpected UnaryOp: {a}"
        raise TypeError(msg)


@_convert_predicate.register(Call)
def _(a: Call) -> Any:
    args = [_convert_predicate(arg) for arg in a.args]
    f = _convert_predicate(a.func)
    if f == "field":
        return args
    elif f in _temporal_conversions:
        # convert from polars-native i64 to ISO8601 string
        return _temporal_conversions[f](*args).isoformat()
    else:
        ref = _convert_predicate(a.func.value)[0]  # type: ignore[attr-defined]
        if f == "isin":
            return pyiceberg.expressions.In(ref, args[0])
        elif f == "is_null":
            return pyiceberg.expressions.IsNull(ref)
        elif f == "is_nan":
            return pyiceberg.expressions.IsNaN(ref)

    msg = f"Unknown call: {f!r}"
    raise ValueError(msg)


@_convert_predicate.register(Attribute)
def _(a: Attribute) -> Any:
    return a.attr


@_convert_predicate.register(BinOp)
def _(a: BinOp) -> Any:
    lhs = _convert_predicate(a.left)
    rhs = _convert_predicate(a.right)

    op = a.op
    if isinstance(op, BitAnd):
        return pyiceberg.expressions.And(lhs, rhs)
    if isinstance(op, BitOr):
        return pyiceberg.expressions.Or(lhs, rhs)
    else:
        msg = f"Unknown: {lhs} {op} {rhs}"
        raise TypeError(msg)


@_convert_predicate.register(Compare)
def _(a: Compare) -> Any:
    op = a.ops[0]
    lhs = _convert_predicate(a.left)[0]
    rhs = _convert_predicate(a.comparators[0])

    if isinstance(op, Gt):
        return pyiceberg.expressions.GreaterThan(lhs, rhs)
    if isinstance(op, GtE):
        return pyiceberg.expressions.GreaterThanOrEqual(lhs, rhs)
    if isinstance(op, Eq):
        return pyiceberg.expressions.EqualTo(lhs, rhs)
    if isinstance(op, Lt):
        return pyiceberg.expressions.LessThan(lhs, rhs)
    if isinstance(op, LtE):
        return pyiceberg.expressions.LessThanOrEqual(lhs, rhs)
    else:
        msg = f"Unknown comparison: {op}"
        raise TypeError(msg)


@_convert_predicate.register(List)
def _(a: List) -> Any:
    return [_convert_predicate(e) for e in a.elts]
