from __future__ import annotations

from typing import TYPE_CHECKING, Any

from narwhals._pandas_like.utils import select_columns_by_name
from narwhals._utils import Implementation, Version, isinstance_or_issubclass
from narwhals.dependencies import get_pyarrow

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

    import dask.dataframe as dd
    import dask.dataframe.dask_expr as dx

    from narwhals._dask.dataframe import DaskLazyFrame, Incomplete
    from narwhals._dask.expr import DaskExpr
    from narwhals.dtypes import DType
    from narwhals.typing import IntoDType
else:
    try:
        import dask.dataframe.dask_expr as dx
    except ModuleNotFoundError:
        import dask_expr as dx


def maybe_evaluate_expr(df: DaskLazyFrame, obj: DaskExpr | object) -> dx.Series | object:
    from narwhals._dask.expr import DaskExpr

    if isinstance(obj, DaskExpr):
        results = obj._call(df)
        assert len(results) == 1  # debug assertion  # noqa: S101
        return results[0]
    return obj


def evaluate_exprs(df: DaskLazyFrame, /, *exprs: DaskExpr) -> list[tuple[str, dx.Series]]:
    native_results: list[tuple[str, dx.Series]] = []
    for expr in exprs:
        native_series_list = expr(df)
        aliases = expr._evaluate_aliases(df)
        if len(aliases) != len(native_series_list):  # pragma: no cover
            msg = f"Internal error: got aliases {aliases}, but only got {len(native_series_list)} results"
            raise AssertionError(msg)
        native_results.extend(zip(aliases, native_series_list))
    return native_results


def align_series_full_broadcast(
    df: DaskLazyFrame, *series: dx.Series | object
) -> Sequence[dx.Series]:
    return [
        s if isinstance(s, dx.Series) else df._native_frame.assign(_tmp=s)["_tmp"]
        for s in series
    ]  # pyright: ignore[reportReturnType]


def add_row_index(frame: dd.DataFrame, name: str) -> dd.DataFrame:
    original_cols = frame.columns
    df: Incomplete = frame.assign(**{name: 1})
    return select_columns_by_name(
        df.assign(**{name: df[name].cumsum(method="blelloch") - 1}),
        [name, *original_cols],
        Implementation.DASK,
    )


def validate_comparand(lhs: dx.Series, rhs: dx.Series) -> None:
    if not dx.expr.are_co_aligned(lhs._expr, rhs._expr):  # pragma: no cover
        # are_co_aligned is a method which cheaply checks if two Dask expressions
        # have the same index, and therefore don't require index alignment.
        # If someone only operates on a Dask DataFrame via expressions, then this
        # should always be the case: expression outputs (by definition) all come from the
        # same input dataframe, and Dask Series does not have any operations which
        # change the index. Nonetheless, we perform this safety check anyway.

        # However, we still need to carefully vet which methods we support for Dask, to
        # avoid issues where `are_co_aligned` doesn't do what we want it to do:
        # https://github.com/dask/dask-expr/issues/1112.
        msg = "Objects are not co-aligned, so this operation is not supported for Dask backend"
        raise RuntimeError(msg)


dtypes = Version.MAIN.dtypes
dtypes_v1 = Version.V1.dtypes
NW_TO_DASK_DTYPES: Mapping[type[DType], str] = {
    dtypes.Float64: "float64",
    dtypes.Float32: "float32",
    dtypes.Boolean: "bool",
    dtypes.Categorical: "category",
    dtypes.Date: "date32[day][pyarrow]",
    dtypes.Int8: "int8",
    dtypes.Int16: "int16",
    dtypes.Int32: "int32",
    dtypes.Int64: "int64",
    dtypes.UInt8: "uint8",
    dtypes.UInt16: "uint16",
    dtypes.UInt32: "uint32",
    dtypes.UInt64: "uint64",
    dtypes.Datetime: "datetime64[us]",
    dtypes.Duration: "timedelta64[ns]",
    dtypes_v1.Datetime: "datetime64[us]",
    dtypes_v1.Duration: "timedelta64[ns]",
}
UNSUPPORTED_DTYPES = (
    dtypes.List,
    dtypes.Struct,
    dtypes.Array,
    dtypes.Time,
    dtypes.Binary,
)


def narwhals_to_native_dtype(dtype: IntoDType, version: Version) -> Any:
    dtypes = version.dtypes
    base_type = dtype.base_type()
    if dask_type := NW_TO_DASK_DTYPES.get(base_type):
        return dask_type
    if isinstance_or_issubclass(dtype, dtypes.String):
        if Implementation.PANDAS._backend_version() >= (2, 0, 0):
            return "string[pyarrow]" if get_pyarrow() else "string[python]"
        return "object"  # pragma: no cover
    if isinstance_or_issubclass(dtype, dtypes.Enum):
        if version is Version.V1:
            msg = "Converting to Enum is not supported in narwhals.stable.v1"
            raise NotImplementedError(msg)
        if isinstance(dtype, dtypes.Enum):
            import pandas as pd

            # NOTE: `pandas-stubs.core.dtypes.dtypes.CategoricalDtype.categories` is too narrow
            # Should be one of the `ListLike*` types
            # https://github.com/pandas-dev/pandas-stubs/blob/8434bde95460b996323cc8c0fea7b0a8bb00ea26/pandas-stubs/_typing.pyi#L497-L505
            return pd.CategoricalDtype(dtype.categories, ordered=True)  # type: ignore[arg-type]
        msg = "Can not cast / initialize Enum without categories present"
        raise ValueError(msg)
    if issubclass(base_type, UNSUPPORTED_DTYPES):  # pragma: no cover
        msg = f"Converting to {base_type.__name__} dtype is not supported for Dask."
        raise NotImplementedError(msg)
    msg = f"Unknown dtype: {dtype}"  # pragma: no cover
    raise AssertionError(msg)
