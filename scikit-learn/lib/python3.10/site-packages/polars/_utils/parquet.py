from typing import Any, Callable

from polars._typing import ParquetMetadataContext, ParquetMetadataFn


def wrap_parquet_metadata_callback(
    fn: ParquetMetadataFn,
) -> Callable[[Any], list[tuple[str, str]]]:
    def pyo3_compatible_callback(ctx: Any) -> list[tuple[str, str]]:
        ctx_py = ParquetMetadataContext(
            arrow_schema=ctx.arrow_schema,
        )
        return list(fn(ctx_py).items())

    return pyo3_compatible_callback
