from __future__ import annotations

import os
import sys
from collections.abc import Iterator
from typing import TYPE_CHECKING, Callable

import polars._reexport as pl
from polars._utils.unstable import unstable

if TYPE_CHECKING:
    from collections.abc import Iterator
    from typing import Callable

    from polars import DataFrame, Expr, LazyFrame
    from polars._typing import SchemaDict


@unstable()
def register_io_source(
    callable: Callable[
        [list[str] | None, Expr | None, int | None, int | None], Iterator[DataFrame]
    ],
    schema: SchemaDict,
) -> LazyFrame:
    """
    Register your IO plugin and initialize a LazyFrame.

    Parameters
    ----------
    callable
        Function that accepts the following arguments:
            with_columns
                Columns that are projected. The reader must
                project these columns if applied
            predicate
                Polars expression. The reader must filter
                their rows accordingly.
            n_rows:
                Materialize only n rows from the source.
                The reader can stop when `n_rows` are read.
            batch_size
                A hint of the ideal batch size the reader's
                generator must produce.
        The function should return a DataFrame batch
        (an iterator over individual DataFrames).
    schema
        Schema that the reader will produce before projection pushdown.

    """

    def wrap(
        with_columns: list[str] | None,
        predicate: bytes | None,
        n_rows: int | None,
        batch_size: int | None,
    ) -> tuple[Iterator[DataFrame], bool]:
        parsed_predicate_success = True
        parsed_predicate = None
        if predicate:
            try:
                parsed_predicate = pl.Expr.deserialize(predicate)
            except Exception as e:
                if os.environ.get("POLARS_VERBOSE"):
                    print(
                        f"failed parsing IO plugin expression\n\nfilter will be handled on Polars' side: {e}",
                        file=sys.stderr,
                    )
                parsed_predicate_success = False

        return callable(
            with_columns, parsed_predicate, n_rows, batch_size
        ), parsed_predicate_success

    return pl.LazyFrame._scan_python_function(
        schema=schema, scan_fn=wrap, pyarrow=False
    )
