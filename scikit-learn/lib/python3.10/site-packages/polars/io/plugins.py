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
    io_source: Callable[
        [list[str] | None, Expr | None, int | None, int | None], Iterator[DataFrame]
    ],
    *,
    schema: Callable[[], SchemaDict] | SchemaDict,
    validate_schema: bool = False,
    is_pure: bool = False,
) -> LazyFrame:
    """
    Register your IO plugin and initialize a LazyFrame.

    See the `user guide <https://docs.pola.rs/user-guide/plugins/io_plugins>`_
    for more information about plugins.

    .. warning::
        This functionality is considered **unstable**. It may be changed
        at any point without it being considered a breaking change.


    Parameters
    ----------
    io_source
        Function that accepts the following arguments:
            with_columns
                Columns that are projected. The reader must
                project these columns if applied
            predicate
                Polars expression. The reader must filter
                their rows accordingly.
            n_rows
                Materialize only n rows from the source.
                The reader can stop when `n_rows` are read.
            batch_size
                A hint of the ideal batch size the reader's
                generator must produce.

        The function should return a an iterator/generator
        that produces DataFrames.
    schema
        Schema or function that when called produces the schema that the reader
        will produce before projection pushdown.
    validate_schema
        Whether the engine should validate if the batches generated match
        the given schema. It's an implementation error if this isn't
        the case and can lead to bugs that are hard to solve.
    is_pure
        Whether the IO source is pure. Repeated occurrences of same IO source in
        a LazyFrame plan can be de-duplicated during optimization if they are
        pure.

    Returns
    -------
    LazyFrame
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

        return io_source(
            with_columns, parsed_predicate, n_rows, batch_size
        ), parsed_predicate_success

    return pl.LazyFrame._scan_python_function(
        schema=schema,
        scan_fn=wrap,
        pyarrow=False,
        validate_schema=validate_schema,
        is_pure=is_pure,
    )


@unstable()
def _defer(
    function: Callable[[], DataFrame],
    *,
    schema: SchemaDict | Callable[[], SchemaDict],
    validate_schema: bool = True,
) -> LazyFrame:
    """
    Deferred execution.

    Takes a function that produces a `DataFrame` but defers execution until the
    `LazyFrame` is collected.

    Parameters
    ----------
    function
        Function that takes no arguments and produces a `DataFrame`.
    schema
        Schema of the `DataFrame` the deferred function will return.
        The caller must ensure this schema is correct.
    validate_schema
        Whether the engine should validate if the batches generated match
        the given schema. It's an implementation error if this isn't
        the case and can lead to bugs that are hard to solve.

    Examples
    --------
    Delay DataFrame execution until query is executed.

    >>> import numpy as np
    >>> np.random.seed(0)
    >>> lf = pl.defer(
    ...     lambda: pl.DataFrame({"a": np.random.randn(3)}), schema={"a": pl.Float64}
    ... )
    >>> lf.collect()
    shape: (3, 1)
    ┌──────────┐
    │ a        │
    │ ---      │
    │ f64      │
    ╞══════════╡
    │ 1.764052 │
    │ 0.400157 │
    │ 0.978738 │
    └──────────┘

     Run an eager source in Polars Cloud

    >>> (
    ...     pl.defer(
    ...         lambda: pl.read_database("select * from tbl"),
    ...         schema={"a": pl.Float64, "b": pl.Boolean},
    ...     )
    ...     .filter("b")
    ...     .sum("a")
    ...     .remote()
    ...     .collect()
    ... )  # doctest: +SKIP


    """

    def source(
        with_columns: list[str] | None,
        predicate: Expr | None,
        n_rows: int | None,
        batch_size: int | None,
    ) -> Iterator[DataFrame]:
        lf = function().lazy()
        if with_columns is not None:
            lf = lf.select(with_columns)
        if predicate is not None:
            lf = lf.filter(predicate)
        if n_rows is not None:
            lf = lf.limit(n_rows)
        yield lf.collect()

    return register_io_source(
        io_source=source, schema=schema, validate_schema=validate_schema
    )
