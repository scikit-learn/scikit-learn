from __future__ import annotations

import contextlib
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING

from polars import DataFrame, col
from polars._typing import PartitioningScheme
from polars._utils.unstable import issue_unstable_warning
from polars.expr import Expr

if TYPE_CHECKING:
    with contextlib.suppress(ImportError):  # Module not available when building docs
        from polars._plr import PyDataFrame, PyExpr

    from typing import IO, Any, Callable

with contextlib.suppress(ImportError):  # Module not available when building docs
    from polars._plr import PyPartitioning


class KeyedPartition:
    """
    A key-value pair for a partition.

    .. warning::
        This functionality is currently considered **unstable**. It may be
        changed at any point without it being considered a breaking change.

    See Also
    --------
    PartitionByKey
    PartitionParted
    KeyedPartitionContext
    """

    def __init__(self, name: str, str_value: str, raw_value: Any) -> None:
        self.name = name
        self.str_value = str_value
        self.raw_value = raw_value

    name: str  #: Name of the key column.
    str_value: str  #: Value of the key as a path and URL safe string.
    raw_value: Any  #: Value of the key for this partition.

    def hive_name(self) -> str:
        """Get the `key=value`."""
        return f"{self.name}={self.str_value}"


class KeyedPartitionContext:
    """
    Callback context for a partition creation using keys.

    .. warning::
        This functionality is currently considered **unstable**. It may be
        changed at any point without it being considered a breaking change.

    See Also
    --------
    PartitionByKey
    PartitionParted
    """

    def __init__(
        self,
        file_idx: int,
        part_idx: int,
        in_part_idx: int,
        keys: list[KeyedPartition],
        file_path: Path,
        full_path: Path,
    ) -> None:
        self.file_idx = file_idx
        self.part_idx = part_idx
        self.in_part_idx = in_part_idx
        self.keys = keys
        self.file_path = file_path
        self.full_path = full_path

    file_idx: int  #: The index of the created file starting from zero.
    part_idx: int  #: The index of the created partition starting from zero.
    in_part_idx: int  #: The index of the file within this partition starting from zero.
    keys: list[KeyedPartition]  #: All the key names and values used for this partition.
    file_path: Path  #: The chosen output path before the callback was called without `base_path`.
    full_path: (
        Path  #: The chosen output path before the callback was called with `base_path`.
    )

    def hive_dirs(self) -> Path:
        """The keys mapped to hive directories."""
        assert len(self.keys) > 0
        p = Path(self.keys[0].hive_name())
        for key in self.keys[1:]:
            p /= Path(key.hive_name())
        return p


class BasePartitionContext:
    """
    Callback context for a partition creation.

    .. warning::
        This functionality is currently considered **unstable**. It may be
        changed at any point without it being considered a breaking change.

    See Also
    --------
    PartitionMaxSize
    """

    def __init__(self, file_idx: int, file_path: Path, full_path: Path) -> None:
        self.file_idx = file_idx
        self.file_path = file_path
        self.full_path = full_path

    file_idx: int  #: The index of the created file starting from zero.
    file_path: Path  #: The chosen output path before the callback was called without `base_path`.
    full_path: (
        Path  #: The chosen output path before the callback was called with `base_path`.
    )


def _cast_base_file_path_cb(
    file_path_cb: Callable[[BasePartitionContext], Path | str | IO[bytes] | IO[str]]
    | None,
) -> Callable[[BasePartitionContext], Path | str | IO[bytes] | IO[str]] | None:
    if file_path_cb is None:
        return None
    return lambda ctx: file_path_cb(
        BasePartitionContext(
            file_idx=ctx.file_idx,
            file_path=Path(ctx.file_path),
            full_path=Path(ctx.full_path),
        )
    )


def _cast_keyed_file_path_cb(
    file_path_cb: Callable[[KeyedPartitionContext], Path | str | IO[bytes] | IO[str]]
    | None,
) -> Callable[[KeyedPartitionContext], Path | str | IO[bytes] | IO[str]] | None:
    if file_path_cb is None:
        return None
    return lambda ctx: file_path_cb(
        KeyedPartitionContext(
            file_idx=ctx.file_idx,
            part_idx=ctx.part_idx,
            in_part_idx=ctx.in_part_idx,
            keys=[
                KeyedPartition(
                    name=kv.name, str_value=kv.str_value, raw_value=kv.raw_value
                )
                for kv in ctx.keys
            ],
            file_path=Path(ctx.file_path),
            full_path=Path(ctx.full_path),
        )
    )


def _prepare_per_partition_sort_by(
    e: str | Expr | Iterable[str | Expr] | None,
) -> list[PyExpr] | None:
    def prepare_one(v: str | Expr) -> PyExpr:
        if isinstance(v, str):
            return col(v)._pyexpr
        elif isinstance(v, Expr):
            return v._pyexpr
        else:
            msg = f"cannot do a per partition sort by for {v!r}"
            raise TypeError(msg)

    if e is None:
        return None
    elif isinstance(e, str):
        return [col(e)._pyexpr]
    elif isinstance(e, Expr):
        return [e._pyexpr]
    elif isinstance(e, Iterable):
        return [prepare_one(v) for v in e]
    else:
        msg = f"cannot do a per partition sort by for {e!r}"
        raise TypeError(msg)


def _prepare_finish_callback(
    f: Callable[[DataFrame], None] | None,
) -> Callable[[PyDataFrame], None] | None:
    if f is None:
        return None

    def cb(pydf: PyDataFrame) -> None:
        nonlocal f
        f(DataFrame._from_pydf(pydf))

    return cb


class PartitionMaxSize(PartitioningScheme):
    """
    Partitioning scheme to write files with a maximum size.

    This partitioning scheme generates files that have a given maximum size. If
    the size reaches the maximum size, it is closed and a new file is opened.

    .. warning::
        This functionality is currently considered **unstable**. It may be
        changed at any point without it being considered a breaking change.

    Parameters
    ----------
    base_path
        The base path for the output files.
    file_path
        A callback to register or modify the output path for each partition
        relative to the `base_path`. The callback provides a
        :class:`polars.io.partition.BasePartitionContext` that contains information
        about the partition.

        If no callback is given, it defaults to `{ctx.file_idx}.{EXT}`.
    max_size : int
        The maximum size in rows of each of the generated files.
    per_partition_sort_by
        Columns or expressions to sort over within each partition.

        Note that this might increase the memory consumption needed for each partition.
    finish_callback
        A callback that gets called when the query finishes successfully.

        For parquet files, the callback is given a dataframe with metrics about all
        files written files.

    Examples
    --------
    Split a parquet file by over smaller CSV files with 100 000 rows each:

    >>> pl.scan_parquet("/path/to/file.parquet").sink_csv(
    ...     pl.PartitionMax("./out/", max_size=100_000),
    ... )  # doctest: +SKIP

    See Also
    --------
    PartitionByKey
    PartitionParted
    polars.io.partition.BasePartitionContext
    """

    def __init__(
        self,
        base_path: str | Path,
        *,
        file_path: Callable[[BasePartitionContext], Path | str | IO[bytes] | IO[str]]
        | None = None,
        max_size: int,
        per_partition_sort_by: str | Expr | Iterable[str | Expr] | None = None,
        finish_callback: Callable[[DataFrame], None] | None = None,
    ) -> None:
        issue_unstable_warning("partitioning strategies are considered unstable.")
        super().__init__(
            PyPartitioning.new_max_size(
                base_path=base_path,
                file_path_cb=_cast_base_file_path_cb(file_path),
                max_size=max_size,
                per_partition_sort_by=_prepare_per_partition_sort_by(
                    per_partition_sort_by
                ),
                finish_callback=_prepare_finish_callback(finish_callback),
            )
        )


def _lower_by(
    by: str | Expr | Sequence[str | Expr] | Mapping[str, Expr],
) -> list[PyExpr]:
    def to_expr(i: str | Expr) -> Expr:
        if isinstance(i, str):
            return col(i)
        else:
            return i

    lowered_by: list[PyExpr]
    if isinstance(by, str):
        lowered_by = [col(by)._pyexpr]
    elif isinstance(by, Expr):
        lowered_by = [by._pyexpr]
    elif isinstance(by, Sequence):
        lowered_by = [to_expr(e)._pyexpr for e in by]
    elif isinstance(by, Mapping):
        lowered_by = [e.alias(n)._pyexpr for n, e in by.items()]
    else:
        msg = "invalid `by` type"
        raise TypeError(msg)

    return lowered_by


class PartitionByKey(PartitioningScheme):
    """
    Partitioning scheme to write files split by the values of keys.

    This partitioning scheme generates an arbitrary amount of files splitting
    the data depending on what the value is of key expressions.

    The amount of files that can be written is not limited. However, when
    writing beyond a certain amount of files, the data for the remaining
    partitions is buffered before writing to the file.

    .. warning::
        This functionality is currently considered **unstable**. It may be
        changed at any point without it being considered a breaking change.

    Parameters
    ----------
    base_path
        The base path for the output files.

        Use the `mkdir` option on the `sink_*` methods to ensure directories in
        the path are created.
    file_path
        A callback to register or modify the output path for each partition
        relative to the `base_path`. The callback provides a
        :class:`polars.io.partition.KeyedPartitionContext` that contains information
        about the partition.

        If no callback is given, it defaults to
        `{ctx.keys.hive_dirs()}/{ctx.in_part_idx}.{EXT}`.
    by
        The expressions to partition by.
    include_key : bool
        Whether to include the key columns in the output files.
    per_partition_sort_by
        Columns or expressions to sort over within each partition.

        Note that this might increase the memory consumption needed for each partition.
    finish_callback
        A callback that gets called when the query finishes successfully.

        For parquet files, the callback is given a dataframe with metrics about all
        files written files.

    Examples
    --------
    Split into a hive-partitioning style partition:

    >>> (
    ...     pl.LazyFrame(
    ...         {"a": [1, 2, 3], "b": [5, 7, 9], "c": ["A", "B", "C"]}
    ...     ).sink_parquet(
    ...         pl.PartitionByKey(
    ...             "./out/",
    ...             by=["a", "b"],
    ...             include_key=False,
    ...         ),
    ...         mkdir=True,
    ...     )
    ... )  # doctest: +SKIP

    Split a parquet file by a column `year` into CSV files:

    >>> pl.scan_parquet("/path/to/file.parquet").sink_csv(
    ...     PartitionByKey(
    ...         "./out/",
    ...         file_path=lambda ctx: f"year={ctx.keys[0].str_value}.csv",
    ...         by="year",
    ...     ),
    ... )  # doctest: +SKIP

    See Also
    --------
    PartitionMaxSize
    PartitionParted
    polars.io.partition.KeyedPartitionContext
    """

    def __init__(
        self,
        base_path: str | Path,
        *,
        file_path: Callable[[KeyedPartitionContext], Path | str | IO[bytes] | IO[str]]
        | None = None,
        by: str | Expr | Sequence[str | Expr] | Mapping[str, Expr],
        include_key: bool = True,
        per_partition_sort_by: str | Expr | Iterable[str | Expr] | None = None,
        finish_callback: Callable[[DataFrame], None] | None = None,
    ) -> None:
        issue_unstable_warning("partitioning strategies are considered unstable.")

        lowered_by = _lower_by(by)
        super().__init__(
            PyPartitioning.new_by_key(
                base_path=base_path,
                file_path_cb=_cast_keyed_file_path_cb(file_path),
                by=lowered_by,
                include_key=include_key,
                per_partition_sort_by=_prepare_per_partition_sort_by(
                    per_partition_sort_by
                ),
                finish_callback=_prepare_finish_callback(finish_callback),
            )
        )


class PartitionParted(PartitioningScheme):
    """
    Partitioning scheme to split parted dataframes.

    This is a specialized version of :class:`PartitionByKey`. Where as
    :class:`PartitionByKey` accepts data in any order, this scheme expects the input
    data to be pre-grouped or pre-sorted. This scheme suffers a lot less overhead than
    :class:`PartitionByKey`, but may not be always applicable.

    Each new value of the key expressions starts a new partition, therefore repeating
    the same value multiple times may overwrite previous partitions.

    .. warning::
        This functionality is currently considered **unstable**. It may be
        changed at any point without it being considered a breaking change.

    Parameters
    ----------
    base_path
        The base path for the output files.

        Use the `mkdir` option on the `sink_*` methods to ensure directories in
        the path are created.
    file_path
        A callback to register or modify the output path for each partition
        relative to the `base_path`.The callback provides a
        :class:`polars.io.partition.KeyedPartitionContext` that contains information
        about the partition.

        If no callback is given, it defaults to
        `{ctx.keys.hive_dirs()}/{ctx.in_part_idx}.{EXT}`.
    by
        The expressions to partition by.
    include_key : bool
        Whether to include the key columns in the output files.
    per_partition_sort_by
        Columns or expressions to sort over within each partition.

        Note that this might increase the memory consumption needed for each partition.
    finish_callback
        A callback that gets called when the query finishes successfully.

        For parquet files, the callback is given a dataframe with metrics about all
        files written files.

    Examples
    --------
    Split a parquet file by a column `year` into CSV files:

    >>> pl.scan_parquet("/path/to/file.parquet").sink_csv(
    ...     pl.PartitionParted("./out/", by="year"),
    ...     mkdir=True,
    ... )  # doctest: +SKIP

    See Also
    --------
    PartitionMaxSize
    PartitionByKey
    polars.io.partition.KeyedPartitionContext
    """

    def __init__(
        self,
        base_path: str | Path,
        *,
        file_path: Callable[[KeyedPartitionContext], Path | str | IO[bytes] | IO[str]]
        | None = None,
        by: str | Expr | Sequence[str | Expr] | Mapping[str, Expr],
        include_key: bool = True,
        per_partition_sort_by: str | Expr | Iterable[str | Expr] | None = None,
        finish_callback: Callable[[DataFrame], None] | None = None,
    ) -> None:
        issue_unstable_warning("partitioning strategies are considered unstable.")

        lowered_by = _lower_by(by)
        super().__init__(
            PyPartitioning.new_by_key(
                base_path=base_path,
                file_path_cb=_cast_keyed_file_path_cb(file_path),
                by=lowered_by,
                include_key=include_key,
                per_partition_sort_by=_prepare_per_partition_sort_by(
                    per_partition_sort_by
                ),
                finish_callback=_prepare_finish_callback(finish_callback),
            )
        )
