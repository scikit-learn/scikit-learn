from __future__ import annotations

from typing import TYPE_CHECKING

from polars._utils.unstable import unstable
from polars.io.pyarrow_dataset.anonymous_scan import _scan_pyarrow_dataset

if TYPE_CHECKING:
    from polars import LazyFrame
    from polars.dependencies import pyarrow as pa


@unstable()
def scan_pyarrow_dataset(
    source: pa.dataset.Dataset,
    *,
    allow_pyarrow_filter: bool = True,
    batch_size: int | None = None,
) -> LazyFrame:
    """
    Scan a pyarrow dataset.

    .. warning::
        This functionality is considered **unstable**. It may be changed
        at any point without it being considered a breaking change.

    This can be useful to connect to cloud or partitioned datasets.

    Parameters
    ----------
    source
        Pyarrow dataset to scan.
    allow_pyarrow_filter
        Allow predicates to be pushed down to pyarrow. This can lead to different
        results if comparisons are done with null values as pyarrow handles this
        different than polars does.
    batch_size
        The maximum row count for scanned pyarrow record batches.

    Warnings
    --------
    This method can only can push down predicates that are allowed by PyArrow
    (e.g. not the full Polars API).

    If :func:`scan_parquet` works for your source, you should use that instead.

    Notes
    -----
    When using partitioning, the appropriate `partitioning` option must be set on
    `pyarrow.dataset.dataset` before passing to Polars or the partitioned-on column(s)
    may not get passed to Polars.

    Examples
    --------
    >>> import pyarrow.dataset as ds
    >>> dset = ds.dataset("s3://my-partitioned-folder/", format="ipc")  # doctest: +SKIP
    >>> (
    ...     pl.scan_pyarrow_dataset(dset)
    ...     .filter("bools")
    ...     .select("bools", "floats", "date")
    ...     .collect()
    ... )  # doctest: +SKIP
    shape: (1, 3)
    ┌───────┬────────┬────────────┐
    │ bools ┆ floats ┆ date       │
    │ ---   ┆ ---    ┆ ---        │
    │ bool  ┆ f64    ┆ date       │
    ╞═══════╪════════╪════════════╡
    │ true  ┆ 2.0    ┆ 1970-05-04 │
    └───────┴────────┴────────────┘
    """
    return _scan_pyarrow_dataset(
        source,
        allow_pyarrow_filter=allow_pyarrow_filter,
        batch_size=batch_size,
    )
