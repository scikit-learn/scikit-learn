from __future__ import annotations

from typing import TYPE_CHECKING

import polars.polars as plr

if TYPE_CHECKING:
    from polars import LazyFrame


def prepare_cloud_plan(
    lf: LazyFrame,
    **optimizations: bool,
) -> bytes:
    """
    Prepare the given LazyFrame for execution on Polars Cloud.

    Parameters
    ----------
    lf
        The LazyFrame to prepare.
    **optimizations
        Optimizations to enable or disable in the query optimizer, e.g.
        `projection_pushdown=False`.

    Raises
    ------
    InvalidOperationError
        If the given LazyFrame is not eligible to be run on Polars Cloud.
        The following conditions will disqualify a LazyFrame from being eligible:

        - Contains a user-defined function
        - Scans or sinks to a local filesystem
    ComputeError
        If the given LazyFrame cannot be serialized.
    """
    pylf = lf._set_sink_optimizations(**optimizations)
    return plr.prepare_cloud_plan(pylf)
