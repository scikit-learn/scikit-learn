from __future__ import annotations

from typing import TYPE_CHECKING

import polars._plr as plr
from polars.lazyframe.opt_flags import DEFAULT_QUERY_OPT_FLAGS

if TYPE_CHECKING:
    from polars import LazyFrame, QueryOptFlags


def prepare_cloud_plan(
    lf: LazyFrame,
    *,
    allow_local_scans: bool,
    optimizations: QueryOptFlags = DEFAULT_QUERY_OPT_FLAGS,
) -> bytes:
    """
    Prepare the given LazyFrame for execution on Polars Cloud.

    Parameters
    ----------
    lf
        The LazyFrame to prepare.
    allow_local_scans
        Whether or not to allow local scans in the plan.
    optimizations
        Optimizations to enable or disable in the query optimizer.

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
    optimizations = optimizations.__copy__()
    pylf = lf._ldf.with_optimizations(optimizations._pyoptflags)
    return plr.prepare_cloud_plan(pylf, allow_local_scans=allow_local_scans)
