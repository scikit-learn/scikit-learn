from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from polars._typing import TimeUnit

# Number of rows to scan by default when inferring datatypes
N_INFER_DEFAULT = 100

DTYPE_TEMPORAL_UNITS: frozenset[TimeUnit] = frozenset(["ns", "us", "ms"])
