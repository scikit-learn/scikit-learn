from datetime import tzinfo

import numpy as np

from pandas._typing import npt

def format_array_from_datetime(
    values: npt.NDArray[np.int64],
    tz: tzinfo | None = ...,
    format: str | None = ...,
    na_rep: object = ...,
) -> npt.NDArray[np.object_]: ...
def array_with_unit_to_datetime(
    values: np.ndarray,
    unit: str,
    errors: str = ...,
) -> tuple[np.ndarray, tzinfo | None]: ...
def array_to_datetime(
    values: npt.NDArray[np.object_],
    errors: str = ...,
    dayfirst: bool = ...,
    yearfirst: bool = ...,
    utc: bool = ...,
    require_iso8601: bool = ...,
    allow_mixed: bool = ...,
) -> tuple[np.ndarray, tzinfo | None]: ...

# returned ndarray may be object dtype or datetime64[ns]
