"""
For cython types that cannot be represented precisely, closest-available
python equivalents are used, and the precise types kept as adjacent comments.
"""
from datetime import tzinfo

import numpy as np

from pandas._libs.tslibs.dtypes import Resolution
from pandas._libs.tslibs.offsets import BaseOffset
from pandas._typing import npt

def dt64arr_to_periodarr(
    stamps: npt.NDArray[np.int64],  # const int64_t[:]
    freq: int,
    tz: tzinfo | None,
) -> npt.NDArray[np.int64]: ...  # np.ndarray[np.int64, ndim=1]
def is_date_array_normalized(
    stamps: npt.NDArray[np.int64],  # const int64_t[:]
    tz: tzinfo | None = ...,
) -> bool: ...
def normalize_i8_timestamps(
    stamps: npt.NDArray[np.int64],  # const int64_t[:]
    tz: tzinfo | None,
) -> npt.NDArray[np.int64]: ...
def get_resolution(
    stamps: npt.NDArray[np.int64],  # const int64_t[:]
    tz: tzinfo | None = ...,
) -> Resolution: ...
def ints_to_pydatetime(
    arr: npt.NDArray[np.int64],  # const int64_t[:}]
    tz: tzinfo | None = ...,
    freq: str | BaseOffset | None = ...,
    fold: bool = ...,
    box: str = ...,
) -> npt.NDArray[np.object_]: ...
