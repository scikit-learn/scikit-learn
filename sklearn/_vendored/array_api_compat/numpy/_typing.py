from __future__ import annotations

__all__ = [
    "ndarray",
    "Device",
    "Dtype",
]

import sys
from typing import (
    Literal,
    Union,
    TYPE_CHECKING,
)

from numpy import (
    ndarray,
    dtype,
    int8,
    int16,
    int32,
    int64,
    uint8,
    uint16,
    uint32,
    uint64,
    float32,
    float64,
)

Device = Literal["cpu"]
if TYPE_CHECKING or sys.version_info >= (3, 9):
    Dtype = dtype[Union[
        int8,
        int16,
        int32,
        int64,
        uint8,
        uint16,
        uint32,
        uint64,
        float32,
        float64,
    ]]
else:
    Dtype = dtype
