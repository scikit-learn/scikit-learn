from __future__ import annotations
import numpy as np

# list of numarray data types
integer_types: list[type] = [
    np.int8, np.uint8, np.int16, np.uint16,
    np.int32, np.uint32, np.int64, np.uint64]

float_types: list[type] = [np.float32, np.float64]

complex_types: list[type] = [np.complex64, np.complex128]

types: list[type] = integer_types + float_types
