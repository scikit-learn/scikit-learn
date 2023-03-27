# _typedefs is a declaration only module
# The functions implemented here are for testing purpose.


import numpy as np


ctypedef fused type_t:
    uint8_t
    intp_t
    float32_t
    float64_t
    int32_t
    int64_t


def make_array_from_typed_val(type_t val):
    cdef type_t[:] val_view = <type_t[:1]>&val
    return np.asarray(val_view)
