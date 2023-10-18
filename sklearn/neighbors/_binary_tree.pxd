from ..utils._typedefs cimport float32_t, float64_t, intp_t

cdef struct NodeData_t:
    intp_t idx_start
    intp_t idx_end
    intp_t is_leaf
    float64_t radius
