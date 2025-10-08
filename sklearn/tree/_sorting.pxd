from sklearn.utils._typedefs cimport float32_t, float64_t, intp_t, uint8_t, int32_t, uint32_t


ctypedef fused floating_t:
    float32_t
    float64_t


cdef inline void sort(floating_t* feature_values, intp_t* samples, intp_t n) noexcept nogil
