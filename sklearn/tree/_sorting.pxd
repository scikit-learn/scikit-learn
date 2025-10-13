from sklearn.utils._typedefs cimport float32_t, float64_t, intp_t, uint8_t, int32_t, uint32_t


from cython cimport floating


cdef void sort(floating* feature_values, intp_t* samples, intp_t n) noexcept nogil
