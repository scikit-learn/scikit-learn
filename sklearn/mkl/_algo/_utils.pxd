from ...utils._typedefs cimport float64_t, intp_t


cdef int BINARY
cdef int MULTICLASS
cdef int REGRESSION
cdef int ONECLASS


cdef float64_t[:, ::1] combine_kernels(
    const float64_t[::1] d,
    object kernels,
    const intp_t n_samples,
)


cdef float64_t[::1] fix_precision(
    float64_t[::1] d,
    const double tol,
)
