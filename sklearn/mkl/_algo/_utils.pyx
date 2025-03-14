import numpy as np
from ...utils._typedefs cimport float64_t, intp_t


cdef int BINARY = 0
cdef int MULTICLASS = 1
cdef int REGRESSION = 2
cdef int ONECLASS = 3


cdef float64_t[:, ::1] combine_kernels(
    const float64_t[::1] d,
    object kernels,
):
    cdef intp_t k
    cdef const float64_t[:, ::1] kernel
    cdef float64_t[:, ::1] result

    for k, kernel in enumerate(kernels):
        if k == 0:
            result = np.multiply(d[k], kernel)
        else:
            result = np.add(result, np.multiply(d[k], kernel))

    return result


cdef float64_t[::1] fix_precision(
    float64_t[::1] d,
    const double tol,
):
    cdef intp_t i

    for i in np.where(np.less(d, tol))[0]:
        d[i] = 0
    d = d / np.sum(d)

    return d
