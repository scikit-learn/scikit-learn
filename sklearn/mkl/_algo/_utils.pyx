import numpy as np
from ...utils._typedefs cimport float64_t, intp_t


cdef int BINARY = 0
cdef int MULTICLASS = 1
cdef int REGRESSION = 2
cdef int ONECLASS = 3


cdef float64_t[:, ::1] combine_kernels(
    const float64_t[::1] d,
    object kernels,
    const intp_t n_samples,
):
    cdef intp_t i, j, k
    cdef const float64_t[:, ::1] kernel
    cdef float64_t[:, ::1] result = np.empty((n_samples, n_samples), dtype=np.float64)

    for k, kernel in enumerate(kernels):
        for i in range(n_samples):
            for j in range(n_samples):
                if k == 0:
                    result[i, j] = d[k] * kernel[i, j]
                else:
                    result[i, j] += d[k] * kernel[i, j]

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
