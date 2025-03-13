import numpy as np

from ...utils._typedefs cimport float64_t, intp_t


def learn(
    const float64_t[:, :, ::1] kernels,
):
    cdef intp_t M = kernels.shape[0]

    return np.full(M, 1.0 / M, dtype=np.float64)
