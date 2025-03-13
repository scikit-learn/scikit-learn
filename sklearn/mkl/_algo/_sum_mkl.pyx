import numpy as np

from ...utils._typedefs cimport float64_t


def learn(
    const float64_t[:, :, ::1] kernels,
):
    return np.full(kernels.shape[0], 1, dtype=np.float64)
