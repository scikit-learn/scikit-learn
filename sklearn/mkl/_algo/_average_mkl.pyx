import numpy as np

from .._utils import kernel_generator

from ...utils._typedefs cimport float64_t, intp_t
from ._utils cimport combine_kernels


def learn(
    object X,
    object y,
    object svm,
    object kernels,
    object kernels_scope,
    object kernels_params,
    const bint precomputed_kernels,
    const intp_t n_kernels,
    const intp_t n_samples,
    const int verbose=0,
):
    cdef const float64_t[::1] d = np.full(
        n_kernels,
        1.0 / n_kernels,
        dtype=np.float64,
    )

    svm.fit(
        combine_kernels(
            d,
            kernel_generator(
                X,
                kernels,
                kernels_scope,
                kernels_params,
                precomputed_kernels,
            ),
            n_samples
        ),
        y,
    )

    return d.base, svm
