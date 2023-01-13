cimport cython
from libc.math cimport pow
from libc.string cimport memset

cimport numpy as cnp
import numpy as np

cnp.import_array()

ctypedef fused INT_DTYPE:
    cnp.int64_t
    cnp.int32_t

def _fit_encoding_fast(
    INT_DTYPE[:, ::1] X_int,
    cnp.float64_t[:] y,
    cnp.int64_t[::1] n_categories,
    double smooth,
    double y_mean,
):
    """Fit a target encoding on X_int and y."""
    cdef:
        list encodings = []
        cnp.int64_t i, j, n_cats
        INT_DTYPE X_int_tmp
        int n_samples = X_int.shape[0]
        int n_features = X_int.shape[1]
        double smooth_sum = smooth * y_mean
        cnp.int64_t max_cats = np.max(n_categories)
        double[::1] sums = np.empty(max_cats, dtype=np.float64)
        double[::1] counts = np.empty(max_cats, dtype=np.float64)

    for j in range(n_features):
        n_cats = n_categories[j]

        for i in range(n_cats):
            sums[i] = smooth_sum
            counts[i] = smooth

        for i in range(n_samples):
            X_int_tmp = X_int[i, j]
            # -1 are unknown categories, which are not counted
            if X_int_tmp == -1:
                continue
            sums[X_int_tmp] += y[i]
            counts[X_int_tmp] += 1.0

        current_encoding = np.empty(shape=n_cats, dtype=np.float64)
        for i in range(n_cats):
            current_encoding[i] = sums[i] / counts[i]
        encodings.append(current_encoding)
    return encodings
