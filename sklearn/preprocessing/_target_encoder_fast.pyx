cimport cython

cimport numpy as cnp
import numpy as np

cnp.import_array()

def _fit_encoding_fast(
    cnp.int64_t[:, :] X_int,
    cnp.float64_t[:] y,
    cnp.int64_t[::1] n_categories,
    double smooth,
    double y_mean,
):
    """Fit a target encoding on X_int and y."""
    cdef:
        list encodings = []
        cnp.int64_t i, j, n_cats
        int n_samples = X_int.shape[0]
        double smooth_sum = smooth * y_mean
        cnp.int64_t max_cats = np.max(n_categories)
        double[::1] current_sum = np.empty(max_cats, dtype=np.float64)
        double[::1] current_cnt = np.empty(max_cats, dtype=np.float64)

    for i in range(n_categories.shape[0]):
        n_cats = n_categories[i]
        current_encoding = np.empty(shape=n_cats, dtype=np.float64)

        for j in range(n_cats):
            current_sum[j] = smooth_sum
            current_cnt[j] = smooth

        for j in range(n_samples):
            if X_int[j, i] == -1:
                continue
            current_sum[X_int[j, i]] += y[j]
            current_cnt[X_int[j, i]] += 1.0

        for j in range(n_cats):
            current_encoding[j] = current_sum[j] / current_cnt[j]
        encodings.append(current_encoding)
    return encodings
