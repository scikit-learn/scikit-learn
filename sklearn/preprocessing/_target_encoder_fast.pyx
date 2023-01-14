cimport cython
from libc.math cimport pow, isnan
from libc.string cimport memset

cimport numpy as cnp
import numpy as np

cnp.import_array()

ctypedef fused INT_DTYPE:
    cnp.int64_t
    cnp.int32_t

ctypedef fused Y_DTYPE:
    cnp.int64_t
    cnp.int32_t
    cnp.float64_t
    cnp.float32_t

def _fit_encoding_fast(
    INT_DTYPE[:, ::1] X_int,
    Y_DTYPE[:] y,
    cnp.int64_t[::1] n_categories,
    double smooth,
    double y_mean,
):
    """Fit a target encoding on X_int and y.

    This implementation uses Equation 7 from [1] to compute the encoding.

    [1]: Micci-Barreca, Daniele. "A preprocessing scheme for high-cardinality
         categorical attributes in classification and prediction problems"
    """
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


def _fit_encoding_fast_auto_smooth(
    INT_DTYPE[:, ::1] X_int,
    Y_DTYPE[:] y,
    cnp.int64_t[::1] n_categories,
    double y_mean,
    double y_variance,
):
    """Fit a target encoding on X_int and y with auto smoothing.

    This implementation uses Equation 6 from [1] to compute `lambda`.

    [1]: Micci-Barreca, Daniele. "A preprocessing scheme for high-cardinality
         categorical attributes in classification and prediction problems"
    """
    cdef:
        list encodings = []
        cnp.int64_t i, j, n_cats
        INT_DTYPE X_int_tmp
        int n_samples = X_int.shape[0]
        int n_features = X_int.shape[1]
        cnp.int64_t max_cats = np.max(n_categories)
        double[::1] means = np.empty(max_cats, dtype=np.float64)
        double[::1] counts = np.empty(max_cats, dtype=np.float64)
        double[::1] sum_of_squared_diffs = np.empty(max_cats, dtype=np.float64)
        double lambda_

    for j in range(n_features):
        n_cats = n_categories[j]

        for i in range(n_cats):
            means[i] = 0.0
            counts[i] = 0.0
            sum_of_squared_diffs[i] = 0.0

        # first pass to compute the mean
        for i in range(n_samples):
            X_int_tmp = X_int[i, j]

            # -1 are unknown categories, which are not counted
            if X_int_tmp == -1:
                continue
            counts[X_int_tmp] += 1.0
            means[X_int_tmp] += (y[i] - means[X_int_tmp]) / counts[X_int_tmp]

        # second pass to compute the sum of squared differences
        for i in range(n_samples):
            X_int_tmp = X_int[i, j]
            if X_int_tmp == -1:
                continue
            sum_of_squared_diffs[X_int_tmp] += pow(y[i] - means[X_int_tmp], 2.0)

        current_encoding = np.empty(shape=n_cats, dtype=np.float64)
        for i in range(n_cats):
            lambda_ = y_variance * counts[i] / (y_variance * counts[i] + sum_of_squared_diffs[i] / counts[i])
            if isnan(lambda_):
                current_encoding[i] = y_mean
            else:
                current_encoding[i] = lambda_ * means[i] + (1 - lambda_) * y_mean
        encodings.append(current_encoding)
    return encodings
