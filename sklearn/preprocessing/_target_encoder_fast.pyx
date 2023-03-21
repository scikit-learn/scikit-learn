from libc.math cimport isnan
from libcpp.vector cimport vector

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

    This implementation uses Eq 7 from [1] to compute the encoding.
    As stated in the paper, Eq 7 is the same as Eq 3.

    [1]: Micci-Barreca, Daniele. "A preprocessing scheme for high-cardinality
         categorical attributes in classification and prediction problems"
    """
    cdef:
        cnp.int64_t sample_idx, feat_idx, cat_idx, n_cats
        INT_DTYPE X_int_tmp
        int n_samples = X_int.shape[0]
        int n_features = X_int.shape[1]
        double smooth_sum = smooth * y_mean
        cnp.int64_t max_n_cats = np.max(n_categories)
        double[::1] sums = np.empty(max_n_cats, dtype=np.float64)
        double[::1] counts = np.empty(max_n_cats, dtype=np.float64)
        list encodings = []
        double[::1] current_encoding
        # Gives access to encodings without gil
        vector[double*] encoding_vec

    encoding_vec.resize(n_features)
    for feat_idx in range(n_features):
        current_encoding = np.empty(shape=n_categories[feat_idx], dtype=np.float64)
        encoding_vec[feat_idx] = &current_encoding[0]
        encodings.append(np.asarray(current_encoding))

    with nogil:
        for feat_idx in range(n_features):
            n_cats = n_categories[feat_idx]

            for cat_idx in range(n_cats):
                sums[cat_idx] = smooth_sum
                counts[cat_idx] = smooth

            for sample_idx in range(n_samples):
                X_int_tmp = X_int[sample_idx, feat_idx]
                # -1 are unknown categories, which are not counted
                if X_int_tmp == -1:
                    continue
                sums[X_int_tmp] += y[sample_idx]
                counts[X_int_tmp] += 1.0

            for cat_idx in range(n_cats):
                if counts[cat_idx] == 0:
                    encoding_vec[feat_idx][cat_idx] = y_mean
                else:
                    encoding_vec[feat_idx][cat_idx] = sums[cat_idx] / counts[cat_idx]

    return encodings


def _fit_encoding_fast_auto_smooth(
    INT_DTYPE[:, ::1] X_int,
    Y_DTYPE[:] y,
    cnp.int64_t[::1] n_categories,
    double y_mean,
    double y_variance,
):
    """Fit a target encoding on X_int and y with auto smoothing.

    This implementation uses Eq 5 and 6 from [1].

    [1]: Micci-Barreca, Daniele. "A preprocessing scheme for high-cardinality
         categorical attributes in classification and prediction problems"
    """
    cdef:
        cnp.int64_t sample_idx, feat_idx, cat_idx, n_cats
        INT_DTYPE X_int_tmp
        double diff
        int n_samples = X_int.shape[0]
        int n_features = X_int.shape[1]
        cnp.int64_t max_n_cats = np.max(n_categories)
        double[::1] means = np.empty(max_n_cats, dtype=np.float64)
        cnp.int64_t[::1] counts = np.empty(max_n_cats, dtype=np.int64)
        double[::1] sum_of_squared_diffs = np.empty(max_n_cats, dtype=np.float64)
        double lambda_
        list encodings = []
        double[::1] current_encoding
        # Gives access to encodings without gil
        vector[double*] encoding_vec

    encoding_vec.resize(n_features)
    for feat_idx in range(n_features):
        current_encoding = np.empty(shape=n_categories[feat_idx], dtype=np.float64)
        encoding_vec[feat_idx] = &current_encoding[0]
        encodings.append(np.asarray(current_encoding))

    # TODO: parallelize this with OpenMP prange. When n_features >= n_threads, it's
    # probably good to parallelize the outer loop. When n_features is too small,
    # then it would probably better to parallelize the nested loops on n_samples and
    # n_cats, but the code to handle thread-local temporary variables might be
    # significantly more complex.
    with nogil:
        for feat_idx in range(n_features):
            n_cats = n_categories[feat_idx]

            for cat_idx in range(n_cats):
                means[cat_idx] = 0.0
                counts[cat_idx] = 0
                sum_of_squared_diffs[cat_idx] = 0.0

            # first pass to compute the mean
            for sample_idx in range(n_samples):
                X_int_tmp = X_int[sample_idx, feat_idx]

                # -1 are unknown categories, which are not counted
                if X_int_tmp == -1:
                    continue
                counts[X_int_tmp] += 1
                means[X_int_tmp] += y[sample_idx]

            for cat_idx in range(n_cats):
                means[cat_idx] /= counts[cat_idx]

            # second pass to compute the sum of squared differences
            for sample_idx in range(n_samples):
                X_int_tmp = X_int[sample_idx, feat_idx]
                if X_int_tmp == -1:
                    continue
                diff = y[sample_idx] - means[X_int_tmp]
                sum_of_squared_diffs[X_int_tmp] += diff * diff

            for cat_idx in range(n_cats):
                lambda_ = (
                    y_variance * counts[cat_idx] /
                    (y_variance * counts[cat_idx] + sum_of_squared_diffs[cat_idx] /
                     counts[cat_idx])
                )
                if isnan(lambda_):
                    # A nan can happen when:
                    # 1. counts[cat_idx] == 0
                    # 2. y_variance == 0 and sum_of_squared_diffs[cat_idx] == 0
                    encoding_vec[feat_idx][cat_idx] = y_mean
                else:
                    encoding_vec[feat_idx][cat_idx] = (
                        lambda_ * means[cat_idx] + (1 - lambda_) * y_mean
                    )

    return encodings
