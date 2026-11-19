from libcpp.vector cimport vector

from sklearn.utils._typedefs cimport (
    float32_t,
    float64_t,
    int32_t,
    int64_t,
    intp_t,
    inlinable_isnan,
)

import numpy as np


ctypedef fused INT_DTYPE:
    int64_t
    int32_t

ctypedef fused Y_DTYPE:
    int64_t
    int32_t
    float64_t
    float32_t


def _fit_encoding_fast(
    INT_DTYPE[::1, :] X_int,
    const Y_DTYPE[:] y,
    int64_t[::1] n_categories,
    double smooth,
    double y_mean,
    const intp_t[::1] X_indices=None,
    const float64_t[:] sample_weight=None,
):
    """Fit a target encoding on X_int and y.

    This implementation uses Eq 7 from [1] to compute the encoding.
    As stated in the paper, Eq 7 is the same as Eq 3.

    When `sample_weight` is given, the per-category sum of targets and the
    per-category count are replaced by their weighted counterparts. `y` and
    `sample_weight` are indexed by position in the (possibly subsampled) training
    set, i.e. the same way.

    [1]: Micci-Barreca, Daniele. "A preprocessing scheme for high-cardinality
         categorical attributes in classification and prediction problems"
    """
    cdef:
        intp_t sample_idx, row_idx
        int64_t feat_idx, cat_idx, n_cats
        INT_DTYPE X_int_tmp
        bint use_X_indices = X_indices is not None
        bint use_sample_weight = sample_weight is not None
        double weight = 1.0
        intp_t n_samples
        int n_features = X_int.shape[1]
        double smooth_sum = smooth * y_mean
        int64_t max_n_cats = np.max(n_categories)
        double[::1] sums = np.empty(max_n_cats, dtype=np.float64)
        double[::1] weighted_counts = np.empty(max_n_cats, dtype=np.float64)
        list encodings = []
        double[::1] current_encoding
        # Gives access to encodings without gil
        vector[double*] encoding_vec

    if use_X_indices:
        n_samples = X_indices.shape[0]
    else:
        n_samples = X_int.shape[0]

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
                weighted_counts[cat_idx] = smooth

            for sample_idx in range(n_samples):
                if use_X_indices:
                    row_idx = X_indices[sample_idx]
                else:
                    row_idx = sample_idx
                X_int_tmp = X_int[row_idx, feat_idx]
                # -1 are unknown categories, which are not counted
                if X_int_tmp == -1:
                    continue
                if use_sample_weight:
                    weight = sample_weight[sample_idx]
                sums[X_int_tmp] += weight * y[sample_idx]
                weighted_counts[X_int_tmp] += weight

            for cat_idx in range(n_cats):
                # A zero count happens when smooth == 0 and the category has no
                # samples (or only zero-weight samples) in the training data.
                if weighted_counts[cat_idx] == 0:
                    encoding_vec[feat_idx][cat_idx] = y_mean
                else:
                    encoding_vec[feat_idx][cat_idx] = (
                        sums[cat_idx] / weighted_counts[cat_idx]
                    )

    return encodings


def _fit_encoding_fast_auto_smooth(
    INT_DTYPE[::1, :] X_int,
    const Y_DTYPE[:] y,
    int64_t[::1] n_categories,
    double y_mean,
    double y_variance,
    const intp_t[::1] X_indices=None,
    const float64_t[:] sample_weight=None,
):
    """Fit a target encoding on X_int and y with auto smoothing.

    This implementation uses Eq 5 and 6 from [1].

    When `sample_weight` is given, the per-category count, mean and sum of
    squared differences are replaced by their weighted counterparts, so that
    integer weights are equivalent to repeating samples. `y_mean` and
    `y_variance` are expected to be the weighted mean and variance of `y`.

    [1]: Micci-Barreca, Daniele. "A preprocessing scheme for high-cardinality
         categorical attributes in classification and prediction problems"
    """
    cdef:
        intp_t sample_idx, row_idx
        int64_t feat_idx, cat_idx, n_cats
        INT_DTYPE X_int_tmp
        double diff
        bint use_X_indices = X_indices is not None
        bint use_sample_weight = sample_weight is not None
        double weight = 1.0
        intp_t n_samples
        int n_features = X_int.shape[1]
        int64_t max_n_cats = np.max(n_categories)
        double[::1] means = np.empty(max_n_cats, dtype=np.float64)
        double[::1] weighted_counts = np.empty(max_n_cats, dtype=np.float64)
        double[::1] sum_of_squared_diffs = np.empty(max_n_cats, dtype=np.float64)
        double lambda_
        list encodings = []
        double[::1] current_encoding
        # Gives access to encodings without gil
        vector[double*] encoding_vec

    if use_X_indices:
        n_samples = X_indices.shape[0]
    else:
        n_samples = X_int.shape[0]

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
                weighted_counts[cat_idx] = 0.0
                sum_of_squared_diffs[cat_idx] = 0.0

            # first pass to compute the (weighted) mean
            for sample_idx in range(n_samples):
                if use_X_indices:
                    row_idx = X_indices[sample_idx]
                else:
                    row_idx = sample_idx
                X_int_tmp = X_int[row_idx, feat_idx]

                # -1 are unknown categories, which are not counted
                if X_int_tmp == -1:
                    continue
                if use_sample_weight:
                    weight = sample_weight[sample_idx]
                weighted_counts[X_int_tmp] += weight
                means[X_int_tmp] += weight * y[sample_idx]

            for cat_idx in range(n_cats):
                means[cat_idx] /= weighted_counts[cat_idx]

            # second pass to compute the (weighted) sum of squared differences
            for sample_idx in range(n_samples):
                if use_X_indices:
                    row_idx = X_indices[sample_idx]
                else:
                    row_idx = sample_idx
                X_int_tmp = X_int[row_idx, feat_idx]
                if X_int_tmp == -1:
                    continue
                if use_sample_weight:
                    weight = sample_weight[sample_idx]
                diff = y[sample_idx] - means[X_int_tmp]
                sum_of_squared_diffs[X_int_tmp] += weight * diff * diff

            for cat_idx in range(n_cats):
                lambda_ = (
                    y_variance * weighted_counts[cat_idx] /
                    (y_variance * weighted_counts[cat_idx] +
                     sum_of_squared_diffs[cat_idx] / weighted_counts[cat_idx])
                )
                if inlinable_isnan(lambda_):
                    # A nan can happen when:
                    # 1. weighted_counts[cat_idx] == 0 (no samples, or only
                    #    zero-weight samples, in the category)
                    # 2. y_variance == 0 and sum_of_squared_diffs[cat_idx] == 0
                    encoding_vec[feat_idx][cat_idx] = y_mean
                else:
                    encoding_vec[feat_idx][cat_idx] = (
                        lambda_ * means[cat_idx] + (1 - lambda_) * y_mean
                    )

    return encodings
