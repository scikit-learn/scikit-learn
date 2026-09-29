"""Single-pass per-feature mean and variance of dense arrays."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from cython cimport floating

import numpy as np

# Number of values of a feature accumulated around the same center before
# being merged into the running statistics (see `_dense_mean_variance_axis0`).
# The merge costs a handful of flops per feature and tile, so it is negligible
# as long as tiles are not too small.
cdef Py_ssize_t TILE_SIZE = 4096

# Maximum number of features processed together for C-contiguous inputs: this
# keeps the per-feature accumulators of a strip in the L1/L2 caches.
cdef Py_ssize_t MAX_STRIP_WIDTH = 1024

# Number of independent accumulators used in the reductions of the
# F-contiguous kernel. They break the dependency chain of the floating point
# additions (which the compiler is not allowed to reorder by itself), letting
# the CPU pipeline them and the compiler vectorize them.
cdef enum:
    N_ACC = 8


def _dense_mean_variance_axis0(X):
    """Compute per-feature mean and variance of a dense 2D array.

    Equivalent to `(np.mean(X, axis=0), np.var(X, axis=0))`, but each value of X
    is read only once and all accumulations are done in float64, whatever the
    input dtype.

    Algorithm: the values of each feature are processed by tiles of at most
    `TILE_SIZE` samples. The values of a tile are accumulated around a center
    `a` (the running mean of the previous tiles, or the first sample for the
    first tile)::

        c = sum(x - a)      q = sum((x - a) ** 2)

    from which the tile mean `a + c / n` and sum of squared deviations
    `M2 = q - c ** 2 / n` are derived. The cancellation in `M2` is governed by
    `(tile_mean - a) ** 2 / tile_variance`, which stays small because `a`
    follows the running mean. Contrary to `E[X²] - E[X]²`, this is therefore
    stable for features with a large offset relative to their standard
    deviation. Tile statistics are then merged into the running statistics
    with the pairwise update of Chan, Golub and LeVeque. The running mean is
    stored relative to the first sample, so that merges operate on quantities
    of the order of the standard deviation rather than of the offset.

    Parameters
    ----------
    X : ndarray of shape (n_samples, n_features), dtype={np.float32, np.float64}
        Input data, without NaNs. C- and F-contiguous arrays are processed without
        any copy; other memory layouts are copied into a C-contiguous array.

    Returns
    -------
    mean : ndarray of shape (n_features,), dtype=np.float64
        Per-feature mean.

    variance : ndarray of shape (n_features,), dtype=np.float64
        Per-feature population variance (ddof=0).

    References
    ----------
    T. Chan, G. Golub, R. LeVeque. Algorithms for computing the sample
    variance: analysis and recommendations, The American Statistician,
    Vol. 37, No. 3, pp. 242-247
    """
    if X.ndim != 2:
        raise ValueError(f"Expected a 2D array, got {X.ndim}D array instead.")
    if X.dtype not in (np.float32, np.float64):
        raise TypeError(f"Expected float32 or float64 data, got {X.dtype} instead.")

    n_samples, n_features = X.shape
    if n_samples == 0:
        return np.full(n_features, np.nan), np.full(n_features, np.nan)

    # First sample of each feature: reference point for the running mean.
    shift = X[0].astype(np.float64)
    # Running mean relative to `shift` and running sum of squared deviations.
    mean = np.zeros(n_features, dtype=np.float64)
    m2 = np.zeros(n_features, dtype=np.float64)

    # F-contiguous check first: a single column array is both C and F
    # contiguous, and the F-contiguous kernel is the one suited to it.
    if X.flags.f_contiguous:
        _mean_m2_f_contiguous(X, shift, mean, m2)
    else:
        _mean_m2_c_contiguous(np.ascontiguousarray(X), shift, mean, m2)

    mean += shift
    # Rounding errors can make M2 slightly negative only when it is almost 0.
    np.maximum(m2, 0, out=m2)
    m2 /= n_samples
    return mean, m2


cdef inline void _merge_tile(
    double* mean,
    double* m2,
    double n_seen,
    double center,
    double dev_sum,
    double sq_dev_sum,
    double n_tile,
) noexcept nogil:
    """Merge the statistics of a tile into the running statistics.

    `mean` and `center` are relative to the same shift. `dev_sum` and
    `sq_dev_sum` are the sums of the deviations of the tile values from
    `center`, and of their squares.
    """
    cdef double tile_dev = dev_sum / n_tile
    cdef double tile_m2 = sq_dev_sum - dev_sum * tile_dev
    cdef double n_total = n_seen + n_tile
    # tile_mean - mean, without cancellation since center ~= mean (up to the
    # rounding of the absolute center).
    cdef double delta = (center - mean[0]) + tile_dev
    mean[0] += delta * (n_tile / n_total)
    m2[0] += tile_m2 + delta * delta * (n_seen * n_tile / n_total)


def _mean_m2_c_contiguous(
    const floating[:, ::1] X,
    const double[::1] shift,
    double[::1] mean,
    double[::1] m2,
):
    cdef:
        Py_ssize_t n_samples = X.shape[0]
        Py_ssize_t n_features = X.shape[1]
        Py_ssize_t strip_width = min(n_features, MAX_STRIP_WIDTH)
        Py_ssize_t n_row_tiles = (n_samples + TILE_SIZE - 1) // TILE_SIZE
        Py_ssize_t n_col_strips = (n_features + strip_width - 1) // strip_width
        Py_ssize_t row_tile, col_strip, row_start, row_stop, col_start, width
        Py_ssize_t i, j
        double n_seen, n_tile, d
        const floating* row
        # Per-feature absolute center, sum of deviations and of their squares.
        double[::1] center_buf = np.empty(strip_width, dtype=np.float64)
        double[::1] dev_buf = np.empty(strip_width, dtype=np.float64)
        double[::1] sq_dev_buf = np.empty(strip_width, dtype=np.float64)
        double* center = &center_buf[0]
        double* dev = &dev_buf[0]
        double* sq_dev = &sq_dev_buf[0]

    with nogil:
        # Column strips as outer loop: a strip is a contiguous block of memory
        # when the array is not wider than MAX_STRIP_WIDTH, the common case.
        # Otherwise, the row segments of a strip are read in sequence, which
        # hardware prefetchers handle well.
        for col_strip in range(n_col_strips):
            col_start = col_strip * strip_width
            width = min(strip_width, n_features - col_start)
            for row_tile in range(n_row_tiles):
                row_start = row_tile * TILE_SIZE
                row_stop = min(row_start + TILE_SIZE, n_samples)
                n_seen = <double> row_start
                n_tile = <double> (row_stop - row_start)

                for j in range(width):
                    center[j] = shift[col_start + j] + mean[col_start + j]
                    dev[j] = 0.0
                    sq_dev[j] = 0.0

                # Hot loop: contiguous accesses over features and no
                # reduction, so that the compiler vectorizes it.
                for i in range(row_start, row_stop):
                    row = &X[i, col_start]
                    for j in range(width):
                        d = row[j] - center[j]
                        dev[j] += d
                        sq_dev[j] += d * d

                for j in range(width):
                    _merge_tile(
                        &mean[col_start + j],
                        &m2[col_start + j],
                        n_seen,
                        # Exact (Sterbenz lemma) when the offset dominates.
                        center[j] - shift[col_start + j],
                        dev[j],
                        sq_dev[j],
                        n_tile,
                    )


def _mean_m2_f_contiguous(
    const floating[::1, :] X,
    const double[::1] shift,
    double[::1] mean,
    double[::1] m2,
):
    cdef:
        Py_ssize_t n_samples = X.shape[0]
        Py_ssize_t n_features = X.shape[1]
        Py_ssize_t n_tiles = (n_samples + TILE_SIZE - 1) // TILE_SIZE
        Py_ssize_t tile, chunk, start, n_tile_int, n_chunks, i, j, k
        double n_seen, n_tile, center, dev, sq_dev, d
        double acc_dev[N_ACC]
        double acc_sq_dev[N_ACC]
        const floating* col

    with nogil:
        for j in range(n_features):
            for tile in range(n_tiles):
                start = tile * TILE_SIZE
                n_tile_int = min(TILE_SIZE, n_samples - start)
                n_chunks = n_tile_int // N_ACC
                n_seen = <double> start
                n_tile = <double> n_tile_int
                col = &X[start, j]
                center = shift[j] + mean[j]

                for k in range(N_ACC):
                    acc_dev[k] = 0.0
                    acc_sq_dev[k] = 0.0
                for chunk in range(n_chunks):
                    i = chunk * N_ACC
                    for k in range(N_ACC):
                        d = col[i + k] - center
                        acc_dev[k] += d
                        acc_sq_dev[k] += d * d
                dev = 0.0
                sq_dev = 0.0
                for k in range(N_ACC):
                    dev += acc_dev[k]
                    sq_dev += acc_sq_dev[k]
                for i in range(n_chunks * N_ACC, n_tile_int):
                    d = col[i] - center
                    dev += d
                    sq_dev += d * d

                _merge_tile(
                    &mean[j], &m2[j], n_seen, center - shift[j], dev, sq_dev, n_tile
                )
