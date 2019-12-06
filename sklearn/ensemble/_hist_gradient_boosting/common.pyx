import numpy as np

# Y_DYTPE is the dtype to which the targets y are converted to. This is also
# dtype for leaf values, gains, and sums of gradients / hessians. The gradients
# and hessians arrays are stored as floats to avoid using too much memory.
Y_DTYPE = np.float64
X_DTYPE = np.float64
X_BINNED_DTYPE = np.uint8  # hence max_bins == 256
# dtype for gradients and hessians arrays
G_H_DTYPE = np.float32

HISTOGRAM_DTYPE = np.dtype([
    ('sum_gradients', Y_DTYPE),  # sum of sample gradients in bin
    ('sum_hessians', Y_DTYPE),  # sum of sample hessians in bin
    ('count', np.uint32),  # number of samples in bin
])

PREDICTOR_RECORD_DTYPE = np.dtype([
    ('value', Y_DTYPE),
    ('count', np.uint32),
    ('feature_idx', np.uint32),
    ('threshold', X_DTYPE),
    ('missing_go_to_left', np.uint8),
    ('left', np.uint32),
    ('right', np.uint32),
    ('gain', Y_DTYPE),
    ('depth', np.uint32),
    ('is_leaf', np.uint8),
    ('bin_threshold', X_BINNED_DTYPE),
])

ALMOST_INF = 1e300  # see LightGBM AvoidInf()
