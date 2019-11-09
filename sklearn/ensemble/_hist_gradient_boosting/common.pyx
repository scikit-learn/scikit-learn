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


cpdef inline Y_DTYPE_C compute_value(
        Y_DTYPE_C sum_gradient,
        Y_DTYPE_C sum_hessian,
        Y_DTYPE_C lower_bound,
        Y_DTYPE_C upper_bound,
        Y_DTYPE_C l2_regularization) nogil:
    """Compute node value.

    The value is capped in the [lower_bound, upper_bound] interval to respect
    monotonic constraints. Shrinkage is ignored.

    See Equation 5 of:
    XGBoost: A Scalable Tree Boosting System, T. Chen, C. Guestrin, 2016
    https://arxiv.org/abs/1603.02754
    """

    cdef:
        Y_DTYPE_C value

    value = -sum_gradient / (sum_hessian + l2_regularization + 1e-15)

    if value < lower_bound:
        value = lower_bound
    elif value > upper_bound:
        value = upper_bound

    return value


ALMOST_INF = 1e300  # see LightGBM AvoidInf()
