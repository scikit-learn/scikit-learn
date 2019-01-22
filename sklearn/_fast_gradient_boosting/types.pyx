import numpy as np

# Y_DYTPE is the dtype to which the targets y are converted to. This is also
# the dtype for gradients, hessians, leaf values, etc. because they are all
# homogeneous to a target.
Y_DTYPE = np.float64
X_DTYPE = np.float64
X_BINNED_DTYPE = np.uint8  # hence max_bins == 256

HISTOGRAM_DTYPE = np.dtype([
    ('sum_gradients', Y_DTYPE),  # sum of sample gradients in bin
    ('sum_hessians', Y_DTYPE),  # sum of sample hessians in bin
    ('count', np.uint32),  # number of samples in bin
])
