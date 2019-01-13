import numpy as np

Y_DTYPE = np.float32
X_DTYPE = np.float64
X_BINNED_DTYPE = np.uint8

HISTOGRAM_DTYPE = np.dtype([
    ('sum_gradients', np.float32),  # sum of sample gradients in bin
    ('sum_hessians', np.float32),  # sum of sample hessians in bin
    ('count', np.uint32),  # number of samples in bin
])
