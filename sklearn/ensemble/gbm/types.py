import numpy as np


Y_DTYPE = np.float32
X_DTYPE = np.float64
X_BINNED_DTYPE = np.uint8

HISTOGRAM_DTYPE = np.dtype([
    ('sum_gradients', np.float32),
    ('sum_hessians', np.float32),
    ('count', np.uint32),
])
