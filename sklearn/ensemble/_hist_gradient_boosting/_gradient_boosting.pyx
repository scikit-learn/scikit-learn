# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3

# Author: Nicolas Hug

cimport cython
from cython.parallel import prange
import numpy as np
cimport numpy as np

from .types import Y_DTYPE
from .types cimport Y_DTYPE_C


def _update_raw_predictions(
        Y_DTYPE_C [::1] raw_predictions,  # OUT
        grower):
    """Update raw_predictions with the predictions of the newest tree.

    This is equivalent to (and much faster than):
        raw_predictions += last_estimator.predict(X_train)

    It's only possible for data X_train that is used to train the trees (it
    isn't usable for e.g. X_val).
    """
    cdef:
        unsigned int [::1] starts  # start of each leaf in partition
        unsigned int [::1] stops  # end of each leaf in partition
        Y_DTYPE_C [::1] values  # value of each leaf
        const unsigned int [::1] partition = grower.splitter.partition
        list leaves

    leaves = grower.finalized_leaves
    starts = np.array([leaf.partition_start for leaf in leaves],
                      dtype=np.uint32)
    stops = np.array([leaf.partition_stop for leaf in leaves],
                     dtype=np.uint32)
    values = np.array([leaf.value for leaf in leaves], dtype=Y_DTYPE)

    _update_raw_predictions_helper(raw_predictions, starts, stops, partition,
                                   values)


cdef inline void _update_raw_predictions_helper(
        Y_DTYPE_C [::1] raw_predictions,  # OUT
        const unsigned int [::1] starts,
        const unsigned int [::1] stops,
        const unsigned int [::1] partition,
        const Y_DTYPE_C [::1] values):

    cdef:
        unsigned int position
        int leaf_idx
        int n_leaves = starts.shape[0]

    for leaf_idx in prange(n_leaves, nogil=True):
        for position in range(starts[leaf_idx], stops[leaf_idx]):
            raw_predictions[partition[position]] += values[leaf_idx]
