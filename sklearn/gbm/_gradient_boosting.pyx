# cython: profile=True
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3

cimport cython
from cython.parallel import prange
import numpy as np
cimport numpy as np

from .types import Y_DTYPE
from .types cimport Y_DTYPE_C


def _update_raw_predictions(Y_DTYPE_C [:] raw_predictions, grower):
    cdef:
        unsigned int [:] starts
        unsigned int [:] stops
        unsigned int [:] partition
        Y_DTYPE_C [:] values
        list leaves

    leaves = grower.finalized_leaves
    starts = np.array([leaf.start for leaf in leaves], dtype=np.uint32)
    stops = np.array([leaf.stop for leaf in leaves], dtype=np.uint32)
    values = np.array([leaf.value for leaf in leaves], dtype=Y_DTYPE)
    partition = grower.splitting_context.partition

    _update_raw_predictions_helper(raw_predictions, starts, stops, partition,
                                   values)

cdef void _update_raw_predictions_helper(
    Y_DTYPE_C [:] raw_predictions,
    unsigned int [:] starts,
    unsigned int [:] stops,
    unsigned int [:] partition,
    Y_DTYPE_C [:] values) nogil:

    cdef:
        int sample_idx
        int leaf_idx
        int n_leaves

    n_leaves = starts.shape[0]
    for leaf_idx in prange(n_leaves):
        for sample_idx in range(starts[leaf_idx], stops[leaf_idx]):
            raw_predictions[sample_idx] += values[leaf_idx]
