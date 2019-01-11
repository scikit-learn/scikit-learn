# cython: profile=True
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3
cimport cython

import numpy as np
cimport numpy as np

from .types import Y_DTYPE

ctypedef np.npy_float32 NPY_Y_DTYPE

def _update_raw_predictions(NPY_Y_DTYPE [:] raw_predictions, grower):
    cdef:
        unsigned int [:] starts
        unsigned int [:] stops
        unsigned int [:] partition
        NPY_Y_DTYPE [:] values
        list leaves

    leaves = grower.finalized_leaves
    starts = np.array([leaf.start for leaf in leaves], dtype=np.uint32)
    stops = np.array([leaf.stop for leaf in leaves], dtype=np.uint32)
    values = np.array([leaf.value for leaf in leaves], dtype=Y_DTYPE)
    partition = grower.splitting_context.partition

    _update_raw_predictions_helper(raw_predictions, starts, stops, partition,
                                   values)

cdef void _update_raw_predictions_helper(
    NPY_Y_DTYPE [:] raw_predictions,
    unsigned int [:] starts,
    unsigned int [:] stops,
    unsigned int [:] partition,
    NPY_Y_DTYPE [:] values) nogil:

    cdef:
        unsigned int sample_idx
        unsigned int n_leaves

    n_leaves = starts.shape[0]
    for leaf_idx in range(n_leaves):
        for sample_idx in range(starts[leaf_idx], stops[leaf_idx]):
            raw_predictions[sample_idx] += values[leaf_idx]
