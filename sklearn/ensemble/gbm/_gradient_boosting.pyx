# cython: profile=True
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3
cimport cython

import numpy as np
cimport numpy as np

ctypedef np.npy_float32 NPY_Y_DTYPE
ctypedef np.npy_uint8 NPY_X_BINNED_DTYPE

def _update_raw_predictions(float [:] leaves_values,
                            list samples_indices_at_leaf,
                            NPY_Y_DTYPE [:] raw_predictions):
    """Update raw_predictions by reading the predictions of the ith tree
    directly form the leaves.

    Can only be used for predicting the training data. raw_predictions
    contains the sum of the tree values from iteration 0 to i - 1. This adds
    the predictions of the ith tree to raw_predictions.

    Parameters
    ----------
    leaves_data: list of tuples (leaf.value, leaf.sample_indices)
        The leaves data used to update raw_predictions.
    raw_predictions : array-like, shape=(n_samples,)
        The raw predictions for the training data.
    """
    cdef:
        int leaf_idx
        float val
        unsigned int [:] sample_indices

    for leaf_idx in range(leaves_values.shape[0]):
        samples_indices = samples_indices_at_leaf[leaf_idx]
        val = leaves_values[leaf_idx]
        blop(samples_indices, raw_predictions, val)

cdef void blop(unsigned int [:] samples_indices, NPY_Y_DTYPE [:] raw_predictions, float
                val):
    cdef:
        unsigned int sample_idx
    for sample_idx in samples_indices:
        raw_predictions[sample_idx] += val