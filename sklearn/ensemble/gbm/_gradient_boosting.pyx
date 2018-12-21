# cython: profile=True
cimport cython

import numpy as np
cimport numpy as np

ctypedef fused float_or_double:
    float
    double

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
def _update_raw_predictions__(float [:] leaves_values, list samples_indices_at_leaf, float_or_double [:] raw_predictions):
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
        unsigned int sample_idx
        unsigned int [:] sample_indices

    for leaf_idx in range(leaves_values.shape[0]):
        samples_indices = samples_indices_at_leaf[leaf_idx]
        for sample_idx in samples_indices:
            raw_predictions[sample_idx] += leaves_values[leaf_idx]