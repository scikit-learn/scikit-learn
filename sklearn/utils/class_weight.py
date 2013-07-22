# Authors: Andreas Mueller
# License: BSD 3 clause

import numpy as np

from .fixes import bincount


def compute_class_weight(class_weight, classes, y_ind):
    """Estimate class weights for unbalanced datasets.

    Parameters
    ----------
    class_weight : dict, 'auto' or None
        If 'auto', class weights will be given inverse proportional
        to the frequency of the class in the data.
        If a dictionary is given, keys are classes and values
        are corresponding class weights.
        If None is given, the class weights will be uniform.

    classes : ndarray
        Array of the classes occurring in the data, as given by
        ``np.unique(y_org)`` with ``y_org`` the original class labels.

    y_ind : array-like, shape=(n_samples,), dtype=int
        Array of class indices per sample;
        0 <= y_ind[i] < n_classes for i in range(n_samples).

    Returns
    -------
    class_weight_vect : ndarray, shape=(n_classes,)
        Array with class_weight_vect[i] the weight for i-th class
        (as determined by sorting).
    """
    if class_weight is None or len(class_weight) == 0:
        # uniform class weights
        weight = np.ones(classes.shape[0], dtype=np.float64, order='C')
    elif class_weight == 'auto':
        # inversely proportional to the number of samples in the class
        counts = bincount(y_ind, minlength=len(classes))
        counts = np.maximum(counts, 1)
        weight = 1. / counts
        weight *= classes.shape[0] / np.sum(weight)
    else:
        # user-defined dictionary
        weight = np.ones(classes.shape[0], dtype=np.float64, order='C')
        if not isinstance(class_weight, dict):
            raise ValueError("class_weight must be dict, 'auto', or None,"
                             " got: %r" % class_weight)
        for c in class_weight:
            i = np.searchsorted(classes, c)
            if classes[i] != c:
                raise ValueError("Class label %d not present." % c)
            else:
                weight[i] = class_weight[c]

    return weight
