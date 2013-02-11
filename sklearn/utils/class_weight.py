# Authors: Andreas Mueller
# License: Simplified BSD

import numpy as np


def compute_class_weight(class_weight, classes, y):
    """Estimate class weights for unbalanced datasets.

    Parameters
    ----------
    class_weight : dict, 'auto' or None
        If 'auto', class weights will be given inverse proportional
        to the frequency of the class in the data.
        If a dictionary is given, keys are classes and values
        are corresponding class weights.
        If None is given, the class weights will be uniform.

    classes : list
        List of the classes occuring in the data, as given by
        ``np.unique(y_org)`` with ``y_org`` the original class labels.

    y : array-like, shape=(n_samples,), dtype=int
        Array of class indices per sample;
        0 <= y[i] < n_classes for i in range(n_samples).

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
        # anti-proportional to the number of samples in the class
        weight = np.array([1.0 / np.sum(y == i) for i in classes],
                          dtype=np.float64, order='C')
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
