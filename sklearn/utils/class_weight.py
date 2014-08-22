# Authors: Andreas Mueller
#          Manoj Kumar
# License: BSD 3 clause

import numpy as np
from ..utils import column_or_1d


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

    classes : ndarray
        Array of the classes occurring in the data, as given by
        ``np.unique(y_org)`` with ``y_org`` the original class labels.

    y : array-like, shape (n_samples,)
        Array of original class labels per sample;

    Returns
    -------
    class_weight_vect : ndarray, shape (n_classes,)
        Array with class_weight_vect[i] the weight for i-th class
    """
    # Import error caused by circular imports.
    from ..preprocessing import LabelEncoder

    if class_weight is None or len(class_weight) == 0:
        # uniform class weights
        weight = np.ones(classes.shape[0], dtype=np.float64, order='C')
    elif class_weight == 'auto':
        if not all(np.in1d(y, classes)):
            raise ValueError("all classes in y should be "
                             "included the classes attribute")

        # Find the weight of each class as present in y.
        le = LabelEncoder()
        class_ind = le.fit_transform(classes)
        y_ind = le.transform(y).ravel()

        # inversely proportional to the number of samples in the class
        recip_freq = 1. / np.bincount(y_ind)
        recip_freq.resize(class_ind.shape[0])
        weight = recip_freq[class_ind] / np.mean(recip_freq)
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
