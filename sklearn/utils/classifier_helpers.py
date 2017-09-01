import numpy as np
from . import _check_classes, _check_y_classes
from .. preprocessing import label_binarize


def _set_classes(classes, y, label_binarize_y=False):
        """Returns the value to be set in classes_ attribute and transformed y

        The classes are selected based on the classes present in y and
        the classes passed in estimator init, if any. The function can also
        return binarized and label encoded y.

        Also performs checks on classes and y.

        Parameters
        ----------
        y : array-like, shape = (n_samples, n_outputs)
            The true labels

        Returns
        -------
        classes : array-like, shape = (n_classes)
            A sorted array of unique classes to be set as class_ attriute

        y_transformed : array-like or None
            The transformed version of y based on function inputs or None if
            no transformations passed as True.
        """
        y_output = None

        # check if classes was passed in init
        if classes is None:
            classes_ = np.unique(y)
        else:
            _check_classes(classes)
            classes_ = np.asarray(classes)
            _check_y_classes(np.unique(y), classes_)

        if label_binarize_y:
            y_output = label_binarize(y, classes=classes_)

        return classes_, y_output
