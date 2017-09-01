import numpy as np
from . import _check_classes, _check_y_classes
from .. preprocessing import LabelEncoder


def _set_classes(classes, y, multioutput=False, return_present_classes=False):
        """Returns the value to be set in classes_ attribute and transformed y

        The classes are selected based on the classes present in y and
        the classes passed in estimator init, if any. The function can also
        return binarized and label encoded y.

        Also performs checks on classes and y.

        Parameters
        ----------
        y : array-like, shape = (n_samples, n_outputs)
            The true labels

        multioutput : boolean
            If True, it will return: (classes, n_classes, y_encoded)
            If False, it will return: (classes, y_encoded)

        return_present_classes : boolean
            This will only work if multioutput is False. Otherwise, this will
            be ignored.
            If True, returns: (classes, present_classes, y_encoded)
            If Falase, returns: (classes, y_encoded)
        """
        if multioutput:
            # multioutput case:
            n_outputs_ = y.shape[1]
            classes_ = []
            n_classes_ = []
            y_encoded = np.zeros(y.shape, dtype=np.int)

            if classes is None:
                for k in range(n_outputs_):
                    classes_k, y_encoded[:, k] = np.unique(
                                                    y[:, k],
                                                    return_inverse=True)
                    classes_.append(classes_k)
                    n_classes_.append(classes_k.shape[0])
            else:
                _check_classes(classes, n_outputs_)
                if n_outputs_ == 1:
                    classes_ = np.atleast_2d(classes)
                else:
                    classes_ = np.asarray(classes)
                # encode y:
                for k in range(n_outputs_):
                    classes_k = np.asarray(classes_[k])
                    # check y has subset of classes:
                    _check_y_classes(np.unique(y[:, k]), classes_k)
                    le = LabelEncoder().fit(classes_k)
                    y_encoded[:, k] = le.transform(y[:, k])
                    n_classes_.append(classes_k.shape[0])

            return classes_, n_classes_, y_encoded

        # non multi-output case:
        # check if classes was passed in init
        present_classes = np.unique(y)
        if classes is None:
            classes_ = present_classes
        else:
            _check_classes(classes)
            classes_ = np.asarray(classes)
            _check_y_classes(present_classes, classes_)

        y_encoded = LabelEncoder().fit(classes_).transform(y)

        if return_present_classes:
            return classes_, present_classes, y_encoded
        else:
            return classes_, y_encoded
