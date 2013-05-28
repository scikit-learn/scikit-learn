# Author: Arnaud Joly
#
# License: BSD 3 clause
"""
Multi-class / multi-label utility function
==========================================

"""
from collections import Sequence

import numpy as np

from ..externals.six import string_types


def unique_labels(*lists_of_labels):
    """Extract an ordered array of unique labels

    Parameters
    ----------
    lists_of_labels : list of labels,
        The supported "list of labels" are:
            - a list / tuple / numpy array of int
            - a list of lists / tuples of int;
            - a binary indicator matrix (2D numpy array)

    Returns
    -------
    out : numpy array of shape [n_unique_labels]
        An ordered array of unique labels.

    Examples
    --------
    >>> from sklearn.utils.multiclass import unique_labels
    >>> unique_labels([3, 5, 5, 5, 7, 7])
    array([3, 5, 7])
    >>> unique_labels([1, 2, 3, 4], [2, 2, 3, 4])
    array([1, 2, 3, 4])
    >>> unique_labels([1, 2, 10], [5, 11])
    array([ 1,  2,  5, 10, 11])
    >>> unique_labels(np.array([[0.0, 1.0], [1.0, 1.0]]), np.zeros((2, 2)))
    array([0, 1])
    >>> unique_labels([(1, 2), (3,)], [(1, 2), tuple()])
    array([1, 2, 3])

    """
    def _unique_labels(y):
        classes = None
        if is_multilabel(y):
            if is_label_indicator_matrix(y):
                classes = np.arange(y.shape[1])
            else:
                classes = np.array(sorted(set.union(*map(set, y))))

        else:
            classes = np.unique(y)

        return classes

    if not lists_of_labels:
        raise ValueError('No list of labels has been passed.')

    return np.unique(np.hstack(_unique_labels(y) for y in lists_of_labels))


def is_label_indicator_matrix(y):
    """ Check if ``y`` is in the label indicator matrix format (multilabel).

    Parameters
    ----------
    y : numpy array of shape [n_samples] or sequence of sequences
        Target values. In the multilabel case the nested sequences can
        have variable lengths.

    Returns
    -------
    out : bool,
        Return ``True``, if ``y`` is in a label indicator matrix format,
        else ``False``.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.utils.multiclass import is_label_indicator_matrix
    >>> is_label_indicator_matrix([0, 1, 0, 1])
    False
    >>> is_label_indicator_matrix([[1], [0, 2], []])
    False
    >>> is_label_indicator_matrix(np.array([[1, 0], [0, 0]]))
    True
    >>> is_label_indicator_matrix(np.array([[1], [0], [0]]))
    False
    >>> is_label_indicator_matrix(np.array([[1, 0, 0]]))
    False

    """
    return (hasattr(y, "shape") and len(y.shape) == 2 and y.shape[1] > 1 and
            y.shape[0] > 1 and np.size(np.unique(y)) <= 2)


def is_multilabel(y):
    """ Check if ``y`` is in a multilabel format.

    Parameters
    ----------
    y : numpy array of shape [n_samples] or sequence of sequences
        Target values. In the multilabel case the nested sequences can
        have variable lengths.

    Returns
    -------
    out : bool,
        Return ``True``, if ``y`` is in a multilabel format, else ```False``.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.utils.multiclass import is_multilabel
    >>> is_multilabel([0, 1, 0, 1])
    False
    >>> is_multilabel([[1], [0, 2], []])
    True
    >>> is_multilabel(np.array([[1, 0], [0, 0]]))
    True
    >>> is_multilabel(np.array([[1], [0], [0]]))
    False
    >>> is_multilabel(np.array([[1, 0, 0]]))
    False

    """
    # the explicit check for ndarray is for forward compatibility; future
    # versions of Numpy might want to register ndarray as a Sequence
    return (not isinstance(y[0], np.ndarray) and isinstance(y[0], Sequence) and
            not isinstance(y[0], string_types) or is_label_indicator_matrix(y))
