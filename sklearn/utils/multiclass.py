# Author: Arnaud Joly, Joel Nothman
#
# License: BSD 3 clause
"""
Multi-class / multi-label utility function
==========================================

"""
from collections import Sequence
from itertools import chain

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
                classes = np.array(sorted(set(chain(*y))))

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
    return (hasattr(y, "shape") and y.ndim == 2 and y.shape[1] > 1 and
            y.shape[0] > 1 and np.size(np.unique(y)) <= 2 and
            issubclass(y.dtype.type, (np.float, np.int)) and
            np.all(y == y.astype(int)))


def is_sequence_of_sequences(y):
    """ Check if ``y`` is in the sequence of sequences format (multilabel).

    Parameters
    ----------
    y : sequence or array.

    Returns
    -------
    out : bool,
        Return ``True``, if ``y`` is a sequence of sequences else ``False``.

    >>> import numpy as np
    >>> from sklearn.utils.multiclass import is_multilabel
    >>> is_sequence_of_sequences([0, 1, 0, 1])
    False
    >>> is_sequence_of_sequences([[1], [0, 2], []])
    True
    >>> is_sequence_of_sequences(np.array([[1], [0, 2], []]))
    True
    >>> is_sequence_of_sequences([(1,), (0, 2), ()])
    True
    >>> is_sequence_of_sequences(np.array([[1, 0], [0, 0]]))
    False
    >>> is_sequence_of_sequences(np.array([[1], [0], [0]]))
    False
    >>> is_sequence_of_sequences(np.array([[1, 0, 0]]))
    False
    """
    # the explicit check for ndarray is for forward compatibility; future
    # versions of Numpy might want to register ndarray as a Sequence
    try:
        return (not isinstance(y[0], np.ndarray) and isinstance(y[0], Sequence)
                and not isinstance(y[0], string_types))
    except IndexError:
        return False


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
    return is_label_indicator_matrix(y) or is_sequence_of_sequences(y)


def type_of_target(y):
    """Determine the type of data indicated by target `y`

    Parameters
    ----------
    y : array-like

    Returns
    -------
    target_type : string
        One of:
        * 'continuous': `y` is an array-like of floats that are not all
          integers, and is a vector in 1 or 2 dimensions.
        * 'continuous-multioutput': `y` is a 2d array of floats that are
          not all integers, and both dimensions are of size > 1.
        * 'binary': `y` contains two values and is not an indicator matrix.
        * 'multiclass': `y` contains more or less than two discrete values,
          is not a sequence of sequences, and is a vector in 1 or 2 dimensions.
        * 'mutliclass-multioutput': `y` is a 2d array that contains more
          than two discrete values, is not a sequence of sequences, and both
          dimensions are of size > 1.
        * 'multilabel-sequences': `y` is a sequence of sequences, a 1d
          array-like of objects that are sequences of labels.
        * 'multilabel-indicator': `y` is a label indicator matrix, an array
          of two dimensions each of size > 1, and at most 2 unique values.

    Examples
    --------
    >>> type_of_target([0.1, 0.6])
    'continuous'
    >>> type_of_target([1, -1, -1, 1])
    'binary'
    >>> type_of_target(['a', 'b', 'a'])
    'binary'
    >>> type_of_target([1, 0, 2])
    'multiclass'
    >>> type_of_target(['a', 'b', 'c'])
    'multiclass'
    >>> type_of_target(['a'])
    'multiclass'
    >>> type_of_target([])
    'multiclass'
    >>> type_of_target([['a', 'b'], ['c'], []])
    'multilabel-sequences'
    >>> type_of_target([[]])
    'multilabel-sequences'
    >>> import numpy as np
    >>> type_of_target(np.array([[0, 1], [1, 1]]))
    'multilabel-indicator'
    """
    if is_sequence_of_sequences(y):
        return 'multilabel-sequences'
    elif is_label_indicator_matrix(y):
        return 'multilabel-indicator'

    y = np.asarray(y)
    if y.ndim > 2 or y.dtype == object:
        return 'unknown'
    if y.ndim == 2 and y.shape[1] > 1 and y.shape[0] > 1:
        suffix = '-multioutput'
    else:
        suffix = ''

    if issubclass(y.dtype.type, np.float) and np.any(y != y.astype(int)):
        return 'continuous' + suffix
    if len(np.unique(y)) == 2:  # XXX: should this be <= 2?
        if suffix:
            # strange case... looks like a label indicator matrix, but did not
            # meet its criteria (e.g. a list of lists, rather than an array)
            return 'unknown'
        return 'binary'
    else:
        return 'multiclass' + suffix
