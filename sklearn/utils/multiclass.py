# Author: Arnaud Joly
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
    return (hasattr(y, "shape") and len(y.shape) == 2 and y.shape[1] > 1 and
            y.shape[0] > 1 and np.size(np.unique(y)) <= 2)


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
    return (not isinstance(y[0], np.ndarray) and isinstance(y[0], Sequence) and
            not isinstance(y[0], string_types))


_mixed_error = ValueError(
    'Received a mix of single-label and multi-label data')


def is_multilabel(*ys):
    """Checks if all of `ys` are multilabel, and what types

    If `ys` is a mix of multilabel and non-multilabel data, raises `ValueError`.

    Parameters
    ----------
    *ys : array-likes

    Returns
    -------
    type_name : string or False
        Returns False if none of `ys` is multilabel. Returns 'matrix' if all
        are label indicator matrices, 'sos' if all are sequence-of-sequences,
        and 'mixed' if both occur.

    Examples:
    >>> import numpy as np
    >>> from sklearn.utils.multiclass import is_multilabel
    >>> from nose.tools import assert_raises
    >>> ind = np.array([[0, 1], [1, 0]])
    >>> seq1 = [[0], [0, 2]]
    >>> seq2 = np.array([[0, 2], [1]])
    >>> other = [1, 2, 3]
    >>> is_multilabel(ind)
    'indicator'
    >>> is_multilabel(seq1)
    'sequences'
    >>> is_multilabel(seq2)
    'sequences'
    >>> is_multilabel(seq1, seq2)
    'sequences'
    >>> is_multilabel(ind, ind)
    'indicator'
    >>> is_multilabel(ind, seq1)
    'mixed'
    >>> is_multilabel(ind, seq1, ind)
    'mixed'
    >>> is_multilabel(other)
    False
    >>> is_multilabel(other, other)
    False
    >>> assert_raises(ValueError, is_multilabel, other, seq1)
    >>> assert_raises(ValueError, is_multilabel, other, ind)
    >>> assert_raises(ValueError, is_multilabel, seq1, other)
    >>> assert_raises(ValueError, is_multilabel, ind, other)
    """
    if not ys:
        raise ValueError('Need at least one set of targets to check')

    ret = None
    for y in ys:
        is_ml = is_label_indicator_matrix(y)
        if is_ml:
            if ret is None:
                ret = 'indicator'
            elif ret == 'sequences':
                ret = 'mixed'
            elif ret == False:
                raise _mixed_error
            continue

        is_ml = is_sequence_of_sequences(y)
        if not is_ml:
            if ret:
                raise _mixed_error
            ret = False
        elif ret is None:
            ret = 'sequences'
        elif ret == 'indicator':
            ret = 'mixed'
        elif ret == False:
            raise _mixed_error
    return ret


    if all(is_label_indicator_matrix(y) for y in ys):
        if len(set(y.shape for y in ys)) > 1:
            raise ValueError('Indicator matrices are not all the same shape')
        return 'matrix'
    
    if all(is_multilabel(y) for y in ys):
        if len(set(getattr(y, 'shape', (len(y),))[0] for y in ys)) > 1:
            raise ValueError('Multilabel targets have different sample sizes')
        if any(is_label_indicator_matrix(y) for y in ys):
            return 'mixed'
        return 'sos'

    if any(is_multilabel(y) for y in ys):
        raise ValueError('Received a mix of single-label and multi-label data')
    return False

