# -*- coding: utf8 -*-
# Author: Alain Pena <alain.pena.r@gmail.com>
#
# License: BSD 3 clause

import numpy as np

from .utils.validation import check_random_state

from .basic_checks import check_is_iterable

__all__ = [
    "convert_probabilities",
    "cycle_permutations",
    "get_most_probable_class",
    "get_random_numbers",
    "iterize",
    "MAX_SEED_INT",
    "next_permutation",
    "next_permutations",
    "probs_to_class",
    "round_nearest",
    "shuffle_array",
]

"""
Maximum value for a random seed.
"""
MAX_SEED_INT = np.iinfo(np.int32).max


def next_permutation(start, compare=cmp):
    """
    Generates the next useful permutation from the start.
    See `next_permutations`.

    Parameters
    ----------
    start : iterable
        The iterable containing the values to rearrange.

    compare : function
        The compare function used to determine the next permutation.

    Returns
    -------
    list of any
        The next permutation.
    """
    for i in next_permutations(start, compare):
        return i
    return None


def next_permutations(start, compare=cmp):
    """
    Generates the next useful permutations from the start
    (start not included).

    Parameters
    ----------
    start : iterable
        The iterable containing the values to rearrange.

    compare : function
        The compare function used to determine the next permutation.

    Returns
    -------
    iterator on a list of any
        The next permutations in order. The iterator will stop when there is
        no greater permutation.

    Reference
    ---------
        http://stackoverflow.com/questions/352203/
        generating-permutations-lazily/353248#353248
    """
    try:
        start_ = list(start[:])
    except TypeError:
        if check_is_iterable(start, string_allowed=False):
            start_ = list(start)
        else:
            # Start is not iterable, so there is no next.
            return

    def find_tail_min_index(lis, thresh):
        a = []
        retv = a
        for x in xrange(1, len(lis) + 1):
            if ((retv is a) and (cmp(lis[- x], thresh) > 0)):
                retv = - x
        return retv

    running = True
    while running:
        for tail_end in xrange(1, len(start_)):
            if (cmp(start_[- tail_end - 1], start_[- tail_end]) < 0):
                # find_tail_min_index WILL return an index
                ind = find_tail_min_index(start_, start_[-tail_end - 1])
                start_[- tail_end - 1], start_[ind] = (start_[ind],
                                                       start_[- tail_end - 1])
                start_[(-tail_end):] = start_[(-tail_end):][::-1]
                yield start_[:]
                break
        else:
            running = False


def cycle_permutations(start, compare=cmp, n_iter=None):
    """
    Generates the next useful permutations from the start
    (start not included).

    Parameters
    ----------
    start : iterable
        The iterable containing the values to rearrange.

    compare : function
        The compare function used to determine the next permutation.

    n_iter : int or None
        The (positive) number of values it should iterate on. If None,
        will forever loop.

    Returns
    -------
    iterator on a list of any
        The next permutations in order. The iterator will loop back
        from the lowest permutation after the highest one has been generated
        until `n_iter` values have been generated.

    Example
    -------
    >>> for p in cycle_permutations("201", n_iter = 10):
    >>>     p
    ['2', '1', '0']
    ['0', '1', '2']
    ['0', '2', '1']
    ['1', '0', '2']
    ['1', '2', '0']
    ['2', '0', '1']
    ['2', '1', '0']
    ['0', '1', '2']
    ['0', '2', '1']
    ['1', '0', '2']
    """
    try:
        start_ = list(start[:])
    except TypeError:
        if check_is_iterable(start, string_allowed=False):
            start_ = list(start)
        else:
            return
    retv = start_
    iter = 10 if (n_iter is None) else n_iter

    def dec_iter(iter_):
        it = iter_
        if n_iter is not None:
            it = it - 1
        return it

    while (iter > 0):
        for result in next_permutations(start_, cmp):
            retv = result
            yield result
            iter = dec_iter(iter)
            if (iter <= 0):
                break
        if (iter > 0):
            start_ = retv[::-1]
            yield start_
            iter = dec_iter(iter)


def round_nearest(array, inplace=False, threshold=0.5):
    """
    Returns an array rounded to 0 or 1 following if the
    each value is [below] or [equal or above] a certain threshold.

    Parameters
    ----------
    array : numpy.array
        The array to round to nearest

    inplace : bool
        Whether or not the "array" argument should be modified.

    threshold : float
        The value of the threshold. Any number < this one will result to 0,
        while others will result to 1.

    Returns
    -------
    numpy.array
        A modified array. If "inplace" is set to True, it returns
        the "array" argument, else it returns a modified copy of it.
    """
    if (inplace):
        retv = array
    else:
        retv = np.empty(array.shape)
    mask = array < threshold
    retv[mask] = 0
    retv[np.logical_not(mask)] = 1
    return retv


def get_most_probable_class(y_proba, axis=1):
    """
    Returns the most probable class from a probabilities array.

    In case of several classes with same probabilities, it selects the first
    one.

    Parameters
    ----------
    y_proba : array_like
        The array of probabilities.

    axis : int
        The axis along which the analysis must be done. By default,
        it selects the most probable class from probabilities displayed
        horizontally.

    Returns
    -------
    numpy.ndarray of ints
        An array containing the indices of the most probable class.

    Raises
    ------
    ValueError
        In case numpy.argmax throws an error.

    Example
    -------
    >>> get_most_probable_class([[0.1, 0.2, 0.5, 0.2],
                                 [0.4, 0.1, 0.1, 0.4],
                                 [0.1, 0.4, 0.2, 0.3]],
                                 axis = 1)
    [2, 0, 1]
    """
    return np.argmax(y_proba, axis=axis)


def convert_probabilities(probas, y_shape):
    """
    Converts the result of a predict_proba classifier to
    the same format of `list(numpy_array)` for multilabel or
    `numpy_array` for single label.

    Where:
        - len(list) = the number of labels

        - numpy_array.shape = (n_samples, 2)

    Parameters
    ---------
    probas : array_like
        The probabilities array to be converted.

        In case this array is a list, the work is supposed to be
        already done, in which case it returns the array unmodified.

    y_shape : array-like
        The shape of the ground-truth.

    Returns
    -------
    (list of) numpy.ndarray
        The list (if multilabel) of labels probability predictions where:
            - len(list) = the number of labels

            - numpy_array.shape = (n_samples, 2)

    TODO
    ----
    May need some tweak to cover multiclass case
    """
    try:
        y_shape[1]
    except:
        return np.array(probas)

    try:
        view = np.reshape(probas, (probas.shape[0], 1))
    except ValueError:  # the array is not unidimensional
        view = np.reshape(probas, (probas.shape[0], probas.shape[1]))
    except AttributeError:  # the array is already a list
        return probas

    retv = [np.empty((view.shape[0], 2), dtype=float)
            for _ in xrange(view.shape[1])]
    for (y, x), value in np.ndenumerate(view):
        retv[x][y, 0] = 1.0 - value
        retv[x][y, 1] = value

    return retv


def shuffle_array(array, inplace=True, random_state=None):
    """
    Shuffles an array following the random state.

    Parameters
    ----------
    array : array_like
        The array to shuffle.

    inplace : bool
        Whether or not the "array" argument should be modified.

    random_state : int or RandomState
        RNG seed.

    Returns
    -------
    array_like (same type as "array" argument)
        A shuffled array. If "inplace" is set to True, it returns
        the "array" argument, else it returns a modified copy of it.
    """
    if (inplace):
        retv = array
    else:
        retv = array.copy()
    check_random_state(seed=random_state).shuffle(retv)
    return retv


def get_random_numbers(shape=None, random_state=None):
    """
    Generates random numbers with a seed.

    Parameters
    ----------
    shape : int or array_like
        The shape of the result.

    random_state : int or RandomState
        RNG seed.

    Returns
    -------
    int or array_like
        Random numbers in [0, sys.maxint].
    """
    r = check_random_state(seed=random_state).randint(MAX_SEED_INT,
                                                      size=shape)
    return r


def iterize(elem):
    """
    Makes an object iterable.
    In the specific case of a string, the object is supposed
    not iterable.

    Parameters
    ----------
    elem : any
        The object to make iterable

    Returns
    -------
    iterable
        Either the object itself if it was already iterable,
        or a tuple containing it.
        In the case of a string, returns a tuple of this string.
    """
    if check_is_iterable(elem, string_allowed=False):
        return elem
    return (elem,)


def probs_to_class(probabilities):
    """
    Converts a list of numpy array of probabilities into their corresponding
    predictions.

    Parameters
    ----------
    probabilities : np.array or list of np.array
        The output from the probabilities predictions from an estimator.

        If it is nor a list nor a numpy array, this method will raise an
        TypeError exception.

    Returns
    -------
    numpy.ndarray
        The predicted classes.

    Raises
    ------
    TypeError
        If the input is neither a list nor a numpy array.
    """
    try:
        probabilities.shape
    except AttributeError:
        if not isinstance(probabilities, list):
            raise TypeError("Incompatible type: {0}.".format(
                type(probabilities)))
    else:
        return get_most_probable_class(probabilities)
    retv = np.asarray([get_most_probable_class(x) for x in probabilities])
    if (retv.size in retv.shape):
        return retv.ravel()
    return np.transpose(retv)
