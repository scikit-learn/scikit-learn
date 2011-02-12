
import numpy as np

__all__ = ['normaliselabels', 'ctransforms']

def normaliselabels(labels):
    '''
    normalised, names = normaliselabels(labels)

    Normalises the labels to be integers from 0 through N-1

    `normalised` is a np.array, while `names` is a list mapping the indices to
    the old names.

    Parameters
    ----------
    labels : any iterable of labels

    Returns
    ------
    normalised : a numpy ndarray of integers 0 .. N-1
    names : list of label names
    '''
    names = sorted(set(labels))
    normalised = map(names.index, labels)
    return np.array(normalised), names
