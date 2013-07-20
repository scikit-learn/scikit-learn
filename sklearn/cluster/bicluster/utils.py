import numpy as np


def check_array_ndim(X):
    if X.ndim != 2:
        raise ValueError("Argument `X` has the wrong dimensionality."
                         " It must have exactly two dimensions, but"
                         " {} != 2".format(X.ndim))


def get_indices(rows, columns):
    """Convert indicator vectors to lists of indices for bicluster `i`."""
    row_idx = np.nonzero(rows)[0]
    col_idx = np.nonzero(columns)[0]
    return row_idx, col_idx


def get_indicators(rows, columns, shape):
    """Convert indices to indicator vectors"""
    row_ind = np.zeros(shape[0], dtype=np.bool)
    col_ind = np.zeros(shape[1], dtype=np.bool)
    row_ind[rows] = True
    col_ind[columns] = True
    return row_ind, col_ind


def get_shape(rows, columns):
    """Returns shape of bicluster from indicator vectors"""
    indices = get_indices(rows, columns)
    return tuple(len(i) for i in indices)


def get_submatrix(rows, columns, data):
    """Returns the submatrix corresponding to bicluster `i`.

    Works with sparse matrices.

    """
    rows, cols = get_indices(rows, columns)
    return data[rows[:, np.newaxis], cols]
