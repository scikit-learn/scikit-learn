from __future__ import division

import numpy as np

from sklearn.utils.linear_assignment_ import linear_assignment
from sklearn.utils.validation import check_arrays


def _check_rows_and_columns(a, b):
    """Unpacks the row and column arrays and checks their shape."""
    a_rows, a_cols = check_arrays(*a)
    b_rows, b_cols = check_arrays(*b)
    a_rows = a_rows.astype(np.bool)
    a_cols = a_cols.astype(np.bool)
    b_rows = b_rows.astype(np.bool)
    b_cols = b_cols.astype(np.bool)
    return a_rows, a_cols, b_rows, b_cols


def _size(rows, cols):
    return rows.sum() * cols.sum()


def _intersection(a_rows, a_cols, b_rows, b_cols):
    return ((a_rows * b_rows).sum() *
            (a_cols * b_cols).sum())


def _jaccard(a_rows, a_cols, b_rows, b_cols):
    """Jaccard coefficient on the elements of the two biclusters."""
    num = _intersection(a_rows, a_cols, b_rows, b_cols)
    denom = (_size(a_rows, a_cols) + _size(b_rows, b_cols) - num)
    return num / denom


def _dice(a_rows, a_cols, b_rows, b_cols):
    num = 2 * _intersection(a_rows, a_cols, b_rows, b_cols)
    denom = _size(a_rows, a_cols) + _size(b_rows, b_cols)
    return num / denom


def _goodness(a_rows, a_cols, b_rows, b_cols):
    inter = _intersection(a_rows, a_cols, b_rows, b_cols)
    prec = inter / _size(a_rows, a_cols)
    recall = inter / _size(b_rows, b_cols)
    return 0.5 * (prec + recall)


def _pairwise_similarity(a, b, similarity):
    """Computes pairwise similarity matrix.

    result[i, j] is the Jaccard coefficient of a's bicluster i and b's
    bicluster j.

    """
    a_rows, a_cols, b_rows, b_cols = _check_rows_and_columns(a, b)
    n_a = a_rows.shape[0]
    n_b = b_rows.shape[0]
    result = np.zeros((n_a, n_b))
    result = np.array(list(list(similarity(a_rows[i], a_cols[i],
                                           b_rows[j], b_cols[j])
                                for j in range(n_b))
                           for i in range(n_a)))
    return result


def consensus_score(a, b, similarity="jaccard"):
    """The similarity of two sets of biclusters.

    Similarity between individual biclusters is computed. Then the
    best matching between sets is found using the Hungarian algorithm.
    The final score is the sum of similarities divided by the size of
    the larger set.

    Parameters
    ----------
    a : (rows, columns)
        Tuple of row and column indicators for a set of biclusters.

    b : (rows, columns)
        Another set of biclusters like ``a``.

    similarity : string or function, optional, default: "jaccard"
        May be the strings "jaccard", "dice", or "goodness", or
        any function that takes four arguments, each of which is a 1d
        indicator vector: (a_rows, a_columns, b_rows, b_columns).

    References
    ----------

    * Hochreiter, Bodenhofer, et. al., 2010. `FABIA: factor analysis
      for bicluster acquisition
      <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2881408/>`__.

    """
    if isinstance(similarity, str):
        if similarity == "jaccard":
            similarity = _jaccard
        elif similarity == "dice":
            similarity = _dice
        elif similarity == "goodness":
            similarity = _goodness
        else:
            raise ValueError("unknown similarity {}".format(similarity))
    if not callable(similarity):
        raise ValueError("'similarity' argument is not callable")
    matrix = _pairwise_similarity(a, b, similarity)
    indices = linear_assignment(1. - matrix)
    n_a = len(a[0])
    n_b = len(b[0])
    return np.trace(matrix[:, indices[:, 1]]) / max(n_a, n_b)
