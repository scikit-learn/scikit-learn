from __future__ import division

import numpy as np

from sklearn.utils.linear_assignment_ import linear_assignment
from sklearn.utils.validation import check_arrays


def replace_nan(f):
    """Transforms nan return value to 0."""
    # ensures that biclusters with no rows/columns have 0 similarity
    # with everything.
    def new_f(*args, **kwargs):
        result = f(*args, **kwargs)
        return 0 if np.isnan(result) else result
    return new_f


def _check_rows_and_columns(a, b):
    """Unpacks the row and column arrays. Checks shape and dtype."""
    a_rows, a_cols = check_arrays(*a)
    b_rows, b_cols = check_arrays(*b)
    a_rows = a_rows.astype(np.bool)
    a_cols = a_cols.astype(np.bool)
    b_rows = b_rows.astype(np.bool)
    b_cols = b_cols.astype(np.bool)
    return a_rows, a_cols, b_rows, b_cols


def _intersection(a_rows, a_cols, b_rows, b_cols):
    """Number of elements in bicluster intersection."""
    return ((a_rows * b_rows).sum() * (a_cols * b_cols).sum())


def _size(rows, cols):
    """Number of elements in bicluster."""
    return rows.sum() * cols.sum()


@replace_nan
def _precision(a_rows, a_cols, b_rows, b_cols):
    """Precision, assuming 'a' is expected bicluster."""
    num = _intersection(a_rows, a_cols, b_rows, b_cols)
    denom = _size(a_rows, a_cols)
    return num / denom


@replace_nan
def _recall(a_rows, a_cols, b_rows, b_cols):
    """Recall, assuming 'a' is expected bicluster."""
    num = _intersection(a_rows, a_cols, b_rows, b_cols)
    denom = _size(b_rows, b_cols)
    return num / denom


@replace_nan
def _precision_corrected(a_rows, a_cols, b_rows, b_cols, dsize):
    """Precision corrected for size bias."""
    p = _precision(a_rows, a_cols, b_rows, b_cols)
    asize = _size(a_rows, a_cols)
    return (dsize * p - asize) / (dsize - asize)


@replace_nan
def _recall_corrected(a_rows, a_cols, b_rows, b_cols, dsize):
    """Recall corrected for size bias."""
    r = _recall(a_rows, a_cols, b_rows, b_cols)
    bsize = _size(b_rows, b_cols)
    return (dsize * r - bsize) / (dsize - bsize)


@replace_nan
def _jaccard(precision, recall):
    """Jaccard coefficient"""
    return (precision * recall) / (precision + recall - precision * recall)


@replace_nan
def _dice(precision, recall):
    """Dice measure. Same as the traditional balanced F-score."""
    return 2 * precision * recall / (precision + recall)


@replace_nan
def _goodness(precision, recall):
    """Goodness measure. Mean of precision and recall."""
    return (precision + recall) / 2.0


def _pairwise_similarity(a, b, similarity):
    """Computes pairwise similarity matrix.

    result[i, j] is the similarity of a's bicluster i and b's
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


def _make_similarity(f, dsize):
    @replace_nan
    def result(a_rows, a_cols, b_rows, b_cols):
        if dsize is None:
            p = _precision(a_rows, a_cols, b_rows, b_cols)
            r = _recall(a_rows, a_cols, b_rows, b_cols)
            return f(p, r)
        else:
            p = _precision_corrected(a_rows, a_cols, b_rows, b_cols, dsize)
            r = _recall_corrected(a_rows, a_cols, b_rows, b_cols, dsize)
            return f(p, r)
    return result


def consensus_score(a, b, similarity="jaccard", correction=None):
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

    correction : int or None, optional, default: None
        If provided, this should be data.size. Used to correct for
        bicluster size bias, as described in Hanczar, et. al. If this
        is used, bicluster similarities may be less than 0, to
        indicate that they are worse than random chance.

    References
    ----------

    * Hochreiter, Bodenhofer, et. al., 2010. `FABIA: factor analysis
      for bicluster acquisition
      <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2881408/>`__.

    * Hanczar, B., & Nadif, M., 2013. `Precision-recall space to
      correct external indices for biclustering.
      <http://jmlr.csail.mit.edu/proceedings/papers/v28/hanczar13.pdf>`__
      136-144).

    """
    if isinstance(similarity, str):
        if similarity == "jaccard":
            similarity = _make_similarity(_jaccard, correction)
        elif similarity == "dice":
            similarity = _make_similarity(_dice, correction)
        elif similarity == "goodness":
            similarity = _make_similarity(_goodness, correction)
        else:
            raise ValueError("unknown similarity {}".format(similarity))
    if not callable(similarity):
        raise ValueError("'similarity' argument is not callable")
    matrix = _pairwise_similarity(a, b, similarity)
    indices = linear_assignment(1.0 - matrix)
    return matrix[indices[:, 0], indices[:, 1]].sum() / max(matrix.shape)
