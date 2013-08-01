"""Implements the Cheng and Church biclustering algorithm.

Authors : Kemal Eren
License: BSD 3 clause

"""
from abc import ABCMeta

import numpy as np

from sklearn.base import BaseEstimator, BiclusterMixin
from sklearn.externals import six

from sklearn.utils.validation import check_arrays
from sklearn.utils.validation import check_random_state

from .utils import check_array_ndim
from ._squared_residue import compute_msr


class EmptyBiclusterException(Exception):
    pass


class IncrementalMSR(object):
    """Incrementally calculates MSR during node deletion."""
    def __init__(self, rows, cols, arr, tol=1e-5):
        assert rows.dtype == np.bool
        assert cols.dtype == np.bool

        self.arr = arr
        self.rows = rows
        self.cols = cols
        self.tol = tol

        self._sum = arr.sum()
        self._row_sum = arr.sum(axis=1)
        self._col_sum = arr.sum(axis=0)
        self._row_idxs = None
        self._col_idxs = None

        self._reset()

    def _reset(self):
        self._msr = None
        self._row_msr = None
        self._col_msr = None

    @property
    def row_idxs(self):
        if self._row_idxs is None:
            self._row_idxs = np.nonzero(self.rows)[0]
        return self._row_idxs

    @property
    def col_idxs(self):
        if self._col_idxs is None:
            self._col_idxs = np.nonzero(self.cols)[0]
        return self._col_idxs

    def remove_row(self, row):
        if not self.rows[row]:
            raise ValueError('cannot remove row {}; it is not in the'
                             ' bicluster'.format(row))
        if len(self.row_idxs) <= 1:
            raise EmptyBiclusterException()
        self._reset()
        vec = self.arr[row, self.col_idxs].ravel()
        self._sum -= vec.sum()
        self._col_sum -= vec

        idx = np.searchsorted(self.row_idxs, row)
        self._row_sum = np.delete(self._row_sum, idx)
        self.rows[row] = False
        self._row_idxs = None

    def remove_col(self, col):
        if not self.cols[col]:
            raise ValueError('cannot remove col {}; it is not in the'
                             ' bicluster'.format(col))
        if len(self.col_idxs) <= 1:
            raise EmptyBiclusterException()
        self._reset()
        vec = self.arr[self.row_idxs, col].ravel()
        self._sum -= vec.sum()
        self._row_sum -= vec

        idx = np.searchsorted(self.col_idxs, col)
        self._col_sum = np.delete(self._col_sum, idx)
        self.cols[col] = False
        self._col_idxs = None

    def remove_rows(self, rows):
        for r in rows:
            self.remove_row(r)

    def remove_cols(self, cols):
        for c in cols:
            self.remove_col(c)

    def _compute(self):
        n_rows = len(self.row_idxs)
        n_cols = len(self.col_idxs)

        row_mean = self._row_sum / n_cols
        col_mean = self._col_sum / n_rows
        mean = self._sum / (n_rows * n_cols)

        self._msr, self._row_msr, self._col_msr = \
            compute_msr(self.row_idxs, self.col_idxs, row_mean,
                        col_mean, mean, self.arr)
        self._msr = 0 if self._msr < self.tol else self._msr
        self._row_msr[self._row_msr < self.tol] = 0
        self.col_msr[self._col_msr < self.tol] = 0

    @property
    def msr(self):
        if self._msr is None:
            self._compute()
        return self._msr

    @property
    def row_msr(self):
        if self._row_msr is None:
            self._compute()
        return self._row_msr

    @property
    def col_msr(self):
        if self._col_msr is None:
            self._compute()
        return self._col_msr


class ChengChurch(six.with_metaclass(ABCMeta, BaseEstimator,
                                     BiclusterMixin)):
    """Algorithm that finds biclusters with small mean squared residue (MSR).

    The residue of an array ``X`` is calculated as ``X -
    X.mean(axis=1, keepdims=True) - X.mean(axis=0) + X.mean()``. It
    measures an element's coherence with the overall array, other
    elements in the same row, and other elements in the same column.
    To get the mean squared residue, the residues are squaredd and their
    mean is calculated.

    ChengChurch tries to maximize bicluser size with the constraint
    that its mean squared residue cannot exceed ``max_msr``.

    Parameters
    -----------
    n_clusters : integer, optional, default: 100
        The number of biclusters to find.

    max_msr : float, default: 1.0
        Maximum mean squared residue of a bicluster. Equivalent to
        'delta` in original paper.

    deletion_threshold : float, optional, default: 1.5
        Multiplier for multiple node deletion. Equivalent to `alpha`
        in original paper.

    row_deletion_cutoff : integer, optional, default: 100
        Number of rows at which to switch to single node deletion.

    column_deletion_cutoff : integer, optional, default: 100
        Number of columns at which to switch to single node deletion.

    inverse_rows : bool, optional, default: False
        If the inverse of a row has a low MSR, add it to the bicluster
        during node addition.

    inverse_columns : bool, optional, default: False
        If the inverse of a column has a low MSR, add it to the bicluster
        during node addition.

    random_state : int seed, RandomState instance, or None (default)
        A pseudo random number generator used by the K-Means
        initialization.

    Attributes
    ----------
    `rows_` : array-like, shape (n_row_clusters, n_rows)
        Results of the clustering. `rows[i, r]` is True if cluster `i`
        contains row `r`. Available only after calling ``fit()``.

    `columns_` : array-like, shape (n_column_clusters, n_columns)
        Results of the clustering, like `rows`.


    References
    ----------

    - Cheng, Y., & Church, G. M. (2000). `Biclustering of
      expression data
      <ftp://samba.ad.sdsc.edu/pub/sdsc/biology/ISMB00/157.pdf>`__.

    """
    def __init__(self, n_clusters=100, max_msr=1.0,
                 deletion_threshold=1.5, row_deletion_cutoff=100,
                 column_deletion_cutoff=100, inverse_rows=False,
                 inverse_columns=False, random_state=None):
        self.n_clusters = n_clusters
        self.max_msr = max_msr
        self.deletion_threshold = deletion_threshold
        self.row_deletion_cutoff = row_deletion_cutoff
        self.column_deletion_cutoff = column_deletion_cutoff
        self.inverse_rows = inverse_rows
        self.inverse_columns = inverse_columns
        self.random_state = random_state

    def _check_parameters(self):
        if self.n_clusters < 1:
            raise ValueError("'n_clusters' must be > 0, but its value"
                             " is {}".format(self.n_clusters))
        if self.max_msr < 0:
            raise ValueError("'max_msr' must be > 0.0, but its value"
                             " is {}".format(self.max_msr))
        if self.deletion_threshold < 1:
            raise ValueError("'deletion_threshold' must be >= 1.0, but its"
                             " value is {}".format(self.deletion_threshold))
        if self.row_deletion_cutoff < 1:
            raise ValueError("'row_deletion_cutoff' must be >= 1, but its"
                             " value is {}".format(self.row_deletion_cutoff))
        if self.column_deletion_cutoff < 1:
            raise ValueError("'column_deletion_cutoff' must be >= 1, but its"
                             " value is {}".format(
                                 self.column_deletion_cutoff))

    def _msr(self, rows, cols, X):
        rows = rows.nonzero()[0][:, np.newaxis]
        cols = cols.nonzero()[0]
        if not rows.size or not cols.size:
            raise EmptyBiclusterException()
        sub = X[rows, cols]
        residue = (sub - sub.mean(axis=1, keepdims=True) -
                   sub.mean(axis=0) + sub.mean())
        return np.power(residue, 2).mean()

    def _row_msr(self, rows, cols, X, inverse=False):
        if not np.count_nonzero(rows) or not np.count_nonzero(cols):
            raise EmptyBiclusterException()
        row_mean = X[:, cols].mean(axis=1, keepdims=True)
        col_mean = X[rows][:, cols].mean(axis=0)
        if inverse:
            arr = (-X[:, cols] + row_mean - col_mean +
                   X[rows][:, cols].mean())
        else:
            arr = (X[:, cols] - row_mean - col_mean +
                   X[rows][:, cols].mean())
        return np.power(arr, 2).mean(axis=1)

    def _col_msr(self, rows, cols, X, inverse=False):
        if not rows.size or not cols.size:
            raise EmptyBiclusterException()
        row_mean = X[rows][:, cols].mean(axis=1, keepdims=True)
        col_mean = X[rows, :].mean(axis=0)
        if inverse:
            arr = (-X[rows, :] - row_mean + col_mean +
                   X[rows][:, cols].mean())
        else:
            arr = X[rows, :] - row_mean - col_mean + X[rows][:, cols].mean()
        return np.power(arr, 2).mean(axis=0)

    def _single_node_deletion(self, rows, cols, X):
        inc = IncrementalMSR(rows, cols, X)
        while inc.msr > self.max_msr:
            row_idx = np.argmax(inc.row_msr)
            col_idx = np.argmax(inc.col_msr)
            if inc.row_msr[row_idx] > inc.col_msr[col_idx]:
                inc.remove_row(inc.row_idxs[row_idx])
            else:
                inc.remove_col(inc.col_idxs[col_idx])
        return inc.rows, inc.cols

    def _multiple_node_deletion(self, rows, cols, X):
        inc = IncrementalMSR(rows, cols, X)
        while inc.msr > self.max_msr:
            n_rows = np.count_nonzero(rows)
            n_cols = np.count_nonzero(cols)
            if n_rows >= self.row_deletion_cutoff:
                to_remove = inc.row_msr > (self.deletion_threshold * inc.msr)
                inc.remove_rows(inc.row_idxs[to_remove])

            if n_cols >= self.column_deletion_cutoff:
                to_remove = inc.col_msr > (self.deletion_threshold *
                                           inc.msr)
                inc.remove_cols(inc.col_idxs[to_remove])

            if (n_rows == np.count_nonzero(rows)) and \
               (n_cols == np.count_nonzero(cols)):
                break
        return inc.rows, inc.cols

    def _node_addition(self, rows, cols, X):
        while True:
            n_rows = np.count_nonzero(rows)
            n_cols = np.count_nonzero(cols)

            old_cols = cols.copy()  # save for inverse
            msr = self._msr(rows, cols, X)
            col_msr = self._col_msr(rows, cols, X)
            to_add = col_msr < msr
            cols = cols + to_add

            if self.inverse_columns:
                col_msr = self._col_msr(rows, old_cols, X,
                                        inverse=True)
                to_add = col_msr < msr
                cols = cols + to_add

            old_rows = rows.copy()  # save for inverse
            msr = self._msr(rows, cols, X)
            row_msr = self._row_msr(rows, cols, X)
            to_add = row_msr < msr
            rows = rows + to_add

            if self.inverse_rows:
                row_msr = self._row_msr(old_rows, cols, X,
                                        inverse=True)
                to_add = row_msr < msr
                rows = rows + to_add

            if (n_rows == np.count_nonzero(rows)) and \
               (n_cols == np.count_nonzero(cols)):
                break
        return rows, cols

    def _mask(self, X, rows, cols, generator, minval, maxval):
        shape = np.count_nonzero(rows), np.count_nonzero(cols)
        mask_vals = generator.uniform(minval, maxval, shape)
        r = rows.nonzero()[0][:, np.newaxis]
        c = cols.nonzero()[0]
        X[r, c] = mask_vals

    def fit(self, X):
        X = X.copy()  # need to modify it in-place
        self._check_parameters()
        X, = check_arrays(X, dtype=np.float64)
        check_array_ndim(X)
        minval, maxval = X.min(), X.max()
        n_rows, n_cols = X.shape

        generator = check_random_state(self.random_state)
        result_rows = []
        result_cols = []

        for i in range(self.n_clusters):
            rows = np.ones(n_rows, dtype=np.bool)
            cols = np.ones(n_cols, dtype=np.bool)
            rows, cols = self._multiple_node_deletion(rows, cols, X)
            rows, cols = self._single_node_deletion(rows, cols, X)
            rows, cols = self._node_addition(rows, cols, X)
            self._mask(X, rows, cols, generator, minval, maxval)
            result_rows.append(rows)
            result_cols.append(cols)

        self.rows_ = np.vstack(result_rows)
        self.columns_ = np.vstack(result_cols)
