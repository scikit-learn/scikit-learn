"""Implements spectral biclustering algorithms.

Authors : Kemal Eren
License: BSD 3 clause

"""
from ..base import BaseEstimator, BiclusterMixin
from ..cluster.k_means_ import k_means
from sklearn.utils.extmath import randomized_svd
from sklearn.utils.arpack import svds
from sklearn.utils.validation import assert_all_finite, check_arrays

import numpy as np


def _make_nonnegative(X, min_value=0):
    """Ensure `X.min()` >= `min_value`."""
    min_ = X.min()
    if min_ < min_value:
        X = X + (min_value - min_)
    return X


def _scale_preprocess(X):
    """Normalize `X` by scaling rows and columns independently.

    Returns the normalized matrix and the row and column scaling
    factors.

    """
    X = _make_nonnegative(X)
    row_diag = 1.0 / np.sqrt(np.sum(X, axis=1))
    col_diag = 1.0 / np.sqrt(np.sum(X, axis=0))
    an = row_diag[:, np.newaxis] * X * col_diag
    return an, row_diag, col_diag


def _bistochastic_preprocess(X, maxiter=1000, tol=1e-5):
    """Normalize rows and columns of `X` simultaneously so that all
    rows sum to one constant and all columns sum to a different
    constant.

    """
    # According to paper, this can also be done more efficiently with
    # deviation reduction and balancing algorithms.
    X = _make_nonnegative(X)
    X_scaled = X
    dist = None
    for _ in range(maxiter):
        X_new, _, _ = _scale_preprocess(X_scaled)
        dist = np.linalg.norm(X_scaled - X_new)
        X_scaled = X_new
        if dist is not None and dist < tol:
            break
    return X_scaled


def _log_preprocess(X):
    """Normalize `X` according to Kluger's log-interactions scheme."""
    X = _make_nonnegative(X, min_value=1)
    L = np.log(X)
    row_avg = np.mean(L, axis=1)[:, np.newaxis]
    col_avg = np.mean(L, axis=0)
    avg = np.mean(L)
    return L - row_avg - col_avg + avg


def _svd(array, n_components, method, kwargs, random_state):
    """Returns first `n_components` left and right singular
    vectors u and v.

    """
    if method == 'randomized':
        u, _, vt = randomized_svd(array, n_components,
                                 random_state=random_state, **kwargs)

    elif method == 'arpack':
        u, _, vt = svds(array, k=n_components, **kwargs)

    assert_all_finite(u)
    assert_all_finite(vt)
    return u, vt.T


def _fit_best_piecewise(vectors, k, n_clusters, random_state,
                        kmeans_kwargs, return_piecewise=False):
    """Find the `k` vectors that are best approximated by piecewise
    constant vectors.

    The piecewise vectors are found by k-means; the best is chosen
    according to Euclidean distance.

    """
    def make_piecewise(v):
        centroid, labels, _ = k_means(v.reshape(-1, 1), n_clusters,
                                      random_state=random_state,
                                      **kmeans_kwargs)
        return centroid[labels].reshape(-1)
    piecewise_vectors = np.apply_along_axis(make_piecewise, 1,
                                            vectors)
    dists = np.apply_along_axis(np.linalg.norm, 1,
                                vectors - piecewise_vectors)
    if return_piecewise:
        result = piecewise_vectors[np.argsort(dists)[:k]]
    else:
        result = vectors[np.argsort(dists)[:k]]
    return result


def _project_and_cluster(data, vectors, n_clusters, random_state,
                         kmeans_kwargs):
    """Project `data` to `vectors` and cluster the result."""
    projected = np.dot(data, vectors)
    _, labels, _ = k_means(projected, n_clusters,
                           random_state=random_state,
                           **kmeans_kwargs)
    return labels


class SpectralBiclustering(BaseEstimator, BiclusterMixin):
    """Spectral biclustering.

    For equivalence with the Spectral Co-Clustering algorithm
    (Dhillon, 2001), use method='dhillon'.

    For the Spectral Biclustering algorithm (Kluger, 2003), use
    one of 'scale', 'bistochastic', or 'log'.

    Parameters
    ----------
    n_clusters : integer or tuple (rows, columns)
        The number of biclusters to find. If method is not 'dhillon',
        the number of row and column clusters may be different.

    method : string
        Method of normalizing and converting singular vectors into
        biclusters. May be one of 'dhillon', 'scale', 'bistochastic',
        or 'log'.

    n_components : integer
        Number of singular vectors to check. Not used if
        `self.method` is 'dhillon'.

    n_best : integer
        Number of best singular vectors to which to project the data
        for clustering. Not used if `self.method` is `dhillon`.

    svd_method : string, optional, default: 'randomized'
        Selects the algorithm for finding singular vectors. May be
        'randomized' or 'arpack'. If 'randomized', uses
        `sklearn.utils.extmath.randomized_svd`, which is faster. If
        'arpack', uses `sklearn.utils.arpack.svds`, which is slower,
        but more accurate.

    svd_kwargs : dictionary, optional, default: {}
        Keyword arguments to pass to the svd function.

    kmeans_kwargs : dictionary, optinal, default: {}
        Keyword arguments to pass to sklearn.cluster.k_means_.

    random_state : int seed, RandomState instance, or None (default)
        A pseudo random number generator used by the K-Means
        initialization.

    Attributes
    ----------
    `rows_` : array-like, shape (n_row_clusters, n_rows)
        Results of the clustering. `rows[i, r]` is True if cluster `i`
        contains row `r`. Available only after calling ``fit``.

    `columns_` : array-like, shape (n_column_clusters, n_columns)
        Results of the clustering, like `rows`.

    `row_labels_` : array-like, shape (n_rows,)
        If `method` is `dhillon`, gives the bicluster label of each
        row. Otherwise, the resulting biclusters have a checkerboard
        structure.

    `column_labels_` : array-like, shape (n_cols,)
        See `row_labels_`.


    References
    ----------

    - Co-clustering documents and words using
      bipartite spectral graph partitioning, 2001
      Dhillon, Inderjit S.
      http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.140.3011

    - Spectral biclustering of microarray data:
      coclustering genes and conditions, 2003
      Kluger, Yuval, et al.
      http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.135.1608

    """
    def __init__(self, n_clusters=3, method='bistochastic',
                 n_components=6, n_best=3,
                 svd_method='randomized', svd_kwargs={}, kmeans_kwargs={},
                 random_state=None):
        legal_methods = ('dhillon', 'bistochastic', 'scale', 'log')
        if method not in legal_methods:
            raise ValueError("Unknown method: '{}'. `method` must be"
                             " one of {}.".format(method, legal_methods))
        legal_svd_methods = ('randomized', 'arpack')
        if svd_method not in legal_svd_methods:
            raise ValueError("Unknown SVD method: '{}'. `svd_method` must be"
                             " one of {}.".format(svd_method,
                                                  legal_svd_methods))

        if hasattr(n_clusters, '__len__'):
            if len(n_clusters) != 2:
                raise ValueError("The parameter `n_clusters` has the"
                                 " wrong length. It should either be an "
                                 " integer or have two entries:"
                                 " `(n_row_clusters, n_column_clusters)`")
            if method == 'dhillon' and n_clusters[0] != n_clusters[1]:
                raise ValueError("Different number of row and column"
                                 " clusters not supported when"
                                 " method=='dhillon'.")
        if n_best > n_components:
            raise ValueError("`n_best` cannot be larger than"
                             " `n_components`.")
        self.n_clusters = n_clusters
        self.method = method
        self.n_components = n_components
        self.n_best = n_best
        self.svd_method = svd_method
        self.svd_kwargs = svd_kwargs
        self.kmeans_kwargs = kmeans_kwargs
        self.random_state = random_state

    def _dhillon(self, X):
        if hasattr(self.n_clusters, '__len__'):
            n_clusters = self.n_clusters[0]
        else:
            n_clusters = self.n_clusters
        normalized_data, row_diag, col_diag = _scale_preprocess(X)
        n_singular_vals = 1 + int(np.ceil(np.log2(n_clusters)))
        u, v = _svd(normalized_data, n_singular_vals, self.svd_method,
                    self.svd_kwargs, self.random_state)
        z = np.vstack((row_diag[:, np.newaxis] * u[:, 1:],
                       col_diag[:, np.newaxis] * v[:, 1:]))
        _, labels, _ = k_means(z, n_clusters,
                               random_state=self.random_state,
                               **self.kmeans_kwargs)

        n_rows = X.shape[0]
        self.row_labels_ = labels[0:n_rows]
        self.column_labels_ = labels[n_rows:]

        self.rows_ = np.vstack(self.row_labels_ == c
                               for c in range(n_clusters))
        self.columns_ = np.vstack(self.column_labels_ == c
                                  for c in range(n_clusters))

    def _kluger(self, X):
        n_sv = self.n_components
        if self.method == 'bistochastic':
            normalized_data = _bistochastic_preprocess(X)
            n_sv += 1
        elif self.method == 'scale':
            normalized_data, _, _ = _scale_preprocess(X)
            n_sv += 1
        elif self.method == 'log':
            normalized_data = _log_preprocess(X)
        u, v = _svd(normalized_data, n_sv, self.svd_method,
                    self.svd_kwargs, self.random_state)
        ut = u.T
        vt = v.T
        if self.method != 'log':
            ut = ut[1:]
            vt = vt[1:]

        if hasattr(self.n_clusters, '__len__'):
            n_row_clusters, n_col_clusters = self.n_clusters
        else:
            n_row_clusters = n_col_clusters = self.n_clusters

        best_ut = _fit_best_piecewise(ut, self.n_best,
                                      n_row_clusters,
                                      self.random_state,
                                      self.kmeans_kwargs)

        best_vt = _fit_best_piecewise(vt, self.n_best,
                                      n_col_clusters,
                                      self.random_state,
                                      self.kmeans_kwargs)

        self.row_labels_ = _project_and_cluster(X, best_vt.T,
                                                n_row_clusters,
                                                self.random_state,
                                                self.kmeans_kwargs)

        self.column_labels_ = _project_and_cluster(X.T, best_ut.T,
                                                   n_col_clusters,
                                                   self.random_state,
                                                   self.kmeans_kwargs)

        self.rows_ = np.vstack(self.row_labels_ == label
                               for label in range(n_row_clusters)
                               for _ in range(n_col_clusters))
        self.columns_ = np.vstack(self.column_labels_ == label
                                  for _ in range(n_row_clusters)
                                  for label in range(n_col_clusters))

    def fit(self, X):
        """Creates a biclustering for X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        """
        X = check_arrays(X, sparse_format='csr', dtype=np.float64)[0]
        if X.ndim != 2:
            raise ValueError("Argument `X` has the wrong dimensionality."
                             " It must have exactly two dimensions.")
        if self.method == 'dhillon':
            self._dhillon(X)
        else:
            self._kluger(X)
