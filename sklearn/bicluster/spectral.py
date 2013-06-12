"""Implements spectral biclustering algorithms.

Authors : Kemal Eren
License: BSD 3 clause

"""
from ..base import BaseEstimator, BiclusterMixin
from ..cluster.k_means_ import k_means

import numpy as np
from scipy.sparse.linalg import svds

# TODO: re-use existing functionality in scikit-learn
# TODO: within-cluster rankings
# TODO: can we use Dhillon's preprocessing but Kluger's postprocessing?

def scaling_preprocess(X):
    raise NotImplementedError()

def bistochastic_preprocess(X):
    raise NotImplementedError()

def log_preprocess(X):
    raise NotImplementedError()


class SpectralBiclustering(BaseEstimator, BiclusterMixin):
    """Spectral biclustering.

    For equivalence with the Spectral Co-Clustering algorithm
    (Dhillon, 2001), use method='dhillon'.

    For the Spectral Biclustering algorithm (Kluger, 2003), use
    one of 'scaling', 'bistochastic', or 'log'.

    Parameters
    -----------
    n_clusters : integer
        The number of biclusters to find.

    method : string
        Method of preparing data matrix for SVD and converting
        singular vectors into biclusters. May be one of 'dhillon',
        'scaling', 'bistochastic', or 'log'.

    n_init : int, optional, default: 10
        Number of time the k-means algorithm will be run with different
        centroid seeds. The final results will be the best output of
        n_init consecutive runs in terms of inertia.

    random_state : int seed, RandomState instance, or None (default)
        A pseudo random number generator used for the initialization
        of the lobpcg eigen vectors decomposition when eigen_solver == 'amg'
        and by the K-Means initialization.

    Attributes
    ----------
    `rows` : array-like, shape (n_clusters, n_rows)
        Results of the clustering. `rows[i, r]` is True if cluster `i`
        contains row `r`. Available only after calling ``fit``.

    `columns` : array-like, shape (n_clusters, n_columns)
        Results of the clustering, like `rows`.


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
    def __init__(self, n_clusters=8, method='dhillon', maxiter=None,
                 n_init=10, random_state=None):
        if method not in ('dhillon', 'bistochastic', 'scaling', 'log'):
            raise Exception('unknown method: {}'.format(method))
        self.n_clusters = n_clusters
        self.method = method
        self.maxiter = maxiter
        self.n_init = n_init
        self.random_state=random_state

    def _dhillon(self, X):
        diag1 = 1.0 / np.sqrt(np.sum(X, axis=1))
        diag2 = 1.0 / np.sqrt(np.sum(X, axis=0))
        an = diag1[:, np.newaxis] * X * diag2
        n_singular_vals = 1 + int(np.ceil(np.log2(self.n_clusters)))
        u, s, vt = svds(an, k=n_singular_vals, maxiter=self.maxiter)

        z = np.vstack((diag1[:, np.newaxis] * u[:, 1:],
                       diag2[:, np.newaxis] * vt.T[:, 1:]))
        _, labels, _ = k_means(z, self.n_clusters,
                               random_state=self.random_state,
                               n_init=self.n_init)

        n_rows = X.shape[0]
        row_labels = labels[0:n_rows]
        col_labels = labels[n_rows:]

        self.rows = np.vstack(row_labels == c for c in range(self.n_clusters))
        self.columns = np.vstack(col_labels == c for c in range(self.n_clusters))

    def _kluger(self, X):
        raise NotImplementedError()

    def fit(self, X):
        """Creates a biclustering for X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        """
        if X.ndim != 2:
            raise Exception('data array must be 2 dimensional')

        if self.method == 'dhillon':
            self._dhillon(X)
        else:
            self._kluger(X)
