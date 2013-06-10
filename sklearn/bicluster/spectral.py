"""Spectral co-clustering and spectral biclustering.

Authors : Kemal Eren
License: BSD 3 clause

"""
from .util import Biclusters
from ..base import BaseEstimator, BiclusterMixin
from ..cluster.k_means_ import k_means

import numpy as np
from scipy.sparse.linalg import svds

# TODO: re-use existing functionality in scikit-learn
# TODO: normalization

def _make_diag_root(m):
    return np.diag(1 / np.sqrt(m))


class SpectralCoclustering(BaseEstimator, BiclusterMixin):
    """Spectral Co-Clustering.

    Parameters
    -----------
    n_clusters : integer
        The number of biclusters to find.

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
    `biclusters` : tuple, (rows, columns)
        Results of the clustering. `rows` has shape (n_clusters,
        n_rows), and similarly for `columns`. Available only after
        calling ``fit``.

    References
    ----------

    - Co-clustering documents and words using
      bipartite spectral graph partitioning, 2001
      Dhillon, Inderjit S.
      http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.140.3011

    """
    def __init__(self, n_clusters, maxiter=None, n_init=10, random_state=None):
        self.n_clusters = n_clusters
        self.maxiter = maxiter
        self.n_init = n_init
        self.random_state=random_state

    def fit(self, X):
        """Creates a biclustering for X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        """
        if X.ndim != 2:
            raise Exception('data array must be 2 dimensional')
        diag1 = _make_diag_root(np.sum(X, axis=1))
        diag2 = _make_diag_root(np.sum(X, axis=0))
        an = diag1.dot(X).dot(diag2)
        n_singular_vals = 1 + int(np.ceil(np.log2(self.n_clusters)))
        u, s, vt = svds(an, k=n_singular_vals, maxiter=self.maxiter)
        z = np.vstack((diag1.dot(u[:, 1:]),
                       diag2.dot(vt.T[:, 1:])))
        _, labels, _ = k_means(z, self.n_clusters,
                               random_state=self.random_state,
                               n_init=self.n_init)

        n_rows = X.shape[0]
        row_labels = labels[0:n_rows]
        col_labels = labels[n_rows:]

        rows = np.vstack(row_labels == c for c in range(self.n_clusters))
        cols = np.vstack(col_labels == c for c in range(self.n_clusters))
        self.biclusters = Biclusters(rows, cols)


class SpectralBiclustering(BaseEstimator, BiclusterMixin):
    """Spectral Biclustering. A different way of optimizing the same
    model as in Spectral Co-Clustering.

    The main differences are:
        * The number of row clusters may differ from the number of
          column clusters.
        * The partitioning scheme does not combine row and column
          eigenvectors.

    References
    ----------

    - Spectral biclustering of microarray data:
      coclustering genes and conditions, 2003
      Kluger, Yuval, et al.
      http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.135.1608

    """
    def __init__(self):
        raise NotImplementedError()