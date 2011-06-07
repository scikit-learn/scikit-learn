"""K-means clustering for sparse data"""

# Authors: Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Thomas Rueckstiess <ruecksti@in.tum.de>
#          James Bergstra <james.bergstra@umontreal.ca>
#          Jan Schlueter <scikit-learn@jan-schlueter.de>
#          Nelle Varoquaux
#          Peter Prettenhofer
# License: BSD

from itertools import cycle, izip

import numpy as np
from scipy import sparse
from math import floor

from ...base import BaseEstimator
from ...utils import check_random_state

from . import _fast_kmeans


def _compute_cache_dot(centers, X, x_squared_norms):
    """Cache the center nearest to each x in X."""
    return (X * centers.T).argmax(axis=1).astype(np.int32)

def _compute_cache_euclid(centers, X, x_squared_norms):
    """Cache the center nearest to each x in X."""
    XX = np.sum(centers * centers, axis=1)[:, np.newaxis]

    # If added same as compute dot
    #XX = np.ones((centers.shape[0],))[:, np.newaxis]
    
    YY = x_squared_norms[np.newaxis, :]
    
    distances = XX + YY
    distances -= 2 * (centers * X.T)
    distances = np.maximum(distances, 0)
    return distances.argmin(axis=0).astype(np.int32)


compute_cache = _compute_cache_euclid


def _mini_batch_step(X, batch, centers, counts, x_squared_norms):
    """Incremental update of the centers for the Minibatch K-Means algorithm

    Parameters
    ----------

    X: csr_matrix, shape (n_samples, n_features)
        The data matrix in sparse CSR format.

    batch: array, shape (chunk_size)
        The array of sample indices belonging to the current mini-batch.

    centers: array, shape (k, n_features)
        The cluster centers

    counts: array, shape (k, )
         The vector in which we keep track of the numbers of elements in a
         cluster

    x_squared_norms: array, shape (n_samples,)
         The squared norms of each sample in `X`.
    """

    cache = compute_cache(centers, X[batch],
                          x_squared_norms[batch])

    _fast_kmeans._mini_batch_update(X.data, X.indices, X.indptr, batch,
                                    centers, counts, cache)


def compute_squared_norms(X):
    """Compute squared row norms of CSR matrix."""
    ## x_squared_norms = np.zeros((X.shape[0],), dtype=np.float64)
##     for i in xrange(X.shape[0]):
##         x_squared_norms[i] = (X[i].data**2.0).sum()
##     return x_squared_norms
    return np.ones((X.shape[0],), dtype=np.float64)


## def compute_cache(centers, X, x_squared_norms):
##     """Cache the center nearest to each x in X."""
##     return (X * centers.T).argmax(axis=1).astype(np.int32)


def _init_centroids(X, k, init, random_state):
    """Compute initial centroids based on `init`."""
    if init == 'random':
        n_samples = X.shape[0]
        seeds = np.argsort(random_state.rand(n_samples))[:k]
        centers = X[seeds].toarray()
    else:
        raise ValueError("Cannot init centroids with %s." % str(init))

    return centers


class MiniBatchKMeans(BaseEstimator):
    """
    Batch K-Means clustering

    Parameters
    ----------

    k : int, optional, default: 8
        The number of clusters to form as well as the number of
        centroids to generate.

    n_iter : int
        Number of iterations of the k-means algorithm for a
        single run.

    chunk_size: int, optional, default: 300
        Size of the mini batches

    init : {'random', ndarray}
        Method for initialization, defaults to 'random':

        'random': choose k observations (rows) at random from data for
        the initial centroids.

    verbose : int
        Verbose output.

    random_state : RandomState
        The random state

    Methods
    -------

    fit(X):
        Compute K-Means clustering

    partial_fit(X):
        Compute a partial K-Means clustering

    Attributes
    ----------

    cluster_centers_: array, [n_clusters, n_features]
        Coordinates of cluster centers

    labels_:
        Labels of each point

    Reference
    ---------
    http://www.eecs.tufts.edu/~dsculley/papers/fastkmeans.pdf
    """

    def __init__(self, k=8, init='random', n_iter=300, chunk_size=300,
                 verbose=0, random_state=None):
        self.k = k
        self.init = init
        self.n_iter = n_iter
        self.verbose = verbose
        self.random_state = random_state
        self.counts = None
        self.cluster_centers_ = None
        self.chunk_size = chunk_size

    def _check_data(self, X, **params):
        """ Set parameters and check the sample given is larger than k. """
        if not sparse.isspmatrix_csr(X):
            X = sparse.csr_matrix(X)

        if X.shape[0] < self.k:
            raise ValueError("n_samples=%d should be larger than k=%d" % (
                X.shape[0], self.k))
        self._set_params(**params)
        return X

    def fit(self, X, y=None, **params):
        """
        Calculates the centroids on a batch X

        params
        ------
        X: CSR matrix, [n_samples, n_features]
            Data matrix in CSR (scipy.sparse.csr_matrix) format.
        """

        self.random_state = check_random_state(self.random_state)

        X = self._check_data(X, **params)

        idx = np.arange(X.shape[0], dtype=np.int32)
        self.random_state.shuffle(idx)

        self.cluster_centers_ = _init_centroids(
                X, self.k, self.init, self.random_state)
        self.counts = np.zeros(self.k, dtype=np.int32)

        x_squared_norms = compute_squared_norms(X)

        try:
            n_batches = floor(float(X.shape[0]) / self.chunk_size)
            batches = np.array_split(idx, n_batches)
            n_batches = len(batches)
        except ValueError:
            batches = [idx]
            n_batches = 1

        batches_cycle = cycle(b for b in batches)
        n_iterations = xrange(self.n_iter * n_batches)
        print "xnorms:\n", np.unique(x_squared_norms)
        for i, batch in izip(n_iterations, batches_cycle):
            # old_centers = self.cluster_centers_.copy()
            _mini_batch_step(X, batch, self.cluster_centers_, self.counts,
                             x_squared_norms)

        self.labels_ = compute_cache(self.cluster_centers_, X,
                                     x_squared_norms)

        return self

    def partial_fit(self, X, y=None, **params):
        """Update k means estimate on a single mini-batch X"""

        self.random_state = check_random_state(self.random_state)

        X = self._check_data(X, **params)

        if len(X) == 0:
            return self

        if self.counts is None:
            # this is the first call partial_fit on this object:
            # initialize the cluster centers
            self.cluster_centers_ = _init_centroids(
                X, self.k, self.init, self.random_state)
            self.counts = np.zeros(self.k, dtype=np.int32)

        x_squared_norms = compute_squared_norms(X)

        batch = np.arange(X.shape[0])

        _mini_batch_step(X, batch, self.cluster_centers_, self.counts,
                         x_squared_norms)

        self.labels_ = compute_cache(self.cluster_centers_, X,
                                     x_squared_norms)

        return self
