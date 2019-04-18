"""Fuzzy C-Means clustering"""

# Authors: Junyi Li (lijy263@mail2.sysu.edu.cn)
# License: BSD 3 clause

import warnings

import numpy as np
import scipy.sparse as sp

from ..base import BaseEstimator, ClusterMixin, TransformerMixin
from ..metrics.pairwise import euclidean_distances
from ..metrics.pairwise import pairwise_distances_argmin_min
from ..utils.extmath import row_norms, squared_norm, stable_cumsum
from ..utils.sparsefuncs_fast import assign_rows_csr
from ..utils.sparsefuncs import mean_variance_axis
from ..utils.validation import _num_samples
from ..utils import check_array
from ..utils import gen_batches
from ..utils import check_random_state
from ..utils.validation import check_is_fitted
from ..utils.validation import FLOAT_DTYPES
from ..utils._joblib import Parallel
from ..utils._joblib import delayed
from ..utils._joblib import effective_n_jobs
from ..exceptions import ConvergenceWarning
from . import _k_means
from ._k_means_elkan import k_means_elkan


def fcm(X, c_clusters, m=2, eps=10, random_state=None, max_iter=300):
    
    if m <= 1:
        raise ValueError("Invalid number of m."
                         " m=%d must be bigger than 1." % m)
    if eps <= 0:
        raise ValueError("Invalid number of eps."
                         " eps=%d must be bigger than 0." % eps)
    random_state = check_random_state(random_state)

    if max_iter <= 0:
        raise ValueError('Number of iterations should be a positive number,'
                         ' got %d instead' % max_iter)
    
    membership_mat = np.random.random((len(X), c_clusters))
    membership_mat = np.divide(membership_mat, np.sum(membership_mat, axis=1)[:, np.newaxis])

    for iter_time in range(max_iter):
        working_membership_mat = membership_mat ** m
        Centroids = np.divide(np.dot(working_membership_mat.T, X), np.sum(working_membership_mat.T, axis=1)[:, np.newaxis])
        
        n_c_distance_mat = np.zeros((len(X), c_clusters))
        for i, x in enumerate(X):
            for j, c in enumerate(Centroids):
                n_c_distance_mat[i][j] = np.linalg.norm(x-c, 2)
        
        new_membership_mat = np.zeros((len(X), c_clusters))
        
        for i, x in enumerate(X):
            for j, c in enumerate(Centroids):
                new_membership_mat[i][j] = 1. / np.sum((n_c_distance_mat[i][j] / n_c_distance_mat[i]) ** (2 / (m - 1)))
        if np.sum(abs(new_membership_mat - membership_mat)) < eps:
            break
        membership_mat =  new_membership_mat
    
    return Centroids, np.argmax(new_membership_mat, axis=1), new_membership_mat

class FCM(BaseEstimator, ClusterMixin, TransformerMixin):
    """K-Means clustering

    Read more in the :ref:`User Guide <k_means>`.

    Parameters
    ----------

    n_clusters : int, optional, default: 8
        The number of clusters to form as well as the number of
        centroids to generate.

    init : {'random'}
        Method for initialization, defaults to 'random':

        'random': choose k observations (rows) at random from data for
        the initial centroids.

        If an ndarray is passed, it should be of shape (n_clusters, n_features)
        and gives the initial centers.

    max_iter : int, default: 300
        Maximum number of iterations of the k-means algorithm for a
        single run.

    random_state : int, RandomState instance or None (default)
        Determines random number generation for centroid initialization. Use
        an int to make the randomness deterministic.
        See :term:`Glossary <random_state>`.

    Attributes
    ----------
    cluster_centers_ : array, [n_clusters, n_features]
        Coordinates of cluster centers. If the algorithm stops before fully
        converging (see ``tol`` and ``max_iter``), these will not be
        consistent with ``labels_``.

    labels_ :
        Labels of each point

    Examples
    --------
    See also
    --------

    """

    def __init__(self, n_clusters=8, m=2, eps=10, init='random', max_iter=300, random_state=None):

        self.n_clusters = n_clusters
        self.init = init
        self.m = m
        self.eps = eps
        self.max_iter = max_iter
        self.random_state = random_state
        

    def _check_test_data(self, X):
        X = check_array(X, accept_sparse='csr', dtype=FLOAT_DTYPES)
        n_samples, n_features = X.shape
        expected_n_features = self.cluster_centers_.shape[1]
        if not n_features == expected_n_features:
            raise ValueError("Incorrect number of features. "
                             "Got %d features, expected %d" % (
                                 n_features, expected_n_features))

        return X

    def fit(self, X, y=None, sample_weight=None):
        """Compute Fuzzy c-means clustering.

        Parameters
        ----------
        X : array-like or sparse matrix, shape=(n_samples, n_features)
            Training instances to cluster. 

        y : Ignored
            not used, present here for API consistency by convention.

        sample_weight : array-like, shape (n_samples,), optional
            The weights for each observation in X. If None, all observations
            are assigned equal weight (default: None)

        """
        random_state = check_random_state(self.random_state)

        self.cluster_centers_, self.labels_, self.membership_mat = \
            fcm(
                X, 
                c_clusters = self.n_clusters,
                m = self.m,
                eps = self.eps,
                random_state = random_state, 
                max_iter=self.max_iter)
        return self

    def fit_predict(self, X, y=None, sample_weight=None):
        """Compute cluster centers and predict cluster index for each sample.

        Convenience method; equivalent to calling fit(X) followed by
        predict(X).

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            New data to transform.

        y : Ignored
            not used, present here for API consistency by convention.

        sample_weight : array-like, shape (n_samples,), optional
            The weights for each observation in X. If None, all observations
            are assigned equal weight (default: None)

        Returns
        -------
        labels : array, shape [n_samples,]
            Index of the cluster each sample belongs to.
        """
        return self.fit(X, sample_weight=sample_weight).labels_

    def fit_transform(self, X, y=None, sample_weight=None):
        """Compute clustering and transform X to cluster-probabiliy space.

        Equivalent to fit(X).transform(X), but more efficiently implemented.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            New data to transform.

        y : Ignored
            not used, present here for API consistency by convention.

        sample_weight : array-like, shape (n_samples,), optional
            The weights for each observation in X. If None, all observations
            are assigned equal weight (default: None)

        Returns
        -------
        X_new : array, shape [n_samples, k]
            X transformed in the new space.
        """
        # Currently, this just skips a copy of the data if it is not in
        # np.array or CSR format already.
        # XXX This skips _check_test_data, which may change the dtype;
        # we should refactor the input validation.
        return self.fit(X, sample_weight=sample_weight)._transform(X)

    def transform(self, X):
        """Transform X to a cluster-probabiliy space.

        In the new space, each dimension is the probality to each cluster
        centers.  Note that even if X is sparse, the array returned by
        `transform` will typically be dense.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            New data to transform.

        Returns
        -------
        X_new : array, shape [n_samples, k]
            X transformed in the new space.
        """
        check_is_fitted(self, 'cluster_centers_')

        X = self._check_test_data(X)
        return self._transform(X)

    def _transform(self, X):
        """guts of transform method; no input validation"""
        return self.membership_mat