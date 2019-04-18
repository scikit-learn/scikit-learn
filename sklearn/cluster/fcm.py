"""Fuzzy C-Means clustering"""

# Authors: Junyi Li (lijy263@mail2.sysu.edu.cn)
# License: BSD 3 clause
# The import part this file is copy from k_means_.py, some functions is not used for now.
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


def fcm(X, c_clusters, m=2, eps=10, random_state=None, max_iter=300, sample_weight=None):
    
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
    
    X = check_array(X, accept_sparse='csr', dtype=FLOAT_DTYPES)

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

class FCM(BaseEstimator, ClusterMixin):
    """Fuzzy CMeans clustering

    Parameters
    ----------

    n_clusters : int, optional, default: 3
        The number of clusters to form as well as the number of
        centroids to generate.

    init : {'random'}
        Method for initialization, defaults to 'random':

        'random': choose k observations (rows) at random from data for
        the initial centroids.

        If an ndarray is passed, it should be of shape (n_clusters, n_features)
        and gives the initial centers.

    max_iter : int, default: 300
        Maximum number of iterations of the Fuzzy c means algorithm for a
        single run.

    random_state : int, RandomState instance or None (default)
        Determines random number generation for centroid initialization. Use
        an int to make the randomness deterministic.
        See :term:`Glossary <random_state>`.

    m  : int, optional, default: 2
        A key number for the update the centroids and the cluster-probability matrix.
    
    eps : float, optional, default: 10
        If the sum of |new_cluster_probability_matrix - old_cluster_probability_matrix| is smaller than eps, then the algorithm will stop.
        

    Attributes
    ----------
    cluster_centers_ : array, [n_clusters, n_features]
        Coordinates of cluster centers. 

    labels_ :
        Labels of each point

    Examples
    --------
    See also
    --------

    Note
    --------

    Now, something remains implementing:
        
        * sample weighted section
        * parallel run the model
        * ... 

    """

    def __init__(self, n_clusters=3, m=2, eps=10, init='random', max_iter=300, random_state=None):

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
        X = check_array(X, accept_sparse='csr')
        self.cluster_centers_, self.labels_, self.membership_mat = \
            fcm(
                X, 
                c_clusters = self.n_clusters,
                m = self.m,
                eps = self.eps,
                random_state = random_state, 
                max_iter=self.max_iter,
                sample_weight=sample_weight
            )
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