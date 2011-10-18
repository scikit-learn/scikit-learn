"""K-means clustering"""

# Authors: Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Thomas Rueckstiess <ruecksti@in.tum.de>
#          James Bergstra <james.bergstra@umontreal.ca>
#          Jan Schlueter <scikit-learn@jan-schlueter.de>
#          Nelle Varoquaux
#          Peter Prettenhofer <peter.prettenhofer@gmail.com>
# License: BSD

import warnings
from itertools import cycle, izip

import numpy as np
import scipy.sparse as sp

from ..base import BaseEstimator
from ..metrics.pairwise import euclidean_distances
from ..utils import check_arrays
from ..utils import check_random_state
from ..utils import gen_even_slices
from ..utils import shuffle
from ..utils import warn_if_not_float

from . import _k_means


###############################################################################
# Initialisation heuristic


def k_init(X, k, n_local_trials=None, random_state=None, x_squared_norms=None):
    """Init k seeds according to kmeans++

    Parameters
    -----------
    X: array, shape (n_samples, n_features)
        The data to pick seeds for

    k: integer
        The number of seeds to choose

    n_local_trials: integer, optional
        The number of seeding trials for each center (except the first),
        of which the one reducing inertia the most is greedily chosen.
        Set to None to make the number of trials depend logarithmically
        on the number of seeds (2+log(k)); this is the default.

    random_state: numpy.RandomState, optional
        The generator used to initialize the centers. Defaults to numpy.random.

    x_squared_norms: array, shape (n_samples,), optional
        Squared euclidean norm of each data point. Pass it if you have it at
        hands already to avoid it being recomputed here. Default: None

    Notes
    ------
    Selects initial cluster centers for k-mean clustering in a smart way
    to speed up convergence. see: Arthur, D. and Vassilvitskii, S.
    "k-means++: the advantages of careful seeding". ACM-SIAM symposium
    on Discrete algorithms. 2007

    Version ported from http://www.stanford.edu/~darthur/kMeansppTest.zip,
    which is the implementation used in the aforementioned paper.
    """
    n_samples, n_features = X.shape
    random_state = check_random_state(random_state)

    centers = np.empty((k, n_features))

    # Set the number of local seeding trials if none is given
    if n_local_trials is None:
        # This is what Arthur/Vassilvitskii tried, but did not report
        # specific results for other than mentioning in the conclusion
        # that it helped.
        n_local_trials = 2 + int(np.log(k))

    # Pick first center randomly
    center_id = random_state.randint(n_samples)
    centers[0] = X[center_id]

    # Initialize list of closest distances and calculate current potential
    if x_squared_norms is None:
        x_squared_norms = X.copy()
        x_squared_norms **= 2
        x_squared_norms = x_squared_norms.sum(axis=1)
    closest_dist_sq = euclidean_distances(
        np.atleast_2d(centers[0]), X, Y_norm_squared=x_squared_norms,
        squared=True)
    current_pot = closest_dist_sq.sum()

    # Pick the remaining k-1 points
    for c in xrange(1, k):
        # Choose center candidates by sampling with probability proportional
        # to the squared distance to the closest existing center
        rand_vals = random_state.random_sample(n_local_trials) * current_pot
        candidate_ids = np.searchsorted(closest_dist_sq.cumsum(), rand_vals)

        # Compute distances to center candidates
        distance_to_candidates = euclidean_distances(
            X[candidate_ids], X, Y_norm_squared=x_squared_norms, squared=True)

        # Decide which candidate is the best
        best_candidate = None
        best_pot = None
        best_dist_sq = None
        for trial in xrange(n_local_trials):
            # Compute potential when including center candidate
            new_dist_sq = np.minimum(closest_dist_sq,
                                     distance_to_candidates[trial])
            new_pot = new_dist_sq.sum()

            # Store result if it is the best local trial so far
            if (best_candidate is None) or (new_pot < best_pot):
                best_candidate = candidate_ids[trial]
                best_pot = new_pot
                best_dist_sq = new_dist_sq

        # Permanently add best center candidate found in local tries
        centers[c] = X[best_candidate]
        current_pot = best_pot
        closest_dist_sq = best_dist_sq

    return centers


###############################################################################
# K-means estimation by EM (expectation maximisation)


def k_means(X, k, init='k-means++', n_init=10, max_iter=300, verbose=0,
                    tol=1e-4, random_state=None, copy_x=True):
    """ K-means clustering algorithm.

    Parameters
    ----------
    X: ndarray
        A M by N array of M observations in N dimensions or a length
        M array of M one-dimensional observations.

    k: int or ndarray
        The number of clusters to form as well as the number of
        centroids to generate.

    max_iter: int, optional, default 300
        Maximum number of iterations of the k-means algorithm to run.

    n_init: int, optional, default: 10
        Number of time the k-means algorithm will be run with different
        centroid seeds. The final results will be the best output of
        n_init consecutive runs in terms of inertia.

    init: {'k-means++', 'random', or ndarray, or a callable}, optional
        Method for initialization, default to 'k-means++':

        'k-means++' : selects initial cluster centers for k-mean
        clustering in a smart way to speed up convergence. See section
        Notes in k_init for more details.

        'random': generate k centroids from a Gaussian with mean and
        variance estimated from the data.

        If an ndarray is passed, it should be of shape (k, p) and gives
        the initial centers.

        If a callable is passed, it should take arguments X, k and
        and a random state and return an initialization.

    tol: float, optional
        The relative increment in the results before declaring convergence.

    verbose: boolean, optional
        Terbosity mode

    random_state: numpy.RandomState, optional
        The generator used to initialize the centers. Defaults to numpy.random.

    copy_x: boolean, optional
        When pre-computing distances it is more numerically accurate to center
        the data first.  If copy_x is True, then the original data is not
        modified.  If False, the original data is modified, and put back before
        the function returns, but small numerical differences may be introduced
        by subtracting and then adding the data mean.

    Returns
    -------
    centroid: ndarray
        A k by N array of centroids found at the last iteration of
        k-means.

    label: ndarray
        label[i] is the code or index of the centroid the
        i'th observation is closest to.

    inertia: float
        The final value of the inertia criterion

    """
    random_state = check_random_state(random_state)

    vdata = np.mean(np.var(X, 0))
    best_inertia = np.infty

    # subtract of mean of x for more accurate distance computations
    X_mean = X.mean(axis=0)
    if copy_x:
        X = X.copy()
    X -= X_mean

    if hasattr(init, '__array__'):
        init = np.asarray(init).copy()
        init -= X_mean
        if not n_init == 1:
            warnings.warn('Explicit initial center position passed: '
                          'performing only one init in the k-means')
            n_init = 1

    # precompute squared norms of data points
    x_squared_norms = X.copy()
    x_squared_norms **= 2
    x_squared_norms = x_squared_norms.sum(axis=1)
    for it in range(n_init):
        # init
        centers = _init_centroids(X, k, init, random_state=random_state,
                                  x_squared_norms=x_squared_norms)
        if verbose:
            print 'Initialization complete'

        # iterations
        for i in range(max_iter):
            centers_old = centers.copy()
            labels, inertia = _e_step(X, centers,
                                      x_squared_norms=x_squared_norms)
            centers = _m_step(X, labels, k)

            if verbose:
                print 'Iteration %i, inertia %s' % (i, inertia)

            if np.sum((centers_old - centers) ** 2) < tol * vdata:
                if verbose:
                    print 'Converged to similar centers at iteration', i
                break

            if inertia < best_inertia:
                best_labels = labels.copy()
                best_centers = centers.copy()
                best_inertia = inertia

    else:
        best_labels = labels
        best_centers = centers
        best_inertia = inertia
    if not copy_x:
        X += X_mean
    return best_centers + X_mean, best_labels, best_inertia


def _calculate_labels_inertia(X, centers, x_squared_norms=None):
    """Compute the inertia and the labels of the given samples and centers"""
    distance = euclidean_distances(centers, X, x_squared_norms, squared=True)
    return distance.min(axis=0).sum(), distance.argmin(axis=0)


def _m_step(X, z, k):
    """M step of the K-means EM algorithm

    Computation of cluster centers/means

    Parameters
    ----------
    X: array, shape (n_samples, n_features)

    z: array, shape (n_samples)
        Current assignment

    k: int
        Number of desired clusters

    Returns
    -------
    centers: array, shape (k, n_features)
        The resulting centers
    """
    dim = X.shape[1]
    centers = np.empty((k, dim))
    X_center = None
    for q in range(k):
        center_mask = (z == q)
        if not np.any(center_mask):
            # The centroid of empty clusters is set to the center of
            # everything
            if X_center is None:
                X_center = X.mean(axis=0)
            centers[q] = X_center
        else:
            centers[q] = np.mean(X[center_mask], axis=0)
    return centers


def _init_centroids(X, k, init, random_state=None, x_squared_norms=None):
    """Compute the initial centroids

    Parameters
    ----------

    X: array, shape (n_samples, n_features)

    k: int
        number of centroids

    init: {'k-means++', 'random' or ndarray or callable} optional
        Method for initialisation

    random_state: numpy.RandomState, optional
        The generator used to initialise the centers. Defaults to numpy.random

    x_squared_norms:  array, shape (n_samples,), optional
        Squared euclidean norm of each data point. Pass it if you have it at
        hands already to avoid it being recomputed here. Default: None

    Returns
    -------
    centers: array, shape(k, n_features)
    """
    random_state = check_random_state(random_state)

    n_samples = X.shape[0]
    if init == 'k-means++':
        if sp.issparse(X):
            raise ValueError("Init method 'k-means++' only for dense X.")
        centers = k_init(X, k,
                        random_state=random_state,
                        x_squared_norms=x_squared_norms)
    elif init == 'random':
        seeds = np.argsort(random_state.rand(n_samples))[:k]
        centers = X[seeds]
    elif hasattr(init, '__array__'):
        centers = init
    elif callable(init):
        centers = init(X, k, random_state=random_state)
    else:
        raise ValueError("the init parameter for the k-means should "
            "be 'k-means++' or 'random' or an ndarray, "
            "'%s' (type '%s') was passed.")

    if sp.issparse(centers):
        centers = centers.toarray()
    return centers


def _e_step(x, centers, precompute_distances=True, x_squared_norms=None):
    """E step of the K-means EM algorithm

    Computation of the input-to-cluster assignment

    Parameters
    ----------
    x: array, shape (n_samples, n_features)

    centers: array, shape (k, n_features)
        The cluster centers

    precompute_distances: bool, optional
        Whether to compute the full distance matrix between centers and data
        points at once for more speed at the cost of memory. Default: True

    x_squared_norms: array, shape (n_samples,), optional
        Squared euclidean norm of each data point, speeds up computations in
        case of precompute_distances == True. Default: None

    Returns
    -------
    z: array of shape(n)
        The resulting assignment

    inertia: float
        The value of the inertia criterion with the assignment
    """

    n_samples = x.shape[0]
    k = centers.shape[0]

    if precompute_distances:
        distances = euclidean_distances(centers, x, x_squared_norms,
                                        squared=True)
    z = np.empty(n_samples, dtype=np.int)
    z.fill(-1)
    mindist = np.empty(n_samples)
    mindist.fill(np.infty)
    for q in range(k):
        if precompute_distances:
            dist = distances[q]
        else:
            dist = np.sum((x - centers[q]) ** 2, axis=1)
        z[dist < mindist] = q
        mindist = np.minimum(dist, mindist)
    inertia = mindist.sum()
    return z, inertia


class KMeans(BaseEstimator):
    """K-Means clustering

    Parameters
    ----------

    k : int, optional, default: 8
        The number of clusters to form as well as the number of
        centroids to generate.

    max_iter : int
        Maximum number of iterations of the k-means algorithm for a
        single run.

    n_init: int, optional, default: 10
        Number of time the k-means algorithm will be run with different
        centroid seeds. The final results will be the best output of
        n_init consecutive runs in terms of inertia.

    init : {'k-means++', 'random' or an ndarray}
        Method for initialization, defaults to 'k-means++':

        'k-means++' : selects initial cluster centers for k-mean
        clustering in a smart way to speed up convergence. See section
        Notes in k_init for more details.

        'random': choose k observations (rows) at random from data for
        the initial centroids.

        if init is an 2d array, it is used as a seed for the centroids

    tol: float, optional default: 1e-4
        Relative tolerance w.r.t. inertia to declare convergence


    Methods
    -------

    fit(X):
        Compute K-Means clustering

    Attributes
    ----------

    cluster_centers_: array, [n_clusters, n_features]
        Coordinates of cluster centers

    labels_:
        Labels of each point

    inertia_: float
        The value of the inertia criterion associated with the chosen
        partition.

    Notes
    ------

    The k-means problem is solved using the Lloyd algorithm.

    The average complexity is given by O(k n T), were n is the number of
    samples and T is the number of iteration.

    The worst case complexity is given by O(n^(k+2/p)) with
    n = n_samples, p = n_features. (D. Arthur and S. Vassilvitskii,
    'How slow is the k-means method?' SoCG2006)

    In practice, the K-means algorithm is very fast (one of the fastest
    clustering algorithms available), but it falls in local minima. That's why
    it can be useful to restart it several times.
    """

    def __init__(self, k=8, init='k-means++', n_init=10, max_iter=300,
                 tol=1e-4, verbose=0, random_state=None, copy_x=True):

        if hasattr(init, '__array__'):
            k = init.shape[0]

        self.k = k
        self.init = init
        self.max_iter = max_iter
        self.tol = tol
        self.n_init = n_init
        self.verbose = verbose
        self.random_state = random_state
        self.copy_x = copy_x

    def _check_data(self, X):
        """Verify that the number of samples given is larger than k"""
        if sp.issparse(X):
            raise ValueError("K-Means does not support sparse input matrices.")
        X = np.asanyarray(X)
        if X.shape[0] < self.k:
            raise ValueError("n_samples=%d should be >= k=%d" % (
                X.shape[0], self.k))
        return X

    def fit(self, X, y=None):
        """Compute k-means"""
        self.random_state = check_random_state(self.random_state)

        X = self._check_data(X)
        warn_if_not_float(X, self)

        self.cluster_centers_, self.labels_, self.inertia_ = k_means(
            X, k=self.k, init=self.init, n_init=self.n_init,
            max_iter=self.max_iter, verbose=self.verbose,
            tol=self.tol, random_state=self.random_state, copy_x=self.copy_x)
        return self

    def transform(self, X, y=None):
        """ Transform the data to a cluster-distance space

        In the new space, each dimension is the distance to the cluster centers.
        Note that even if X is sparse, the array returned by `transform` will
        typically be dense.

        Parameters
        ----------
        X: {array-like, sparse matrix}, shape = [n_samples, n_features]
            New data to transform.

        Returns
        -------
        X_new : array, shape [n_samples, k]
            X transformed in the new space.
        """
        if not hasattr(self, "cluster_centers_"):
            raise AttributeError("Model has not been trained. "
                                 "Train k-means before using transform.")
        cluster_shape = self.cluster_centers_.shape[1]
        if not X.shape[1] == cluster_shape:
            raise ValueError("Incorrect number of features for points. "
                             "Got %d features, expected %d" % (X.shape[1],
                                                               cluster_shape))
        return euclidean_distances(X, self.cluster_centers_)

    def predict(self, X):
        """Predict the closest cluster each sample in X belongs to.

        In the vector quantization literature, `cluster_centers_` is called
        the code book and each value returned by `predict` is the index of
        the closest code in the code book.

        Parameters
        ----------
        X: {array-like, sparse matrix}, shape = [n_samples, n_features]
            New data to predict.

        Returns
        -------
        Y : array, shape [n_samples,]
            Index of the closest center each sample belongs to.
        """
        if not hasattr(self, "cluster_centers_"):
            raise AttributeError("Model has not been trained yet. "
                                 "Fit k-means before using predict.")
        expected_n_features = self.cluster_centers_.shape[1]
        if not X.shape[1] == expected_n_features:
            raise ValueError("Incorrect number of features. "
                             "Got %d features, expected %d" % (
                                 X.shape[1], expected_n_features))
        return _e_step(X, self.cluster_centers_)[0]


def _mini_batch_step_dense(X, batch_slice, centers, counts, x_squared_norms):
    """Incremental update of the centers for the Minibatch K-Means algorithm

    Parameters
    ----------

    X: array, shape (n_samples, n_features)
        The original data array.

    batch_slice: slice
        The row slice of the mini batch.

    centers: array, shape (k, n_features)
        The cluster centers. This array is MODIFIED IN PLACE

    counts: array, shape (k, )
         The vector in which we keep track of the numbers of elements in a
         cluster. This array is MODIFIED IN PLACE

    x_squared_norms: array, shape (n_samples,)
        Squared euclidean norm of each data point.
    """
    # This is inefficient but saves mem and fits to sparse matrices.
    X = X[batch_slice]
    x_squared_norms = x_squared_norms[batch_slice]

    cache = euclidean_distances(centers, X, Y_norm_squared=x_squared_norms,
                                squared=True).argmin(axis=0)

    k = centers.shape[0]
    for q in range(k):
        center_mask = (cache == q)
        c = center_mask.sum()
        if np.any(center_mask):
            centers[q] = (1. / (counts[q] + c)) * (
                counts[q] * centers[q] + np.sum(X[center_mask], axis=0))
            counts[q] += c
    return counts, centers


def _mini_batch_step_sparse(X, batch_slice, centers, counts, x_squared_norms):
    """Incremental update of the centers for the Minibatch K-Means algorithm

    Parameters
    ----------

    X: csr_matrix, shape (n_samples, n_features)
        The data matrix in sparse CSR format.

    batch_slice: slice
        The row slice of the mini batch.

    centers: array, shape (k, n_features)
        The cluster centers. This array is MODIFIED IN PLACE

    counts: array, shape (k, )
         The vector in which we keep track of the numbers of elements in a
         cluster. This array is MODIFIED IN PLACE

    x_squared_norms: array, shape (n_samples,)
         The squared norms of each sample in `X`.
    """
    cache = euclidean_distances(centers, X[batch_slice],
              x_squared_norms[batch_slice]).argmin(axis=0).astype(np.int32)

    _k_means._mini_batch_update_sparse(X.data, X.indices, X.indptr,
                                       batch_slice, centers, counts, cache)
    return counts, centers


class MiniBatchKMeans(KMeans):
    """Mini-Batch K-Means clustering

    Parameters
    ----------

    k : int, optional, default: 8
        The number of clusters to form as well as the number of
        centroids to generate.

    max_iter : int
        Maximum number of iterations of the k-means algorithm for a
        single run.

    chunk_size: int, optional, default: 1000
        Size of the mini batches

    init : {'k-means++', 'random' or an ndarray}
        Method for initialization, defaults to 'random':

        'k-means++' : selects initial cluster centers for k-mean
        clustering in a smart way to speed up convergence. See section
        Notes in k_init for more details. Only for dense `X`.

        'random': choose k observations (rows) at random from data for
        the initial centroids.

        if init is an 2d array, it is used as a seed for the centroids

    tol: float, optional default: 1e-4
        Relative tolerance w.r.t. inertia to declare convergence

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

    inertia_: float
        The value of the inertia criterion associated with the chosen
        partition.

    References
    ----------
    http://www.eecs.tufts.edu/~dsculley/papers/fastkmeans.pdf
    """

    def __init__(self, k=8, init='random', max_iter=100,
                 chunk_size=1000, tol=1e-4, verbose=0, random_state=None):

        super(MiniBatchKMeans, self).__init__(k, init, 1,
              max_iter, tol, verbose, random_state)

        self.counts = None
        self.cluster_centers_ = None
        self.chunk_size = chunk_size

    def fit(self, X, y=None):
        """Compute the centroids on X by chunking it into mini-batches.

        Parameters
        ----------
        X: array-like, shape = [n_samples, n_features]
            Coordinates of the data points to cluster
        """
        self.random_state = check_random_state(self.random_state)
        X = check_arrays(X, sparse_format="csr", copy=False)[0]
        warn_if_not_float(X, self)
        n_samples, n_features = X.shape
        if n_samples < self.k:
            raise ValueError("Number of samples smaller than number "\
                             "of clusters.")

        if hasattr(self.init, '__array__'):
            self.init = np.asarray(self.init)

        X_shuffled = shuffle(X, random_state=self.random_state)

        if sp.issparse(X_shuffled):
            x_squared_norms = _k_means.csr_row_norm_l2(X)
        else:
            x_squared_norms = np.sum(X ** 2.0, axis=1)

        self.cluster_centers_ = _init_centroids(
            X_shuffled, self.k, self.init, random_state=self.random_state,
            x_squared_norms=x_squared_norms)
        self.counts = np.zeros(self.k, dtype=np.int32)

        n_batches = int(np.ceil(float(n_samples) / self.chunk_size))
        batch_slices = list(gen_even_slices(n_samples, n_batches))
        n_iterations = xrange(int(self.max_iter * n_batches))
        if sp.issparse(X_shuffled):
            _mini_batch_step = _mini_batch_step_sparse
            tol = self.tol
        else:
            _mini_batch_step = _mini_batch_step_dense
            tol = np.mean(np.var(X_shuffled, axis=0)) * self.tol

        for i, batch_slice in izip(n_iterations, cycle(batch_slices)):
            old_centers = self.cluster_centers_.copy()
            self.counts, self.cluster_centers_ = _mini_batch_step(
                            X_shuffled, batch_slice, 
                            self.cluster_centers_, self.counts, 
                            x_squared_norms=x_squared_norms)

            if np.sum((old_centers - self.cluster_centers_) ** 2) < tol:
                if self.verbose:
                    print 'Converged to similar centers at iteration', i
                break

        self.inertia_ = 0
        self.labels_ = np.empty((n_samples,), dtype=np.int)
        for batch_slice in batch_slices:
            batch_inertia, batch_labels = _calculate_labels_inertia(
            X[batch_slice], self.cluster_centers_)
            self.inertia_ += batch_inertia
            self.labels_[batch_slice] = batch_labels

        return self

    def partial_fit(self, X, y=None):
        """Update k means estimate on a single mini-batch X.

        Parameters
        ----------
        X: array-like, shape = [n_samples, n_features]
            Coordinates of the data points to cluster.
        """
        self.random_state = check_random_state(self.random_state)

        X = check_arrays(X, sparse_format="csr", copy=False)[0]
        n_samples, n_features = X.shape
        if hasattr(self.init, '__array__'):
            self.init = np.asarray(self.init)

        if n_samples == 0:
            return self

        if sp.issparse(X):
            x_squared_norms = _k_means.csr_row_norm_l2(X)
        else:
            x_squared_norms = (X ** 2).sum(axis=1)

        if self.counts is None:
            # this is the first call partial_fit on this object:
            # initialize the cluster centers
            self.cluster_centers_ = _init_centroids(
                X, self.k, self.init, random_state=self.random_state,
                x_squared_norms=x_squared_norms)

            self.counts = np.zeros(self.k, dtype=np.int32)

        batch_slice = slice(0, n_samples, None)
        if sp.issparse(X):
            _mini_batch_step = _mini_batch_step_sparse
        else:
            _mini_batch_step = _mini_batch_step_dense

        self.counts, self.cluster_centers_ = _mini_batch_step(X, 
                        batch_slice, self.cluster_centers_, self.counts,
                        x_squared_norms=x_squared_norms)

        self.inertia_, self.labels_ = _calculate_labels_inertia(
            X, self.cluster_centers_, x_squared_norms)

        return self
