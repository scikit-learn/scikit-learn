""" K-means clustering
"""

# Authors: Gael Varoquaux <gael.xaroquaux@normalesup.org>
#          Thomas Rueckstiess <ruecksti@in.tum.de>
#          James Bergstra <james.bergstra@umontreal.ca>
# License: BSD

import warnings

import numpy as np

from ..base import BaseEstimator

def all_pairs_l2_distance_squared(A, B, B_norm_squared=None):
    """
    Returns the squared l2 norms of the differences between rows of A and B.

    Parameters
    ----------
    A: array, [n_rows_A, n_cols]

    B: array, [n_rows_B, n_cols]

    B_norm_squared: array [n_rows_B], or None
        pre-computed (B**2).sum(axis=1)

    Returns
    -------

    array [n_rows_A, n_rows_B]
        entry [i,j] is ((A[i] - B[i])**2).sum(axis=1)

    """
    if B_norm_squared is None:
        B_norm_squared = (B**2).sum(axis=1)
    if A is B:
        A_norm_squared = B_norm_squared
    else:
        A_norm_squared = (A**2).sum(axis=1)
    return (B_norm_squared + A_norm_squared.reshape((A.shape[0],1)) - 2*np.dot(A, B.T))

################################################################################
# Initialisation heuristic

# kinit originaly from pybrain:
# http://github.com/pybrain/pybrain/raw/master/pybrain/auxiliary/kmeans.py
def k_init(X, k, n_samples_max=500, rng=None):
    """Init k seeds according to kmeans++

    Parameters
    -----------
    X: array, shape (n_samples, n_features)
        The data

    k: integer
        The number of seeds to choose

    n_samples_max: integer, optional
        The maximum number of samples to use: the complexity of the
        algorithm is n_samples**2, if n_samples > n_samples_max,
        we use the Niquist strategy, and choose our centers in the
        n_samples_max samples randomly choosen.

    Notes
    ------
    Selects initial cluster centers for k-mean clustering in a smart way
    to speed up convergence. see: Arthur, D. and Vassilvitskii, S.
    "k-means++: the advantages of careful seeding". ACM-SIAM symposium
    on Discrete algorithms. 2007

    Implementation from Yong Sun's website
    http://blogs.sun.com/yongsun/entry/k_means_and_k_means
    """
    n_samples = X.shape[0]
    if rng is None:
        rng = np.random

    if n_samples >= n_samples_max:
        X = X[rng.randint(n_samples, size=n_samples_max)]
        n_samples = n_samples_max

    distances = all_pairs_l2_distance_squared(X, X)

    'choose the 1st seed randomly, and store D(x)^2 in D[]'
    first_idx =rng.randint(n_samples)
    centers = [X[first_idx]]
    D = distances[first_idx]

    for _ in range(k - 1):
        bestDsum = bestIdx = -1

        for i in range(n_samples):
            'Dsum = sum_{x in X} min(D(x)^2,||x-xi||^2)'
            Dsum = np.minimum(D, distances[i]).sum()

            if bestDsum < 0 or Dsum < bestDsum:
                bestDsum, bestIdx = Dsum, i

        centers.append(X[bestIdx])
        D = np.minimum(D, distances[bestIdx])

    return np.array(centers)


################################################################################
# K-means estimation by EM (expectation maximisation)

def k_means(X, k, init='k-means++', n_init=10, max_iter=300, verbose=0,
                    delta=1e-4, rng=None):
    """ K-means clustering algorithm.

    Parameters
    ----------
    X: ndarray
        A M by N array of M observations in N dimensions or a length
        M array of M one-dimensional observations.

    k: int or ndarray
        The number of clusters to form as well as the number of
        centroids to generate. If minit initialization string is
        'matrix', or if a ndarray is given instead, it is
        interpreted as initial cluster to use instead.

    max_iter: int, optional, default 300
        Maximum number of iterations of the k-means algorithm to run.

    n_init: int, optional, default: 10
        Number of time the k-means algorithm will be run with different centroid
        seeds. The final results will be the best output of n_init consecutive
        runs in terms of inertia.

    init: {'k-means++', 'random', or ndarray, or a callable}, optional
        Method for initialization, default to 'k-means++':

        'k-means++' : selects initial cluster centers for k-mean
        clustering in a smart way to speed up convergence. See section
        Notes in k_init for more details.

        'random': generate k centroids from a Gaussian with mean and
        variance estimated from the data.

        If an ndarray is passed, it should be of shape (k, p) and gives
        the initial centers.

    delta: float, optional
        The relative increment in the results before declaring convergence.

    verbose: boolean, optional
        Terbosity mode

    rng: numpy.RandomState, optional
        The generator used to initialize the centers

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
    if rng is None:
        rng = np.random
    n_samples = X.shape[0]

    vdata = np.mean(np.var(X, 0))
    best_inertia = np.infty
    if hasattr(init, '__array__'):
        init = np.asarray(init)
        if not n_init == 1:
            warnings.warn('Explicit initial center position passed: '
                          'performing only one init in the k-means')
            n_init = 1
    'subtract of mean of x for more accurate distance computations'
    Xmean = X.mean(axis=0)
    X = X-Xmean # TODO: offer an argument to allow doing this inplace
    for it in range(n_init):
        # init
        if init == 'k-means++':
            centers = k_init(X, k, rng=rng)
        elif init == 'random':
            seeds = np.argsort(rng.rand(n_samples))[:k]
            centers = X[seeds]
        elif hasattr(init, '__array__'):
            centers = np.asanyarray(init).copy()
        elif callable(init):
            centers = init(X, k, rng=rng)
        else:
            raise ValueError("the init parameter for the k-means should "
                "be 'k-means++' or 'random' or an ndarray, "
                "'%s' (type '%s') was passed.")

        if verbose:
            print 'Initialization complete' 
        # iterations
        for i in range(max_iter):
            centers_old = centers.copy()
            labels, inertia = _e_step(X, centers)
            centers = _m_step(X, labels, k)
            if verbose:
                print 'Iteration %i, intertia %s' % (i, inertia)
            if np.sum((centers_old - centers) ** 2) < delta * vdata:
                if verbose:
                    print 'Converged to similar centers at iteration', i
                break

            if inertia < best_inertia:
                best_centers = centers.copy()
                best_labels  = labels.copy()
                best_inertia = inertia
    else:
        best_centers = centers
        best_labels  = labels
        best_inertia = inertia
    return best_centers+Xmean, best_labels, best_inertia


def _m_step(x, z ,k):
    """ M step of the K-means EM algorithm

    Computation of cluster centers/means

    Parameters
    ----------
    x array of shape (n,p)
      n = number of samples, p = number of features
    z, array of shape (x.shape[0])
        Current assignment
    k, int
        Number of desired clusters

    Returns
    -------
    centers, array of shape (k, p)
        The resulting centers
    """
    dim = x.shape[1]
    centers = np.repeat(np.reshape(x.mean(0), (1, dim)), k, 0)
    for q in range(k):
        if np.sum(z==q)==0:
            pass
        else:
            centers[q] = np.mean(x[z==q], axis=0)
    return centers


def _e_step(x, centers):
    """E step of the K-means EM algorithm

    Computation of the input-to-cluster assignment

    Parameters
    ----------
    x: array of shape (n, p)
      n = number of samples, p = number of features

    centers: array of shape (k, p)
        The cluster centers

    Returns
    -------
    z: array of shape(n)
        The resulting assignment

    inertia: float
        The value of the inertia criterion with the assignment
    """

    n_samples = x.shape[0]
    k = centers.shape[0]

    there_is_memory_to_compute_distances_matrix = True

    if there_is_memory_to_compute_distances_matrix:
        distances = (
                (x**2).sum(axis=1)
                + (centers**2).sum(axis=1).reshape((k,1))
                - 2*np.dot(centers, x.T))
    # distances is a matrix of shape (k, n_samples) 

    z = -np.ones(n_samples).astype(np.int)
    mindist = np.infty * np.ones(n_samples)
    for q in range(k):
        if there_is_memory_to_compute_distances_matrix:
            dist = distances[q]
        else:
            dist = np.sum((x - centers[q]) ** 2, 1)
        z[dist<mindist] = q
        mindist = np.minimum(dist, mindist)
    inertia = mindist.sum()
    return z, inertia



################################################################################

class KMeans(BaseEstimator):
    """ K-Means clustering

    Parameters
    ----------

    data : ndarray
        A M by N array of M observations in N dimensions or a length
        M array of M one-dimensional observations.

    k : int or ndarray
        The number of clusters to form as well as the number of
        centroids to generate. If init initialization string is
        'matrix', or if a ndarray is given instead, it is
        interpreted as initial cluster to use instead.

    max_iter : int
        Maximum number of iterations of the k-means algorithm for a single run.

    n_init: int, optional, default: 10
        Number of time the k-means algorithm will be run with different centroid
        seeds. The final results will be the best output of n_init consecutive
        runs in terms of inertia.

    init : {'k-means++', 'random', 'points', 'matrix'}
        Method for initialization, defaults to 'k-means++':

        'k-means++' : selects initial cluster centers for k-mean
        clustering in a smart way to speed up convergence. See section
        Notes in k_init for more details.

        'random': generate k centroids from a Gaussian with mean and
        variance estimated from the data.

        'points': choose k observations (rows) at random from data for
        the initial centroids.

        'matrix': interpret the k parameter as a k by M (or length k
        array for one-dimensional data) array of initial centroids.

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


    def __init__(self, k=8, init='random', n_init=10, max_iter=300):
        self.k = k
        self.init = init
        self.max_iter = max_iter
        self.n_init = n_init

    def fit(self, X, **params):
        """ Compute k-means"""
        X = np.asanyarray(X)
        self._set_params(**params)
        self.cluster_centers_, self.labels_, self.inertia_ = k_means(X,
                    k=self.k, init=self.init, n_init=self.n_init,
                    max_iter=self.max_iter)
        return self

