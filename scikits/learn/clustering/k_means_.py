""" K-means clustering
"""

# Authors: Gael Varoquaux <gael.xaroquaux@normalesup.org>
#          Thomas Rueckstiess <ruecksti@in.tum.de>
# License: BSD


import numpy as np

from scipy import cluster

from ..base import BaseEstimator

# kinit originaly from pybrain:
# http://github.com/pybrain/pybrain/raw/master/pybrain/auxiliary/kmeans.py
def kinit(X, k, n_samples_max=500):
    """ Init k seeds according to kmeans++

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
    if n_samples >= n_samples_max:
        X = X[np.random.randint(n_samples, size=n_samples_max)]
        n_samples = n_samples_max

    'choose the 1st seed randomly, and store D(x)^2 in D[]'
    centers = [X[np.random.randint(n_samples)]]
    D = ((X - centers[0])**2).sum(axis=-1)

    for _ in range(k - 1):
        bestDsum = bestIdx = -1

        for i in range(n_samples):
            'Dsum = sum_{x in X} min(D(x)^2,||x-xi||^2)'
            Dsum = np.minimum(D, ((X - X[i])**2).sum(axis=-1)
                              ).sum()

            if bestDsum < 0 or Dsum < bestDsum:
                bestDsum, bestIdx = Dsum, i

        centers.append(X[bestIdx])
        D = np.minimum(D, ((X - X[bestIdx])**2).sum(axis=-1))

    return np.array(centers)


def k_means(X, k, init='k-means++', n_iter=300, 
                        thresh=1e-5, missing='warn'):
    """ K-means clustering alorithm.

    Parameters
    ----------
    data : ndarray
        A M by N array of M observations in N dimensions or a length
        M array of M one-dimensional observations.
    k : int or ndarray
        The number of clusters to form as well as the number of
        centroids to generate. If minit initialization string is
        'matrix', or if a ndarray is given instead, it is
        interpreted as initial cluster to use instead.
    iter : int
        Number of iterations of the k-means algrithm to run. Note
        that this differs in meaning from the iters parameter to
        the kmeans function.
    thresh : float
        (not used yet).
    minit : {'k-means++', 'random', 'points', 'matrix'}
        Method for initialization:

        'random': generate k centroids from a Gaussian with mean and
        variance estimated from the data.

        'points': choose k observations (rows) at random from data for
        the initial centroids.

        'matrix': interpret the k parameter as a k by M (or length k
        array for one-dimensional data) array of initial centroids.

    Returns
    -------
    centroid : ndarray
        A k by N array of centroids found at the last iteration of
        k-means.
    label : ndarray
        label[i] is the code or index of the centroid the
        i'th observation is closest to.

    Notes
    ------
    This is currently scipy.cluster.vq.kmeans2 with the
    additional 'k-means++' initialization.

    """
    if init == 'k-means++':
        k = kinit(X, k)
        init='points'
    return cluster.vq.kmeans2(X, k, minit=init, missing=missing,
                                iter=n_iter)


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
        centroids to generate. If minit initialization string is
        'matrix', or if a ndarray is given instead, it is
        interpreted as initial cluster to use instead.
    iter : int
        Number of iterations of the k-means algrithm to run. Note
        that this differs in meaning from the iters parameter to
        the kmeans function.
    thresh : float
        (not used yet).
    minit : {'k-means++', 'random', 'points', 'matrix'}
        Method for initialization:

        'random': generate k centroids from a Gaussian with mean and
        variance estimated from the data.

        'points': choose k observations (rows) at random from data for
        the initial centroids.

        'matrix': interpret the k parameter as a k by M (or length k
        array for one-dimensional data) array of initial centroids.

    Methods
    -------

    fit(X):
        Compute MeanShift clustering

    Attributes
    ----------

    cluster_centers_: array, [n_clusters, n_features]
        Coordinates of cluster centers

    labels_:
        Labels of each point
    """


    def __init__(self, k=8, init='k-means++', n_iter=300, missing='warn'):
        self.k = k
        self.init = init
        self.n_iter = n_iter
        self.missing = missing

    def fit(self, X, **params):
        """ Compute k-means"""
        self._set_params(**params)
        self.cluster_centers_, self.labels_ = k_means(X, 
                    k=self.k, init=self.init, missing=self.missing,
                    n_iter=self.n_iter)
        return self
 
