import numpy as np
from sklearn.metrics.pairwise import euclidean_distances, check_pairwise_arrays


class Probabilistic:

    @staticmethod
    def distances(X, Y=None,  **kwargs):
        return euclidean_distances(X, Y, **kwargs)

    @staticmethod
    def memberships(distances, m=2.):
        """Calculate probabilistic memberships based on distances.

        Parameters
        ----------
        distances : array-like, shape (n_samples, n_samples)
            The distances from each point to every other.
        m : float, optional, default 2.0
            The fuzziness factor. Larger values dilute membership.

        Returns
        -------
        memberships : float ndarray, shape (n_samples, n_clusters)
            Updated memberships

        """
        distance_ratio = np.divide(distances[:, :, np.newaxis],
                                   distances[:, np.newaxis, :])
        distance_ratio[np.isnan(distance_ratio)] = 1.
        memberships = np.sum(np.power(distance_ratio, 2 / (m - 1)),
                             axis=2) ** -1
        return memberships

    @staticmethod
    def centers(X, memberships, m=2.):
        """Calculate probabilistic centers based on memberships.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            The data points to cluster.
        memberships : float ndarray, shape (n_samples, n_clusters)
            Membership of each data point to the cluster centers.
        m : float, optional, default 2.0
            The fuzziness factor. Larger values dilute membership.

        Returns
        -------
        centers : float ndarray, shape (n_clusters, n_features)
            The positions of the cluster centers.

        """
        centers = np.dot(memberships.T ** m, X) / \
            np.sum(memberships.T ** m, axis=1)[..., np.newaxis]
        return centers


class Possibilistic:

    def __init__(self, weights):
        self.weights = weights

    @staticmethod
    def distances(X, Y=None, **kwargs):
        return euclidean_distances(X, Y, **kwargs)

    def memberships(self, distances, m=2.):
        """Calculate possibilistic memberships based on distances.

        Parameters
        ----------
        distances : array-like, shape (n_samples, n_samples)
            The distances from each point to every other.
        weights : array-like, shape (n_clusters, )
            The relative weighting of the clusters.
        m : float, optional, default 2.0
            The fuzziness factor. Larger values dilute membership.

        Returns
        -------
        memberships : float ndarray, shape (n_samples, n_clusters)
            Updated memberships

        """
        memberships = (1. + (distances / self.weights) ** (1. / (m - 1))) ** -1.
        return memberships

    @staticmethod
    def centers(X, memberships, m=2.):
        """Calculate possibilistic centers based on memberships.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            The data points to cluster.
        memberships : float ndarray, shape (n_samples, n_clusters)
            Membership of each data point to the cluster centers.
        m : float, optional, default 2.0
            The fuzziness factor. Larger values dilute membership.

        Returns
        -------
        centers : float ndarray, shape (n_clusters, n_features)
            The positions of the cluster centers.

        """
        return np.divide(np.dot((memberships ** m).T, X),
                         np.sum(memberships ** m, axis=0)[..., np.newaxis])


class GustafsonKessel(Probabilistic):

    def __init__(self, covariances=None):
        self.covariances_ = covariances

    def distances(self, X, Y, **kwargs):
        if self.covariances_ is None:
            raise ValueError('Cannot compute Gustafson-Kessel'
                             ' distance without covariance matrices')
        v = X - Y[:, np.newaxis]
        m_inv = (np.linalg.det(self.covariances_) ** (-1 / X.shape[1]))[
                    ..., np.newaxis, np.newaxis] * self.covariances_
        m = np.linalg.inv(m_inv)
        return np.einsum('...ki,...ij,...kj->...k', v, m, v).T

    @staticmethod
    def covariances(X, centers, memberships, m=2.):
        v = X - centers[:, np.newaxis]
        numerator = np.einsum('k...,...ki,...kj->...ij', memberships ** m, v, v)
        dnominator = np.sum(memberships ** m, axis=0)[..., np.newaxis, np.newaxis]
        return numerator / dnominator

