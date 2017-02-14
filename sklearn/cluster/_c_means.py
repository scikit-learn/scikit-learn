import numpy as np


def _memberships_probabilistic(distances, m=2.):
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
    memberships = np.sum(np.power(distance_ratio, 2 / (m - 1)), axis=2) ** -1
    return memberships


def _centers_probabilistic(X, memberships, m=2.):
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

    """
    centers = np.dot(memberships.T ** m, X) / \
              np.sum(memberships.T ** m, axis=1)[..., np.newaxis]
    return centers