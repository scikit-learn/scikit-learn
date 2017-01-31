import numpy as np


def _calculate_memberships(distances, m):
    distance_ratio = np.divide(distances[:, :, np.newaxis], distances[:, np.newaxis, :])
    distance_ratio[np.isnan(distance_ratio)] = 1.
    memberships = np.sum(np.power(distance_ratio, 2 / (m - 1)), axis=2) ** -1
    return memberships


def _calculate_centers(X, memberships, m):
    centers = np.dot(memberships.T ** m, X) / \
              np.sum(memberships.T ** m, axis=1)[..., np.newaxis]
    return centers