""" Algorithms for clustering : Meanshift,  Affinity propagation and spectral
clustering.

"""
# Author: Alexandre Gramfort alexandre.gramfort@inria.fr
#        Gael Varoquaux gael.varoquaux@normalesup.org

# License: BSD

import numpy as np
import warnings

from ..base import BaseEstimator
from ..utils import as_float_array
from ..metrics import euclidean_distances


def affinity_propagation(S, p=None, convit=30, max_iter=200, damping=0.5,
            copy=True, verbose=False):
    """Perform Affinity Propagation Clustering of data

    Parameters
    ----------

    S: array [n_points, n_points]
        Matrix of similarities between points

    p: array [n_points,] or float, optional
        Preferences for each point - points with larger values of
        preferences are more likely to be chosen as exemplars. The number of
        exemplars, ie of clusters, is influenced by the input preferences
        value. If the preferences are not passed as arguments, they will be
        set to the median of the input similarities (resulting in a moderate
        number of clusters). For a smaller amount of clusters, this can be set
        to the minimum value of the similarities.

    damping : float, optional
        Damping factor

    copy: boolean, optional
        If copy is False, the affinity matrix is modified inplace by the
        algorithm, for memory efficiency

    verbose: boolean, optional
        The verbosity level

    Returns
    -------

    cluster_centers_indices: array [n_clusters]
        index of clusters centers

    labels : array [n_points]
        cluster labels for each point

    Notes
    -----
    See examples/plot_affinity_propagation.py for an example.

    References
    ----------
    Brendan J. Frey and Delbert Dueck, "Clustering by Passing Messages
    Between Data Points", Science Feb. 2007
    """
    S = as_float_array(S, copy=copy)

    n_points = S.shape[0]

    if S.shape[0] != S.shape[1]:
        raise ValueError("S must be a square array (shape=%r)" % S.shape)

    if p is None:
        p = np.median(S)

    if damping < 0.5 or damping >= 1:
        raise ValueError('damping must be >= 0.5 and < 1')

    random_state = np.random.RandomState(0)

    # Place preferences on the diagonal of S
    S.flat[::(n_points + 1)] = p

    A = np.zeros((n_points, n_points))
    R = np.zeros((n_points, n_points))  # Initialize messages

    # Remove degeneracies
    S += (np.finfo(np.double).eps * S + np.finfo(np.double).tiny * 100) * \
         random_state.randn(n_points, n_points)

    # Execute parallel affinity propagation updates
    e = np.zeros((n_points, convit))

    ind = np.arange(n_points)

    for it in range(max_iter):
        # Compute responsibilities
        Rold = R.copy()
        AS = A + S

        I = np.argmax(AS, axis=1)
        Y = AS[np.arange(n_points), I]  # np.max(AS, axis=1)

        AS[ind, I[ind]] = - np.finfo(np.double).max

        Y2 = np.max(AS, axis=1)
        R = S - Y[:, np.newaxis]

        R[ind, I[ind]] = S[ind, I[ind]] - Y2[ind]

        R = (1 - damping) * R + damping * Rold  # Damping

        # Compute availabilities
        Aold = A
        Rp = np.maximum(R, 0)
        Rp.flat[::n_points + 1] = R.flat[::n_points + 1]

        A = np.sum(Rp, axis=0)[np.newaxis, :] - Rp

        dA = np.diag(A)
        A = np.minimum(A, 0)

        A.flat[::n_points + 1] = dA

        A = (1 - damping) * A + damping * Aold  # Damping

        # Check for convergence
        E = (np.diag(A) + np.diag(R)) > 0
        e[:, it % convit] = E
        K = np.sum(E, axis=0)

        if it >= convit:
            se = np.sum(e, axis=1)
            unconverged = np.sum((se == convit) + (se == 0)) != n_points
            if (not unconverged and (K > 0)) or (it == max_iter):
                if verbose:
                    print "Converged after %d iterations." % it
                break
    else:
        if verbose:
            print "Did not converged"

    I = np.where(np.diag(A + R) > 0)[0]
    K = I.size  # Identify exemplars

    if K > 0:
        c = np.argmax(S[:, I], axis=1)
        c[I] = np.arange(K)  # Identify clusters
        # Refine the final set of exemplars and clusters and return results
        for k in range(K):
            ii = np.where(c == k)[0]
            j = np.argmax(np.sum(S[ii[:, np.newaxis], ii], axis=0))
            I[k] = ii[j]

        c = np.argmax(S[:, I], axis=1)
        c[I] = np.arange(K)
        labels = I[c]
        # Reduce labels to a sorted, gapless, list
        cluster_centers_indices = np.unique(labels)
        labels = np.searchsorted(cluster_centers_indices, labels)
    else:
        labels = np.empty((n_points, 1))
        cluster_centers_indices = None
        labels.fill(np.nan)

    return cluster_centers_indices, labels


###############################################################################

class AffinityPropagation(BaseEstimator):
    """Perform Affinity Propagation Clustering of data

    Parameters
    ----------
    damping : float, optional
        Damping factor

    max_iter : int, optional
        Maximum number of iterations

    convit : int, optional
        Number of iterations with no change in the number
        of estimated clusters that stops the convergence.

    copy : boolean, optional
        Make a copy of input data. True by default.

    p : array [n_points,] or float, optional
        Preferences for each point - points with larger values of
        preferences are more likely to be chosen as exemplars. The number
        of exemplars, ie of clusters, is influenced by the input
        preferences value. If the preferences are not passed as arguments,
        they will be set to the median of the input similarities.

    precomputed : bool, optional
        Whether a precomputed affinity matrix is used for X, instead
        of the original data.


    Attributes
    ----------
    `cluster_centers_indices_` : array, [n_clusters]
        Indices of cluster centers

    `labels_` : array, [n_samples]
        Labels of each point

    `affinity_matrix_` : array-like, [n_samples, n_samples]
        Stores the affinity matrix used in ``fit``.

    Notes
    -----
    See examples/plot_affinity_propagation.py for an example.

    The algorithmic complexity of affinity propagation is quadratic
    in the number of points.

    References
    ----------

    Brendan J. Frey and Delbert Dueck, "Clustering by Passing Messages
    Between Data Points", Science Feb. 2007
    """

    def __init__(self, damping=.5, max_iter=200, convit=30, copy=True, p=None,
            precomputed=False):
        self.damping = damping
        self.max_iter = max_iter
        self.convit = convit
        self.copy = copy
        self.p = p
        self.precomputed = precomputed

    @property
    def _pairwise(self):
        return self.precomputed

    def fit(self, X):
        """ Create affinity matrix from negative euclidean distances, then
        apply affinity propagation clustering.

        Parameters
        ----------

        X: array [n_points, n_points]
            Matrix of similarities between points
        """

        if X.shape[0] == X.shape[1] and not self.precomputed:
            warnings.warn("The API of AffinityPropagation has changed."
                "Now ``fit`` constructs an affinity matrix from the data."
                "To use a custom affinity matrix, set ``precomputed=True``.")
        if self.precomputed:
            self.affinity_matrix_ = X
        else:
            self.affinity_matrix_ = -euclidean_distances(X, squared=True)

        self.cluster_centers_indices_, self.labels_ = affinity_propagation(
                self.affinity_matrix_, self.p,
                max_iter=self.max_iter, convit=self.convit,
                damping=self.damping, copy=self.copy)
        return self
