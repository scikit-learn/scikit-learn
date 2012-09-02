""" Algorithms for clustering : Meanshift,  Affinity propagation and spectral
clustering.

"""
# Author: Alexandre Gramfort alexandre.gramfort@inria.fr
#        Gael Varoquaux gael.varoquaux@normalesup.org

# License: BSD

import numpy as np
import warnings

from ..base import BaseEstimator, ClusterMixin
from ..utils import as_float_array
from ..metrics import euclidean_distances


def affinity_propagation(S, preference=None, p=None, convergence_iter=15,
        convit=None, max_iter=200,
        damping=0.5, copy=True, verbose=False):
    """Perform Affinity Propagation Clustering of data

    Parameters
    ----------

    S: array [n_samples, n_samples]
        Matrix of similarities between points

    preference: array [n_samples,] or float, optional, default: None
        Preferences for each point - points with larger values of
        preferences are more likely to be chosen as exemplars. The number of
        exemplars, i.e. of clusters, is influenced by the input preferences
        value. If the preferences are not passed as arguments, they will be
        set to the median of the input similarities (resulting in a moderate
        number of clusters). For a smaller amount of clusters, this can be set
        to the minimum value of the similarities.

    convergence_iter: int, optional, default: 15
        Number of iterations with no change in the number
        of estimated clusters that stops the convergence.

    max_iter: int, optional, default: 200
        Maximum number of iterations

    damping: float, optional, default: 200
        Damping factor between 0.5 and 1.

    copy: boolean, optional, default: True
        If copy is False, the affinity matrix is modified inplace by the
        algorithm, for memory efficiency

    verbose: boolean, optional, default: False
        The verbosity level

    Returns
    -------

    cluster_centers_indices: array [n_clusters]
        index of clusters centers

    labels : array [n_samples]
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
    if convit is not None:
        warnings.warn("``convit`` is deprecated and will be removed in"
                      "version 0.14. Use ``convergence_iter`` instead",
                      DeprecationWarning)
        convergence_iter = convit

    n_samples = S.shape[0]

    if S.shape[0] != S.shape[1]:
        raise ValueError("S must be a square array (shape=%r)" % S.shape)

    if not p is None:
        warnings.warn("p is deprecated and will be removed in version 0.14."
                " Use ``preference`` instead.", DeprecationWarning)
        preference = p

    if preference is None:
        preference = np.median(S)

    if damping < 0.5 or damping >= 1:
        raise ValueError('damping must be >= 0.5 and < 1')

    random_state = np.random.RandomState(0)

    # Place preference on the diagonal of S
    S.flat[::(n_samples + 1)] = preference

    A = np.zeros((n_samples, n_samples))
    R = np.zeros((n_samples, n_samples))  # Initialize messages

    # Remove degeneracies
    S += (np.finfo(np.double).eps * S + np.finfo(np.double).tiny * 100) * \
         random_state.randn(n_samples, n_samples)

    # Execute parallel affinity propagation updates
    e = np.zeros((n_samples, convergence_iter))

    ind = np.arange(n_samples)

    for it in range(max_iter):
        # Compute responsibilities
        Rold = R.copy()
        AS = A + S

        I = np.argmax(AS, axis=1)
        Y = AS[np.arange(n_samples), I]  # np.max(AS, axis=1)

        AS[ind, I[ind]] = - np.finfo(np.double).max

        Y2 = np.max(AS, axis=1)
        R = S - Y[:, np.newaxis]

        R[ind, I[ind]] = S[ind, I[ind]] - Y2[ind]

        R = (1 - damping) * R + damping * Rold  # Damping

        # Compute availabilities
        Aold = A
        Rp = np.maximum(R, 0)
        Rp.flat[::n_samples + 1] = R.flat[::n_samples + 1]

        A = np.sum(Rp, axis=0)[np.newaxis, :] - Rp

        dA = np.diag(A)
        A = np.minimum(A, 0)

        A.flat[::n_samples + 1] = dA

        A = (1 - damping) * A + damping * Aold  # Damping

        # Check for convergence
        E = (np.diag(A) + np.diag(R)) > 0
        e[:, it % convergence_iter] = E
        K = np.sum(E, axis=0)

        if it >= convergence_iter:
            se = np.sum(e, axis=1)
            unconverged = np.sum((se == convergence_iter) +\
                                 (se == 0)) != n_samples
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
        labels = np.empty((n_samples, 1))
        cluster_centers_indices = None
        labels.fill(np.nan)

    return cluster_centers_indices, labels


###############################################################################

class AffinityPropagation(BaseEstimator, ClusterMixin):
    """Perform Affinity Propagation Clustering of data

    Parameters
    ----------
    damping: float, optional, default: 0.5
        Damping factor between 0.5 and 1.

    convergence_iter: int, optional, default: 15
        Number of iterations with no change in the number
        of estimated clusters that stops the convergence.

    max_iter: int, optional, default: 200
        Maximum number of iterations

    copy: boolean, optional, default: True
        Make a copy of input data.

    preference: array [n_samples,] or float, optional, default: None
        Preferences for each point - points with larger values of
        preferences are more likely to be chosen as exemplars. The number
        of exemplars, ie of clusters, is influenced by the input
        preferences value. If the preferences are not passed as arguments,
        they will be set to the median of the input similarities.

    affinity: string, optional, default=``euclidean``
        Which affinity to use. At the moment ``precomputed`` and
        ``euclidean`` are supported. ``euclidean`` uses the
        negative squared euclidean distance between points.

    verbose: boolean, optional, default: False
        Whether to be verbose.


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

    def __init__(self, damping=.5, max_iter=200, convergence_iter=15,
            convit=None, copy=True,
            preference=None, p=None, affinity='euclidean', verbose=False):

        if convit is not None:
            warnings.warn("``convit`` is deprectaed and will be removed in "
                          "version 0.14. Use ``convergence_iter`` "
                          "instead", DeprecationWarning)
            convergence_iter = convit

        self.damping = damping
        self.max_iter = max_iter
        self.convergence_iter = convergence_iter
        self.copy = copy
        self.verbose = verbose
        if not p is None:
            warnings.warn("p is deprecated and will be removed in version 0.14"
                    ". Use ``preference`` instead.", DeprecationWarning)
            preference = p

        self.preference = preference
        self.affinity = affinity

    @property
    def _pairwise(self):
        return self.affinity is "precomputed"

    def fit(self, X):
        """ Create affinity matrix from negative euclidean distances, then
        apply affinity propagation clustering.

        Parameters
        ----------

        X: array [n_samples, n_features] or [n_samples, n_samples]
            Data matrix or, if affinity is ``precomputed``, matrix of
            similarities / affinities.
        """

        if X.shape[0] == X.shape[1] and not self._pairwise:
            warnings.warn("The API of AffinityPropagation has changed."
                "Now ``fit`` constructs an affinity matrix from the data."
                "To use a custom affinity matrix, set "
                "``affinity=precomputed``.")
        if self.affinity is "precomputed":
            self.affinity_matrix_ = X
        elif self.affinity is "euclidean":
            self.affinity_matrix_ = -euclidean_distances(X, squared=True)
        else:
            raise ValueError("Affinity must be 'precomputed' or "
                "'euclidean'. Got %s instead" % str(self.affinity))

        self.cluster_centers_indices_, self.labels_ = affinity_propagation(
                self.affinity_matrix_, self.preference,
                max_iter=self.max_iter,
                convergence_iter=self.convergence_iter,
                damping=self.damping, copy=self.copy, verbose=self.verbose)
        return self
