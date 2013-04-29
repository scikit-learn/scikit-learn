# -*- coding: utf-8 -*-
"""
EAC: Evidence Accumulation Clustering
"""

# Author: Robert Layton <robertlayton@gmail.com>
#
# License: 3-clause BSD.

import numpy as np

from sklearn.cluster import KMeans, SpectralClustering

from ..base import BaseEstimator, ClusterMixin
from ..utils import check_random_state, atleast2d_or_csr


def eac(X, initial_clusterers=None, final_clusterer=None, use_distance=False,
        random_state=None):
    """Perform EAC clustering from vector array or distance matrix.

    Parameters
    ----------
    X: array [n_samples, n_samples] or [n_samples, n_features]
        Array of distances between samples, or a feature array.
        The array is treated as a feature array unless the metric is given as
        'precomputed'.
    initial_clusterers: iterable, or None
        The clusterers used in the first step of the process. If an iterable is
        given, then each one is called. If None is given (default), 150 runs of
        k-means with k randomly selected between 10 and 30 are used.
    final_clusterer: model (extends ClusterMixin), or None
        The clusterer to apply to the final clustering matrix. The method must
        be able to take a coassociation matrix as input, which is an array of
        size [n_samples, n_samples].
        If None, the default model is used, which is MST.
    use_distance: boolean, or callable TODO: default should be false
        If True, convert the coassociation matrix to distance using
        `D=1./(C + 1)`. If callable, the function is called with the
        coassociation matrix as input. If False (default), then the matrix is
        given as input to the `final_clusterer`.
    random_state: numpy.RandomState, optional
        The generator used to initialize the initial_clusterers.
        Defaults to numpy.random.

    Returns
    -------
    final_model: model (extends ClusterMixin)
        The model given as `final_clusterer`, fitted with the evidence
        accumulated through this process.

    Notes
    -----
    See examples/plot_eac.py for an example.

    References
    ----------
    Fred, Ana LN, and Anil K. Jain. "Data clustering using evidence
    accumulation." Pattern Recognition, 2002. Proceedings. 16th International
    Conference on. Vol. 4. IEEE, 2002.
    """
    X = atleast2d_or_csr(X)
    n_samples = X.shape[0]
    # If index order not given, create random order.
    random_state = check_random_state(random_state)
    # If initial_clusterers is None, it is k-means 150 times with randomly
    # initialised k values (as per original paper).
    if initial_clusterers is None:
        initial_clusterers = _kmeans_random_k(n_samples, random_state)
    # If the final_clusterer is None, create the default model
    if final_clusterer is None:
        final_clusterer = create_default_final_clusterer(random_state)
    # Co-association matrix, originally zeros everywhere
    C = np.zeros((n_samples, n_samples), dtype='float')
    num_initial_clusterers = 0
    for model in initial_clusterers:
        num_initial_clusterers += 1
        # Update random state 
        # Fit model to X
        model.fit(X)
        # Calculate new coassociation matrix and add that to the tally
        C = update_coassociation_matrix(C, model.labels_)
    C /= num_initial_clusterers
    if use_distance:
        if use_distance is True:
            # Turn into a distance matrix
            C = 1. - C
        elif callable(use_distance):  # If a callable
            C = use_distance(C)
    final_clusterer.fit(C)
    return final_clusterer


def update_coassociation_matrix(C, labels):
    """Updates a co-association matrix from an array of labels.
    """
    labels = np.asarray(labels)
    for i in range(len(labels)):
        indices = np.where(labels[i:] == labels[i])[0] + i
        C[i, indices] += 1.
    return C


def _kmeans_random_k(n_samples, random_state=None, **kmeans_args):
    """Returns a generator for the default initial clustering for EAC
    
    This initial clustering is k-means, initialised randomly with k values
    chosen randomly between 10 and 30 inclusive.
    Parameters
    ----------
    random_state: numpy.RandomState, optional
        The generator used to initialize the initial_clusterers.
        Defaults to numpy.random.
    
    kmeans_args: other keywords
        Any additional arguments to this function are passed onto the
        initialiser for `sklearn.cluster.KMeans` with the exception of the
        `n_clusters` argument, which cannot be given (an error is raised if
        it is).
    
    Returns
    -------
    models: generator of KMeans instances
        Length will be 150, each instance initialised randomly using the
        supplied random_state.
    
    References
    ----------
    Fred, Ana LN, and Anil K. Jain. "Data clustering using evidence
    accumulation." Pattern Recognition, 2002. Proceedings. 16th International
    Conference on. Vol. 4. IEEE, 2002.
    """
    if 'n_clusters' in kmeans_args:
        error_msg = "n_clusters cannot be assigned for the default clusterers."
        raise ValueError(error_msg)
    random_state = check_random_state(random_state)
    num_iterations = 150
    k_low, k_high = (10, 30)
    if n_samples < k_high:
        k_high = n_samples
        k_low = min(k_low, int(k_high / 2))
    k_values = random_state.randint(k_low, high=k_high, size=num_iterations)
    return (KMeans(n_clusters=k, **kmeans_args) for k in k_values)


def create_default_final_clusterer(random_state=None):
    random_state = check_random_state(random_state)
    # TODO: MST clustering is the default in the paper, should use that, but it hasn't been implemented yet.
    return SpectralClustering(n_clusters=3, affinity='precomputed')


class EAC(BaseEstimator, ClusterMixin):
    """Perform EAC clustering from vector array or distance matrix.

    Evidence Accumulation Clustering (EAC) is an ensemble cluster that uses
    many iterations of k-means with randomly chosen k values (``n_clusters``)
    each time. The number of times two instances are clustered together is
    given in a co-association matrix, which is then clustered a final time to
    produce the 'final clustering'. In practice, this gives a more easily
    separable set of attributes that the original attributes.
    

    Parameters
    ----------
    X: array [n_samples, n_samples] or [n_samples, n_features]
        Array of distances between samples, or a feature array.
        The array is treated as a feature array unless the metric is given as
        'precomputed'.
    initial_clusterers: iterable, or None
        The clusterers used in the first step of the process. If an iterable is
        given, then each one is called. If None is given (default), 150 runs of
        k-means with k randomly selected between 10 and 30 are used.
    final_clusterer: model (extends ClusterMixin), or None
        The clusterer to apply to the final clustering matrix. The method must
        be able to take a coassociation matrix as input, which is an array of
        size [n_samples, n_samples].
        If None, the default model is used, which is MST.
    use_distance: boolean, or callable
        If True, convert the coassociation matrix to distance using
        `D=1./(C + 1)`. If callable, the function is called with the
        coassication matrix as input. If False (default), then the matrix is
        given as input to the `final_clusterer`.
    random_state: numpy.RandomState, optional
        The generator used to initialize the initial_clusterers.
        Defaults to numpy.random.

    Attributes
    ----------
    final_clusterer: model (same as the parameter final_clusterer)
        The clusterer given as a parameter, now fit to the evidence
        accumulation from learning.

    `labels_` : array, shape = [n_samples]
        Cluster labels for each point in the dataset given to fit().
        Same as the self.final_clusterer.labels_

    Notes
    -----
    See examples/plot_eac.py for an example.

    References
    ----------
    Fred, Ana LN, and Anil K. Jain. "Data clustering using evidence
    accumulation." Pattern Recognition, 2002. Proceedings. 16th International
    Conference on. Vol. 4. IEEE, 2002.
    """

    def __init__(self, initial_clusterers=None, final_clusterer=None,
                 use_distance=False, random_state=None):
        self.initial_clusterers = initial_clusterers
        self.final_clusterer = final_clusterer
        self.use_distance = use_distance
        self.random_state = random_state

    def fit(self, X):
        """Perform EAC clustering from vector array or distance matrix.

        Parameters
        ----------
        X: array [n_samples, n_samples] or [n_samples, n_features]
            Array of distances between samples, or a feature array.
            The array is treated as a feature array unless the metric is
            given as 'precomputed'.
        """
        final_clusterer = eac(X, final_clusterer=self.final_clusterer)
        self.labels_ = final_clusterer.labels_
        return self
