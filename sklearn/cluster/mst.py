# -*- coding: utf-8 -*-
"""
MST Cluster: Compute Minimum Spanning Tree, cut weak links using a threshold.
"""

# Author: Robert Layton <robertlayton@gmail.com>
#
# License: 3-clause BSD.

from scipy.sparse import csr_matrix

from ..base import BaseEstimator, ClusterMixin
from ..metrics import pairwise_distances
from ..utils import atleast2d_or_csr
from ..utils.sparsetools import minimum_spanning_tree
from ..utils import connected_components


class MSTCluster(BaseEstimator, ClusterMixin):
    """Perform clustering using the minimum spanning tree and a threshold cut.

    A minimum spanning tree is created from the input data (interpreted as a
    graph) and then weak links are cut. The threshold parameter determines
    whether a link is cut.

    Parameters
    ----------
    X: array [n_samples, n_samples]
        Array of similarities between samples.

    threshold: float (default 0.5)
        The threshold to cut weak links. This algorithm is based on similarity,
        and therefore any links with distance lower than this values are cut.

    metric: string, or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string or callable, it must be one of
        the options allowed by metrics.pairwise.calculate_distance for its
        metric parameter.
        If metric is "precomputed", X is assumed to be a distance matrix and
        must be square.

    Attributes
    ----------
    `labels_` : array, shape = [n_samples]
        Cluster labels for each point in the dataset given to fit().

    Notes
    -----
    See examples/plot_eac.py for an example. The Evidence Accumulation
    Clustering (EAC) algorithm uses this clusterer in its final clustering
    step.

    References
    ----------
    Fred, Ana LN, and Anil K. Jain. "Data clustering using evidence
    accumulation." Pattern Recognition, 2002. Proceedings. 16th International
    Conference on. Vol. 4. IEEE, 2002.
    """

    def __init__(self, threshold=0.85, metric='euclidean'):
        self.threshold = threshold
        self.metric = metric

    def fit(self, X):
        """Perform MST clustering from a similarity matrix.

        Parameters
        ----------
        X: array (dense or sparse) [n_samples, n_samples]
            Array of similarities between samples.
        """
        X = atleast2d_or_csr(X)
        X = pairwise_distances(X, metric=self.metric)
        assert X.shape[0] == X.shape[1]
        span_tree = minimum_spanning_tree(X)
        idx = span_tree.data < self.threshold
        data = span_tree.data[idx]
        rows, cols = span_tree.nonzero()
        rows = rows[idx]
        cols = cols[idx]
        # Compute clusters by finding connected subgraphs in self.span_tree
        new_data = (data, (rows, cols))
        self.span_tree = csr_matrix(new_data, shape=span_tree.shape)
        n_components, self.labels_ = connected_components(self.span_tree,
                                                          directed=False)
        return self
