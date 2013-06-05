# -*- coding: utf-8 -*-
"""
MST Cluster: Compute Minimum Spanning Tree, cut weak links using a threshold.
"""

# Author: Robert Layton <robertlayton@gmail.com>
#
# License: 3-clause BSD.

import numpy as np

from scipy.sparse import cs_graph_components

from ..base import BaseEstimator, ClusterMixin
from ..utils import check_random_state, atleast2d_or_csr
from ..utils.mst import minimum_spanning_tree



class MSTCluster(BaseEstimator):
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

    def __init__(self, threshold=0.001):
        self.threshold = threshold

    def fit(self, X):
        """Perform MST clustering from a similarity matrix.

        Parameters
        ----------
        X: array (dense or sparse) [n_samples, n_samples]
            Array of similarities between samples.
        """
        X = atleast2d_or_csr(X)
        assert X.shape[0] == X.shape[1]
        # Compute spanning tree from X
        self.span_tree = minimum_spanning_tree(X)
        # Remote any link with a distance more than self.threshold
        print sorted(self.span_tree.data)
        self.span_tree.data[self.span_tree.data < self.threshold] = 0
        #self.span_tree = self.span_tree.trunc()
        print sorted(self.span_tree.data, reverse=True)[:5]
        # Compute clusters by finding connected subgraphs in self.span_tree
        n_components, self.labels_ = cs_graph_components(self.span_tree)
        return self
