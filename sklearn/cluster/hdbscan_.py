# -*- coding: utf-8 -*-
"""
HDBSCAN: Hierarchical Density-Based Spatial Clustering 
         of Applications with Noise
"""
# Author: Leland McInnes <leland.mcinnes@gmail.com>
#         Steve Astels <sastels@gmail.com>
#         John Healy <jchealy@gmail.com>
#
# License: BSD 3 clause

import warnings

import numpy as np
import pandas as pd

from ..base import BaseEstimator, ClusterMixin
from ..metrics import pairwise_distances
from ..utils import check_array, check_consistent_length

from ._hdbscan_linkage import single_linkage
from ._hdbscan_tree import igraph_to_tree, \
                           igraph_condense_tree, \
                           igraph_tree_to_dataframe, \
                           igraph_get_clusters
                           
def mutual_reachability(distance_matrix, min_points=5):
    """Compute the weighted adjacency matrix of the mutual reachability
    graph of a distance matrix.
    
    Parameters
    ----------
    distance_matrix : array [n_samples, n_samples]
        Array of distances between samples.
        
    min_points : int optional
        The number of points in a neighbourhood for a point to be considered
        a core point. (defaults to 5)

    Returns
    -------
    mututal_reachability: array [n_samples, n_samples]
        Weighted adjacency matrix of the mutual reachability graph.
        
    References
    ----------
    R. Campello, D. Moulavi, and J. Sander, "Density-Based Clustering Based on
    Hierarchical Density Estimates"
    In: Advances in Knowledge Discovery and Data Mining, Springer, pp 160-172.
    2013
    """
    dim = distance_matrix.shape[0]
    core_distances = np.partition(distance_matrix, 
                                  min_points, 
                                  axis=0)[min_points]
                                  
    stage1 = np.where(core_distances > distance_matrix, 
                      core_distances, distance_matrix)
    result = np.where(core_distances > stage1.T,
                      core_distances.T, stage1.T).T
    return result
    
def compute_stability(cluster_tree):
    """Compute the stability of clusters.
    
    Parameters
    ----------
    cluster_tree : dataframe
        A dataframe giving the cluster tree specifying the parent id, child id, 
        and lambda value and size of splits in the tree.
        
    Returns
    -------
    stability : dataframe
        A dataframe indexed by clusters giving their stability values.
        
    References
    ----------
    R. Campello, D. Moulavi, and J. Sander, "Density-Based Clustering Based on
    Hierarchical Density Estimates"
    In: Advances in Knowledge Discovery and Data Mining, Springer, pp 160-172.
    2013
    """
    births = cluster_tree.groupby('child').min()[['lambda']]
    biths_and_deaths = cluster_tree.join(births,
                                         on='parent'
                                         lsuffix='_death',
                                         rsuffix='_birth')
    births_and_deaths['stability'] = (births_and_deaths['lambda_death'] -
                                      births_and_deaths['lambd_birth']) \
                                     * births_and_deaths['child_size']
    raw_stability = births_and_deaths.groupby('parent')[['stability']].sum()
    normalization = pd.DataFrame(births_and_deaths.parent.value_counts(),
                                 columns=['stability'])
    normalized_stability = raw_stability / normalization
    return normalized_stability
    
def hdbscan(X, min_cluster_size=5, min_samples=None, metric='minkowski', p=2):
    """Perform HDBSCAN clustering from a vector array or distance matrix.
    
    Parameters
    ----------
    X : array or sparse (CSR) matrix of shape (n_samples, n_features), or \
            array of shape (n_samples, n_samples)
        A feature array, or array of distances between samples if
        ``metric='precomputed'``.
        
    min_cluster_size : int optional
        The minimum number of samples in a groupo for that group to be 
        considered a cluster; groupings smaller than this size will be left
        as noise.

    min_samples : int, optional
        The number of samples in a neighborhood for a point
        to be considered as a core point. This includes the point itself.
        defaults to the min_cluster_size.

    metric : string, or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string or callable, it must be one of
        the options allowed by metrics.pairwise.pairwise_distances for its
        metric parameter.
        If metric is "precomputed", X is assumed to be a distance matrix and
        must be square.

    Returns
    -------
    labels : array [n_samples]
        Cluster labels for each point.  Noisy samples are given the label -1.
        
    References
    ----------
    R. Campello, D. Moulavi, and J. Sander, "Density-Based Clustering Based on
    Hierarchical Density Estimates"
    In: Advances in Knowledge Discovery and Data Mining, Springer, pp 160-172.
    2013
    """
    if min_samples is None:
        min_samples = min_cluster_size
        
    X = check_array(X, accept_sparse='csr')
    
    if metric == 'minkowski':
        if p is None:
            raise TypeError('Minkowski metrix given but no p value supplied!')
        distance_matrix = pairwise_distances(X, metric=metric, p=p)
    else:
        distance_matrix = pairwise_distances(X, metric=metric)
        
    mutual_reachability_graph = mutual_reachability(distance_matrix,
                                                    min_samples)
    raw_tree = single_linkage(mutual_reachability_graph)
    igraph_tree = igraph_to_tree(raw_tree)
    condensed_tree = igraph_condense_tree(igraph_tree, min_cluster_size)
    tree_dataframe = igraph_tree_to_dataframe(condensed_tree)
    stability_dict = compute_stability(tree_dataframe).to_dict()["stability"]
    cluster_list = igraph.get_clusters(condensed_tree, stability_dict)
    
    labels = -1 * np.ones(distance_matrix.shape[0])
    for index, cluster in enumerate(cluster_list):
        labels[cluster] = index
    return labels

class HDBSCAN(BaseEstimator, ClusterMixin):
    """Perform HDBSCAN clustering from vector array or distance matrix.
    
    HDBSCAN - Hierarchical Desnity-Based Spatial Clustering of Applications
    with Nopise. Performs DBSCAN over vvarying epsilon values and integrates 
    the result to find a clustering that gives the best stability over epsilon.
    This allows HDBSCAN to find clusters of varying densities (unlike DBSCAN),
    and be more robust to parameter selection.
    
    Parameters
    ----------
    min_cluster_size : int, optional
        The minumum size of clusters; single linkage splits that contain
        fewer points than this will be considered points "falling out" of a
        cluster rather than a cluster splitting into two new clusters.
        
    min_samples : int, optional
        The number of samples in a neighbourhood for a point to be
        considered a core point. (defaults to min_cluster_size)
        
    
    metric : string, or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string or callable, it must be one of
        the options allowed by metrics.pairwise.pairwise_distances for its
        metric parameter.
        If metric is "precomputed", X is assumed to be a distance matrix and
        must be square.
        
    Attributes
    ----------
    labels_ : array, shape = [n_samples]
        Cluster labels for each point in the dataset given to fit().
        Noisy samples are given the label -1.
        
    References
    ----------
    R. Campello, D. Moulavi, and J. Sander, "Density-Based Clustering Based on
    Hierarchical Density Estimates"
    In: Advances in Knowledge Discovery and Data Mining, Springer, pp 160-172.
    2013
    """
    
    def __init__(self, min_cluster_size=5, min_samples=None, 
                 metric='euclidean', p=None):
        self.min_cluster_size = min_cluster_size
        if min_samples is None:
            self.min_samples = min_cluster_size
        else:
            self.min_samples = min_samples
            
        self.metric = metric
        self.p = p
        
    def fit(self, X, y=None):
        """Perform HDBSCAN clustering from features or distance matrix.

        Parameters
        ----------
        X : array or sparse (CSR) matrix of shape (n_samples, n_features), or \
                array of shape (n_samples, n_samples)
            A feature array, or array of distances between samples if
            ``metric='precomputed'``.
        """
        X = check_array(X, accept_sparse='csr')
        self.labels_ = hdbscan(X, **self.get_params())
        return self
        
    def fit_predict(self, X, y=None, sample_weight=None):
        """Performs clustering on X and returns cluster labels.

        Parameters
        ----------
        X : array or sparse (CSR) matrix of shape (n_samples, n_features), or \
                array of shape (n_samples, n_samples)
            A feature array, or array of distances between samples if
            ``metric='precomputed'``.

        Returns
        -------
        y : ndarray, shape (n_samples,)
            cluster labels
        """
        self.fit(X, sample_weight=sample_weight)
        return self.labels_
