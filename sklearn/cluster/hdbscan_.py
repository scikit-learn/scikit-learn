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

from ..base import BaseEstimator, ClusterMixin
from ..metrics import pairwise_distances
from ..utils import check_array, check_consistent_length

from ._hdbscan_linkage import single_linkage
from ._hdbscan_tree import (get_points,
                            condense_tree, 
                            compute_stability, 
                            get_clusters)
                           
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
    min_points = min(dim - 1, min_points)
    try:
        core_distances = np.partition(distance_matrix, 
                                      min_points, 
                                      axis=0)[min_points]
    except AttributeError:
        core_distances = np.sort(distance_matrix,
                                 axis=0)[min_points]        
                                  
    stage1 = np.where(core_distances > distance_matrix, 
                      core_distances, distance_matrix)
    result = np.where(core_distances > stage1.T,
                      core_distances.T, stage1.T).T
    return result.astype(np.double)

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

    if type(min_samples) is not int or type(min_cluster_size) is not int:
        raise ValueError('Min samples and min cluster size must be integers!')

    if min_samples <= 0 or min_cluster_size <= 0:
        raise ValueError('Min samples and Min cluster size must be positive integers')
    
    X = check_array(X, accept_sparse='csr')
    
    if metric == 'minkowski':
        if p is None:
            raise TypeError('Minkowski metric given but no p value supplied!')
        if p < 0:
            raise ValueError('Minkowski metric with negative p value is not defined!')
        distance_matrix = pairwise_distances(X, metric=metric, p=p)
    else:
        distance_matrix = pairwise_distances(X, metric=metric)
        
    mutual_reachability_graph = mutual_reachability(distance_matrix,
                                                    min_samples)
    raw_tree = single_linkage(mutual_reachability_graph)
    condensed_tree, new_points = condense_tree(raw_tree, get_points(raw_tree), min_cluster_size)
    stability_dict = compute_stability(condensed_tree)
    cluster_list = get_clusters(condensed_tree, stability_dict, new_points)
    
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
        
    def fit_predict(self, X, y=None):
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
        self.fit(X)
        return self.labels_
