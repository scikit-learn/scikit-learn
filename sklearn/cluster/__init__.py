"""
The :mod:`sklearn.cluster` module gathers popular unsupervised clustering
algorithms.
"""

from .spectral import spectral_clustering, SpectralClustering
from .mean_shift_ import (mean_shift, MeanShift,
                          estimate_bandwidth, get_bin_seeds)
from .affinity_propagation_ import affinity_propagation, AffinityPropagation
from .hierarchical import (ward_tree, AgglomerativeClustering, linkage_tree,
                           FeatureAgglomeration)
from .k_means_ import k_means, KMeans, MiniBatchKMeans
from .dbscan_ import dbscan, DBSCAN
from .bicluster import SpectralBiclustering, SpectralCoclustering
from .birch import Birch
from .som_ import SelfOrganizingMap

import numpy as np

__all__ = ['AffinityPropagation',
           'AgglomerativeClustering',
           'Birch',
           'DBSCAN',
           'KMeans',
           'FeatureAgglomeration',
           'MeanShift',
           'MiniBatchKMeans',
           'SpectralClustering',
           'affinity_propagation',
           'dbscan',
           'estimate_bandwidth',
           'get_bin_seeds',
           'k_means',
           'linkage_tree',
           'mean_shift',
           'spectral_clustering',
           'ward_tree',
           'SpectralBiclustering',
           'SpectralCoclustering']

def calinski_index(X,labels,centroids):
    mean = np.mean(X,axis=0) 
    B = np.sum([ (c - mean)**2 for c in centroids])
    W = np.sum([ (x-centroids[labels[i]])**2 
                 for i,x in enumerate(X)])
    c = len(centroids)
    n = len(X)
    return (B /(c-1))/(W/ (n-c))
    
