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
from .optics_ import OPTICS, cluster_optics_dbscan, compute_optics_graph
from .bicluster import SpectralBiclustering, SpectralCoclustering
from .birch import Birch

__all__ = ['AffinityPropagation',
           'AgglomerativeClustering',
           'Birch',
           'DBSCAN',
           'OPTICS',
           'cluster_optics_dbscan',
           'compute_optics_graph',
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
