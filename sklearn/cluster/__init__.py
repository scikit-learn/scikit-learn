"""
The :mod:`sklearn.cluster` module gathers popular unsupervised clustering
algorithms.
"""

from .spectral import spectral_clustering, SpectralClustering
from .mean_shift_ import mean_shift, MeanShift, estimate_bandwidth, \
     get_bin_seeds
from .affinity_propagation_ import affinity_propagation, AffinityPropagation
from .hierarchical import ward_tree, Ward, WardAgglomeration
from .k_means_ import k_means, KMeans, MiniBatchKMeans
from .dbscan_ import dbscan, DBSCAN

__all__ = ['AffinityPropagation',
           'DBSCAN',
           'KMeans',
           'MeanShift',
           'MiniBatchKMeans',
           'SpectralClustering',
           'Ward',
           'WardAgglomeration',
           'affinity_propagation',
           'dbscan',
           'estimate_bandwidth',
           'get_bin_seeds',
           'k_means',
           'mean_shift',
           'spectral_clustering',
           'ward_tree']
