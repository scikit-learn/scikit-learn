"""
Clustering algorithms
"""
from .spectral import spectral_clustering, SpectralClustering
from .mean_shift_ import mean_shift, MeanShift, estimate_bandwidth
from .affinity_propagation_ import affinity_propagation, AffinityPropagation
from .k_means_ import k_means, KMeans
from .hierarchical import ward_tree, Ward, WardAgglomeration

