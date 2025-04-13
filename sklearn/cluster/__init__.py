"""Popular unsupervised clustering algorithms."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ._affinity_propagation import AffinityPropagation, affinity_propagation
from ._agglomerative import (
    AgglomerativeClustering,
    FeatureAgglomeration,
    linkage_tree,
    ward_tree,
)
from ._bicluster import SpectralBiclustering, SpectralCoclustering
from ._birch import Birch
from ._bisect_k_means import BisectingKMeans
from ._dbscan import DBSCAN, dbscan
from ._hdbscan.hdbscan import HDBSCAN
from ._kmeans import KMeans, MiniBatchKMeans, k_means, kmeans_plusplus
from ._mean_shift import MeanShift, estimate_bandwidth, get_bin_seeds, mean_shift
from ._optics import (
    OPTICS,
    cluster_optics_dbscan,
    cluster_optics_xi,
    compute_optics_graph,
)
from ._spectral import SpectralClustering, spectral_clustering

__all__ = [
    "DBSCAN",
    "HDBSCAN",
    "OPTICS",
    "AffinityPropagation",
    "AgglomerativeClustering",
    "Birch",
    "BisectingKMeans",
    "FeatureAgglomeration",
    "KMeans",
    "MeanShift",
    "MiniBatchKMeans",
    "SpectralBiclustering",
    "SpectralClustering",
    "SpectralCoclustering",
    "affinity_propagation",
    "cluster_optics_dbscan",
    "cluster_optics_xi",
    "compute_optics_graph",
    "dbscan",
    "estimate_bandwidth",
    "get_bin_seeds",
    "k_means",
    "kmeans_plusplus",
    "linkage_tree",
    "mean_shift",
    "spectral_clustering",
    "ward_tree",
]
