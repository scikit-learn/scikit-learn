"""Popular unsupervised clustering algorithms."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from sklearn.cluster._affinity_propagation import (
    AffinityPropagation,
    affinity_propagation,
)
from sklearn.cluster._agglomerative import (
    AgglomerativeClustering,
    FeatureAgglomeration,
    linkage_tree,
    ward_tree,
)
from sklearn.cluster._bicluster import SpectralBiclustering, SpectralCoclustering
from sklearn.cluster._birch import Birch
from sklearn.cluster._bisect_k_means import BisectingKMeans
from sklearn.cluster._dbscan import DBSCAN, dbscan
from sklearn.cluster._hdbscan.hdbscan import HDBSCAN
from sklearn.cluster._kmeans import KMeans, MiniBatchKMeans, k_means, kmeans_plusplus
from sklearn.cluster._mean_shift import (
    MeanShift,
    estimate_bandwidth,
    get_bin_seeds,
    mean_shift,
)
from sklearn.cluster._optics import (
    OPTICS,
    cluster_optics_dbscan,
    cluster_optics_xi,
    compute_optics_graph,
)
from sklearn.cluster._spectral import SpectralClustering, spectral_clustering

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
