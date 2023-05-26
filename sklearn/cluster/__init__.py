"""
The :mod:`sklearn.cluster` module gathers popular unsupervised clustering
algorithms.
"""
from ..externals import _lazy_loader

__getattr__, __dir__, __all__ = _lazy_loader.attach(
    __name__,
    submod_attrs={
        "_affinity_propagation": ["AffinityPropagation", "affinity_propagation"],
        "_agglomerative": [
            "linkage_tree",
            "AgglomerativeClustering",
            "FeatureAgglomeration",
            "ward_tree",
        ],
        "_bicluster": ["SpectralBiclustering", "SpectralCoclustering"],
        "_birch": ["Birch"],
        "_bisect_k_means": ["BisectingKMeans"],
        "_dbscan": ["dbscan", "DBSCAN"],
        "_kmeans": ["KMeans", "MiniBatchKMeans", "k_means", "kmeans_plusplus"],
        "_mean_shift": [
            "estimate_bandwidth",
            "get_bin_seeds",
            "MeanShift",
            "mean_shift",
        ],
        "_optics": [
            "OPTICS",
            "cluster_optics_dbscan",
            "compute_optics_graph",
            "cluster_optics_xi",
        ],
        "_spectral": ["spectral_clustering", "SpectralClustering"],
    },
)
