"""
The :mod:`sklearn.metrics.cluster` submodule contains evaluation metrics for
cluster analysis results. There are two forms of evaluation:

- supervised, which uses a ground truth class values for each sample.
- unsupervised, which does not and measures the 'quality' of the model itself.
"""
from .supervised import adjusted_mutual_info_score
from .supervised import normalized_mutual_info_score
from .supervised import adjusted_rand_score
from .supervised import completeness_score
from .supervised import contingency_matrix
from .supervised import expected_mutual_information
from .supervised import homogeneity_completeness_v_measure
from .supervised import homogeneity_score
from .supervised import mutual_info_score
from .supervised import v_measure_score
from .supervised import fowlkes_mallows_score
from .supervised import entropy
from .unsupervised import silhouette_samples
from .unsupervised import silhouette_score
from .unsupervised import calinski_harabasz_score
from .unsupervised import calinski_harabaz_score
from .unsupervised import davies_bouldin_score
from .bicluster import consensus_score

__all__ = ["adjusted_mutual_info_score", "normalized_mutual_info_score",
           "adjusted_rand_score", "completeness_score", "contingency_matrix",
           "expected_mutual_information", "homogeneity_completeness_v_measure",
           "homogeneity_score", "mutual_info_score", "v_measure_score",
           "fowlkes_mallows_score", "entropy", "silhouette_samples",
           "silhouette_score", "calinski_harabaz_score",
           "calinski_harabasz_score", "davies_bouldin_score",
           "consensus_score"]
