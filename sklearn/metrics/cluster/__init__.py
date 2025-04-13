"""Evaluation metrics for cluster analysis results.

- Supervised evaluation uses a ground truth class values for each sample.
- Unsupervised evaluation does not use ground truths and measures the "quality" of the
  model itself.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ._bicluster import consensus_score
from ._supervised import (
    adjusted_mutual_info_score,
    adjusted_rand_score,
    completeness_score,
    contingency_matrix,
    entropy,
    expected_mutual_information,
    fowlkes_mallows_score,
    homogeneity_completeness_v_measure,
    homogeneity_score,
    mutual_info_score,
    normalized_mutual_info_score,
    pair_confusion_matrix,
    rand_score,
    v_measure_score,
)
from ._unsupervised import (
    calinski_harabasz_score,
    davies_bouldin_score,
    silhouette_samples,
    silhouette_score,
)

__all__ = [
    "adjusted_mutual_info_score",
    "adjusted_rand_score",
    "calinski_harabasz_score",
    "completeness_score",
    "consensus_score",
    "contingency_matrix",
    "davies_bouldin_score",
    "entropy",
    "expected_mutual_information",
    "fowlkes_mallows_score",
    "homogeneity_completeness_v_measure",
    "homogeneity_score",
    "mutual_info_score",
    "normalized_mutual_info_score",
    "pair_confusion_matrix",
    "rand_score",
    "silhouette_samples",
    "silhouette_score",
    "v_measure_score",
]
