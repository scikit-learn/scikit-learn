"""
The :mod:`sklearn.metrics` module includes score functions, performance metrics
and pairwise metrics and distance computations.
"""

from .metrics import (accuracy_score,
                      average_precision_score,
                      auc,
                      roc_auc_score,
                      classification_report,
                      confusion_matrix,
                      explained_variance_score,
                      f1_score,
                      fbeta_score,
                      hamming_loss,
                      hinge_loss,
                      jaccard_similarity_score,
                      log_loss,
                      matthews_corrcoef,
                      mean_squared_error,
                      mean_absolute_error,
                      precision_recall_curve,
                      precision_recall_fscore_support,
                      precision_score,
                      recall_score,
                      r2_score,
                      roc_curve,
                      zero_one_loss)


# Deprecated in 0.16
from .metrics import auc_score

from .scorer import make_scorer, SCORERS

from . import cluster
from .cluster import (adjusted_rand_score,
                      adjusted_mutual_info_score,
                      completeness_score,
                      homogeneity_completeness_v_measure,
                      homogeneity_score,
                      mutual_info_score,
                      normalized_mutual_info_score,
                      silhouette_score,
                      silhouette_samples,
                      v_measure_score,
                      consensus_score)

from .pairwise import (euclidean_distances,
                       pairwise_distances,
                       pairwise_distances_argmin_min,
                       pairwise_distances_argmin,
                       pairwise_kernels)

__all__ = ['accuracy_score',
           'adjusted_mutual_info_score',
           'adjusted_rand_score',
           'auc',
           'roc_auc_score',
           'average_precision_score',
           'classification_report',
           'cluster',
           'completeness_score',
           'confusion_matrix',
           'euclidean_distances',
           'pairwise_distances_argmin_min',
           'explained_variance_score',
           'f1_score',
           'fbeta_score',
           'hamming_loss',
           'hinge_loss',
           'homogeneity_completeness_v_measure',
           'homogeneity_score',
           'jaccard_similarity_score',
           'log_loss',
           'matthews_corrcoef',
           'mean_squared_error',
           'mean_absolute_error',
           'mutual_info_score',
           'normalized_mutual_info_score',
           'pairwise_distances',
           'pairwise_distances_argmin',
           'pairwise_distances_argmin_min',
           'pairwise_kernels',
           'precision_recall_curve',
           'precision_recall_fscore_support',
           'precision_score',
           'r2_score',
           'recall_score',
           'roc_curve',
           'silhouette_score',
           'silhouette_samples',
           'v_measure_score',
           'consensus_score',
           'zero_one_loss',
           'make_scorer',
           'SCORERS']
