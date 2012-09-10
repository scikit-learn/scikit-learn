"""
The :mod:`sklearn.metrics` module includes score functions, performance metrics
and pairwise metrics and distance computations.
"""

from .metrics import confusion_matrix, roc_curve, auc, precision_score, \
                recall_score, fbeta_score, f1_score, zero_one_score, \
                precision_recall_fscore_support, classification_report, \
                precision_recall_curve, explained_variance_score, r2_score, \
                zero_one, mean_square_error, hinge_loss, matthews_corrcoef, \
                mean_squared_error, average_precision_score, auc_score

from . import cluster
from .cluster import adjusted_rand_score
from .cluster import homogeneity_completeness_v_measure
from .cluster import homogeneity_score
from .cluster import completeness_score
from .cluster import v_measure_score
from .cluster import silhouette_score
from .cluster import mutual_info_score
from .cluster import adjusted_mutual_info_score
from .cluster import normalized_mutual_info_score
from .pairwise import euclidean_distances, pairwise_distances, pairwise_kernels

__all__ = ['adjusted_mutual_info_score',
           'adjusted_rand_score',
           'auc',
           'auc_score',
           'average_precision_score',
           'classification_report',
           'cluster',
           'completeness_score',
           'confusion_matrix',
           'euclidean_distances',
           'explained_variance_score',
           'f1_score',
           'fbeta_score',
           'hinge_loss',
           'homogeneity_completeness_v_measure',
           'homogeneity_score',
           'matthews_corrcoef',
           'mean_square_error',
           'mean_squared_error',
           'mutual_info_score',
           'normalized_mutual_info_score',
           'pairwise_distances',
           'pairwise_kernels',
           'precision_recall_curve',
           'precision_recall_fscore_support',
           'precision_score',
           'r2_score',
           'recall_score',
           'roc_curve',
           'silhouette_score',
           'v_measure_score',
           'zero_one',
           'zero_one_score']
