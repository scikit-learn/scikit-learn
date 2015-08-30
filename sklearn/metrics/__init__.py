"""
The :mod:`sklearn.metrics` module includes score functions, performance metrics
and pairwise metrics and distance computations.
"""


from .ranking import auc
from .ranking import average_precision_score
from .ranking import coverage_error
from .ranking import label_ranking_average_precision_score
from .ranking import label_ranking_loss
from .ranking import precision_recall_curve
from .ranking import roc_auc_score
from .ranking import roc_curve

from .classification import accuracy_score
from .classification import classification_report
from .classification import cohen_kappa_score
from .classification import confusion_matrix
from .classification import f1_score
from .classification import fbeta_score
from .classification import hamming_loss
from .classification import hinge_loss
from .classification import jaccard_similarity_score
from .classification import log_loss
from .classification import matthews_corrcoef
from .classification import precision_recall_fscore_support
from .classification import precision_score
from .classification import recall_score
from .classification import zero_one_loss
from .classification import brier_score_loss

from . import cluster
from .cluster import adjusted_mutual_info_score
from .cluster import adjusted_rand_score
from .cluster import completeness_score
from .cluster import consensus_score
from .cluster import homogeneity_completeness_v_measure
from .cluster import homogeneity_score
from .cluster import mutual_info_score
from .cluster import normalized_mutual_info_score
from .cluster import silhouette_samples
from .cluster import silhouette_score
from .cluster import v_measure_score

from .pairwise import euclidean_distances
from .pairwise import pairwise_distances
from .pairwise import pairwise_distances_argmin
from .pairwise import pairwise_distances_argmin_min
from .pairwise import pairwise_kernels

from .regression import explained_variance_score
from .regression import mean_absolute_error
from .regression import mean_squared_error
from .regression import median_absolute_error
from .regression import r2_score

from .scorer import make_scorer
from .scorer import SCORERS
from .scorer import get_scorer

__all__ = [
    'accuracy_score',
    'adjusted_mutual_info_score',
    'adjusted_rand_score',
    'auc',
    'average_precision_score',
    'classification_report',
    'cluster',
    'completeness_score',
    'confusion_matrix',
    'consensus_score',
    'coverage_error',
    'euclidean_distances',
    'explained_variance_score',
    'f1_score',
    'fbeta_score',
    'get_scorer',
    'hamming_loss',
    'hinge_loss',
    'homogeneity_completeness_v_measure',
    'homogeneity_score',
    'jaccard_similarity_score',
    'label_ranking_average_precision_score',
    'label_ranking_loss',
    'log_loss',
    'make_scorer',
    'matthews_corrcoef',
    'mean_absolute_error',
    'mean_squared_error',
    'median_absolute_error',
    'mutual_info_score',
    'normalized_mutual_info_score',
    'pairwise_distances',
    'pairwise_distances_argmin',
    'pairwise_distances_argmin_min',
    'pairwise_distances_argmin_min',
    'pairwise_kernels',
    'precision_recall_curve',
    'precision_recall_fscore_support',
    'precision_score',
    'r2_score',
    'recall_score',
    'roc_auc_score',
    'roc_curve',
    'SCORERS',
    'silhouette_samples',
    'silhouette_score',
    'v_measure_score',
    'zero_one_loss',
    'brier_score_loss',
]
