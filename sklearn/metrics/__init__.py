"""
The :mod:`sklearn.metrics` module includes score functions, performance metrics
and pairwise metrics and distance computations.
"""

from .metrics import confusion_matrix, roc_curve, auc, precision_score, \
                recall_score, fbeta_score, f1_score, zero_one_score, \
                precision_recall_fscore_support, classification_report, \
                precision_recall_curve, explained_variance_score, r2_score, \
                zero_one, mean_square_error, hinge_loss

from . import cluster
from .cluster import adjusted_rand_score
from .cluster import homogeneity_completeness_v_measure
from .cluster import homogeneity_score
from .cluster import completeness_score
from .cluster import v_measure_score
from .cluster import silhouette_score
from .cluster import mutual_info_score
from .cluster import adjusted_mutual_info_score
from .pairwise import euclidean_distances, pairwise_distances, pairwise_kernels
