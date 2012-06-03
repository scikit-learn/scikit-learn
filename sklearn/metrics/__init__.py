"""
The :mod:`sklearn.metrics` module includes score functions, performance metrics
and pairwise metrics and distance computations.
"""

from .metrics import confusion_matrix
from .metrics import roc_curve
from .metrics import auc
from .metrics import precision_score
from .metrics import recall_score
from .metrics import fbeta_score
from .metrics import f1_score
from .metrics import zero_one_score
from .metrics import precision_recall_fscore_support
from .metrics import classification_report
from .metrics import precision_recall_curve
from .metrics import explained_variance_score
from .metrics import r2_score
from .metrics import zero_one
from .metrics import mean_square_error
from .metrics import hinge_loss
from .metrics import matthews_corrcoef
from .metrics import mean_squared_error
from .metrics import brier_score
from .metrics import calibration_plot

from . import cluster
from .cluster import adjusted_rand_score
from .cluster import homogeneity_completeness_v_measure
from .cluster import homogeneity_score
from .cluster import completeness_score
from .cluster import v_measure_score
from .cluster import silhouette_score
from .cluster import mutual_info_score
from .cluster import adjusted_mutual_info_score

from .pairwise import euclidean_distances
from .pairwise import pairwise_distances
from .pairwise import pairwise_kernels
