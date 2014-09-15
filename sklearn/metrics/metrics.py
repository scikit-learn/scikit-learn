import warnings
warnings.warn("sklearn.metrics.metrics is deprecated and will be removed in "
              "0.18. Please import from sklearn.metrics",
              DeprecationWarning)


from .ranking import auc
from .ranking import average_precision_score
from .ranking import label_ranking_average_precision_score
from .ranking import precision_recall_curve
from .ranking import roc_auc_score
from .ranking import roc_curve

from .classification import accuracy_score
from .classification import classification_report
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

from .regression import explained_variance_score
from .regression import mean_absolute_error
from .regression import mean_squared_error
from .regression import r2_score

# Deprecated in 0.16
from .ranking import auc_score
