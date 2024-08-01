"""Semi-supervised learning algorithms.

These algorithms utilize small amounts of labeled data and large amounts of unlabeled
data for classification tasks.
"""

from ._label_propagation import LabelPropagation, LabelSpreading
from ._self_training import SelfTrainingClassifier

__all__ = ["SelfTrainingClassifier", "LabelPropagation", "LabelSpreading"]
