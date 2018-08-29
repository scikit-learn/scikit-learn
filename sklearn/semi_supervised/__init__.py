"""
The :mod:`sklearn.semi_supervised` module implements semi-supervised learning
algorithms. These algorithms utilize small amounts of labeled data and large
amounts of unlabeled data for classification tasks. This module includes Label
Propagation.
"""

from .label_propagation import LabelPropagation, LabelSpreading
from .self_training import SelfTraining

__all__ = ['SelfTraining', 'LabelPropagation', 'LabelSpreading']
