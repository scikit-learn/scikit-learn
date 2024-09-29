"""Semi-supervised learning algorithms.

These algorithms utilize small amounts of labeled data and large amounts of unlabeled
data for classification tasks.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ._label_propagation import LabelPropagation, LabelSpreading
from ._self_training import SelfTrainingClassifier

__all__ = ["SelfTrainingClassifier", "LabelPropagation", "LabelSpreading"]
