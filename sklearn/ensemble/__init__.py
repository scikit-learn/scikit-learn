"""Ensemble-based methods for classification, regression and anomaly detection."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ._bagging import BaggingClassifier, BaggingRegressor
from ._base import BaseEnsemble
from ._bagging import BaggingClassifier
from ._bagging import BaggingRegressor
from ._iforest import IsolationForest
from ._weight_boosting import AdaBoostClassifier
from ._weight_boosting import AdaBoostRegressor
from ._voting import VotingClassifier
from ._voting import VotingRegressor
from ._stacking import StackingClassifier
from ._stacking import StackingRegressor
from ._ensemble_selection import EnsembleSelection
from ._forest import (
    ExtraTreesClassifier,
    ExtraTreesRegressor,
    RandomForestClassifier,
    RandomForestRegressor,
    RandomTreesEmbedding,
)
from ._gb import GradientBoostingClassifier, GradientBoostingRegressor

from ._hist_gradient_boosting.gradient_boosting import (
    HistGradientBoostingClassifier,
    HistGradientBoostingRegressor,
)
from ._iforest import IsolationForest
from ._stacking import StackingClassifier, StackingRegressor
from ._voting import VotingClassifier, VotingRegressor
from ._weight_boosting import AdaBoostClassifier, AdaBoostRegressor

__all__ = [
    "BaseEnsemble",
    "RandomForestClassifier",
    "RandomForestRegressor",
    "RandomTreesEmbedding",
    "ExtraTreesClassifier",
    "ExtraTreesRegressor",
    "BaggingClassifier",
    "BaggingRegressor",
    "IsolationForest",
    "GradientBoostingClassifier",
    "GradientBoostingRegressor",
    "AdaBoostClassifier",
    "AdaBoostRegressor",
    "VotingClassifier",
    "VotingRegressor",
    "StackingClassifier",
    "StackingRegressor",
    "EnsembleSelection",
    "HistGradientBoostingClassifier",
    "HistGradientBoostingRegressor",
]
