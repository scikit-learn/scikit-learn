"""
The :mod:`sklearn.ensemble` module includes ensemble-based methods for
classification, regression and anomaly detection.
"""

from ._base import BaseEnsemble
from ._forest import RandomForestClassifier
from ._forest import RandomForestRegressor
from ._forest import RandomTreesEmbedding
from ._forest import ExtraTreesClassifier
from ._forest import ExtraTreesRegressor
from ._bagging import BaggingClassifier
from ._bagging import BaggingRegressor
from ._iforest import IsolationForest
from ._weight_boosting import AdaBoostClassifier
from ._weight_boosting import AdaBoostRegressor
from ._gb import GradientBoostingClassifier
from ._gb import GradientBoostingRegressor
from ._voting import VotingClassifier
from ._voting import VotingRegressor
from ._partial_dependence import partial_dependence
from ._partial_dependence import plot_partial_dependence

__all__ = ["BaseEnsemble",
           "RandomForestClassifier", "RandomForestRegressor",
           "RandomTreesEmbedding", "ExtraTreesClassifier",
           "ExtraTreesRegressor", "BaggingClassifier",
           "BaggingRegressor", "IsolationForest", "GradientBoostingClassifier",
           "GradientBoostingRegressor", "AdaBoostClassifier",
           "AdaBoostRegressor", "VotingClassifier", "VotingRegressor",
           "bagging", "forest", "gradient_boosting",
           "partial_dependence", "plot_partial_dependence", "weight_boosting"]
