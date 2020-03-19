"""
The :mod:`sklearn.ensemble` module includes ensemble-based methods for
classification, regression and anomaly detection.
"""
import typing

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
from ._stacking import StackingClassifier
from ._stacking import StackingRegressor

if typing.TYPE_CHECKING:
    # Workaround for type checkers (e.g. mypy) to avoid
    # import errors for experimenal estimators.
    # TODO: remove the above check once the estimator is no longer
    #       experimental.
    from ._hist_gradient_boosting.gradient_boosting import (
        HistGradientBoostingRegressor, HistGradientBoostingClassifier
    )

__all__ = ["BaseEnsemble",
           "RandomForestClassifier", "RandomForestRegressor",
           "RandomTreesEmbedding", "ExtraTreesClassifier",
           "ExtraTreesRegressor", "BaggingClassifier",
           "BaggingRegressor", "IsolationForest", "GradientBoostingClassifier",
           "GradientBoostingRegressor", "AdaBoostClassifier",
           "AdaBoostRegressor", "VotingClassifier", "VotingRegressor",
           "StackingClassifier", "StackingRegressor",
           ]
