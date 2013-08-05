"""
The :mod:`sklearn.ensemble` module includes ensemble-based methods for
classification and regression.
"""

from .base import BaseEnsemble
from .forest import RandomForestClassifier
from .forest import RandomForestRegressor
from .forest import RandomTreesEmbedding
from .forest import ExtraTreesClassifier
from .forest import ExtraTreesRegressor
from .weight_boosting import AdaBoostClassifier
from .weight_boosting import AdaBoostRegressor
from .gradient_boosting import GradientBoostingClassifier
from .gradient_boosting import GradientBoostingRegressor

from . import forest
from . import weight_boosting
from . import gradient_boosting
from . import partial_dependence

__all__ = ["BaseEnsemble", "RandomForestClassifier", "RandomForestRegressor",
           "RandomTreesEmbedding", "ExtraTreesClassifier",
           "ExtraTreesRegressor", "GradientBoostingClassifier",
           "GradientBoostingRegressor", "AdaBoostClassifier",
           "AdaBoostRegressor", "forest", "gradient_boosting",
           "partial_dependence", "weight_boosting"]
