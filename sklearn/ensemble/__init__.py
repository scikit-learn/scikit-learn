"""
The :mod:`sklearn.ensemble` module includes ensemble-based methods for
classification and regression.
"""

from .base import BaseEnsemble
from .bagging import BaggedClassifier
from .bagging import BaggedRegressor
from .forest import RandomForestClassifier
from .forest import RandomForestRegressor
from .forest import ExtraTreesClassifier
from .forest import ExtraTreesRegressor
