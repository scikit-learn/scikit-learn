"""This module implements the 'fast' gradient boosting estimators.

The implementation is a port from pygbm which is itself strongly inspired
from LightGBM.
"""
from .gradient_boosting import HistGradientBoostingClassifier
from .gradient_boosting import HistGradientBoostingRegressor

__all__ = ["HistGradientBoostingClassifier", "HistGradientBoostingRegressor"]
