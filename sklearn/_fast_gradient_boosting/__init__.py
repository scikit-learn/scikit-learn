"""This module implements the 'fast' gradient boosting estimators.

The implementation is a port from pygbm which is itself strongly inspired
from LightGBM.
"""
from .gradient_boosting import FastGradientBoostingClassifier
from .gradient_boosting import FastGradientBoostingRegressor

__all__ = ["FastGradientBoostingClassifier", "FastGradientBoostingRegressor"]
