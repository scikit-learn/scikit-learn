"""
The :mod:`sklearn.linear_model.sparse` submodule is the sparse counterpart of
the :mod:`sklearn.linear_model` module.
"""

from .coordinate_descent import Lasso, ElasticNet
from .logistic import LogisticRegression
from .stochastic_gradient import SGDClassifier, SGDRegressor
