"""
:mod:`sklearn.linear_model.sparse` is the sparse counterpart
of :mod:`sklearn.linear_model`.

"""

from .coordinate_descent import Lasso, ElasticNet
from .logistic import LogisticRegression
from .stochastic_gradient import SGDClassifier, SGDRegressor
