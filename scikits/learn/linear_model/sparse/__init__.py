"""
Generalized linear models with sparse data
==========================================

scikits.learn.glm.sparse is the sparse counterpart
of scikits.learn.glm

"""

from .coordinate_descent import Lasso, ElasticNet
from .logistic import LogisticRegression
from .stochastic_gradient import SGDClassifier, SGDRegressor
from .ridge import Ridge, RidgeClassifier
