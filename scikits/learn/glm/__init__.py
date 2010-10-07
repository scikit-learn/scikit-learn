"""
Generalized linear models
=========================

scikits.learn.glm is a module to fit genelarized linear models.
It includes Ridge regression, Bayesian Regression, Lasso and
Elastic Net estimators computed with Least Angle Regression
and coordinate descent.

"""

from .bayes import BayesianRidge, ARDRegression
from .base import LinearRegression
from .bayes import BayesianRidge, ARDRegression
from .lars import LARS, LassoLARS, lars_path
from .coordinate_descent import Lasso, ElasticNet, LassoCV, ElasticNetCV, \
                                lasso_path, enet_path
from .ridge import Ridge

