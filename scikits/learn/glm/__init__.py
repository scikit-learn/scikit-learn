"""
Generalized linear models
=========================

scikits.learn.glm is a module to fit genelarized linear models.
It includes Ridge regression, Bayesian Regression, Lasso and
Elastic Net estimators computed with Least Angle Regression
and coordinate descent.

"""

from .base import LinearRegression
from .lars import LARS, LassoLARS, lars_path, LeastAngleRegression
from .coordinate_descent import Lasso, ElasticNet, LassoCV, ElasticNetCV, \
                                lasso_path, enet_path
from .bayes import BayesianRidge, ARDRegression
from .bayes_classification import ARDClassification
from .ridge import Ridge

