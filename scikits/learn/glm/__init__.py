"""
Generalized linear models
=========================

scikits.learn.glm is a module to fit genelarized linear models.
It includes Ridge regression, Bayesian Regression, Lasso and
Elastic Net estimators computed with Least Angle Regression
and coordinate descent.

"""

from .lars import LARS, LassoLARS, lars_path, LeastAngleRegression
from .coordinate_descent import Lasso, ElasticNet, LassoCV, ElasticNetCV
from .bayes import Ridge, BayesianRidge, ARDRegression

