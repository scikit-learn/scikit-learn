"""
Generalized linear models
=========================

scikits.learn.linear_model is a module to fit genelarized linear
models.  It includes Ridge regression, Bayesian Regression, Lasso and
Elastic Net estimators computed with Least Angle Regression and
coordinate descent.

It also implements Stochastic Gradient Descent related algorithms.

See http://scikit-learn.sourceforge.net/modules/sgd.html and
http://scikit-learn.sourceforge.net/modules/linear_model.html for
complete documentation.
"""

from .base import LinearRegression, Log, ModifiedHuber, Hinge, \
     SquaredLoss, Huber

from .bayes import BayesianRidge, ARDRegression
from .least_angle import LARS, LassoLARS, lars_path
from .coordinate_descent import Lasso, ElasticNet, LassoCV, ElasticNetCV, \
                                lasso_path, enet_path
from .stochastic_gradient import SGDClassifier, SGDRegressor
from .ridge import Ridge, RidgeCV, RidgeClassifier, RidgeClassifierCV
from .logistic import LogisticRegression

from . import sparse

