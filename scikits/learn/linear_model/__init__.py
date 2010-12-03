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

from .bayes import BayesianRidge, ARDRegression
from .base import LinearRegression
from .bayes import BayesianRidge, ARDRegression
from .least_angle import LARS, LassoLARS, lars_path
from .coordinate_descent import Lasso, ElasticNet, LassoCV, ElasticNetCV, \
                                lasso_path, enet_path
from .ridge import Ridge
from .logistic import LogisticRegression

from . import sparse
from .stochastic_gradient import SGDClassifier, SGDRegressor
from .base import Log, ModifiedHuber, Hinge, SquaredLoss, Huber
