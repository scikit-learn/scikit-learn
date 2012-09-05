"""
The :mod:`sklearn.linear_model` module implements genelarized linear models. It
includes Ridge regression, Bayesian Regression, Lasso and Elastic Net
estimators computed with Least Angle Regression and coordinate descent. It also
implements Stochastic Gradient Descent related algorithms.
"""

# See http://scikit-learn.sourceforge.net/modules/sgd.html and
# http://scikit-learn.sourceforge.net/modules/linear_model.html for
# complete documentation.

from .base import LinearRegression

from .bayes import BayesianRidge, ARDRegression
from .least_angle import Lars, LassoLars, lars_path, LarsCV, LassoLarsCV, \
                         LassoLarsIC
from .coordinate_descent import Lasso, ElasticNet, LassoCV, ElasticNetCV, \
                                lasso_path, enet_path, MultiTaskLasso, \
                                MultiTaskElasticNet
from .sgd_fast import Hinge, Log, ModifiedHuber, SquaredLoss, Huber
from .stochastic_gradient import SGDClassifier, SGDRegressor
from .ridge import Ridge, RidgeCV, RidgeClassifier, RidgeClassifierCV, \
                   ridge_regression
from .logistic import LogisticRegression
from .omp import orthogonal_mp, orthogonal_mp_gram, OrthogonalMatchingPursuit
from .perceptron import Perceptron
from .randomized_l1 import RandomizedLasso, RandomizedLogisticRegression, \
                           lasso_stability_path

from . import sparse

__all__ = ['ARDRegression',
           'BayesianRidge',
           'ElasticNet',
           'ElasticNetCV',
           'Hinge',
           'Huber',
           'LARS',
           'Lars',
           'LarsCV',
           'Lasso',
           'LassoCV',
           'LassoLARS',
           'LassoLars',
           'LassoLarsCV',
           'LassoLarsIC',
           'LinearRegression',
           'Log',
           'LogisticRegression',
           'ModifiedHuber',
           'MultiTaskElasticNet',
           'MultiTaskLasso',
           'OrthogonalMatchingPursuit',
           'Perceptron',
           'RandomizedLasso',
           'RandomizedLogisticRegression',
           'Ridge',
           'RidgeCV',
           'RidgeClassifier',
           'RidgeClassifierCV',
           'SGDClassifier',
           'SGDRegressor',
           'SquaredLoss',
           'enet_path',
           'lars_path',
           'lasso_path',
           'lasso_stability_path',
           'orthogonal_mp',
           'orthogonal_mp_gram',
           'ridge_regression',
           'sparse']
