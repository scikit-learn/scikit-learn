"""
The :mod:`sklearn.linear_model` module implements generalized linear models. It
includes Ridge regression, Bayesian Regression, Lasso and Elastic Net
estimators computed with Least Angle Regression and coordinate descent. It also
implements Stochastic Gradient Descent related algorithms.
"""

# See http://scikit-learn.sourceforge.net/modules/sgd.html and
# http://scikit-learn.sourceforge.net/modules/linear_model.html for
# complete documentation.

from .base import LinearRegression

from .bayes import BayesianRidge, ARDRegression
from .least_angle import (Lars, LassoLars, lars_path, lars_path_gram, LarsCV,
                          LassoLarsCV, LassoLarsIC)
from .coordinate_descent import (Lasso, ElasticNet, LassoCV, ElasticNetCV,
                                 lasso_path, enet_path, MultiTaskLasso,
                                 MultiTaskElasticNet, MultiTaskElasticNetCV,
                                 MultiTaskLassoCV)
from .huber import HuberRegressor
from .sgd_fast import Hinge, Log, ModifiedHuber, SquaredLoss, Huber
from .stochastic_gradient import SGDClassifier, SGDRegressor
from .ridge import (Ridge, RidgeCV, RidgeClassifier, RidgeClassifierCV,
                    ridge_regression)
from .logistic import (LogisticRegression, LogisticRegressionCV,
                       logistic_regression_path)
from .omp import (orthogonal_mp, orthogonal_mp_gram, OrthogonalMatchingPursuit,
                  OrthogonalMatchingPursuitCV)
from .passive_aggressive import PassiveAggressiveClassifier
from .passive_aggressive import PassiveAggressiveRegressor
from .perceptron import Perceptron

from .ransac import RANSACRegressor
from .theil_sen import TheilSenRegressor

__all__ = ['ARDRegression',
           'BayesianRidge',
           'ElasticNet',
           'ElasticNetCV',
           'Hinge',
           'Huber',
           'HuberRegressor',
           'Lars',
           'LarsCV',
           'Lasso',
           'LassoCV',
           'LassoLars',
           'LassoLarsCV',
           'LassoLarsIC',
           'LinearRegression',
           'Log',
           'LogisticRegression',
           'LogisticRegressionCV',
           'ModifiedHuber',
           'MultiTaskElasticNet',
           'MultiTaskElasticNetCV',
           'MultiTaskLasso',
           'MultiTaskLassoCV',
           'OrthogonalMatchingPursuit',
           'OrthogonalMatchingPursuitCV',
           'PassiveAggressiveClassifier',
           'PassiveAggressiveRegressor',
           'Perceptron',
           'Ridge',
           'RidgeCV',
           'RidgeClassifier',
           'RidgeClassifierCV',
           'SGDClassifier',
           'SGDRegressor',
           'SquaredLoss',
           'TheilSenRegressor',
           'enet_path',
           'lars_path',
           'lars_path_gram',
           'lasso_path',
           'logistic_regression_path',
           'orthogonal_mp',
           'orthogonal_mp_gram',
           'ridge_regression',
           'RANSACRegressor']
