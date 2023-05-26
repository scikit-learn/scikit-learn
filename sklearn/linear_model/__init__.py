"""
The :mod:`sklearn.linear_model` module implements a variety of linear models.
"""

# See http://scikit-learn.sourceforge.net/modules/sgd.html and
# http://scikit-learn.sourceforge.net/modules/linear_model.html for
# complete documentation.
from ..externals import _lazy_loader

__getattr__, __dir__, __all__ = _lazy_loader.attach(
    __name__,
    submod_attrs={
        "_base": ["LinearRegression"],
        "_bayes": ["ARDRegression", "BayesianRidge"],
        "_coordinate_descent": [
            "ElasticNet",
            "ElasticNetCV",
            "Lasso",
            "LassoCV",
            "MultiTaskElasticNet",
            "MultiTaskElasticNetCV",
            "MultiTaskLasso",
            "MultiTaskLassoCV",
            "enet_path",
            "lasso_path",
        ],
        "_glm.glm": ["PoissonRegressor", "GammaRegressor", "TweedieRegressor"],
        "_huber": ["HuberRegressor"],
        "_least_angle": [
            "Lars",
            "LarsCV",
            "LassoLars",
            "LassoLarsCV",
            "LassoLarsIC",
            "lars_path",
            "lars_path_gram",
        ],
        "_logistic": ["LogisticRegression", "LogisticRegressionCV"],
        "_omp": [
            "OrthogonalMatchingPursuit",
            "OrthogonalMatchingPursuitCV",
            "orthogonal_mp",
            "orthogonal_mp_gram",
        ],
        "_passive_aggressive": [
            "PassiveAggressiveClassifier",
            "PassiveAggressiveRegressor",
        ],
        "_perceptron": ["Perceptron"],
        "_quantile": ["QuantileRegressor"],
        "_ransac": ["RANSACRegressor"],
        "_ridge": [
            "Ridge",
            "RidgeCV",
            "RidgeClassifier",
            "RidgeClassifierCV",
            "ridge_regression",
        ],
        "_sgd_fast": ["Hinge", "Huber", "Log", "ModifiedHuber", "SquaredLoss"],
        "_stochastic_gradient": ["SGDClassifier", "SGDRegressor", "SGDOneClassSVM"],
        "_theil_sen": ["TheilSenRegressor"],
    },
)
