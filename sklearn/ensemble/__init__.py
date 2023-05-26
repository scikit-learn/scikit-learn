"""
The :mod:`sklearn.ensemble` module includes ensemble-based methods for
classification, regression and anomaly detection.
"""
from ..externals import _lazy_loader

__getattr__, __dir__, __all__ = _lazy_loader.attach(
    __name__,
    submod_attrs={
        "_bagging": ["BaggingClassifier", "BaggingRegressor"],
        "_base": ["BaseEnsemble"],
        "_forest": [
            "RandomForestClassifier",
            "RandomForestRegressor",
            "RandomTreesEmbedding",
            "ExtraTreesClassifier",
            "ExtraTreesRegressor",
        ],
        "_gb": ["GradientBoostingClassifier", "GradientBoostingRegressor"],
        "_hist_gradient_boosting.gradient_boosting": [
            "HistGradientBoostingClassifier",
            "HistGradientBoostingRegressor",
        ],
        "_iforest": ["IsolationForest"],
        "_stacking": ["StackingClassifier", "StackingRegressor"],
        "_voting": ["VotingClassifier", "VotingRegressor"],
        "_weight_boosting": ["AdaBoostClassifier", "AdaBoostRegressor"],
    },
)
