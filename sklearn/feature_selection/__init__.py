"""
The :mod:`sklearn.feature_selection` module implements feature selection
algorithms. It currently includes univariate filter selection methods and the
recursive feature elimination algorithm.
"""
from ..externals import _lazy_loader

__getattr__, __dir__, __all__ = _lazy_loader.attach(
    __name__,
    submod_attrs={
        "_base": ["SelectorMixin"],
        "_from_model": ["SelectFromModel"],
        "_mutual_info": ["mutual_info_classif", "mutual_info_regression"],
        "_rfe": ["RFE", "RFECV"],
        "_sequential": ["SequentialFeatureSelector"],
        "_univariate_selection": [
            "f_oneway",
            "SelectFwe",
            "SelectFpr",
            "f_regression",
            "SelectFdr",
            "chi2",
            "SelectKBest",
            "GenericUnivariateSelect",
            "SelectPercentile",
            "f_classif",
            "r_regression",
        ],
        "_variance_threshold": ["VarianceThreshold"],
    },
)
