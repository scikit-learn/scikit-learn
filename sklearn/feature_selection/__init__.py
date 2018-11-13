"""
The :mod:`sklearn.feature_selection` module implements feature selection
algorithms. It currently includes univariate filter selection methods and the
recursive feature elimination algorithm.
"""

from .univariate_selection import (
    chi2, f_classif, f_oneway, f_regression, SelectPercentile,
    SelectKBest, SelectFpr, SelectFdr, SelectFwe, GenericUnivariateSelect,
    ANOVAFScorerClassification, ANOVAFScorerRegression, Chi2Scorer,
    MutualInfoScorerClassification, MutualInfoScorerRegression)
from .variance_threshold import VarianceThreshold
from .rfe import RFE, RFECV
from .from_model import SelectFromModel
from .mutual_info_ import mutual_info_regression, mutual_info_classif


__all__ = ['GenericUnivariateSelect',
           'RFE',
           'RFECV',
           'SelectFdr',
           'SelectFpr',
           'SelectFwe',
           'SelectKBest',
           'SelectFromModel',
           'SelectPercentile',
           'VarianceThreshold',
           'chi2',
           'f_classif',
           'f_oneway',
           'f_regression',
           'mutual_info_classif',
           'mutual_info_regression',
           'ANOVAFScorerClassification',
           'ANOVAFScorerRegression',
           'Chi2Scorer',
           'MutualInfoScorerClassification',
           'MutualInfoScorerRegression']
