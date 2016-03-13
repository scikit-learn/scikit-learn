"""
The :mod:`sklearn.feature_selection` module implements feature selection
algorithms. It currently includes univariate filter selection methods and the
recursive feature elimination algorithm.
"""

from .univariate_selection import chi2
from .univariate_selection import ig
from .univariate_selection import igr
from .univariate_selection import f_classif
from .univariate_selection import f_oneway
from .univariate_selection import f_regression
from .univariate_selection import SelectPercentile
from .univariate_selection import SelectKBest
from .univariate_selection import SelectFpr
from .univariate_selection import SelectFdr
from .univariate_selection import SelectFwe
from .univariate_selection import GenericUnivariateSelect
from .mutual_info_ import mutual_info_classif
from .mutual_info_ import mutual_info_regression

from .variance_threshold import VarianceThreshold

from .rfe import RFE
from .rfe import RFECV

from .from_model import SelectFromModel

__all__ = ['GenericUnivariateSelect',
           'RFE',
           'RFECV',
           'SelectFdr',
           'SelectFpr',
           'SelectFwe',
           'SelectKBest',
           'SelectPercentile',
           'VarianceThreshold',
           'chi2',
           'ig',
           'igr',
           'f_classif',
           'f_oneway',
           'f_regression',
           'mutual_info_classif',
           'mutual_info_regression',
           'SelectFromModel']
