"""Feature selection algorithms.

These include univariate filter selection methods and the recursive feature elimination
algorithm.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from sklearn.feature_selection._base import SelectorMixin
from sklearn.feature_selection._from_model import SelectFromModel
from sklearn.feature_selection._mutual_info import (
    mutual_info_classif,
    mutual_info_regression,
)
from sklearn.feature_selection._rfe import RFE, RFECV
from sklearn.feature_selection._sequential import SequentialFeatureSelector
from sklearn.feature_selection._univariate_selection import (
    GenericUnivariateSelect,
    SelectFdr,
    SelectFpr,
    SelectFwe,
    SelectKBest,
    SelectPercentile,
    chi2,
    f_classif,
    f_oneway,
    f_regression,
    r_regression,
)
from sklearn.feature_selection._variance_threshold import VarianceThreshold

__all__ = [
    "RFE",
    "RFECV",
    "GenericUnivariateSelect",
    "SelectFdr",
    "SelectFpr",
    "SelectFromModel",
    "SelectFwe",
    "SelectKBest",
    "SelectPercentile",
    "SelectorMixin",
    "SequentialFeatureSelector",
    "VarianceThreshold",
    "chi2",
    "f_classif",
    "f_oneway",
    "f_regression",
    "mutual_info_classif",
    "mutual_info_regression",
    "r_regression",
]
