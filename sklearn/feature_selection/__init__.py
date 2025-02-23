"""Feature selection algorithms.

These include univariate filter selection methods and the recursive feature elimination
algorithm.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ._base import SelectorMixin
from ._from_model import SelectFromModel
from ._mutual_info import mutual_info_classif, mutual_info_regression
from ._rfe import RFE, RFECV
from ._sequential import SequentialFeatureSelector
from ._univariate_selection import (
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
from ._variance_threshold import VarianceThreshold

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
