"""Support vector machine algorithms."""

# See http://scikit-learn.sourceforge.net/modules/svm.html for complete
# documentation.

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from sklearn.svm._bounds import l1_min_c
from sklearn.svm._classes import (
    SVC,
    SVR,
    LinearSVC,
    LinearSVR,
    NuSVC,
    NuSVR,
    OneClassSVM,
)

__all__ = [
    "SVC",
    "SVR",
    "LinearSVC",
    "LinearSVR",
    "NuSVC",
    "NuSVR",
    "OneClassSVM",
    "l1_min_c",
]
