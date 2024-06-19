"""Support vector machine algorithms."""

# See http://scikit-learn.sourceforge.net/modules/svm.html for complete
# documentation.

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ._bounds import l1_min_c
from ._classes import SVC, SVR, LinearSVC, LinearSVR, NuSVC, NuSVR, OneClassSVM

__all__ = [
    "LinearSVC",
    "LinearSVR",
    "NuSVC",
    "NuSVR",
    "OneClassSVM",
    "SVC",
    "SVR",
    "l1_min_c",
]
