"""Support vector machine algorithms with alpha seeding and other small changes."""

# See http://scikit-learn.sourceforge.net/modules/svm.html for complete
# documentation.

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ._classes import SVC, SVR, NuSVC, NuSVR, OneClassSVM

__all__ = [
    "SVC",
    "SVR",
    "NuSVC",
    "NuSVR",
    "OneClassSVM",
]
