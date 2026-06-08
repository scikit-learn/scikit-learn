"""This is now a no-op and can be safely removed from your code.

It used to enable the use of :class:`~sklearn.impute.IterativeImputer` when it
was still :term:`experimental`, but this estimator is now stable and can be
imported normally from ``sklearn.impute``.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# Don't remove this file, we don't want to break users code just because the
# feature isn't experimental anymore.

import warnings

warnings.warn(
    "Since version 1.10, "
    "it is not needed to import enable_iterative_imputer anymore. "
    "IterativeImputer is now stable and can be normally imported from "
    "sklearn.impute."
)
