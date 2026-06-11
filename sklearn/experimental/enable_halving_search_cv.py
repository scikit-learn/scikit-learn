"""This is now a no-op and can be safely removed from your code.

It used to enable the use of
:class:`~sklearn.model_selection.HalvingGridSearchCV` and
:class:`~sklearn.model_selection.HalvingRandomSearchCV` when they were still
:term:`experimental`, but these estimators are now stable and can be imported
normally from `sklearn.model_selection`.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# TODO(1.14): remove this file

import warnings

warnings.warn(
    "Since version 1.10, it is not needed to import enable_halving_search_cv "
    "anymore. HalvingGridSearchCV and HalvingRandomSearchCV are now stable and "
    "can be imported normally from sklearn.model_selection.",
    UserWarning,
)
