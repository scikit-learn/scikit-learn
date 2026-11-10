"""This is now a no-op and can be safely removed from your code.

It used to enable the use of
:class:`~sklearn.ensemble.HistGradientBoostingClassifier` and
:class:`~sklearn.ensemble.HistGradientBoostingRegressor` when they were still
:term:`experimental`, but these estimators are now stable and can be imported
normally from `sklearn.ensemble`.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# TODO(1.12): remove this file

import warnings

warnings.warn(
    "The sklearn.experimental.enable_hist_gradient_boosting module is "
    "deprecated since version 1.10 and will be removed in 1.12. It is no "
    "longer needed: HistGradientBoostingClassifier and "
    "HistGradientBoostingRegressor are stable and can be imported normally "
    "from sklearn.ensemble.",
    FutureWarning,
)
