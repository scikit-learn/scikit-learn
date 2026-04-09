"""Enables IterativeImputer

The API and results of this estimator might change without any deprecation
cycle.

Importing this file dynamically sets :class:`~sklearn.impute.IterativeImputer`
as an attribute of the impute module::

    >>> # explicitly require this experimental feature
    >>> from sklearn.experimental import enable_iterative_imputer  # noqa
    >>> # now you can import normally from impute
    >>> from sklearn.impute import IterativeImputer
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from sklearn import impute
from sklearn.impute._iterative import IterativeImputer

# use settattr to avoid mypy errors when monkeypatching
setattr(impute, "IterativeImputer", IterativeImputer)
impute.__all__ += ["IterativeImputer"]
