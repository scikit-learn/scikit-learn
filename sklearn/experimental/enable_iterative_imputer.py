"""Enables IterativeImputer

The API and results of this estimators might change without any deprecation
cycle.

Importing this file dynamically sets :class:`sklearn.impute.IterativeImputer`
as an attribute of the impute module::

    >>> # explicitly require this experimental feature
    >>> from sklearn.experimental import enable_iterative_imputer  # noqa
    >>> # now you can import normally from ensemble
    >>> from sklearn.impute import IterativeImputer
"""

from ..impute._iterative import IterativeImputer
from .. import impute

impute.IterativeImputer = IterativeImputer
impute.__all__ += ['IterativeImputer']
