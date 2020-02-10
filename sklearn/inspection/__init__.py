"""The :mod:`sklearn.inspection` module includes tools for model inspection."""

# TODO: remove me in 0.24 (as well as the noqa markers) and
# import the partial_dependence func directly from the
# ._partial_dependence module instead.
# Pre-cache the import of the deprecated module so that import
# sklearn.inspection.partial_dependence returns the function as in
# 0.21, instead of the module
# https://github.com/scikit-learn/scikit-learn/issues/15842
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=FutureWarning)
    from .partial_dependence import partial_dependence

from ._partial_dependence import plot_partial_dependence  # noqa
from ._partial_dependence import PartialDependenceDisplay  # noqa
from ._permutation_importance import permutation_importance  # noqa


__all__ = [
    'partial_dependence',
    'plot_partial_dependence',
    'permutation_importance',
    'PartialDependenceDisplay'
]
