"""The :mod:`sklearn.inspection` module includes tools for model inspection."""

# TODO: remove me in 0.24:
# pre-cache the import of the deprecated module to avoid showing
# the partial_dependence function.
# https://github.com/scikit-learn/scikit-learn/issues/15842
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=FutureWarning)
    from . import partial_dependence as _ignored
    del _ignored

from ._partial_dependence import partial_dependence  # noqa
from ._partial_dependence import plot_partial_dependence  # noqa
from ._partial_dependence import PartialDependenceDisplay  # noqa
from ._permutation_importance import permutation_importance  # noqa


__all__ = [
    'partial_dependence',
    'plot_partial_dependence',
    'permutation_importance',
    'PartialDependenceDisplay'
]
