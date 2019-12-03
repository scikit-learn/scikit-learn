"""The :mod:`sklearn.inspection` module includes tools for model inspection."""
from ._partial_dependence import partial_dependence
from ._partial_dependence import plot_partial_dependence
from ._partial_dependence import PartialDependenceDisplay
from ._permutation_importance import permutation_importance

__all__ = [
    'partial_dependence',
    'plot_partial_dependence',
    'permutation_importance',
    'PartialDependenceDisplay'
]
