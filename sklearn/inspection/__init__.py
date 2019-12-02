"""The :mod:`sklearn.inspection` module includes tools for model inspection."""
from ._partial_dependence import partial_dependence
from ._partial_dependence import plot_partial_dependence
from ._partial_dependence import PartialDependenceDisplay
from ._permutation_importance import permutation_importance
from ._plot.decision_boundary import plot_decision_boundary

__all__ = [
    'partial_dependence',
    'plot_partial_dependence',
    'permutation_importance',
    'PartialDependenceDisplay',
    'plot_decision_boundary'
]
