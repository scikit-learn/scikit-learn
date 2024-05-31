"""Tools for model inspection."""

from ._h_statistic import h_statistic
from ._partial_dependence import partial_dependence
from ._permutation_importance import permutation_importance
from ._plot.decision_boundary import DecisionBoundaryDisplay
from ._plot.partial_dependence import PartialDependenceDisplay

__all__ = [
    "h_statistic",
    "partial_dependence",
    "permutation_importance",
    "PartialDependenceDisplay",
    "DecisionBoundaryDisplay",
]
