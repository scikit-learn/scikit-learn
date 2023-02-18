"""The :mod:`sklearn.inspection` module includes tools for model inspection."""


from ._permutation_importance import permutation_importance
from ._plot.decision_boundary import DecisionBoundaryDisplay

from ._partial_dependence import partial_dependence
from ._plot.partial_dependence import PartialDependenceDisplay

from ._metric_threshold_curve import metric_threshold_curve


__all__ = [
    "partial_dependence",
    "permutation_importance",
    "PartialDependenceDisplay",
    "DecisionBoundaryDisplay",
    "metric_threshold_curve",
]
