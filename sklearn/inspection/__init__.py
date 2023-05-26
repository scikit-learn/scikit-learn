"""The :mod:`sklearn.inspection` module includes tools for model inspection."""
from ..externals import _lazy_loader

__getattr__, __dir__, __all__ = _lazy_loader.attach(
    __name__,
    submod_attrs={
        "_partial_dependence": ["partial_dependence"],
        "_permutation_importance": ["permutation_importance"],
        "_plot.decision_boundary": ["DecisionBoundaryDisplay"],
        "_plot.partial_dependence": ["PartialDependenceDisplay"],
    },
)
