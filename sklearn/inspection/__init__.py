"""The :mod:`sklearn.inspection` module includes tools for model inspection."""
from .partial_dependence import partial_dependence
from .partial_dependence import plot_partial_dependence
from ._display_estimator import display_estimator


__all__ = [
    'partial_dependence',
    'plot_partial_dependence',
    'display_estimator'
]
