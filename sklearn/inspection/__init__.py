"""The :mod:`sklearn.inspection` module includes tools for model inspection."""
from .partial_dependence import partial_dependence
from .partial_dependence import plot_partial_dependence


__all__ = [
    'partial_dependence',
    'plot_partial_dependence',
]
