"""The :mod:`sklearn.inspection` module includes tools for model inspection."""
from .partial_dependence import partial_dependence
from .partial_dependence import plot_partial_dependence
from ._plot_estimators import export_html


__all__ = [
    'partial_dependence',
    'plot_partial_dependence',
    'export_html'
]
