"""
The :mod:`sklearn.metric_learning` module includes models that find
a metric that makes samples from the same class closer than samples
from different classes.
"""

from .nca import NCATransformer

__all__ = ['NCATransformer']
