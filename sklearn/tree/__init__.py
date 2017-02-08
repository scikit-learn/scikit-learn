"""
The :mod:`sklearn.tree` module includes decision tree-based models for
classification and regression.
"""

from .tree import DecisionTreeClassifier
from .export import export_graphviz

__all__ = ["DecisionTreeClassifier", "export_graphviz"]
