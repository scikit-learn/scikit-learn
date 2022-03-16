"""
The :mod:`sklearn.tree` module includes decision tree-based models for
classification and regression.
"""

from ._classes import BaseDecisionTree
from ._classes import DecisionTreeClassifier
from ._classes import DecisionTreeRegressor
from ._classes import ExtraTreeClassifier
from ._classes import ExtraTreeRegressor
from ._classes import ObliqueDecisionTreeClassifier
from ._export import export_graphviz, plot_tree, export_text

__all__ = [
    "BaseDecisionTree",
    "DecisionTreeClassifier",
    "DecisionTreeRegressor",
    "ExtraTreeClassifier",
    "ExtraTreeRegressor",
    "ObliqueDecisionTreeClassifier",
    "export_graphviz",
    "plot_tree",
    "export_text",
]
