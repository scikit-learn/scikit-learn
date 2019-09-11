"""
The :mod:`sklearn.tree` module includes decision tree-based models for
classification and regression.
"""

from ._tree_base import DecisionTreeClassifier
from ._tree_base import DecisionTreeRegressor
from ._tree_base import ExtraTreeClassifier
from ._tree_base import ExtraTreeRegressor
from ._export import export_graphviz, plot_tree, export_text

__all__ = ["DecisionTreeClassifier", "DecisionTreeRegressor",
           "ExtraTreeClassifier", "ExtraTreeRegressor", "export_graphviz",
           "plot_tree", "export_text"]
