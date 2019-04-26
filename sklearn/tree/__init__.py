"""
The :mod:`sklearn.tree` module includes decision tree-based models for
classification and regression.
"""

from .tree import DecisionTreeClassifier
from .tree import DecisionTreeRegressor
from .tree import ExtraTreeClassifier
from .tree import ExtraTreeRegressor
from .tree_m5 import M5Base
from .tree_m5 import M5Prime
from .export import export_graphviz, plot_tree, export_text

__all__ = ["DecisionTreeClassifier", "DecisionTreeRegressor",
           "ExtraTreeClassifier", "ExtraTreeRegressor",
           "M5Base", "M5Prime",
           "export_graphviz", "plot_tree",  "export_text"]
