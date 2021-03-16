"""
The :mod:`sklearn.tree` module includes decision tree-based models for
classification and regression.
"""

from ._classes import BaseDecisionTree
from ._classes import DecisionTreeClassifier
from ._classes import DecisionTreeRegressor
from ._classes import ExtraTreeClassifier
from ._classes import ExtraTreeRegressor
from ._classes_m5 import M5Base
from ._classes_m5 import M5Prime
from ._export import export_graphviz, plot_tree, export_text
from ._export_m5p import export_text_m5  # temporary

__all__ = ["BaseDecisionTree",
           "DecisionTreeClassifier", "DecisionTreeRegressor",
           "ExtraTreeClassifier", "ExtraTreeRegressor", "M5Base", "M5Prime",
           "export_text_m5",
           "export_graphviz", "plot_tree", "export_text"]
