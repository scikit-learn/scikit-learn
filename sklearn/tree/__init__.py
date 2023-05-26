"""
The :mod:`sklearn.tree` module includes decision tree-based models for
classification and regression.
"""
from ..externals import _lazy_loader

__getattr__, __dir__, __all__ = _lazy_loader.attach(
    __name__,
    submod_attrs={
        "_classes": [
            "BaseDecisionTree",
            "DecisionTreeClassifier",
            "DecisionTreeRegressor",
            "ExtraTreeClassifier",
            "ExtraTreeRegressor",
        ],
        "_export": ["export_graphviz", "plot_tree", "export_text"],
    },
)
