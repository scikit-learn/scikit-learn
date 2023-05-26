"""
The :mod:`sklearn.svm` module includes Support Vector Machine algorithms.
"""

# See http://scikit-learn.sourceforge.net/modules/svm.html for complete
# documentation.

# Author: Fabian Pedregosa <fabian.pedregosa@inria.fr> with help from
#         the scikit-learn community. LibSVM and LibLinear are copyright
#         of their respective owners.
# License: BSD 3 clause (C) INRIA 2010
from ..externals import _lazy_loader

__getattr__, __dir__, __all__ = _lazy_loader.attach(
    __name__,
    submod_attrs={
        "_classes": [
            "SVC",
            "NuSVC",
            "SVR",
            "NuSVR",
            "OneClassSVM",
            "LinearSVC",
            "LinearSVR",
        ],
        "_bounds": ["l1_min_c"],
    },
)
