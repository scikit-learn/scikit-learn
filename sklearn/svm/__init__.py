"""
The :mod:`sklearn.svm` module includes Support Vector Machine algorithms.
"""

# See http://scikit-learn.sourceforge.net/modules/svm.html for complete
# documentation.

# Author: Fabian Pedregosa <fabian.pedregosa@inria.fr> with help from
#         the scikit-learn community. LibSVM and LibLinear are copyright
#         of their respective owners.
# License: New BSD, (C) INRIA 2010

from .classes import SVC, NuSVC, SVR, NuSVR, OneClassSVM, LinearSVC
from .bounds import l1_min_c
from . import sparse, libsvm, liblinear
