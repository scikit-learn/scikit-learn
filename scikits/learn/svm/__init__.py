"""
Module that implements Support Vector Machnine related algorithms.

See http://scikit-learn.sourceforge.net/modules/svm.html for complete
documentation.
"""

from .libsvm import SVC, NuSVC, SVR, NuSVR, OneClassSVM
from .liblinear import LinearSVC
from . import sparse
