"""
Support Vector Machnine algorithms.

See http://scikit-learn.sourceforge.net/modules/svm.html for complete
documentation.

Author: Fabian Pedregosa <fabian.pedregosa@inria.fr> with help from
        the scikit-learn community. LibSVM and LibLinear are copyright
        of their respective owners.
License: New BSD, (C) INRIA 2010
"""

from .classes import SVC, NuSVC, SVR, NuSVR, OneClassSVM, LinearSVC
from . import sparse, libsvm, liblinear
