"""
The :mod:`sklearn.svm.sparse` module includes Support Vector Machine algorithms
for sparse matrices.

This module should have the same API as :mod:`sklearn.svm`, except
that matrices are expected to be in some sparse format supported by
scipy.sparse.

.. note::

    Some fields, like `dual_coef_` are not sparse matrices strictly speaking.
    However, they are converted to a sparse matrix for consistency and
    efficiency when multiplying to other sparse matrices.
"""

# see http://scikit-learn.sourceforge.net/modules/svm.html
# Author: Fabian Pedregosa <fabian.pedregosa@inria.fr> with help from
#         the scikit-learn community.
# License: New BSD, (C) INRIA 2010

from .classes import SVC, NuSVC, SVR, NuSVR, OneClassSVM, LinearSVC
from .. import libsvm_sparse as libsvm
