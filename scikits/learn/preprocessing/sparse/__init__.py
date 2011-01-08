
# Authors: Mathieu Blondel <mathieu@mblondel.org>
# License: BSD

import numpy as np
import scipy.sparse as sp

from .. import Normalizer as DenseNormalizer
from .. import LengthNormalizer as DenseLengthNormalizer
from .. import Binarizer as DenseBinarizer

from ._preprocessing import normalize_axis1_sparse, \
                            normalize_length_axis1_sparse

class Normalizer(DenseNormalizer):

    def transform(self, X, y=None, copy=True):
        if not sp.isspmatrix_csr(X):
            X = sp.csr_matrix(X)
        elif copy:
            X = X.copy()

        normalize_axis1_sparse(X)

        return X

class LengthNormalizer(DenseNormalizer):

    def transform(self, X, y=None, copy=True):
        if not sp.isspmatrix_csr(X):
            X = sp.csr_matrix(X)
        elif copy:
            X = X.copy()

        normalize_length_axis1_sparse(X)

        return X

class Binarizer(DenseBinarizer):
    """
    Binarize data according to a threshold.
    """

    def __init__(self, threshold=0.0):
        if threshold < 0:
            # FIXME: sparsity structure changed
            raise NotImplementedError
        self.threshold = threshold

    def transform(self, X, y=None, copy=True):
        if not sp.isspmatrix_csr(X) and not sp.isspmatrix_csc(X):
            X = sp.csr_matrix(X)
        elif copy:
            X = X.copy()

        cond = X.data > self.threshold
        not_cond = np.logical_not(cond)

        X.data[cond] = 1
        # FIXME: if enough values became 0, it may be worth changing
        #        the sparsity structure
        X.data[not_cond] = 0

        return X

