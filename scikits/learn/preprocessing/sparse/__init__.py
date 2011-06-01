
# Authors: Mathieu Blondel <mathieu@mblondel.org>
# License: BSD

import numpy as np
import scipy.sparse as sp

from .. import SampleNormalizer as DenseNormalizer
from .. import Binarizer as DenseBinarizer


class Binarizer(DenseBinarizer):
    """Binarize data according to a threshold"""

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
