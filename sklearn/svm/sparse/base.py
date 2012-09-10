import numpy as np
import scipy.sparse
from abc import ABCMeta, abstractmethod

from ..base import BaseLibSVM


class SparseBaseLibSVM(BaseLibSVM):
    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self, impl, kernel, degree, gamma, coef0,
                 tol, C, nu, epsilon, shrinking, probability, cache_size,
                 class_weight, verbose):

        super(SparseBaseLibSVM, self).__init__(impl, kernel, degree, gamma,
                coef0, tol, C, nu, epsilon, shrinking, probability, cache_size,
                True, class_weight, verbose)

    def fit(self, X, y, sample_weight=None):
        X = scipy.sparse.csr_matrix(X, dtype=np.float64)
        return super(SparseBaseLibSVM, self).fit(X, y,
                                                 sample_weight=sample_weight)
