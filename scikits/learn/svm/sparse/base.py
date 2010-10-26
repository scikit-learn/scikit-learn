"""
Support Vector Machine algorithms for sparse matrices.

Warning: this module is a work in progress. It is not tested and surely
contains bugs.

Notes
-----

Some fields, like dual_coef_ are not sparse matrices strictly speaking.
However, they are converted to a sparse matrix for consistency and
efficiency when multiplying to other sparse matrices.

Author: Fabian Pedregosa <fabian.pedregosa@inria.fr>
License: New BSD
"""

import numpy as np
from scipy import sparse

from ...base import ClassifierMixin
from ..base import BaseLibSVM, BaseLibLinear

from ._libsvm_sparse import libsvm_sparse_train, \
     libsvm_sparse_predict, set_verbosity_wrap

from .. import _liblinear

class SparseBaseLibSVM(BaseLibSVM):

    _kernel_types = ['linear', 'poly', 'rbf', 'sigmoid', 'precomputed']
    _svm_types = ['c_svc', 'nu_svc', 'one_class', 'epsilon_svr', 'nu_svr']

    def __init__(self, impl, kernel, degree, gamma, coef0, cache_size,
                 eps, C, nu, p, shrinking, probability):
        assert impl in self._svm_types, \
            "impl should be one of %s, %s was given" % (
                self._svm_types, impl)
        assert kernel in self._kernel_types or callable(kernel), \
            "kernel should be one of %s or a callable, %s was given." % (
                self._kernel_types, kernel)
        self.kernel = kernel
        self.impl = impl
        self.degree = degree
        self.gamma = gamma
        self.coef0 = coef0
        self.cache_size = cache_size
        self.eps = eps
        self.C = C
        self.nu = nu
        self.p = p
        self.shrinking = shrinking
        self.probability = probability

        # container for when we call fit
        self._support_data    = np.empty (0, dtype=np.float64, order='C')
        self._support_indices = np.empty (0, dtype=np.int32, order='C')
        self._support_indptr  = np.empty (0, dtype=np.int32, order='C')

        # strictly speaking, dual_coef is not sparse (see Notes above)
        self._dual_coef_data    = np.empty (0, dtype=np.float64, order='C')
        self._dual_coef_indices = np.empty (0, dtype=np.int32,   order='C')
        self._dual_coef_indptr  = np.empty (0, dtype=np.int32,   order='C')
        self.intercept_         = np.empty (0, dtype=np.float64, order='C')

        # only used in classification
        self.n_support = np.empty(0, dtype=np.int32, order='C')


    def fit(self, X, Y, class_weight={}):
        """
        X is expected to be a sparse matrix. For maximum effiency, use a
        sparse matrix in csr format (scipy.sparse.csr_matrix)
        """

        X = sparse.csr_matrix(X)
        X.data = np.asanyarray(X.data, dtype=np.float64, order='C')
        Y      = np.asanyarray(Y,      dtype=np.float64, order='C')

        solver_type = self._svm_types.index(self.impl)
        kernel_type = self._kernel_types.index(self.kernel)

        self.weight       = np.asarray(class_weight.values(),
                                      dtype=np.float64, order='C')
        self.weight_label = np.asarray(class_weight.keys(),
                                       dtype=np.int32, order='C')

        if (kernel_type == 2) and (self.gamma == 0):
            # if custom gamma is not provided ...
            self.gamma = 1.0/X.shape[0]

        self.label_, self.probA_, self.probB_ = libsvm_sparse_train (
                 X.shape[1], X.data, X.indices, X.indptr, Y,
                 solver_type, kernel_type, self.degree,
                 self.gamma, self.coef0, self.eps, self.C,
                 self._support_data, self._support_indices,
                 self._support_indptr, self._dual_coef_data,
                 self.intercept_, self.weight_label, self.weight,
                 self.n_support, self.nu, self.cache_size, self.p,
                 int(self.shrinking),
                 int(self.probability))

        n_class = len(self.label_) - 1
        n_SV = self._support_indptr.size - 1

        dual_coef_indices = np.tile(np.arange(n_SV), n_class)
        dual_coef_indptr = np.arange(0, dual_coef_indices.size + 1,
                                     dual_coef_indices.size / n_class)

        # this will fail if n_SV is zero. This is a limitation
        # in scipy.sparse, which does not permit empty matrices
        self.support_vectors_ = sparse.csr_matrix((self._support_data,
                                           self._support_indices,
                                           self._support_indptr),
                                           (n_SV, X.shape[1]) )

        self.dual_coef_ = sparse.csr_matrix((self._dual_coef_data,
                                             dual_coef_indices,
                                             dual_coef_indptr),
                                            (n_class, n_SV)
                                            )
        return self


    def predict(self, T):
        """
        This function does classification or regression on an array of
        test vectors T.

        For a classification model, the predicted class for each
        sample in T is returned.  For a regression model, the function
        value of T calculated is returned.

        For an one-class model, +1 or -1 is returned.

        Parameters
        ----------
        T : scipy.sparse.csr, shape = [n_samples, n_features]

        Returns
        -------
        C : array, shape = [nsample]
        """
        T = sparse.csr_matrix(T)
        T.data = np.asanyarray(T.data, dtype=np.float64, order='C')
        kernel_type = self._kernel_types.index(self.kernel)

        return libsvm_sparse_predict (T.data, T.indices, T.indptr,
                      self.support_vectors_.data,
                      self.support_vectors_.indices,
                      self.support_vectors_.indptr,
                      self.dual_coef_.data, self.intercept_,
                      self._svm_types.index(self.impl), kernel_type,
                      self.degree, self.gamma, self.coef0, self.eps,
                      self.C, self.weight_label, self.weight, self.nu,
                      self.cache_size, self.p, self.shrinking,
                      self.probability, self.n_support, self.label_,
                      self.probA_, self.probB_)


class SparseBaseLibLinear(BaseLibLinear):

    def fit(self, X, Y, **params):
        """
        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.
        Y : array, shape = [n_samples]
            Target vector relative to X
        """
        self._set_params(**params)
        X = sparse.csr_matrix(X)
        X.data = np.asanyarray(X.data, dtype=np.float64, order='C')
        Y = np.asanyarray(Y, dtype=np.int32, order='C')

        self.raw_coef_, self.label_ = \
                       _liblinear.csr_train_wrap(X.shape[1], X.data, X.indices,
                       X.indptr, Y,
                       self._get_solver_type(),
                       self.eps, self._get_bias(), self.C, self._weight_label,
                       self._weight)
        return self

    def predict(self, T):
        T = sparse.csr_matrix(T)
        T.data = np.asanyarray(T.data, dtype=np.float64, order='C')
        return _liblinear.csr_predict_wrap(T.shape[1],
                                      T.data, T.indices, T.indptr,
                                      self.raw_coef_,
                                      self._get_solver_type(),
                                      self.eps, self.C,
                                      self._weight_label,
                                      self._weight, self.label_,
                                      self._get_bias())


set_verbosity_wrap(0)
