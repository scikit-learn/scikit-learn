import numpy as np

from ..base import BaseLibSVM, BaseLibLinear, _get_class_weight
from . import libsvm
from .. import liblinear


class SparseBaseLibSVM(BaseLibSVM):

    _kernel_types = ['linear', 'poly', 'rbf', 'sigmoid', 'precomputed']
    _svm_types = ['c_svc', 'nu_svc', 'one_class', 'epsilon_svr', 'nu_svr']

    def __init__(self, impl, kernel, degree, gamma, coef0, cache_size,
                 tol, C, nu, p, shrinking, probability):

        assert impl in self._svm_types, \
            "impl should be one of %s, %s was given" % (
                self._svm_types, impl)

        assert kernel in self._kernel_types, \
               "kernel should be one of %s, "\
               "%s was given." % (self._kernel_types, kernel)

        self.kernel = kernel
        self.impl = impl
        self.degree = degree
        self.gamma = gamma
        self.coef0 = coef0
        self.cache_size = cache_size
        self.tol = tol
        self.C = C
        self.nu = nu
        self.p = p
        self.shrinking = shrinking
        self.probability = probability

        # container for when we call fit
        self._support_data = np.empty(0, dtype=np.float64, order='C')
        self._support_indices = np.empty(0, dtype=np.int32, order='C')
        self._support_indptr = np.empty(0, dtype=np.int32, order='C')

        # strictly speaking, dual_coef is not sparse (see Notes above)
        self._dual_coef_data = np.empty(0, dtype=np.float64, order='C')
        self._dual_coef_indices = np.empty(0, dtype=np.int32,   order='C')
        self._dual_coef_indptr = np.empty(0, dtype=np.int32,   order='C')
        self.intercept_ = np.empty(0, dtype=np.float64, order='C')

        # only used in classification
        self.n_support = np.empty(0, dtype=np.int32, order='C')

    def fit(self, X, y, class_weight={}, sample_weight=[], **params):
        """
        Fit the SVM model according to the given training data and
        parameters.

        Parameters
        ----------
        X : sparse matrix, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples and
            n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target values (integers in classification, real numbers in
            regression)

        class_weight : {dict, 'auto'}, optional
            Weights associated with classes in the form
            {class_label : weight}. If not given, all classes are
            supposed to have weight one.

            The 'auto' mode uses the values of y to automatically adjust
            weights inversely proportional to class frequencies.

        sample_weight : array-like, shape = [n_samples], optional
            Weights applied to individual samples (1. for unweighted).

        Returns
        -------
        self : object
            Returns an instance of self.

        Notes
        -----
        For maximum effiency, use a sparse matrix in csr format
        (scipy.sparse.csr_matrix)
        """
        self._set_params(**params)

        import scipy.sparse
        X = scipy.sparse.csr_matrix(X)
        X.data = np.asanyarray(X.data, dtype=np.float64, order='C')
        y = np.asanyarray(y, dtype=np.float64, order='C')
        sample_weight = np.asanyarray(sample_weight, dtype=np.float64,
                                      order='C')

        solver_type = self._svm_types.index(self.impl)
        kernel_type = self._kernel_types.index(self.kernel)

        self.class_weight, self.class_weight_label = \
                     _get_class_weight(class_weight, y)

        if (kernel_type in [1, 2]) and (self.gamma == 0):
            # if custom gamma is not provided ...
            self.gamma = 1.0 / X.shape[0]

        self.label_, self.probA_, self.probB_ = libsvm.libsvm_sparse_train(
                 X.shape[1], X.data, X.indices, X.indptr, y,
                 solver_type, kernel_type, self.degree, self.gamma,
                 self.coef0, self.tol, self.C, self._support_data,
                 self._support_indices, self._support_indptr,
                 self._dual_coef_data, self.intercept_,
                 self.class_weight_label, self.class_weight, sample_weight,
                 self.n_support, self.nu, self.cache_size, self.p,
                 int(self.shrinking), int(self.probability))

        n_class = len(self.label_) - 1
        n_SV = self._support_indptr.size - 1

        dual_coef_indices = np.tile(np.arange(n_SV), n_class)
        dual_coef_indptr = np.arange(0, dual_coef_indices.size + 1,
                                     dual_coef_indices.size / n_class)

        # this will fail if n_SV is zero. This is a limitation
        # in scipy.sparse, which does not permit empty matrices
        self.support_vectors_ = scipy.sparse.csr_matrix((self._support_data,
                                           self._support_indices,
                                           self._support_indptr),
                                           (n_SV, X.shape[1]))

        self.dual_coef_ = scipy.sparse.csr_matrix((self._dual_coef_data,
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
        C : array, shape = [n_samples]
        """
        import scipy.sparse
        T = scipy.sparse.csr_matrix(T)
        T.data = np.asanyarray(T.data, dtype=np.float64, order='C')
        kernel_type = self._kernel_types.index(self.kernel)

        return libsvm.libsvm_sparse_predict(T.data, T.indices, T.indptr,
                      self.support_vectors_.data,
                      self.support_vectors_.indices,
                      self.support_vectors_.indptr,
                      self.dual_coef_.data, self.intercept_,
                      self._svm_types.index(self.impl), kernel_type,
                      self.degree, self.gamma, self.coef0, self.tol,
                      self.C, self.class_weight_label, self.class_weight,
                      self.nu, self.cache_size, self.p, self.shrinking,
                      self.probability, self.n_support, self.label_,
                      self.probA_, self.probB_)


class SparseBaseLibLinear(BaseLibLinear):

    def fit(self, X, y, class_weight={}, **params):
        """
        Fit the model using X, y as training data.

        Parameters
        ----------
        X : sparse matrix, shape = [n_samples, n_features]
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.
        y : array, shape = [n_samples]
            Target vector relative to X

        Returns
        -------
        self : object
            Returns an instance of self.
        """
        self._set_params(**params)

        import scipy.sparse
        X = scipy.sparse.csr_matrix(X)
        X.data = np.asanyarray(X.data, dtype=np.float64, order='C')
        y = np.asanyarray(y, dtype=np.int32, order='C')

        self.class_weight, self.class_weight_label = \
                     _get_class_weight(class_weight, y)

        self.raw_coef_, self.label_ = \
                       liblinear.csr_train_wrap(X.shape[1], X.data, X.indices,
                       X.indptr, y,
                       self._get_solver_type(),
                       self.tol, self._get_bias(), self.C,
                       self.class_weight_label, self.class_weight)

        return self

    def predict(self, X):
        """
        Predict target values of X according to the fitted model.

        Parameters
        ----------
        X : sparse matrix, shape = [n_samples, n_features]

        Returns
        -------
        C : array, shape = [n_samples]
        """
        import scipy.sparse
        X = scipy.sparse.csr_matrix(X)
        self._check_n_features(X)
        X.data = np.asanyarray(X.data, dtype=np.float64, order='C')

        return liblinear.csr_predict_wrap(X.shape[1], X.data,
                                      X.indices, X.indptr,
                                      self.raw_coef_,
                                      self._get_solver_type(),
                                      self.tol, self.C,
                                      self.class_weight_label,
                                      self.class_weight, self.label_,
                                      self._get_bias())

    def decision_function(self, X):
        """
        Return the decision function of X according to the trained
        model.

        Parameters
        ----------
        X : sparse matrix, shape = [n_samples, n_features]

        Returns
        -------
        T : array-like, shape = [n_samples, n_class]
            Returns the decision function of the sample for each class
            in the model.
        """
        import scipy.sparse
        X = scipy.sparse.csr_matrix(X)
        self._check_n_features(X)
        X.data = np.asanyarray(X.data, dtype=np.float64, order='C')

        dec_func = liblinear.csr_decision_function_wrap(
            X.shape[1], X.data, X.indices, X.indptr, self.raw_coef_,
            self._get_solver_type(), self.tol, self.C,
            self.class_weight_label, self.class_weight, self.label_,
            self._get_bias())

        if len(self.label_) <= 2:
            # in the two-class case, the decision sign needs be flipped
            # due to liblinear's design
            return -dec_func
        else:
            return dec_func

libsvm.set_verbosity_wrap(0)
