import numpy as np
import scipy.sparse

from ..base import BaseLibSVM, BaseLibLinear, LIBSVM_IMPL, _get_class_weight
from . import libsvm
from .. import liblinear


class SparseBaseLibSVM(BaseLibSVM):

    _kernel_types = ['linear', 'poly', 'rbf', 'sigmoid', 'precomputed']

    def __init__(self, impl, kernel, degree, gamma, coef0,
                 tol, C, nu, epsilon, shrinking, probability, cache_size,
                 scale_C):

        assert kernel in self._kernel_types, \
               "kernel should be one of %s, "\
               "%s was given." % (self._kernel_types, kernel)

        super(SparseBaseLibSVM, self).__init__(impl, kernel, degree, gamma,
                coef0, tol, C, nu, epsilon, shrinking, probability, cache_size,
                scale_C)

    def fit(self, X, y, class_weight=None, sample_weight=None):
        """
        Fit the SVM model according to the given training data and parameters.

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

        X = scipy.sparse.csr_matrix(X)
        X.data = np.asarray(X.data, dtype=np.float64, order='C')
        y = np.asarray(y, dtype=np.float64, order='C')
        sample_weight = np.asarray([] if sample_weight is None
                                      else sample_weight, dtype=np.float64)

        if X.shape[0] != y.shape[0]:
            raise ValueError("X and y have incompatible shapes: %r vs %r\n"
                             "Note: Sparse matrices cannot be indexed w/"
                             "boolean masks (use `indices=True` in CV)."
                             % (X.shape, y.shape))

        if sample_weight.shape[0] > 0 and sample_weight.shape[0] != X.shape[0]:
            raise ValueError("sample_weight and X have incompatible shapes:"
                             "%r vs %r\n"
                             "Note: Sparse matrices cannot be indexed w/"
                             "boolean masks (use `indices=True` in CV)."
                             % (sample_weight.shape, X.shape))

        solver_type = LIBSVM_IMPL.index(self.impl)
        kernel_type = self._kernel_types.index(self.kernel)

        self.class_weight, self.class_weight_label = \
                     _get_class_weight(class_weight, y)

        if (kernel_type in [1, 2]) and (self.gamma == 0):
            # if custom gamma is not provided ...
            self.gamma = 1.0 / X.shape[1]

        C = self.C
        if self.scale_C:
            C = C / float(X.shape[0])

        self.support_vectors_, dual_coef_data, self.intercept_, self.label_, \
            self.n_support_, self.probA_, self.probB_ = \
            libsvm.libsvm_sparse_train(
                 X.shape[1], X.data, X.indices, X.indptr, y, solver_type,\
                 kernel_type, self.degree, self.gamma, self.coef0, self.tol,\
                 C, self.class_weight_label, self.class_weight,\
                 sample_weight, self.nu, self.cache_size, self.epsilon,\
                 int(self.shrinking), int(self.probability))

        n_class = len(self.label_) - 1
        n_SV = self.support_vectors_.shape[0]

        dual_coef_indices = np.tile(np.arange(n_SV), n_class)
        dual_coef_indptr = np.arange(0, dual_coef_indices.size + 1,
                                     dual_coef_indices.size / n_class)
        self.dual_coef_ = scipy.sparse.csr_matrix(
            (dual_coef_data, dual_coef_indices, dual_coef_indptr),
            (n_class, n_SV))
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
        T = scipy.sparse.csr_matrix(T)
        T.data = np.asarray(T.data, dtype=np.float64, order='C')
        kernel_type = self._kernel_types.index(self.kernel)

        return libsvm.libsvm_sparse_predict(T.data, T.indices, T.indptr,
                      self.support_vectors_.data,
                      self.support_vectors_.indices,
                      self.support_vectors_.indptr,
                      self.dual_coef_.data, self.intercept_,
                      LIBSVM_IMPL.index(self.impl), kernel_type,
                      self.degree, self.gamma, self.coef0, self.tol,
                      self.C, self.class_weight_label, self.class_weight,
                      self.nu, self.epsilon, self.shrinking,
                      self.probability, self.n_support_, self.label_,
                      self.probA_, self.probB_)

    def predict_proba(self, X):
        """
        This function does classification or regression on a test vector X
        given a model with probability information.

        Parameters
        ----------
        X : scipy.sparse.csr, shape = [n_samples, n_features]

        Returns
        -------
        X : array-like, shape = [n_samples, n_classes]
            Returns the probability of the sample for each class in
            the model, where classes are ordered by arithmetical
            order.

        Notes
        -----
        The probability model is created using cross validation, so
        the results can be slightly different than those obtained by
        predict. Also, it will meaningless results on very small
        datasets.
        """

        if not self.probability:
            raise ValueError(
                    "probability estimates must be enabled to use this method")

        if self.impl not in ('c_svc', 'nu_svc'):
            raise NotImplementedError("predict_proba only implemented for " +
                                      "SVC and NuSVC")

        X = scipy.sparse.csr_matrix(X)
        X.data = np.asarray(X.data, dtype=np.float64, order='C')
        kernel_type = self._kernel_types.index(self.kernel)

        return libsvm.libsvm_sparse_predict_proba(
            X.data, X.indices, X.indptr,
            self.support_vectors_.data,
            self.support_vectors_.indices,
            self.support_vectors_.indptr,
            self.dual_coef_.data, self.intercept_,
            LIBSVM_IMPL.index(self.impl), kernel_type,
            self.degree, self.gamma, self.coef0, self.tol,
            self.C, self.class_weight_label, self.class_weight,
            self.nu, self.epsilon, self.shrinking,
            self.probability, self.n_support_, self.label_,
            self.probA_, self.probB_)


class SparseBaseLibLinear(BaseLibLinear):
    def fit(self, X, y, class_weight=None):
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
        X = scipy.sparse.csr_matrix(X)
        y = np.asarray(y, dtype=np.int32, order='C')
        if X.shape[0] != y.shape[0]:
            raise ValueError("X and y have incompatible shapes.\n" +
                             "Note: Sparse matrices cannot be indexed w/" +
                             "boolean masks (use `indices=True` in CV).")

        X.data = np.asarray(X.data, dtype=np.float64, order='C')

        self.class_weight, self.class_weight_label = \
                     _get_class_weight(class_weight, y)

        C = self.C
        if self.scale_C:
            C = C / float(X.shape[0])

        self.raw_coef_, self.label_ = \
                       liblinear.csr_train_wrap(X.shape[1], X.data, X.indices,
                       X.indptr, y,
                       self._get_solver_type(),
                       self.tol, self._get_bias(), C,
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
        X = scipy.sparse.csr_matrix(X)
        self._check_n_features(X)
        X.data = np.asarray(X.data, dtype=np.float64, order='C')

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
        X = scipy.sparse.csr_matrix(X)
        self._check_n_features(X)
        X.data = np.asarray(X.data, dtype=np.float64, order='C')

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
