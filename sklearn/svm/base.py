import numpy as np
import scipy.sparse as sp
import warnings

from . import libsvm, liblinear
from . import libsvm_sparse
from ..base import BaseEstimator
from ..utils import array2d, atleast2d_or_csr
from ..utils.extmath import safe_sparse_dot


LIBSVM_IMPL = ['c_svc', 'nu_svc', 'one_class', 'epsilon_svr', 'nu_svr']


def _get_class_weight(class_weight, y):
    """Estimate class weights for unbalanced datasets."""
    if class_weight == 'auto':
        uy = np.unique(y)
        weight_label = np.asarray(uy, dtype=np.int32, order='C')
        weight = np.array([1.0 / np.sum(y == i) for i in uy],
                          dtype=np.float64, order='C')
        weight *= uy.shape[0] / np.sum(weight)
    else:
        if class_weight is None:
            keys = values = []
        else:
            keys = class_weight.keys()
            values = class_weight.values()
        weight = np.asarray(values, dtype=np.float64, order='C')
        weight_label = np.asarray(keys, dtype=np.int32, order='C')

    return weight, weight_label


def _one_vs_one_coef(dual_coef, n_support, support_vectors):
    """Generate primal coefficients from dual coefficients
    for the one-vs-one multi class LibSVM in the case
    of a linear kernel."""

    # get 1vs1 weights for all n*(n-1) classifiers.
    # this is somewhat messy.
    # shape of dual_coef_ is nSV * (n_classes -1)
    # see docs for details
    n_class = dual_coef.shape[0] + 1

    # XXX we could do preallocation of coef but
    # would have to take care in the sparse case
    coef = []
    sv_locs = np.cumsum(np.hstack([[0], n_support]))
    for class1 in xrange(n_class):
        # SVs for class1:
        sv1 = support_vectors[sv_locs[class1]:sv_locs[class1 + 1], :]
        for class2 in xrange(class1 + 1, n_class):
            # SVs for class1:
            sv2 = support_vectors[sv_locs[class2]:sv_locs[class2 + 1], :]

            # dual coef for class1 SVs:
            alpha1 = dual_coef[class2 - 1, sv_locs[class1]:sv_locs[class1 + 1]]
            # dual coef for class2 SVs:
            alpha2 = dual_coef[class1, sv_locs[class2]:sv_locs[class2 + 1]]
            # build weight for class1 vs class2

            coef.append(safe_sparse_dot(alpha1, sv1)
                    + safe_sparse_dot(alpha2, sv2))
    return coef


class BaseLibSVM(BaseEstimator):
    """Base class for estimators that use libsvm as backing library

    This implements support vector machine classification and regression.
    """

    _sparse_kernels = ["linear", "poly", "rbf", "sigmoid", "precomputed"]

    def __init__(self, impl, kernel, degree, gamma, coef0,
                 tol, C, nu, epsilon, shrinking, probability, cache_size,
                 scale_C, sparse, class_weight):

        if not impl in LIBSVM_IMPL:
            raise ValueError("impl should be one of %s, %s was given" % (
                LIBSVM_IMPL, impl))
        if hasattr(kernel, '__call__'):
            self.kernel_function = kernel
            self.kernel = 'precomputed'
        else:
            self.kernel = kernel
        if not scale_C:
            warnings.warn('SVM: scale_C will disappear and be assumed to be '
                          'True in scikit-learn 0.12', FutureWarning,
                          stacklevel=2)

        self.impl = impl
        self.degree = degree
        self.gamma = gamma
        self.coef0 = coef0
        self.tol = tol
        self.C = C
        self.nu = nu
        self.epsilon = epsilon
        self.shrinking = shrinking
        self.probability = probability
        self.cache_size = cache_size
        self.scale_C = scale_C
        self.sparse = sparse
        self.class_weight = class_weight

    def fit(self, X, y, class_weight=None, sample_weight=None):
        """Fit the SVM model according to the given training data.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target values (integers in classification, real numbers in
            regression)

        sample_weight : array-like, shape = [n_samples], optional
            Weights applied to individual samples (1. for unweighted).

        Returns
        -------
        self : object
            Returns self.

        Notes
        ------
        If X and y are not C-ordered and contiguous arrays of np.float64 and
        X is not a scipy.sparse.csr_matrix, X and/or y may be copied.

        If X is a dense array, then the other methods will not support sparse
        matrices as input.
        """
        self._sparse = sp.isspmatrix(X) if self.sparse == "auto" else self.sparse
        if class_weight != None:
            warnings.warn("'class_weight' is now an initialization parameter."
                    "Using it in the 'fit' method is deprecated.",
                    DeprecationWarning)
            self.class_weight = class_weight
        fit = self._sparse_fit if self._sparse else self._dense_fit
        fit(X, y, sample_weight)
        return self

    def _dense_fit(self, X, y, sample_weight=None):
        X = np.asarray(X, dtype=np.float64, order='C')
        y = np.asarray(y, dtype=np.float64, order='C')
        sample_weight = np.asarray([] if sample_weight is None
                                      else sample_weight, dtype=np.float64)

        if hasattr(self, 'kernel_function'):
            # you must store a reference to X to compute the kernel in predict
            # TODO: add keyword copy to copy on demand
            self.__Xfit = X
            X = self._compute_kernel(X)

        self.class_weight_, self.class_weight_label_ = \
                     _get_class_weight(self.class_weight, y)

        # check dimensions
        solver_type = LIBSVM_IMPL.index(self.impl)
        if solver_type != 2 and X.shape[0] != y.shape[0]:
            raise ValueError("X and y have incompatible shapes.\n" +
                             "X has %s samples, but y has %s." %
                             (X.shape[0], y.shape[0]))

        if self.kernel == "precomputed" and X.shape[0] != X.shape[1]:
            raise ValueError("X.shape[0] should be equal to X.shape[1]")

        if (self.kernel in ['poly', 'rbf']) and (self.gamma == 0):
            # if custom gamma is not provided ...
            self.gamma = 1.0 / X.shape[1]
        self.shape_fit_ = X.shape

        # set default parameters
        C = self.C
        if getattr(self, 'scale_C', False):
            C = self.C / float(X.shape[0])
        epsilon = self.epsilon
        if epsilon is None:
            epsilon = 0.1

        # we don't pass **self.get_params() to allow subclasses to
        # add other parameters to __init__
        self.support_, self.support_vectors_, self.n_support_, \
        self.dual_coef_, self.intercept_, self.label_, self.probA_, \
        self.probB_ = libsvm.fit(X, y,
            svm_type=solver_type, sample_weight=sample_weight,
            class_weight=self.class_weight_,
            class_weight_label=self.class_weight_label_,
            kernel=self.kernel, C=C, nu=self.nu,
            probability=self.probability, degree=self.degree,
            shrinking=self.shrinking, tol=self.tol, cache_size=self.cache_size,
            coef0=self.coef0, gamma=self.gamma, epsilon=epsilon)

        # In binary case, we need to flip the sign of coef, intercept and
        # decision function. Use self._intercept_ internally.
        self._intercept_ = self.intercept_.copy()
        if len(self.label_) == 2 and self.impl != 'one_class':
            self.intercept_ *= -1

    def _sparse_fit(self, X, y, sample_weight=None):
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

        X = sp.csr_matrix(X)
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
        kernel_type = self._sparse_kernels.index(self.kernel)

        self.class_weight_, self.class_weight_label_ = \
                     _get_class_weight(self.class_weight, y)

        if (kernel_type in [1, 2]) and (self.gamma == 0):
            # if custom gamma is not provided ...
            self.gamma = 1.0 / X.shape[1]

        C = self.C
        if self.scale_C:
            C /= float(X.shape[0])

        self.support_vectors_, dual_coef_data, self.intercept_, self.label_, \
            self.n_support_, self.probA_, self.probB_ = \
            libsvm_sparse.libsvm_sparse_train(
                 X.shape[1], X.data, X.indices, X.indptr, y, solver_type,
                 kernel_type, self.degree, self.gamma, self.coef0, self.tol,
                 C, self.class_weight_label_, self.class_weight_,
                 sample_weight, self.nu, self.cache_size, self.epsilon,
                 int(self.shrinking), int(self.probability))

        # In binary case, we need to flip the sign of coef, intercept and
        # decision function. Use self._intercept_ internally.
        self._intercept_ = self.intercept_.copy()
        if len(self.label_) == 2 and self.impl != 'one_class':
            self.intercept_ *= -1

        n_class = len(self.label_) - 1
        n_SV = self.support_vectors_.shape[0]

        dual_coef_indices = np.tile(np.arange(n_SV), n_class)
        dual_coef_indptr = np.arange(0, dual_coef_indices.size + 1,
                                     dual_coef_indices.size / n_class)
        self.dual_coef_ = sp.csr_matrix(
            (dual_coef_data, dual_coef_indices, dual_coef_indptr),
            (n_class, n_SV))

    def predict(self, X):
        """Perform classification or regression samples in X.

        For a classification model, the predicted class for each
        sample in X is returned.  For a regression model, the function
        value of X calculated is returned.

        For an one-class model, +1 or -1 is returned.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]

        Returns
        -------
        C : array, shape = [n_samples]
        """
        X = self._validate_for_predict(X)
        predict = self._sparse_predict if self._sparse else self._dense_predict
        return predict(X)

    def _dense_predict(self, X):
        X = np.asarray(X, dtype=np.float64, order='C')
        if X.ndim == 1:
            # don't use np.atleast_2d, it doesn't guarantee C-contiguity
            X = np.reshape(X, (1, -1), order='C')
        n_samples, n_features = X.shape
        X = self._compute_kernel(X)

        if self.kernel == "precomputed":
            if X.shape[1] != self.shape_fit_[0]:
                raise ValueError("X.shape[1] = %d should be equal to %d, "
                                 "the number of samples at training time" %
                                 (X.shape[1], self.shape_fit_[0]))
        elif n_features != self.shape_fit_[1]:
            raise ValueError("X.shape[1] = %d should be equal to %d, "
                             "the number of features at training time" %
                             (n_features, self.shape_fit_[1]))

        epsilon = self.epsilon
        if epsilon == None:
            epsilon = 0.1

        svm_type = LIBSVM_IMPL.index(self.impl)
        return libsvm.predict(
            X, self.support_, self.support_vectors_, self.n_support_,
            self.dual_coef_, self._intercept_,
            self.label_, self.probA_, self.probB_,
            svm_type=svm_type,
            kernel=self.kernel, C=self.C, nu=self.nu,
            probability=self.probability, degree=self.degree,
            shrinking=self.shrinking, tol=self.tol, cache_size=self.cache_size,
            coef0=self.coef0, gamma=self.gamma, epsilon=epsilon)

    def _sparse_predict(self, X):
        X = sp.csr_matrix(X, dtype=np.float64)
        kernel_type = self._sparse_kernels.index(self.kernel)

        return libsvm_sparse.libsvm_sparse_predict(
                      X.data, X.indices, X.indptr,
                      self.support_vectors_.data,
                      self.support_vectors_.indices,
                      self.support_vectors_.indptr,
                      self.dual_coef_.data, self._intercept_,
                      LIBSVM_IMPL.index(self.impl), kernel_type,
                      self.degree, self.gamma, self.coef0, self.tol,
                      self.C, self.class_weight_label_, self.class_weight_,
                      self.nu, self.epsilon, self.shrinking,
                      self.probability, self.n_support_, self.label_,
                      self.probA_, self.probB_)

    def predict_proba(self, X):
        """Compute the likehoods each possible outcomes of samples in T.

        The model need to have probability information computed at training
        time: fit with attribute `probability` set to True.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

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
            raise NotImplementedError("predict_proba only implemented for SVC "
                                      "and NuSVC")

        X = self._validate_for_predict(X)
        pred_proba = self._sparse_predict_proba if self._sparse \
                                                else self._dense_predict_proba
        return pred_proba(X)

    def _dense_predict_proba(self, X):
        X = self._compute_kernel(X)

        epsilon = self.epsilon
        if epsilon == None:
            epsilon = 0.1

        svm_type = LIBSVM_IMPL.index(self.impl)
        pprob = libsvm.predict_proba(
            X, self.support_, self.support_vectors_, self.n_support_,
            self.dual_coef_, self._intercept_, self.label_,
            self.probA_, self.probB_,
            svm_type=svm_type, kernel=self.kernel, C=self.C, nu=self.nu,
            probability=self.probability, degree=self.degree,
            shrinking=self.shrinking, tol=self.tol, cache_size=self.cache_size,
            coef0=self.coef0, gamma=self.gamma, epsilon=epsilon)

        return pprob

    def _compute_kernel(self, X):
        """Return the data transformed by a callable kernel"""
        if hasattr(self, 'kernel_function'):
            # in the case of precomputed kernel given as a function, we
            # have to compute explicitly the kernel matrix
            X = np.asarray(self.kernel_function(X, self.__Xfit),
                           dtype=np.float64, order='C')
        return X

    def _sparse_predict_proba(self, X):
        X.data = np.asarray(X.data, dtype=np.float64, order='C')
        kernel_type = self._sparse_kernels.index(self.kernel)

        return libsvm_sparse.libsvm_sparse_predict_proba(
            X.data, X.indices, X.indptr,
            self.support_vectors_.data,
            self.support_vectors_.indices,
            self.support_vectors_.indptr,
            self.dual_coef_.data, self._intercept_,
            LIBSVM_IMPL.index(self.impl), kernel_type,
            self.degree, self.gamma, self.coef0, self.tol,
            self.C, self.class_weight_label_, self.class_weight_,
            self.nu, self.epsilon, self.shrinking,
            self.probability, self.n_support_, self.label_,
            self.probA_, self.probB_)

    def predict_log_proba(self, X):
        """Compute the log likehoods each possible outcomes of samples in X.

        The model need to have probability information computed at training
        time: fit with attribute `probability` set to True.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        X : array-like, shape = [n_samples, n_classes]
            Returns the log-probabilities of the sample for each class in
            the model, where classes are ordered by arithmetical
            order.

        Notes
        -----
        The probability model is created using cross validation, so
        the results can be slightly different than those obtained by
        predict. Also, it will meaningless results on very small
        datasets.
        """
        return np.log(self.predict_proba(X))

    def decision_function(self, X):
        """Distance of the samples X to the separating hyperplane.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        X : array-like, shape = [n_samples, n_class * (n_class-1) / 2]
            Returns the decision function of the sample for each class
            in the model.
        """
        if self._sparse:
            raise ValueError("decision_function not supported for sparse SVM")

        X = array2d(X, dtype=np.float64, order="C")

        epsilon = self.epsilon
        if epsilon == None:
            epsilon = 0.1
        dec_func = libsvm.decision_function(
            X, self.support_, self.support_vectors_, self.n_support_,
            self.dual_coef_, self._intercept_, self.label_,
            self.probA_, self.probB_,
            svm_type=LIBSVM_IMPL.index(self.impl),
            kernel=self.kernel, C=self.C, nu=self.nu,
            probability=self.probability, degree=self.degree,
            shrinking=self.shrinking, tol=self.tol, cache_size=self.cache_size,
            coef0=self.coef0, gamma=self.gamma, epsilon=epsilon)

        # In binary case, we need to flip the sign of coef, intercept and
        # decision function.
        if len(self.label_) == 2 and self.impl != 'one_class':
            return -dec_func

        return dec_func

    def _validate_for_predict(self, X):
        X = atleast2d_or_csr(X, dtype=np.float64, order="C")
        if self._sparse and not sp.isspmatrix(X):
            X = sp.csr_matrix(X)
        if sp.issparse(X) and not self._sparse:
            raise ValueError(
                "cannot use sparse input in %r trained on dense data"
                % type(self).__name__)
        return X

    @property
    def coef_(self):
        if self.kernel != 'linear':
            raise ValueError('coef_ is only available when using a '
                             'linear kernel')

        if self.dual_coef_.shape[0] == 1:
            # binary classifier
            coef = -safe_sparse_dot(self.dual_coef_, self.support_vectors_)
        else:
            # 1vs1 classifier
            coef = _one_vs_one_coef(self.dual_coef_, self.n_support_,
                    self.support_vectors_)
            if sp.issparse(coef[0]):
                coef = sp.vstack(coef).tocsr()
            else:
                coef = np.vstack(coef)

        # coef_ being a read-only property it's better to mark the value as
        # immutable to avoid hiding potential bugs for the unsuspecting user
        if sp.issparse(coef):
            # sparse matrix do not have global flags
            coef.data.flags.writeable = False
        else:
            # regular dense array
            coef.flags.writeable = False
        return coef


class BaseLibLinear(BaseEstimator):
    """Base for classes binding liblinear (dense and sparse versions)"""

    _solver_type_dict = {
        'PL2_LLR_D0': 0,  # L2 penalty, logistic regression
        'PL2_LL2_D1': 1,  # L2 penalty, L2 loss, dual form
        'PL2_LL2_D0': 2,  # L2 penalty, L2 loss, primal form
        'PL2_LL1_D1': 3,  # L2 penalty, L1 Loss, dual form
        'MC_SVC': 4,      # Multi-class Support Vector Classification
        'PL1_LL2_D0': 5,  # L1 penalty, L2 Loss, primal form
        'PL1_LLR_D0': 6,  # L1 penalty, logistic regression
        'PL2_LLR_D1': 7,  # L2 penalty, logistic regression, dual form
        }

    def __init__(self, penalty='l2', loss='l2', dual=True, tol=1e-4, C=1.0,
                 multi_class=False, fit_intercept=True, intercept_scaling=1,
                 scale_C=True, class_weight=None):
        self.penalty = penalty
        self.loss = loss
        self.dual = dual
        self.tol = tol
        self.C = C
        self.fit_intercept = fit_intercept
        self.intercept_scaling = intercept_scaling
        self.multi_class = multi_class
        self.scale_C = scale_C
        self.class_weight = class_weight

        if not scale_C:
            warnings.warn('SVM: scale_C will disappear and be assumed to be '
                          'True in scikit-learn 0.12', FutureWarning,
                          stacklevel=2)

        # Check that the arguments given are valid:
        self._get_solver_type()

    def _get_solver_type(self):
        """Find the liblinear magic number for the solver.

        This number depends on the values of the following attributes:
          - multi_class
          - penalty
          - loss
          - dual
        """
        if self.multi_class:
            solver_type = 'MC_SVC'
        else:
            solver_type = "P%s_L%s_D%d" % (
                self.penalty.upper(), self.loss.upper(), int(self.dual))
        if not solver_type in self._solver_type_dict:
            if self.penalty.upper() == 'L1' and self.loss.upper() == 'L1':
                error_string = ("The combination of penalty='l1' "
                    "and loss='l1' is not supported.")
            elif self.penalty.upper() == 'L2' and self.loss.upper() == 'L1':
                # this has to be in primal
                error_string = ("loss='l2' and penalty='l1' is "
                    "only supported when dual='true'.")
            else:
                # only PL1 in dual remains
                error_string = ("penalty='l1' is only supported "
                    "when dual='false'.")
            raise ValueError('Not supported set of arguments: '
                             + error_string)
        return self._solver_type_dict[solver_type]

    def fit(self, X, y, class_weight=None):
        """Fit the model according to the given training data.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target vector relative to X

        class_weight : {dict, 'auto'}, optional
            Weights associated with classes. If not given, all classes
            are supposed to have weight one.

        Returns
        -------
        self : object
            Returns self.
        """

        if class_weight != None:
            warnings.warn("'class_weight' is now an initialization parameter."
                    "Using it in the 'fit' method is deprecated.",
                    DeprecationWarning)
            self.class_weight = class_weight

        X = atleast2d_or_csr(X, dtype=np.float64, order="C")
        y = np.asarray(y, dtype=np.int32).ravel()
        self._sparse = sp.isspmatrix(X)

        self.class_weight_, self.class_weight_label_ = \
                     _get_class_weight(self.class_weight, y)

        if X.shape[0] != y.shape[0]:
            raise ValueError("X and y have incompatible shapes.\n" +
                             "X has %s samples, but y has %s." % \
                             (X.shape[0], y.shape[0]))

        C = self.C
        if self.scale_C:
            C = C / float(X.shape[0])

        train = liblinear.csr_train_wrap if self._sparse \
                                         else liblinear.train_wrap
        self.raw_coef_, self.label_ = train(X, y, self._get_solver_type(),
                                            self.tol, self._get_bias(), C,
                                            self.class_weight_label_,
                                            self.class_weight_)

        return self

    def predict(self, X):
        """Predict target values of X according to the fitted model.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]

        Returns
        -------
        C : array, shape = [n_samples]
        """
        X = self._validate_for_predict(X)
        self._check_n_features(X)

        predict = liblinear.csr_predict_wrap if self._sparse \
                                             else liblinear.predict_wrap
        return predict(X, self.raw_coef_, self._get_solver_type(), self.tol,
                       self.C, self.class_weight_label_, self.class_weight_,
                       self.label_, self._get_bias())

    def decision_function(self, X):
        """Decision function value for X according to the trained model.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        T : array-like, shape = [n_samples, n_class]
            Returns the decision function of the sample for each class
            in the model.
        """
        X = self._validate_for_predict(X)
        self._check_n_features(X)

        dfunc_wrap = liblinear.csr_decision_function_wrap \
                       if self._sparse \
                       else liblinear.decision_function_wrap

        dec_func = dfunc_wrap(X, self.raw_coef_, self._get_solver_type(),
                self.tol, self.C, self.class_weight_label_, self.class_weight_,
                self.label_, self._get_bias())

        return dec_func

    def _check_n_features(self, X):
        n_features = self.raw_coef_.shape[1]
        if self.fit_intercept:
            n_features -= 1
        if X.shape[1] != n_features:
            raise ValueError("X.shape[1] should be %d, not %d." % (n_features,
                                                                   X.shape[1]))

    def _validate_for_predict(self, X):
        X = atleast2d_or_csr(X, dtype=np.float64, order="C")
        if self._sparse and not sp.isspmatrix(X):
            X = sp.csr_matrix(X)
        elif sp.isspmatrix(X) and not self._sparse:
            raise ValueError(
                "cannot use sparse input in %r trained on dense data"
                % type(self).__name__)
        return X

    def _get_intercept_(self):
        if self.fit_intercept:
            ret = self.intercept_scaling * self.raw_coef_[:, -1]
            return ret
        return 0.0

    def _set_intercept_(self, intercept):
        self.fit_intercept = True

        intercept /= self.intercept_scaling
        intercept = intercept.reshape(-1, 1)

        self.raw_coef_ = np.hstack((self.raw_coef_[:, : -1], intercept))
        # We need fortran ordered arrays for the predict
        self.raw_coef_ = np.asfortranarray(self.raw_coef_)

    intercept_ = property(_get_intercept_, _set_intercept_)

    def _get_coef_(self):
        if self.fit_intercept:
            ret = self.raw_coef_[:, : -1].copy()
        else:
            ret = self.raw_coef_.copy()

        # mark the returned value as immutable
        # to avoid silencing potential bugs
        ret.flags.writeable = False
        return ret

    def _set_coef_(self, coef):
        raw_intercept = self.raw_coef_[:, -1].reshape(-1, 1)

        self.raw_coef_ = coef

        if self.fit_intercept:
            self.raw_coef_ = np.hstack((self.raw_coef_, raw_intercept))

        # We need fortran ordered arrays for the predict
        self.raw_coef_ = np.asfortranarray(self.raw_coef_)

    coef_ = property(_get_coef_, _set_coef_)

    def _get_bias(self):
        if self.fit_intercept:
            return self.intercept_scaling
        else:
            return -1.0


libsvm.set_verbosity_wrap(0)
libsvm_sparse.set_verbosity_wrap(0)
