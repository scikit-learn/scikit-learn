import numpy as np
import scipy.sparse as sp
import warnings
from abc import ABCMeta, abstractmethod

from . import libsvm, liblinear
from . import libsvm_sparse
from ..base import BaseEstimator, ClassifierMixin
from ..utils import atleast2d_or_csr, array2d
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

    __metaclass__ = ABCMeta
    _sparse_kernels = ["linear", "poly", "rbf", "sigmoid", "precomputed"]

    @abstractmethod
    def __init__(self, impl, kernel, degree, gamma, coef0,
                 tol, C, nu, epsilon, shrinking, probability, cache_size,
                 sparse, class_weight, verbose):

        if not impl in LIBSVM_IMPL:
            raise ValueError("impl should be one of %s, %s was given" % (
                LIBSVM_IMPL, impl))

        if C is None:
            warnings.warn("Using 'None' for C of BaseLibSVM is deprecated "
                    "since version 0.12, and backward compatibility "
                    "won't be maintained from version 0.14 onward. "
                    "Setting C=1.0.", DeprecationWarning, stacklevel=2)
            C = 1.0

        self.impl = impl
        self.kernel = kernel
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
        self.sparse = sparse
        self.class_weight = class_weight
        self.verbose = verbose

    @property
    def _pairwise(self):
        kernel = self.kernel
        return kernel == "precomputed" or hasattr(kernel, "__call__")

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

        if self.sparse == "auto":
            self._sparse = sp.isspmatrix(X) and not self._pairwise
        else:
            self._sparse = self.sparse

        if self._sparse and self._pairwise:
            raise ValueError("Sparse precomputed kernels are not supported. "
                    "Using sparse data and dense kernels is possible by not "
                    "using the ``sparse`` parameter")

        X = atleast2d_or_csr(X, dtype=np.float64, order='C')
        y = np.asarray(y, dtype=np.float64, order='C')

        if self.impl != "one_class" and len(np.unique(y)) < 2:
            raise ValueError("The number of classes has to be greater than"
                    " one.")

        if class_weight != None:
            warnings.warn("'class_weight' is now an initialization parameter."
                          "Using it in the 'fit' method is deprecated and "
                          "will be removed in 0.13.", DeprecationWarning,
                          stacklevel=2)
            self.class_weight = class_weight

        sample_weight = np.asarray([] if sample_weight is None
                                      else sample_weight, dtype=np.float64)
        solver_type = LIBSVM_IMPL.index(self.impl)
        self.class_weight_, self.class_weight_label_ = \
                     _get_class_weight(self.class_weight, y)

        # input validation
        if solver_type != 2 and X.shape[0] != y.shape[0]:
            raise ValueError("X and y have incompatible shapes.\n" +
                             "X has %s samples, but y has %s." %
                             (X.shape[0], y.shape[0]))

        if self.kernel == "precomputed" and X.shape[0] != X.shape[1]:
            raise ValueError("X.shape[0] should be equal to X.shape[1]")

        if sample_weight.shape[0] > 0 and sample_weight.shape[0] != X.shape[0]:
            raise ValueError("sample_weight and X have incompatible shapes:"
                             "%r vs %r\n"
                             "Note: Sparse matrices cannot be indexed w/"
                             "boolean masks (use `indices=True` in CV)."
                             % (sample_weight.shape, X.shape))

        if (self.kernel in ['poly', 'rbf']) and (self.gamma == 0):
            # if custom gamma is not provided ...
            self._gamma = 1.0 / X.shape[1]
        else:
            self._gamma = self.gamma

        kernel = self.kernel
        if hasattr(kernel, '__call__'):
            kernel = 'precomputed'

        fit = self._sparse_fit if self._sparse else self._dense_fit
        if self.verbose:
            print '[LibSVM]',
        fit(X, y, sample_weight, solver_type, kernel)

        self.shape_fit_ = X.shape

        # In binary case, we need to flip the sign of coef, intercept and
        # decision function. Use self._intercept_ internally.
        self._intercept_ = self.intercept_.copy()
        if len(self.label_) == 2 and self.impl != 'one_class':
            self.intercept_ *= -1
        return self

    def _dense_fit(self, X, y, sample_weight, solver_type, kernel):

        if hasattr(self.kernel, '__call__'):
            # you must store a reference to X to compute the kernel in predict
            # TODO: add keyword copy to copy on demand
            self.__Xfit = X
            X = self._compute_kernel(X)

        if hasattr(self.kernel, '__call__') and X.shape[0] != X.shape[1]:
            raise ValueError("X.shape[0] should be equal to X.shape[1]")

        libsvm.set_verbosity_wrap(self.verbose)

        # we don't pass **self.get_params() to allow subclasses to
        # add other parameters to __init__
        self.support_, self.support_vectors_, self.n_support_, \
        self.dual_coef_, self.intercept_, self.label_, self.probA_, \
        self.probB_ = libsvm.fit(X, y,
            svm_type=solver_type, sample_weight=sample_weight,
            class_weight=self.class_weight_,
            class_weight_label=self.class_weight_label_,
            kernel=kernel, C=self.C, nu=self.nu,
            probability=self.probability, degree=self.degree,
            shrinking=self.shrinking, tol=self.tol, cache_size=self.cache_size,
            coef0=self.coef0, gamma=self._gamma, epsilon=self.epsilon)

    def _sparse_fit(self, X, y, sample_weight, solver_type, kernel):
        X.data = np.asarray(X.data, dtype=np.float64, order='C')

        kernel_type = self._sparse_kernels.index(kernel)

        libsvm_sparse.set_verbosity_wrap(self.verbose)

        self.support_vectors_, dual_coef_data, self.intercept_, self.label_, \
            self.n_support_, self.probA_, self.probB_ = \
            libsvm_sparse.libsvm_sparse_train(
                 X.shape[1], X.data, X.indices, X.indptr, y, solver_type,
                 kernel_type, self.degree, self._gamma, self.coef0, self.tol,
                 self.C, self.class_weight_label_, self.class_weight_,
                 sample_weight, self.nu, self.cache_size, self.epsilon,
                 int(self.shrinking), int(self.probability))

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
        y_pred : array, shape = [n_samples]
        """
        X = self._validate_for_predict(X)
        predict = self._sparse_predict if self._sparse else self._dense_predict
        return predict(X)

    def _dense_predict(self, X):
        n_samples, n_features = X.shape
        X = self._compute_kernel(X)
        if X.ndim == 1:
            X = array2d(X, order='C')

        kernel = self.kernel
        if hasattr(self.kernel, "__call__"):
            kernel = 'precomputed'
            if X.shape[1] != self.shape_fit_[0]:
                raise ValueError("X.shape[1] = %d should be equal to %d, "
                                 "the number of samples at training time" %
                                 (X.shape[1], self.shape_fit_[0]))

        C = 0.0  # C is not useful here

        svm_type = LIBSVM_IMPL.index(self.impl)

        return libsvm.predict(
            X, self.support_, self.support_vectors_, self.n_support_,
            self.dual_coef_, self._intercept_,
            self.label_, self.probA_, self.probB_,
            svm_type=svm_type,
            kernel=kernel, C=C, nu=self.nu,
            probability=self.probability, degree=self.degree,
            shrinking=self.shrinking, tol=self.tol, cache_size=self.cache_size,
            coef0=self.coef0, gamma=self._gamma, epsilon=self.epsilon)

    def _sparse_predict(self, X):
        X = sp.csr_matrix(X, dtype=np.float64)

        kernel = self.kernel
        if hasattr(kernel, '__call__'):
            kernel = 'precomputed'

        kernel_type = self._sparse_kernels.index(kernel)

        C = 0.0  # C is not useful here

        return libsvm_sparse.libsvm_sparse_predict(
                      X.data, X.indices, X.indptr,
                      self.support_vectors_.data,
                      self.support_vectors_.indices,
                      self.support_vectors_.indptr,
                      self.dual_coef_.data, self._intercept_,
                      LIBSVM_IMPL.index(self.impl), kernel_type,
                      self.degree, self._gamma, self.coef0, self.tol,
                      C, self.class_weight_label_, self.class_weight_,
                      self.nu, self.epsilon, self.shrinking,
                      self.probability, self.n_support_, self.label_,
                      self.probA_, self.probB_)

    def _compute_kernel(self, X):
        """Return the data transformed by a callable kernel"""
        if hasattr(self.kernel, '__call__'):
            # in the case of precomputed kernel given as a function, we
            # have to compute explicitly the kernel matrix
            kernel = self.kernel(X, self.__Xfit)
            if sp.issparse(kernel):
                kernel = kernel.toarray()
            X = np.asarray(kernel, dtype=np.float64, order='C')
        return X

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
            raise NotImplementedError("decision_function not supported for"
                    " sparse SVM")

        X = self._validate_for_predict(X)

        C = 0.0  # C is not useful here

        kernel = self.kernel
        if hasattr(kernel, '__call__'):
            kernel = 'precomputed'

        dec_func = libsvm.decision_function(
            X, self.support_, self.support_vectors_, self.n_support_,
            self.dual_coef_, self._intercept_, self.label_,
            self.probA_, self.probB_,
            svm_type=LIBSVM_IMPL.index(self.impl),
            kernel=kernel, C=C, nu=self.nu,
            probability=self.probability, degree=self.degree,
            shrinking=self.shrinking, tol=self.tol, cache_size=self.cache_size,
            coef0=self.coef0, gamma=self._gamma, epsilon=self.epsilon)

        # In binary case, we need to flip the sign of coef, intercept and
        # decision function.
        if len(self.label_) == 2 and self.impl != 'one_class':
            return -dec_func

        return dec_func

    def _validate_for_predict(self, X):
        X = atleast2d_or_csr(X, dtype=np.float64, order="C")
        if self._sparse and not sp.isspmatrix(X):
            X = sp.csr_matrix(X)
        if (sp.issparse(X) and not self._sparse and
                not hasattr(self.kernel, '__call__')):
            raise ValueError(
                "cannot use sparse input in %r trained on dense data"
                % type(self).__name__)
        n_samples, n_features = X.shape

        if self.kernel == "precomputed":
            if X.shape[1] != self.shape_fit_[0]:
                raise ValueError("X.shape[1] = %d should be equal to %d, "
                                 "the number of samples at training time" %
                                 (X.shape[1], self.shape_fit_[0]))
        elif n_features != self.shape_fit_[1]:
            raise ValueError("X.shape[1] = %d should be equal to %d, "
                             "the number of features at training time" %
                             (n_features, self.shape_fit_[1]))
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


class BaseSVC(BaseLibSVM, ClassifierMixin):
    """ABC for LibSVM-based classifiers."""

    def predict_proba(self, X):
        """Compute probabilities of possible outcomes for samples in X.

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
        predict. Also, it will produce meaningless results on very small
        datasets.
        """
        if not self.probability:
            raise NotImplementedError(
                    "probability estimates must be enabled to use this method")

        if self.impl not in ('c_svc', 'nu_svc'):
            raise NotImplementedError("predict_proba only implemented for SVC "
                                      "and NuSVC")

        X = self._validate_for_predict(X)
        pred_proba = self._sparse_predict_proba if self._sparse \
                                                else self._dense_predict_proba
        return pred_proba(X)

    def predict_log_proba(self, X):
        """Compute log probabilities of possible outcomes for samples in X.

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
        predict. Also, it will produce meaningless results on very small
        datasets.
        """
        return np.log(self.predict_proba(X))

    def _dense_predict_proba(self, X):
        X = self._compute_kernel(X)

        C = 0.0  # C is not useful here

        kernel = self.kernel
        if hasattr(kernel, '__call__'):
            kernel = 'precomputed'

        svm_type = LIBSVM_IMPL.index(self.impl)
        pprob = libsvm.predict_proba(
            X, self.support_, self.support_vectors_, self.n_support_,
            self.dual_coef_, self._intercept_, self.label_,
            self.probA_, self.probB_,
            svm_type=svm_type, kernel=kernel, C=C, nu=self.nu,
            probability=self.probability, degree=self.degree,
            shrinking=self.shrinking, tol=self.tol, cache_size=self.cache_size,
            coef0=self.coef0, gamma=self._gamma, epsilon=self.epsilon)

        return pprob

    def _sparse_predict_proba(self, X):
        X.data = np.asarray(X.data, dtype=np.float64, order='C')

        kernel = self.kernel
        if hasattr(kernel, '__call__'):
            kernel = 'precomputed'

        kernel_type = self._sparse_kernels.index(kernel)

        return libsvm_sparse.libsvm_sparse_predict_proba(
            X.data, X.indices, X.indptr,
            self.support_vectors_.data,
            self.support_vectors_.indices,
            self.support_vectors_.indptr,
            self.dual_coef_.data, self._intercept_,
            LIBSVM_IMPL.index(self.impl), kernel_type,
            self.degree, self._gamma, self.coef0, self.tol,
            self.C, self.class_weight_label_, self.class_weight_,
            self.nu, self.epsilon, self.shrinking,
            self.probability, self.n_support_, self.label_,
            self.probA_, self.probB_)


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
            multi_class='ovr', fit_intercept=True, intercept_scaling=1,
            class_weight=None, verbose=0):

        if C is None:
            warnings.warn("Using 'None' for C of BaseLibLinear is deprecated "
                    "since version 0.12, and backward compatibility "
                    "won't be maintained from version 0.14 onward. "
                    "Setting C=1.0.", DeprecationWarning, stacklevel=2)
            C = 1.0

        self.penalty = penalty
        self.loss = loss
        self.dual = dual
        self.tol = tol
        self.C = C
        self.fit_intercept = fit_intercept
        self.intercept_scaling = intercept_scaling
        self.multi_class = multi_class
        self.class_weight = class_weight
        self.verbose = verbose

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
        if self.multi_class == 'crammer_singer':
            solver_type = 'MC_SVC'
        else:
            if self.multi_class != 'ovr':
                raise ValueError("`multi_class` must be one of `ovr`, "
                        "`crammer_singer`")
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
        if len(np.unique(y)) < 2:
            raise ValueError("The number of classes has to be greater than"
                    " one.")

        if class_weight != None:
            warnings.warn("'class_weight' is now an initialization parameter."
                          "Using it in the 'fit' method is deprecated and "
                          "will be removed in 0.13.", DeprecationWarning,
                          stacklevel=2)
            self.class_weight = class_weight

        X = atleast2d_or_csr(X, dtype=np.float64, order="C")
        y = np.asarray(y, dtype=np.float64).ravel()
        self._sparse = sp.isspmatrix(X)

        self.class_weight_, self.class_weight_label_ = \
                     _get_class_weight(self.class_weight, y)

        if X.shape[0] != y.shape[0]:
            raise ValueError("X and y have incompatible shapes.\n" +
                             "X has %s samples, but y has %s." % \
                             (X.shape[0], y.shape[0]))

        liblinear.set_verbosity_wrap(self.verbose)

        if self._sparse:
            train = liblinear.csr_train_wrap
        else:
            train = liblinear.train_wrap

        if self.verbose:
            print '[LibLinear]',
        self.raw_coef_, self.label_ = train(X, y, self._get_solver_type(),
                                            self.tol, self._get_bias(), self.C,
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

        C = 0.0  # C is not useful here

        predict = liblinear.csr_predict_wrap if self._sparse \
                                             else liblinear.predict_wrap
        return predict(X, self.raw_coef_, self._get_solver_type(), self.tol,
                       C, self.class_weight_label_, self.class_weight_,
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

        C = 0.0  # C is not useful here

        dfunc_wrap = liblinear.csr_decision_function_wrap \
                       if self._sparse \
                       else liblinear.decision_function_wrap

        dec_func = dfunc_wrap(X, self.raw_coef_, self._get_solver_type(),
                self.tol, C, self.class_weight_label_, self.class_weight_,
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
        if not self.raw_coef_.flags['F_CONTIGUOUS']:
            warnings.warn('Coefficients are the fortran-contiguous. '
                          'Copying them.', RuntimeWarning,
                          stacklevel=3)
            self.raw_coef_ = np.asfortranarray(self.raw_coef_)
        self._check_n_features(X)
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
liblinear.set_verbosity_wrap(0)
