from __future__ import print_function

import numpy as np
import scipy.sparse as sp
import warnings
from abc import ABCMeta, abstractmethod

from . import libsvm, liblinear
from . import libsvm_sparse
from ..base import BaseEstimator, ClassifierMixin
from ..preprocessing import LabelEncoder
from ..utils import atleast2d_or_csr, array2d, check_random_state, column_or_1d
from ..utils import ConvergenceWarning, compute_class_weight
from ..utils.extmath import safe_sparse_dot
from ..externals import six


LIBSVM_IMPL = ['c_svc', 'nu_svc', 'one_class', 'epsilon_svr', 'nu_svr']


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
    for class1 in range(n_class):
        # SVs for class1:
        sv1 = support_vectors[sv_locs[class1]:sv_locs[class1 + 1], :]
        for class2 in range(class1 + 1, n_class):
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


class BaseLibSVM(six.with_metaclass(ABCMeta, BaseEstimator)):
    """Base class for estimators that use libsvm as backing library

    This implements support vector machine classification and regression.

    Parameter documentation is in the derived `SVC` class.
    """

    # The order of these must match the integer values in LibSVM.
    # XXX These are actually the same in the dense case. Need to factor
    # this out.
    _sparse_kernels = ["linear", "poly", "rbf", "sigmoid", "precomputed"]

    @abstractmethod
    def __init__(self, impl, kernel, degree, gamma, coef0,
                 tol, C, nu, epsilon, shrinking, probability, cache_size,
                 class_weight, verbose, max_iter, random_state):

        if not impl in LIBSVM_IMPL:  # pragma: no cover
            raise ValueError("impl should be one of %s, %s was given" % (
                LIBSVM_IMPL, impl))

        self._impl = impl
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
        self.class_weight = class_weight
        self.verbose = verbose
        self.max_iter = max_iter
        self.random_state = random_state

    @property
    def _pairwise(self):
        # Used by cross_val_score.
        kernel = self.kernel
        return kernel == "precomputed" or callable(kernel)

    def fit(self, X, y, sample_weight=None):
        """Fit the SVM model according to the given training data.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape (n_samples,)
            Target values (class labels in classification, real numbers in
            regression)

        sample_weight : array-like, shape (n_samples,)
            Per-sample weights. Rescale C per sample. Higher weights
            force the classifier to put more emphasis on these points.

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

        rnd = check_random_state(self.random_state)

        sparse = sp.isspmatrix(X)
        if sparse and self.kernel == "precomputed":
            raise TypeError("Sparse precomputed kernels are not supported.")
        self._sparse = sparse and not callable(self.kernel)

        X = atleast2d_or_csr(X, dtype=np.float64, order='C')
        y = self._validate_targets(y)

        sample_weight = np.asarray([]
                                   if sample_weight is None
                                   else sample_weight, dtype=np.float64)
        solver_type = LIBSVM_IMPL.index(self._impl)

        # input validation
        if solver_type != 2 and X.shape[0] != y.shape[0]:
            raise ValueError("X and y have incompatible shapes.\n" +
                             "X has %s samples, but y has %s." %
                             (X.shape[0], y.shape[0]))

        if self.kernel == "precomputed" and X.shape[0] != X.shape[1]:
            raise ValueError("X.shape[0] should be equal to X.shape[1]")

        if sample_weight.shape[0] > 0 and sample_weight.shape[0] != X.shape[0]:
            raise ValueError("sample_weight and X have incompatible shapes: "
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
        if callable(kernel):
            kernel = 'precomputed'

        fit = self._sparse_fit if self._sparse else self._dense_fit
        if self.verbose:  # pragma: no cover
            print('[LibSVM]', end='')

        seed = rnd.randint(np.iinfo('i').max)
        fit(X, y, sample_weight, solver_type, kernel, random_seed=seed)
        # see comment on the other call to np.iinfo in this file

        self.shape_fit_ = X.shape

        # In binary case, we need to flip the sign of coef, intercept and
        # decision function. Use self._intercept_ internally.
        self._intercept_ = self.intercept_.copy()
        if self._impl in ['c_svc', 'nu_svc'] and len(self.classes_) == 2:
            self.intercept_ *= -1
        return self

    def _validate_targets(self, y):
        """Validation of y and class_weight.

        Default implementation for SVR and one-class; overridden in BaseSVC.
        """
        # XXX this is ugly.
        # Regression models should not have a class_weight_ attribute.
        self.class_weight_ = np.empty(0)
        return np.asarray(y, dtype=np.float64, order='C')

    def _warn_from_fit_status(self):
        assert self.fit_status_ in (0, 1)
        if self.fit_status_ == 1:
            warnings.warn('Solver terminated early (max_iter=%i).'
                          '  Consider pre-processing your data with'
                          ' StandardScaler or MinMaxScaler.'
                          % self.max_iter, ConvergenceWarning)

    def _dense_fit(self, X, y, sample_weight, solver_type, kernel,
                   random_seed):
        if callable(self.kernel):
            # you must store a reference to X to compute the kernel in predict
            # TODO: add keyword copy to copy on demand
            self.__Xfit = X
            X = self._compute_kernel(X)

            if X.shape[0] != X.shape[1]:
                raise ValueError("X.shape[0] should be equal to X.shape[1]")

        libsvm.set_verbosity_wrap(self.verbose)

        # we don't pass **self.get_params() to allow subclasses to
        # add other parameters to __init__
        self.support_, self.support_vectors_, self.n_support_, \
            self.dual_coef_, self.intercept_, self.probA_, \
            self.probB_, self.fit_status_ = libsvm.fit(
                X, y,
                svm_type=solver_type, sample_weight=sample_weight,
                class_weight=self.class_weight_, kernel=kernel, C=self.C,
                nu=self.nu, probability=self.probability, degree=self.degree,
                shrinking=self.shrinking, tol=self.tol,
                cache_size=self.cache_size, coef0=self.coef0,
                gamma=self._gamma, epsilon=self.epsilon,
                max_iter=self.max_iter, random_seed=random_seed)

        self._warn_from_fit_status()

    def _sparse_fit(self, X, y, sample_weight, solver_type, kernel,
                    random_seed):
        X.data = np.asarray(X.data, dtype=np.float64, order='C')
        X.sort_indices()

        kernel_type = self._sparse_kernels.index(kernel)

        libsvm_sparse.set_verbosity_wrap(self.verbose)

        self.support_, self.support_vectors_, dual_coef_data, \
            self.intercept_, self.n_support_, \
            self.probA_, self.probB_, self.fit_status_ = \
            libsvm_sparse.libsvm_sparse_train(
                X.shape[1], X.data, X.indices, X.indptr, y, solver_type,
                kernel_type, self.degree, self._gamma, self.coef0, self.tol,
                self.C, self.class_weight_,
                sample_weight, self.nu, self.cache_size, self.epsilon,
                int(self.shrinking), int(self.probability), self.max_iter,
                random_seed)

        self._warn_from_fit_status()

        if hasattr(self, "classes_"):
            n_class = len(self.classes_) - 1
        else:   # regression
            n_class = 1
        n_SV = self.support_vectors_.shape[0]

        dual_coef_indices = np.tile(np.arange(n_SV), n_class)
        dual_coef_indptr = np.arange(0, dual_coef_indices.size + 1,
                                     dual_coef_indices.size / n_class)
        self.dual_coef_ = sp.csr_matrix(
            (dual_coef_data, dual_coef_indices, dual_coef_indptr),
            (n_class, n_SV))

    def predict(self, X):
        """Perform regression on samples in X.

        For an one-class model, +1 or -1 is returned.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)

        Returns
        -------
        y_pred : array, shape (n_samples,)
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
        if callable(self.kernel):
            kernel = 'precomputed'
            if X.shape[1] != self.shape_fit_[0]:
                raise ValueError("X.shape[1] = %d should be equal to %d, "
                                 "the number of samples at training time" %
                                 (X.shape[1], self.shape_fit_[0]))

        svm_type = LIBSVM_IMPL.index(self._impl)

        return libsvm.predict(
            X, self.support_, self.support_vectors_, self.n_support_,
            self.dual_coef_, self._intercept_,
            self.probA_, self.probB_, svm_type=svm_type, kernel=kernel,
            degree=self.degree, coef0=self.coef0, gamma=self._gamma,
            cache_size=self.cache_size)

    def _sparse_predict(self, X):
        X = sp.csr_matrix(X, dtype=np.float64)

        kernel = self.kernel
        if callable(kernel):
            kernel = 'precomputed'

        kernel_type = self._sparse_kernels.index(kernel)

        C = 0.0  # C is not useful here

        return libsvm_sparse.libsvm_sparse_predict(
            X.data, X.indices, X.indptr,
            self.support_vectors_.data,
            self.support_vectors_.indices,
            self.support_vectors_.indptr,
            self.dual_coef_.data, self._intercept_,
            LIBSVM_IMPL.index(self._impl), kernel_type,
            self.degree, self._gamma, self.coef0, self.tol,
            C, self.class_weight_,
            self.nu, self.epsilon, self.shrinking,
            self.probability, self.n_support_,
            self.probA_, self.probB_)

    def _compute_kernel(self, X):
        """Return the data transformed by a callable kernel"""
        if callable(self.kernel):
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
            raise NotImplementedError("Decision_function not supported for"
                                      " sparse SVM.")

        X = self._validate_for_predict(X)
        X = self._compute_kernel(X)

        kernel = self.kernel
        if callable(kernel):
            kernel = 'precomputed'

        dec_func = libsvm.decision_function(
            X, self.support_, self.support_vectors_, self.n_support_,
            self.dual_coef_, self._intercept_,
            self.probA_, self.probB_,
            svm_type=LIBSVM_IMPL.index(self._impl),
            kernel=kernel, degree=self.degree, cache_size=self.cache_size,
            coef0=self.coef0, gamma=self._gamma)

        # In binary case, we need to flip the sign of coef, intercept and
        # decision function.
        if self._impl in ['c_svc', 'nu_svc'] and len(self.classes_) == 2:
            return -dec_func

        return dec_func

    def _validate_for_predict(self, X):
        X = atleast2d_or_csr(X, dtype=np.float64, order="C")
        if self._sparse and not sp.isspmatrix(X):
            X = sp.csr_matrix(X)
        if self._sparse:
            X.sort_indices()

        if sp.issparse(X) and not self._sparse and not callable(self.kernel):
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

    def _validate_targets(self, y):
        y_ = column_or_1d(y, warn=True)
        cls, y = np.unique(y_, return_inverse=True)
        self.class_weight_ = compute_class_weight(self.class_weight, cls, y_)
        if len(cls) < 2:
            raise ValueError(
                "The number of classes has to be greater than one; got %d"
                % len(cls))

        self.classes_ = cls

        return np.asarray(y, dtype=np.float64, order='C')

    def predict(self, X):
        """Perform classification on samples in X.

        For an one-class model, +1 or -1 is returned.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]

        Returns
        -------
        y_pred : array, shape = [n_samples]
            Class labels for samples in X.
        """
        y = super(BaseSVC, self).predict(X)
        return self.classes_.take(np.asarray(y, dtype=np.intp))

    # Hacky way of getting predict_proba to raise an AttributeError when
    # probability=False using properties. Do not use this in new code; when
    # probabilities are not available depending on a setting, introduce two
    # estimators.
    def _check_proba(self):
        if not self.probability:
            raise AttributeError("predict_proba is not available when"
                                 " probability=%r" % self.probability)
        if self._impl not in ('c_svc', 'nu_svc'):
            raise AttributeError("predict_proba only implemented for SVC"
                                 " and NuSVC")

    @property
    def predict_proba(self):
        """Compute probabilities of possible outcomes for samples in X.

        The model need to have probability information computed at training
        time: fit with attribute `probability` set to True.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        T : array-like, shape = [n_samples, n_classes]
            Returns the probability of the sample for each class in
            the model. The columns correspond to the classes in sorted
            order, as they appear in the attribute `classes_`.

        Notes
        -----
        The probability model is created using cross validation, so
        the results can be slightly different than those obtained by
        predict. Also, it will produce meaningless results on very small
        datasets.
        """
        self._check_proba()
        return self._predict_proba

    def _predict_proba(self, X):
        X = self._validate_for_predict(X)
        pred_proba = (self._sparse_predict_proba
                      if self._sparse else self._dense_predict_proba)
        return pred_proba(X)

    @property
    def predict_log_proba(self):
        """Compute log probabilities of possible outcomes for samples in X.

        The model need to have probability information computed at training
        time: fit with attribute `probability` set to True.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        T : array-like, shape = [n_samples, n_classes]
            Returns the log-probabilities of the sample for each class in
            the model. The columns correspond to the classes in sorted
            order, as they appear in the attribute `classes_`.

        Notes
        -----
        The probability model is created using cross validation, so
        the results can be slightly different than those obtained by
        predict. Also, it will produce meaningless results on very small
        datasets.
        """
        self._check_proba()
        return self._predict_log_proba

    def _predict_log_proba(self, X):
        return np.log(self.predict_proba(X))

    def _dense_predict_proba(self, X):
        X = self._compute_kernel(X)

        kernel = self.kernel
        if callable(kernel):
            kernel = 'precomputed'

        svm_type = LIBSVM_IMPL.index(self._impl)
        pprob = libsvm.predict_proba(
            X, self.support_, self.support_vectors_, self.n_support_,
            self.dual_coef_, self._intercept_,
            self.probA_, self.probB_,
            svm_type=svm_type, kernel=kernel, degree=self.degree,
            cache_size=self.cache_size, coef0=self.coef0, gamma=self._gamma)

        return pprob

    def _sparse_predict_proba(self, X):
        X.data = np.asarray(X.data, dtype=np.float64, order='C')

        kernel = self.kernel
        if callable(kernel):
            kernel = 'precomputed'

        kernel_type = self._sparse_kernels.index(kernel)

        return libsvm_sparse.libsvm_sparse_predict_proba(
            X.data, X.indices, X.indptr,
            self.support_vectors_.data,
            self.support_vectors_.indices,
            self.support_vectors_.indptr,
            self.dual_coef_.data, self._intercept_,
            LIBSVM_IMPL.index(self._impl), kernel_type,
            self.degree, self._gamma, self.coef0, self.tol,
            self.C, self.class_weight_,
            self.nu, self.epsilon, self.shrinking,
            self.probability, self.n_support_,
            self.probA_, self.probB_)


class BaseLibLinear(six.with_metaclass(ABCMeta, BaseEstimator)):
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

    @abstractmethod
    def __init__(self, penalty='l2', loss='l2', dual=True, tol=1e-4, C=1.0,
                 multi_class='ovr', fit_intercept=True, intercept_scaling=1,
                 class_weight=None, verbose=0, random_state=None):

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
        self.random_state = random_state

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
        elif self.multi_class == 'ovr':
            solver_type = "P%s_L%s_D%d" % (
                self.penalty.upper(), self.loss.upper(), int(self.dual))
        else:
            raise ValueError("`multi_class` must be one of `ovr`, "
                             "`crammer_singer`, got %r" % self.multi_class)
        if not solver_type in self._solver_type_dict:
            if self.penalty.upper() == 'L1' and self.loss.upper() == 'L1':
                error_string = ("The combination of penalty='l1' "
                                "and loss='l1' is not supported.")
            elif self.penalty.upper() == 'L2' and self.loss.upper() == 'L1':
                # this has to be in primal
                error_string = ("penalty='l2' and loss='l1' is "
                                "only supported when dual='true'.")
            else:
                # only PL1 in dual remains
                error_string = ("penalty='l1' is only supported "
                                "when dual='false'.")
            raise ValueError('Unsupported set of arguments: %s, '
                             'Parameters: penalty=%r, loss=%r, dual=%r'
                             % (error_string, self.penalty,
                                self.loss, self.dual))
        return self._solver_type_dict[solver_type]

    def fit(self, X, y):
        """Fit the model according to the given training data.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target vector relative to X

        Returns
        -------
        self : object
            Returns self.
        """
        self._enc = LabelEncoder()
        y_ind = self._enc.fit_transform(y)
        if len(self.classes_) < 2:
            raise ValueError("The number of classes has to be greater than"
                             " one.")

        X = atleast2d_or_csr(X, dtype=np.float64, order="C")

        self.class_weight_ = compute_class_weight(self.class_weight,
                                                  self.classes_, y)

        if X.shape[0] != y_ind.shape[0]:
            raise ValueError("X and y have incompatible shapes.\n"
                             "X has %s samples, but y has %s." %
                             (X.shape[0], y_ind.shape[0]))

        liblinear.set_verbosity_wrap(self.verbose)

        rnd = check_random_state(self.random_state)
        if self.verbose:
            print('[LibLinear]', end='')

        # LibLinear wants targets as doubles, even for classification
        y_ind = np.asarray(y_ind, dtype=np.float64).ravel()
        self.raw_coef_ = liblinear.train_wrap(X, y_ind,
                                              sp.isspmatrix(X),
                                              self._get_solver_type(),
                                              self.tol, self._get_bias(),
                                              self.C,
                                              self.class_weight_,
                                              rnd.randint(np.iinfo('i').max))
        # Regarding rnd.randint(..) in the above signature:
        # seed for srand in range [0..INT_MAX); due to limitations in Numpy
        # on 32-bit platforms, we can't get to the UINT_MAX limit that
        # srand supports

        if self.fit_intercept:
            self.coef_ = self.raw_coef_[:, :-1]
            self.intercept_ = self.intercept_scaling * self.raw_coef_[:, -1]
        else:
            self.coef_ = self.raw_coef_
            self.intercept_ = 0.

        if self.multi_class == "crammer_singer" and len(self.classes_) == 2:
            self.coef_ = (self.coef_[1] - self.coef_[0]).reshape(1, -1)
            if self.fit_intercept:
                intercept = self.intercept_[1] - self.intercept_[0]
                self.intercept_ = np.array([intercept])

        return self

    @property
    def classes_(self):
        return self._enc.classes_

    def _get_bias(self):
        if self.fit_intercept:
            return self.intercept_scaling
        else:
            return -1.0


libsvm.set_verbosity_wrap(0)
libsvm_sparse.set_verbosity_wrap(0)
liblinear.set_verbosity_wrap(0)
