from abc import ABCMeta, abstractmethod

import numpy as np
import scipy.sparse as sp

from . import libsvm, liblinear
from ..base import BaseEstimator
from ..utils import array2d
from ..utils import safe_asarray
from ..utils.extmath import safe_sparse_dot
import warnings


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


class BaseLibSVM(BaseEstimator):
    """Base class for estimators that use libsvm as backing library

    This implements support vector machine classification and regression.
    """

    __metaclass__ = ABCMeta

    def __init__(self, impl, kernel, degree, gamma, coef0,
                 tol, C, nu, epsilon, shrinking, probability, cache_size,
                 scale_C):

        if not impl in LIBSVM_IMPL:
            raise ValueError("impl should be one of %s, %s was given" % (
                LIBSVM_IMPL, impl))
        if hasattr(kernel, '__call__'):
            self.kernel_function = kernel
            self.kernel = 'precomputed'
        else:
            self.kernel = kernel
        if not scale_C:
            warnings.warn('SVM: scale_C will be True by default in '
                          'scikit-learn 0.11', FutureWarning,
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

    @abstractmethod
    def fit(self, X, y, class_weight=None, sample_weight=None):
        pass

    @abstractmethod
    def predict(self, X):
        pass

    @abstractmethod
    def predict_proba(self, X):
        pass

    def predict_log_proba(self, T):
        """Compute the log likehoods each possible outcomes of samples in T.

        The model need to have probability information computed at training
        time: fit with attribute `probability` set to True.

        Parameters
        ----------
        T : array-like, shape = [n_samples, n_features]

        Returns
        -------
        T : array-like, shape = [n_samples, n_classes]
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
        return np.log(self.predict_proba(T))

    @property
    def coef_(self):
        if self.kernel != 'linear':
            raise ValueError('coef_ is only available when using a '
                             'linear kernel')
        coef = safe_sparse_dot(self.dual_coef_, self.support_vectors_)
        # coef_ being a read-only property it's better to mark the value as
        # immutable to avoid hiding potential bugs for the unsuspecting user
        if sp.issparse(coef):
            # sparse matrix do not have global flags
            coef.data.flags.writeable = False
        else:
            # regular dense array
            coef.flags.writeable = False
        return coef


class DenseBaseLibSVM(BaseLibSVM):
    def _compute_kernel(self, X):
        """Return the data transformed by a callable kernel"""
        if hasattr(self, 'kernel_function'):
            # in the case of precomputed kernel given as a function, we
            # have to compute explicitly the kernel matrix
            X = np.asarray(self.kernel_function(X, self.__Xfit),
                           dtype=np.float64, order='C')
        return X

    def fit(self, X, y, class_weight=None, sample_weight=None):
        """Fit the SVM model according to the given training data.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target values (integers in classification, real numbers in
            regression)

        class_weight : {dict, 'auto'}, optional
            Set the parameter C of class i to class_weight[i]*C for
            SVC. If not given, all classes are supposed to have
            weight one. The 'auto' mode uses the values of y to
            automatically adjust weights inversely proportional to
            class frequencies.

        sample_weight : array-like, shape = [n_samples], optional
            Weights applied to individual samples (1. for unweighted).

        Returns
        -------
        self : object
            Returns self.

        Notes
        ------
        If X and y are not C-ordered and contiguous arrays, they are
        copied.

        """

        X = np.asarray(X, dtype=np.float64, order='C')
        y = np.asarray(y, dtype=np.float64, order='C')
        sample_weight = np.asarray([] if sample_weight is None
                                      else sample_weight, dtype=np.float64)

        if hasattr(self, 'kernel_function'):
            # you must store a reference to X to compute the kernel in predict
            # TODO: add keyword copy to copy on demand
            self.__Xfit = X
            X = self._compute_kernel(X)

        class_weight, class_weight_label = \
                     _get_class_weight(class_weight, y)

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

        params = self._get_params()
        if 'scale_C' in params:
            if params['scale_C']:
                params['C'] = params['C'] / float(X.shape[0])
            del params['scale_C']

        self.support_, self.support_vectors_, self.n_support_, \
        self.dual_coef_, self.intercept_, self.label_, self.probA_, \
        self.probB_ = libsvm.fit(X, y,
            svm_type=solver_type, sample_weight=sample_weight,
            class_weight=class_weight,
            class_weight_label=class_weight_label,
            **params)

        return self

    def predict(self, X):
        """Perform classification or regression samples in X.

        For a classification model, the predicted class for each
        sample in X is returned.  For a regression model, the function
        value of X calculated is returned.

        For an one-class model, +1 or -1 is returned.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        C : array, shape = [n_samples]
        """
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

        params = self._get_params()
        if 'scale_C' in params:
            del params['scale_C']

        svm_type = LIBSVM_IMPL.index(self.impl)
        return libsvm.predict(
            X, self.support_, self.support_vectors_, self.n_support_,
            self.dual_coef_, self.intercept_,
            self.label_, self.probA_, self.probB_,
            svm_type=svm_type, **params)

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
        X = np.asarray(X, dtype=np.float64, order='C')
        if X.ndim == 1:
            # don't use np.atleast_2d, it doesn't guarantee C-contiguity
            X = np.reshape(X, (1, -1), order='C')
        X = self._compute_kernel(X)
        if self.impl not in ('c_svc', 'nu_svc'):
            raise NotImplementedError("predict_proba only implemented for SVC "
                                      "and NuSVC")

        params = self._get_params()
        if 'scale_C' in params:
            del params['scale_C']

        svm_type = LIBSVM_IMPL.index(self.impl)
        pprob = libsvm.predict_proba(
            X, self.support_, self.support_vectors_, self.n_support_,
            self.dual_coef_, self.intercept_, self.label_,
            self.probA_, self.probB_,
            svm_type=svm_type, **params)

        return pprob

    def decision_function(self, X):
        """Distance of the samples T to the separating hyperplane.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        X : array-like, shape = [n_samples, n_class * (n_class-1) / 2]
            Returns the decision function of the sample for each class
            in the model.
        """
        X = np.asarray(X, dtype=np.float64, order='C')
        if X.ndim == 1:
            # don't use np.atleast_2d, it doesn't guarantee C-contiguity
            X = np.reshape(X, (1, -1), order='C')
        X = self._compute_kernel(X)

        params = self._get_params()
        if 'scale_C' in params:
            del params['scale_C']

        dec_func = libsvm.decision_function(
            X, self.support_, self.support_vectors_, self.n_support_,
            self.dual_coef_, self.intercept_, self.label_,
            self.probA_, self.probB_,
            svm_type=LIBSVM_IMPL.index(self.impl),
            **params)

        if self.impl != 'one_class':
            # libsvm has the convention of returning negative values for
            # rightmost labels, so we invert the sign since our label_ is
            # sorted by increasing order
            return - dec_func
        else:
            return dec_func


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
                 scale_C=False):
        self.penalty = penalty
        self.loss = loss
        self.dual = dual
        self.tol = tol
        self.C = C
        self.fit_intercept = fit_intercept
        self.intercept_scaling = intercept_scaling
        self.multi_class = multi_class
        self.scale_C = scale_C

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
        X : array-like, shape = [n_samples, n_features]
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

        self.class_weight, self.class_weight_label = \
                     _get_class_weight(class_weight, y)

        X = safe_asarray(X, dtype=np.float64, order='C')
        if not isinstance(X, np.ndarray):   # sparse X passed in by user
            raise ValueError("Training vectors should be array-like, not %s"
                             % type(X))
        y = np.asarray(y, dtype=np.int32, order='C')

        if X.shape[0] != y.shape[0]:
            raise ValueError("X and y have incompatible shapes.\n" +
                             "X has %s samples, but y has %s." % \
                             (X.shape[0], y.shape[0]))

        C = self.C
        if self.scale_C:
            C = C / float(X.shape[0])

        self.raw_coef_, self.label_ = liblinear.train_wrap(X, y,
                       self._get_solver_type(), self.tol,
                       self._get_bias(), C,
                       self.class_weight_label, self.class_weight)

        return self

    def predict(self, X):
        """Predict target values of X according to the fitted model.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        C : array, shape = [n_samples]
        """
        X = array2d(X, dtype=np.float64, order='C')
        self._check_n_features(X)

        coef = self.raw_coef_

        return liblinear.predict_wrap(X, coef,
                                      self._get_solver_type(),
                                      self.tol, self.C,
                                      self.class_weight_label,
                                      self.class_weight, self.label_,
                                      self._get_bias())

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
        X = np.asarray(X, dtype=np.float64, order='C')
        if X.ndim == 1:
            # don't use np.atleast_2d, it doesn't guarantee C-contiguity
            X = np.reshape(X, (1, -1), order='C')
        self._check_n_features(X)

        dec_func = liblinear.decision_function_wrap(
            X, self.raw_coef_, self._get_solver_type(), self.tol,
            self.C, self.class_weight_label, self.class_weight,
            self.label_, self._get_bias())

        if len(self.label_) <= 2:
            # in the two-class case, the decision sign needs be flipped
            # due to liblinear's design
            return -dec_func
        else:
            return dec_func

    def _check_n_features(self, X):
        n_features = self.raw_coef_.shape[1]
        if self.fit_intercept:
            n_features -= 1
        if X.shape[1] != n_features:
            raise ValueError("X.shape[1] should be %d, not %d." % (n_features,
                                                                   X.shape[1]))

    @property
    def intercept_(self):
        if self.fit_intercept:
            ret = self.intercept_scaling * self.raw_coef_[:, -1]
            if len(self.label_) <= 2:
                ret *= -1
            return ret
        return 0.0

    @property
    def coef_(self):
        if self.fit_intercept:
            ret = self.raw_coef_[:, : -1].copy()
        else:
            ret = self.raw_coef_.copy()

        # as coef_ is readonly property, mark the returned value as immutable
        # to avoid silencing potential bugs
        if len(self.label_) <= 2:
            ret *= -1
            ret.flags.writeable = False
            return ret
        else:
            ret.flags.writeable = False
            return ret

    def _get_bias(self):
        if self.fit_intercept:
            return self.intercept_scaling
        else:
            return -1.0


libsvm.set_verbosity_wrap(0)
