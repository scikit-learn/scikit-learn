import numpy as np

from . import libsvm, liblinear
from ..base import BaseEstimator


def _get_class_weight(class_weight, y):
    """
    Estimate class weights for unbalanced datasets.
    """
    if class_weight == 'auto':
        uy = np.unique(y)
        weight_label = np.asarray(uy, dtype=np.int32, order='C')
        weight = np.array([1.0 / np.sum(y == i) for i in uy],
                          dtype=np.float64, order='C')
        weight *= uy.shape[0] / np.sum(weight)
    else:
        weight = np.asarray(class_weight.values(),
                            dtype=np.float64, order='C')
        weight_label = np.asarray(class_weight.keys(),
                                  dtype=np.int32, order='C')

    return weight, weight_label


class BaseLibSVM(BaseEstimator):
    """
    Base class for classifiers that use libsvm as library for
    support vector machine classification and regression.

    Should not be used directly, use derived classes instead
    """

    _kernel_types = ['linear', 'poly', 'rbf', 'sigmoid', 'precomputed']
    _svm_types = ['c_svc', 'nu_svc', 'one_class', 'epsilon_svr', 'nu_svr']

    def __init__(self, impl, kernel, degree, gamma, coef0, cache_size,
                 tol, C, nu, p, shrinking, probability):

        if not impl in self._svm_types:
            raise ValueError("impl should be one of %s, %s was given" % (
                self._svm_types, impl))

        if not (kernel in self._kernel_types or hasattr(kernel, '__call__')):
            raise ValueError("kernel should be one of %s or a callable, " \
                             "%s was given." % (self._kernel_types, kernel))

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

    def _get_kernel(self, X):
        """ Get the kernel type code as well as the data transformed by
            the kernel (if the kernel is a callable.
        """
        if hasattr(self.kernel, '__call__'):
            # in the case of precomputed kernel given as a function, we
            # have to compute explicitly the kernel matrix
            _X = np.asanyarray(self.kernel(X, self.__Xfit),
                               dtype=np.float64, order='C')
            kernel_type = 4
        else:
            kernel_type = self._kernel_types.index(self.kernel)
            _X = X
        return kernel_type, _X

    def fit(self, X, y, class_weight={}, sample_weight=[], **params):
        """
        Fit the SVM model according to the given training data and
        parameters.

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
        """
        self._set_params(**params)

        X = np.asanyarray(X, dtype=np.float64, order='C')
        y = np.asanyarray(y, dtype=np.float64, order='C')
        sample_weight = np.asanyarray(sample_weight, dtype=np.float64,
                                      order='C')

        if hasattr(self.kernel, '__call__'):
            # you must store a reference to X to compute the kernel in predict
            # there's a way around this, but it involves patching libsvm
            # TODO: put keyword copy to copy on demand
            self.__Xfit = X
        kernel_type, _X = self._get_kernel(X)

        self.class_weight, self.class_weight_label = \
                     _get_class_weight(class_weight, y)

        # check dimensions
        solver_type = self._svm_types.index(self.impl)
        if solver_type != 2 and _X.shape[0] != y.shape[0]:
            raise ValueError("X and y have incompatible shapes.\n" +
                             "X has %s features, but y has %s." % \
                             (_X.shape[0], y.shape[0]))

        if self.kernel == "precomputed" and X.shape[0] != X.shape[1]:
            raise ValueError("X.shape[0] should be equal to X.shape[1]")

        if (kernel_type in [1, 2]) and (self.gamma == 0):
            # if custom gamma is not provided ...
            self.gamma = 1.0 / _X.shape[0]

        self.shape_fit_ = X.shape

        self.support_, self.support_vectors_, self.n_support_, \
        self.dual_coef_, self.intercept_, self.label_, self.probA_, \
        self.probB_ = \
        libsvm.train(_X, y, solver_type, kernel_type, self.degree,
                      self.gamma, self.coef0, self.tol, self.C,
                      self.nu, self.cache_size, self.p,
                      self.class_weight_label, self.class_weight,
                      sample_weight, int(self.shrinking),
                      int(self.probability))

        return self

    def predict(self, X):
        """
        This function does classification or regression on an array of
        test vectors X.

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
        X = np.atleast_2d(np.asanyarray(X, dtype=np.float64, order='C'))
        n_samples, n_features = X.shape
        kernel_type, X = self._get_kernel(X)

        if self.kernel == "precomputed":
            if X.shape[1] != self.shape_fit_[0]:
                raise ValueError("X.shape[1] should be equal to the number of "
                                 "samples at training time!")
        elif n_features != self.shape_fit_[1]:
            raise ValueError("X.shape[1] should be equal to the number of "
                             "features at training time!")

        return libsvm.predict( X, self.support_vectors_,
            self.dual_coef_, self.intercept_,
            self._svm_types.index(self.impl), kernel_type,
            self.degree, self.gamma, self.coef0, self.tol, self.C,
            self.nu, self.cache_size, self.p, self.n_support_,
            self.support_, self.label_, self.class_weight_label,
            self.class_weight, self.probA_, self.probB_,
            int(self.shrinking), int(self.probability))

    def predict_proba(self, T):
        """
        This function does classification or regression on a test vector T
        given a model with probability information.

        Parameters
        ----------
        T : array-like, shape = [n_samples, n_features]

        Returns
        -------
        T : array-like, shape = [n_samples, n_classes]
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
        T = np.atleast_2d(np.asanyarray(T, dtype=np.float64, order='C'))
        kernel_type, T = self._get_kernel(T)
        if self.impl not in ('c_svc', 'nu_svc'):
            raise NotImplementedError

        pprob = libsvm.predict_proba(T, self.support_vectors_,
                      self.dual_coef_, self.intercept_,
                      self._svm_types.index(self.impl), kernel_type,
                      self.degree, self.gamma, self.coef0, self.tol,
                      self.C, self.nu, self.cache_size,
                      self.p, self.n_support_,
                      self.support_, self.label_,
                      self.class_weight_label,
                      self.class_weight, 
                      self.probA_, self.probB_, int(self.shrinking),
                      int(self.probability))

        return pprob

    def predict_log_proba(self, T):
        """
        This function does classification or regression on a test vector T
        given a model with probability information.

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

    def decision_function(self, T):
        """
        Calculate the distance of the samples T to the separating hyperplane.

        Parameters
        ----------
        T : array-like, shape = [n_samples, n_features]

        Returns
        -------
        T : array-like, shape = [n_samples, n_class * (n_class-1) / 2]
            Returns the decision function of the sample for each class
            in the model.
        """
        T = np.atleast_2d(np.asanyarray(T, dtype=np.float64, order='C'))
        kernel_type, T = self._get_kernel(T)

        dec_func = libsvm.decision_function(T, self.support_vectors_,
                      self.dual_coef_, self.intercept_,
                      self._svm_types.index(self.impl), kernel_type,
                      self.degree, self.gamma, self.coef0, self.tol,
                      self.C, self.class_weight_label,
                      self.class_weight, self.nu, self.cache_size,
                      self.p, int(self.shrinking),
                      int(self.probability), self.n_support_,
                      self.support_, self.label_, self.probA_,
                      self.probB_)

        if self.impl != 'one_class':
            # libsvm has the convention of returning negative values for
            # rightmost labels, so we invert the sign since our label_ is
            # sorted by increasing order
            return -dec_func
        else:
            return dec_func

    @property
    def coef_(self):
        if self.kernel != 'linear':
            raise NotImplementedError('coef_ is only available when using a linear kernel')
        return np.dot(self.dual_coef_, self.support_vectors_)


class BaseLibLinear(BaseEstimator):
    """
    Base for classes binding liblinear (dense and sparse versions)
    """

    _solver_type_dict = {
        'PL2_LLR_D0' : 0,  # L2 penalty, logistic regression
        'PL2_LL2_D1' : 1,  # L2 penalty, L2 loss, dual form
        'PL2_LL2_D0' : 2,  # L2 penalty, L2 loss, primal form
        'PL2_LL1_D1' : 3,  # L2 penalty, L1 Loss, dual form
        'MC_SVC'     : 4,  # Multi-class Support Vector Classification
        'PL1_LL2_D0' : 5,  # L1 penalty, L2 Loss, primal form
        'PL1_LLR_D0' : 6,  # L1 penalty, logistic regression
        'PL2_LLR_D1' : 7,  # L2 penalty, logistic regression, dual form
        }

    def __init__(self, penalty='l2', loss='l2', dual=True, tol=1e-4, C=1.0,
                 multi_class=False, fit_intercept=True, intercept_scaling=1):
        self.penalty = penalty
        self.loss = loss
        self.dual = dual
        self.tol = tol
        self.C = C
        self.fit_intercept = fit_intercept
        self.intercept_scaling = intercept_scaling
        self.multi_class = multi_class

        # Check that the arguments given are valid:
        self._get_solver_type()

    def _get_solver_type(self):
        """ Return the magic number for the solver described by the
            settings.
        """
        if self.multi_class:
            solver_type = 'MC_SVC'
        else:
            solver_type = "P%s_L%s_D%d" % (
                self.penalty.upper(), self.loss.upper(), int(self.dual))
        if not solver_type in self._solver_type_dict:
            raise ValueError('Not supported set of arguments: '
                             + solver_type)
        return self._solver_type_dict[solver_type]

    def fit(self, X, y, class_weight={}, **params):
        """
        Fit the model according to the given training data and
        parameters.

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
        self._set_params(**params)

        self.class_weight, self.class_weight_label = \
                     _get_class_weight(class_weight, y)

        X = np.asanyarray(X, dtype=np.float64, order='C')
        y = np.asanyarray(y, dtype=np.int32, order='C')

        self.raw_coef_, self.label_ = liblinear.train_wrap(X, y,
                       self._get_solver_type(), self.tol,
                       self._get_bias(), self.C,
                       self.class_weight_label, self.class_weight)

        return self

    def predict(self, X):
        """
        Predict target values of X according to the fitted model.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        C : array, shape = [n_samples]
        """
        X = np.asanyarray(X, dtype=np.float64, order='C')
        self._check_n_features(X)

        coef = self.raw_coef_

        return liblinear.predict_wrap(X, coef,
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
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        T : array-like, shape = [n_samples, n_class]
            Returns the decision function of the sample for each class
            in the model.
        """
        X = np.atleast_2d(np.asanyarray(X, dtype=np.float64, order='C'))
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
            ret = self.raw_coef_[:, : -1]
        else:
            ret = self.raw_coef_
        if len(self.label_) <= 2:
            return -ret
        else:
            return ret

    def predict_proba(self, T):
        # only available for logistic regression
        raise NotImplementedError(
                'liblinear does not provide this functionality')

    def _get_bias(self):
        if self.fit_intercept:
            return self.intercept_scaling
        else:
            return -1.0


libsvm.set_verbosity_wrap(0)
