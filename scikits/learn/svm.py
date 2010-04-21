import numpy as np
from . import libsvm, liblinear

_kernel_types = ['linear', 'poly', 'rbf', 'sigmoid', 'precomputed']
_svm_types = ['c_svc', 'nu_svc', 'one_class', 'epsilon_svr', 'nu_svr']


class BaseLibsvm(object):
    """
    Base class for classifiers that use support vector machine.

    Should not be used directly, use derived classes instead

    Parameters
    ----------
    X : array-like, shape = [N, D]
        It will be converted to a floating-point array.
    y : array, shape = [N]
        target vector relative to X
        It will be converted to a floating-point array.
    """
    support_ = None
    coef_ = None
    rho_ = None

    weight = np.empty(0, dtype=np.float64)
    weight_label = np.empty(0, dtype=np.int32)

    def __init__(self, impl, kernel, degree, gamma, coef0, cache_size,
                 eps, C, nr_weight, nu, p, shrinking, probability):
        self.svm = _svm_types.index(impl)
        self.kernel = _kernel_types.index(kernel)
        self.degree = degree
        self.gamma = gamma
        self.coef0 = coef0
        self.cache_size = cache_size
        self.eps = eps
        self.C = C
        self.nr_weight = 0
        self.nu = nu
        self.p = p
        self.shrinking = int(shrinking)
        self.probability = int(probability)

    def fit(self, X, y):
        """
        Fit the model with vectors X, Y.

        Parameters
        ----------
        X : array-like, shape = [nsamples, nfeatures]
            Training vector, where nsamples in the number of samples and
            nfeatures is the number of features.
        Y : array, shape = [nsamples]
            Target vector relative to X

        """
        X = np.asanyarray(X, dtype=np.float, order='C')
        y = np.asanyarray(y, dtype=np.float, order='C')

        # check dimensions
        if X.shape[0] != y.shape[0]: raise ValueError("Incompatible shapes")

        if (self.gamma == 0): self.gamma = 1.0/X.shape[0]
        self.coef_, self.rho_, self.support_, self.nclass_, self.nSV_, \
                 self.label_, self.probA_, self.probB_ = libsvm.train_wrap(X, y,
                 self.svm, self.kernel, self.degree, self.gamma,
                 self.coef0, self.eps, self.C, self.nr_weight,
                 self.weight_label, self.weight, self.nu, self.cache_size, self.p,
                 self.shrinking, int(self.probability))
        return self

    def predict(self, T):
        T = np.asanyarray(T, dtype=np.float, order='C')
        return libsvm.predict_from_model_wrap(T, self.support_,
                      self.coef_, self.rho_, self.svm,
                      self.kernel, self.degree, self.gamma,
                      self.coef0, self.eps, self.C, self.nr_weight,
                      np.empty(0, dtype=np.int), np.empty(0,
                      dtype=np.float), self.nu, self.cache_size,
                      self.p, self.shrinking, self.probability,
                      self.nclass_, self.nSV_, self.label_,
                      self.probA_, self.probB_)

    def predict_proba(self, T):
        if not self.probability:
            raise ValueError("probability estimates must be enabled to use this method")
        T = np.asanyarray(T, dtype=np.float, order='C')
        return libsvm.predict_prob_from_model_wrap(T, self.support_,
                      self.coef_, self.rho_, self.svm,
                      self.kernel, self.degree, self.gamma,
                      self.coef0, self.eps, self.C, self.nr_weight,
                      np.empty(0, dtype=np.int), np.empty(0,
                      dtype=np.float), self.nu, self.cache_size,
                      self.p, self.shrinking, self.probability,
                      self.nclass_, self.nSV_, self.label_,
                      self.probA_, self.probB_)

###
# Public API
# No processing should go into these classes

class SVC(BaseLibsvm):
    """
    Support Vector Classification

    Implements C-SVC, Nu-SVC

    Parameters
    ----------
    X : array-like, shape = [nsamples, nfeatures]
        Training vector, where nsamples in the number of samples and
        nfeatures is the number of features.
    Y : array, shape = [nsamples]
        Target vector relative to X

    impl : string, optional
        SVM implementation to choose from. This refers to different
        formulations of the SVM optimization problem.
        Can be one of 'c_svc', 'nu_svc'. By default 'c_svc' will be chosen.

    nu : float, optional
        An upper bound on the fraction of training errors and a lower
        bound of the fraction of support vectors. Should be in the
        interval (0, 1].  By default 0.5 will be taken.  Only
        available if impl='nu_svc'

    kernel : string, optional
         Specifies the kernel type to be used in the algorithm.
         one of 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed'.
         If none is given 'rbf' will be used.

    degree : int, optional
        degree of kernel function
        is significant only in poly, rbf, sigmoid

    probability: boolean, optional (False by default)
        especify if probability estimates must be enabled
        must be enabled prior to calling prob_predict

    coef0 : float, optional

    Attributes
    ----------
    `support_` : array-like, shape = [nSV, nfeatures]
        Support vectors

    `coef_` : array, shape = [nclasses-1, nfeatures]
        Coefficient of the support vector in the decision function.

    `rho_` : array, shape = [nclasses-1]
        constants in decision function

    Methods
    -------
    fit(X, Y) : self
        Fit the model

    predict(X) : array
        Predict using the model.

    Examples
    --------
    >>> X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    >>> Y = np.array([1, 1, 2, 2])
    >>> clf = SVC()
    >>> clf.fit(X, Y)    #doctest: +ELLIPSIS
    <scikits.learn.svm.SVC object at 0x...>
    >>> print clf.predict([[-0.8, -1]])
    [ 1.]

    See also
    --------
    SVR
    """

    def __init__(self, impl='c_svc', kernel='rbf', degree=3,
                 gamma=0.0, coef0=0.0, cache_size=100.0, eps=1e-3,
                 C=1.0, nr_weight=0, nu=0.5, p=0.1, shrinking=True,
                 probability=False):
        BaseLibsvm.__init__(self, impl, kernel, degree, gamma, coef0,
                         cache_size, eps, C, nr_weight, nu, p,
                         shrinking, probability)


class SVR(BaseLibsvm):
    """
    Support Vector Regression.

    Attributes
    ----------
    `support_` : array-like, shape = [nSV, nfeatures]
        Support vectors

    `coef_` : array, shape = [nclasses-1, nfeatures]
        Coefficient of the support vector in the decision function.

    `rho_` : array, shape = [nclasses-1]
        constants in decision function

    Methods
    -------
    fit(X, Y) : self
        Fit the model

    predict(X) : array
        Predict using the model.

    See also
    --------
    SVC
    """
    def __init__(self, impl='epsilon_svr', kernel='rbf', degree=3,
                 gamma=0.0, coef0=0.0, cache_size=100.0, eps=1e-3,
                 C=1.0, nr_weight=0, nu=0.5, p=0.1, shrinking=True,
                 probability=False):
        BaseLibsvm.__init__(self, impl, kernel, degree, gamma, coef0,
                         cache_size, eps, C, nr_weight, nu, p,
                         shrinking, probability)

class OneClassSVM(BaseLibsvm):
    """
    Outlayer detection

    Methods
    -------
    fit(X, Y) : self
        Fit the model

    predict(X) : array
        Predict using the model.
    """
    def __init__(self, kernel='rbf', degree=3,
                 gamma=0.0, coef0=0.0, cache_size=100.0, eps=1e-3,
                 C=1.0, nr_weight=0, nu=0.5, p=0.1, shrinking=True,
                 probability=False):
        impl = 'one_class'
        BaseLibsvm.__init__(self, impl, kernel, degree, gamma, coef0,
                         cache_size, eps, C, nr_weight, nu, p,
                         shrinking, probability)


class LinearSVC(object):
    """
    Linear Support Vector Classification.


    Parameters
    ----------
    
    Also accepts parameter penalty, that can have values 'l1' or 'l2'

    Similar to SVC with parameter kernel='linear', but uses internally
    liblinear rather than libsvm, so it has more flexibility in the
    choice of penalties and loss functions and should be faster for
    huge datasets.

    TODO: wrap Cramer & Singer
    """
    _solver_type = {'l2l2_1': 1, 'l2l2_0' : 2, 'l2l1_1' : 3, 'l1l2_0' : 5}

    def __init__(self, penalty='l2', loss='l2', dual=False, eps=1e-4, C=1.0):
        s = penalty + loss + '_' + str(int(dual))
        try: self.solver_type = self._solver_type[s]
        except KeyError:
            raise ValueError('Not supported set of arguments')
        self.eps = eps
        self.C = C

    _penalties = {'l2': 0, 'l1' : 6}
    _weight_label = np.empty(0, dtype=np.int)
    _weight = np.empty(0, dtype=np.float64)

    def fit(self, X, Y):
        X = np.asanyarray(X, dtype=np.float64, order='C')
        Y = np.asanyarray(Y, dtype=np.int, order='C')
        self.coef_, self.label_, self.bias_ = liblinear.train_wrap(X,
                                          Y, self.solver_type, self.eps, 1.0,
                                          self.C, 0,
                                          self._weight_label,
                                          self._weight)

    def predict(self, T):
        T = np.asanyarray(T, dtype=np.float64, order='C')
        return liblinear.predict_wrap(T, self.coef_, self.solver_type,
                                      self.eps, self.C,
                                      self._weight_label,
                                      self._weight, self.label_,
                                      self.bias_)


    def predict_proba(self, T):
        raise NotImplementedError('liblinear does not provide this functionality')

