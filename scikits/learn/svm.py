import numpy as np
import libsvm

_kernel_types = ['linear', 'poly', 'rbf', 'sigmoid', 'precomputed']
_svm_types = ['c_svc', 'nu_svc', 'one_class', 'epsilon_svr', 'nu_svr']


class BaseSVM(object):
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

    def __init__(self, svm, kernel, degree, gamma, coef0, cache_size,
                 eps, C, nr_weight, nu, p, shrinking, probability,
                 scale):
        self.svm = _svm_types.index(svm)
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
        self.shrinking = shrinking
        self.probability = probability
        self.scale = scale

    def fit(self, X, y):
        """
        should empty arrays created be order='C' ?
        """
        X = np.asanyarray(X, dtype=np.float, order='C')
        y = np.asanyarray(y, dtype=np.float, order='C')

        if self.scale:
            self.mean = X.mean(0)
            self.std = X.std(0)
            X = (X - self.mean) / self.std

        # check dimensions
        if X.shape[0] != y.shape[0]: raise ValueError("Incompatible shapes")

        if (self.gamma == 0): self.gamma = 1.0/X.shape[0]
        self.coef_, self.rho_, self.support_, self.nclass_, self.nSV_, self.label_  = \
             libsvm.train_wrap(X, y, self.svm, self.kernel, self.degree,
                 self.gamma, self.coef0, self.eps, self.C, self.nr_weight,
                 np.empty(0, dtype=np.int), np.empty(0, dtype=np.float), self.nu,
                 self.cache_size, self.p, self.shrinking, self.probability)
        return self

    def predict(self, T):
        T = np.asanyarray(T, dtype=np.float, order='C')
        if self.scale: T = (T - self.mean) / self.std
        return libsvm.predict_from_model_wrap(T, self.support_,
                      self.coef_, self.rho_, self.svm,
                      self.kernel, self.degree, self.gamma,
                      self.coef0, self.eps, self.C, self.nr_weight,
                      np.empty(0, dtype=np.int), np.empty(0,
                      dtype=np.float), self.nu, self.cache_size,
                      self.p, self.shrinking, self.probability,
                      self.nclass_, self.nSV_, self.label_)


def predict(X, y, T, svm='c_svc', kernel='rbf', degree=3,
                 gamma=0.0, coef0=0.0, cache_size=100.0, eps=1e-3,
                 C=1.0, nr_weight=0, nu=0.5, p=0.1, shrinking=1,
                 probability=0):
    """
    Shortcut that does fit and predict in a single step.

    Should be faster than instatating the object, since less copying is done.

    Parameters
    ----------
    X : array-like
        data points
    y : array
        targets
    T : array
        test points

    Optional Parameters
    -------------------
    TODO

    Examples
    --------
    """
    X = np.asanyarray(X, dtype=np.float, order='C')
    y = np.asanyarray(y, dtype=np.float, order='C')
    T = np.asanyarray(T, dtype=np.float, order='C')
    if X.shape[0] != y.shape[0]: raise ValueError("Incompatible shapes")
    return libsvm.predict_wrap(X, y, T, _svm_types.index(svm),
                               _kernel_types.index(kernel),
                               degree, gamma, coef0, eps, C,
                               nr_weight, np.empty(0, dtype=np.int),
                               np.empty(0, dtype=np.float), nu,
                               cache_size, p, shrinking, probability)


###
# Public API
# No processing should go into these classes

class SVC(BaseSVM):
    """
    Support Vector Classification

    Implementats C-SVC, nu-SVC

    Parameters
    ----------
    X : array-like, shape = [nsamples, nfeatures]
        Training vector, where nsamples in the number of samples and
        nfeatures is the number of features.
    Y : array, shape = [nsamples]
        Target vector relative to X

    Optional Parameters
    -------------------
    impl : string
        SVM implementation to choose from. This refers to different
        formulations of the SVM optimization problem.
        Can be one of 'c_svc', 'nu_svc'. By default 'c_svc' will be chosen.

    scale : boolean
        Scale data before fitting the model. Caution: if activated,
        estimated parameters will also be scaled.
        Default is True.

    nu: float
        An upper bound on the fraction of training errors and a lower
        bound of the fraction of support vectors. Should be in the
        interval (0, 1].
        By default 0.5 will be taken.
        Only available is impl is set to 'nu_svc'


    kernel : string
         Specifies the kernel type to be used in the algorithm.
         one of 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed'.
         If none is given 'rbf' will be used.

    degree : int
        degree of kernel function
        is significant only in POLY, RBF, SIGMOID

    Members
    -------
    support_ : array-like, shape = [nSV, D]
        estimated support vectors.
        where nSV is the number of support vectors, D is the dimension
        of the underlying space.

    coef_ : array
        coefficient of the support vector in the decission function.

    rho_ : array
        constants in decision function


    Examples
    --------
    >>> X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    >>> y = np.array([1, 1, 2, 2])
    >>> clf = SVM()
    >>> clf.fit(X, y)    #doctest: +ELLIPSIS
    <scikits.learn.svm.svm.SVM object at 0x...>
    >>> print clf.predict([[-0.8, -1]])
    [ 1.]

    References
    ----------
    http://scikit-learn.sourceforge.net/doc/modules/svm.html
    http://www.csie.ntu.edu.tw/~cjlin/papers/libsvm.pdf
    """
    def __init__(self, impl='c_svc', kernel='rbf', degree=3,
                 gamma=0.0, coef0=0.0, cache_size=100.0, eps=1e-3,
                 C=1.0, nr_weight=0, nu=0.5, p=0.1, shrinking=1,
                 probability=0, scale=True):

        BaseSVM.__init__(self, impl, kernel, degree, gamma, coef0,
                         cache_size, eps, C, nr_weight, nu, p,
                         shrinking, probability, scale)    


class SVR(BaseSVM):
    """
    Support Vector Regression.

    Parameters
    ----------
    X : array-like, shape = [N, D]
        Training vector
    Y : array, shape = [N]
        Target vector relative to X
    
    """
    def __init__(self, svm='epsilon_svr', kernel='rbf', degree=3,
                 gamma=0.0, coef0=0.0, cache_size=100.0, eps=1e-3,
                 C=1.0, nr_weight=0, nu=0.5, p=0.1, shrinking=1,
                 probability=0, scale=True):
        BaseSVM.__init__(self, svm, kernel, degree, gamma, coef0,
                         cache_size, eps, C, nr_weight, nu, p,
                         shrinking, probability, scale)

class OneClassSVM(BaseSVM):
    """
    Outlayer detection
    """
    def __init__(self, kernel='rbf', degree=3,
                 gamma=0.0, coef0=0.0, cache_size=100.0, eps=1e-3,
                 C=1.0, nr_weight=0, nu=0.5, p=0.1, shrinking=1,
                 probability=0, scale=True):
        svm = 'one_class'
        BaseSVM.__init__(self, svm, kernel, degree, gamma, coef0,
                         cache_size, eps, C, nr_weight, nu, p,
                         shrinking, probability, scale)
