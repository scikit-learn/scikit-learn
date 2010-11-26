
from ..base import ClassifierMixin, RegressorMixin
from .base import BaseLibSVM


class SVC(BaseLibSVM, ClassifierMixin):
    """
    C-Support Vector Classification.


    Parameters
    ----------
    C : float, optional (default=1.0)
        penalty parameter C of the error term.

    kernel : string, optional
         Specifies the kernel type to be used in the algorithm.
         one of 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed'.
         If none is given 'rbf' will be used.

    degree : int, optional
        degree of kernel function
        is significant only in poly, rbf, sigmoid

    gamma : float, optional
        kernel coefficient for rbf and poly, by default 1/n_features
        will be taken.

    coef0 : float, optional
        independent term in kernel function. It is only significant
        in poly/sigmoid.

    probability: boolean, optional (False by default)
        enable probability estimates. This must be enabled prior
        to calling prob_predict.

    shrinking: boolean, optional
         wether to use the shrinking heuristic.

    eps: float, optional
         precision for stopping criteria

    cache_size: float, optional
         specify the size of the cache (in MB)


    Attributes
    ----------
    `support_` : array-like, shape = [n_SV]
        Index of support vectors.

    `support_vectors_` : array-like, shape = [n_SV, n_features]
        Support vectors.

    `n_support_` : array-like, dtype=int32, shape = [n_class]
        number of support vector for each class.

    `dual_coef_` : array, shape = [n_class-1, n_SV]
        Coefficients of the support vector in the decision function.

    `coef_` : array, shape = [n_class-1, n_features]
        Weights asigned to the features (coefficients in the primal
        problem). This is only available in the case of linear kernel.

    `intercept_` : array, shape = [n_class * (n_class-1) / 2]
        Constants in decision function.


    Examples
    --------
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    >>> y = np.array([1, 1, 2, 2])
    >>> from scikits.learn.svm import SVC
    >>> clf = SVC()
    >>> clf.fit(X, y)
    SVC(kernel='rbf', C=1.0, probability=False, degree=3, coef0=0.0, eps=0.001,
      cache_size=100.0, shrinking=True, gamma=0.25)
    >>> print clf.predict([[-0.8, -1]])
    [ 1.]

    See also
    --------
    SVR, LinearSVC
    """

    def __init__(self, C=1.0, kernel='rbf', degree=3, gamma=0.0,
                 coef0=0.0, shrinking=True, probability=False,
                 eps=1e-3, cache_size=100.0):

        BaseLibSVM.__init__(self, 'c_svc', kernel, degree, gamma, coef0,
                         cache_size, eps, C, 0., 0.,
                         shrinking, probability)


class NuSVC(BaseLibSVM, ClassifierMixin):
    """
    Nu-Support Vector Classification.


    Parameters
    ----------
    nu : float, optional
        An upper bound on the fraction of training errors and a lower
        bound of the fraction of support vectors. Should be in the
        interval (0, 1].  By default 0.5 will be taken.

    kernel : string, optional
         Specifies the kernel type to be used in the algorithm.
         one of 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed'.
         If none is given 'rbf' will be used.

    degree : int, optional
        degree of kernel function
        is significant only in poly, rbf, sigmoid

    gamma : float, optional
        kernel coefficient for rbf and poly, by default 1/n_features
        will be taken.

    probability: boolean, optional (False by default)
        enable probability estimates. This must be enabled prior
        to calling prob_predict.

    coef0 : float, optional
        independent term in kernel function. It is only significant
        in poly/sigmoid.

    shrinking: boolean, optional
         wether to use the shrinking heuristic.

    eps: float, optional
         precision for stopping criteria

    cache_size: float, optional
         specify the size of the cache (in MB)


    Attributes
    ----------
    `support_` : array-like, shape = [n_SV]
        Index of support vectors.

    `support_vectors_` : array-like, shape = [n_SV, n_features]
        Support vectors.

    `n_support_` : array-like, dtype=int32, shape = [n_class]
        number of support vector for each class.

    `dual_coef_` : array, shape = [n_classes-1, n_SV]
        Coefficients of the support vector in the decision function.

    `coef_` : array, shape = [n_classes-1, n_features]
        Weights asigned to the features (coefficients in the primal
        problem). This is only available in the case of linear kernel.

    `intercept_` : array, shape = [n_class * (n_class-1) / 2]
        Constants in decision function.


    Methods
    -------
    fit(X, y) : self
        Fit the model

    predict(X) : array
        Predict using the model.

    predict_proba(X) : array
        Return probability estimates.

    decision_function(X) : array
        Return distance to predicted margin.

    Examples
    --------
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    >>> y = np.array([1, 1, 2, 2])
    >>> from scikits.learn.svm import NuSVC
    >>> clf = NuSVC()
    >>> clf.fit(X, y)
    NuSVC(kernel='rbf', probability=False, degree=3, coef0=0.0, eps=0.001,
       cache_size=100.0, shrinking=True, nu=0.5, gamma=0.25)
    >>> print clf.predict([[-0.8, -1]])
    [ 1.]

    See also
    --------
    SVC, LinearSVC, SVR
    """

    def __init__(self, nu=0.5, kernel='rbf', degree=3, gamma=0.0,
                 coef0=0.0, shrinking=True, probability=False,
                 eps=1e-3, cache_size=100.0):

        BaseLibSVM.__init__(self, 'nu_svc', kernel, degree, gamma, coef0,
                         cache_size, eps, 0., nu, 0.,
                         shrinking, probability)


class SVR(BaseLibSVM, RegressorMixin):
    """
    Support Vector Regression.


    Parameters
    ----------
    nu : float, optional
        An upper bound on the fraction of training errors and a lower bound of
        the fraction of support vectors. Should be in the interval (0, 1].  By
        default 0.5 will be taken.  Only available if impl='nu_svc'

    kernel : string, optional
         Specifies the kernel type to be used in the algorithm.
         one of 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed'.
         If none is given 'rbf' will be used.

    p : float
        epsilon in the epsilon-SVR model.

    degree : int, optional
        degree of kernel function
        is significant only in poly, rbf, sigmoid

    gamma : float, optional
        kernel coefficient for rbf and poly, by default 1/n_features
        will be taken.

    C : float, optional (default=1.0)
        penalty parameter C of the error term.

    probability: boolean, optional (False by default)
        enable probability estimates. This must be enabled prior
        to calling prob_predict.

    eps: float, optional
         precision for stopping criteria

    coef0 : float, optional
        independent term in kernel function. It is only significant
        in poly/sigmoid.

    cache_size: float, optional
         specify the size of the cache (in MB)

    shrinking: boolean, optional
         wether to use the shrinking heuristic.

    Attributes
    ----------
    `support_` : array-like, shape = [n_SV]
        Index of support vectors.

    `support_vectors_` : array-like, shape = [nSV, n_features]
        Support vectors.

    `dual_coef_` : array, shape = [n_classes-1, n_SV]
        Coefficients of the support vector in the decision function.

    `coef_` : array, shape = [n_classes-1, n_features]
        Weights asigned to the features (coefficients in the primal
        problem). This is only available in the case of linear kernel.

    `intercept_` : array, shape = [n_class * (n_class-1) / 2]
        Constants in decision function.

    See also
    --------
    NuSVR
    """
    def __init__(self, kernel='rbf', degree=3, gamma=0.0, coef0=0.0,
                 cache_size=100.0, eps=1e-3, C=1.0, nu=0.5, p=0.1,
                 shrinking=True, probability=False):

        BaseLibSVM.__init__(self, 'epsilon_svr', kernel, degree, gamma, coef0,
                         cache_size, eps, C, nu, p,
                         shrinking, probability)

    def fit(self, X, y, sample_weight=[]):
        """
        Fit the SVM model according to the given training data and parameters.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vector, where n_samples is the number of samples and
            n_features is the number of features.
        y : array, shape = [n_samples]
            Target values. Array of floating-point numbers.

        Returns
        -------
        self : object
            Returns self.
        """
        # we copy this method because SVR does not accept class_weight
        return BaseLibSVM.fit(self, X, y, sample_weight=sample_weight)


class NuSVR(BaseLibSVM, RegressorMixin):
    """
    Nu Support Vector Regression. Similar to NuSVC, for regression,
    uses a paramter nu to control the number of support
    vectors. However, unlike NuSVC, where nu replaces with C, here nu
    replaces with the parameter p of SVR.

    Parameters
    ----------
    nu : float, optional
        An upper bound on the fraction of training errors and a lower bound of
        the fraction of support vectors. Should be in the interval (0, 1].  By
        default 0.5 will be taken.  Only available if impl='nu_svc'

    C : float, optional (default=1.0)
        penalty parameter C of the error term.

    kernel : string, optional
         Specifies the kernel type to be used in the algorithm.
         one of 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed'.
         If none is given 'rbf' will be used.

    degree : int, optional
        degree of kernel function
        is significant only in poly, rbf, sigmoid

    gamma : float, optional
        kernel coefficient for rbf and poly, by default 1/n_features
        will be taken.

    eps: float, optional
         precision for stopping criteria

    probability: boolean, optional (False by default)
        enable probability estimates. This must be enabled prior
        to calling prob_predict.

    coef0 : float, optional
        independent term in kernel function. It is only significant
        in poly/sigmoid.

    shrinking: boolean, optional
         wether to use the shrinking heuristic.

    cache_size: float, optional
         specify the size of the cache (in MB)

    Attributes
    ----------
    `support_` : array-like, shape = [n_SV]
        Index of support vectors.

    `support_vectors_` : array-like, shape = [nSV, n_features]
        Support vectors.

    `dual_coef_` : array, shape = [n_classes-1, n_SV]
        Coefficients of the support vector in the decision function.

    `coef_` : array, shape = [n_classes-1, n_features]
        Weights asigned to the features (coefficients in the primal
        problem). This is only available in the case of linear kernel.

    `intercept_` : array, shape = [n_class * (n_class-1) / 2]
        Constants in decision function.

    See also
    --------
    NuSVR
    """

    def __init__(self, nu=0.5, C=1.0, kernel='rbf', degree=3,
                 gamma=0.0, coef0=0.0, shrinking=True,
                 probability=False, cache_size=100.0, eps=1e-3):

        BaseLibSVM.__init__(self, 'epsilon_svr', kernel, degree, gamma, coef0,
                         cache_size, eps, C, nu, 0.,
                         shrinking, probability)

    def fit(self, X, y):
        """
        Fit the SVM model according to the given training data and parameters.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vector, where n_samples is the number of samples and
            n_features is the number of features.
        y : array, shape = [n_samples]
            Target values. Array of floating-point numbers.

        Returns
        -------
        self : object
            Returns self.
        """
        # we copy this method because SVR does not accept class_weight
        return BaseLibSVM.fit(self, X, y)


class OneClassSVM(BaseLibSVM):
    """Unsupervised outliers detection

    Estimate the support of a high-dimensional distribution.

    Parameters
    ----------
    kernel : string, optional
        Specifies the kernel type to be used in
        the algorithm. Can be one of 'linear', 'poly', 'rbf', 'sigmoid',
        'precomputed'. If none is given 'rbf' will be used.

    nu : float, optional
        An upper bound on the fraction of training
        errors and a lower bound of the fraction of support
        vectors. Should be in the interval (0, 1]. By default 0.5
        will be taken.

    degree : int, optional
        Degree of kernel function. Significant only in poly, rbf, sigmoid.

    gamma : float, optional
        kernel coefficient for rbf and poly, by default 1/n_features
        will be taken.

    coef0 : float, optional
        Independent term in kernel function. It is only significant in
        poly/sigmoid.

    eps: float, optional
         precision for stopping criteria

    shrinking: boolean, optional
         wether to use the shrinking heuristic.

    cache_size: float, optional
         specify the size of the cache (in MB)

    Attributes
    ----------
    `support_` : array-like, shape = [n_SV]
        Index of support vectors.

    `support_vectors_` : array-like, shape = [nSV, n_features]
        Support vectors.

    `dual_coef_` : array, shape = [n_classes-1, n_SV]
        Coefficient of the support vector in the decision function.

    `coef_` : array, shape = [n_classes-1, n_features]
        Weights asigned to the features (coefficients in the primal
        problem). This is only available in the case of linear kernel.

    `intercept_` : array, shape = [n_classes-1]
        Constants in decision function.

    """
    def __init__(self, kernel='rbf', degree=3, gamma=0.0, coef0=0.0,
                 cache_size=100.0, eps=1e-3, nu=0.5, shrinking=True):
        BaseLibSVM.__init__(self, 'one_class', kernel, degree, gamma, coef0,
                             cache_size, eps, 0.0, nu, 0.0, shrinking, False)

    def fit(self, X):
        """
        Detects the soft boundary (aka soft boundary) of the set of samples X.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Set of samples, where n_samples is the number of samples and
            n_features is the number of features.

        """
        super(OneClassSVM, self).fit(X, [])
