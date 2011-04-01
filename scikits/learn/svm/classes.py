from ..base import ClassifierMixin, RegressorMixin
from ..linear_model.base import CoefSelectTransformerMixin
from .base import BaseLibLinear, BaseLibSVM


class LinearSVC(BaseLibLinear, ClassifierMixin, CoefSelectTransformerMixin):
    """Linear Support Vector Classification.

    Similar to SVC with parameter kernel='linear', but uses internally
    liblinear rather than libsvm, so it has more flexibility in the
    choice of penalties and loss functions and should be faster for
    huge datasets.

    Parameters
    ----------
    loss : string, 'l1' or 'l2' (default 'l2')
        Specifies the loss function. With 'l1' it is the standard SVM
        loss (a.k.a. hinge Loss) while with 'l2' it is the squared loss.
        (a.k.a. squared hinge Loss)

    penalty : string, 'l1' or 'l2' (default 'l2')
        Specifies the norm used in the penalization. The 'l2'
        penalty is the standard used in SVC. The 'l1' leads to coef_
        vectors that are sparse.

    dual : bool, (default True)
        Select the algorithm to either solve the dual or primal
        optimization problem.

    tol: float, optional
         tolerance for stopping criteria

    multi_class: boolean, optional
         perform multi-class SVM by Cramer and Singer. If active,
         options loss, penalty and dual will be ignored.

    intercept_scaling : float, default: 1
        when self.fit_intercept is True, instance vector x becomes
        [x, self.intercept_scaling],
        i.e. a "synthetic" feature with constant value equals to
        intercept_scaling is appended to the instance vector.
        The intercept becomes intercept_scaling * synthetic feature weight
        Note! the synthetic feature weight is subject to l1/l2 regularization
        as all other features.
        To lessen the effect of regularization on synthetic feature weight
        (and therefore on the intercept) intercept_scaling has to be increased

    Attributes
    ----------
    `coef_` : array, shape = [n_features] if n_classes == 2 else [n_classes, n_features]
        Weights asigned to the features (coefficients in the primal
        problem). This is only available in the case of linear kernel.

    `intercept_` : array, shape = [1] if n_classes == 2 else [n_classes]
        Constants in decision function.

    Notes
    -----
    The underlying C implementation uses a random number generator to
    select features when fitting the model. It is thus not uncommon,
    to have slightly different results for the same input data. If
    that happens, try with a smaller tol parameter.

    See also
    --------
    SVC

    References
    ----------
    LIBLINEAR -- A Library for Large Linear Classification
    http://www.csie.ntu.edu.tw/~cjlin/liblinear/

    """

    # all the implementation is provided by the mixins
    pass


class SVC(BaseLibSVM, ClassifierMixin):
    """C-Support Vector Classification.

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

    tol: float, optional
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
    SVC(kernel='rbf', C=1.0, probability=False, degree=3, coef0=0.0, tol=0.001,
      cache_size=100.0, shrinking=True, gamma=0.25)
    >>> print clf.predict([[-0.8, -1]])
    [ 1.]

    See also
    --------
    SVR, LinearSVC
    """

    def __init__(self, C=1.0, kernel='rbf', degree=3, gamma=0.0,
                 coef0=0.0, shrinking=True, probability=False,
                 tol=1e-3, cache_size=100.0):

        BaseLibSVM.__init__(self, 'c_svc', kernel, degree, gamma, coef0,
                         cache_size, tol, C, 0., 0.,
                         shrinking, probability)


class NuSVC(BaseLibSVM, ClassifierMixin):
    """Nu-Support Vector Classification.

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

    tol: float, optional
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

    predict_log_proba(X) : array
        Return log-probability estimates.

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
    NuSVC(kernel='rbf', probability=False, degree=3, coef0=0.0, tol=0.001,
       cache_size=100.0, shrinking=True, nu=0.5, gamma=0.25)
    >>> print clf.predict([[-0.8, -1]])
    [ 1.]

    See also
    --------
    SVC, LinearSVC, SVR
    """

    def __init__(self, nu=0.5, kernel='rbf', degree=3, gamma=0.0,
                 coef0=0.0, shrinking=True, probability=False,
                 tol=1e-3, cache_size=100.0):

        BaseLibSVM.__init__(self, 'nu_svc', kernel, degree, gamma,
                         coef0, cache_size, tol, 0., nu, 0.,
                         shrinking, probability)


class SVR(BaseLibSVM, RegressorMixin):
    """Support Vector Regression.

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

    tol: float, optional
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

    Examples
    --------
    >>> from scikits.learn.svm import SVR
    >>> import numpy as np
    >>> n_samples, n_features = 10, 5
    >>> np.random.seed(0)
    >>> y = np.random.randn(n_samples)
    >>> X = np.random.randn(n_samples, n_features)
    >>> clf = SVR(C=1.0, p=0.2)
    >>> clf.fit(X, y)
    SVR(kernel='rbf', C=1.0, probability=False, degree=3, shrinking=True, p=0.2,
      tol=0.001, cache_size=100.0, coef0=0.0, nu=0.5, gamma=0.1)

    See also
    --------
    NuSVR
    """
    def __init__(self, kernel='rbf', degree=3, gamma=0.0, coef0=0.0,
                 cache_size=100.0, tol=1e-3, C=1.0, nu=0.5, p=0.1,
                 shrinking=True, probability=False):

        BaseLibSVM.__init__(self, 'epsilon_svr', kernel, degree, gamma, coef0,
                         cache_size, tol, C, nu, p,
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
    """Nu Support Vector Regression.

    Similar to NuSVC, for regression, uses a paramter nu to control
    the number of support vectors. However, unlike NuSVC, where nu
    replaces with C, here nu replaces with the parameter p of SVR.

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

    tol: float, optional
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

    Examples
    --------
    >>> from scikits.learn.svm import NuSVR
    >>> import numpy as np
    >>> n_samples, n_features = 10, 5
    >>> np.random.seed(0)
    >>> y = np.random.randn(n_samples)
    >>> X = np.random.randn(n_samples, n_features)
    >>> clf = NuSVR(nu=0.1, C=1.0)
    >>> clf.fit(X, y)
    NuSVR(kernel='rbf', C=1.0, probability=False, degree=3, shrinking=True,
       tol=0.001, cache_size=100.0, coef0=0.0, nu=0.1, gamma=0.1)

    See also
    --------
    NuSVR
    """

    def __init__(self, nu=0.5, C=1.0, kernel='rbf', degree=3,
                 gamma=0.0, coef0=0.0, shrinking=True,
                 probability=False, cache_size=100.0, tol=1e-3):

        BaseLibSVM.__init__(self, 'epsilon_svr', kernel, degree, gamma, coef0,
                         cache_size, tol, C, nu, 0.,
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
    """Unsupervised Outliers Detection.

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

    tol: float, optional
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
                 cache_size=100.0, tol=1e-3, nu=0.5, shrinking=True):
        BaseLibSVM.__init__(self, 'one_class', kernel, degree, gamma, coef0,
                             cache_size, tol, 0.0, nu, 0.0, shrinking, False)

    def fit(self, X, class_weight={}, sample_weight=[], **params):
        """
        Detects the soft boundary of the set of samples X.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Set of samples, where n_samples is the number of samples and
            n_features is the number of features.

        Returns
        -------
        self : object
            Returns self.
        """
        super(OneClassSVM, self).fit(
            X, [], class_weight=class_weight, sample_weight=sample_weight,
            **params)
