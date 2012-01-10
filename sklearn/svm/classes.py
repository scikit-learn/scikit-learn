from ..base import ClassifierMixin, RegressorMixin
from ..linear_model.base import CoefSelectTransformerMixin
from .base import BaseLibLinear, DenseBaseLibSVM


class LinearSVC(BaseLibLinear, ClassifierMixin, CoefSelectTransformerMixin):
    """Linear Support Vector Classification.

    Similar to SVC with parameter kernel='linear', but uses internally
    liblinear rather than libsvm, so it has more flexibility in the
    choice of penalties and loss functions and should be faster for
    huge datasets.

    Parameters
    ----------
    C : float, optional (default=1.0)
        Penalty parameter C of the error term.

    loss : string, 'l1' or 'l2' (default='l2')
        Specifies the loss function. 'l1' is the hinge loss (standard SVM)
        while 'l2' is the squared hinge loss.

    penalty : string, 'l1' or 'l2' (default='l2')
        Specifies the norm used in the penalization. The 'l2'
        penalty is the standard used in SVC. The 'l1' leads to `coef_`
        vectors that are sparse.

    dual : bool, (default=True)
        Select the algorithm to either solve the dual or primal
        optimization problem. Prefer dual=False when n_samples > n_features.

    tol: float, optional (default=1e-4)
        Tolerance for stopping criteria

    multi_class: boolean, optional (default=False)
        Perform multi-class SVM as per Cramer and Singer. If active,
        the options loss, penalty and dual will be ignored.

    fit_intercept : boolean, optional (default=True)
        Whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    intercept_scaling : float, optional (default=1)
        when self.fit_intercept is True, instance vector x becomes
        [x, self.intercept_scaling],
        i.e. a "synthetic" feature with constant value equals to
        intercept_scaling is appended to the instance vector.
        The intercept becomes intercept_scaling * synthetic feature weight
        Note! the synthetic feature weight is subject to l1/l2 regularization
        as all other features.
        To lessen the effect of regularization on synthetic feature weight
        (and therefore on the intercept) intercept_scaling has to be increased

    scale_C : bool
        Scale C with number of samples. It makes the setting of C independent
        of the number of samples.

    Attributes
    ----------
    `coef_` : array, shape = [n_features] if n_classes == 2 \
            else [n_classes, n_features]
        Weights asigned to the features (coefficients in the primal
        problem). This is only available in the case of linear kernel.

        `coef_` is readonly property derived from `raw_coef_` that \
        follows the internal memory layout of liblinear.

    `intercept_` : array, shape = [1] if n_classes == 2 else [n_classes]
        Constants in decision function.

    Notes
    -----
    The underlying C implementation uses a random number generator to
    select features when fitting the model. It is thus not uncommon,
    to have slightly different results for the same input data. If
    that happens, try with a smaller tol parameter.

    **References:**
    `LIBLINEAR: A Library for Large Linear Classification
    <http://www.csie.ntu.edu.tw/~cjlin/liblinear/>`__

    See also
    --------
    SVC


    """

    # all the implementation is provided by the mixins
    pass


class SVC(DenseBaseLibSVM, ClassifierMixin):
    """C-Support Vector Classification.

    Parameters
    ----------
    C : float, optional (default=1.0)
        Penalty parameter C of the error term.

    kernel : string, optional (default='rbf')
         Specifies the kernel type to be used in the algorithm.
         It must be one of 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed'.
         If none is given, 'rbf' will be used.

    degree : int, optional (default=3)
        Degree of kernel function.
        It is significant only in 'poly' and 'sigmoid'.

    gamma : float, optional (default=0.0)
        Kernel coefficient for 'rbf' and 'poly'.
        If gamma is 0.0 then 1/n_features will be used instead.

    coef0 : float, optional (default=0.0)
        Independent term in kernel function.
        It is only significant in 'poly' and 'sigmoid'.

    probability: boolean, optional (default=False)
        Whether to enable probability estimates. This must be enabled prior
        to calling predict_proba.

    shrinking: boolean, optional (default=True)
        Whether to use the shrinking heuristic.

    tol: float, optional (default=1e-3)
        Tolerance for stopping criterion.

    cache_size: float, optional
        Specify the size of the kernel cache (in MB)

    scale_C : bool
        Scale C with number of samples. It makes the setting of C independant
        of the number of samples.

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

        `coef_` is readonly property derived from `dual_coef_` and
        `support_vectors_`

    `intercept_` : array, shape = [n_class * (n_class-1) / 2]
        Constants in decision function.

    Examples
    --------
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    >>> y = np.array([1, 1, 2, 2])
    >>> from sklearn.svm import SVC
    >>> clf = SVC()
    >>> clf.fit(X, y)
    SVC(C=1.0, cache_size=200, coef0=0.0, degree=3, gamma=0.5, kernel='rbf',
      probability=False, scale_C=False, shrinking=True, tol=0.001)
    >>> print clf.predict([[-0.8, -1]])
    [ 1.]

    See also
    --------
    SVR, LinearSVC
    """

    def __init__(self, C=1.0, kernel='rbf', degree=3, gamma=0.0,
                 coef0=0.0, shrinking=True, probability=False,
                 tol=1e-3, cache_size=200, scale_C=False):

        super(SVC, self).__init__('c_svc', kernel, degree, gamma, coef0, tol,
                                  C, 0., 0., shrinking, probability,
                                  cache_size, scale_C)


class NuSVC(DenseBaseLibSVM, ClassifierMixin):
    """Nu-Support Vector Classification.

    Parameters
    ----------
    nu : float, optional (default=0.5)
        An upper bound on the fraction of training errors and a lower
        bound of the fraction of support vectors. Should be in the
        interval (0, 1].

    kernel : string, optional (default='rbf')
         Specifies the kernel type to be used in the algorithm.
         one of 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed'.
         If none is given 'rbf' will be used.

    degree : int, optional (default=3)
        degree of kernel function
        is significant only in poly, rbf, sigmoid

    gamma : float, optional (default=0.0)
        kernel coefficient for rbf and poly, if gamma is 0.0 then 1/n_features
        will be taken.

    coef0 : float, optional (default=0.0)
        independent term in kernel function. It is only significant
        in poly/sigmoid.

    probability: boolean, optional (default=False)
        Whether to enable probability estimates. This must be enabled prior
        to calling predict_proba.

    shrinking: boolean, optional (default=True)
        Whether to use the shrinking heuristic.

    tol: float, optional (default=1e-3)
        Tolerance for stopping criterion.

    cache_size: float, optional
        Specify the size of the kernel cache (in MB)

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

        `coef_` is readonly property derived from `dual_coef_` and
        `support_vectors_`

    `intercept_` : array, shape = [n_class * (n_class-1) / 2]
        Constants in decision function.

    Examples
    --------
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    >>> y = np.array([1, 1, 2, 2])
    >>> from sklearn.svm import NuSVC
    >>> clf = NuSVC()
    >>> clf.fit(X, y)
    NuSVC(cache_size=200, coef0=0.0, degree=3, gamma=0.5, kernel='rbf', nu=0.5,
       probability=False, shrinking=True, tol=0.001)
    >>> print clf.predict([[-0.8, -1]])
    [ 1.]

    See also
    --------
    SVC, LinearSVC, SVR
    """

    def __init__(self, nu=0.5, kernel='rbf', degree=3, gamma=0.0,
                 coef0=0.0, shrinking=True, probability=False,
                 tol=1e-3, cache_size=200):

        super(NuSVC, self).__init__('nu_svc', kernel, degree, gamma, coef0,
                                    tol, 0., nu, 0., shrinking, probability,
                                    cache_size, scale_C=None)


class SVR(DenseBaseLibSVM, RegressorMixin):
    """epsilon-Support Vector Regression.

    The free parameters in the model are C and epsilon.

    Parameters
    ----------
    C : float, optional (default=1.0)
        penalty parameter C of the error term.

    epsilon : float, optional (default=0.1)
         epsilon in the epsilon-SVR model. It specifies the epsilon-tube
         within which no penalty is associated in the training loss function
         with points predicted within a distance epsilon from the actual
         value.

    kernel : string, optional (default='rbf')
         Specifies the kernel type to be used in the algorithm.
         one of 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed'.
         If none is given 'rbf' will be used.

    degree : int, optional (default=3)
        degree of kernel function
        is significant only in poly, rbf, sigmoid

    gamma : float, optional (default=0.0)
        kernel coefficient for rbf and poly, if gamma is 0.0 then 1/n_features
        will be taken.

    coef0 : float, optional (default=0.0)
        independent term in kernel function. It is only significant
        in poly/sigmoid.

    probability: boolean, optional (default=False)
        Whether to enable probability estimates. This must be enabled prior
        to calling predict_proba.

    shrinking: boolean, optional (default=True)
        Whether to use the shrinking heuristic.

    tol: float, optional (default=1e-3)
        Tolerance for stopping criterion.

    cache_size: float, optional
        Specify the size of the kernel cache (in MB)

    scale_C : bool
        Scale C with number of samples. It makes the setting of C independant
        of the number of samples.

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

        `coef_` is readonly property derived from `dual_coef_` and
        `support_vectors_`

    `intercept_` : array, shape = [n_class * (n_class-1) / 2]
        Constants in decision function.

    Examples
    --------
    >>> from sklearn.svm import SVR
    >>> import numpy as np
    >>> n_samples, n_features = 10, 5
    >>> np.random.seed(0)
    >>> y = np.random.randn(n_samples)
    >>> X = np.random.randn(n_samples, n_features)
    >>> clf = SVR(C=1.0, epsilon=0.2)
    >>> clf.fit(X, y)
    SVR(C=1.0, cache_size=200, coef0=0.0, degree=3, epsilon=0.2, gamma=0.2,
      kernel='rbf', probability=False, scale_C=False, shrinking=True,
      tol=0.001)

    See also
    --------
    NuSVR
    """
    def __init__(self, kernel='rbf', degree=3, gamma=0.0, coef0=0.0,
                 tol=1e-3, C=1.0, epsilon=0.1, shrinking=True,
                 probability=False, cache_size=200, scale_C=False):

        super(SVR, self).__init__('epsilon_svr', kernel, degree, gamma, coef0,
                                  tol, C, 0., epsilon, shrinking, probability,
                                  cache_size, scale_C)

    def fit(self, X, y, sample_weight=None, **params):
        """
        Fit the SVM model according to the given training data and parameters.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vector, where n_samples is the number of samples and
            n_features is the number of features.
        y : array, shape = [n_samples]
            Target values. Array of floating-point numbers.
        cache_size: float, optional
            Specify the size of the cache (in MB)


        Returns
        -------
        self : object
            Returns self.
        """
        # we copy this method because SVR does not accept class_weight
        return super(SVR, self).fit(X, y, sample_weight=sample_weight,
                                    **params)


class NuSVR(DenseBaseLibSVM, RegressorMixin):
    """Nu Support Vector Regression.

    Similar to NuSVC, for regression, uses a parameter nu to control
    the number of support vectors. However, unlike NuSVC, where nu
    replaces C, here nu replaces with the parameter epsilon of SVR.

    Parameters
    ----------
    C : float, optional (default=1.0)
        penalty parameter C of the error term.

    nu : float, optional
        An upper bound on the fraction of training errors and a lower bound of
        the fraction of support vectors. Should be in the interval (0, 1].  By
        default 0.5 will be taken.  Only available if impl='nu_svc'.

    kernel : string, optional (default='rbf')
         Specifies the kernel type to be used in the algorithm.
         one of 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed'.
         If none is given 'rbf' will be used.

    degree : int, optional (default=3)
        degree of kernel function
        is significant only in poly, rbf, sigmoid

    gamma : float, optional (default=0.0)
        kernel coefficient for rbf and poly, if gamma is 0.0 then 1/n_features
        will be taken.

    coef0 : float, optional (default=0.0)
        independent term in kernel function. It is only significant
        in poly/sigmoid.

    probability: boolean, optional (default=False)
        Whether to enable probability estimates. This must be enabled prior
        to calling predict_proba.

    shrinking: boolean, optional (default=True)
        Whether to use the shrinking heuristic.

    tol: float, optional (default=1e-3)
        Tolerance for stopping criterion.

    cache_size: float, optional
        Specify the size of the kernel cache (in MB)

    scale_C : bool
        Scale C with number of samples. It makes the setting of C independant
        of the number of samples.

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

        `coef_` is readonly property derived from `dual_coef_` and
        `support_vectors_`

    `intercept_` : array, shape = [n_class * (n_class-1) / 2]
        Constants in decision function.

    Examples
    --------
    >>> from sklearn.svm import NuSVR
    >>> import numpy as np
    >>> n_samples, n_features = 10, 5
    >>> np.random.seed(0)
    >>> y = np.random.randn(n_samples)
    >>> X = np.random.randn(n_samples, n_features)
    >>> clf = NuSVR(C=1.0, nu=0.1)
    >>> clf.fit(X, y)
    NuSVR(C=1.0, cache_size=200, coef0=0.0, degree=3, gamma=0.2, kernel='rbf',
       nu=0.1, probability=False, scale_C=False, shrinking=True, tol=0.001)

    See also
    --------
    NuSVC, SVR
    """

    def __init__(self, nu=0.5, C=1.0, kernel='rbf', degree=3,
                 gamma=0.0, coef0=0.0, shrinking=True,
                 probability=False, tol=1e-3, cache_size=200,
                 scale_C=False):

        super(NuSVR, self).__init__('nu_svr', kernel, degree, gamma, coef0,
                                    tol, C, nu, None, shrinking, probability,
                                    cache_size, scale_C=scale_C)

    def fit(self, X, y, sample_weight=None, **params):
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
        return super(NuSVR, self).fit(X, y, sample_weight=[], **params)


class OneClassSVM(DenseBaseLibSVM):
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

    gamma : float, optional (default=0.0)
        kernel coefficient for rbf and poly, if gamma is 0.0 then 1/n_features
        will be taken.

    coef0 : float, optional
        Independent term in kernel function. It is only significant in
        poly/sigmoid.

    tol: float, optional
        Tolerance for stopping criterion.

    shrinking: boolean, optional
        Whether to use the shrinking heuristic.

    cache_size: float, optional
        Specify the size of the kernel cache (in MB)

    scale_C : bool
        Scale C with number of samples. It makes the setting of C independant
        of the number of samples.


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

        `coef_` is readonly property derived from `dual_coef_` and
        `support_vectors_`

    `intercept_` : array, shape = [n_classes-1]
        Constants in decision function.

    """
    def __init__(self, kernel='rbf', degree=3, gamma=0.0, coef0=0.0,
                 tol=1e-3, nu=0.5, shrinking=True, cache_size=200):

        super(OneClassSVM, self).__init__('one_class', kernel, degree, gamma,
                                          coef0, tol, 0., nu, 0., shrinking,
                                          False, cache_size, scale_C=None)

    def fit(self, X, class_weight={}, sample_weight=None, **params):
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

        Notes
        -----
        If X is not a C-ordered contiguous array, it is copied.

        """
        super(OneClassSVM, self).fit(
            X, [], class_weight=class_weight, sample_weight=sample_weight,
            **params)
        return self
