import warnings
import numpy as np

from .base import _fit_liblinear, BaseSVC, BaseLibSVM
from ..base import BaseEstimator, RegressorMixin, OutlierMixin
from ..linear_model.base import LinearClassifierMixin, SparseCoefMixin, \
    LinearModel
from ..utils import check_X_y
from ..utils.validation import _num_samples
from ..utils.multiclass import check_classification_targets


class LinearSVC(BaseEstimator, LinearClassifierMixin,
                SparseCoefMixin):
    """Linear Support Vector Classification.

    Similar to SVC with parameter kernel='linear', but implemented in terms of
    liblinear rather than libsvm, so it has more flexibility in the choice of
    penalties and loss functions and should scale better to large numbers of
    samples.

    This class supports both dense and sparse input and the multiclass support
    is handled according to a one-vs-the-rest scheme.

    Read more in the :ref:`User Guide <svm_classification>`.

    Parameters
    ----------
    penalty : string, 'l1' or 'l2' (default='l2')
        Specifies the norm used in the penalization. The 'l2'
        penalty is the standard used in SVC. The 'l1' leads to ``coef_``
        vectors that are sparse.

    loss : string, 'hinge' or 'squared_hinge' (default='squared_hinge')
        Specifies the loss function. 'hinge' is the standard SVM loss
        (used e.g. by the SVC class) while 'squared_hinge' is the
        square of the hinge loss.

    dual : bool, (default=True)
        Select the algorithm to either solve the dual or primal
        optimization problem. Prefer dual=False when n_samples > n_features.

    tol : float, optional (default=1e-4)
        Tolerance for stopping criteria.

    C : float, optional (default=1.0)
        Penalty parameter C of the error term.

    multi_class : string, 'ovr' or 'crammer_singer' (default='ovr')
        Determines the multi-class strategy if `y` contains more than
        two classes.
        ``"ovr"`` trains n_classes one-vs-rest classifiers, while
        ``"crammer_singer"`` optimizes a joint objective over all classes.
        While `crammer_singer` is interesting from a theoretical perspective
        as it is consistent, it is seldom used in practice as it rarely leads
        to better accuracy and is more expensive to compute.
        If ``"crammer_singer"`` is chosen, the options loss, penalty and dual
        will be ignored.

    fit_intercept : boolean, optional (default=True)
        Whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (i.e. data is expected to be already centered).

    intercept_scaling : float, optional (default=1)
        When self.fit_intercept is True, instance vector x becomes
        ``[x, self.intercept_scaling]``,
        i.e. a "synthetic" feature with constant value equals to
        intercept_scaling is appended to the instance vector.
        The intercept becomes intercept_scaling * synthetic feature weight
        Note! the synthetic feature weight is subject to l1/l2 regularization
        as all other features.
        To lessen the effect of regularization on synthetic feature weight
        (and therefore on the intercept) intercept_scaling has to be increased.

    class_weight : {dict, 'balanced'}, optional
        Set the parameter C of class i to ``class_weight[i]*C`` for
        SVC. If not given, all classes are supposed to have
        weight one.
        The "balanced" mode uses the values of y to automatically adjust
        weights inversely proportional to class frequencies in the input data
        as ``n_samples / (n_classes * np.bincount(y))``

    verbose : int, (default=0)
        Enable verbose output. Note that this setting takes advantage of a
        per-process runtime setting in liblinear that, if enabled, may not work
        properly in a multithreaded context.

    random_state : int, RandomState instance or None, optional (default=None)
        The seed of the pseudo random number generator to use when shuffling
        the data for the dual coordinate descent (if ``dual=True``). When
        ``dual=False`` the underlying implementation of :class:`LinearSVC`
        is not random and ``random_state`` has no effect on the results. If
        int, random_state is the seed used by the random number generator; If
        RandomState instance, random_state is the random number generator; If
        None, the random number generator is the RandomState instance used by
        `np.random`.

    max_iter : int, (default=1000)
        The maximum number of iterations to be run.

    Attributes
    ----------
    coef_ : array, shape = [n_features] if n_classes == 2 else [n_classes, n_features]
        Weights assigned to the features (coefficients in the primal
        problem). This is only available in the case of a linear kernel.

        ``coef_`` is a readonly property derived from ``raw_coef_`` that
        follows the internal memory layout of liblinear.

    intercept_ : array, shape = [1] if n_classes == 2 else [n_classes]
        Constants in decision function.

    Examples
    --------
    >>> from sklearn.svm import LinearSVC
    >>> from sklearn.datasets import make_classification
    >>> X, y = make_classification(n_features=4, random_state=0)
    >>> clf = LinearSVC(random_state=0, tol=1e-5)
    >>> clf.fit(X, y)  # doctest: +NORMALIZE_WHITESPACE
    LinearSVC(C=1.0, class_weight=None, dual=True, fit_intercept=True,
         intercept_scaling=1, loss='squared_hinge', max_iter=1000,
         multi_class='ovr', penalty='l2', random_state=0, tol=1e-05, verbose=0)
    >>> print(clf.coef_)
    [[0.085... 0.394... 0.498... 0.375...]]
    >>> print(clf.intercept_)
    [0.284...]
    >>> print(clf.predict([[0, 0, 0, 0]]))
    [1]

    Notes
    -----
    The underlying C implementation uses a random number generator to
    select features when fitting the model. It is thus not uncommon
    to have slightly different results for the same input data. If
    that happens, try with a smaller ``tol`` parameter.

    The underlying implementation, liblinear, uses a sparse internal
    representation for the data that will incur a memory copy.

    Predict output may not match that of standalone liblinear in certain
    cases. See :ref:`differences from liblinear <liblinear_differences>`
    in the narrative documentation.

    References
    ----------
    `LIBLINEAR: A Library for Large Linear Classification
    <https://www.csie.ntu.edu.tw/~cjlin/liblinear/>`__

    See also
    --------
    SVC
        Implementation of Support Vector Machine classifier using libsvm:
        the kernel can be non-linear but its SMO algorithm does not
        scale to large number of samples as LinearSVC does.

        Furthermore SVC multi-class mode is implemented using one
        vs one scheme while LinearSVC uses one vs the rest. It is
        possible to implement one vs the rest with SVC by using the
        :class:`sklearn.multiclass.OneVsRestClassifier` wrapper.

        Finally SVC can fit dense data without memory copy if the input
        is C-contiguous. Sparse data will still incur memory copy though.

    sklearn.linear_model.SGDClassifier
        SGDClassifier can optimize the same cost function as LinearSVC
        by adjusting the penalty and loss parameters. In addition it requires
        less memory, allows incremental (online) learning, and implements
        various loss functions and regularization regimes.

    """

    def __init__(self, penalty='l2', loss='squared_hinge', dual=True, tol=1e-4,
                 C=1.0, multi_class='ovr', fit_intercept=True,
                 intercept_scaling=1, class_weight=None, verbose=0,
                 random_state=None, max_iter=1000):
        self.dual = dual
        self.tol = tol
        self.C = C
        self.multi_class = multi_class
        self.fit_intercept = fit_intercept
        self.intercept_scaling = intercept_scaling
        self.class_weight = class_weight
        self.verbose = verbose
        self.random_state = random_state
        self.max_iter = max_iter
        self.penalty = penalty
        self.loss = loss

    def fit(self, X, y, sample_weight=None):
        """Fit the model according to the given training data.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target vector relative to X

        sample_weight : array-like, shape = [n_samples], optional
            Array of weights that are assigned to individual
            samples. If not provided,
            then each sample is given unit weight.

        Returns
        -------
        self : object
        """
        # FIXME Remove l1/l2 support in 1.0 -----------------------------------
        msg = ("loss='%s' has been deprecated in favor of "
               "loss='%s' as of 0.16. Backward compatibility"
               " for the loss='%s' will be removed in %s")

        if self.loss in ('l1', 'l2'):
            old_loss = self.loss
            self.loss = {'l1': 'hinge', 'l2': 'squared_hinge'}.get(self.loss)
            warnings.warn(msg % (old_loss, self.loss, old_loss, '1.0'),
                          DeprecationWarning)
        # ---------------------------------------------------------------------

        if self.C < 0:
            raise ValueError("Penalty term must be positive; got (C=%r)"
                             % self.C)

        X, y = check_X_y(X, y, accept_sparse='csr',
                         dtype=np.float64, order="C",
                         accept_large_sparse=False)
        check_classification_targets(y)
        self.classes_ = np.unique(y)

        self.coef_, self.intercept_, self.n_iter_ = _fit_liblinear(
            X, y, self.C, self.fit_intercept, self.intercept_scaling,
            self.class_weight, self.penalty, self.dual, self.verbose,
            self.max_iter, self.tol, self.random_state, self.multi_class,
            self.loss, sample_weight=sample_weight)

        if self.multi_class == "crammer_singer" and len(self.classes_) == 2:
            self.coef_ = (self.coef_[1] - self.coef_[0]).reshape(1, -1)
            if self.fit_intercept:
                intercept = self.intercept_[1] - self.intercept_[0]
                self.intercept_ = np.array([intercept])

        return self


class LinearSVR(LinearModel, RegressorMixin):
    """Linear Support Vector Regression.

    Similar to SVR with parameter kernel='linear', but implemented in terms of
    liblinear rather than libsvm, so it has more flexibility in the choice of
    penalties and loss functions and should scale better to large numbers of
    samples.

    This class supports both dense and sparse input.

    Read more in the :ref:`User Guide <svm_regression>`.

    Parameters
    ----------
    epsilon : float, optional (default=0.0)
        Epsilon parameter in the epsilon-insensitive loss function. Note
        that the value of this parameter depends on the scale of the target
        variable y. If unsure, set ``epsilon=0``.

    tol : float, optional (default=1e-4)
        Tolerance for stopping criteria.

    C : float, optional (default=1.0)
        Penalty parameter C of the error term. The penalty is a squared
        l2 penalty. The bigger this parameter, the less regularization is used.

    loss : string, optional (default='epsilon_insensitive')
        Specifies the loss function. The epsilon-insensitive loss
        (standard SVR) is the L1 loss, while the squared epsilon-insensitive
        loss ('squared_epsilon_insensitive') is the L2 loss.

    fit_intercept : boolean, optional (default=True)
        Whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (i.e. data is expected to be already centered).

    intercept_scaling : float, optional (default=1)
        When self.fit_intercept is True, instance vector x becomes
        [x, self.intercept_scaling],
        i.e. a "synthetic" feature with constant value equals to
        intercept_scaling is appended to the instance vector.
        The intercept becomes intercept_scaling * synthetic feature weight
        Note! the synthetic feature weight is subject to l1/l2 regularization
        as all other features.
        To lessen the effect of regularization on synthetic feature weight
        (and therefore on the intercept) intercept_scaling has to be increased.

    dual : bool, (default=True)
        Select the algorithm to either solve the dual or primal
        optimization problem. Prefer dual=False when n_samples > n_features.

    verbose : int, (default=0)
        Enable verbose output. Note that this setting takes advantage of a
        per-process runtime setting in liblinear that, if enabled, may not work
        properly in a multithreaded context.

    random_state : int, RandomState instance or None, optional (default=None)
        The seed of the pseudo random number generator to use when shuffling
        the data.  If int, random_state is the seed used by the random number
        generator; If RandomState instance, random_state is the random number
        generator; If None, the random number generator is the RandomState
        instance used by `np.random`.

    max_iter : int, (default=1000)
        The maximum number of iterations to be run.

    Attributes
    ----------
    coef_ : array, shape = [n_features] if n_classes == 2 else [n_classes, n_features]
        Weights assigned to the features (coefficients in the primal
        problem). This is only available in the case of a linear kernel.

        `coef_` is a readonly property derived from `raw_coef_` that
        follows the internal memory layout of liblinear.

    intercept_ : array, shape = [1] if n_classes == 2 else [n_classes]
        Constants in decision function.

    Examples
    --------
    >>> from sklearn.svm import LinearSVR
    >>> from sklearn.datasets import make_regression
    >>> X, y = make_regression(n_features=4, random_state=0)
    >>> regr = LinearSVR(random_state=0, tol=1e-5)
    >>> regr.fit(X, y)  # doctest: +NORMALIZE_WHITESPACE
    LinearSVR(C=1.0, dual=True, epsilon=0.0, fit_intercept=True,
         intercept_scaling=1.0, loss='epsilon_insensitive', max_iter=1000,
         random_state=0, tol=1e-05, verbose=0)
    >>> print(regr.coef_)
    [16.35... 26.91... 42.30... 60.47...]
    >>> print(regr.intercept_)
    [-4.29...]
    >>> print(regr.predict([[0, 0, 0, 0]]))
    [-4.29...]

    See also
    --------
    LinearSVC
        Implementation of Support Vector Machine classifier using the
        same library as this class (liblinear).

    SVR
        Implementation of Support Vector Machine regression using libsvm:
        the kernel can be non-linear but its SMO algorithm does not
        scale to large number of samples as LinearSVC does.

    sklearn.linear_model.SGDRegressor
        SGDRegressor can optimize the same cost function as LinearSVR
        by adjusting the penalty and loss parameters. In addition it requires
        less memory, allows incremental (online) learning, and implements
        various loss functions and regularization regimes.
    """

    def __init__(self, epsilon=0.0, tol=1e-4, C=1.0,
                 loss='epsilon_insensitive', fit_intercept=True,
                 intercept_scaling=1., dual=True, verbose=0,
                 random_state=None, max_iter=1000):
        self.tol = tol
        self.C = C
        self.epsilon = epsilon
        self.fit_intercept = fit_intercept
        self.intercept_scaling = intercept_scaling
        self.verbose = verbose
        self.random_state = random_state
        self.max_iter = max_iter
        self.dual = dual
        self.loss = loss

    def fit(self, X, y, sample_weight=None):
        """Fit the model according to the given training data.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target vector relative to X

        sample_weight : array-like, shape = [n_samples], optional
            Array of weights that are assigned to individual
            samples. If not provided,
            then each sample is given unit weight.

        Returns
        -------
        self : object
        """
        # FIXME Remove l1/l2 support in 1.0 -----------------------------------
        msg = ("loss='%s' has been deprecated in favor of "
               "loss='%s' as of 0.16. Backward compatibility"
               " for the loss='%s' will be removed in %s")

        if self.loss in ('l1', 'l2'):
            old_loss = self.loss
            self.loss = {'l1': 'epsilon_insensitive',
                         'l2': 'squared_epsilon_insensitive'
                         }.get(self.loss)
            warnings.warn(msg % (old_loss, self.loss, old_loss, '1.0'),
                          DeprecationWarning)
        # ---------------------------------------------------------------------

        if self.C < 0:
            raise ValueError("Penalty term must be positive; got (C=%r)"
                             % self.C)

        X, y = check_X_y(X, y, accept_sparse='csr',
                         dtype=np.float64, order="C",
                         accept_large_sparse=False)
        penalty = 'l2'  # SVR only accepts l2 penalty
        self.coef_, self.intercept_, self.n_iter_ = _fit_liblinear(
            X, y, self.C, self.fit_intercept, self.intercept_scaling,
            None, penalty, self.dual, self.verbose,
            self.max_iter, self.tol, self.random_state, loss=self.loss,
            epsilon=self.epsilon, sample_weight=sample_weight)
        self.coef_ = self.coef_.ravel()

        return self


class SVC(BaseSVC):
    """C-Support Vector Classification.

    The implementation is based on libsvm. The fit time scales at least
    quadratically with the number of samples and may be impractical
    beyond tens of thousands of samples. For large datasets
    consider using :class:`sklearn.linear_model.LinearSVC` or
    :class:`sklearn.linear_model.SGDClassifier` instead, possibly after a
    :class:`sklearn.kernel_approximation.Nystroem` transformer.

    The multiclass support is handled according to a one-vs-one scheme.

    For details on the precise mathematical formulation of the provided
    kernel functions and how `gamma`, `coef0` and `degree` affect each
    other, see the corresponding section in the narrative documentation:
    :ref:`svm_kernels`.

    Read more in the :ref:`User Guide <svm_classification>`.

    Parameters
    ----------
    C : float, optional (default=1.0)
        Penalty parameter C of the error term.

    kernel : string, optional (default='rbf')
        Specifies the kernel type to be used in the algorithm.
        It must be one of 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed' or
        a callable.
        If none is given, 'rbf' will be used. If a callable is given it is
        used to pre-compute the kernel matrix from data matrices; that matrix
        should be an array of shape ``(n_samples, n_samples)``.

    degree : int, optional (default=3)
        Degree of the polynomial kernel function ('poly').
        Ignored by all other kernels.

    gamma : float, optional (default='auto')
        Kernel coefficient for 'rbf', 'poly' and 'sigmoid'.

        Current default is 'auto' which uses 1 / n_features,
        if ``gamma='scale'`` is passed then it uses 1 / (n_features * X.var())
        as value of gamma. The current default of gamma, 'auto', will change
        to 'scale' in version 0.22. 'auto_deprecated', a deprecated version of
        'auto' is used as a default indicating that no explicit value of gamma
        was passed.

    coef0 : float, optional (default=0.0)
        Independent term in kernel function.
        It is only significant in 'poly' and 'sigmoid'.

    shrinking : boolean, optional (default=True)
        Whether to use the shrinking heuristic.

    probability : boolean, optional (default=False)
        Whether to enable probability estimates. This must be enabled prior
        to calling `fit`, and will slow down that method.

    tol : float, optional (default=1e-3)
        Tolerance for stopping criterion.

    cache_size : float, optional
        Specify the size of the kernel cache (in MB).

    class_weight : {dict, 'balanced'}, optional
        Set the parameter C of class i to class_weight[i]*C for
        SVC. If not given, all classes are supposed to have
        weight one.
        The "balanced" mode uses the values of y to automatically adjust
        weights inversely proportional to class frequencies in the input data
        as ``n_samples / (n_classes * np.bincount(y))``

    verbose : bool, default: False
        Enable verbose output. Note that this setting takes advantage of a
        per-process runtime setting in libsvm that, if enabled, may not work
        properly in a multithreaded context.

    max_iter : int, optional (default=-1)
        Hard limit on iterations within solver, or -1 for no limit.

    decision_function_shape : 'ovo', 'ovr', default='ovr'
        Whether to return a one-vs-rest ('ovr') decision function of shape
        (n_samples, n_classes) as all other classifiers, or the original
        one-vs-one ('ovo') decision function of libsvm which has shape
        (n_samples, n_classes * (n_classes - 1) / 2). However, one-vs-one
        ('ovo') is always used as multi-class strategy.

        .. versionchanged:: 0.19
            decision_function_shape is 'ovr' by default.

        .. versionadded:: 0.17
           *decision_function_shape='ovr'* is recommended.

        .. versionchanged:: 0.17
           Deprecated *decision_function_shape='ovo' and None*.

    random_state : int, RandomState instance or None, optional (default=None)
        The seed of the pseudo random number generator used when shuffling
        the data for probability estimates. If int, random_state is the
        seed used by the random number generator; If RandomState instance,
        random_state is the random number generator; If None, the random
        number generator is the RandomState instance used by `np.random`.

    Attributes
    ----------
    support_ : array-like, shape = [n_SV]
        Indices of support vectors.

    support_vectors_ : array-like, shape = [n_SV, n_features]
        Support vectors.

    n_support_ : array-like, dtype=int32, shape = [n_class]
        Number of support vectors for each class.

    dual_coef_ : array, shape = [n_class-1, n_SV]
        Coefficients of the support vector in the decision function.
        For multiclass, coefficient for all 1-vs-1 classifiers.
        The layout of the coefficients in the multiclass case is somewhat
        non-trivial. See the section about multi-class classification in the
        SVM section of the User Guide for details.

    coef_ : array, shape = [n_class * (n_class-1) / 2, n_features]
        Weights assigned to the features (coefficients in the primal
        problem). This is only available in the case of a linear kernel.

        `coef_` is a readonly property derived from `dual_coef_` and
        `support_vectors_`.

    intercept_ : array, shape = [n_class * (n_class-1) / 2]
        Constants in decision function.

    fit_status_ : int
        0 if correctly fitted, 1 otherwise (will raise warning)

    probA_ : array, shape = [n_class * (n_class-1) / 2]
    probB_ : array, shape = [n_class * (n_class-1) / 2]
        If probability=True, the parameters learned in Platt scaling to
        produce probability estimates from decision values. If
        probability=False, an empty array. Platt scaling uses the logistic
        function
        ``1 / (1 + exp(decision_value * probA_ + probB_))``
        where ``probA_`` and ``probB_`` are learned from the dataset [2]_. For
        more information on the multiclass case and training procedure see
        section 8 of [1]_.

    Examples
    --------
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    >>> y = np.array([1, 1, 2, 2])
    >>> from sklearn.svm import SVC
    >>> clf = SVC(gamma='auto')
    >>> clf.fit(X, y) #doctest: +NORMALIZE_WHITESPACE
    SVC(C=1.0, cache_size=200, class_weight=None, coef0=0.0,
        decision_function_shape='ovr', degree=3, gamma='auto', kernel='rbf',
        max_iter=-1, probability=False, random_state=None, shrinking=True,
        tol=0.001, verbose=False)
    >>> print(clf.predict([[-0.8, -1]]))
    [1]

    See also
    --------
    SVR
        Support Vector Machine for Regression implemented using libsvm.

    LinearSVC
        Scalable Linear Support Vector Machine for classification
        implemented using liblinear. Check the See also section of
        LinearSVC for more comparison element.

    References
    ----------
    .. [1] `LIBSVM: A Library for Support Vector Machines
        <http://www.csie.ntu.edu.tw/~cjlin/papers/libsvm.pdf>`_

    .. [2] `Platt, John (1999). "Probabilistic outputs for support vector
        machines and comparison to regularizedlikelihood methods."
        <http://citeseer.ist.psu.edu/viewdoc/summary?doi=10.1.1.41.1639>`_
    """

    _impl = 'c_svc'

    def __init__(self, C=1.0, kernel='rbf', degree=3, gamma='auto_deprecated',
                 coef0=0.0, shrinking=True, probability=False,
                 tol=1e-3, cache_size=200, class_weight=None,
                 verbose=False, max_iter=-1, decision_function_shape='ovr',
                 random_state=None):

        super().__init__(
            kernel=kernel, degree=degree, gamma=gamma,
            coef0=coef0, tol=tol, C=C, nu=0., shrinking=shrinking,
            probability=probability, cache_size=cache_size,
            class_weight=class_weight, verbose=verbose, max_iter=max_iter,
            decision_function_shape=decision_function_shape,
            random_state=random_state)


class NuSVC(BaseSVC):
    """Nu-Support Vector Classification.

    Similar to SVC but uses a parameter to control the number of support
    vectors.

    The implementation is based on libsvm.

    Read more in the :ref:`User Guide <svm_classification>`.

    Parameters
    ----------
    nu : float, optional (default=0.5)
        An upper bound on the fraction of training errors and a lower
        bound of the fraction of support vectors. Should be in the
        interval (0, 1].

    kernel : string, optional (default='rbf')
         Specifies the kernel type to be used in the algorithm.
         It must be one of 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed' or
         a callable.
         If none is given, 'rbf' will be used. If a callable is given it is
         used to precompute the kernel matrix.

    degree : int, optional (default=3)
        Degree of the polynomial kernel function ('poly').
        Ignored by all other kernels.

    gamma : float, optional (default='auto')
        Kernel coefficient for 'rbf', 'poly' and 'sigmoid'.

        Current default is 'auto' which uses 1 / n_features,
        if ``gamma='scale'`` is passed then it uses 1 / (n_features * X.var())
        as value of gamma. The current default of gamma, 'auto', will change
        to 'scale' in version 0.22. 'auto_deprecated', a deprecated version of
        'auto' is used as a default indicating that no explicit value of gamma
        was passed.

    coef0 : float, optional (default=0.0)
        Independent term in kernel function.
        It is only significant in 'poly' and 'sigmoid'.

    shrinking : boolean, optional (default=True)
        Whether to use the shrinking heuristic.

    probability : boolean, optional (default=False)
        Whether to enable probability estimates. This must be enabled prior
        to calling `fit`, and will slow down that method.

    tol : float, optional (default=1e-3)
        Tolerance for stopping criterion.

    cache_size : float, optional
        Specify the size of the kernel cache (in MB).

    class_weight : {dict, 'balanced'}, optional
        Set the parameter C of class i to class_weight[i]*C for
        SVC. If not given, all classes are supposed to have
        weight one. The "balanced" mode uses the values of y to automatically
        adjust weights inversely proportional to class frequencies as
        ``n_samples / (n_classes * np.bincount(y))``

    verbose : bool, default: False
        Enable verbose output. Note that this setting takes advantage of a
        per-process runtime setting in libsvm that, if enabled, may not work
        properly in a multithreaded context.

    max_iter : int, optional (default=-1)
        Hard limit on iterations within solver, or -1 for no limit.

    decision_function_shape : 'ovo', 'ovr', default='ovr'
        Whether to return a one-vs-rest ('ovr') decision function of shape
        (n_samples, n_classes) as all other classifiers, or the original
        one-vs-one ('ovo') decision function of libsvm which has shape
        (n_samples, n_classes * (n_classes - 1) / 2).

        .. versionchanged:: 0.19
            decision_function_shape is 'ovr' by default.

        .. versionadded:: 0.17
           *decision_function_shape='ovr'* is recommended.

        .. versionchanged:: 0.17
           Deprecated *decision_function_shape='ovo' and None*.

    random_state : int, RandomState instance or None, optional (default=None)
        The seed of the pseudo random number generator used when shuffling
        the data for probability estimates. If int, random_state is the seed
        used by the random number generator; If RandomState instance,
        random_state is the random number generator; If None, the random
        number generator is the RandomState instance used by `np.random`.

    Attributes
    ----------
    support_ : array-like, shape = [n_SV]
        Indices of support vectors.

    support_vectors_ : array-like, shape = [n_SV, n_features]
        Support vectors.

    n_support_ : array-like, dtype=int32, shape = [n_class]
        Number of support vectors for each class.

    dual_coef_ : array, shape = [n_class-1, n_SV]
        Coefficients of the support vector in the decision function.
        For multiclass, coefficient for all 1-vs-1 classifiers.
        The layout of the coefficients in the multiclass case is somewhat
        non-trivial. See the section about multi-class classification in
        the SVM section of the User Guide for details.

    coef_ : array, shape = [n_class * (n_class-1) / 2, n_features]
        Weights assigned to the features (coefficients in the primal
        problem). This is only available in the case of a linear kernel.

        `coef_` is readonly property derived from `dual_coef_` and
        `support_vectors_`.

    intercept_ : array, shape = [n_class * (n_class-1) / 2]
        Constants in decision function.

    Examples
    --------
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    >>> y = np.array([1, 1, 2, 2])
    >>> from sklearn.svm import NuSVC
    >>> clf = NuSVC(gamma='scale')
    >>> clf.fit(X, y) #doctest: +NORMALIZE_WHITESPACE
    NuSVC(cache_size=200, class_weight=None, coef0=0.0,
          decision_function_shape='ovr', degree=3, gamma='scale', kernel='rbf',
          max_iter=-1, nu=0.5, probability=False, random_state=None,
          shrinking=True, tol=0.001, verbose=False)
    >>> print(clf.predict([[-0.8, -1]]))
    [1]

    See also
    --------
    SVC
        Support Vector Machine for classification using libsvm.

    LinearSVC
        Scalable linear Support Vector Machine for classification using
        liblinear.

    Notes
    -----
    **References:**
    `LIBSVM: A Library for Support Vector Machines
    <http://www.csie.ntu.edu.tw/~cjlin/papers/libsvm.pdf>`__
    """

    _impl = 'nu_svc'

    def __init__(self, nu=0.5, kernel='rbf', degree=3, gamma='auto_deprecated',
                 coef0=0.0, shrinking=True, probability=False, tol=1e-3,
                 cache_size=200, class_weight=None, verbose=False, max_iter=-1,
                 decision_function_shape='ovr', random_state=None):

        super().__init__(
            kernel=kernel, degree=degree, gamma=gamma,
            coef0=coef0, tol=tol, C=0., nu=nu, shrinking=shrinking,
            probability=probability, cache_size=cache_size,
            class_weight=class_weight, verbose=verbose, max_iter=max_iter,
            decision_function_shape=decision_function_shape,
            random_state=random_state)


class SVR(BaseLibSVM, RegressorMixin):
    """Epsilon-Support Vector Regression.

    The free parameters in the model are C and epsilon.

    The implementation is based on libsvm. The fit time complexity
    is more than quadratic with the number of samples which makes it hard
    to scale to datasets with more than a couple of 10000 samples. For large
    datasets consider using :class:`sklearn.linear_model.LinearSVR` or
    :class:`sklearn.linear_model.SGDRegressor` instead, possibly after a
    :class:`sklearn.kernel_approximation.Nystroem` transformer.

    Read more in the :ref:`User Guide <svm_regression>`.

    Parameters
    ----------
    kernel : string, optional (default='rbf')
         Specifies the kernel type to be used in the algorithm.
         It must be one of 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed' or
         a callable.
         If none is given, 'rbf' will be used. If a callable is given it is
         used to precompute the kernel matrix.

    degree : int, optional (default=3)
        Degree of the polynomial kernel function ('poly').
        Ignored by all other kernels.

    gamma : float, optional (default='auto')
        Kernel coefficient for 'rbf', 'poly' and 'sigmoid'.

        Current default is 'auto' which uses 1 / n_features,
        if ``gamma='scale'`` is passed then it uses 1 / (n_features * X.var())
        as value of gamma. The current default of gamma, 'auto', will change
        to 'scale' in version 0.22. 'auto_deprecated', a deprecated version of
        'auto' is used as a default indicating that no explicit value of gamma
        was passed.

    coef0 : float, optional (default=0.0)
        Independent term in kernel function.
        It is only significant in 'poly' and 'sigmoid'.

    tol : float, optional (default=1e-3)
        Tolerance for stopping criterion.

    C : float, optional (default=1.0)
        Penalty parameter C of the error term.

    epsilon : float, optional (default=0.1)
         Epsilon in the epsilon-SVR model. It specifies the epsilon-tube
         within which no penalty is associated in the training loss function
         with points predicted within a distance epsilon from the actual
         value.

    shrinking : boolean, optional (default=True)
        Whether to use the shrinking heuristic.

    cache_size : float, optional
        Specify the size of the kernel cache (in MB).

    verbose : bool, default: False
        Enable verbose output. Note that this setting takes advantage of a
        per-process runtime setting in libsvm that, if enabled, may not work
        properly in a multithreaded context.

    max_iter : int, optional (default=-1)
        Hard limit on iterations within solver, or -1 for no limit.

    Attributes
    ----------
    support_ : array-like, shape = [n_SV]
        Indices of support vectors.

    support_vectors_ : array-like, shape = [nSV, n_features]
        Support vectors.

    dual_coef_ : array, shape = [1, n_SV]
        Coefficients of the support vector in the decision function.

    coef_ : array, shape = [1, n_features]
        Weights assigned to the features (coefficients in the primal
        problem). This is only available in the case of a linear kernel.

        `coef_` is readonly property derived from `dual_coef_` and
        `support_vectors_`.

    intercept_ : array, shape = [1]
        Constants in decision function.

    Examples
    --------
    >>> from sklearn.svm import SVR
    >>> import numpy as np
    >>> n_samples, n_features = 10, 5
    >>> rng = np.random.RandomState(0)
    >>> y = rng.randn(n_samples)
    >>> X = rng.randn(n_samples, n_features)
    >>> clf = SVR(gamma='scale', C=1.0, epsilon=0.2)
    >>> clf.fit(X, y) #doctest: +NORMALIZE_WHITESPACE
    SVR(C=1.0, cache_size=200, coef0=0.0, degree=3, epsilon=0.2, gamma='scale',
        kernel='rbf', max_iter=-1, shrinking=True, tol=0.001, verbose=False)

    See also
    --------
    NuSVR
        Support Vector Machine for regression implemented using libsvm
        using a parameter to control the number of support vectors.

    LinearSVR
        Scalable Linear Support Vector Machine for regression
        implemented using liblinear.

    Notes
    -----
    **References:**
    `LIBSVM: A Library for Support Vector Machines
    <http://www.csie.ntu.edu.tw/~cjlin/papers/libsvm.pdf>`__
    """

    _impl = 'epsilon_svr'

    def __init__(self, kernel='rbf', degree=3, gamma='auto_deprecated',
                 coef0=0.0, tol=1e-3, C=1.0, epsilon=0.1, shrinking=True,
                 cache_size=200, verbose=False, max_iter=-1):

        super().__init__(
            kernel=kernel, degree=degree, gamma=gamma,
            coef0=coef0, tol=tol, C=C, nu=0., epsilon=epsilon, verbose=verbose,
            shrinking=shrinking, probability=False, cache_size=cache_size,
            class_weight=None, max_iter=max_iter, random_state=None)


class NuSVR(BaseLibSVM, RegressorMixin):
    """Nu Support Vector Regression.

    Similar to NuSVC, for regression, uses a parameter nu to control
    the number of support vectors. However, unlike NuSVC, where nu
    replaces C, here nu replaces the parameter epsilon of epsilon-SVR.

    The implementation is based on libsvm.

    Read more in the :ref:`User Guide <svm_regression>`.

    Parameters
    ----------
    nu : float, optional
        An upper bound on the fraction of training errors and a lower bound of
        the fraction of support vectors. Should be in the interval (0, 1].  By
        default 0.5 will be taken.

    C : float, optional (default=1.0)
        Penalty parameter C of the error term.

    kernel : string, optional (default='rbf')
         Specifies the kernel type to be used in the algorithm.
         It must be one of 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed' or
         a callable.
         If none is given, 'rbf' will be used. If a callable is given it is
         used to precompute the kernel matrix.

    degree : int, optional (default=3)
        Degree of the polynomial kernel function ('poly').
        Ignored by all other kernels.

    gamma : float, optional (default='auto')
        Kernel coefficient for 'rbf', 'poly' and 'sigmoid'.

        Current default is 'auto' which uses 1 / n_features,
        if ``gamma='scale'`` is passed then it uses 1 / (n_features * X.var())
        as value of gamma. The current default of gamma, 'auto', will change
        to 'scale' in version 0.22. 'auto_deprecated', a deprecated version of
        'auto' is used as a default indicating that no explicit value of gamma
        was passed.

    coef0 : float, optional (default=0.0)
        Independent term in kernel function.
        It is only significant in 'poly' and 'sigmoid'.

    shrinking : boolean, optional (default=True)
        Whether to use the shrinking heuristic.

    tol : float, optional (default=1e-3)
        Tolerance for stopping criterion.

    cache_size : float, optional
        Specify the size of the kernel cache (in MB).

    verbose : bool, default: False
        Enable verbose output. Note that this setting takes advantage of a
        per-process runtime setting in libsvm that, if enabled, may not work
        properly in a multithreaded context.

    max_iter : int, optional (default=-1)
        Hard limit on iterations within solver, or -1 for no limit.

    Attributes
    ----------
    support_ : array-like, shape = [n_SV]
        Indices of support vectors.

    support_vectors_ : array-like, shape = [nSV, n_features]
        Support vectors.

    dual_coef_ : array, shape = [1, n_SV]
        Coefficients of the support vector in the decision function.

    coef_ : array, shape = [1, n_features]
        Weights assigned to the features (coefficients in the primal
        problem). This is only available in the case of a linear kernel.

        `coef_` is readonly property derived from `dual_coef_` and
        `support_vectors_`.

    intercept_ : array, shape = [1]
        Constants in decision function.

    Examples
    --------
    >>> from sklearn.svm import NuSVR
    >>> import numpy as np
    >>> n_samples, n_features = 10, 5
    >>> np.random.seed(0)
    >>> y = np.random.randn(n_samples)
    >>> X = np.random.randn(n_samples, n_features)
    >>> clf = NuSVR(gamma='scale', C=1.0, nu=0.1)
    >>> clf.fit(X, y)  #doctest: +NORMALIZE_WHITESPACE
    NuSVR(C=1.0, cache_size=200, coef0=0.0, degree=3, gamma='scale',
          kernel='rbf', max_iter=-1, nu=0.1, shrinking=True, tol=0.001,
          verbose=False)

    See also
    --------
    NuSVC
        Support Vector Machine for classification implemented with libsvm
        with a parameter to control the number of support vectors.

    SVR
        epsilon Support Vector Machine for regression implemented with libsvm.

    Notes
    -----
    **References:**
    `LIBSVM: A Library for Support Vector Machines
    <http://www.csie.ntu.edu.tw/~cjlin/papers/libsvm.pdf>`__
    """

    _impl = 'nu_svr'

    def __init__(self, nu=0.5, C=1.0, kernel='rbf', degree=3,
                 gamma='auto_deprecated', coef0=0.0, shrinking=True,
                 tol=1e-3, cache_size=200, verbose=False, max_iter=-1):

        super().__init__(
            kernel=kernel, degree=degree, gamma=gamma, coef0=coef0,
            tol=tol, C=C, nu=nu, epsilon=0., shrinking=shrinking,
            probability=False, cache_size=cache_size, class_weight=None,
            verbose=verbose, max_iter=max_iter, random_state=None)


class OneClassSVM(BaseLibSVM, OutlierMixin):
    """Unsupervised Outlier Detection.

    Estimate the support of a high-dimensional distribution.

    The implementation is based on libsvm.

    Read more in the :ref:`User Guide <outlier_detection>`.

    Parameters
    ----------
    kernel : string, optional (default='rbf')
         Specifies the kernel type to be used in the algorithm.
         It must be one of 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed' or
         a callable.
         If none is given, 'rbf' will be used. If a callable is given it is
         used to precompute the kernel matrix.

    degree : int, optional (default=3)
        Degree of the polynomial kernel function ('poly').
        Ignored by all other kernels.

    gamma : float, optional (default='auto')
        Kernel coefficient for 'rbf', 'poly' and 'sigmoid'.

        Current default is 'auto' which uses 1 / n_features,
        if ``gamma='scale'`` is passed then it uses 1 / (n_features * X.var())
        as value of gamma. The current default of gamma, 'auto', will change
        to 'scale' in version 0.22. 'auto_deprecated', a deprecated version of
        'auto' is used as a default indicating that no explicit value of gamma
        was passed.

    coef0 : float, optional (default=0.0)
        Independent term in kernel function.
        It is only significant in 'poly' and 'sigmoid'.

    tol : float, optional
        Tolerance for stopping criterion.

    nu : float, optional
        An upper bound on the fraction of training
        errors and a lower bound of the fraction of support
        vectors. Should be in the interval (0, 1]. By default 0.5
        will be taken.

    shrinking : boolean, optional
        Whether to use the shrinking heuristic.

    cache_size : float, optional
        Specify the size of the kernel cache (in MB).

    verbose : bool, default: False
        Enable verbose output. Note that this setting takes advantage of a
        per-process runtime setting in libsvm that, if enabled, may not work
        properly in a multithreaded context.

    max_iter : int, optional (default=-1)
        Hard limit on iterations within solver, or -1 for no limit.

    random_state : int, RandomState instance or None, optional (default=None)
        Ignored.

        .. deprecated:: 0.20
           ``random_state`` has been deprecated in 0.20 and will be removed in
           0.22.

    Attributes
    ----------
    support_ : array-like, shape = [n_SV]
        Indices of support vectors.

    support_vectors_ : array-like, shape = [nSV, n_features]
        Support vectors.

    dual_coef_ : array, shape = [1, n_SV]
        Coefficients of the support vectors in the decision function.

    coef_ : array, shape = [1, n_features]
        Weights assigned to the features (coefficients in the primal
        problem). This is only available in the case of a linear kernel.

        `coef_` is readonly property derived from `dual_coef_` and
        `support_vectors_`

    intercept_ : array, shape = [1,]
        Constant in the decision function.

    offset_ : float
        Offset used to define the decision function from the raw scores.
        We have the relation: decision_function = score_samples - `offset_`.
        The offset is the opposite of `intercept_` and is provided for
        consistency with other outlier detection algorithms.

    Examples
    --------
    >>> from sklearn.svm import OneClassSVM
    >>> X = [[0], [0.44], [0.45], [0.46], [1]]
    >>> clf = OneClassSVM(gamma='auto').fit(X)
    >>> clf.predict(X)
    array([-1,  1,  1,  1, -1])
    >>> clf.score_samples(X)  # doctest: +ELLIPSIS
    array([1.7798..., 2.0547..., 2.0556..., 2.0561..., 1.7332...])
    """

    _impl = 'one_class'

    def __init__(self, kernel='rbf', degree=3, gamma='auto_deprecated',
                 coef0=0.0, tol=1e-3, nu=0.5, shrinking=True, cache_size=200,
                 verbose=False, max_iter=-1, random_state=None):

        super().__init__(
            kernel, degree, gamma, coef0, tol, 0., nu, 0.,
            shrinking, False, cache_size, None, verbose, max_iter,
            random_state)

    def fit(self, X, y=None, sample_weight=None, **params):
        """
        Detects the soft boundary of the set of samples X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Set of samples, where n_samples is the number of samples and
            n_features is the number of features.

        sample_weight : array-like, shape (n_samples,)
            Per-sample weights. Rescale C per sample. Higher weights
            force the classifier to put more emphasis on these points.

        y : Ignored
            not used, present for API consistency by convention.

        Returns
        -------
        self : object

        Notes
        -----
        If X is not a C-ordered contiguous array it is copied.

        """

        if self.random_state is not None:
            warnings.warn("The random_state parameter is deprecated and will"
                          " be removed in version 0.22.", DeprecationWarning)

        super().fit(X, np.ones(_num_samples(X)),
                    sample_weight=sample_weight, **params)
        self.offset_ = -self._intercept_
        return self

    def decision_function(self, X):
        """Signed distance to the separating hyperplane.

        Signed distance is positive for an inlier and negative for an outlier.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        Returns
        -------
        dec : array-like, shape (n_samples,)
            Returns the decision function of the samples.
        """
        dec = self._decision_function(X).ravel()
        return dec

    def score_samples(self, X):
        """Raw scoring function of the samples.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        Returns
        -------
        score_samples : array-like, shape (n_samples,)
            Returns the (unshifted) scoring function of the samples.
        """
        return self.decision_function(X) + self.offset_

    def predict(self, X):
        """
        Perform classification on samples in X.

        For a one-class model, +1 or -1 is returned.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            For kernel="precomputed", the expected shape of X is
            [n_samples_test, n_samples_train]

        Returns
        -------
        y_pred : array, shape (n_samples,)
            Class labels for samples in X.
        """
        y = super().predict(X)
        return np.asarray(y, dtype=np.intp)
