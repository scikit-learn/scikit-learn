import warnings
import numpy as np

from .base import _fit_liblinear, BaseSVC, BaseLibSVM
from ..base import BaseEstimator, RegressorMixin
from ..linear_model.base import LinearClassifierMixin, SparseCoefMixin, \
    LinearModel
from ..feature_selection.from_model import _LearntSelectorMixin
from ..utils import check_X_y


class LinearSVC(BaseEstimator, LinearClassifierMixin,
                _LearntSelectorMixin, SparseCoefMixin):
    """Linear Support Vector Classification.

    Similar to SVC with parameter kernel='linear', but implemented in terms of
    liblinear rather than libsvm, so it has more flexibility in the choice of
    penalties and loss functions and should scale better (to large numbers of
    samples).

    This class supports both dense and sparse input and the multiclass support
    is handled according to a one-vs-the-rest scheme.

    Parameters
    ----------
    C : float, optional (default=1.0)
        Penalty parameter C of the error term.

    loss : string, 'hinge' or 'squared_hinge' (default='squared_hinge')
        Specifies the loss function. 'hinge' is the standard SVM loss
        (used e.g. by the SVC class) while 'squared_hinge' is the
        square of the hinge loss.

    penalty : string, 'l1' or 'l2' (default='l2')
        Specifies the norm used in the penalization. The 'l2'
        penalty is the standard used in SVC. The 'l1' leads to `coef_`
        vectors that are sparse.

    dual : bool, (default=True)
        Select the algorithm to either solve the dual or primal
        optimization problem. Prefer dual=False when n_samples > n_features.

    tol : float, optional (default=1e-4)
        Tolerance for stopping criteria.

    multi_class: string, 'ovr' or 'crammer_singer' (default='ovr')
        Determines the multi-class strategy if `y` contains more than
        two classes.
        `ovr` trains n_classes one-vs-rest classifiers, while `crammer_singer`
        optimizes a joint objective over all classes.
        While `crammer_singer` is interesting from an theoretical perspective
        as it is consistent it is seldom used in practice and rarely leads to
        better accuracy and is more expensive to compute.
        If `crammer_singer` is chosen, the options loss, penalty and dual will
        be ignored.

    fit_intercept : boolean, optional (default=True)
        Whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    intercept_scaling : float, optional (default=1)
        When self.fit_intercept is True, instance vector x becomes
        [x, self.intercept_scaling],
        i.e. a "synthetic" feature with constant value equals to
        intercept_scaling is appended to the instance vector.
        The intercept becomes intercept_scaling * synthetic feature weight
        Note! the synthetic feature weight is subject to l1/l2 regularization
        as all other features.
        To lessen the effect of regularization on synthetic feature weight
        (and therefore on the intercept) intercept_scaling has to be increased

    class_weight : {dict, 'auto'}, optional
        Set the parameter C of class i to class_weight[i]*C for
        SVC. If not given, all classes are supposed to have
        weight one. The 'auto' mode uses the values of y to
        automatically adjust weights inversely proportional to
        class frequencies.

    verbose : int, (default=0)
        Enable verbose output. Note that this setting takes advantage of a
        per-process runtime setting in liblinear that, if enabled, may not work
        properly in a multithreaded context.

    random_state : int seed, RandomState instance, or None (default=None)
        The seed of the pseudo random number generator to use when
        shuffling the data.

    max_iter : int, (default=1000)
        The maximum number of iterations to be run.

    Attributes
    ----------
    coef_ : array, shape = [n_features] if n_classes == 2 \
            else [n_classes, n_features]
        Weights assigned to the features (coefficients in the primal
        problem). This is only available in the case of linear kernel.

        `coef_` is a readonly property derived from `raw_coef_` that \
        follows the internal memory layout of liblinear.

    intercept_ : array, shape = [1] if n_classes == 2 else [n_classes]
        Constants in decision function.

    Notes
    -----
    The underlying C implementation uses a random number generator to
    select features when fitting the model. It is thus not uncommon,
    to have slightly different results for the same input data. If
    that happens, try with a smaller tol parameter.

    The underlying implementation (liblinear) uses a sparse internal
    representation for the data that will incur a memory copy.

    Predict output may not match that of standalone liblinear in certain
    cases. See :ref:`differences from liblinear <liblinear_differences>`
    in the narrative documentation.

    **References:**
    `LIBLINEAR: A Library for Large Linear Classification
    <http://www.csie.ntu.edu.tw/~cjlin/liblinear/>`__

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
        # FIXME Remove l1/l2 support in 1.0 -----------------------------------
        loss_l = self.loss.lower()

        msg = ("loss='%s' has been deprecated in favor of "
               "loss='%s' as of 0.16. Backward compatibility"
               " for the loss='%s' will be removed in %s")

        # FIXME change loss_l --> self.loss after 0.18
        if loss_l in ('l1', 'l2'):
            old_loss = self.loss
            self.loss = {'l1': 'hinge', 'l2': 'squared_hinge'}.get(loss_l)
            warnings.warn(msg % (old_loss, self.loss, old_loss, '1.0'),
                          DeprecationWarning)
        # ---------------------------------------------------------------------

        if self.C < 0:
            raise ValueError("Penalty term must be positive; got (C=%r)"
                             % self.C)

        X, y = check_X_y(X, y, accept_sparse='csr',
                         dtype=np.float64, order="C")
        self.classes_ = np.unique(y)

        self.coef_, self.intercept_, self.n_iter_ = _fit_liblinear(
            X, y, self.C, self.fit_intercept, self.intercept_scaling,
            self.class_weight, self.penalty, self.dual, self.verbose,
            self.max_iter, self.tol, self.random_state, self.multi_class,
            self.loss
            )

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
    penalties and loss functions and should scale better (to large numbers of
    samples).

    This class supports both dense and sparse input.

    Parameters
    ----------
    C : float, optional (default=1.0)
        Penalty parameter C of the error term. The penalty is a squared
        l2 penalty. The bigger this parameter, the less regularization is used.

    loss : string, 'epsilon_insensitive' or 'squared_epsilon_insensitive' \
           (default='epsilon_insensitive')
        Specifies the loss function. 'l1' is the epsilon-insensitive loss
        (standard SVR) while 'l2' is the squared epsilon-insensitive loss.

    epsilon : float, optional (default=0.1)
        Epsilon parameter in the epsilon-insensitive loss function. Note
        that the value of this parameter depends on the scale of the target
        variable y. If unsure, set epsilon=0.

    dual : bool, (default=True)
        Select the algorithm to either solve the dual or primal
        optimization problem. Prefer dual=False when n_samples > n_features.

    tol : float, optional (default=1e-4)
        Tolerance for stopping criteria.

    fit_intercept : boolean, optional (default=True)
        Whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

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

    verbose : int, (default=0)
        Enable verbose output. Note that this setting takes advantage of a
        per-process runtime setting in liblinear that, if enabled, may not work
        properly in a multithreaded context.

    random_state : int seed, RandomState instance, or None (default=None)
        The seed of the pseudo random number generator to use when
        shuffling the data.

    max_iter : int, (default=1000)
        The maximum number of iterations to be run.

    Attributes
    ----------
    coef_ : array, shape = [n_features] if n_classes == 2 \
            else [n_classes, n_features]
        Weights assigned to the features (coefficients in the primal
        problem). This is only available in the case of linear kernel.

        `coef_` is a readonly property derived from `raw_coef_` that \
        follows the internal memory layout of liblinear.

    intercept_ : array, shape = [1] if n_classes == 2 else [n_classes]
        Constants in decision function.


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
        # FIXME Remove l1/l2 support in 1.0 -----------------------------------
        loss_l = self.loss.lower()

        msg = ("loss='%s' has been deprecated in favor of "
               "loss='%s' as of 0.16. Backward compatibility"
               " for the loss='%s' will be removed in %s")

        # FIXME change loss_l --> self.loss after 0.18
        if loss_l in ('l1', 'l2'):
            old_loss = self.loss
            self.loss = {'l1': 'epsilon_insensitive',
                         'l2': 'squared_epsilon_insensitive'
                         }.get(loss_l)
            warnings.warn(msg % (old_loss, self.loss, old_loss, '1.0'),
                          DeprecationWarning)
        # ---------------------------------------------------------------------

        if self.C < 0:
            raise ValueError("Penalty term must be positive; got (C=%r)"
                             % self.C)

        X, y = check_X_y(X, y, accept_sparse='csr',
                         dtype=np.float64, order="C")
        penalty = 'l2'  # SVR only accepts l2 penalty
        self.coef_, self.intercept_, self.n_iter_ = _fit_liblinear(
            X, y, self.C, self.fit_intercept, self.intercept_scaling,
            None, penalty, self.dual, self.verbose,
            self.max_iter, self.tol, self.random_state, loss=self.loss,
            epsilon=self.epsilon)
        self.coef_ = self.coef_.ravel()

        return self


class SVC(BaseSVC):
    """C-Support Vector Classification.

    The implementation is based on libsvm. The fit time complexity
    is more than quadratic with the number of samples which makes it hard
    to scale to dataset with more than a couple of 10000 samples.

    The multiclass support is handled according to a one-vs-one scheme.

    For details on the precise mathematical formulation of the provided
    kernel functions and how `gamma`, `coef0` and `degree` affect each
    other, see the corresponding section in the narrative documentation:
    :ref:`svm_kernels`.

    .. The narrative documentation is available at http://scikit-learn.org/

    Parameters
    ----------
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

    gamma : float, optional (default=0.0)
        Kernel coefficient for 'rbf', 'poly' and 'sigmoid'.
        If gamma is 0.0 then 1/n_features will be used instead.

    coef0 : float, optional (default=0.0)
        Independent term in kernel function.
        It is only significant in 'poly' and 'sigmoid'.

    probability: boolean, optional (default=False)
        Whether to enable probability estimates. This must be enabled prior
        to calling `fit`, and will slow down that method.

    shrinking: boolean, optional (default=True)
        Whether to use the shrinking heuristic.

    tol : float, optional (default=1e-3)
        Tolerance for stopping criterion.

    cache_size : float, optional
        Specify the size of the kernel cache (in MB).

    class_weight : {dict, 'auto'}, optional
        Set the parameter C of class i to class_weight[i]*C for
        SVC. If not given, all classes are supposed to have
        weight one. The 'auto' mode uses the values of y to
        automatically adjust weights inversely proportional to
        class frequencies.

    verbose : bool, default: False
        Enable verbose output. Note that this setting takes advantage of a
        per-process runtime setting in libsvm that, if enabled, may not work
        properly in a multithreaded context.

    max_iter : int, optional (default=-1)
        Hard limit on iterations within solver, or -1 for no limit.

    random_state : int seed, RandomState instance, or None (default)
        The seed of the pseudo random number generator to use when
        shuffling the data for probability estimation.

    Attributes
    ----------
    support_ : array-like, shape = [n_SV]
        Indices of support vectors.

    support_vectors_ : array-like, shape = [n_SV, n_features]
        Support vectors.

    n_support_ : array-like, dtype=int32, shape = [n_class]
        Number of support vectors for each class.

    dual_coef_ : array, shape = [n_class-1, n_SV]
        Coefficients of the support vector in the decision function. \
        For multiclass, coefficient for all 1-vs-1 classifiers. \
        The layout of the coefficients in the multiclass case is somewhat \
        non-trivial. See the section about multi-class classification in the \
        SVM section of the User Guide for details.

    coef_ : array, shape = [n_class-1, n_features]
        Weights assigned to the features (coefficients in the primal
        problem). This is only available in the case of linear kernel.

        `coef_` is a readonly property derived from `dual_coef_` and
        `support_vectors_`.

    intercept_ : array, shape = [n_class * (n_class-1) / 2]
        Constants in decision function.

    Examples
    --------
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    >>> y = np.array([1, 1, 2, 2])
    >>> from sklearn.svm import SVC
    >>> clf = SVC()
    >>> clf.fit(X, y) #doctest: +NORMALIZE_WHITESPACE
    SVC(C=1.0, cache_size=200, class_weight=None, coef0=0.0, degree=3,
        gamma=0.0, kernel='rbf', max_iter=-1, probability=False,
        random_state=None, shrinking=True, tol=0.001, verbose=False)
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

    """

    def __init__(self, C=1.0, kernel='rbf', degree=3, gamma=0.0,
                 coef0=0.0, shrinking=True, probability=False,
                 tol=1e-3, cache_size=200, class_weight=None,
                 verbose=False, max_iter=-1, random_state=None):

        super(SVC, self).__init__(
            'c_svc', kernel, degree, gamma, coef0, tol, C, 0., 0., shrinking,
            probability, cache_size, class_weight, verbose, max_iter,
            random_state)


class NuSVC(BaseSVC):
    """Nu-Support Vector Classification.

    Similar to SVC but uses a parameter to control the number of support
    vectors.

    The implementation is based on libsvm.

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
        Degree of kernel function
        is significant only in poly, rbf, sigmoid.

    gamma : float, optional (default=0.0)
        Kernel coefficient for rbf and poly, if gamma is 0.0 then 1/n_features
        will be taken.

    coef0 : float, optional (default=0.0)
        Independent term in kernel function. It is only significant
        in poly/sigmoid.

    probability: boolean, optional (default=False)
        Whether to enable probability estimates. This must be enabled prior
        to calling `fit`, and will slow down that method.

    shrinking: boolean, optional (default=True)
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

    random_state : int seed, RandomState instance, or None (default)
        The seed of the pseudo random number generator to use when
        shuffling the data for probability estimation.

    Attributes
    ----------
    support_ : array-like, shape = [n_SV]
        Indices of support vectors.

    support_vectors_ : array-like, shape = [n_SV, n_features]
        Support vectors.

    n_support_ : array-like, dtype=int32, shape = [n_class]
        Number of support vector for each class.

    dual_coef_ : array, shape = [n_class-1, n_SV]
        Coefficients of the support vector in the decision function. \
        For multiclass, coefficient for all 1-vs-1 classifiers. \
        The layout of the coefficients in the multiclass case is somewhat \
        non-trivial. See the section about multi-class classification in \
        the SVM section of the User Guide for details.

    coef_ : array, shape = [n_class-1, n_features]
        Weights assigned to the features (coefficients in the primal
        problem). This is only available in the case of linear kernel.

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
    >>> clf = NuSVC()
    >>> clf.fit(X, y) #doctest: +NORMALIZE_WHITESPACE
    NuSVC(cache_size=200, coef0=0.0, degree=3, gamma=0.0, kernel='rbf',
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
    """

    def __init__(self, nu=0.5, kernel='rbf', degree=3, gamma=0.0,
                 coef0=0.0, shrinking=True, probability=False,
                 tol=1e-3, cache_size=200, verbose=False, max_iter=-1,
                 random_state=None):

        super(NuSVC, self).__init__(
            'nu_svc', kernel, degree, gamma, coef0, tol, 0., nu, 0., shrinking,
            probability, cache_size, None, verbose, max_iter, random_state)


class SVR(BaseLibSVM, RegressorMixin):
    """Epsilon-Support Vector Regression.

    The free parameters in the model are C and epsilon.

    The implementation is based on libsvm.

    Parameters
    ----------
    C : float, optional (default=1.0)
        Penalty parameter C of the error term.

    epsilon : float, optional (default=0.1)
         Epsilon in the epsilon-SVR model. It specifies the epsilon-tube
         within which no penalty is associated in the training loss function
         with points predicted within a distance epsilon from the actual
         value.

    kernel : string, optional (default='rbf')
         Specifies the kernel type to be used in the algorithm.
         It must be one of 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed' or
         a callable.
         If none is given, 'rbf' will be used. If a callable is given it is
         used to precompute the kernel matrix.

    degree : int, optional (default=3)
        Degree of kernel function
        is significant only in poly, rbf, sigmoid.

    gamma : float, optional (default=0.0)
        Kernel coefficient for rbf and poly, if gamma is 0.0 then 1/n_features
        will be taken.

    coef0 : float, optional (default=0.0)
        independent term in kernel function. It is only significant
        in poly/sigmoid.

    shrinking: boolean, optional (default=True)
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
        problem). This is only available in the case of linear kernel.

        `coef_` is readonly property derived from `dual_coef_` and
        `support_vectors_`.

    intercept_ : array, shape = [1]
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
    >>> clf.fit(X, y) #doctest: +NORMALIZE_WHITESPACE
    SVR(C=1.0, cache_size=200, coef0=0.0, degree=3, epsilon=0.2, gamma=0.0,
        kernel='rbf', max_iter=-1, shrinking=True, tol=0.001, verbose=False)

    See also
    --------
    NuSVR
        Support Vector Machine for regression implemented using libsvm
        using a parameter to control the number of support vectors.

    LinearSVR
        Scalable Linear Support Vector Machine for regression
        implemented using liblinear.
    """
    def __init__(self, kernel='rbf', degree=3, gamma=0.0, coef0=0.0, tol=1e-3,
                 C=1.0, epsilon=0.1, shrinking=True, cache_size=200,
                 verbose=False, max_iter=-1):

        super(SVR, self).__init__(
            'epsilon_svr', kernel=kernel, degree=degree, gamma=gamma,
            coef0=coef0, tol=tol, C=C, nu=0., epsilon=epsilon, verbose=verbose,
            shrinking=shrinking, probability=False, cache_size=cache_size,
            class_weight=None, max_iter=max_iter, random_state=None)


class NuSVR(BaseLibSVM, RegressorMixin):
    """Nu Support Vector Regression.

    Similar to NuSVC, for regression, uses a parameter nu to control
    the number of support vectors. However, unlike NuSVC, where nu
    replaces C, here nu replaces the parameter epsilon of epsilon-SVR.

    The implementation is based on libsvm.

    Parameters
    ----------
    C : float, optional (default=1.0)
        Penalty parameter C of the error term.

    nu : float, optional
        An upper bound on the fraction of training errors and a lower bound of
        the fraction of support vectors. Should be in the interval (0, 1].  By
        default 0.5 will be taken.

    kernel : string, optional (default='rbf')
         Specifies the kernel type to be used in the algorithm.
         It must be one of 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed' or
         a callable.
         If none is given, 'rbf' will be used. If a callable is given it is
         used to precompute the kernel matrix.

    degree : int, optional (default=3)
        Degree of kernel function
        is significant only in poly, rbf, sigmoid.

    gamma : float, optional (default=0.0)
        Kernel coefficient for rbf and poly, if gamma is 0.0 then 1/n_features
        will be taken.

    coef0 : float, optional (default=0.0)
        Independent term in kernel function. It is only significant
        in poly/sigmoid.

    shrinking: boolean, optional (default=True)
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
        problem). This is only available in the case of linear kernel.

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
    >>> clf = NuSVR(C=1.0, nu=0.1)
    >>> clf.fit(X, y)  #doctest: +NORMALIZE_WHITESPACE
    NuSVR(C=1.0, cache_size=200, coef0=0.0, degree=3, gamma=0.0, kernel='rbf',
          max_iter=-1, nu=0.1, shrinking=True, tol=0.001, verbose=False)

    See also
    --------
    NuSVC
        Support Vector Machine for classification implemented with libsvm
        with a parameter to control the number of support vectors.

    SVR
        epsilon Support Vector Machine for regression implemented with libsvm.
    """

    def __init__(self, nu=0.5, C=1.0, kernel='rbf', degree=3,
                 gamma=0.0, coef0=0.0, shrinking=True, tol=1e-3,
                 cache_size=200, verbose=False, max_iter=-1):

        super(NuSVR, self).__init__(
            'nu_svr', kernel=kernel, degree=degree, gamma=gamma, coef0=coef0,
            tol=tol, C=C, nu=nu, epsilon=0., shrinking=shrinking,
            probability=False, cache_size=cache_size, class_weight=None,
            verbose=verbose, max_iter=max_iter, random_state=None)


class OneClassSVM(BaseLibSVM):
    """Unsupervised Outlier Detection.

    Estimate the support of a high-dimensional distribution.

    The implementation is based on libsvm.

    Parameters
    ----------
    kernel : string, optional (default='rbf')
         Specifies the kernel type to be used in the algorithm.
         It must be one of 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed' or
         a callable.
         If none is given, 'rbf' will be used. If a callable is given it is
         used to precompute the kernel matrix.

    nu : float, optional
        An upper bound on the fraction of training
        errors and a lower bound of the fraction of support
        vectors. Should be in the interval (0, 1]. By default 0.5
        will be taken.

    degree : int, optional (default=3)
        Degree of the polynomial kernel function ('poly').
        Ignored by all other kernels.

    gamma : float, optional (default=0.0)
        Kernel coefficient for 'rbf', 'poly' and 'sigmoid'.
        If gamma is 0.0 then 1/n_features will be used instead.

    coef0 : float, optional (default=0.0)
        Independent term in kernel function.
        It is only significant in 'poly' and 'sigmoid'.

    tol : float, optional
        Tolerance for stopping criterion.

    shrinking: boolean, optional
        Whether to use the shrinking heuristic.

    cache_size : float, optional
        Specify the size of the kernel cache (in MB).

    verbose : bool, default: False
        Enable verbose output. Note that this setting takes advantage of a
        per-process runtime setting in libsvm that, if enabled, may not work
        properly in a multithreaded context.

    max_iter : int, optional (default=-1)
        Hard limit on iterations within solver, or -1 for no limit.

    random_state : int seed, RandomState instance, or None (default)
        The seed of the pseudo random number generator to use when
        shuffling the data for probability estimation.

    Attributes
    ----------
    support_ : array-like, shape = [n_SV]
        Indices of support vectors.

    support_vectors_ : array-like, shape = [nSV, n_features]
        Support vectors.

    dual_coef_ : array, shape = [n_classes-1, n_SV]
        Coefficients of the support vectors in the decision function.

    coef_ : array, shape = [n_classes-1, n_features]
        Weights assigned to the features (coefficients in the primal
        problem). This is only available in the case of linear kernel.

        `coef_` is readonly property derived from `dual_coef_` and
        `support_vectors_`

    intercept_ : array, shape = [n_classes-1]
        Constants in decision function.

    """
    def __init__(self, kernel='rbf', degree=3, gamma=0.0, coef0=0.0, tol=1e-3,
                 nu=0.5, shrinking=True, cache_size=200, verbose=False,
                 max_iter=-1, random_state=None):

        super(OneClassSVM, self).__init__(
            'one_class', kernel, degree, gamma, coef0, tol, 0., nu, 0.,
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

        Returns
        -------
        self : object
            Returns self.

        Notes
        -----
        If X is not a C-ordered contiguous array it is copied.

        """
        super(OneClassSVM, self).fit(X, [], sample_weight=sample_weight,
                                     **params)
        return self
