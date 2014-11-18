import numpy as np
import scipy.sparse as sp

from abc import ABCMeta

from .base import LinearClassifierMixin, LinearModel
from ..base import RegressorMixin
from sklearn.feature_selection.from_model import _LearntSelectorMixin
from ..utils import (check_array, check_random_state, check_X_y)
from ..externals import six
from .sgd_fast import Log, SquaredLoss
from ..externals.joblib import Parallel, delayed
from .sag_fast import fast_fit_sparse, get_auto_eta
from ..utils.seq_dataset import ArrayDataset, CSRDataset


class BaseSAG(six.with_metaclass(ABCMeta)):
    def __init__(self, alpha=0.0001,
                 fit_intercept=True, n_iter=5, verbose=0,
                 random_state=None, eta0=0.0,
                 class_weight=None, warm_start=False):

        self.gradient_memory = None
        self.alpha = alpha
        self.fit_intercept = fit_intercept
        self.n_iter = n_iter
        self.verbose = verbose
        self.eta0 = eta0
        self.random_state = random_state
        self.class_weight = class_weight
        self.warm_start = warm_start

    def partial_fit(self, X, y, sample_weight=None):
        raise ValueError("partial fit not supported for SAG")

    def _fit(self, X, y, coef_init=None, intercept_init=None,
             sample_weight=None):
        X, y = check_X_y(X, y, "csr", copy=False, order='C', dtype=np.float64)
        n_samples, n_features = X.shape[0], X.shape[1]
        sample_weight = np.ones(n_samples, dtype=np.float64, order='C')
        if sp.issparse(X):
            dataset = CSRDataset(X.data, X.indptr, X.indices,
                                 y, sample_weight, seed=self.random_state)
        else:
            dataset = ArrayDataset(X, y, sample_weight, seed=self.random_state)

        if self.eta0 == 'auto':
            step_size = get_auto_eta(dataset, self.alpha, n_samples)
        else:
            step_size = self.eta0
        coef_ = np.zeros(n_features, dtype=np.float64, order='C')

        # intercept_ = fast_fit(dataset, coef_, n_samples, n_features,
        #                       self.n_iter, loss_function,
        #                       self.eta0, self.alpha)
        intercept_ = fast_fit_sparse(dataset, coef_, n_samples, n_features,
                                     self.n_iter, self.loss_function,
                                     step_size, self.alpha)
        return coef_.reshape(1, -1), intercept_


class BaseSAGClassifier(six.with_metaclass(ABCMeta, BaseSAG,
                                           LinearClassifierMixin)):
    def __init__(self, alpha=0.0001,
                 fit_intercept=True, n_iter=5, verbose=0,
                 n_jobs=1, random_state=None,
                 eta0=0.0, class_weight=None, warm_start=False):
        self.n_jobs = n_jobs
        self.loss_function = Log()
        super(BaseSAGClassifier, self).__init__(alpha=alpha,
                                                fit_intercept=fit_intercept,
                                                n_iter=n_iter,
                                                verbose=verbose,
                                                random_state=random_state,
                                                eta0=eta0,
                                                warm_start=warm_start)

    def _fit(self, X, y, coef_init=None, intercept_init=None,
             sample_weight=None):
        n_samples, n_features = X.shape[0], X.shape[1]
        self.classes_ = np.unique(y)

        if len(self.classes_) == 2:
            coef, intercept = self._fit_target_class(X, y, self.classes_[1])
        else:
            coef = []
            intercept = []

            results = Parallel(n_jobs=self.n_jobs,
                               backend="threading",
                               verbose=self.verbose)(
                delayed(self._fit_target_class)(X, y, cl)
                for cl in self.classes_)

            for coef_cl, intercept_cl in results:
                coef.append(coef_cl)
                intercept.append(intercept_cl)

            coef = np.vstack(coef)
            intercept = np.array(intercept)

        self.coef_ = coef
        self.intercept_ = intercept

    def _fit_target_class(self, X, y, target_class):
        n_samples, n_features = X.shape[0], X.shape[1]

        y_encoded = np.ones(n_samples)
        y_encoded[y != target_class] = -1.0

        return super(BaseSAGClassifier, self)._fit(X, y_encoded,
                                                   coef_init=None,
                                                   intercept_init=None,
                                                   sample_weight=None)


class SAGClassifier(BaseSAGClassifier, _LearntSelectorMixin):
    def __init__(self, alpha=0.0001, fit_intercept=True, n_iter=5,
                 verbose=0, n_jobs=1, random_state=None,
                 eta0=0.0, class_weight=None, warm_start=False):

        super(SAGClassifier, self).__init__(alpha=alpha,
                                            fit_intercept=fit_intercept,
                                            n_iter=n_iter,
                                            verbose=verbose,
                                            n_jobs=n_jobs,
                                            random_state=random_state,
                                            eta0=eta0,
                                            warm_start=warm_start)

    def fit(self, X, y, coef_init=None, intercept_init=None,
            sample_weight=None):
        super(SAGClassifier, self)._fit(X, y, coef_init=None,
                                        intercept_init=None,
                                        sample_weight=None)


class BaseSAGRegressor(six.with_metaclass(ABCMeta, BaseSAG,
                                          LinearModel, RegressorMixin)):
    def __init__(self, alpha=0.0001,
                 fit_intercept=True, n_iter=5, verbose=0,
                 n_jobs=1, random_state=None,
                 eta0=0.001, class_weight=None, warm_start=False):

        self.loss_function = SquaredLoss()
        super(BaseSAGRegressor, self).__init__(alpha=alpha,
                                               fit_intercept=fit_intercept,
                                               n_iter=n_iter,
                                               verbose=verbose,
                                               random_state=random_state,
                                               eta0=eta0,
                                               warm_start=warm_start)

    def _fit(self, X, y, coef_init=None, intercept_init=None,
             sample_weight=None):
        self.coef_, self.intercept_ = \
            super(BaseSAGRegressor, self)._fit(X, y,
                                               coef_init=None,
                                               intercept_init=None,
                                               sample_weight=None)


class SAGRegressor(BaseSAGRegressor, _LearntSelectorMixin):
    """Linear model fitted by minimizing a regularized empirical loss with SAG

    SAG stands for Stochastic Average Gradient: the gradient of the loss is
    estimated each sample at a time and the model is updated along the way with
    a constant learning rate.

    The regularizer is a penalty added to the loss function that shrinks model
    parameters towards the zero vector using the squared euclidean norm
    L2.

    This implementation works with data represented as dense or sparse numpy
    arrays of floating point values for the features.

    Parameters
    ----------
    alpha : float
        Constant that multiplies the regularization term. Defaults to 0.0001

    fit_intercept: bool
        Whether the intercept should be estimated or not. If False, the
        data is assumed to be already centered. Defaults to True.

    n_iter : int, optional
        The number of passes over the training data (aka epochs). The number
        of iterations is set to 1 if using partial_fit.
        Defaults to 5.

    verbose: integer, optional
        The verbosity level.

    eta0 : double, optional
        The initial learning rate [default 0.01].

    power_t : double, optional
        The exponent for inverse scaling learning rate [default 0.25].

    warm_start : bool, optional
        When set to True, reuse the solution of the previous call to fit as
        initialization, otherwise, just erase the previous solution.

    Attributes
    ----------
    coef_ : array, shape (n_features,)
        Weights asigned to the features.

    intercept_ : array, shape (1,)
        The intercept term.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn import linear_model
    >>> n_samples, n_features = 10, 5
    >>> np.random.seed(0)
    >>> y = np.random.randn(n_samples)
    >>> X = np.random.randn(n_samples, n_features)
    >>> clf = linear_model.SAGRegressor()
    >>> clf.fit(X, y)
    ... #doctest: +NORMALIZE_WHITESPACE
    SAGRegressor(alpha=0.0001, eta0=0.01,
                 fit_intercept=True, n_iter=5, random_state=None,
                 verbose=0, warm_start=False)

    See also
    --------
    SGDRegressor, Ridge, ElasticNet, Lasso, SVR

    """
    def __init__(self, alpha=0.0001, fit_intercept=True, n_iter=5, verbose=0,
                 n_jobs=1, random_state=None,
                 eta0=0.001, class_weight=None, warm_start=False):

        super(SAGRegressor, self).__init__(alpha=alpha,
                                           fit_intercept=fit_intercept,
                                           n_iter=n_iter,
                                           verbose=verbose,
                                           random_state=random_state,
                                           eta0=eta0, warm_start=warm_start)

    def fit(self, X, y, coef_init=None, intercept_init=None,
            sample_weight=None):
        super(SAGRegressor, self)._fit(X, y, coef_init=None,
                                       intercept_init=None,
                                       sample_weight=None)

