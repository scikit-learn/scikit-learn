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

        self.coef_ = None
        self.intercept_ = None

        self.num_seen_ = None
        self.seen_ = None
        self.sum_gradient_ = None
        self.gradient_memory_ = None

    def partial_fit(self, X, y, sample_weight=None):
        raise ValueError("partial fit not supported for SAG")

    def _fit(self, X, y, coef_init=None, intercept_init=None,
             sample_weight=None, sum_gradient_init=None,
             gradient_memory_init=None, seen_init=None, num_seen_init=None):
        X, y = check_X_y(X, y, "csr", copy=False, order='C', dtype=np.float64)
        n_samples, n_features = X.shape[0], X.shape[1]

        # initialize all parameters if there is no init
        if sample_weight is None:
            sample_weight = np.ones(n_samples, dtype=np.float64, order='C')

        if intercept_init is None:
            intercept_init = 0.0

        if coef_init is None:
            coef_init = np.zeros(n_features, dtype=np.float64, order='C')

        if sum_gradient_init is None:
            sum_gradient_init = np.zeros(n_features, dtype=np.float64,
                                         order='C')

        if gradient_memory_init is None:
            gradient_memory_init = np.zeros(n_samples, dtype=np.float64,
                                            order='C')

        if seen_init is None:
            seen_init = np.zeros(n_samples, dtype=np.int32, order='C')

        if num_seen_init is None:
            num_seen_init = 0

        # check which type of Sequential Dataset is needed
        if sp.issparse(X):
            dataset = CSRDataset(X.data, X.indptr, X.indices,
                                 y, sample_weight, seed=self.random_state)
        else:
            dataset = ArrayDataset(X, y, sample_weight, seed=self.random_state)

        # set the eta0 if needed, 'auto' is 1 / 4L where L is the max sum of
        # squares for over all samples
        if self.eta0 == 'auto':
            step_size = get_auto_eta(dataset, self.alpha, n_samples)
        else:
            step_size = self.eta0

        # intercept_ = fast_fit(dataset, coef_, n_samples, n_features,
        #                       self.n_iter, loss_function,
        #                       self.eta0, self.alpha)
        intercept_, num_seen = fast_fit_sparse(dataset, coef_init.ravel(),
                                               intercept_init, n_samples,
                                               n_features, self.n_iter,
                                               self.loss_function,
                                               step_size, self.alpha,
                                               sum_gradient_init.ravel(),
                                               gradient_memory_init.ravel(),
                                               seen_init.ravel(),
                                               num_seen_init)
        return (coef_init.reshape(1, -1), intercept_,
                sum_gradient_init.reshape(1, -1),
                gradient_memory_init.reshape(1, -1),
                seen_init.reshape(1, -1),
                num_seen)


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
             sample_weight=None, sum_gradient_init=None,
             gradient_memory_init=None, seen_init=None,
             num_seen_init=None):
        n_samples, n_features = X.shape[0], X.shape[1]
        self.classes_ = np.unique(y)

        if len(self.classes_) == 2:
            if self.warm_start:
                coef_init = self.coef_
                intercept_init = self.intercept_
                sum_gradient_init = self.sum_gradient_
                gradient_memory_init = self.gradient_memory_
                seen_init = self.seen_
                num_seen_init = self.num_seen_

            coef, intercept, sum_gradient, gradient_memory, seen, num_seen = \
                self._fit_target_class(X, y, self.classes_[1],
                                       coef_init, intercept_init,
                                       sample_weight, sum_gradient_init,
                                       gradient_memory_init,
                                       seen_init, num_seen_init)
        else:
            coef = []
            intercept = []

            # TODO do something nice here...
            sum_gradient = None
            gradient_memory = None
            seen = None
            num_seen = None

            results = Parallel(n_jobs=self.n_jobs,
                               backend="threading",
                               verbose=self.verbose)(
                delayed(self._fit_target_class)(X, y, cl)
                for cl in self.classes_)

            # TODO: do something with the warm start params
            for coef_cl, intercept_cl, _, _, _, _ in results:
                coef.append(coef_cl)
                intercept.append(intercept_cl)

            coef = np.vstack(coef)
            intercept = np.array(intercept)

        self.coef_ = coef
        self.intercept_ = intercept
        self.sum_gradient_ = sum_gradient
        self.gradient_memory_ = gradient_memory
        self.seen_ = seen
        self.num_seen_ = num_seen

    def _fit_target_class(self, X, y, target_class, coef_init=None,
                          intercept_init=None, sample_weight=None,
                          sum_gradient_init=None, gradient_memory_init=None,
                          seen_init=None, num_seen_init=None):

        n_samples, n_features = X.shape[0], X.shape[1]

        y_encoded = np.ones(n_samples)
        y_encoded[y != target_class] = -1.0

        return super(BaseSAGClassifier, self)._fit(X, y_encoded,
                                                   coef_init, intercept_init,
                                                   sample_weight,
                                                   sum_gradient_init,
                                                   gradient_memory_init,
                                                   seen_init, num_seen_init)


class SAGClassifier(BaseSAGClassifier, _LearntSelectorMixin):
    def __init__(self, alpha=0.0001, fit_intercept=True, n_iter=5,
                 verbose=0, n_jobs=1, random_state=1.0,
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
            sample_weight=None, sum_gradient_init=None,
            gradient_memory_init=None, seen_init=None,
            num_seen_init=None):
        super(SAGClassifier, self)._fit(X, y, coef_init, intercept_init,
                                        sample_weight, sum_gradient_init,
                                        gradient_memory_init, seen_init,
                                        num_seen_init)


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
             sample_weight=None, sum_gradient_init=None,
             gradient_memory_init=None, seen_init=None,
             num_seen_init=None):
        if self.warm_start:
            coef_init = self.coef_
            intercept_init = self.intercept_
            sum_gradient_init = self.sum_gradient_
            gradient_memory_init = self.gradient_memory_
            seen_init = self.seen_
            num_seen_init = self.num_seen_

        (self.coef_, self.intercept_, self.sum_gradient_,
         self.gradient_memory_, self.seen_, self.num_seen_) = \
            super(BaseSAGRegressor, self)._fit(X, y, coef_init,
                                               intercept_init,
                                               sample_weight,
                                               sum_gradient_init,
                                               gradient_memory_init,
                                               seen_init, num_seen_init)


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
                 n_jobs=1, random_state=1.0,
                 eta0=0.001, class_weight=None, warm_start=False):

        super(SAGRegressor, self).__init__(alpha=alpha,
                                           fit_intercept=fit_intercept,
                                           n_iter=n_iter,
                                           verbose=verbose,
                                           random_state=random_state,
                                           eta0=eta0, warm_start=warm_start)

    def fit(self, X, y, coef_init=None, intercept_init=None,
            sample_weight=None, sum_gradient_init=None,
            gradient_memory_init=None, seen_init=None,
            num_seen_init=None):
        super(SAGRegressor, self)._fit(X, y, coef_init,
                                       intercept_init,
                                       sample_weight,
                                       sum_gradient_init,
                                       gradient_memory_init,
                                       seen_init, num_seen_init)

