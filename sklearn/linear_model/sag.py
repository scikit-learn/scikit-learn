import numpy as np
import scipy.sparse as sp

from abc import ABCMeta

from .stochastic_gradient import BaseSGD
from .base import LinearClassifierMixin, LinearModel
from ..base import RegressorMixin
from sklearn.feature_selection.from_model import _LearntSelectorMixin
from ..utils import check_random_state
from ..externals import six
from .sgd_fast import Log, SquaredLoss
from ..externals.joblib import Parallel, delayed
from .sag_fast import fast_fit, fast_fit_sparse
from ..utils.seq_dataset import ArrayDataset, CSRDataset


class BaseSAG(six.with_metaclass(ABCMeta, BaseSGD)):
    def __init__(self, loss=None, penalty='l2', alpha=0.0001, l1_ratio=0.0,
                 fit_intercept=True, n_iter=5, shuffle=False, verbose=0,
                 epsilon=.1, n_jobs=1, random_state=None,
                 learning_rate="optimal", eta0=0.0, power_t=0.5,
                 class_weight=None, warm_start=False):

        self.gradient_memory = None

        super(BaseSAG, self).__init__(loss=loss, penalty=penalty,
                                      alpha=alpha, l1_ratio=l1_ratio,
                                      fit_intercept=fit_intercept,
                                      n_iter=n_iter, shuffle=shuffle,
                                      verbose=verbose,
                                      epsilon=epsilon,
                                      random_state=random_state,
                                      learning_rate=learning_rate,
                                      eta0=eta0, power_t=power_t,
                                      warm_start=warm_start,
                                      average=False)

    def partial_fit(self, X, y, sample_weight=None):
        raise ValueError("partial fit not supported for SAG")

    def _fit(self, X, y, coef_init=None, intercept_init=None,
             sample_weight=None):

        n_samples, n_features = X.shape[0], X.shape[1]
        sample_weight = np.ones(n_samples, dtype=np.float64, order='C')
        if sp.issparse(X):
            dataset = CSRDataset(X.data, X.indptr, X.indices,
                                 y, sample_weight, seed=self.random_state)
        else:
            dataset = ArrayDataset(X, y, sample_weight, seed=self.random_state)

        coef_ = np.zeros(n_features, dtype=np.float64, order='C')
        loss_function = self._get_loss_function(self.loss)

        # intercept_ = fast_fit(dataset, coef_, n_samples, n_features,
        #                       self.n_iter, loss_function,
        #                       self.eta0, self.alpha)
        intercept_ = fast_fit_sparse(dataset, coef_, n_samples, n_features,
                                     self.n_iter, loss_function,
                                     self.eta0, self.alpha)
        return coef_.reshape(1, -1), intercept_


class BaseSAGClassifier(six.with_metaclass(ABCMeta, BaseSAG,
                                           LinearClassifierMixin)):
    def __init__(self, penalty='l2', alpha=0.0001,
                 fit_intercept=True, n_iter=5, shuffle=False, verbose=0,
                 epsilon=.1, n_jobs=1, random_state=None,
                 learning_rate="optimal", eta0=0.0, power_t=0.5,
                 class_weight=None, warm_start=False):
        self.n_jobs = n_jobs
        self.loss_functions = {"log": (Log, )}
        super(BaseSAGClassifier, self).__init__(penalty=penalty,
                                                loss="log",
                                                alpha=alpha,
                                                fit_intercept=fit_intercept,
                                                n_iter=n_iter, shuffle=shuffle,
                                                verbose=verbose,
                                                epsilon=epsilon,
                                                random_state=random_state,
                                                learning_rate=learning_rate,
                                                eta0=eta0, power_t=power_t,
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
    def __init__(self, penalty='l2', alpha=0.0001,
                 fit_intercept=True, n_iter=5, shuffle=False, verbose=0,
                 epsilon=.1, n_jobs=1, random_state=None,
                 learning_rate="optimal", eta0=0.0, power_t=0.5,
                 class_weight=None, warm_start=False):

        super(SAGClassifier, self).__init__(penalty=penalty,
                                            alpha=alpha,
                                            fit_intercept=fit_intercept,
                                            n_iter=n_iter, shuffle=shuffle,
                                            verbose=verbose,
                                            n_jobs=n_jobs,
                                            epsilon=epsilon,
                                            random_state=random_state,
                                            learning_rate=learning_rate,
                                            eta0=eta0, power_t=power_t,
                                            warm_start=warm_start)

    def fit(self, X, y, coef_init=None, intercept_init=None,
            sample_weight=None):
        super(SAGClassifier, self)._fit(X, y, coef_init=None,
                                        intercept_init=None,
                                        sample_weight=None)


class BaseSAGRegressor(six.with_metaclass(ABCMeta, BaseSAG,
                                          LinearModel, RegressorMixin)):
    def __init__(self, penalty='l2', alpha=0.0001,
                 fit_intercept=True, n_iter=5, shuffle=False, verbose=0,
                 epsilon=.1, n_jobs=1, random_state=None,
                 learning_rate="optimal", eta0=0.001, power_t=0.5,
                 class_weight=None, warm_start=False):

        self.loss_functions = {"squared_loss": (SquaredLoss, )}
        super(BaseSAGRegressor, self).__init__(penalty=penalty,
                                               alpha=alpha,
                                               loss="squared_loss",
                                               fit_intercept=fit_intercept,
                                               n_iter=n_iter, shuffle=shuffle,
                                               verbose=verbose,
                                               epsilon=epsilon,
                                               random_state=random_state,
                                               learning_rate=learning_rate,
                                               eta0=eta0, power_t=power_t,
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
    a decreasing strength schedule (aka learning rate).

    The regularizer is a penalty added to the loss function that shrinks model
    parameters towards the zero vector using either the squared euclidean norm
    L2 or the absolute norm L1 or a combination of both (Elastic Net). If the
    parameter update crosses the 0.0 value because of the regularizer, the
    update is truncated to 0.0 to allow for learning sparse models and achieve
    online feature selection.

    This implementation works with data represented as dense numpy arrays of
    floating point values for the features.

    Parameters
    ----------
    loss : str, 'squared_loss', 'huber', 'epsilon_insensitive', \
                or 'squared_epsilon_insensitive'
        The loss function to be used. Defaults to 'squared_loss' which refers
        to the ordinary least squares fit. 'huber' modifies 'squared_loss' to
        focus less on getting outliers correct by switching from squared to
        linear loss past a distance of epsilon. 'epsilon_insensitive' ignores
        errors less than epsilon and is linear past that; this is the loss
        function used in SVR. 'squared_epsilon_insensitive' is the same but
        becomes squared loss past a tolerance of epsilon.

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

    epsilon: float
        Epsilon in the epsilon-insensitive loss functions; only if `loss` is
        'huber', 'epsilon_insensitive', or 'squared_epsilon_insensitive'.
        For 'huber', determines the threshold at which it becomes less
        important to get the prediction exactly right.
        For epsilon-insensitive, any differences between the current prediction
        and the correct label are ignored if they are less than this threshold.

    learning_rate : string, optional
        The learning rate:
        constant: eta = eta0
        optimal: eta = 1.0/(alpha * t)
        invscaling: eta = eta0 / pow(t, power_t) [default]

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
    SAGRegressor(alpha=0.0001, epsilon=0.1, eta0=0.01,
                 fit_intercept=True, l1_ratio=0.15, learning_rate='invscaling',
                 n_iter=5, power_t=0.25, random_state=None,
                 shuffle=False, verbose=0, warm_start=False)

    See also
    --------
    SGDRegressor, Ridge, ElasticNet, Lasso, SVR

    """
    def __init__(self, penalty='l2', alpha=0.0001,
                 fit_intercept=True, n_iter=5, shuffle=False, verbose=0,
                 epsilon=.1, n_jobs=1, random_state=None,
                 learning_rate="invscaling", eta0=0.001, power_t=0.5,
                 class_weight=None, warm_start=False):

        super(SAGRegressor, self).__init__(penalty=penalty,
                                           alpha=alpha,
                                           fit_intercept=fit_intercept,
                                           n_iter=n_iter, shuffle=shuffle,
                                           verbose=verbose,
                                           epsilon=epsilon,
                                           random_state=random_state,
                                           learning_rate=learning_rate,
                                           eta0=eta0, power_t=power_t,
                                           warm_start=warm_start)

    def fit(self, X, y, coef_init=None, intercept_init=None,
            sample_weight=None):
        super(SAGRegressor, self)._fit(X, y, coef_init=None,
                                       intercept_init=None,
                                       sample_weight=None)

