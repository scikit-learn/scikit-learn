import numpy as np
from abc import ABCMeta
from .stochastic_gradient import BaseSGD
from .base import LinearClassifierMixin, LinearModel
from ..base import RegressorMixin
from sklearn.feature_selection.from_model import _LearntSelectorMixin
from ..utils import check_random_state
from ..externals import six
from .sgd_fast import Log, SquaredLoss
from ..externals.joblib import Parallel, delayed
from .sag_fast import fast_fit
from ..utils.seq_dataset import ArrayDataset


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
        dataset = ArrayDataset(X, y, sample_weight, seed=self.random_state)
        coef_ = np.zeros(n_features, dtype=np.float64, order='C')
        loss_function = self._get_loss_function(self.loss)

        intercept_ = fast_fit(dataset, coef_, n_samples, n_features,
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

