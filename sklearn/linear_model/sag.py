import numpy as np
from abc import ABCMeta
from .stochastic_gradient import BaseSGD
from .base import LinearClassifierMixin
from sklearn.feature_selection.from_model \
    import _LearntSelectorMixin
from ..utils import check_random_state
from ..externals import six
from .sgd_fast import Log
from ..externals.joblib import Parallel, delayed


class BaseSAG(six.with_metaclass(ABCMeta, BaseSGD)):
    def __init__(self, loss="log", penalty='l2', alpha=0.0001, l1_ratio=0.0,
                 fit_intercept=True, n_iter=5, shuffle=False, verbose=0,
                 epsilon=.1, n_jobs=1, random_state=None,
                 learning_rate="optimal", eta0=0.0, power_t=0.5,
                 class_weight=None, warm_start=False, average=False):

        self.gradient_memory = None
        self.loss_functions = {"log": (Log, )}

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
                                      average=average)

    def _fit(self, X, y, coef_init=None, intercept_init=None,
             sample_weight=None):

        n_samples, n_features = X.shape[0], X.shape[1]
        gradient_memory = np.zeros((n_samples, n_features))
        coef_ = np.zeros(n_features)
        intercept_ = 0.0
        sum_gradient = np.zeros(n_features)
        rng = check_random_state(self.random_state)
        seen = set()
        loss_function = self._get_loss_function(self.loss)
        eta = self.eta0

        for i in range(self.n_iter * n_samples):
            idx = int(rng.rand(1) * n_samples)
            sample = X[idx]
            seen.add(idx)

            p = np.dot(sample, coef_) + intercept_
            gradient = loss_function.dloss(p, y[idx])
            # eta = eta

            update = X[idx] * gradient + self.alpha * coef_
            sum_gradient += update - gradient_memory[idx]
            gradient_memory[idx] = update

            coef_ -= eta * sum_gradient / len(seen)
            intercept_ -= eta * gradient

        return coef_.reshape(1, -1), intercept_


class BaseSAGClassifier(six.with_metaclass(ABCMeta, BaseSAG,
                                           LinearClassifierMixin)):
    def __init__(self, penalty='l2', alpha=0.0001,
                 fit_intercept=True, n_iter=5, shuffle=False, verbose=0,
                 epsilon=.1, n_jobs=1, random_state=None,
                 learning_rate="optimal", eta0=0.0, power_t=0.5,
                 class_weight=None, warm_start=False, average=False):
        self.n_jobs = n_jobs
        super(BaseSAGClassifier, self).__init__(penalty=penalty,
                                                alpha=alpha,
                                                fit_intercept=fit_intercept,
                                                n_iter=n_iter, shuffle=shuffle,
                                                verbose=verbose,
                                                epsilon=epsilon,
                                                random_state=random_state,
                                                learning_rate=learning_rate,
                                                eta0=eta0, power_t=power_t,
                                                warm_start=warm_start,
                                                average=average)

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
                 class_weight=None, warm_start=False, average=False):

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
                                            warm_start=warm_start,
                                            average=average)

    def fit(self, X, y, coef_init=None, intercept_init=None,
            sample_weight=None):
        super(SAGClassifier, self)._fit(X, y, coef_init=None,
                                        intercept_init=None,
                                        sample_weight=None)

    def partial_fit(self, X, y, sample_weight=None):
        raise ValueError("partial fit not supported for SAG")

