import numpy as np
from .stochastic_gradient import BaseSGDClassifier
from sklearn.feature_selection.from_model \
    import _LearntSelectorMixin
from ..utils import check_random_state

DEFAULT_EPSILON = 0.1


class BaseSAGClassifier(BaseSGDClassifier):
    def __init__(self, loss="hinge", penalty='l2', alpha=0.0001, l1_ratio=0.15,
                 fit_intercept=True, n_iter=5, shuffle=False, verbose=0,
                 epsilon=DEFAULT_EPSILON, n_jobs=1, random_state=None,
                 learning_rate="optimal", eta0=0.0, power_t=0.5,
                 class_weight=None, warm_start=False, average=False,
                 sag=False):

        self.gradient_memory = None
        super(BaseSGDClassifier, self).__init__(loss=loss, penalty=penalty,
                                                alpha=alpha, l1_ratio=l1_ratio,
                                                fit_intercept=fit_intercept,
                                                n_iter=n_iter, shuffle=shuffle,
                                                verbose=verbose,
                                                epsilon=epsilon,
                                                random_state=random_state,
                                                learning_rate=learning_rate,
                                                eta0=eta0, power_t=power_t,
                                                warm_start=warm_start,
                                                average=average, sag=sag)

    def partial_fit(self, X, y, sample_weight=None):
        raise ValueError("partial fit not supported for SAG")

    def _fit(self, X, y, coef_init=None, intercept_init=None,
             sample_weight=None):

        n_samples, n_features = X.shape[0], X.shape[1]
        gradient_memory = np.zeros((n_samples, n_features))
        self.coef_ = np.zeros(n_features)
        self.intercept_ = 0.0
        self.classes_ = np.unique(y)
        sum_gradient = np.zeros(n_features)
        rng = check_random_state(self.random_state)
        seen = set()
        loss_function = self._get_loss_function(self.loss)

        for i in range(self.n_iter * n_samples):
            idx = int(rng.rand(1) * n_samples)
            sample = X[idx]
            seen.add(idx)

            p = np.dot(sample, self.coef_) + self.intercept_
            gradient = loss_function.dloss(p, y[idx])
            eta = self.eta0

            update = X[idx] * gradient + self.alpha * self.coef_
            sum_gradient += update - gradient_memory[idx]
            gradient_memory[idx] = update

            self.coef_ -= eta * sum_gradient / len(seen)
            self.intercept_ -= eta * gradient

        self.coef_ = self.coef_.reshape(1, -1)


class SAGClassifier(BaseSAGClassifier, _LearntSelectorMixin):
    def __init__(self, loss="hinge", penalty='l2', alpha=0.0001, l1_ratio=0.15,
                 fit_intercept=True, n_iter=5, shuffle=False, verbose=0,
                 epsilon=DEFAULT_EPSILON, n_jobs=1, random_state=None,
                 learning_rate="optimal", eta0=0.0, power_t=0.5,
                 class_weight=None, warm_start=False, average=False,
                 sag=False):

        self.gradient_memory = None
        super(BaseSAGClassifier, self).__init__(loss=loss, penalty=penalty,
                                                alpha=alpha, l1_ratio=l1_ratio,
                                                fit_intercept=fit_intercept,
                                                n_iter=n_iter, shuffle=shuffle,
                                                verbose=verbose,
                                                epsilon=epsilon,
                                                random_state=random_state,
                                                learning_rate=learning_rate,
                                                eta0=eta0, power_t=power_t,
                                                warm_start=warm_start,
                                                average=average, sag=sag)

    def fit(self, X, y, coef_init=None, intercept_init=None,
            sample_weight=None):
        return self._fit(X, y, coef_init=None, intercept_init=None,
                         sample_weight=None)

