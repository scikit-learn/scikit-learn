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

    def fit(self, X, y, coef_init=None, intercept_init=None,
            sample_weight=None):
        return self._fit(X, y, coef_init=None, intercept_init=None,
            sample_weight=None)

    def _fit(X, y, coef_init=None, intercept_init=None,
            sample_weight=None):

        n_samples, n_features = X.shape[0], X.shape[1]
        gradient_memory = np.zeros(n_samples)
        rng = check_random_state(self.random_state)

        for _ in range(n_iter * X.shape[0]):
            idx = int(rng.rand(1) * n_samples))
            sample = X[]




class SAGClassifier(BaseSAGClassifier, _LearntSelectorMixin):

