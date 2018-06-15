# Authors: Rob Zinkov, Mathieu Blondel
# License: BSD 3 clause

from .stochastic_gradient import BaseSGDClassifier
from .stochastic_gradient import BaseSGDRegressor
from .stochastic_gradient import DEFAULT_EPSILON


class PassiveAggressiveClassifier(BaseSGDClassifier):
    """Passive Aggressive Classifier

    Read more in the :ref:`User Guide <passive_aggressive>`.

    Parameters
    ----------

    C : float
        Maximum step size (regularization). Defaults to 1.0.

    fit_intercept : bool, default=False
        Whether the intercept should be estimated or not. If False, the
        data is assumed to be already centered.

    max_iter : int, optional
        The maximum number of passes over the training data (aka epochs).
        It only impacts the behavior in the ``fit`` method, and not the
        `partial_fit`.
        Defaults to 5. Defaults to 1000 from 0.21, or if tol is not None.

        .. versionadded:: 0.19

    tol : float or None, optional
        The stopping criterion. If it is not None, the iterations will stop
        when (loss > previous_loss - tol). Defaults to None.
        Defaults to 1e-3 from 0.21.

        .. versionadded:: 0.19

    shuffle : bool, default=True
        Whether or not the training data should be shuffled after each epoch.

    verbose : integer, optional
        The verbosity level

    loss : string, optional
        The loss function to be used:
        hinge: equivalent to PA-I in the reference paper.
        squared_hinge: equivalent to PA-II in the reference paper.

    n_jobs : integer, optional
        The number of CPUs to use to do the OVA (One Versus All, for
        multi-class problems) computation. -1 means 'all CPUs'. Defaults
        to 1.

    random_state : int, RandomState instance or None, optional, default=None
        The seed of the pseudo random number generator to use when shuffling
        the data.  If int, random_state is the seed used by the random number
        generator; If RandomState instance, random_state is the random number
        generator; If None, the random number generator is the RandomState
        instance used by `np.random`.

    warm_start : bool, optional
        When set to True, reuse the solution of the previous call to fit as
        initialization, otherwise, just erase the previous solution.
        See :term:`the Glossary <warm_start>`.

        Repeatedly calling fit or partial_fit when warm_start is True can
        result in a different solution than when calling fit a single time
        because of the way the data is shuffled.

    class_weight : dict, {class_label: weight} or "balanced" or None, optional
        Preset for the class_weight fit parameter.

        Weights associated with classes. If not given, all classes
        are supposed to have weight one.

        The "balanced" mode uses the values of y to automatically adjust
        weights inversely proportional to class frequencies in the input data
        as ``n_samples / (n_classes * np.bincount(y))``

        .. versionadded:: 0.17
           parameter *class_weight* to automatically weight samples.

    average : bool or int, optional
        When set to True, computes the averaged SGD weights and stores the
        result in the ``coef_`` attribute. If set to an int greater than 1,
        averaging will begin once the total number of samples seen reaches
        average. So average=10 will begin averaging after seeing 10 samples.

        .. versionadded:: 0.19
           parameter *average* to use weights averaging in SGD

    n_iter : int, optional
        The number of passes over the training data (aka epochs).
        Defaults to None. Deprecated, will be removed in 0.21.

        .. versionchanged:: 0.19
            Deprecated

    Attributes
    ----------
    coef_ : array, shape = [1, n_features] if n_classes == 2 else [n_classes,\
            n_features]
        Weights assigned to the features.

    intercept_ : array, shape = [1] if n_classes == 2 else [n_classes]
        Constants in decision function.

    n_iter_ : int
        The actual number of iterations to reach the stopping criterion.
        For multiclass fits, it is the maximum over every binary fit.

    Examples
    --------
    >>> from sklearn.linear_model import PassiveAggressiveClassifier
    >>> from sklearn.datasets import make_classification
    >>>
    >>> X, y = make_classification(n_features=4, random_state=0)
    >>> clf = PassiveAggressiveClassifier(random_state=0)
    >>> clf.fit(X, y)
    PassiveAggressiveClassifier(C=1.0, average=False, class_weight=None,
                  fit_intercept=True, loss='hinge', max_iter=None, n_iter=None,
                  n_jobs=1, random_state=0, shuffle=True, tol=None, verbose=0,
                  warm_start=False)
    >>> print(clf.coef_)
    [[0.49324685 1.0552176  1.49519589 1.33798314]]
    >>> print(clf.intercept_)
    [2.18438388]
    >>> print(clf.predict([[0, 0, 0, 0]]))
    [1]

    See also
    --------

    SGDClassifier
    Perceptron

    References
    ----------
    Online Passive-Aggressive Algorithms
    <http://jmlr.csail.mit.edu/papers/volume7/crammer06a/crammer06a.pdf>
    K. Crammer, O. Dekel, J. Keshat, S. Shalev-Shwartz, Y. Singer - JMLR (2006)

    """
    def __init__(self, C=1.0, fit_intercept=True, max_iter=None, tol=None,
                 shuffle=True, verbose=0, loss="hinge", n_jobs=1,
                 random_state=None, warm_start=False, class_weight=None,
                 average=False, n_iter=None, iters=1):
        super(PassiveAggressiveClassifier, self).__init__(
            penalty=None,
            fit_intercept=fit_intercept,
            max_iter=max_iter,
            tol=tol,
            shuffle=shuffle,
            verbose=verbose,
            random_state=random_state,
            eta0=1.0,
            warm_start=warm_start,
            class_weight=class_weight,
            average=average,
            n_jobs=n_jobs,
            n_iter=n_iter,
            iters=iters)

        self.C = C
        self.loss = loss

    def partial_fit(self, X, y, classes=None):
        """Fit linear model with Passive Aggressive algorithm.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Subset of the training data

        y : numpy array of shape [n_samples]
            Subset of the target values

        classes : array, shape = [n_classes]
            Classes across all calls to partial_fit.
            Can be obtained by via `np.unique(y_all)`, where y_all is the
            target vector of the entire dataset.
            This argument is required for the first call to partial_fit
            and can be omitted in the subsequent calls.
            Note that y doesn't need to contain all labels in `classes`.

        Returns
        -------
        self : returns an instance of self.
        """
        if self.class_weight == 'balanced':
            raise ValueError("class_weight 'balanced' is not supported for "
                             "partial_fit. For 'balanced' weights, use "
                             "`sklearn.utils.compute_class_weight` with "
                             "`class_weight='balanced'`. In place of y you "
                             "can use a large enough subset of the full "
                             "training set target to properly estimate the "
                             "class frequency distributions. Pass the "
                             "resulting weights as the class_weight "
                             "parameter.")
        lr = "pa1" if self.loss == "hinge" else "pa2"
        return self._partial_fit(X, y, alpha=1.0, C=self.C,
                                 loss="hinge", learning_rate=lr,
                                 max_iter=self.iters, classes=classes,
                                 sample_weight=None, coef_init=None,
                                 intercept_init=None)

    def fit(self, X, y, coef_init=None, intercept_init=None):
        """Fit linear model with Passive Aggressive algorithm.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training data

        y : numpy array of shape [n_samples]
            Target values

        coef_init : array, shape = [n_classes,n_features]
            The initial coefficients to warm-start the optimization.

        intercept_init : array, shape = [n_classes]
            The initial intercept to warm-start the optimization.

        Returns
        -------
        self : returns an instance of self.
        """
        lr = "pa1" if self.loss == "hinge" else "pa2"
        return self._fit(X, y, alpha=1.0, C=self.C,
                         loss="hinge", learning_rate=lr,
                         coef_init=coef_init, intercept_init=intercept_init)


class PassiveAggressiveRegressor(BaseSGDRegressor):
    """Passive Aggressive Regressor

    Read more in the :ref:`User Guide <passive_aggressive>`.

    Parameters
    ----------

    C : float
        Maximum step size (regularization). Defaults to 1.0.

    fit_intercept : bool
        Whether the intercept should be estimated or not. If False, the
        data is assumed to be already centered. Defaults to True.

    max_iter : int, optional
        The maximum number of passes over the training data (aka epochs).
        It only impacts the behavior in the ``fit`` method, and not the
        `partial_fit`.
        Defaults to 5. Defaults to 1000 from 0.21, or if tol is not None.

        .. versionadded:: 0.19

    tol : float or None, optional
        The stopping criterion. If it is not None, the iterations will stop
        when (loss > previous_loss - tol). Defaults to None.
        Defaults to 1e-3 from 0.21.

        .. versionadded:: 0.19

    shuffle : bool, default=True
        Whether or not the training data should be shuffled after each epoch.

    verbose : integer, optional
        The verbosity level

    loss : string, optional
        The loss function to be used:
        epsilon_insensitive: equivalent to PA-I in the reference paper.
        squared_epsilon_insensitive: equivalent to PA-II in the reference
        paper.

    epsilon : float
        If the difference between the current prediction and the correct label
        is below this threshold, the model is not updated.

    random_state : int, RandomState instance or None, optional, default=None
        The seed of the pseudo random number generator to use when shuffling
        the data.  If int, random_state is the seed used by the random number
        generator; If RandomState instance, random_state is the random number
        generator; If None, the random number generator is the RandomState
        instance used by `np.random`.

    warm_start : bool, optional
        When set to True, reuse the solution of the previous call to fit as
        initialization, otherwise, just erase the previous solution.
        See :term:`the Glossary <warm_start>`.

        Repeatedly calling fit or partial_fit when warm_start is True can
        result in a different solution than when calling fit a single time
        because of the way the data is shuffled.

    average : bool or int, optional
        When set to True, computes the averaged SGD weights and stores the
        result in the ``coef_`` attribute. If set to an int greater than 1,
        averaging will begin once the total number of samples seen reaches
        average. So average=10 will begin averaging after seeing 10 samples.

        .. versionadded:: 0.19
           parameter *average* to use weights averaging in SGD

    n_iter : int, optional
        The number of passes over the training data (aka epochs).
        Defaults to None. Deprecated, will be removed in 0.21.

        .. versionchanged:: 0.19
            Deprecated

    Attributes
    ----------
    coef_ : array, shape = [1, n_features] if n_classes == 2 else [n_classes,\
            n_features]
        Weights assigned to the features.

    intercept_ : array, shape = [1] if n_classes == 2 else [n_classes]
        Constants in decision function.

    n_iter_ : int
        The actual number of iterations to reach the stopping criterion.

    Examples
    --------
    >>> from sklearn.linear_model import PassiveAggressiveRegressor
    >>> from sklearn.datasets import make_regression
    >>>
    >>> X, y = make_regression(n_features=4, random_state=0)
    >>> regr = PassiveAggressiveRegressor(random_state=0)
    >>> regr.fit(X, y)
    PassiveAggressiveRegressor(C=1.0, average=False, epsilon=0.1,
                  fit_intercept=True, iters=1, loss='epsilon_insensitive',
                  max_iter=None, n_iter=None, random_state=0, shuffle=True,
                  tol=None, verbose=0, warm_start=False)
    >>> print(regr.coef_)
    [20.48736655 34.18818427 67.59122734 87.94731329]
    >>> print(regr.intercept_)
    [-0.02306214]
    >>> print(regr.predict([[0, 0, 0, 0]]))
    [-0.02306214]

    See also
    --------

    SGDRegressor

    References
    ----------
    Online Passive-Aggressive Algorithms
    <http://jmlr.csail.mit.edu/papers/volume7/crammer06a/crammer06a.pdf>
    K. Crammer, O. Dekel, J. Keshat, S. Shalev-Shwartz, Y. Singer - JMLR (2006)

    """
    def __init__(self, C=1.0, fit_intercept=True, max_iter=None, tol=None,
                 shuffle=True, verbose=0, loss="epsilon_insensitive",
                 epsilon=DEFAULT_EPSILON, random_state=None, warm_start=False,
                 average=False, n_iter=None, iters=1):
        super(PassiveAggressiveRegressor, self).__init__(
            penalty=None,
            l1_ratio=0,
            epsilon=epsilon,
            eta0=1.0,
            fit_intercept=fit_intercept,
            max_iter=max_iter,
            tol=tol,
            shuffle=shuffle,
            verbose=verbose,
            random_state=random_state,
            warm_start=warm_start,
            average=average,
            n_iter=n_iter,
            iters=iters)
        self.C = C
        self.loss = loss

    def partial_fit(self, X, y):
        """Fit linear model with Passive Aggressive algorithm.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Subset of training data

        y : numpy array of shape [n_samples]
            Subset of target values

        Returns
        -------
        self : returns an instance of self.
        """
        self._validate_params()
        lr = "pa1" if self.loss == "epsilon_insensitive" else "pa2"
        return self._partial_fit(X, y, alpha=1.0, C=self.C,
                                 loss="epsilon_insensitive",
                                 learning_rate=lr, max_iter=self.iters,
                                 sample_weight=None,
                                 coef_init=None, intercept_init=None)

    def fit(self, X, y, coef_init=None, intercept_init=None):
        """Fit linear model with Passive Aggressive algorithm.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training data

        y : numpy array of shape [n_samples]
            Target values

        coef_init : array, shape = [n_features]
            The initial coefficients to warm-start the optimization.

        intercept_init : array, shape = [1]
            The initial intercept to warm-start the optimization.

        Returns
        -------
        self : returns an instance of self.
        """
        lr = "pa1" if self.loss == "epsilon_insensitive" else "pa2"
        return self._fit(X, y, alpha=1.0, C=self.C,
                         loss="epsilon_insensitive",
                         learning_rate=lr,
                         coef_init=coef_init,
                         intercept_init=intercept_init)
