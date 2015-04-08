# Authors: Rob Zinkov, Mathieu Blondel
# License: BSD 3 clause

from .stochastic_gradient import BaseSGDClassifier
from .stochastic_gradient import BaseSGDRegressor
from .stochastic_gradient import DEFAULT_EPSILON


class PassiveAggressiveClassifier(BaseSGDClassifier):
    """Passive Aggressive Classifier

    Parameters
    ----------

    C : float
        Maximum step size (regularization). Defaults to 1.0.

    fit_intercept : bool, default=False
        Whether the intercept should be estimated or not. If False, the
        data is assumed to be already centered.

    n_iter : int, optional
        The number of passes over the training data (aka epochs).
        Defaults to 5.

    shuffle : bool, default=True
        Whether or not the training data should be shuffled after each epoch.

    random_state : int seed, RandomState instance, or None (default)
        The seed of the pseudo random number generator to use when
        shuffling the data.

    verbose : integer, optional
        The verbosity level

    n_jobs : integer, optional
        The number of CPUs to use to do the OVA (One Versus All, for
        multi-class problems) computation. -1 means 'all CPUs'. Defaults
        to 1.

    loss : string, optional
        The loss function to be used:
        hinge: equivalent to PA-I in the reference paper.
        squared_hinge: equivalent to PA-II in the reference paper.

    warm_start : bool, optional
        When set to True, reuse the solution of the previous call to fit as
        initialization, otherwise, just erase the previous solution.

    Attributes
    ----------
    coef_ : array, shape = [1, n_features] if n_classes == 2 else [n_classes,\
            n_features]
        Weights assigned to the features.

    intercept_ : array, shape = [1] if n_classes == 2 else [n_classes]
        Constants in decision function.

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
    def __init__(self, C=1.0, fit_intercept=True,
                 n_iter=5, shuffle=True, verbose=0, loss="hinge",
                 n_jobs=1, random_state=None, warm_start=False):
        BaseSGDClassifier.__init__(self,
                                   penalty=None,
                                   fit_intercept=fit_intercept,
                                   n_iter=n_iter,
                                   shuffle=shuffle,
                                   verbose=verbose,
                                   random_state=random_state,
                                   eta0=1.0,
                                   warm_start=warm_start,
                                   n_jobs=n_jobs)
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
        lr = "pa1" if self.loss == "hinge" else "pa2"
        return self._partial_fit(X, y, alpha=1.0, C=self.C,
                                 loss="hinge", learning_rate=lr, n_iter=1,
                                 classes=classes, sample_weight=None,
                                 coef_init=None, intercept_init=None)

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

    Parameters
    ----------

    C : float
        Maximum step size (regularization). Defaults to 1.0.

    epsilon : float
        If the difference between the current prediction and the correct label
        is below this threshold, the model is not updated.

    fit_intercept : bool
        Whether the intercept should be estimated or not. If False, the
        data is assumed to be already centered. Defaults to True.

    n_iter : int, optional
        The number of passes over the training data (aka epochs).
        Defaults to 5.

    shuffle : bool, default=True
        Whether or not the training data should be shuffled after each epoch.

    random_state : int seed, RandomState instance, or None (default)
        The seed of the pseudo random number generator to use when
        shuffling the data.

    verbose : integer, optional
        The verbosity level

    loss : string, optional
        The loss function to be used:
        epsilon_insensitive: equivalent to PA-I in the reference paper.
        squared_epsilon_insensitive: equivalent to PA-II in the reference
        paper.

    warm_start : bool, optional
        When set to True, reuse the solution of the previous call to fit as
        initialization, otherwise, just erase the previous solution.

    Attributes
    ----------
    coef_ : array, shape = [1, n_features] if n_classes == 2 else [n_classes,\
            n_features]
        Weights assigned to the features.

    intercept_ : array, shape = [1] if n_classes == 2 else [n_classes]
        Constants in decision function.

    See also
    --------

    SGDRegressor

    References
    ----------
    Online Passive-Aggressive Algorithms
    <http://jmlr.csail.mit.edu/papers/volume7/crammer06a/crammer06a.pdf>
    K. Crammer, O. Dekel, J. Keshat, S. Shalev-Shwartz, Y. Singer - JMLR (2006)

    """
    def __init__(self, C=1.0, fit_intercept=True, n_iter=5, shuffle=True,
                 verbose=0, loss="epsilon_insensitive",
                 epsilon=DEFAULT_EPSILON, random_state=None, class_weight=None,
                 warm_start=False):
        BaseSGDRegressor.__init__(self,
                                  penalty=None,
                                  l1_ratio=0,
                                  epsilon=epsilon,
                                  eta0=1.0,
                                  fit_intercept=fit_intercept,
                                  n_iter=n_iter,
                                  shuffle=shuffle,
                                  verbose=verbose,
                                  random_state=random_state,
                                  warm_start=warm_start)
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
        lr = "pa1" if self.loss == "epsilon_insensitive" else "pa2"
        return self._partial_fit(X, y, alpha=1.0, C=self.C,
                                 loss="epsilon_insensitive",
                                 learning_rate=lr, n_iter=1,
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
