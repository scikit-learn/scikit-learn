# Author: Rob Zinkov
# License: BSD Style.

from .stochastic_gradient import SGDClassifier
from .stochastic_gradient import SGDRegressor
from .stochastic_gradient import DEFAULT_EPSILON


class PassiveAggressiveClassifier(SGDClassifier):
    """Passive Aggressive Classifier

    Parameters
    ----------

    C : float
        C parameter that scales the regularization term. Defaults to 0.001

    fit_intercept: bool
        Whether the intercept should be estimated or not. If False, the
        data is assumed to be already centered. Defaults to True.

    n_iter: int, optional
        The number of passes over the training data (aka epochs).
        Defaults to 5.

    shuffle: bool, optional
        Whether or not the training data should be shuffled after each epoch.
        Defaults to False.

    random_state: int seed, RandomState instance, or None (default)
        The seed of the pseudo random number generator to use when
        shuffling the data.

    verbose: integer, optional
        The verbosity level

    n_jobs: integer, optional
        The number of CPUs to use to do the OVA (One Versus All, for
        multi-class problems) computation. -1 means 'all CPUs'. Defaults
        to 1.

    loss : string, optional
        The loss function to be used:
        hinge: equivalent to PA-I in the reference paper.
        squared_hinge: equivalent to PA-II in the reference paper.

    class_weight : dict, {class_label : weight} or "auto" or None, optional
        Preset for the class_weight fit parameter.

        Weights associated with classes. If not given, all classes
        are supposed to have weight one.

        The "auto" mode uses the values of y to automatically adjust
        weights inversely proportional to class frequencies.

    warm_start : bool, optional
        When set to True, reuse the solution of the previous call to fit as
        initialization, otherwise, just erase the previous solution.

    Attributes
    ----------
    `coef_` : array, shape = [1, n_features] if n_classes == 2 else [n_classes,
    n_features]
        Weights assigned to the features.

    `intercept_` : array, shape = [1] if n_classes == 2 else [n_classes]
        Constants in decision function.

    See also
    --------

    SGDClassifier
    Perceptron

    References
    ----------
    Online Passive-Aggressive Algorithms
    <http://jmlr.csail.mit.edu/papers/volume7/crammer06a/crammer06a.pdf>
    K. Crammer, O. Dekel, J. Keshat, S. Shalev-Shwartz, Y. Singer - JMLR 7 (2006)

    """
    def __init__(self, C=1.0, fit_intercept=True,
                 n_iter=5, shuffle=False, verbose=0, loss="hinge",
                 n_jobs=1, random_state=0, class_weight=None,
                 warm_start=False):
        self.C = C
        self.loss = loss
        SGDClassifier.__init__(self,
                               penalty=None,
                               fit_intercept=fit_intercept,
                               n_iter=n_iter,
                               shuffle=shuffle,
                               verbose=verbose,
                               random_state=random_state,
                               eta0=1.0,
                               warm_start=warm_start,
                               class_weight=class_weight,
                               n_jobs=n_jobs)

    def partial_fit(self, X, y, classes=None, sample_weight=None):
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

        sample_weight : array-like, shape = [n_samples], optional
            Weights applied to individual samples.
            If not provided, uniform weights are assumed.

        Returns
        -------
        self : returns an instance of self.
        """
        lr = "pa1" if self.loss == "hinge" else "pa2"
        return self._partial_fit(X, y, alpha=1.0, C=self.C,
                                 loss="hinge", learning_rate=lr, n_iter=1,
                                 classes=classes, sample_weight=sample_weight,
                                 coef_init=None, intercept_init=None)

    def fit(self, X, y, coef_init=None, intercept_init=None,
            class_weight=None, sample_weight=None):
        """Fit linear model with Passive Aggressive algorithm.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training data

        y : numpy array of shape [n_samples]
            Target values

        coef_init : array, shape = [n_classes,n_features]
            The initial coeffients to warm-start the optimization.

        intercept_init : array, shape = [n_classes]
            The initial intercept to warm-start the optimization.

        sample_weight : array-like, shape = [n_samples], optional
            Weights applied to individual samples.
            If not provided, uniform weights are assumed.

        Returns
        -------
        self : returns an instance of self.
        """
        lr = "pa1" if self.loss == "hinge" else "pa2"
        return self._fit(X, y, alpha=1.0, C=self.C,
                         loss="hinge", learning_rate=lr,
                         coef_init=coef_init, intercept_init=intercept_init,
                         class_weight=class_weight,
                         sample_weight=sample_weight)


class PassiveAggressiveRegressor(SGDRegressor):
    """Passive Aggressive Regressor

    Parameters
    ----------

    C : float
        C parameter that scales the regularization term. Defaults to 0.001

    fit_intercept: bool
        Whether the intercept should be estimated or not. If False, the
        data is assumed to be already centered. Defaults to True.

    n_iter: int, optional
        The number of passes over the training data (aka epochs).
        Defaults to 5.

    shuffle: bool, optional
        Whether or not the training data should be shuffled after each epoch.
        Defaults to False.

    random_state: int seed, RandomState instance, or None (default)
        The seed of the pseudo random number generator to use when
        shuffling the data.

    verbose: integer, optional
        The verbosity level

    loss : string, optional
        The loss function to be used:
        pa1: eta = min(alpha, loss/norm(x))
        pa2: eta = 1.0 / (norm(x) + 0.5*alpha)

    warm_start : bool, optional
        When set to True, reuse the solution of the previous call to fit as
        initialization, otherwise, just erase the previous solution.

    Attributes
    ----------
    `coef_` : array, shape = [1, n_features] if n_classes == 2 else [n_classes,
    n_features]
        Weights assigned to the features.

    `intercept_` : array, shape = [1] if n_classes == 2 else [n_classes]
        Constants in decision function.

    See also
    --------

    SGDClassifier
    Perceptron

    References
    ----------
    Online Passive-Aggressive Algorithms
    <http://jmlr.csail.mit.edu/papers/volume7/crammer06a/crammer06a.pdf>
    K. Crammer, O. Dekel, J. Keshat, S. Shalev-Shwartz, Y. Singer - JMLR 7 (2006)

    """
    def __init__(self, C=1.0, fit_intercept=True,
                 n_iter=5, shuffle=False, verbose=0, loss="pa1",
                 epsilon=DEFAULT_EPSILON,
                 random_state=0, class_weight=None, warm_start=False):
        self.C = C
        SGDRegressor.__init__(self,
                              loss="epsilon_insensitive",
                              penalty=None,
                              l1_ratio=0,
                              epsilon=epsilon,
                              eta0=1.0,
                              fit_intercept=fit_intercept,
                              n_iter=n_iter,
                              shuffle=shuffle,
                              verbose=verbose,
                              random_state=random_state,
                              learning_rate=loss,
                              warm_start=warm_start)

    def partial_fit(self, X, y, sample_weight=None):
        """Fit linear model with Passive Aggressive algorithm.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Subset of training data

        y : numpy array of shape [n_samples]
            Subset of target values

        sample_weight : array-like, shape = [n_samples], optional
            Weights applied to individual samples.
            If not provided, uniform weights are assumed.

        Returns
        -------
        self : returns an instance of self.
        """
        return self._partial_fit(X, y, alpha=1.0, C=self.C, n_iter=1,
                                 sample_weight=sample_weight)

    def fit(self, X, y, coef_init=None, intercept_init=None,
            sample_weight=None):
        """Fit linear model with Passive Aggressive algorithm.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training data

        y : numpy array of shape [n_samples]
            Target values

        coef_init : array, shape = [n_features]
            The initial coeffients to warm-start the optimization.

        intercept_init : array, shape = [1]
            The initial intercept to warm-start the optimization.

        sample_weight : array-like, shape = [n_samples], optional
            Weights applied to individual samples (1. for unweighted).

        Returns
        -------
        self : returns an instance of self.
        """
        return self._fit(X, y, alpha=1.0, C=self.C,
                         coef_init=coef_init,
                         intercept_init=intercept_init,
                         sample_weight=sample_weight)
