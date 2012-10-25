# Author: Rob Zinkov
# License: BSD Style.

from .stochastic_gradient import SGDClassifier, SGDRegressor


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

    learning_rate : string, optional
        The learning rate:
        pa1: eta = min(alpha, loss/norm(x))
        pa2: eta = 1.0 / (norm(x) + 0.5*alpha)

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
    def __init__(self, C=0.0001, fit_intercept=True,
                 n_iter=5, shuffle=False, verbose=0, learning_rate="pa1",
                 n_jobs=1, random_state=0, class_weight=None, warm_start=False):
        super(PassiveAggressiveClassifier, self).__init__(loss="hinge",
                                                          penalty=None,
                                                          C=C, l1_ratio=0,
                                                          fit_intercept=fit_intercept,
                                                          n_iter=n_iter,
                                                          shuffle=shuffle,
                                                          verbose=verbose,
                                                          random_state=random_state,
                                                          learning_rate=learning_rate,
                                                          warm_start=warm_start,
                                                          class_weight=class_weight,
                                                          n_jobs=n_jobs)


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

    n_jobs: integer, optional
        The number of CPUs to use to do the OVA (One Versus All, for
        multi-class problems) computation. -1 means 'all CPUs'. Defaults
        to 1.

    learning_rate : string, optional
        The learning rate:
        pa1: eta = min(alpha, loss/norm(x))
        pa2: eta = 1.0 / (norm(x) + 0.5*alpha)

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
    def __init__(self, C=0.0001, fit_intercept=True,
                 n_iter=5, shuffle=False, verbose=0, learning_rate="pa1",
                 n_jobs=1, random_state=0, class_weight=None, warm_start=False):
        super(
            PassiveAggressiveRegressor, self).__init__(loss="epsilon_insensitive",
                                                       penalty=None,
                                                       C=C, l1_ratio=0,
                                                       fit_intercept=fit_intercept,
                                                       n_iter=n_iter,
                                                       shuffle=shuffle,
                                                       verbose=verbose,
                                                       random_state=random_state,
                                                       learning_rate=learning_rate,
                                                       warm_start=warm_start,
                                                       class_weight=class_weight,
                                                       n_jobs=n_jobs)
