# Author: Mathieu Blondel
# License: BSD Style.

from .stochastic_gradient import SGDClassifier

class Perceptron(SGDClassifier):
    """Perceptron

    Parameters
    ----------

    fit_intercept: bool
        Whether the intercept should be estimated or not. If False, the
        data is assumed to be already centered. Defaults to True.

    n_iter: int, optional
        The number of passes over the training data (aka epochs).
        Defaults to 5.

    shuffle: bool, optional
        Whether or not the training data should be shuffled after each epoch.
        Defaults to False.

    seed: int, optional
        The seed of the pseudo random number generator to use when
        shuffling the data.

    verbose: integer, optional
        The verbosity level

    n_jobs: integer, optional
        The number of CPUs to use to do the OVA (One Versus All, for
        multi-class problems) computation. -1 means 'all CPUs'. Defaults
        to 1.

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

    """
    def __init__(self, fit_intercept=True, n_iter=5, shuffle=False, verbose=0,
                 n_jobs=1, seed=0, class_weight=None, warm_start=False):
        super(Perceptron, self).__init__(loss="perceptron",
                                         penalty=None,
                                         alpha=0, rho=0,
                                         fit_intercept=fit_intercept,
                                         n_iter=n_iter,
                                         shuffle=shuffle,
                                         verbose=verbose,
                                         seed=seed,
                                         learning_rate="constant",
                                         eta0=1.0,
                                         power_t=0.5,
                                         warm_start=warm_start)
