# Author: Mathieu Blondel
# License: BSD 3 clause

from .stochastic_gradient import BaseSGDClassifier
from ..feature_selection.from_model import _LearntSelectorMixin


class Perceptron(BaseSGDClassifier, _LearntSelectorMixin):
    """Perceptron

    Read more in the :ref:`User Guide <perceptron>`.

    Parameters
    ----------

    penalty : None, 'l2' or 'l1' or 'elasticnet'
        The penalty (aka regularization term) to be used. Defaults to None.

    alpha : float
        Constant that multiplies the regularization term if regularization is
        used. Defaults to 0.0001

    fit_intercept : bool
        Whether the intercept should be estimated or not. If False, the
        data is assumed to be already centered. Defaults to True.

    n_iter : int, optional
        The number of passes over the training data (aka epochs).
        Defaults to 5.

    shuffle : bool, optional, default True
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

    eta0 : double
        Constant by which the updates are multiplied. Defaults to 1.

    class_weight : dict, {class_label: weight} or "balanced" or None, optional
        Preset for the class_weight fit parameter.

        Weights associated with classes. If not given, all classes
        are supposed to have weight one.

        The "balanced" mode uses the values of y to automatically adjust
        weights inversely proportional to class frequencies in the input data
        as ``n_samples / (n_classes * np.bincount(y))``

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

    Notes
    -----

    `Perceptron` and `SGDClassifier` share the same underlying implementation.
    In fact, `Perceptron()` is equivalent to `SGDClassifier(loss="perceptron",
    eta0=1, learning_rate="constant", penalty=None)`.

    See also
    --------

    SGDClassifier

    References
    ----------

    https://en.wikipedia.org/wiki/Perceptron and references therein.
    """
    def __init__(self, penalty=None, alpha=0.0001, fit_intercept=True,
                 n_iter=5, shuffle=True, verbose=0, eta0=1.0, n_jobs=1,
                 random_state=0, class_weight=None, warm_start=False):
        super(Perceptron, self).__init__(loss="perceptron",
                                         penalty=penalty,
                                         alpha=alpha, l1_ratio=0,
                                         fit_intercept=fit_intercept,
                                         n_iter=n_iter,
                                         shuffle=shuffle,
                                         verbose=verbose,
                                         random_state=random_state,
                                         learning_rate="constant",
                                         eta0=eta0,
                                         power_t=0.5,
                                         warm_start=warm_start,
                                         class_weight=class_weight,
                                         n_jobs=n_jobs)
