# Author: Mathieu Blondel
# License: BSD 3 clause

from .stochastic_gradient import BaseSGDClassifier


class Perceptron(BaseSGDClassifier):
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

    shuffle : bool, optional, default True
        Whether or not the training data should be shuffled after each epoch.

    verbose : integer, optional
        The verbosity level

    eta0 : double
        Constant by which the updates are multiplied. Defaults to 1.

    n_jobs : integer, optional
        The number of CPUs to use to do the OVA (One Versus All, for
        multi-class problems) computation. -1 means 'all CPUs'. Defaults
        to 1.

    random_state : int, RandomState instance or None, optional, default None
        The seed of the pseudo random number generator to use when shuffling
        the data.  If int, random_state is the seed used by the random number
        generator; If RandomState instance, random_state is the random number
        generator; If None, the random number generator is the RandomState
        instance used by `np.random`.

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
                 max_iter=None, tol=None, shuffle=True, verbose=0, eta0=1.0,
                 n_jobs=1, random_state=0, class_weight=None,
                 warm_start=False, n_iter=None):
        super(Perceptron, self).__init__(loss="perceptron",
                                         penalty=penalty,
                                         alpha=alpha, l1_ratio=0,
                                         fit_intercept=fit_intercept,
                                         max_iter=max_iter,
                                         tol=tol,
                                         shuffle=shuffle,
                                         verbose=verbose,
                                         random_state=random_state,
                                         learning_rate="constant",
                                         eta0=eta0,
                                         power_t=0.5,
                                         warm_start=warm_start,
                                         class_weight=class_weight,
                                         n_jobs=n_jobs,
                                         n_iter=n_iter)
