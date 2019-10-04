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

    max_iter : int, optional (default=1000)
        The maximum number of passes over the training data (aka epochs).
        It only impacts the behavior in the ``fit`` method, and not the
        :meth:`partial_fit` method.

        .. versionadded:: 0.19

    tol : float or None, optional (default=1e-3)
        The stopping criterion. If it is not None, the iterations will stop
        when (loss > previous_loss - tol).

        .. versionadded:: 0.19

    shuffle : bool, default=True
        Whether or not the training data should be shuffled after each epoch.

    verbose : integer, default=0
        The verbosity level

    eta0 : double
        Constant by which the updates are multiplied. Defaults to 1.

    n_jobs : int or None, optional (default=None)
        The number of CPUs to use to do the OVA (One Versus All, for
        multi-class problems) computation.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    random_state : int, RandomState instance or None, optional, default None
        The seed of the pseudo random number generator to use when shuffling
        the data.  If int, random_state is the seed used by the random number
        generator; If RandomState instance, random_state is the random number
        generator; If None, the random number generator is the RandomState
        instance used by `np.random`.

    early_stopping : bool, default=False
        Whether to use early stopping to terminate training when validation.
        score is not improving. If set to True, it will automatically set aside
        a stratified fraction of training data as validation and terminate
        training when validation score is not improving by at least tol for
        n_iter_no_change consecutive epochs.

        .. versionadded:: 0.20

    validation_fraction : float, default=0.1
        The proportion of training data to set aside as validation set for
        early stopping. Must be between 0 and 1.
        Only used if early_stopping is True.

        .. versionadded:: 0.20

    n_iter_no_change : int, default=5
        Number of iterations with no improvement to wait before early stopping.

        .. versionadded:: 0.20

    class_weight : dict, {class_label: weight} or "balanced" or None, optional
        Preset for the class_weight fit parameter.

        Weights associated with classes. If not given, all classes
        are supposed to have weight one.

        The "balanced" mode uses the values of y to automatically adjust
        weights inversely proportional to class frequencies in the input data
        as ``n_samples / (n_classes * np.bincount(y))``

    warm_start : bool, default=False
        When set to True, reuse the solution of the previous call to fit as
        initialization, otherwise, just erase the previous solution. See
        :term:`the Glossary <warm_start>`.

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

    classes_ : array of shape (n_classes,)
        The unique classes labels.

    t_ : int
        Number of weight updates performed during training.
        Same as ``(n_iter_ * n_samples)``.

    Notes
    -----

    ``Perceptron`` is a classification algorithm which shares the same
    underlying implementation with ``SGDClassifier``. In fact,
    ``Perceptron()`` is equivalent to `SGDClassifier(loss="perceptron",
    eta0=1, learning_rate="constant", penalty=None)`.

    Examples
    --------
    >>> from sklearn.datasets import load_digits
    >>> from sklearn.linear_model import Perceptron
    >>> X, y = load_digits(return_X_y=True)
    >>> clf = Perceptron(tol=1e-3, random_state=0)
    >>> clf.fit(X, y)
    Perceptron()
    >>> clf.score(X, y)
    0.939...

    See also
    --------

    SGDClassifier

    References
    ----------

    https://en.wikipedia.org/wiki/Perceptron and references therein.
    """
    def __init__(self, penalty=None, alpha=0.0001, fit_intercept=True,
                 max_iter=1000, tol=1e-3, shuffle=True, verbose=0, eta0=1.0,
                 n_jobs=None, random_state=0, early_stopping=False,
                 validation_fraction=0.1, n_iter_no_change=5,
                 class_weight=None, warm_start=False):
        super().__init__(
            loss="perceptron", penalty=penalty, alpha=alpha, l1_ratio=0,
            fit_intercept=fit_intercept, max_iter=max_iter, tol=tol,
            shuffle=shuffle, verbose=verbose, random_state=random_state,
            learning_rate="constant", eta0=eta0, early_stopping=early_stopping,
            validation_fraction=validation_fraction,
            n_iter_no_change=n_iter_no_change, power_t=0.5,
            warm_start=warm_start, class_weight=class_weight, n_jobs=n_jobs)
