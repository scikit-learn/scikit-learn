"""Permutation importance for estimators"""
from contextlib import contextmanager

import numpy as np

from ..base import is_classifier, clone
from ..utils import check_random_state
from ..utils._joblib import Parallel, delayed
from ..model_selection import check_cv
from ..metrics import check_scoring
from ..utils.metaestimators import _safe_split


@contextmanager
def _permute_column(X, column, random_state):
    """Context manager to permute a column"""
    original_feature = X[:, column].copy()
    X[:, column] = random_state.permutation(X[:, column])
    yield X
    X[:, column] = original_feature


def _fit_and_calcuate_permutation_importance(estimator, X, y, train_indices,
                                             test_indices, columns, scoring,
                                             random_state):
    """Fits and calculates permutation importance

    Fits ``estimator`` on ``X`` and ``y``

    Parameters
    ----------
    estimator : object
        A supervised learning estimator with a `fit` and is compatible with
        ``scorer``.

    X : array-like, shape = (n_samples, n_features)
        Training data.

    y : array-like, shape = (n_samples, ...)
        Target relative to ``X``.

    train_indices : array of int
        Train indicies.

    test_indices : array of int
        Test indices.

    columns : list of integers
        A list of columns to calculate the permutation importance. If `None`,
        all columns will be used.

    scoring : string, callable or None
        A string (see model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.

    random_state: : RandomState instance
        Random number generator.

    Returns
    -------
    permutation_importance_scores : list
        Permutation importance scores for each column on the validation set
        defined by ``test_indices``.
    """
    X_train, y_train = _safe_split(estimator, X, y, train_indices)
    X_test, y_test = _safe_split(estimator, X, y, test_indices, train_indices)

    estimator.fit(X_train, y_train)
    baseline_score = scoring(estimator, X_test, y_test)

    permutation_importance_scores = []
    for column in columns:
        with _permute_column(X_test, column, random_state) as X_perm:
            feature_score = scoring(estimator, X_perm, y_test)
            permutation_importance_scores.append(baseline_score -
                                                 feature_score)

    return permutation_importance_scores


def permutation_importance(estimator, X, y, columns=None, scoring=None, cv=5,
                           n_jobs=None, pre_dispatch='2*n_jobs',
                           random_state=None):
    """Permutation importance for feature evaluation.

    The permutation importance of a feature is calculated as follows. First,
    the estimator is trained on a training set. Then a baseline metric, defined
    by ``scoring``, is evaluated on a validation set. Next, a feature column
    from the validation set is permuted and evaluated again. The permutation
    importance is defined to be the difference between the baseline metric
    and metric from permutating the feature column.

    Parameters
    ----------
    estimator : object
        A supervised learning estimator with a `fit` and is compatible with
        ``scorer``.

    X : array-like, shape = (n_samples, n_features)
        Training data.

    y : array-like, shape = (n_samples, ...)
        Target relative to ``X``.

    columns : list of integers, optional (default=None)
        A list of columns to calculate the permutation importance. If `None`,
        all columns will be used

    scorer : string, callable or None, optional (default=None)
        A string (see model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.

    cv : int, cross-validation generator or an iterable, optional (default=5)
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:

        - integer, to specify the number of folds.
        - :term:`CV splitter`,
        - An iterable yielding (train, test) splits as arrays of indices.

        For integer/None inputs, :class:`KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

    n_jobs : int or None, optional (default=None)
        Number of CPUs to use during the cross validation.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    random_state : int, RandomState instance or None, optional, default None
        The seed of the pseudo random number generator that selects a random
        feature to update.  If int, random_state is the seed used by the random
        number generator; If RandomState instance, random_state is the random
        number generator; If None, the random number generator is the
        RandomState instance used by `np.random`.

    Returns
    -------

    permutation_importance_scores : array, shape (n_columns, n_cv)
        Permutation importance scores where the rows are ordered corresponding
        to the ``columns`` argument.
    """

    cv = check_cv(cv, y, classifier=is_classifier(estimator))
    random_state = check_random_state(random_state)
    scoring = check_scoring(estimator, scoring=scoring)

    parallel = Parallel(n_jobs=n_jobs, pre_dispatch=pre_dispatch)

    if columns is None:
        columns = range(0, X.shape[1])

    with parallel:
        permutation_importance_scores = parallel(
            delayed(_fit_and_calcuate_permutation_importance)(
                clone(estimator), X, y, train_indices,
                test_indices, columns, scoring, random_state
            ) for train_indices, test_indices in cv.split(X, y))

    return np.array(permutation_importance_scores).T
