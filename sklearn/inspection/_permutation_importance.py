"""Permutation importance for estimators"""
import numpy as np
from joblib import Parallel
from joblib import delayed

from ..metrics import check_scoring
from ..utils import Bunch
from ..utils import check_random_state
from ..utils import check_array
from ..utils.validation import _deprecate_positional_args


def _calculate_permutation_scores(estimator, X, y, col_idx, random_state,
                                  n_repeats, scorer):
    """Calculate score when `col_idx` is permuted."""
    random_state = check_random_state(random_state)

    # Work on a copy of X to to ensure thread-safety in case of threading based
    # parallelism. Furthermore, making a copy is also useful when the joblib
    # backend is 'loky' (default) or the old 'multiprocessing': in those cases,
    # if X is large it will be automatically be backed by a readonly memory map
    # (memmap). X.copy() on the other hand is always guaranteed to return a
    # writable data-structure whose columns can be shuffled inplace.
    X_permuted = X.copy()
    scores = np.zeros(n_repeats)
    shuffling_idx = np.arange(X.shape[0])
    for n_round in range(n_repeats):
        random_state.shuffle(shuffling_idx)
        if hasattr(X_permuted, "iloc"):
            col = X_permuted.iloc[shuffling_idx, col_idx]
            col.index = X_permuted.index
            X_permuted.iloc[:, col_idx] = col
        else:
            X_permuted[:, col_idx] = X_permuted[shuffling_idx, col_idx]
        feature_score = scorer(estimator, X_permuted, y)
        scores[n_round] = feature_score

    return scores


@_deprecate_positional_args
def permutation_importance(estimator, X, y, *, scoring=None, n_repeats=5,
                           n_jobs=None, random_state=None):
    """Permutation importance for feature evaluation [BRE]_.

    The :term:`estimator` is required to be a fitted estimator. `X` can be the
    data set used to train the estimator or a hold-out set. The permutation
    importance of a feature is calculated as follows. First, a baseline metric,
    defined by :term:`scoring`, is evaluated on a (potentially different)
    dataset defined by the `X`. Next, a feature column from the validation set
    is permuted and the metric is evaluated again. The permutation importance
    is defined to be the difference between the baseline metric and metric from
    permutating the feature column.

    Read more in the :ref:`User Guide <permutation_importance>`.

    Parameters
    ----------
    estimator : object
        An estimator that has already been :term:`fitted` and is compatible
        with :term:`scorer`.

    X : ndarray or DataFrame, shape (n_samples, n_features)
        Data on which permutation importance will be computed.

    y : array-like or None, shape (n_samples, ) or (n_samples, n_classes)
        Targets for supervised or `None` for unsupervised.

    scoring : string, callable or None, default=None
        Scorer to use. It can be a single
        string (see :ref:`scoring_parameter`) or a callable (see
        :ref:`scoring`). If None, the estimator's default scorer is used.

    n_repeats : int, default=5
        Number of times to permute a feature.

    n_jobs : int or None, default=None
        The number of jobs to use for the computation.
        `None` means 1 unless in a :obj:`joblib.parallel_backend` context.
        `-1` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    random_state : int, RandomState instance, default=None
        Pseudo-random number generator to control the permutations of each
        feature.
        Pass an int to get reproducible results across function calls.
        See :term: `Glossary <random_state>`.

    Returns
    -------
    result : :class:`~sklearn.utils.Bunch`
        Dictionary-like object, with the following attributes.

        importances_mean : ndarray, shape (n_features, )
            Mean of feature importance over `n_repeats`.
        importances_std : ndarray, shape (n_features, )
            Standard deviation over `n_repeats`.
        importances : ndarray, shape (n_features, n_repeats)
            Raw permutation importance scores.

    References
    ----------
    .. [BRE] L. Breiman, "Random Forests", Machine Learning, 45(1), 5-32,
             2001. https://doi.org/10.1023/A:1010933404324

    Examples
    --------
    >>> from sklearn.linear_model import LogisticRegression
    >>> from sklearn.inspection import permutation_importance
    >>> X = [[1, 9, 9],[1, 9, 9],[1, 9, 9],
    ...      [0, 9, 9],[0, 9, 9],[0, 9, 9]]
    >>> y = [1, 1, 1, 0, 0, 0]
    >>> clf = LogisticRegression().fit(X, y)
    >>> result = permutation_importance(clf, X, y, n_repeats=10,
    ...                                 random_state=0)
    >>> result.importances_mean
    array([0.4666..., 0.       , 0.       ])
    >>> result.importances_std
    array([0.2211..., 0.       , 0.       ])
    """
    if not hasattr(X, "iloc"):
        X = check_array(X, force_all_finite='allow-nan', dtype=None)

    # Precompute random seed from the random state to be used
    # to get a fresh independent RandomState instance for each
    # parallel call to _calculate_permutation_scores, irrespective of
    # the fact that variables are shared or not depending on the active
    # joblib backend (sequential, thread-based or process-based).
    random_state = check_random_state(random_state)
    random_seed = random_state.randint(np.iinfo(np.int32).max + 1)

    scorer = check_scoring(estimator, scoring=scoring)
    baseline_score = scorer(estimator, X, y)

    scores = Parallel(n_jobs=n_jobs)(delayed(_calculate_permutation_scores)(
        estimator, X, y, col_idx, random_seed, n_repeats, scorer
    ) for col_idx in range(X.shape[1]))

    importances = baseline_score - np.array(scores)
    return Bunch(importances_mean=np.mean(importances, axis=1),
                 importances_std=np.std(importances, axis=1),
                 importances=importances)
