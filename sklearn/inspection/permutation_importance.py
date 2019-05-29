"""Permutation importance for estimators"""
import numpy as np
import scipy.sparse as sp

from ..utils import check_random_state
from ..utils.validation import check_array
from ..metrics import check_scoring


def permutation_importance(estimator, X, y, scoring=None, n_rounds=1,
                           random_state=None):
    """Permutation importance for feature evaluation. [BRE]_

    The permutation importance of a feature is calculated as follows. First,
    the estimator is trained on a training set. Then a baseline metric, defined
    by ``scoring``, is evaluated on a (potentially different) dataset defined
    by the ``X`` parameter. Next, a feature column from the validation set is
    permuted and the metric is evaluated again. The permutation importance is
    defined to be the difference between the baseline metric and metric from
    permutating the feature column.

    Read more in the :ref:`User Guide <permutation_importance>`.

    Parameters
    ----------
    estimator : object
        An estimator that has already been `fit` and is compatible with
        ``scorer``.

    X : array-like or DataFrame, shape = (n_samples, n_features)
        Data on which permutaiton importance will be computed.

    y : array-like, shape = (n_samples, ...)
        Targets for supervised learning.

    scoring : string, callable or None, optional (default=None)
        Scorer to use. It can be a single
        string (see :ref:`scoring_parameter`) or a callable (see
        :ref:`scoring`). If None, the estimator's default scorer is used.

    n_rounds : int, optional (default=1)
        Number of times to permute a feature.

    random_state : int, RandomState instance or None, optional, default None
        Pseudo-random number generator to control the permutations.
        See :term:`random_state`.

    Returns
    -------
    importances : array, shape (n_features, n_rounds)
        Permutation importance scores.

    References
    ----------
    .. [BRE] L. Breiman, "Random Forests", Machine Learning, 45(1), 5-32,
        2001.https://doi.org/10.1023/A:1010933404324
    """
    if not (hasattr(X, '__array__') or sp.issparse(X)):
        X = check_array(X, force_all_finite='allow-nan', dtype=np.object)

    n_features = X.shape[1]

    random_state = check_random_state(random_state)
    scorer = check_scoring(estimator, scoring=scoring)
    scores = np.empty(shape=(n_features, n_rounds), dtype=np.float)

    if hasattr(X, 'iloc'):  # pandas dataframe
        # Need to copy in case where the dataframe is a view
        X = X.copy()
        X_iloc = X.iloc
    else:
        X_iloc = X

    baseline_score = scorer(estimator, X, y)
    for col_idx in range(n_features):
        original_feature = X_iloc[:, col_idx].copy()

        for n_round in range(n_rounds):
            X_perm = random_state.permutation(original_feature)
            X_iloc[:, col_idx] = X_perm

            feature_score = scorer(estimator, X, y)
            scores[col_idx, n_round] = baseline_score - feature_score

        X_iloc[:, col_idx] = original_feature

    return scores
