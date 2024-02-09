"""Friedman and Popescu's H-Statistic"""

import itertools

import numpy as np

from ..utils import Bunch, _get_column_indices, _safe_assign, _safe_indexing
from ..utils.random import sample_without_replacement
from ..utils.validation import _check_sample_weight, check_is_fitted


def _calculate_pd_over_data(estimator, X, feature_indices, sample_weight=None):
    """Calculates partial dependence over the data distribution.

    It returns a 1D or 2D numpy array of the same length as X.
    """

    # Select grid columns and remove duplicates (will compensate below)
    grid = _safe_indexing(X, feature_indices, axis=1)

    # Unfortunately, the next line does not work in all cases, especially not in
    # the important case of discrete pandas Dataframes.
    try:
        compressed = True
        grid, ix_reconstruct = np.unique(grid, return_inverse=True, axis=0)
    except TypeError:  # TODO Better solution
        compressed = False
    n_grid = grid.shape[0]

    # X is stacked n_grid times, and grid columns are replaced by replicated grid
    n = X.shape[0]
    X_stacked = _safe_indexing(X, np.tile(np.arange(n), n_grid), axis=0)
    _safe_assign(
        X_stacked, values=np.repeat(grid, n, axis=0), column_indexer=feature_indices
    )

    # Predict on stacked data
    if hasattr(estimator, "predict_proba"):
        preds = estimator.predict_proba(X_stacked)
    else:
        preds = estimator.predict(X_stacked)

    # Predictions are averaged by grid value, mapped to original row order, and centered
    averaged_predictions = np.array(
        [np.average(Z, axis=0, weights=sample_weight) for Z in np.split(preds, n_grid)]
    )
    if compressed:
        averaged_predictions = averaged_predictions[ix_reconstruct]
    column_means = np.average(averaged_predictions, axis=0, weights=sample_weight)
    return averaged_predictions - column_means


def h_statistic(
    estimator,
    X,
    *,
    features=None,
    n_max=500,
    random_state=None,
    sample_weight=None,
    eps=1e-10,
):
    """Friedman and Popescu's H-statistic of pairwise interaction strength.

    For each feature pair, Friedman and Popescu's H-statistic of
    interaction strength [FRI]_ is calculated. It equals the proportion of
    effect variability of the two features unexplained by their main effects.
    Besides the (normalized) H-squared statistic, also the unnormalized statistic
    is returned. Its root is on the scale of the predictions, and can directly
    be compared between feature pairs.

    The complexity of the function is :math:`O(n^2 n^2)`, where :math:`n` is
    the number of observations and :math:`p` is the number of features considered.
    The size of `n` is automatically controlled via `n_max=500`, while it is
    the user's responsibility to select some *important* features.

    Parameters
    ----------
    estimator : object
        An estimator that has already been :term:`fitted`.

    X : ndarray or DataFrame, shape (n_observations, n_features)
        Data for which :term:`estimator` is able to calculate predictions.

    features : array-like of {int, str}, default=None
        List of feature names or column indices used to calculate pairwise statistics.
        The default, None, will use all column indices of X.

    n_max : int, default=500
        The number of rows to draw without replacement from X (and `sample_weight`).

    random_state : int, RandomState instance, default=None
        Pseudo-random number generator used for subsampling via `n_max`.
        See :term:`Glossary <random_state>`.

    sample_weight : array-like of shape (n_samples,), default=None
        Sample weights used in calculating partial dependences.

    eps : float, default=1e-10
        Threshold below which numerator values are set to 0.

    Returns
    -------
    result : :class:`~sklearn.utils.Bunch`
        Dictionary-like object, with the following attributes.

        feature_pair : list of length n_feature_pairs
            The list contains tuples of feature pairs (indices) in the same order
            as all pairwise statistics.

        numerator_pairwise : ndarray of shape (n_pairs, ) or (n_pairs, output_dim)
            Numerator of pairwise H-squared statistic.
            Useful to see which feature pair has strongest absolute interaction.
            Take square-root to get values on the scale of the predictions.

        denominator_pairwise : ndarray of shape (n_pairs, ) or (n_pairs, output_dim)
            Denominator of pairwise H-squared statistic.

        hsquared_pairwise : ndarray of shape (n_pairs, ) or (n_pairs, output_dim)
            Pairwise H-squared statistic. Useful to see which feature pair has
            strongest relative interation (relative with respect to joint effect).
            Calculated as numerator_pairwise / denominator_pairwise.

    References
    ----------
    .. [FRI] :doi:`J. H. Friedman and B. E. Popescu,
            "Predictive Learning via Rule Ensembles",
            The Annals of Applied Statistics, 2(3), 916-954,
            2008. <10.1214/07-AOAS148>`

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.ensemble import HistGradientBoostingRegressor
    >>> from sklearn.inspection import permutation_importance, h_statistic
    >>> from sklearn.datasets import load_diabetes

    >>> X, y = load_diabetes(return_X_y=True)
    >>> est = HistGradientBoostingRegressor().fit(X, y)

    >>> # Get Friedman's H-squared for top three predictors
    >>> imp = permutation_importance(est, X, y, n_repeats=10, random_state=0)
    >>> top_3 = np.argsort(imp.importances_mean)[-3:]
    >>> h_statistic(est, X=X, features=top_3, random_state=4)

    >>> # For feature pair (3, 2), about 4% of joint effect variability comes from
    >>> # their interaction. Unnormalized statistics are highest for pair (2, 8):
    >>> # {'feature_pair': [(3, 2), (3, 8), (2, 8)],
    >>> # 'numerator_pairwise': array([50.42733969, 20.34764278, 56.77828891]),
    >>> # 'denominator_pairwise': array([1163.19050558, 1300.27629882, 2198.72605431]),
    >>> # 'h_squared_pairwise': array([0.04335261, 0.01564871, 0.02582327])}

    """
    check_is_fitted(estimator)

    if sample_weight is not None:
        sample_weight = _check_sample_weight(sample_weight, X)

    # Usually, the data is too large and we need subsampling
    if X.shape[0] > n_max:
        row_indices = sample_without_replacement(
            n_population=X.shape[0], n_samples=n_max, random_state=random_state
        )
        X = _safe_indexing(X, row_indices, axis=0)
        if sample_weight is not None:
            sample_weight = _safe_indexing(sample_weight, row_indices, axis=0)
    else:
        X = X.copy()

    # TODO: Improve logic, e.g., use column names if there are some
    if features is None:
        features = feature_indices = np.arange(X.shape[1])
    else:
        feature_indices = np.asarray(
            _get_column_indices(X, features), dtype=np.int32, order="C"
        ).ravel()

    # CALCULATIONS
    pd_univariate = []
    for ind in feature_indices:
        pd_univariate.append(
            _calculate_pd_over_data(
                estimator, X=X, feature_indices=[ind], sample_weight=sample_weight
            )
        )

    hstat_results = Bunch()
    hstat_results["feature_pair"] = list(itertools.combinations(features, 2))
    num = []
    denom = []

    for i, j in itertools.combinations(range(len(feature_indices)), 2):
        pd_bivariate = _calculate_pd_over_data(
            estimator,
            X=X,
            feature_indices=feature_indices[[i, j]],
            sample_weight=sample_weight,
        )
        num.append(
            np.average(
                (pd_bivariate - pd_univariate[i] - pd_univariate[j]) ** 2,
                axis=0,
                weights=sample_weight,
            )
        )
        denom.append(np.average(pd_bivariate**2, axis=0, weights=sample_weight))

    # Round small numerators to 0
    num = np.array(num)
    num[np.abs(num) < eps] = 0

    hstat_results["numerator_pairwise"] = num
    hstat_results["denominator_pairwise"] = np.array(denom)
    hstat_results["h_squared_pairwise"] = (
        hstat_results["numerator_pairwise"] / hstat_results["denominator_pairwise"]
    )

    return hstat_results
