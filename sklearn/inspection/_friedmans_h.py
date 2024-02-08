"""Friedman and Popescu's H-squared statistics"""

import itertools

import numpy as np

from ..ensemble._bagging import _generate_indices
from ..utils import _get_column_indices, _safe_assign, _safe_indexing
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
    except:
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


def hstatistics(
    estimator, X, *, features=None, n_max=500, random_state=None, sample_weight=None
):
    """Friedman and Popescu's H-statistics of pairwise interaction strength.

    For each feature pair, Friedman and Popescu's H-squared statistic [FRI]_ of
    interaction is calculated. It equals the proportion of effect variability
    of the two features unexplained by the main effects. Besides the "official"
    H-squared statistic, also the unnormalized statistic is returned. Its root
    is on the scale of the predictions and can directly compared between feature
    pairs.

    The complexity of the function is of O(n^2 p^2), where n is the number of
    data rows and p is the number of features considered.
    The size of n is automatically controlled via `n_max=500`, while it is
    your responsibility to pass only 2-5 *important* features or features of
    special interest.

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

    Returns
    -------
    result : A list with a tuple per feature pair. Each element contains the two
       feature names/indices, the unnormalized H-squared statistic, and the usual
       (normalized) H-squared statistic.

    References
    ----------
    .. [FRI] :doi:`J. H. Friedman and B. E. Popescu,
            "Predictive Learning via Rule Ensembles",
            The Annals of Applied Statistics, 2(3), 916-954,
            2008. <10.1214/07-AOAS148>`

    Examples
    ----------

    >>> from sklearn.ensemble import HistGradientBoostingRegressor
    >>> from sklearn.inspection import permutation_importance
    >>> from sklearn.datasets import load_diabetes

    >>> X, y = load_diabetes(return_X_y=True)
    >>> est = HistGradientBoostingRegressor().fit(X, y)

    >>> # Get Friedman's H-squared for top three predictors
    >>> imp = permutation_importance(est, X, y, n_repeats=10, random_state=0)
    >>> top_3 = np.argsort(imp.importances_mean)[-3:]
    >>> hstatistics(est, X=X, features=top_3, random_state=4)

    >>> For features 3 and 2, ~4% of joint effect variability comes from interaction:
    >>> [(3, 2, 50.42733968603663, 0.043352605995576346),
    >>> (3, 8, 20.34764277579143, 0.015648706966531766),
    >>> (2, 8, 56.778288912326026, 0.025823266523399196)]

    """
    check_is_fitted(estimator)

    if sample_weight is not None:
        sample_weight = _check_sample_weight(sample_weight, X)

    # Usually, the data is too large and we need subsampling
    if X.shape[0] > n_max:
        row_indices = _generate_indices(
            random_state=random_state,
            bootstrap=False,
            n_population=X.shape[0],
            n_samples=n_max,
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
            _get_column_indices(X, features), dtype="int32", order="C"
        ).ravel()

    # CALCULATIONS
    pd_univariate = []
    for ind in feature_indices:
        pd_univariate.append(
            _calculate_pd_over_data(
                estimator, X=X, feature_indices=[ind], sample_weight=sample_weight
            )
        )

    stats = []
    for i, j in itertools.combinations(range(len(feature_indices)), 2):
        pd_bivariate = _calculate_pd_over_data(
            estimator,
            X=X,
            feature_indices=feature_indices[[i, j]],
            sample_weight=sample_weight,
        )
        numerator = np.average(
            (pd_bivariate - pd_univariate[i] - pd_univariate[j]) ** 2,
            axis=0,
            weights=sample_weight,
        )
        denominator = np.average(pd_bivariate**2, axis=0, weights=sample_weight)
        stats.append((features[i], features[j], numerator, numerator / denominator))

    return stats
