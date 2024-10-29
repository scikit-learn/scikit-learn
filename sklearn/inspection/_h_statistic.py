"""Friedman and Popescu's H-Statistic"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import itertools

import numpy as np
from scipy import sparse

from ..base import is_classifier, is_regressor
from ..utils import Bunch, check_array
from ..utils._indexing import (
    _get_column_indices,
    _safe_assign,
    _safe_indexing,
    resample,
)
from ..utils._param_validation import (
    HasMethods,
    Integral,
    Interval,
    Real,
    validate_params,
)
from ..utils.validation import _check_sample_weight, _num_samples, check_is_fitted


def _calculate_pd_brute_fast(
    pred_fun, X, feature_indices, grid, sample_weight=None, reduce_binary=False
):
    """Fast version of _calculate_partial_dependence_brute()

    Returns np.array of size (n_grid, output_dim).
    """

    # X is stacked n_grid times, and grid columns are replaced by replicated grid
    n = X.shape[0]
    n_grid = grid.shape[0]

    X_stacked = _safe_indexing(X, np.tile(np.arange(n), n_grid), axis=0)
    grid_stacked = _safe_indexing(grid, np.repeat(np.arange(n_grid), n), axis=0)

    # TODO: remove via https://github.com/scikit-learn/scikit-learn/issues/28931
    if hasattr(X, "iloc"):  # pandas<2 does not allow "values" to have repeated indices
        grid_stacked = grid_stacked.reset_index(drop=True)
        X_stacked = X_stacked.reset_index(drop=True)
    _safe_assign(X_stacked, values=grid_stacked, column_indexer=feature_indices)

    preds = pred_fun(X_stacked)

    # Drop reduntant 0 class in binary classification
    if reduce_binary:
        preds = preds[:, 1]

    # Partial dependences are averages per grid block
    pd_values = np.average(
        preds.reshape(n_grid, preds.shape[0] // n_grid, -1),
        axis=1,
        weights=sample_weight,
    )

    return pd_values


def _calculate_pd_over_data(
    pred_fun, X, feature_indices, sample_weight=None, reduce_binary=False
):
    """Calculates centered partial dependence over the data distribution.

    It returns a numpy array of size (n, output_dim).
    """

    # Select grid columns and remove duplicates (will compensate below)
    grid = _safe_indexing(X, feature_indices, axis=1)

    # np.unique() fails for mixed type and sparse objects
    # TODO: check closer the AxisError and when it is available
    possible_error_types = (
        (TypeError, np.AxisError) if hasattr(np, "AxisError") else (TypeError,)
    )
    try:
        ax = 0 if grid.shape[1] > 1 else None  # np.unique works better in 1 dim
        _, ix, ix_reconstruct = np.unique(
            grid, return_index=True, return_inverse=True, axis=ax
        )
        ix, ix_reconstruct = ix.squeeze(), ix_reconstruct.squeeze()  # squeeze to 1D
        grid = _safe_indexing(grid, ix, axis=0)
        compressed_grid = True
    except possible_error_types:
        compressed_grid = False

    pd_values = _calculate_pd_brute_fast(
        pred_fun,
        X=X,
        feature_indices=feature_indices,
        grid=grid,
        sample_weight=sample_weight,
        reduce_binary=reduce_binary,
    )

    if compressed_grid:
        pd_values = pd_values[ix_reconstruct]

    # H-statistics are based on *centered* partial dependences
    column_means = np.average(pd_values, axis=0, weights=sample_weight)

    return pd_values - column_means


@validate_params(
    {
        "estimator": [
            HasMethods(["fit", "predict"]),
            HasMethods(["fit", "predict_proba"]),
        ],
        "X": ["array-like", "sparse matrix"],
        "features": ["array-like", list, None],
        "sample_weight": ["array-like", None],
        "subsample": [Interval(Integral, 1, None, closed="left")],
        "eps": [Interval(Real, 0, None, closed="left")],
        "random_state": ["random_state"],
    },
    prefer_skip_nested_validation=True,
)
def h_statistic(
    estimator,
    X,
    features=None,
    *,
    sample_weight=None,
    subsample=500,
    eps=1e-10,
    random_state=None,
):
    """Friedman and Popescu's H-statistic of pairwise interaction strength.

    Calculates Friedman and Popescu's H-statistic of interaction strength
    for each feature pair j, k, see [FRI]_. The statistic is defined as:

        H_jk^2 = Numerator_jk / Denominator_jk, where

        - Numerator_jk = 1/n * sum(PD_{jk}(x_ij, x_ik) - PD_j(x_ij) - PD_k(x_ik)^2,
        - Denominator_jk = 1/n * sum(PD_{jk}(x_ij, x_ik)^2),
        - PD_j and PD_jk are the one- and two-dimensional partial dependence
          functions centered to mean 0,
        - and the sums run over 1 <= i <= n, where n is the sample size.

    It equals the proportion of effect variability between two features unexplained
    by their main effects. When there is no interaction, the value is
    exactly 0. The numerator (or its square root) provides an absolute measure
    of interaction strength, enabling direct comparison across feature pairs.

    The computational complexity of the function is `O(p^2 n^2)`,
    where `p` denotes the number of features considered. The size of `n` is
    automatically controlled via `subsample=500`, while it is the user's responsibility
    to select only a subset of *important* features. It is crucial to focus on important
    features because for weak predictors, the denominator might be small, and
    even a weak interaction could result in a high Friedman's H, sometimes exceeding 1.

    Read more in the :ref:`User Guide <h_statistic>`.

    .. versionadded:: 1.6

    Parameters
    ----------
    estimator : object
        A fitted estimator.

    X : {array-like or dataframe} of shape (n_samples, n_features)
        Data for which :term:`estimator` is able to calculate predictions.

    features : array-like of {int, str}, default=None
        List of feature names or column indices used to calculate pairwise statistics.
        The default, None, will use all column indices of X.

    sample_weight : array-like of shape (n_samples,), default=None
        Sample weights used in calculating partial dependencies.

    subsample : int, default=500
        Maximum number of samples drawn without replacement for computational
        efficiency. Reducing the number of samples improves speed but increases
        variance of the estimate.

    eps : float, default=1e-10
        Threshold below which numerator values are set to 0.

    random_state : int, RandomState instance, default=None
        Pseudo-random number generator used for subsampling via `subsample`.
        See :term:`Glossary <random_state>`.

    Returns
    -------
    result : :class:`~sklearn.utils.Bunch`
        Dictionary-like object, with the attributes listed below. Note that
        `output_dim` equals the number of values predicted per observation.
        For single-output regression and binary classification, `output_dim` is 1.

        feature_pairs : list of length n_feature_pairs
            The list contains tuples of feature pairs (indices) in the same order
            as all pairwise statistics.

        h_squared_pairwise : ndarray of shape (n_pairs, output_dim)
            Pairwise H-squared statistic. Useful to see which feature pair has
            strongest relative interaction (relative with respect to joint effect).
            Calculated as numerator_pairwise / denominator_pairwise.

        numerator_pairwise : ndarray of shape (n_pairs, output_dim)
            Numerator of pairwise H-squared statistic.
            Useful to see which feature pair has strongest absolute interaction.
            Take square-root to get values on the scale of the predictions.

        denominator_pairwise : ndarray of shape (n_pairs, output_dim)
            Denominator of pairwise H-squared statistic. Used for appropriate
            normalization of H.

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
    >>> est = HistGradientBoostingRegressor(max_iter=5, max_depth=2).fit(X, y)
    >>> m = 3
    >>> imp = permutation_importance(est, X, y, random_state=0)
    >>> top_m = np.argsort(imp.importances_mean)[-m:]
    >>> H = h_statistic(est, X=X, features=top_m, random_state=4)
    >>> # {'feature_pairs': [(3, 8), (3, 2), (8, 2)],
    >>> #  'h_squared_pairwise': array([[0.00985985],
    >>> #         [0.00927104],
    >>> #         [0.03439926]]),
    >>> #  'numerator_pairwise': array([[ 1.2955532 ],
    >>> #         [ 1.2419687 ],
    >>> #         [11.13358385]]),
    >>> #  'denominator_pairwise': array([[131.39690331],
    >>> #         [133.96210997],
    >>> #         [323.6576595 ]])}
    >>> # Interpretation: For features (8, 2), 3.4% of the joint effect variability
    >>> # comes from their interaction. These two features also have strongest absolute
    >>> # interaction across pairs ("numerator_pairwise").
    """
    reduce_binary = False  # In binary classification, we only retain second class probs

    check_is_fitted(estimator)

    if is_regressor(estimator):
        pred_fun = getattr(estimator, "predict", None)
        if pred_fun is None:
            raise ValueError("The regressor has no predict method")
    elif is_classifier(estimator):
        if isinstance(estimator.classes_[0], np.ndarray):
            raise ValueError("Multiclass-multioutput estimators are not supported")
        reduce_binary = len(estimator.classes_) == 2
        pred_fun = getattr(estimator, "predict_proba", None)
        if pred_fun is None:
            raise ValueError("The classifier has no predict_proba method")
    else:
        raise ValueError("'estimator' must be a regressor or classifier")

    # Use check_array only on lists and other non-array-likes / sparse. Do not
    # convert DataFrame into a NumPy array.
    if not (hasattr(X, "__array__") or sparse.issparse(X)):
        X = check_array(X, force_all_finite="allow-nan", dtype=object)

    sample_weight = _check_sample_weight(sample_weight, X)

    # Usually, the data is too large and we need subsampling
    if _num_samples(X) > subsample:
        X, sample_weight = resample(
            *(X, sample_weight),
            replace=False,
            n_samples=subsample,
            random_state=random_state,
        )

    if features is None:
        features = feature_indices = np.arange(X.shape[1])
    else:
        feature_indices = np.asarray(
            _get_column_indices(X, features), dtype=np.intp, order="C"
        ).ravel()

    # CALCULATIONS
    pd_univariate = [
        _calculate_pd_over_data(
            pred_fun,
            X=X,
            feature_indices=[idx],
            sample_weight=sample_weight,
            reduce_binary=reduce_binary,
        )
        for idx in feature_indices
    ]

    n_features = len(features)
    n_pairs = int(n_features * (n_features - 1) / 2)
    output_dim = pd_univariate[0].shape[1]
    output_shape = (n_pairs, output_dim)

    num = np.empty(output_shape)
    denom = np.empty(output_shape)

    for i, (j, k) in enumerate(itertools.combinations(range(n_features), 2)):
        pd_bivariate = _calculate_pd_over_data(
            pred_fun,
            X=X,
            feature_indices=feature_indices[[j, k]],
            sample_weight=sample_weight,
            reduce_binary=reduce_binary,
        )
        num[i] = np.average(
            (pd_bivariate - pd_univariate[j] - pd_univariate[k]) ** 2,
            axis=0,
            weights=sample_weight,
        )
        denom[i] = np.average(pd_bivariate**2, axis=0, weights=sample_weight)

    num[np.abs(num) < eps] = 0  # Round small numerators to 0
    h2_stat = np.divide(num, denom, out=np.zeros_like(num), where=denom > 0)

    return Bunch(
        feature_pairs=list(itertools.combinations(features, 2)),
        h_squared_pairwise=h2_stat,
        numerator_pairwise=num,
        denominator_pairwise=denom,
    )
