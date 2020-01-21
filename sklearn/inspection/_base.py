import numbers
import warnings
from collections.abc import Iterable
from itertools import chain

import numpy as np
from joblib import Parallel, delayed
from scipy import sparse
from scipy.stats.mstats import mquantiles

from sklearn.ensemble._hist_gradient_boosting.gradient_boosting import (
    BaseHistGradientBoosting)
from ..base import is_classifier
from ..base import is_regressor
from ..ensemble._gb import BaseGradientBoosting
from ..exceptions import NotFittedError
from ..pipeline import Pipeline
from ..utils import _determine_key_type
from ..utils import _safe_indexing
from ..utils import check_array
from ..utils import check_matplotlib_support  # noqa
from ..utils.extmath import cartesian
from ..utils.validation import check_is_fitted


def _grid_from_X(X, percentiles, grid_resolution):
    """Generate a grid of points based on the percentiles of X.

    The grid is a cartesian product between the columns of ``values``. The
    ith column of ``values`` consists in ``grid_resolution`` equally-spaced
    points between the percentiles of the jth column of X.
    If ``grid_resolution`` is bigger than the number of unique values in the
    jth column of X, then those unique values will be used instead.

    Parameters
    ----------
    X : ndarray, shape (n_samples, n_target_features)
        The data

    percentiles : tuple of floats
        The percentiles which are used to construct the extreme values of
        the grid. Must be in [0, 1].

    grid_resolution : int
        The number of equally spaced points to be placed on the grid for each
        feature.

    Returns
    -------
    grid : ndarray, shape (n_points, n_target_features)
        A value for each feature at each point in the grid. ``n_points`` is
        always ``<= grid_resolution ** X.shape[1]``.

    values : list of 1d ndarrays
        The values with which the grid has been created. The size of each
        array ``values[j]`` is either ``grid_resolution``, or the number of
        unique values in ``X[:, j]``, whichever is smaller.
    """
    if not isinstance(percentiles, Iterable) or len(percentiles) != 2:
        raise ValueError("'percentiles' must be a sequence of 2 elements.")
    if not all(0 <= x <= 1 for x in percentiles):
        raise ValueError("'percentiles' values must be in [0, 1].")
    if percentiles[0] >= percentiles[1]:
        raise ValueError('percentiles[0] must be strictly less '
                         'than percentiles[1].')

    if grid_resolution <= 1:
        raise ValueError("'grid_resolution' must be strictly greater than 1.")

    values = []
    for feature in range(X.shape[1]):
        uniques = np.unique(_safe_indexing(X, feature, axis=1))
        if uniques.shape[0] < grid_resolution:
            # feature has low resolution use unique vals
            axis = uniques
        else:
            # create axis based on percentiles and grid resolution
            emp_percentiles = mquantiles(
                _safe_indexing(X, feature, axis=1), prob=percentiles, axis=0
            )
            if np.allclose(emp_percentiles[0], emp_percentiles[1]):
                raise ValueError(
                    'percentiles are too close to each other, '
                    'unable to build the grid. Please choose percentiles '
                    'that are further apart.')
            axis = np.linspace(emp_percentiles[0],
                               emp_percentiles[1],
                               num=grid_resolution, endpoint=True)
        values.append(axis)

    return cartesian(values), values


def _get_predictions(est, grid, features_indices, X, response_method):
    predictions = []

    # define the prediction_method (predict, predict_proba, decision_function).
    if is_regressor(est):
        prediction_method = est.predict
    else:
        predict_proba = getattr(est, 'predict_proba', None)
        decision_function = getattr(est, 'decision_function', None)
        if response_method == 'auto':
            # try predict_proba, then decision_function if it doesn't exist
            prediction_method = predict_proba or decision_function
        else:
            prediction_method = (predict_proba if response_method ==
                                 'predict_proba' else decision_function)
        if prediction_method is None:
            if response_method == 'auto':
                raise ValueError(
                    'The estimator has no predict_proba and no '
                    'decision_function method.'
                )
            elif response_method == 'predict_proba':
                raise ValueError('The estimator has no predict_proba method.')
            else:
                raise ValueError(
                    'The estimator has no decision_function method.')

    for new_values in grid:
        X_eval = X.copy()
        for i, variable in enumerate(features_indices):
            if hasattr(X_eval, 'iloc'):
                X_eval.iloc[:, variable] = new_values[i]
            else:
                X_eval[:, variable] = new_values[i]

        try:
            # Note: predictions is of shape
            # (n_points,) for non-multioutput regressors
            # (n_points, n_tasks) for multioutput regressors
            # (n_points, 1) for the regressors in cross_decomposition (I think)
            # (n_points, 2) for binary classification
            # (n_points, n_classes) for multiclass classification
            predictions.append(prediction_method(X_eval))
        except NotFittedError:
            raise ValueError(
                "'estimator' parameter must be a fitted estimator")

    return predictions


def _validate_pdp_ice_parameters(estimator, X, features, response_method,
                                 method):
    if not (is_classifier(estimator) or is_regressor(estimator)):
        raise ValueError(
            "'estimator' must be a fitted regressor or classifier."
        )

    if isinstance(estimator, Pipeline):
        # TODO: to be removed if/when pipeline get a `steps_` attributes
        # assuming Pipeline is the only estimator that does not store a new
        # attribute
        for est in estimator:
            # FIXME: remove the None option when it will be deprecated
            if est not in (None, 'drop'):
                check_is_fitted(est)
    else:
        check_is_fitted(estimator)

    if (is_classifier(estimator) and
            isinstance(estimator.classes_[0], np.ndarray)):
        raise ValueError(
            'Multiclass-multioutput estimators are not supported'
        )

    # Use check_array only on lists and other non-array-likes / sparse. Do not
    # convert DataFrame into a NumPy array.
    if not(hasattr(X, '__array__') or sparse.issparse(X)):
        X = check_array(X, force_all_finite='allow-nan', dtype=np.object)

    accepted_responses = ('auto', 'predict_proba', 'decision_function')
    if response_method not in accepted_responses:
        raise ValueError(
            'response_method {} is invalid. Accepted response_method names '
            'are {}.'.format(response_method, ', '.join(accepted_responses)))

    if is_regressor(estimator) and response_method != 'auto':
        raise ValueError(
            "The response_method parameter is ignored for regressors and "
            "must be 'auto'."
        )

    accepted_methods = ('brute', 'recursion', 'auto')
    if method not in accepted_methods:
        raise ValueError(
            'method {} is invalid. Accepted method names are {}.'.format(
                method, ', '.join(accepted_methods)))

    if method == 'auto':
        if (isinstance(estimator, BaseGradientBoosting) and
                estimator.init is None):
            method = 'recursion'
        elif isinstance(estimator, BaseHistGradientBoosting):
            method = 'recursion'
        else:
            method = 'brute'

    if method == 'recursion':
        if not isinstance(estimator,
                          (BaseGradientBoosting, BaseHistGradientBoosting)):
            supported_classes_recursion = (
                'GradientBoostingClassifier',
                'GradientBoostingRegressor',
                'HistGradientBoostingClassifier',
                'HistGradientBoostingRegressor',
            )
            raise ValueError(
                "Only the following estimators support the 'recursion' "
                "method: {}. Try using method='brute'."
                .format(', '.join(supported_classes_recursion)))
        if response_method == 'auto':
            response_method = 'decision_function'

        if response_method != 'decision_function':
            raise ValueError(
                "With the 'recursion' method, the response_method must be "
                "'decision_function'. Got {}.".format(response_method)
            )

    if _determine_key_type(features, accept_slice=False) == 'int':
        # _get_column_indices() supports negative indexing. Here, we limit
        # the indexing to be positive. The upper bound will be checked
        # by _get_column_indices()
        if np.any(np.less(features, 0)):
            raise ValueError(
                'all features must be in [0, {}]'.format(X.shape[1] - 1)
            )

    return X, response_method, method


def _plot_pdp_ice(estimator, X, features, calc_method, display_class,
                  feature_names=None, target=None, response_method='auto',
                  n_cols=3, grid_resolution=100, percentiles=(0.05, 0.95),
                  method='auto', n_jobs=None, verbose=0, fig=None,
                  line_kw=None, contour_kw=None, ax=None):

    check_matplotlib_support('_plot_pdp_ice')  # noqa
    import matplotlib.pyplot as plt  # noqa

    # set target_idx for multi-class estimators
    if hasattr(estimator, 'classes_') and np.size(estimator.classes_) > 2:
        if target is None:
            raise ValueError('target must be specified for multi-class')
        target_idx = np.searchsorted(estimator.classes_, target)
        if (not (0 <= target_idx < len(estimator.classes_)) or
                estimator.classes_[target_idx] != target):
            raise ValueError('target not in est.classes_, got {}'.format(
                target))
    else:
        # regression and binary classification
        target_idx = 0

    # Use check_array only on lists and other non-array-likes / sparse. Do not
    # convert DataFrame into a NumPy array.
    if not(hasattr(X, '__array__') or sparse.issparse(X)):
        X = check_array(X, force_all_finite='allow-nan', dtype=np.object)
    n_features = X.shape[1]

    # convert feature_names to list
    if feature_names is None:
        if hasattr(X, "loc"):
            # get the column names for a pandas dataframe
            feature_names = X.columns.tolist()
        else:
            # define a list of numbered indices for a numpy array
            feature_names = [str(i) for i in range(n_features)]
    elif hasattr(feature_names, "tolist"):
        # convert numpy array or pandas index to a list
        feature_names = feature_names.tolist()
    if len(set(feature_names)) != len(feature_names):
        raise ValueError('feature_names should not contain duplicates.')

    def convert_feature(fx):
        if isinstance(fx, str):
            try:
                fx = feature_names.index(fx)
            except ValueError:
                raise ValueError('Feature %s not in feature_names' % fx)
        return int(fx)

    # convert features into a seq of int tuples
    tmp_features = []
    for fxs in features:
        if isinstance(fxs, (numbers.Integral, str)):
            fxs = (fxs,)
        try:
            fxs = tuple(convert_feature(fx) for fx in fxs)
        except TypeError:
            raise ValueError('Each entry in features must be either an int, '
                             'a string, or an iterable of size at most 2.')
        if not 1 <= np.size(fxs) <= 2:
            raise ValueError('Each entry in features must be either an int, '
                             'a string, or an iterable of size at most 2.')

        tmp_features.append(fxs)

    features = tmp_features

    # Early exit if the axes does not have the correct number of axes
    if ax is not None and not isinstance(ax, plt.Axes):
        axes = np.asarray(ax, dtype=object)
        if axes.size != len(features):
            raise ValueError("Expected ax to have {} axes, got {}".format(
                             len(features), axes.size))

    for i in chain.from_iterable(features):
        if i >= len(feature_names):
            raise ValueError('All entries of features must be less than '
                             'len(feature_names) = {0}, got {1}.'
                             .format(len(feature_names), i))

    # compute averaged predictions
    pd_results = Parallel(n_jobs=n_jobs, verbose=verbose)(
        delayed(calc_method)(estimator, X, fxs,
                             response_method=response_method, method=method,
                             grid_resolution=grid_resolution,
                             percentiles=percentiles)
        for fxs in features)

    # For multioutput regression, we can only check the validity of target
    # now that we have the predictions.
    # Also note: as multiclass-multioutput classifiers are not supported,
    # multiclass and multioutput scenario are mutually exclusive. So there is
    # no risk of overwriting target_idx here.
    avg_preds, _ = pd_results[0]  # checking the first result is enough
    if is_regressor(estimator) and avg_preds.shape[0] > 1:
        if target is None:
            raise ValueError(
                'target must be specified for multi-output regressors')
        if not 0 <= target <= avg_preds.shape[0]:
            raise ValueError(
                'target must be in [0, n_tasks], got {}.'.format(target))
        target_idx = target

    # get global min and max average predictions of PD grouped by plot type
    pdp_lim = {}
    for avg_preds, values in pd_results:
        min_pd = avg_preds[target_idx].min()
        max_pd = avg_preds[target_idx].max()
        n_fx = len(values)
        old_min_pd, old_max_pd = pdp_lim.get(n_fx, (min_pd, max_pd))
        min_pd = min(min_pd, old_min_pd)
        max_pd = max(max_pd, old_max_pd)
        pdp_lim[n_fx] = (min_pd, max_pd)

    deciles = {}
    for fx in chain.from_iterable(features):
        if fx not in deciles:
            X_col = _safe_indexing(X, fx, axis=1)
            deciles[fx] = mquantiles(X_col, prob=np.arange(0.1, 1.0, 0.1))

    if fig is not None:
        warnings.warn("The fig parameter is deprecated in version "
                      "0.22 and will be removed in version 0.24",
                      FutureWarning)
        fig.clear()
        ax = fig.gca()

    display = display_class(pd_results, features, feature_names, target_idx,
                            pdp_lim, deciles)
    return display.plot(ax=ax, n_cols=n_cols, line_kw=line_kw,
                        contour_kw=contour_kw)
