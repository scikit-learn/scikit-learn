"""Individual Conditional Expectation (ICE) plots for regression and
classification models."""

# Authors: Peter Prettenhofer
#          Trevor Stephens
#          Nicolas Hug
#          Madhura Jayaratne
# License: BSD 3 clause

from collections.abc import Iterable
from itertools import count
from itertools import chain
import numbers
import warnings

import numpy as np
from scipy import sparse
from scipy.stats.mstats import mquantiles
from joblib import Parallel, delayed

from sklearn.ensemble._hist_gradient_boosting.gradient_boosting import (
    BaseHistGradientBoosting)
from ..base import is_classifier
from ..base import is_regressor
from ..ensemble._gb import BaseGradientBoosting
from ..exceptions import NotFittedError
from ..pipeline import Pipeline
from ..utils import _determine_key_type
from ..utils import _get_column_indices
from ..utils import _safe_indexing
from ..utils import check_array
from ..utils import check_matplotlib_support  # noqa
from ..utils.extmath import cartesian
from ..utils.validation import check_is_fitted


__all__ = ['individual_conditional_expectation',
           'plot_individual_conditional_expectation',
           'IndividualConditionalExpectationDisplay']


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


def _get_predictions(est, grid, features, X, response_method):
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
        for i, variable in enumerate(features):
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


def _ice_brute(estimator, grid, features, X, response_method):
    predictions = _get_predictions(estimator, grid, features, X,
                                   response_method)
    instances = X.shape[0]

    # reshape to (n_targets, n_instances, n_points) where n_targets is:
    # - 1 for non-multioutput regression and binary classification (shape is
    #   already correct in those cases)
    # - n_tasks for multi-output regression
    # - n_classes for multiclass classification.
    predictions = np.array(predictions).T
    if is_regressor(estimator) and predictions.ndim == 2:
        # non-multioutput regression, shape is (n_instances, n_points,)
        predictions = predictions.reshape(instances, -1)
    elif is_classifier(estimator) and predictions.shape[0] == 2:
        # Binary classification, shape is (2, n_instances, n_points).
        # we output the effect of **positive** class
        predictions = predictions[1]
        predictions = predictions.reshape(instances, -1)

    return predictions


def individual_conditional_expectation(estimator, X, features,
                                       response_method='auto',
                                       percentiles=(0.05, 0.95),
                                       grid_resolution=100, **kwargs):
    """Individual Conditional Expectation (ICE) of ``features``.

    ICE of a feature (or a set of features) corresponds to the responses of an
    estimator for each possible value of the feature for all instance in ``X``.

    Read more in the :ref:`User Guide <individual_conditional_expectation>`.

    Parameters
    ----------
    estimator : BaseEstimator
        A fitted estimator object implementing :term:`predict`,
        :term:`predict_proba`, or :term:`decision_function`.
        Multioutput-multiclass classifiers are not supported.

    X : {array-like or dataframe} of shape (n_samples, n_features)
        ``X`` is used both to generate a grid of values for the
        ``features``, and to compute the averaged predictions.

    features : array-like of {int, str}
        The feature (e.g. `[0]`) or pair of interacting features
        (e.g. `[(0, 1)]`) for which the ICE should be computed.

    response_method : 'auto', 'predict_proba' or 'decision_function', \
            optional (default='auto')
        Specifies whether to use :term:`predict_proba` or
        :term:`decision_function` as the target response. For regressors
        this parameter is ignored and the response is always the output of
        :term:`predict`. By default, :term:`predict_proba` is tried first
        and we revert to :term:`decision_function` if it doesn't exist.

    percentiles : tuple of float, optional (default=(0.05, 0.95))
        The lower and upper percentile used to create the extreme values
        for the grid. Must be in [0, 1].

    grid_resolution : int, optional (default=100)
        The number of equally spaced points on the grid, for each target
        feature.

    Returns
    -------
    predictions : ndarray, \
            shape (n_outputs, n_instances, len(values[0]), len(values[1]),...)
        The predictions for all the points in the grid, for all
        samples in X. ``n_outputs`` corresponds to the number of classes in
        a multi-class setting, or to the number of tasks for multi-output
        regression. For classical regression and binary classification
        ``n_outputs==1``. ``n_values_feature_j`` corresponds to the size
        ``values[j]``.

    values : seq of 1d ndarrays
        The values with which the grid has been created. The generated grid
        is a cartesian product of the arrays in ``values``. ``len(values) ==
        len(features)``. The size of each array ``values[j]`` is either
        ``grid_resolution``, or the number of unique values in ``X[:, j]``,
        whichever is smaller.

    Examples
    --------
    >>> X = [[0, 0, 2], [1, 0, 0]]
    >>> y = [0, 1]
    >>> from sklearn.ensemble import GradientBoostingClassifier
    >>> gb = GradientBoostingClassifier(random_state=0).fit(X, y)
    >>> individual_conditional_expectation(gb, features=[0], X=X,
    ...       percentiles=(0, 1), grid_resolution=2) # doctest: +SKIP
    (array([[[2.19017313e-05, 5.16319395e-01],
        [4.83680605e-01, 9.99978098e-01]]]), [array([0., 1.])])
    """
    X, response_method, _ = _validate_pdp_ice_parameters(estimator, X,
                                                         features,
                                                         response_method,
                                                         'brute')

    features_indices = np.asarray(
        _get_column_indices(X, features), dtype=np.int32, order='C'
    ).ravel()

    grid, values = _grid_from_X(
        _safe_indexing(X, features_indices, axis=1), percentiles,
        grid_resolution
    )

    predictions = _ice_brute(
        estimator, grid, features_indices, X, response_method
    )

    # reshape predictions to
    # (n_outputs, n_instances, n_values_feature_0, n_values_feature_1, ...)
    predictions = predictions.reshape(
        -1, X.shape[0], *[val.shape[0] for val in values])

    return predictions, values


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


def plot_individual_conditional_expectation(estimator, X, features,
                                            feature_names=None, target=None,
                                            response_method='auto', n_cols=3,
                                            grid_resolution=100,
                                            percentiles=(0.05, 0.95),
                                            n_jobs=None, verbose=0,
                                            line_kw=None, ax=None):
    """Individual Conditional Expectation (ICE) plots.

    The ``len(features)`` plots are arranged in a grid with ``n_cols``
    columns. The deciles of the feature values will be shown with tick marks
    on the x-axes.

    .. note::

        :func:`plot_individual_conditional_expectation` does not support using
        the same axes with multiple calls. To plot the ICE curves for multiple
        estimators, please pass the axes created by the first call to the
        second call::

          >>> from sklearn.inspection import (
          ...                   plot_individual_conditional_expectation)
          >>> from sklearn.datasets import make_friedman1
          >>> from sklearn.linear_model import LinearRegression
          >>> X, y = make_friedman1()
          >>> est = LinearRegression().fit(X, y)
          >>> disp1 = plot_individual_conditional_expectation(est,
          ...                                 X)  # doctest: +SKIP
          >>> disp2 = plot_individual_conditional_expectation(est, X,
          ...                                 ax=disp1.axes_)  # doctest: +SKIP

    Read more in the :ref:`User Guide <individual_conditional_expectation>`.

    Parameters
    ----------
    estimator : BaseEstimator
        A fitted estimator object implementing :term:`predict`,
        :term:`predict_proba`, or :term:`decision_function`.
        Multioutput-multiclass classifiers are not supported.

    X : {array-like or dataframe} of shape (n_samples, n_features)
        The data to use to build the grid of values on which the individual
        conditional expectation will be evaluated.
        This is usually the training data.

    features : list of {int, str}
        The target features for which to create the ICEs.
        if any entry is a string, then it must be in ``feature_names``.

    feature_names : array-like of shape (n_features,), dtype=str, default=None
        Name of each feature; feature_names[i] holds the name of the feature
        with index i.
        By default, the name of the feature corresponds to their numerical
        index for NumPy array and their column name for pandas dataframe.

    target : int, optional (default=None)
        - In a multiclass setting, specifies the class for which the ICEs
          should be computed. Note that for binary classification, the
          positive class (index 1) is always used.
        - In a multioutput setting, specifies the task for which the ICEs
          should be computed.

        Ignored in binary classification or classical regression settings.

    response_method : 'auto', 'predict_proba' or 'decision_function', \
            optional (default='auto')
        Specifies whether to use :term:`predict_proba` or
        :term:`decision_function` as the target response. For regressors
        this parameter is ignored and the response is always the output of
        :term:`predict`. By default, :term:`predict_proba` is tried first
        and we revert to :term:`decision_function` if it doesn't exist.

    n_cols : int, optional (default=3)
        The maximum number of columns in the grid plot. Only active when `ax`
        is a single axis or `None`.

    grid_resolution : int, optional (default=100)
        The number of equally spaced points on the axes of the plots, for each
        target feature.

    percentiles : tuple of float, optional (default=(0.05, 0.95))
        The lower and upper percentile used to create the extreme values
        for the ICE axes. Must be in [0, 1].

    n_jobs : int, optional (default=None)
        The number of CPUs to use to compute the ICEs.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    verbose : int, optional (default=0)
        Verbose output during ICE computations.

    line_kw : dict, optional
        Dict with keywords passed to the ``matplotlib.pyplot.plot`` call.

    ax : Matplotlib axes or array-like of Matplotlib axes, default=None
        - If a single axis is passed in, it is treated as a bounding axes
            and a grid of ICE plots will be drawn within these bounds.
            The `n_cols` parameter controls the number of columns in the grid.
        - If an array-like of axes are passed in, the ICE plots will be drawn
            directly into these axes.
        - If `None`, a figure and a bounding axes is created and treated
            as the single axes case.

        .. versionadded:: 0.24

    Returns
    -------
    display:
    :class:`~sklearn.inspection.IndividualConditionalExpectationDisplay`

    Examples
    --------
    >>> from sklearn.datasets import make_friedman1
    >>> from sklearn.ensemble import GradientBoostingRegressor
    >>> X, y = make_friedman1()
    >>> clf = GradientBoostingRegressor(n_estimators=10).fit(X, y)
    >>> plot_individual_conditional_expectation(
    ...                         clf, X, [0, (0, 1)]) #doctest: +SKIP

    See also
    --------
    sklearn.inspection.individual_conditional_expectation: Return raw
    individual conditional expectation values
    """
    for feature in features:
        if not isinstance(feature, (numbers.Integral, str)):
            raise ValueError('Each entry in features must be either an int '
                             'or a string.')
    return _plot_pdp_ice(estimator, X, features,
                         individual_conditional_expectation,
                         IndividualConditionalExpectationDisplay,
                         feature_names=feature_names, target=target,
                         response_method=response_method, n_cols=n_cols,
                         grid_resolution=grid_resolution,
                         percentiles=percentiles, n_jobs=n_jobs,
                         verbose=verbose, line_kw=line_kw, ax=ax)


class IndividualConditionalExpectationDisplay:
    """Individual Conditional Expectation (ICE) visualization.

    It is recommended to use
    :func:`~sklearn.inspection.individual_conditional_expectation` to create a
    :class:`~sklearn.inspection.IndividualConditionalExpectation`.
    All parameters are stored as attributes.

    Read more in
    :ref:`sphx_glr_auto_examples_plot_individual_conditional_expectation_api.py`
    and the :ref:`User Guide <visualizations>`.

        .. versionadded:: 0.24

    Parameters
    ----------
    ice_results : list of (ndarray, ndarray)
        Results of
        :func:`~sklearn.inspection.individual_conditional_expectation` for
        ``features``. Each tuple corresponds to a (predictions, grid).

    features : list of (int,)
        Indices of features wrapped in tuples for a given plot.

    feature_names : list of str
        Feature names corresponding to the indices in ``features``.

    target_idx : int
        - In a multiclass setting, specifies the class for which the ICEs
          should be computed. Note that for binary classification, the
          positive class (index 1) is always used.
        - In a multioutput setting, specifies the task for which the ICEs
          should be computed.

        Ignored in binary classification or classical regression settings.

    ice_lim : dict
        Global min and max of predictions, such that all plots will have
        the same scale and y limits. `ice_lim[1]` is the global min and max for
        ICE curves.

    deciles : dict
        Deciles for feature indices in ``features``.

    Attributes
    ----------
    bounding_ax_ : matplotlib Axes or None
        If `ax` is an axes or None, the `bounding_ax_` is the axes where the
        grid of ICE plots are drawn. If `ax` is a list of axes or a numpy
        array of axes, `bounding_ax_` is None.

    axes_ : ndarray of matplotlib Axes
        If `ax` is an axes or None, `axes_[i, j]` is the axes on the i-th row
        and j-th column. If `ax` is a list of axes, `axes_[i]` is the i-th item
        in `ax`. Elements that are None corresponds to a nonexisting axes in
        that position.

    lines_ : ndarray of matplotlib Artists
        If `ax` is an axes or None, `line_[i, j, k]` is the k-th ICE curve
        on the i-th row and j-th column. If `ax` is a list of axes,
        `lines_[i, j]` is the j-th ICE curve corresponding to the i-th
        item in `ax`. Elements that are None corresponds to a nonexisting axes
        or an axes that does not include a line plot.

    figure_ : matplotlib Figure
        Figure containing ICE plots.

    """
    def __init__(self, ice_results, features, feature_names, target_idx,
                 ice_lim, deciles):
        self.ice_results = ice_results
        self.features = features
        self.feature_names = feature_names
        self.target_idx = target_idx
        self.ice_lim = ice_lim
        self.deciles = deciles

    def plot(self, ax=None, n_cols=3, line_kw=None, **kwargs):
        """Plot ICE plots.

        Parameters
        ----------
        ax : Matplotlib axes or array-like of Matplotlib axes, default=None
            - If a single axis is passed in, it is treated as a bounding axes
                and a grid of ICE plots will be drawn within these bounds.
                The `n_cols` parameter controls the number of columns in the
                grid.
            - If an array-like of axes are passed in, the ICE plots will be
                drawn directly into these axes.
            - If `None`, a figure and a bounding axes is created and treated
                as the single axes case.

        n_cols : int, default=3
            The maximum number of columns in the grid plot. Only active when
            `ax` is a single axes or `None`.

        line_kw : dict, default=None
            Dict with keywords passed to the `matplotlib.pyplot.plot` call.

        Returns
        -------
        display:
        :class:`~sklearn.inspection.IndividualConditionalExpectationDisplay`
        """

        check_matplotlib_support(
            "IndividualConditionalExpectationDisplay.plot")
        import matplotlib.pyplot as plt  # noqa
        from matplotlib import transforms  # noqa
        from matplotlib.gridspec import GridSpecFromSubplotSpec  # noqa

        if line_kw is None:
            line_kw = {}

        if ax is None:
            _, ax = plt.subplots()

        n_features = len(self.features)
        n_instances = len(self.ice_results[0][0][0])

        if isinstance(ax, plt.Axes):
            # If ax was set off, it has most likely been set to off
            # by a previous call to plot.
            if not ax.axison:
                raise ValueError("The ax was already used in another plot "
                                 "function, please set ax=display.axes_ "
                                 "instead")

            ax.set_axis_off()
            self.bounding_ax_ = ax
            self.figure_ = ax.figure

            n_cols = min(n_cols, n_features)
            n_rows = int(np.ceil(n_features / float(n_cols)))

            self.axes_ = np.empty((n_rows, n_cols), dtype=np.object)
            self.lines_ = np.empty((n_rows, n_cols, n_instances),
                                   dtype=np.object)

            axes_ravel = self.axes_.ravel()

            gs = GridSpecFromSubplotSpec(n_rows, n_cols,
                                         subplot_spec=ax.get_subplotspec())
            for i, spec in zip(range(n_features), gs):
                axes_ravel[i] = self.figure_.add_subplot(spec)

        else:  # array-like
            ax = np.asarray(ax, dtype=object)
            if ax.size != n_features:
                raise ValueError("Expected ax to have {} axes, got {}"
                                 .format(n_features, ax.size))

            if ax.ndim == 2:
                n_cols = ax.shape[1]
            else:
                n_cols = None

            self.bounding_ax_ = None
            self.figure_ = ax.ravel()[0].figure
            self.axes_ = ax
            self.lines_ = np.empty((len(ax), n_instances), dtype=np.object)

        lines_ravel = self.lines_.ravel(order='C')
        for i, axi, fx, (avg_preds, values) in zip(count(),
                                                   self.axes_.ravel(),
                                                   self.features,
                                                   self.ice_results):

            for j, ins in enumerate(avg_preds[self.target_idx]):
                lines_ravel[i * j + j] = axi.plot(values[0], ins.ravel(),
                                                  **line_kw)[0]

            trans = transforms.blended_transform_factory(axi.transData,
                                                         axi.transAxes)
            ylim = axi.get_ylim()
            axi.vlines(self.deciles[fx[0]], 0, 0.05, transform=trans,
                       color='k')
            axi.set_ylim(ylim)

            # Set xlabel if it is not already set
            if not axi.get_xlabel():
                axi.set_xlabel(self.feature_names[fx[0]])

            if n_cols is None or i % n_cols == 0:
                axi.set_ylabel('ICE')
            else:
                axi.set_yticklabels([])
            axi.set_ylim(self.ice_lim[1])
        return self
