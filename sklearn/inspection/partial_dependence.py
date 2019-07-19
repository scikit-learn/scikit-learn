"""Partial dependence plots for regression and classification models."""

# Authors: Peter Prettenhofer
#          Trevor Stephens
#          Nicolas Hug
# License: BSD 3 clause

import warnings
from itertools import count
import numbers
from collections.abc import Iterable

import numpy as np
from scipy.stats.mstats import mquantiles
from joblib import Parallel, delayed

from ..base import is_classifier, is_regressor
from ..utils.extmath import cartesian
from ..utils import check_array
from ..utils import check_matplotlib_support  # noqa
from ..utils.validation import check_is_fitted
from ..tree._tree import DTYPE
from ..exceptions import NotFittedError
from ..ensemble.gradient_boosting import BaseGradientBoosting


__all__ = ['partial_dependence', 'plot_partial_dependence']


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
        uniques = np.unique(X[:, feature])
        if uniques.shape[0] < grid_resolution:
            # feature has low resolution use unique vals
            axis = uniques
        else:
            # create axis based on percentiles and grid resolution
            emp_percentiles = mquantiles(X[:, feature], prob=percentiles,
                                         axis=0)
            if np.allclose(emp_percentiles[0],
                           emp_percentiles[1]):
                raise ValueError(
                    'percentiles are too close to each other, '
                    'unable to build the grid. Please choose percentiles '
                    'that are further apart.')
            axis = np.linspace(emp_percentiles[0],
                               emp_percentiles[1],
                               num=grid_resolution, endpoint=True)
        values.append(axis)

    return cartesian(values), values


def _partial_dependence_recursion(est, grid, features):
    if est.init is not None:
        warnings.warn(
            'Using recursion method with a non-constant init predictor will '
            'lead to incorrect partial dependence values.',
            UserWarning
        )

    # grid needs to be DTYPE
    grid = np.asarray(grid, dtype=DTYPE, order='C')

    n_estimators, n_trees_per_stage = est.estimators_.shape
    averaged_predictions = np.zeros((n_trees_per_stage, grid.shape[0]),
                                    dtype=np.float64, order='C')
    for stage in range(n_estimators):
        for k in range(n_trees_per_stage):
            tree = est.estimators_[stage, k].tree_
            tree.compute_partial_dependence(grid, features,
                                            averaged_predictions[k])
    averaged_predictions *= est.learning_rate

    return averaged_predictions


def _partial_dependence_brute(est, grid, features, X, response_method):
    averaged_predictions = []

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
            X_eval[:, variable] = new_values[i]

        try:
            predictions = prediction_method(X_eval)
        except NotFittedError:
            raise ValueError(
                "'estimator' parameter must be a fitted estimator")

        # Note: predictions is of shape
        # (n_points,) for non-multioutput regressors
        # (n_points, n_tasks) for multioutput regressors
        # (n_points, 1) for the regressors in cross_decomposition (I think)
        # (n_points, 2)  for binary classifaction
        # (n_points, n_classes) for multiclass classification

        # average over samples
        averaged_predictions.append(np.mean(predictions, axis=0))

    # reshape to (n_targets, n_points) where n_targets is:
    # - 1 for non-multioutput regression and binary classification (shape is
    #   already correct in those cases)
    # - n_tasks for multi-output regression
    # - n_classes for multiclass classification.
    averaged_predictions = np.array(averaged_predictions).T
    if is_regressor(est) and averaged_predictions.ndim == 1:
        # non-multioutput regression, shape is (n_points,)
        averaged_predictions = averaged_predictions.reshape(1, -1)
    elif is_classifier(est) and averaged_predictions.shape[0] == 2:
        # Binary classification, shape is (2, n_points).
        # we output the effect of **positive** class
        averaged_predictions = averaged_predictions[1]
        averaged_predictions = averaged_predictions.reshape(1, -1)

    return averaged_predictions


def partial_dependence(estimator, X, features, response_method='auto',
                       percentiles=(0.05, 0.95), grid_resolution=100,
                       method='auto'):
    """Partial dependence of ``features``.

    Partial dependence of a feature (or a set of features) corresponds to
    the average response of an estimator for each possible value of the
    feature.

    Read more in the :ref:`User Guide <partial_dependence>`.

    Parameters
    ----------
    estimator : BaseEstimator
        A fitted estimator object implementing `predict`, `predict_proba`,
        or `decision_function`. Multioutput-multiclass classifiers are not
        supported.
    X : array-like, shape (n_samples, n_features)
        ``X`` is used both to generate a grid of values for the
        ``features``, and to compute the averaged predictions when
        method is 'brute'.
    features : list or array-like of int
        The target features for which the partial dependency should be
        computed.
    response_method : 'auto', 'predict_proba' or 'decision_function', \
            optional (default='auto')
        Specifies whether to use :term:`predict_proba` or
        :term:`decision_function` as the target response. For regressors
        this parameter is ignored and the response is always the output of
        :term:`predict`. By default, :term:`predict_proba` is tried first
        and we revert to :term:`decision_function` if it doesn't exist. If
        ``method`` is 'recursion', the response is always the output of
        :term:`decision_function`.
    percentiles : tuple of float, optional (default=(0.05, 0.95))
        The lower and upper percentile used to create the extreme values
        for the grid. Must be in [0, 1].
    grid_resolution : int, optional (default=100)
        The number of equally spaced points on the grid, for each target
        feature.
    method : str, optional (default='auto')
        The method used to calculate the averaged predictions:

        - 'recursion' is only supported for objects inheriting from
          `BaseGradientBoosting`, but is more efficient in terms of speed.
          With this method, ``X`` is only used to build the
          grid and the partial dependences are computed using the training
          data. This method does not account for the ``init`` predicor of
          the boosting process, which may lead to incorrect values (see
          warning below). With this method, the target response of a
          classifier is always the decision function, not the predicted
          probabilities.

        - 'brute' is supported for any estimator, but is more
          computationally intensive.

        - If 'auto', then 'recursion' will be used for
          ``BaseGradientBoosting`` estimators with ``init=None``, and 'brute'
          for all other.

    Returns
    -------
    averaged_predictions : ndarray, \
            shape (n_outputs, len(values[0]), len(values[1]), ...)
        The predictions for all the points in the grid, averaged over all
        samples in X (or over the training data if ``method`` is
        'recursion'). ``n_outputs`` corresponds to the number of classes in
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
    >>> partial_dependence(gb, features=[0], X=X, percentiles=(0, 1),
    ...                    grid_resolution=2) # doctest: +SKIP
    (array([[-4.52...,  4.52...]]), [array([ 0.,  1.])])

    See also
    --------
    sklearn.inspection.plot_partial_dependence: Plot partial dependence

    Warnings
    --------
    The 'recursion' method only works for gradient boosting estimators, and
    unlike the 'brute' method, it does not account for the ``init``
    predictor of the boosting process. In practice this will produce the
    same values as 'brute' up to a constant offset in the target response,
    provided that ``init`` is a consant estimator (which is the default).
    However, as soon as ``init`` is not a constant estimator, the partial
    dependence values are incorrect for 'recursion'.

    """

    if not (is_classifier(estimator) or is_regressor(estimator)):
        raise ValueError(
            "'estimator' must be a fitted regressor or classifier.")

    if is_classifier(estimator):
        if not hasattr(estimator, 'classes_'):
            raise ValueError(
                "'estimator' parameter must be a fitted estimator"
            )
        if isinstance(estimator.classes_[0], np.ndarray):
            raise ValueError(
                'Multiclass-multioutput estimators are not supported'
            )

    X = check_array(X)

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
        else:
            method = 'brute'

    if method == 'recursion':
        if not isinstance(estimator, BaseGradientBoosting):
            raise ValueError(
                "'estimator' must be an instance of BaseGradientBoosting "
                "for the 'recursion' method. Try using method='brute'.")
        if response_method == 'auto':
            response_method = 'decision_function'

        if response_method != 'decision_function':
            raise ValueError(
                "With the 'recursion' method, the response_method must be "
                "'decision_function'. Got {}.".format(response_method)
            )
        check_is_fitted(estimator, 'estimators_',
                        msg="'estimator' parameter must be a fitted estimator")
        # Note: if method is brute, this check is done at prediction time
        n_features = estimator.n_features_
    else:
        n_features = X.shape[1]

    features = np.asarray(features, dtype=np.int32, order='C').ravel()
    if any(not (0 <= f < n_features) for f in features):
        raise ValueError('all features must be in [0, %d]'
                         % (n_features - 1))

    grid, values = _grid_from_X(X[:, features], percentiles,
                                grid_resolution)
    if method == 'brute':
        averaged_predictions = _partial_dependence_brute(estimator, grid,
                                                         features, X,
                                                         response_method)
    else:
        averaged_predictions = _partial_dependence_recursion(estimator, grid,
                                                             features)

    # reshape averaged_predictions to
    # (n_outputs, n_values_feature_0, n_values_feature_1, ...)
    averaged_predictions = averaged_predictions.reshape(
        -1, *[val.shape[0] for val in values])

    return averaged_predictions, values


def plot_partial_dependence(estimator, X, features, feature_names=None,
                            target=None, response_method='auto', n_cols=3,
                            grid_resolution=100, percentiles=(0.05, 0.95),
                            method='auto', n_jobs=None, verbose=0, fig=None,
                            line_kw=None, contour_kw=None):
    """Partial dependence plots.

    The ``len(features)`` plots are arranged in a grid with ``n_cols``
    columns. Two-way partial dependence plots are plotted as contour plots.

    Read more in the :ref:`User Guide <partial_dependence>`.

    Parameters
    ----------
    estimator : BaseEstimator
        A fitted estimator object implementing `predict`, `predict_proba`,
        or `decision_function`. Multioutput-multiclass classifiers are not
        supported.
    X : array-like, shape (n_samples, n_features)
        The data to use to build the grid of values on which the dependence
        will be evaluated. This is usually the training data.
    features : list of {int, str, pair of int, pair of str}
        The target features for which to create the PDPs.
        If features[i] is an int or a string, a one-way PDP is created; if
        features[i] is a tuple, a two-way PDP is created. Each tuple must be
        of size 2.
        if any entry is a string, then it must be in ``feature_names``.
    feature_names : seq of str, shape (n_features,), optional
        Name of each feature; feature_names[i] holds the name of the feature
        with index i. By default, the name of the feature corresponds to
        their numerical index.
    target : int, optional (default=None)
        - In a multiclass setting, specifies the class for which the PDPs
          should be computed. Note that for binary classification, the
          positive class (index 1) is always used.
        - In a multioutput setting, specifies the task for which the PDPs
          should be computed
        Ignored in binary classification or classical regression settings.
    response_method : 'auto', 'predict_proba' or 'decision_function', \
            optional (default='auto') :
        Specifies whether to use :term:`predict_proba` or
        :term:`decision_function` as the target response. For regressors
        this parameter is ignored and the response is always the output of
        :term:`predict`. By default, :term:`predict_proba` is tried first
        and we revert to :term:`decision_function` if it doesn't exist. If
        ``method`` is 'recursion', the response is always the output of
        :term:`decision_function`.
    n_cols : int, optional (default=3)
        The maximum number of columns in the grid plot.
    grid_resolution : int, optional (default=100)
        The number of equally spaced points on the axes of the plots, for each
        target feature.
    percentiles : tuple of float, optional (default=(0.05, 0.95))
        The lower and upper percentile used to create the extreme values
        for the PDP axes. Must be in [0, 1].
    method : str, optional (default='auto')
        The method to use to calculate the partial dependence predictions:

        - 'recursion' is only supported for objects inheriting from
          `BaseGradientBoosting`, but is more efficient in terms of speed.
          With this method, ``X`` is optional and is only used to build the
          grid and the partial dependences are computed using the training
          data. This method does not account for the ``init`` predicor of
          the boosting process, which may lead to incorrect values (see
          warning below. With this method, the target response of a
          classifier is always the decision function, not the predicted
          probabilities.

        - 'brute' is supported for any estimator, but is more
          computationally intensive.

        - If 'auto', then 'recursion' will be used for
          ``BaseGradientBoosting`` estimators with ``init=None``, and
          'brute' for all other.

        Unlike the 'brute' method, 'recursion' does not account for the
        ``init`` predictor of the boosting process. In practice this still
        produces the same plots, up to a constant offset in the target
        response.
    n_jobs : int, optional (default=None)
        The number of CPUs to use to compute the partial dependences.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.
    verbose : int, optional (default=0)
        Verbose output during PD computations.
    fig : Matplotlib figure object, optional (default=None)
        A figure object onto which the plots will be drawn, after the figure
        has been cleared. By default, a new one is created.
    line_kw : dict, optional
        Dict with keywords passed to the ``matplotlib.pyplot.plot`` call.
        For one-way partial dependence plots.
    contour_kw : dict, optional
        Dict with keywords passed to the ``matplotlib.pyplot.plot`` call.
        For two-way partial dependence plots.

    Examples
    --------
    >>> from sklearn.datasets import make_friedman1
    >>> from sklearn.ensemble import GradientBoostingRegressor
    >>> X, y = make_friedman1()
    >>> clf = GradientBoostingRegressor(n_estimators=10).fit(X, y)
    >>> plot_partial_dependence(clf, X, [0, (0, 1)]) #doctest: +SKIP

    See also
    --------
    sklearn.inspection.partial_dependence: Return raw partial
      dependence values

    Warnings
    --------
    The 'recursion' method only works for gradient boosting estimators, and
    unlike the 'brute' method, it does not account for the ``init``
    predictor of the boosting process. In practice this will produce the
    same values as 'brute' up to a constant offset in the target response,
    provided that ``init`` is a consant estimator (which is the default).
    However, as soon as ``init`` is not a constant estimator, the partial
    dependence values are incorrect for 'recursion'.
    """
    check_matplotlib_support('plot_partial_dependence')  # noqa
    import matplotlib.pyplot as plt  # noqa
    from matplotlib import transforms  # noqa
    from matplotlib.ticker import MaxNLocator  # noqa
    from matplotlib.ticker import ScalarFormatter  # noqa

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

    X = check_array(X)
    n_features = X.shape[1]

    # convert feature_names to list
    if feature_names is None:
        # if feature_names is None, use feature indices as name
        feature_names = [str(i) for i in range(n_features)]
    elif isinstance(feature_names, np.ndarray):
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
            fxs = [convert_feature(fx) for fx in fxs]
        except TypeError:
            raise ValueError('Each entry in features must be either an int, '
                             'a string, or an iterable of size at most 2.')
        if not (1 <= np.size(fxs) <= 2):
            raise ValueError('Each entry in features must be either an int, '
                             'a string, or an iterable of size at most 2.')

        tmp_features.append(fxs)

    features = tmp_features

    names = []
    try:
        for fxs in features:
            names_ = []
            # explicit loop so "i" is bound for exception below
            for i in fxs:
                names_.append(feature_names[i])
            names.append(names_)
    except IndexError:
        raise ValueError('All entries of features must be less than '
                         'len(feature_names) = {0}, got {1}.'
                         .format(len(feature_names), i))

    # compute averaged predictions
    pd_result = Parallel(n_jobs=n_jobs, verbose=verbose)(
        delayed(partial_dependence)(estimator, X, fxs,
                                    response_method=response_method,
                                    method=method,
                                    grid_resolution=grid_resolution,
                                    percentiles=percentiles)
        for fxs in features)

    # For multioutput regression, we can only check the validity of target
    # now that we have the predictions.
    # Also note: as multiclass-multioutput classifiers are not supported,
    # multiclass and multioutput scenario are mutually exclusive. So there is
    # no risk of overwriting target_idx here.
    avg_preds, _ = pd_result[0]  # checking the first result is enough
    if is_regressor(estimator) and avg_preds.shape[0] > 1:
        if target is None:
            raise ValueError(
                'target must be specified for multi-output regressors')
        if not 0 <= target <= avg_preds.shape[0]:
            raise ValueError(
                'target must be in [0, n_tasks], got {}.'.format(target))
        target_idx = target

    # get global min and max values of PD grouped by plot type
    pdp_lim = {}
    for avg_preds, values in pd_result:
        min_pd = avg_preds[target_idx].min()
        max_pd = avg_preds[target_idx].max()
        n_fx = len(values)
        old_min_pd, old_max_pd = pdp_lim.get(n_fx, (min_pd, max_pd))
        min_pd = min(min_pd, old_min_pd)
        max_pd = max(max_pd, old_max_pd)
        pdp_lim[n_fx] = (min_pd, max_pd)

    # create contour levels for two-way plots
    if 2 in pdp_lim:
        Z_level = np.linspace(*pdp_lim[2], num=8)

    if fig is None:
        fig = plt.figure()
    else:
        fig.clear()

    if line_kw is None:
        line_kw = {'color': 'green'}
    if contour_kw is None:
        contour_kw = {}

    n_cols = min(n_cols, len(features))
    n_rows = int(np.ceil(len(features) / float(n_cols)))
    axs = []
    for i, fx, name, (avg_preds, values) in zip(
            count(), features, names, pd_result):
        ax = fig.add_subplot(n_rows, n_cols, i + 1)

        if len(values) == 1:
            ax.plot(values[0], avg_preds[target_idx].ravel(), **line_kw)
        else:
            # make contour plot
            assert len(values) == 2
            XX, YY = np.meshgrid(values[0], values[1])
            Z = avg_preds[target_idx].T
            CS = ax.contour(XX, YY, Z, levels=Z_level, linewidths=0.5,
                            colors='k')
            ax.contourf(XX, YY, Z, levels=Z_level, vmax=Z_level[-1],
                        vmin=Z_level[0], alpha=0.75, **contour_kw)
            ax.clabel(CS, fmt='%2.2f', colors='k', fontsize=10, inline=True)

        # plot data deciles + axes labels
        deciles = mquantiles(X[:, fx[0]], prob=np.arange(0.1, 1.0, 0.1))
        trans = transforms.blended_transform_factory(ax.transData,
                                                     ax.transAxes)
        ylim = ax.get_ylim()
        ax.vlines(deciles, [0], 0.05, transform=trans, color='k')
        ax.set_xlabel(name[0])
        ax.set_ylim(ylim)

        # prevent x-axis ticks from overlapping
        ax.xaxis.set_major_locator(MaxNLocator(nbins=6, prune='lower'))
        tick_formatter = ScalarFormatter()
        tick_formatter.set_powerlimits((-3, 4))
        ax.xaxis.set_major_formatter(tick_formatter)

        if len(values) > 1:
            # two-way PDP - y-axis deciles + labels
            deciles = mquantiles(X[:, fx[1]], prob=np.arange(0.1, 1.0, 0.1))
            trans = transforms.blended_transform_factory(ax.transAxes,
                                                         ax.transData)
            xlim = ax.get_xlim()
            ax.hlines(deciles, [0], 0.05, transform=trans, color='k')
            ax.set_ylabel(name[1])
            # hline erases xlim
            ax.set_xlim(xlim)
        else:
            ax.set_ylabel('Partial dependence')

        if len(values) == 1:
            ax.set_ylim(pdp_lim[1])
        axs.append(ax)

    fig.subplots_adjust(bottom=0.15, top=0.7, left=0.1, right=0.95, wspace=0.4,
                        hspace=0.3)
