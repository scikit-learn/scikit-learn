"""Partial dependence plots for regression and classification models."""

# Authors: Peter Prettenhofer
#          Trevor Stephens
#          Nicolas Hug
# License: BSD 3 clause

import warnings

import numpy as np
from scipy.stats.mstats import mquantiles

from ..base import is_classifier, is_regressor
from ..utils.extmath import cartesian
from ..utils import check_array
from ..utils.validation import check_is_fitted
from ..tree._tree import DTYPE
from ..exceptions import NotFittedError
from ..ensemble.gradient_boosting import BaseGradientBoosting
from ..ensemble._gradient_boosting import _partial_dependence_tree


__all__ = ['partial_dependence']


def _grid_from_X(X, percentiles=(0.05, 0.95), grid_resolution=100):
    """Generate a grid of points based on the percentiles of X.

    The grid is a cartesian product between the columns of ``values``. The
    ith column of ``values`` consists in ``grid_resolution`` equally-spaced
    points between the percentiles of the jth column of X.
    If ``grid_resolution`` is bigger than the number of unique values in the
    jth column of X, then those unique values will be used instead.

    Parameters
    ----------
    X : ndarray, shape=(n_samples, n_target_features)
        The data
    percentiles : tuple of floats
        The percentiles which are used to construct the extreme values of
        the grid. Must be in [0, 1].
    grid_resolution : int
        The number of equally spaced points to be placed on the grid for each
        feature.

    Returns
    -------
    grid : ndarray, shape=(n_points, X.shape[1])
        A value for each feature at each point in the grid. ``n_points`` is
        always ``<= grid_resolution ** X.shape[1]``.
    values : list of 1d ndarrays
        The values with which the grid has been created. The size of each
        array ``values[j]`` is either ``grid_resolution``, or the number of
        unique values in ``X[:, j]``, whichever is smaller.
    """
    try:
        assert len(percentiles) == 2
    except (AssertionError, TypeError):
        raise ValueError('percentiles must be a sequence of 2 elements.')
    if not all(0. <= x <= 1. for x in percentiles):
        raise ValueError('percentiles values must be in [0, 1].')
    if percentiles[0] >= percentiles[1]:
        raise ValueError('percentiles[0] must be strictly less '
                         'than percentiles[1].')

    if grid_resolution <= 1:
        raise ValueError('grid_resolution must be strictly greater than 1.')

    values = []
    for feature in range(X.shape[1]):
        uniques = np.unique(X[:, feature])
        if uniques.shape[0] < grid_resolution:
            # feature has low resolution use unique vals
            axis = uniques
        else:
            # create axis based on percentiles and grid resolution
            emp_percentiles = mquantiles(X, prob=percentiles, axis=0)
            if np.allclose(emp_percentiles[0, feature],
                           emp_percentiles[1, feature]):
                raise ValueError('percentiles are too close to each other, '
                                 'unable to build the grid.')
            axis = np.linspace(emp_percentiles[0, feature],
                               emp_percentiles[1, feature],
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

    n_trees_per_stage = est.estimators_.shape[1]
    n_estimators = est.estimators_.shape[0]
    learning_rate = est.learning_rate
    averaged_predictions = np.zeros((n_trees_per_stage, grid.shape[0]),
                                    dtype=np.float64, order='C')
    for stage in range(n_estimators):
        for k in range(n_trees_per_stage):
            tree = est.estimators_[stage, k].tree_
            _partial_dependence_tree(tree, grid, features,
                                     learning_rate, averaged_predictions[k])

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
            raise ValueError('est parameter must be a fitted estimator')

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
    X : array-like, shape=(n_samples, n_features)
        ``X`` is used both to generate a grid of values for the
        ``features``, and to compute the averaged predictions when
        method is 'brute'.
    features : list or array-like of int
        The target features for which the partial dependency should be
        computed.
    response_method : 'auto', 'predict_proba' or 'decision_function', \
            optional (default='auto') :
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
          :ref:`this warning<warning_recursion_init_plot>`). With this
          method, the target response of a classifier is always the decision
          function, not the predicted probabilities.

        - 'brute' is supported for any estimator, but is more
          computationally intensive.

        - If 'auto', then 'recursion' will be used for
          ``BaseGradientBoosting`` estimators with ``init=None``, and 'brute'
          for all other.

    Returns
    -------
    averaged_predictions : array, \
            shape=(n_outputs, len(values[0]), len(values[1]), ...)
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
    sklearn.plot.plot_partial_dependence: Plot partial dependence

    .. _warning_recursion_init:

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
        raise ValueError('est must be a fitted regressor or classifier.')

    if (hasattr(estimator, 'classes_') and
            isinstance(estimator.classes_[0], np.ndarray)):
        raise ValueError('Multiclass-multioutput estimators are not supported')

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
        if isinstance(estimator, BaseGradientBoosting) and estimator.init is None:
            method = 'recursion'
        else:
            method = 'brute'

    if method == 'recursion':
        if not isinstance(estimator, BaseGradientBoosting):
            raise ValueError(
                'est must be an instance of BaseGradientBoosting '
                'for the "recursion" method. Try using method="brute".')
        if response_method == 'auto':
            response_method = 'decision_function'

        if response_method != 'decision_function':
            raise ValueError(
                "With the 'recursion' method, the response_method must be "
                "'decision_function'. Got {}.".format(response_method)
            )
        check_is_fitted(estimator, 'estimators_',
                        msg='est parameter must be a fitted estimator')
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
    # (n_outputs, n_values_feature_0, # n_values_feature_1, ...)
    averaged_predictions = averaged_predictions.reshape(
        -1, *[val.shape[0] for val in values])

    return averaged_predictions, values
