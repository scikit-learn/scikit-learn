"""
Testing for the partial dependence module.
"""

import numpy as np
from numpy.testing import assert_array_almost_equal
import pytest

import sklearn
from sklearn.model_inspection import partial_dependence
from sklearn.model_inspection.partial_dependence import (
    _grid_from_X,
    _partial_dependence_brute,
    _partial_dependence_recursion
)
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import MultiTaskLasso
from sklearn.datasets import make_classification, make_regression
from sklearn.cluster import KMeans
from sklearn.metrics import r2_score
from sklearn.preprocessing import PolynomialFeatures
from sklearn.dummy import DummyClassifier
from sklearn.base import BaseEstimator, ClassifierMixin


# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
y = [-1, -1, -1, 1, 1, 1]


# (X, y), n_targets  <-- as expected in the output of partial_dep()
binary_classification_data = (make_classification(random_state=0), 1)
multiclass_classification_data = (make_classification(n_classes=3,
                                                      n_clusters_per_class=1,
                                                      random_state=0), 3)
regression_data = (make_regression(random_state=0), 1)
multioutput_regression_data = (make_regression(n_targets=2, random_state=0), 2)


@pytest.mark.filterwarnings('ignore:Default solver will be changed ')  # 0.22
@pytest.mark.filterwarnings('ignore:Default multi_class will be')  # 0.22
@pytest.mark.parametrize('Estimator, method, data', [
    (GradientBoostingClassifier, 'recursion', binary_classification_data),
    (GradientBoostingClassifier, 'recursion', multiclass_classification_data),
    (GradientBoostingClassifier, 'brute', binary_classification_data),
    (GradientBoostingClassifier, 'brute', multiclass_classification_data),
    (GradientBoostingRegressor, 'recursion', regression_data),
    (GradientBoostingRegressor, 'brute', regression_data),
    (LinearRegression, 'brute', regression_data),
    (LinearRegression, 'brute', multioutput_regression_data),
    (LogisticRegression, 'brute', binary_classification_data),
    (LogisticRegression, 'brute', multiclass_classification_data),
    (MultiTaskLasso, 'brute', multioutput_regression_data),
    ])
@pytest.mark.parametrize('grid_resolution', (5, 10))
@pytest.mark.parametrize('features', ([1], [1, 2]))
def test_output_shape(Estimator, method, data, grid_resolution,
                      features):
    # Check that partial_dependence has consistent output shape for different
    # kinds of estimators:
    # - classifiers with binary and multiclass settings
    # - regressors
    # - multi-task regressors

    est = Estimator()

    # n_target corresponds to the number of classes (1 for binary classif) or
    # the number of tasks / outputs in multi task settings. It's equal to 1 for
    # classical regression_data.
    (X, y), n_targets = data

    est.fit(X, y)
    pdp, axes = partial_dependence(est, X=X, features=features,
                                   method=method,
                                   grid_resolution=grid_resolution)

    expected_pdp_shape = (n_targets, *[grid_resolution
                                       for _ in range(len(features))])
    expected_axes_shape = (len(features), grid_resolution)

    assert pdp.shape == expected_pdp_shape
    assert axes is not None
    assert np.asarray(axes).shape == expected_axes_shape


def test_grid_from_X():
    # tests for _grid_from_X: sanity check for output, and for shapes.

    # Make sure that the grid is a cartesian product of the input (it will use
    # the unique values instead of the percentiles)
    X = np.asarray([[1, 2],
                    [3, 4]])
    grid, axes = _grid_from_X(X)
    assert_array_almost_equal(grid, [[1, 2],
                                     [1, 4],
                                     [3, 2],
                                     [3, 4]])
    assert_array_almost_equal(axes, X.T)

    # test shapes of returned objects depending on the number of unique values
    # for a feature.
    rng = np.random.RandomState(0)
    grid_resolution = 15

    # n_unique_values > grid_resolution
    X = rng.normal(size=(20, 2))
    grid, axes = _grid_from_X(X, grid_resolution=grid_resolution)
    assert grid.shape == (grid_resolution * grid_resolution, X.shape[1])
    assert np.asarray(axes).shape == (2, grid_resolution)

    # n_unique_values < grid_resolution, will use actual values
    n_unique_values = 12
    X[n_unique_values - 1:, 0] = 12345
    rng.shuffle(X)  # just to make sure the order is irrelevant
    grid, axes = _grid_from_X(X, grid_resolution=grid_resolution)
    assert grid.shape == (n_unique_values * grid_resolution, X.shape[1])
    # axes is a list of arrays of different shapes
    assert axes[0].shape == (n_unique_values,)
    assert axes[1].shape == (grid_resolution,)

    with pytest.raises(
            ValueError,
            match='percentiles are too close'):
        _grid_from_X(X, grid_resolution=2, percentiles=(0, 0.0001))

    for percentiles in ((1, 2, 3, 4), 12345):
        with pytest.raises(
                ValueError,
                match="percentiles must be a sequence"):
            _grid_from_X(X, percentiles=percentiles)

    for percentiles in ((-1, .95), (.05, 2)):
        with pytest.raises(
                ValueError,
                match="percentiles values must be in"):
            _grid_from_X(X, percentiles=percentiles)

    with pytest.raises(
            ValueError,
            match=r"percentiles\[0\] must be strictly less than"):
        _grid_from_X(X, percentiles=(.9, .1))

    with pytest.raises(
            ValueError,
            match='grid_resolution must be strictly greater than 1.'):
        _grid_from_X(X, grid_resolution=1)


@pytest.mark.parametrize('target_feature', (0, 3))
@pytest.mark.parametrize('est, method',
                         [(LinearRegression(), 'brute'),
                          (GradientBoostingRegressor(random_state=0),
                          'recursion'),
                          (GradientBoostingRegressor(random_state=0),
                          'brute')])
def test_partial_dependence_helpers(est, method, target_feature):
    # Check that what is returned by _partial_dependence_brute or
    # _partial_dependece_recursion is equivalent to manually setting a target
    # feature to a given value, and computing the average prediction over all
    # samples.
    # This also checks that the brute and recursion methods give the same
    # output.

    X, y = make_regression(random_state=0)
    # The 'init' estimator for GBDT (here the average prediction) isn't taken
    # into account with the recursion method, for technical reasons. We set
    # the mean to 0 to that this 'bug' doesn't have any effect.
    y = y - y.mean()
    est.fit(X, y)

    # target feature will be set to .5 and then to 123
    features = np.array([target_feature], dtype=np.int32)
    grid = np.array([[.5],
                     [123]])

    if method == 'brute':
        pdp = _partial_dependence_brute(est, grid, features, X,
                                        response_method='auto')
    else:
        pdp = _partial_dependence_recursion(est, grid, features)

    mean_predictions = []
    for val in (.5, 123):
        X_ = X.copy()
        X_[:, target_feature] = val
        mean_predictions.append(est.predict(X_).mean())

    pdp = pdp[0]  # (shape is (1, 2) so make it (2,))
    assert_array_almost_equal(pdp, mean_predictions, decimal=3)


@pytest.mark.parametrize('target_feature', (0, 1, 2, 3, 4, 5))
def test_recursion_decision_function(target_feature):
    # Make sure the recursion method (implicitely uses decision_function) has
    # the same result as using brute method with response_method=decision

    X, y = make_classification(n_classes=2, n_clusters_per_class=1,
                               random_state=1)
    assert np.mean(y) == .5  # make sure the init estimator predicts 0 anyway

    est = GradientBoostingClassifier(random_state=0, loss='deviance')
    est.fit(X, y)

    preds_1, _ = partial_dependence(est, X, [target_feature],
                                    response_method='decision_function',
                                    method='recursion')
    preds_2, _ = partial_dependence(est, X, [target_feature],
                                    response_method='decision_function',
                                    method='brute')

    assert_array_almost_equal(preds_1, preds_2, decimal=5)


@pytest.mark.parametrize('est', (LinearRegression(),
                                 GradientBoostingRegressor(random_state=0)))
@pytest.mark.parametrize('power', (1, 2))
def test_partial_dependence_easy_target(est, power):
    # If the target y only depends on one feature in an obvious way (linear or
    # quadratic) then the partial dependence for that feature should reflect
    # it.
    # We here fit a linear regression_data model (with polynomial features if
    # needed) and compute r_squared to check that the partial dependence
    # correctly reflects the target.

    rng = np.random.RandomState(0)
    n_samples = 100
    target_variable = 2
    X = rng.normal(size=(n_samples, 5))
    y = X[:, target_variable]**power

    est.fit(X, y)

    averaged_predictions, values = partial_dependence(
        est, features=[target_variable], X=X, grid_resolution=1000)

    new_X = values[0].reshape(-1, 1)
    new_y = averaged_predictions[0]
    # add polynomial features if needed
    new_X = PolynomialFeatures(degree=power).fit_transform(new_X)

    lr = LinearRegression().fit(new_X, new_y)
    r2 = r2_score(new_y, lr.predict(new_X))

    assert r2 > .99


@pytest.mark.filterwarnings('ignore:The default value of ')  # 0.22
@pytest.mark.parametrize('Estimator',
                         (sklearn.tree.DecisionTreeClassifier,
                          sklearn.tree.ExtraTreeClassifier,
                          sklearn.ensemble.ExtraTreesClassifier,
                          sklearn.neighbors.KNeighborsClassifier,
                          sklearn.neighbors.RadiusNeighborsClassifier,
                          sklearn.ensemble.RandomForestClassifier))
def test_multiclass_multioutput(Estimator):
    # Make sure error is raised for multiclass-multioutput classifiers

    # make multiclass-multioutput dataset
    X, y = make_classification(n_classes=3, n_clusters_per_class=1,
                               random_state=0)
    y = np.array([y, y]).T

    est = Estimator()
    est.fit(X, y)

    with pytest.raises(
            ValueError,
            match="Multiclass-multioutput estimators are not supported"):
        partial_dependence(est, X, [0])


def test_partial_dependence_input():
    # Test input validation of partial_dependence.

    X, y = make_classification(random_state=0)

    lr = LinearRegression()
    lr.fit(X, y)
    gbc = GradientBoostingClassifier(random_state=0)
    gbc.fit(X, y)

    with pytest.raises(
            ValueError,
            match="est must be a fitted regressor or classifier"):
        partial_dependence(KMeans(), X, [0])

    with pytest.raises(
            ValueError,
            match='The response_method parameter is ignored for regressors'):
        partial_dependence(lr, X, [0], response_method='predict_proba')

    with pytest.raises(
            ValueError,
            match="With the 'recursion' method, the response_method must be "
                  "'decision_function'."):
        partial_dependence(gbc, X, [0], response_method='predict_proba',
                           method='recursion')

    # for GBDTs, if users want to use predict_proba then they're forced to set
    # 'method' to brute.
    with pytest.raises(
            ValueError,
            match="With the 'recursion' method, the response_method must be "
                  "'decision_function"):
        partial_dependence(gbc, X, [0], response_method='predict_proba',
                           method='auto')

    with pytest.raises(
            ValueError,
            match="response_method blahblah is invalid. Accepted response"):
        partial_dependence(gbc, X, [0], response_method='blahblah')

    class NoPredictProbaNoDecisionFunction(BaseEstimator, ClassifierMixin):
        pass
    bad_clf = NoPredictProbaNoDecisionFunction()

    with pytest.raises(
            ValueError,
            match='The estimator has no predict_proba and no '
                  'decision_function method.'):
        partial_dependence(bad_clf, X, [0], response_method='auto')

    with pytest.raises(
            ValueError,
            match='The estimator has no predict_proba method.'):
        partial_dependence(bad_clf, X, [0], response_method='predict_proba')

    with pytest.raises(
            ValueError,
            match='The estimator has no decision_function method.'):
        partial_dependence(bad_clf, X, [0],
                           response_method='decision_function')

    with pytest.raises(
            ValueError,
            match="method blahblah is invalid. Accepted method names "
                  "are brute, recursion, auto."):
        partial_dependence(lr, X, [0], method='blahblah')

    with pytest.raises(
            ValueError,
            match='est must be an instance of BaseGradientBoosting '
                  'for the "recursion" method'):
        partial_dependence(lr, X, [0], method='recursion')

    for feature in (-1, 1000000):
        for est in (lr, gbc):
            with pytest.raises(
                    ValueError,
                    match="all features must be in"):
                partial_dependence(est, X, [feature])

    for unfitted_est in (LinearRegression(), GradientBoostingRegressor()):
        with pytest.raises(
                ValueError,
                match='est parameter must be a fitted estimator'):
            partial_dependence(unfitted_est, X, [0])

    # check that array-like objects are accepted
    for est in (lr, gbc):
        partial_dependence(est, list(X), [0])


def test_warning_recursion_non_constant_init():
    # make sure that passing a non-constant init parameter to a GBDT and using
    # recursion method yields a warning.

    gbc = GradientBoostingClassifier(init=DummyClassifier(), random_state=0)
    gbc.fit(X, y)

    with pytest.warns(
            UserWarning,
            match='Using recursion method with a non-constant init predictor'):
        partial_dependence(gbc, X, [0], method='recursion')

    with pytest.warns(
            UserWarning,
            match='Using recursion method with a non-constant init predictor'):
        partial_dependence(gbc, X, [0], method='recursion')
