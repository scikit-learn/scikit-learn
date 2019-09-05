"""
Testing for the partial dependence module.
"""

import numpy as np
import pytest

import sklearn
from sklearn.inspection import partial_dependence
from sklearn.inspection.partial_dependence import (
    _grid_from_X,
    _partial_dependence_brute,
    _partial_dependence_recursion
)
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.experimental import enable_hist_gradient_boosting  # noqa
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import MultiTaskLasso
from sklearn.tree import DecisionTreeRegressor
from sklearn.datasets import load_iris
from sklearn.datasets import make_classification, make_regression
from sklearn.cluster import KMeans
from sklearn.metrics import r2_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures
from sklearn.preprocessing import StandardScaler
from sklearn.dummy import DummyClassifier
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.utils.testing import assert_allclose
from sklearn.utils.testing import assert_array_equal


# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
y = [-1, -1, -1, 1, 1, 1]


# (X, y), n_targets  <-- as expected in the output of partial_dep()
binary_classification_data = (make_classification(n_samples=50,
                                                  random_state=0), 1)
multiclass_classification_data = (make_classification(n_samples=50,
                                                      n_classes=3,
                                                      n_clusters_per_class=1,
                                                      random_state=0), 3)
regression_data = (make_regression(n_samples=50, random_state=0), 1)
multioutput_regression_data = (make_regression(n_samples=50, n_targets=2,
                                               random_state=0), 2)


@pytest.mark.parametrize('Estimator, method, data', [
    (GradientBoostingClassifier, 'recursion', binary_classification_data),
    (GradientBoostingClassifier, 'recursion', multiclass_classification_data),
    (GradientBoostingClassifier, 'brute', binary_classification_data),
    (GradientBoostingClassifier, 'brute', multiclass_classification_data),
    (GradientBoostingRegressor, 'recursion', regression_data),
    (GradientBoostingRegressor, 'brute', regression_data),
    (DecisionTreeRegressor, 'brute', regression_data),
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
    percentiles = (.05, .95)
    grid_resolution = 100
    X = np.asarray([[1, 2],
                    [3, 4]])
    grid, axes = _grid_from_X(X, percentiles, grid_resolution)
    assert_array_equal(grid, [[1, 2],
                              [1, 4],
                              [3, 2],
                              [3, 4]])
    assert_array_equal(axes, X.T)

    # test shapes of returned objects depending on the number of unique values
    # for a feature.
    rng = np.random.RandomState(0)
    grid_resolution = 15

    # n_unique_values > grid_resolution
    X = rng.normal(size=(20, 2))
    grid, axes = _grid_from_X(X, percentiles, grid_resolution=grid_resolution)
    assert grid.shape == (grid_resolution * grid_resolution, X.shape[1])
    assert np.asarray(axes).shape == (2, grid_resolution)

    # n_unique_values < grid_resolution, will use actual values
    n_unique_values = 12
    X[n_unique_values - 1:, 0] = 12345
    rng.shuffle(X)  # just to make sure the order is irrelevant
    grid, axes = _grid_from_X(X, percentiles, grid_resolution=grid_resolution)
    assert grid.shape == (n_unique_values * grid_resolution, X.shape[1])
    # axes is a list of arrays of different shapes
    assert axes[0].shape == (n_unique_values,)
    assert axes[1].shape == (grid_resolution,)


@pytest.mark.parametrize(
    "grid_resolution, percentiles, err_msg",
    [(2, (0, 0.0001), "percentiles are too close"),
     (100, (1, 2, 3, 4), "'percentiles' must be a sequence of 2 elements"),
     (100, 12345, "'percentiles' must be a sequence of 2 elements"),
     (100, (-1, .95), r"'percentiles' values must be in \[0, 1\]"),
     (100, (.05, 2), r"'percentiles' values must be in \[0, 1\]"),
     (100, (.9, .1), r"percentiles\[0\] must be strictly less than"),
     (1, (0.05, 0.95), "'grid_resolution' must be strictly greater than 1")]
)
def test_grid_from_X_error(grid_resolution, percentiles, err_msg):
    X = np.asarray([[1, 2], [3, 4]])
    with pytest.raises(ValueError, match=err_msg):
        _grid_from_X(
            X, grid_resolution=grid_resolution, percentiles=percentiles
        )


@pytest.mark.parametrize('target_feature', range(5))
@pytest.mark.parametrize('est, method', [
    (LinearRegression(), 'brute'),
    (GradientBoostingRegressor(random_state=0), 'brute'),
    (GradientBoostingRegressor(random_state=0), 'recursion'),
    (HistGradientBoostingRegressor(random_state=0), 'brute'),
    (HistGradientBoostingRegressor(random_state=0), 'recursion')]
)
def test_partial_dependence_helpers(est, method, target_feature):
    # Check that what is returned by _partial_dependence_brute or
    # _partial_dependence_recursion is equivalent to manually setting a target
    # feature to a given value, and computing the average prediction over all
    # samples.
    # This also checks that the brute and recursion methods give the same
    # output.

    X, y = make_regression(random_state=0, n_features=5, n_informative=5)
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

    # allow for greater margin for error with recursion method
    rtol = 1e-1 if method == 'recursion' else 1e-3
    assert np.allclose(pdp, mean_predictions, rtol=rtol)


@pytest.mark.parametrize('est', (
    GradientBoostingClassifier(random_state=0),
    HistGradientBoostingClassifier(random_state=0),
))
@pytest.mark.parametrize('target_feature', (0, 1, 2, 3, 4, 5))
def test_recursion_decision_function(est, target_feature):
    # Make sure the recursion method (implicitly uses decision_function) has
    # the same result as using brute method with
    # response_method=decision_function

    X, y = make_classification(n_classes=2, n_clusters_per_class=1,
                               random_state=1)
    assert np.mean(y) == .5  # make sure the init estimator predicts 0 anyway

    est.fit(X, y)

    preds_1, _ = partial_dependence(est, X, [target_feature],
                                    response_method='decision_function',
                                    method='recursion')
    preds_2, _ = partial_dependence(est, X, [target_feature],
                                    response_method='decision_function',
                                    method='brute')

    assert_allclose(preds_1, preds_2, atol=1e-7)


@pytest.mark.parametrize('est', (
    LinearRegression(),
    GradientBoostingRegressor(random_state=0),
    HistGradientBoostingRegressor(random_state=0, min_samples_leaf=1,
                                  max_leaf_nodes=None, max_iter=1))
)
@pytest.mark.parametrize('power', (1, 2))
def test_partial_dependence_easy_target(est, power):
    # If the target y only depends on one feature in an obvious way (linear or
    # quadratic) then the partial dependence for that feature should reflect
    # it.
    # We here fit a linear regression_data model (with polynomial features if
    # needed) and compute r_squared to check that the partial dependence
    # correctly reflects the target.

    rng = np.random.RandomState(0)
    n_samples = 200
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


class NoPredictProbaNoDecisionFunction(ClassifierMixin, BaseEstimator):
    def fit(self, X, y):
        # simulate that we have some classes
        self.classes_ = [0, 1]
        return self


@pytest.mark.parametrize(
    "estimator, params, err_msg",
    [(KMeans(),
      {'features': [0]},
      "'estimator' must be a fitted regressor or classifier"),
     (LinearRegression(),
      {'features': [0], 'response_method': 'predict_proba'},
      'The response_method parameter is ignored for regressors'),
     (GradientBoostingClassifier(random_state=0),
      {'features': [0], 'response_method': 'predict_proba',
       'method': 'recursion'},
      "'recursion' method, the response_method must be 'decision_function'"),
     (GradientBoostingClassifier(random_state=0),
      {'features': [0], 'response_method': 'predict_proba', 'method': 'auto'},
      "'recursion' method, the response_method must be 'decision_function'"),
     (GradientBoostingClassifier(random_state=0),
      {'features': [0], 'response_method': 'blahblah'},
      'response_method blahblah is invalid. Accepted response_method'),
     (NoPredictProbaNoDecisionFunction(),
      {'features': [0], 'response_method': 'auto'},
      'The estimator has no predict_proba and no decision_function method'),
     (NoPredictProbaNoDecisionFunction(),
      {'features': [0], 'response_method': 'predict_proba'},
      'The estimator has no predict_proba method.'),
     (NoPredictProbaNoDecisionFunction(),
      {'features': [0], 'response_method': 'decision_function'},
      'The estimator has no decision_function method.'),
     (LinearRegression(),
      {'features': [0], 'method': 'blahblah'},
      'blahblah is invalid. Accepted method names are brute, recursion, auto'),
     (LinearRegression(),
      {'features': [0], 'method': 'recursion'},
      "Only the following estimators support the 'recursion' method:")]
)
def test_partial_dependence_error(estimator, params, err_msg):
    X, y = make_classification(random_state=0)
    estimator.fit(X, y)

    with pytest.raises(ValueError, match=err_msg):
        partial_dependence(estimator, X, **params)


@pytest.mark.parametrize(
    'estimator',
    [LinearRegression(), GradientBoostingClassifier(random_state=0)]
)
@pytest.mark.parametrize('features', [-1, 1000000])
def test_partial_dependence_unknown_feature(estimator, features):
    X, y = make_classification(random_state=0)
    estimator.fit(X, y)

    err_msg = 'all features must be in'
    with pytest.raises(ValueError, match=err_msg):
        partial_dependence(estimator, X, [features])


@pytest.mark.parametrize(
    'estimator',
    [LinearRegression(), GradientBoostingClassifier(random_state=0)]
)
def test_partial_dependence_unfitted_estimator(estimator):
    err_msg = "'estimator' parameter must be a fitted estimator"
    with pytest.raises(ValueError, match=err_msg):
        partial_dependence(estimator, X, [0])


@pytest.mark.parametrize(
    'estimator',
    [LinearRegression(), GradientBoostingClassifier(random_state=0)]
)
def test_partial_dependence_X_list(estimator):
    # check that array-like objects are accepted
    X, y = make_classification(random_state=0)
    estimator.fit(X, y)
    partial_dependence(estimator, list(X), [0])


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


def test_partial_dependence_sample_weight():
    # Test near perfect correlation between partial dependence and diagonal
    # when sample weights emphasize y = x predictions
    # non-regression test for #13193
    # TODO: extend to HistGradientBoosting once sample_weight is supported
    N = 1000
    rng = np.random.RandomState(123456)
    mask = rng.randint(2, size=N, dtype=bool)

    x = rng.rand(N)
    # set y = x on mask and y = -x outside
    y = x.copy()
    y[~mask] = -y[~mask]
    X = np.c_[mask, x]
    # sample weights to emphasize data points where y = x
    sample_weight = np.ones(N)
    sample_weight[mask] = 1000.

    clf = GradientBoostingRegressor(n_estimators=10, random_state=1)
    clf.fit(X, y, sample_weight=sample_weight)

    pdp, values = partial_dependence(clf, X, features=[1])

    assert np.corrcoef(pdp, values)[0, 1] > 0.99


def test_partial_dependence_pipeline():
    # check that the partial dependence support pipeline
    iris = load_iris()

    scaler = StandardScaler()
    clf = DummyClassifier(random_state=42)
    pipe = make_pipeline(scaler, clf)

    clf.fit(scaler.fit_transform(iris.data), iris.target)
    pipe.fit(iris.data, iris.target)

    features = 0
    pdp_pipe, values_pipe = partial_dependence(
        pipe, iris.data, features=[features]
    )
    pdp_clf, values_clf = partial_dependence(
        clf, scaler.transform(iris.data), features=[features]
    )
    assert_allclose(pdp_pipe, pdp_clf)
    assert_allclose(
        values_pipe[0],
        values_clf[0] * scaler.scale_[features] + scaler.mean_[features]
    )
