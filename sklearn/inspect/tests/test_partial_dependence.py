"""
Testing for the partial dependence module.
"""

import numpy as np
from numpy.testing import assert_array_almost_equal
import pytest

import sklearn
from sklearn.utils.testing import assert_raises_regex
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import if_matplotlib
from sklearn.inspect import partial_dependence
from sklearn.inspect import plot_partial_dependence
from sklearn.inspect.partial_dependence import _grid_from_X
from sklearn.inspect.partial_dependence import _partial_dependence_brute
from sklearn.inspect.partial_dependence import _partial_dependence_recursion
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import MultiTaskLasso
from sklearn.datasets import load_boston, load_iris
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
    pdp, axes = partial_dependence(est, features=features,
                                   X=X, method=method,
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

    assert_raises_regex(ValueError, 'percentiles are too close',
                        _grid_from_X, X, grid_resolution=2,
                        percentiles=(0, 0.0001))

    for percentiles in ((1, 2, 3, 4), 12345):
        assert_raises_regex(ValueError, "percentiles must be a sequence",
                            _grid_from_X, X, percentiles=percentiles)

    for percentiles in ((-1, .95), (.05, 2)):
        assert_raises_regex(ValueError, "percentiles values must be in",
                            _grid_from_X, X, percentiles=percentiles)

    assert_raises_regex(ValueError,
                        r"percentiles\[0\] must be strictly less than",
                        _grid_from_X, X, percentiles=(.9, .1))

    assert_raises_regex(ValueError,
                        'grid_resolution must be strictly greater than 1.',
                        _grid_from_X, X, grid_resolution=1)


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

    preds_1, _ = partial_dependence(est, target_feature, X,
                                    response_method='decision_function',
                                    method='recursion')
    preds_2, _ = partial_dependence(est, target_feature, X,
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

    assert_raises_regex(ValueError,
                        "Multiclass-multioutput estimators are not supported",
                        partial_dependence, est, [0], X=X)


def test_partial_dependence_input():
    # Test input validation of partial_dependence.

    X, y = make_classification(random_state=0)

    lr = LinearRegression()
    lr.fit(X, y)
    gbc = GradientBoostingClassifier(random_state=0)
    gbc.fit(X, y)

    assert_raises_regex(ValueError,
                        "est must be a fitted regressor or classifier",
                        partial_dependence, KMeans(), [0], X)

    with pytest.raises(
            ValueError,
            match='The response_method parameter is ignored for regressors'):
        partial_dependence(lr, [0], X, response_method='predict_proba')

    assert_raises_regex(
        ValueError,
        "With the 'recursion' method, the response_method must be "
        "'decision_function'.",
        partial_dependence, gbc, [0], X, response_method='predict_proba',
        method='recursion')

    # for GBDTs, if users want to use predict_proba then they're forced to set
    # 'method' to brute.
    with pytest.raises(
            ValueError,
            match="With the 'recursion' method, the response_method must be "
                  "'decision_function"):
        partial_dependence(gbc, [0], X, response_method='predict_proba',
                           method='auto')

    assert_raises_regex(ValueError,
                        "response_method blahblah is invalid. "
                        "Accepted response",
                        partial_dependence, gbc, [0], X,
                        response_method='blahblah')

    class NoPredictProbaNoDecisionFunction(BaseEstimator, ClassifierMixin):
        pass
    bad_clf = NoPredictProbaNoDecisionFunction()

    assert_raises_regex(
        ValueError,
        'The estimator has no predict_proba and no decision_function method.',
        partial_dependence, bad_clf, [0], X, response_method='auto')

    assert_raises_regex(
        ValueError,
        'The estimator has no predict_proba method.',
        partial_dependence, bad_clf, [0], X, response_method='predict_proba')

    assert_raises_regex(
        ValueError,
        'The estimator has no decision_function method.',
        partial_dependence, bad_clf, [0], X,
        response_method='decision_function')

    assert_raises_regex(ValueError,
                        "method blahblah is invalid. Accepted method names "
                        "are brute, recursion, auto.",
                        partial_dependence, lr, [0], X, method='blahblah')

    assert_raises_regex(ValueError,
                        'est must be an instance of BaseGradientBoosting '
                        'for the "recursion" method',
                        partial_dependence, lr, [0], X, method='recursion')

    for feature in (-1, 1000000):
        for est in (lr, gbc):
            assert_raises_regex(ValueError,
                                "all features must be in",
                                partial_dependence, est, [feature], X=X)

    for unfitted_est in (LinearRegression(), GradientBoostingRegressor()):
        assert_raises_regex(ValueError,
                            'est parameter must be a fitted estimator',
                            partial_dependence, unfitted_est, [0], X=X)

    # check that array-like objects are accepted
    for est in (lr, gbc):
        partial_dependence(est, [0],  X=list(X))


@if_matplotlib
def test_plot_partial_dependence():
    # Test partial dependence plot function.
    boston = load_boston()
    clf = GradientBoostingRegressor(n_estimators=10, random_state=1)
    clf.fit(boston.data, boston.target)

    grid_resolution = 25
    fig, axs = plot_partial_dependence(clf, boston.data, [0, 1, (0, 1)],
                                       grid_resolution=grid_resolution,
                                       feature_names=boston.feature_names)
    assert len(axs) == 3
    assert all(ax.has_data for ax in axs)

    # check with str features and array feature names
    fig, axs = plot_partial_dependence(clf, boston.data, ['CRIM', 'ZN',
                                                          ('CRIM', 'ZN')],
                                       grid_resolution=grid_resolution,
                                       feature_names=boston.feature_names)

    assert len(axs) == 3
    assert all(ax.has_data for ax in axs)

    # check with list feature_names
    feature_names = boston.feature_names.tolist()
    fig, axs = plot_partial_dependence(clf, boston.data, ['CRIM', 'ZN',
                                                          ('CRIM', 'ZN')],
                                       grid_resolution=grid_resolution,
                                       feature_names=feature_names)
    assert len(axs) == 3
    assert all(ax.has_data for ax in axs)


@if_matplotlib
def test_plot_partial_dependence_multiclass():
    # Test partial dependence plot function on multi-class input.
    iris = load_iris()
    clf = GradientBoostingClassifier(n_estimators=10, random_state=1)
    clf.fit(iris.data, iris.target)

    grid_resolution = 25
    fig, axs = plot_partial_dependence(clf, iris.data, [0, 1],
                                       target=0,
                                       grid_resolution=grid_resolution)
    assert len(axs) == 2
    assert all(ax.has_data for ax in axs)

    # now with symbol labels
    target = iris.target_names[iris.target]
    clf = GradientBoostingClassifier(n_estimators=10, random_state=1)
    clf.fit(iris.data, target)

    grid_resolution = 25
    fig, axs = plot_partial_dependence(clf, iris.data, [0, 1],
                                       target='setosa',
                                       grid_resolution=grid_resolution)
    assert len(axs) == 2
    assert all(ax.has_data for ax in axs)


@if_matplotlib
def test_plot_partial_dependence_multioutput():
    # Test partial dependence plot function on multi-output input.
    (X, y), _ = multioutput_regression_data
    clf = LinearRegression()
    clf.fit(X, y)

    grid_resolution = 25
    fig, axs = plot_partial_dependence(clf, X, [0, 1],
                                       target=0,
                                       grid_resolution=grid_resolution)
    assert len(axs) == 2
    assert all(ax.has_data for ax in axs)

    fig, axs = plot_partial_dependence(clf, X, [0, 1],
                                       target=1,
                                       grid_resolution=grid_resolution)
    assert len(axs) == 2
    assert all(ax.has_data for ax in axs)


@if_matplotlib
@pytest.mark.filterwarnings('ignore:Default solver will be changed ')  # 0.22
@pytest.mark.filterwarnings('ignore:Default multi_class will be')  # 0.22
def test_plot_partial_dependence_input():
    X, y = make_classification(random_state=0)

    lr = LinearRegression()
    lr.fit(X, y)
    gbc = GradientBoostingClassifier(random_state=0)
    gbc.fit(X, y)

    # check target param for multiclass
    (X_m, y_m), _ = multiclass_classification_data
    lr_m = LogisticRegression()
    lr_m.fit(X_m, y_m)
    assert_raises_regex(ValueError,
                        'target must be specified for multi-class',
                        plot_partial_dependence, lr_m, X_m, [0],
                        target=None)
    for target in (-1, 100):
        assert_raises_regex(ValueError,
                            'target not in est.classes_',
                            plot_partial_dependence, lr_m, X_m, [0],
                            target=target)

    # check target param for multioutput
    (X_m, y_m), _ = multioutput_regression_data
    lr_m = LinearRegression()
    lr_m.fit(X_m, y_m)
    assert_raises_regex(ValueError,
                        'target must be specified for multi-output',
                        plot_partial_dependence, lr_m, X_m, [0],
                        target=None)
    for target in (-1, 100):
        assert_raises_regex(ValueError,
                            r'target must be in \[0, n_tasks\]',
                            plot_partial_dependence, lr_m, X_m, [0],
                            target=target)

    for feature_names in (None, ['abcd', 'def']):
        assert_raises_regex(ValueError,
                            'Feature foobar not in feature_names',
                            plot_partial_dependence, lr, X,
                            features=['foobar'],
                            feature_names=feature_names)

    for features in([(1, 2, 3)], [1, {}], [tuple()]):
        assert_raises_regex(ValueError,
                            'Each entry in features must be either an int, ',
                            plot_partial_dependence, lr, X,
                            features=features)

    assert_raises_regex(ValueError,
                        'All entries of features must be less than ',
                        plot_partial_dependence, lr, X,
                        features=[123],
                        feature_names=['blah'])

    assert_raises_regex(ValueError,
                        'feature_names should not contain duplicates',
                        plot_partial_dependence, lr, X,
                        features=[0, 1, 2],
                        feature_names=['a', 'b', 'a'])


@pytest.mark.skip('Passing non-constant init fails. Wait for PR #12436 '
                  'to be merged to un-skip this test')
def test_warning_recursion_non_constant_init():
    # make sure that passing a non-constant init parameter to a GBDT and using
    # recursion method yields a warning.

    gbc = GradientBoostingClassifier(init=DummyClassifier(), random_state=0)
    gbc.fit(X, y)

    assert_warns_message(
        UserWarning,
        'Using recursion method with a non-constant init predictor',
        plot_partial_dependence,
        gbc, X, [0], method='recursion'
    )

    assert_warns_message(
        UserWarning,
        'Using recursion method with a non-constant init predictor',
        partial_dependence,
        gbc, [0], X=X, method='recursion'
    )


@if_matplotlib
def test_plot_partial_dependence_fig():
    # Make sure fig object is correctly used if not None

    import matplotlib.pyplot as plt

    (X, y), _ = regression_data
    clf = LinearRegression()
    clf.fit(X, y)

    fig = plt.figure()
    grid_resolution = 25
    returned_fig, axs = plot_partial_dependence(
        clf, X, [0, 1], target=0, grid_resolution=grid_resolution, fig=fig)

    assert returned_fig is fig
