import pytest
import scipy
import numpy as np
from numpy.testing import assert_array_equal

from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.feature_selection import SequentialFeatureSelector
from sklearn.datasets import make_regression
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import cross_val_score


@pytest.mark.parametrize('n_features_to_select', (0, 5, 0., -1, 1.1))
def test_bad_n_features_to_select(n_features_to_select):
    X, y = make_regression(n_features=5)
    sfs = SequentialFeatureSelector(LinearRegression(),
                                    n_features_to_select=n_features_to_select)
    with pytest.raises(ValueError, match="must be either None"):
        sfs.fit(X, y)


def test_bad_direction():
    X, y = make_regression(n_features=5)
    sfs = SequentialFeatureSelector(LinearRegression(), direction='bad')
    with pytest.raises(ValueError, match="must be either 'forward' or"):
        sfs.fit(X, y)


@pytest.mark.parametrize('direction', ('forward', 'backward'))
@pytest.mark.parametrize('n_features_to_select', (1, 5, 9, 'auto'))
def test_n_features_to_select(direction, n_features_to_select):
    # Make sure n_features_to_select is respected

    X, y = make_regression(n_features=10)
    sfs = SequentialFeatureSelector(LinearRegression(),
                                    n_features_to_select=n_features_to_select,
                                    direction=direction, cv=2)
    sfs.fit(X, y)
    if n_features_to_select == 'auto':
        n_features_to_select = 9  # n_features - 1

        assert sfs.get_support(indices=True).shape[0] <= n_features_to_select
        assert sfs.n_features_to_select_ <= n_features_to_select
        assert sfs.transform(X).shape[1] <= n_features_to_select
        assert sfs.get_support(
            indices=True).shape[0] == sfs.n_features_to_select_

    else:
        assert sfs.get_support(indices=True).shape[0] == n_features_to_select
        assert sfs.n_features_to_select_ == n_features_to_select
        assert sfs.transform(X).shape[1] == n_features_to_select


@pytest.mark.parametrize('direction', ('forward', 'backward'))
def test_stopping_criterion(direction):
    # Make sure n_features_to_select is respected

    X, y = make_regression(n_features=50, n_informative=10, random_state=0)

    tol = 1e-3

    sfs = SequentialFeatureSelector(LinearRegression(),
                                    n_features_to_select='auto',
                                    tol=tol,
                                    direction=direction, cv=2)
    sfs.fit(X, y)
    selected_X = sfs.transform(X)

    added_candidates = list(
        set(range(X.shape[1])) - set(sfs.get_support(indices=True)))
    added_X = np.hstack([
        selected_X,
        (X[:, np.random.choice(added_candidates)])[:, np.newaxis],
    ])

    removed_candidate = np.random.choice(
        list(range(sfs.n_features_to_select_)))
    removed_X = np.delete(selected_X, removed_candidate, axis=1)

    plain_cv_score = cross_val_score(
        LinearRegression(), X, y, cv=2).mean()
    sfs_cv_score = cross_val_score(
        LinearRegression(), selected_X, y, cv=2).mean()
    added_cv_score = cross_val_score(
        LinearRegression(), added_X, y, cv=2).mean()
    removed_cv_score = cross_val_score(
        LinearRegression(), removed_X, y, cv=2).mean()

    assert sfs_cv_score >= plain_cv_score

    if direction == 'forward':
        assert (sfs_cv_score - added_cv_score) <= tol
        assert (sfs_cv_score - removed_cv_score) >= tol

    else:
        assert (added_cv_score - sfs_cv_score) <= tol
        assert (removed_cv_score - sfs_cv_score) <= tol


@pytest.mark.parametrize('direction', ('forward', 'backward'))
@pytest.mark.parametrize('n_features_to_select, expected', (
    (.1, 1),
    (1., 10),
    (.5, 5),
    (None, 5),  # just to make sure .5 is equivalent to passing None
))
def test_n_features_to_select_float(direction, n_features_to_select, expected):
    # Test passing a float as n_features_to_select
    X, y = make_regression(n_features=10)
    sfs = SequentialFeatureSelector(LinearRegression(),
                                    n_features_to_select=n_features_to_select,
                                    direction=direction, cv=2)
    sfs.fit(X, y)
    assert sfs.n_features_to_select_ == expected


@pytest.mark.parametrize('seed', range(10))
@pytest.mark.parametrize('direction', ('forward', 'backward'))
@pytest.mark.parametrize('n_features_to_select, expected_selected_features', [
    (2, [0, 2]),  # f1 is dropped since it has no predictive power
    (1, [2]),  # f2 is more predictive than f0 so it's kept
])
def test_sanity(seed, direction, n_features_to_select,
                expected_selected_features):
    # Basic sanity check: 3 features, only f0 and f2 are correlated with the
    # target, f2 having a stronger correlation than f0. We expect f1 to be
    # dropped, and f2 to always be selected.

    rng = np.random.RandomState(seed)
    n_samples = 100
    X = rng.randn(n_samples, 3)
    y = 3 * X[:, 0] - 10 * X[:, 2]

    sfs = SequentialFeatureSelector(LinearRegression(),
                                    n_features_to_select=n_features_to_select,
                                    direction=direction, cv=2)
    sfs.fit(X, y)
    assert_array_equal(sfs.get_support(indices=True),
                       expected_selected_features)


def test_sparse_support():
    # Make sure sparse data is supported

    X, y = make_regression(n_features=10)
    X = scipy.sparse.csr_matrix(X)
    sfs = SequentialFeatureSelector(LinearRegression(), cv=2)
    sfs.fit(X, y)
    sfs.transform(X)


def test_nan_support():
    # Make sure nans are OK if the underlying estimator supports nans

    rng = np.random.RandomState(0)
    n_samples, n_features = 100, 10
    X, y = make_regression(n_samples, n_features, random_state=0)
    nan_mask = rng.randint(0, 2, size=(n_samples, n_features), dtype=bool)
    X[nan_mask] = np.nan
    sfs = SequentialFeatureSelector(HistGradientBoostingRegressor(), cv=2)
    sfs.fit(X, y)
    sfs.transform(X)

    with pytest.raises(ValueError, match='Input contains NaN'):
        # LinearRegression does not support nans
        SequentialFeatureSelector(LinearRegression(), cv=2).fit(X, y)


def test_pipeline_support():
    # Make sure that pipelines can be passed into SFS and that SFS can be
    # passed into a pipeline

    n_samples, n_features = 50, 3
    X, y = make_regression(n_samples, n_features, random_state=0)

    # pipeline in SFS
    pipe = make_pipeline(StandardScaler(), LinearRegression())
    sfs = SequentialFeatureSelector(pipe, cv=2)
    sfs.fit(X, y)
    sfs.transform(X)

    # SFS in pipeline
    sfs = SequentialFeatureSelector(LinearRegression(), cv=2)
    pipe = make_pipeline(StandardScaler(), sfs)
    pipe.fit(X, y)
    pipe.transform(X)
