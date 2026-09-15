import numpy as np
import pytest

from sklearn.datasets import make_classification, make_regression
from sklearn.exceptions import NotFittedError
from sklearn.inspection import accumulated_local_effects
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.utils import Bunch


def test_ale_linear_regression_correlated_features():
    """Verify that ALE recovers true linear coefficients despite high feature correlation."""
    rng = np.random.RandomState(42)
    n = 1500
    # Create features with 0.95 correlation
    cov = [[1.0, 0.95], [0.95, 1.0]]
    X = rng.multivariate_normal([0, 0], cov, size=n)
    # True weights: 2.0 for X0, 3.0 for X1
    y = 2.0 * X[:, 0] + 3.0 * X[:, 1] + rng.normal(0, 0.05, size=n)

    model = LinearRegression().fit(X, y)
    res = accumulated_local_effects(model, X, features=[0, 1], grid_resolution=25)

    assert isinstance(res, Bunch)
    assert len(res.ale) == 2
    assert len(res["values"]) == 2
    assert len(res.grid_values) == 2
    assert res.feature_names == ["x0", "x1"]

    # Compute slope of ALE curves
    slope_0, _ = np.polyfit(res.grid_values[0], res.ale[0][0], 1)
    slope_1, _ = np.polyfit(res.grid_values[1], res.ale[1][0], 1)

    np.testing.assert_allclose(slope_0, 2.0, rtol=0.08)
    np.testing.assert_allclose(slope_1, 3.0, rtol=0.08)


def test_ale_centering_property():
    """Verify that ALE is centered around 0 across the sample distribution."""
    rng = np.random.RandomState(42)
    X, y = make_regression(n_samples=500, n_features=4, random_state=rng)
    model = RandomForestRegressor(n_estimators=10, random_state=rng).fit(X, y)

    res = accumulated_local_effects(model, X, features=[0, 1, 2], grid_resolution=20)
    for i in range(3):
        # The grid values should encompass the feature min and max
        assert res.grid_values[i].min() <= X[:, i].min() + 1e-5
        assert res.grid_values[i].max() >= X[:, i].max() - 1e-5
        # The mean of grid ALE should be close to 0
        assert np.abs(np.mean(res.ale[i][0])) < 0.5


def test_ale_binary_classification():
    """Verify ALE works on binary classifiers with probability and decision_function outputs."""
    rng = np.random.RandomState(42)
    X, y = make_classification(n_samples=300, n_features=4, n_informative=2, n_redundant=1, random_state=rng)
    clf = LogisticRegression().fit(X, y)

    # Auto response_method defaults to predict_proba
    res_proba = accumulated_local_effects(clf, X, features=[0], response_method="auto")
    assert res_proba.ale[0].shape[0] == 1  # positive class probability
    assert np.all(res_proba.ale[0] >= -1.0)
    assert np.all(res_proba.ale[0] <= 1.0)

    # Explicit decision_function
    res_df = accumulated_local_effects(clf, X, features=[0], response_method="decision_function")
    assert res_df.ale[0].shape[0] == 1


def test_ale_multiclass_classification():
    """Verify ALE outputs separate effects per class for multiclass classification."""
    rng = np.random.RandomState(42)
    X, y = make_classification(
        n_samples=300, n_features=6, n_informative=3, n_redundant=1, n_classes=3, n_clusters_per_class=1, random_state=rng
    )
    clf = RandomForestClassifier(n_estimators=10, random_state=rng).fit(X, y)

    res = accumulated_local_effects(clf, X, features=[0, 1])
    assert len(res.ale) == 2
    # For 3 classes, each feature ALE has shape (3, n_grid_points)
    assert res.ale[0].shape[0] == 3
    assert res.ale[1].shape[0] == 3


def test_ale_multi_target_regression():
    """Verify ALE works for multi-output regression."""
    rng = np.random.RandomState(42)
    X = rng.normal(size=(200, 3))
    y = np.column_stack([X[:, 0] + X[:, 1], 2 * X[:, 1] - X[:, 2]])

    model = LinearRegression().fit(X, y)
    res = accumulated_local_effects(model, X, features=[1])
    assert res.ale[0].shape[0] == 2  # 2 targets


def test_ale_pandas_dataframe():
    """Verify ALE properly handles pandas DataFrames with string column names."""
    pd = pytest.importorskip("pandas")
    rng = np.random.RandomState(42)
    df = pd.DataFrame({
        "feature_a": rng.normal(size=200),
        "feature_b": rng.uniform(0, 10, size=200),
        "feature_c": rng.exponential(1.0, size=200),
    })
    y = df["feature_a"] * 2 + df["feature_b"]

    model = LinearRegression().fit(df, y)
    res = accumulated_local_effects(model, df, features=["feature_a", "feature_b"])

    assert res.feature_names == ["feature_a", "feature_b"]
    assert len(res.ale) == 2


def test_ale_sample_weight():
    """Verify sample_weight parameter operates correctly without errors."""
    rng = np.random.RandomState(42)
    X, y = make_regression(n_samples=200, n_features=2, random_state=rng)
    sw = rng.uniform(0.5, 2.0, size=200)

    model = LinearRegression().fit(X, y, sample_weight=sw)
    res = accumulated_local_effects(model, X, features=[0], sample_weight=sw)
    assert res.ale[0].shape[0] == 1


def test_ale_constant_feature_warning():
    """Verify that constant feature triggers a warning and returns zero ALE."""
    X = np.ones((100, 2))
    X[:, 1] = np.linspace(0, 10, 100)
    y = X[:, 1] * 2

    model = LinearRegression().fit(X, y)
    with pytest.warns(UserWarning, match="too few unique values"):
        res = accumulated_local_effects(model, X, features=[0])
    np.testing.assert_array_equal(res.ale[0], 0.0)


def test_ale_unfitted_estimator():
    """Verify that an unfitted estimator raises NotFittedError."""
    model = LinearRegression()
    X = np.ones((50, 2))
    with pytest.raises(NotFittedError):
        accumulated_local_effects(model, X, features=[0])


def test_ale_invalid_percentiles():
    """Verify invalid percentiles raise ValueError."""
    model = LinearRegression().fit(np.ones((50, 2)), np.ones(50))
    with pytest.raises(ValueError, match="percentiles"):
        accumulated_local_effects(model, np.ones((50, 2)), features=[0], percentiles=(0.9, 0.1))
