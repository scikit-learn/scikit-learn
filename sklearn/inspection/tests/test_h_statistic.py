import numpy as np
import pytest
from numpy.testing import assert_allclose

from sklearn.base import RegressorMixin, clone
from sklearn.datasets import make_classification, make_regression
from sklearn.ensemble import (
    HistGradientBoostingRegressor,
    RandomForestClassifier,
    RandomForestRegressor,
)
from sklearn.inspection import h_statistic
from sklearn.linear_model import LinearRegression


@pytest.mark.parametrize(
    "est",
    (
        LinearRegression(),
        HistGradientBoostingRegressor(random_state=0, max_depth=1, max_iter=20),
    ),
)
@pytest.mark.parametrize("features", [None, [0, 1, 2]])
@pytest.mark.parametrize("n_max", [500, 50])
@pytest.mark.parametrize("sample_weight", [None, "ones"])
def test_h_statistic_additive_regression(est, features, n_max, sample_weight):
    # Checks that additive regressions get statistics of 0.
    M_FEAT = 4
    N = 200
    X, y = make_regression(n_samples=N, n_features=M_FEAT, random_state=0)
    w = np.ones(N) if sample_weight == "ones" else None

    model = est.fit(X, y, sample_weight=w)
    result = h_statistic(
        model, X, features=features, random_state=1, n_max=n_max, sample_weight=w
    )

    m = M_FEAT if features is None else len(features)
    expected_length = m * (m - 1) / 2

    assert result.h_squared_pairwise.shape == (expected_length,)
    assert_allclose(result.h_squared_pairwise, 0)


def test_h_statistic_additive_classification():
    # Checks that additive classification gets statistics of 0. The presence
    # of link functions (especially for GradientBoosting with more than 2 classes)
    # would make further tests tricky.
    M_FEAT = 4
    N_CLASSES = 4
    X, y = make_classification(
        n_samples=200,
        n_features=M_FEAT,
        n_informative=M_FEAT,
        n_redundant=0,
        n_classes=N_CLASSES,
        random_state=0,
    )
    model = RandomForestClassifier(n_estimators=20, random_state=0, max_depth=1)
    model.fit(X, y)
    result = h_statistic(model, X, random_state=1)

    expected_length = M_FEAT * (M_FEAT - 1) / 2

    assert result.h_squared_pairwise.shape == (expected_length, N_CLASSES)
    assert_allclose(result.h_squared_pairwise, 0)


@pytest.mark.parametrize(
    "est",
    (
        RandomForestRegressor(random_state=0, max_depth=4, n_estimators=20),
        HistGradientBoostingRegressor(random_state=0, max_iter=20),
    ),
)
@pytest.mark.parametrize("features", [None, [0, 1, 2]])
@pytest.mark.parametrize("n_max", [500, 50])
@pytest.mark.parametrize("sample_weight", [None, "ones"])
def test_h_statistic_regression(est, features, n_max, sample_weight):
    # Tests that models with interactions will produce (some) positive values.
    M_FEAT = 4
    N = 200
    X, y = make_regression(n_samples=N, n_features=M_FEAT, random_state=0)
    w = np.ones(N) if sample_weight == "ones" else None

    model = est.fit(X, y, sample_weight=w)
    result = h_statistic(
        model, X, features=features, random_state=1, n_max=n_max, sample_weight=w
    )

    assert any(result.h_squared_pairwise > 1e-5)


@pytest.mark.parametrize("sample_weight", [None, "ones"])
@pytest.mark.parametrize("n_max", [500, 50])
def test_h_statistic_equivalence_array_dataframe(n_max, sample_weight):
    # Checks that index operations give the same for numpy arrays or dataframes.
    pd = pytest.importorskip("pandas")

    N = 200
    X, y = make_regression(n_samples=N, n_features=4, random_state=0)
    X_df = pd.DataFrame(X)
    w = np.ones(N) if sample_weight == "ones" else None

    # Numpy
    model_np = RandomForestRegressor(random_state=0, max_depth=4, n_estimators=20)
    model_np.fit(X, y, sample_weight=w)
    result_np = h_statistic(model_np, X, random_state=1, n_max=n_max, sample_weight=w)

    # Pandas
    model_pd = clone(model_np)
    model_pd.fit(X_df, y, sample_weight=w)
    result_pd = h_statistic(
        model_pd, X_df, random_state=1, n_max=n_max, sample_weight=w
    )

    assert_allclose(result_np.h_squared_pairwise, result_pd.h_squared_pairwise)


def test_h_statistic_reflect_interaction_constraints():
    # Checks that interaction constraints are reflected.
    X, y = make_regression(n_samples=200, n_features=4, random_state=0)

    model = HistGradientBoostingRegressor(
        random_state=0, max_iter=20, interaction_cst=[[0, 1], [2], [3]]
    )
    model.fit(X, y)
    result = h_statistic(model, X, random_state=1)

    assert result.h_squared_pairwise[0] > 1e-5
    assert_allclose(result.h_squared_pairwise[1:], 0)


def test_h_statistic_matches_other_implementations():
    # Checks for one specific example that the result matches the one from two
    # independent R packages ({iml} of Christoph Molnar, {hstats} of Michael Mayer).
    # The latter also allows for sample weights.
    # X <- data.frame(x1 = c(1, 1, 1, 2), x2 = c(0, 0, 2, 1), x3 = c(2, 2, 1, 4))
    # pred_fun <- function(model, newdata)
    #   newdata[, 1] + newdata[, 2] * sqrt(newdata[, 3])
    #
    # library(iml)     # 0.11.1
    # mod <- Predictor$new(data = X, predict.function = pred_fun)
    # Interaction$new(mod, grid.size = 5, feature = "x3")$results$.interaction^2
    # -> x2:x3 -> 0.078125
    #
    # library(hstats)  # 1.1.2
    # h2_pairwise(hstats(X = X, pred_fun = pred_fun))           # x2:x3 -> 0.078125
    # h2_pairwise(hstats(X = X, pred_fun = pred_fun, w = 0:3))  # weights -> 0.1010093
    class Deterministic(RegressorMixin):
        def __init__(self, pred_fun):
            self._pred_fun = pred_fun
            self.is_fitted_ = True

        def fit(self, *args, **kwargs):
            return self

        def predict(self, X):
            return self._pred_fun(X)

    X = np.array([[1, 1, 1, 2], [0, 0, 2, 1], [2, 2, 1, 4]]).T

    # Model with interaction in last feature pair
    model = Deterministic(lambda X: X[:, 0] + X[:, 1] * np.sqrt(X[:, 2]))

    assert_allclose(h_statistic(model, X=X).h_squared_pairwise[2], 0.078125)
    assert_allclose(
        h_statistic(model, X=X, sample_weight=range(4)).h_squared_pairwise[2],
        0.1010093,
    )
