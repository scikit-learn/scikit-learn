"""Tests for the accumulated local effects."""

from itertools import pairwise

import numpy as np
import pytest

from sklearn.base import BaseEstimator, RegressorMixin, clone, is_classifier
from sklearn.cluster import KMeans
from sklearn.compose import make_column_transformer
from sklearn.datasets import make_classification, make_regression
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.exceptions import NotFittedError
from sklearn.inspection import accumulated_local_effects
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor
from sklearn.utils._testing import assert_allclose

# Reference values from ALEPlot, the R package by D. W. Apley, author of ALE.
# The data come from `_aleplot_data()` and are passed to R with
# `to_csv("data.csv", index=False, float_format="%.17g")`, so that both languages
# use identical inputs. The model is a fixed formula, `_ALEPlotFormula`. In R:
#
#   library(ALEPlot)  # version 1.1, archived on CRAN; R 4.6.1, yaImpute 1.0.36
#   X <- read.csv("data.csv",
#                 colClasses = c("numeric", "numeric", "integer", "factor", "numeric"))
#   pred <- function(X.model, newdata) {
#     effect <- ifelse(newdata$c == "a", 0, ifelse(newdata$c == "b", 1, 3))
#     newdata$x1 + 2 * newdata$x1 * newdata$x2 + newdata$x2^2 + 0.5 * newdata$x3 +
#       newdata$x1 * newdata$x3 + effect + newdata$x2 * (newdata$c == "c") +
#       1.5 * newdata$x1 * newdata$x4
#   }
#   ALEPlot(X, NULL, pred, J = 1, K = 5)        # continuous x1
#   ALEPlot(X, NULL, pred, J = 3, K = 10)       # integer x3
#   ALEPlot(X, NULL, pred, J = 4)               # factor c
#   ALEPlot(X, NULL, pred, J = c(1, 2), K = 4)  # pair x1, x2
#   ALEPlot(X, NULL, pred, J = c(1, 5), K = 4)  # pair x1, x4, with empty cells
#
# ALEPlot returns `x.values` (the grid: interval bounds, or levels for a factor)
# and `f.values` (the ALE at each grid value, a matrix for a pair). They are
# hard-coded below as `expected_grid` and `expected_effects`. ALEPlot's `K` is the
# number of intervals, so it corresponds to `grid_resolution = K + 1`.


def _aleplot_data():
    """Deterministic data with 61 samples, the input to ALEPlot in the R code above."""
    pd = pytest.importorskip("pandas")
    i = np.arange(1, 62)
    x1 = np.mod(i * 0.6180339887498949, 1.0)
    x2 = np.mod(i * 0.7548776662466927, 1.0)
    x3 = np.mod(i * 7, 5) + 1
    score = x1 + 0.6 * x2
    c = np.where(score < 0.55, "a", np.where(score < 1.05, "b", "c"))
    x4 = x1 + 0.1 * (x2 - 0.5)
    return pd.DataFrame({"x1": x1, "x2": x2, "x3": x3, "c": c, "x4": x4})


class _ALEPlotFormula(RegressorMixin, BaseEstimator):
    """Python version of `pred`, the formula used as the model in the R code above."""

    def fit(self, X, y=None):
        self.is_fitted_ = True
        return self

    def predict(self, X):
        effect = np.select([X["c"] == "a", X["c"] == "b"], [0.0, 1.0], 3.0)
        return (
            X["x1"]
            + 2 * X["x1"] * X["x2"]
            + X["x2"] ** 2
            + 0.5 * X["x3"]
            + X["x1"] * X["x3"]
            + effect
            + X["x2"] * (X["c"] == "c")
            + 1.5 * X["x1"] * X["x4"]
        ).to_numpy()


@pytest.mark.parametrize(
    "feature, grid_resolution, expected_grid, expected_effects",
    [
        (
            "x1",
            6,
            [
                0.013155617496426686,
                0.214781741247582758,
                0.416407864998738775,
                0.609903369994114541,
                0.798373876248844039,
                0.991869381244221415,
            ],
            [
                -2.69815306011557565,
                -1.68862250227735022,
                -0.54584864967921209,
                0.62101287133147443,
                1.69314578560799189,
                2.90434268034917986,
            ],
        ),
        (
            "x3",
            11,
            [1, 2, 3, 4, 5],
            [
                -1.71518559878436272,
                -0.70544828649600677,
                0.33822074936611979,
                1.30629139498780966,
                2.28829406310970995,
            ],
        ),
    ],
    ids=["continuous", "integer"],
)
def test_matches_r_aleplot(feature, grid_resolution, expected_grid, expected_effects):
    X = _aleplot_data()
    model = _ALEPlotFormula().fit(X)
    result = accumulated_local_effects(
        model, X, feature, grid_resolution=grid_resolution, categorical_features=None
    )
    assert_allclose(result.grid_values[0], expected_grid, rtol=0, atol=1e-12)
    assert_allclose(result.effects[0], expected_effects, rtol=0, atol=1e-12)


def test_matches_r_aleplot_categorical():
    # ALEPlot orders the levels by multidimensional scaling, here c, b, a. That
    # is the reverse of our sorted order, which gives the same centered effects.
    X = _aleplot_data()
    model = _ALEPlotFormula().fit(X)
    result = accumulated_local_effects(model, X, "c", categorical_features=["c"])
    expected = {
        "c": 2.521322762565242570,
        "b": -0.032536631982445829,
        "a": -1.032536631982445829,
    }
    assert list(result.grid_values[0]) == ["a", "b", "c"]
    assert_allclose(
        result.effects[0], [expected[level] for level in "abc"], rtol=0, atol=1e-12
    )


@pytest.mark.parametrize(
    "features, expected_grid, expected_effects",
    [
        (
            ["x1", "x2"],
            [
                [
                    0.013155617496426686,
                    0.249223594996216491,
                    0.506577808748213343,
                    0.742645786248003148,
                    0.991869381244221415,
                ],
                [
                    0.008516311074714622,
                    0.253638644828022564,
                    0.498760978581330505,
                    0.743883312334638447,
                    0.989005646087946610,
                ],
            ],
            [
                [
                    0.455048135573065782,
                    0.208822914821934102,
                    -0.010041951308413122,
                    -0.228773951022475552,
                    -0.487556484252403610,
                ],
                [
                    0.245285576384817156,
                    0.114791422772029916,
                    0.011657623780026327,
                    -0.091343308795692080,
                    -0.234394774887276558,
                ],
                [
                    -0.0364863427419286346,
                    -0.0408139654024416987,
                    -0.0177812334421711404,
                    0.0053843649343849886,
                    -0.0115005702049256753,
                ],
                [
                    -0.2486599658288918269,
                    -0.1372565213510614779,
                    0.0015072777475531041,
                    0.1404039432624520911,
                    0.2392500752614853399,
                ],
                [
                    -0.4970937273717983551,
                    -0.2635097444302436642,
                    -0.0025654068679046294,
                    0.2585117971107195878,
                    0.4795384675734768454,
                ],
            ],
        ),
        (
            ["x1", "x4"],
            [
                [
                    0.013155617496426686,
                    0.249223594996216491,
                    0.506577808748213343,
                    0.742645786248003148,
                    0.991869381244221415,
                ],
                [
                    0.029739682735181815,
                    0.237043510749128356,
                    0.519609524102774833,
                    0.711130532247384983,
                    1.013956862865847519,
                ],
            ],
            [
                [
                    0.0319952171203266439,
                    -0.0047081294350857991,
                    -0.1396257052697704726,
                    -0.3213230442399233611,
                    -0.6256763055214680236,
                ],
                [
                    -0.0047081294350857991,
                    0.0319952171203266439,
                    -0.0295156656035331227,
                    -0.1372798968716180679,
                    -0.3677000504510947732,
                ],
                [
                    -0.174427541017715015,
                    -0.028644863162806698,
                    0.018923585412829388,
                    -0.014907538153187641,
                    -0.171394584030596320,
                ],
                [
                    -0.408690026252510552,
                    -0.162850167473741869,
                    -0.015224537974245375,
                    0.018762304021866225,
                    -0.030493342751367350,
                ],
                [
                    -0.702385318712688056,
                    -0.356488279010058917,
                    -0.137265318259912728,
                    -0.031681145013151404,
                    0.032270408377985046,
                ],
            ],
        ),
    ],
    ids=["independent", "empty-cells"],
)
def test_matches_r_aleplot_two_way(features, expected_grid, expected_effects):
    X = _aleplot_data()
    model = _ALEPlotFormula().fit(X)
    result = accumulated_local_effects(
        model, X, features, grid_resolution=5, categorical_features=None
    )
    for grid, expected in zip(result.grid_values, expected_grid):
        assert_allclose(grid, expected, rtol=0, atol=1e-12)

    # ALEPlot fills empty cells with the nearest non-empty cell. That only changes
    # the effects at grid points that are not a corner of any non-empty cell, so
    # only the grid points that touch data are compared.
    occupied = result.bin_sizes > 0
    touches_data = np.zeros(result.effects.shape[1:], dtype=bool)
    for i in (0, 1):
        for j in (0, 1):
            touches_data[i : i + occupied.shape[0], j : j + occupied.shape[1]] |= (
                occupied
            )
    assert (~occupied).any() == (features == ["x1", "x4"])
    assert_allclose(
        result.effects[0][touches_data],
        np.asarray(expected_effects)[touches_data],
        rtol=0,
        atol=1e-12,
    )


class _AdditiveRegressor(RegressorMixin, BaseEstimator):
    """Regressor predicting x0 ** 2 + x0 * x1 + x2 without fitting."""

    def fit(self, X, y=None):
        self.is_fitted_ = True
        return self

    def predict(self, X):
        X = np.asarray(X, dtype=np.float64)
        return X[:, 0] ** 2 + X[:, 0] * X[:, 1] + X[:, 2]


@pytest.mark.parametrize(
    "estimator, X, y, n_outputs",
    [
        (LogisticRegression(), *make_classification(n_samples=50, random_state=0), 1),
        (
            LogisticRegression(),
            *make_classification(
                n_samples=50, n_classes=3, n_informative=3, random_state=0
            ),
            3,
        ),
        (LinearRegression(), *make_regression(n_samples=50, random_state=0), 1),
        (
            LinearRegression(),
            *make_regression(n_samples=50, n_targets=2, random_state=0),
            2,
        ),
    ],
    ids=["binary", "multiclass", "regression", "multioutput-regression"],
)
@pytest.mark.parametrize("grid_resolution", [2, 10, 100])
def test_output_shape(estimator, X, y, n_outputs, grid_resolution):
    estimator.fit(X, y)
    result = accumulated_local_effects(
        estimator, X, features=0, grid_resolution=grid_resolution
    )
    n_grid_values = min(grid_resolution, X.shape[0])
    assert result.effects.shape == (n_outputs, n_grid_values)
    assert len(result.grid_values) == 1
    assert result.grid_values[0].shape == (n_grid_values,)
    assert result.bin_sizes.shape == (n_grid_values - 1,)


def test_linear_model_effects_have_the_coefficient_as_slope():
    # from ALE’s definition: for a linear model, the ALE of each feature
    # is a straight line whose slope is that feature’s coefficient.
    X, y = make_regression(n_samples=100, n_features=3, random_state=0)
    reg = LinearRegression().fit(X, y)
    for feature in range(X.shape[1]):
        result = accumulated_local_effects(reg, X, feature, grid_resolution=10)
        slopes = np.diff(result.effects[0]) / np.diff(result.grid_values[0])
        assert_allclose(slopes, reg.coef_[feature])


def test_additive_effect_is_exact_at_grid_values():
    # x2 is added to the model (f = x0**2 + x0 * x1 + x2), so its ALE must equal x2
    # plus a constant, whatever the non-linear terms and interactions of x0 and x1.
    X = np.random.RandomState(0).uniform(-1, 1, size=(200, 3))
    reg = _AdditiveRegressor().fit(X)
    result = accumulated_local_effects(reg, X, features=2, grid_resolution=20)
    assert_allclose(np.diff(result.effects[0]), np.diff(result.grid_values[0]))


# Reference implementation: a slow, loop-based reproduction of the continuous ALE,
# used to check the vectorised code in `accumulated_local_effects`.
def _ale_reference_continuous(predict, X, feature, edges):
    """Loop implementation of Apley and Zhu (2020), Eq. 15 and ALEPlot."""
    x = X[:, feature]
    n_bins = len(edges) - 1
    local_effects, counts = np.zeros(n_bins), np.zeros(n_bins)
    for k in range(n_bins):
        in_bin = (x > edges[k]) & (x <= edges[k + 1])
        if k == 0:
            in_bin |= x == edges[0]
        X_upper, X_lower = X[in_bin].copy(), X[in_bin].copy()
        X_upper[:, feature] = edges[k + 1]
        X_lower[:, feature] = edges[k]
        local_effects[k] = np.mean(predict(X_upper) - predict(X_lower))
        counts[k] = in_bin.sum()
    ale = np.concatenate([[0], np.cumsum(local_effects)])
    ale -= np.sum((ale[:-1] + ale[1:]) / 2 * counts) / counts.sum()
    return ale


@pytest.mark.parametrize("grid_resolution", [3, 20, 100])
@pytest.mark.parametrize("with_ties", [False, True], ids=["continuous", "ties"])
def test_continuous_matches_reference(grid_resolution, with_ties):
    rng = np.random.RandomState(0)
    X = rng.normal(size=(300, 3))
    if with_ties:
        X[:, 0] = rng.randint(0, 7, size=300)
    y = X[:, 0] ** 2 + X[:, 0] * X[:, 1] + rng.normal(scale=0.1, size=300)
    reg = DecisionTreeRegressor(max_depth=5, random_state=0).fit(X, y)

    result = accumulated_local_effects(reg, X, 0, grid_resolution=grid_resolution)
    # Recompute the ALE with the loop-based reference, on the same edges.
    edges = result.grid_values[0]
    expected = _ale_reference_continuous(reg.predict, X, 0, edges)
    assert_allclose(result.effects[0], expected)


# Reference implementation: a slow, loop-based reproduction of the categorical ALE,
# used to check the vectorised code in `accumulated_local_effects`.
def _ale_reference_categorical(predict, X, feature, levels):
    """Loop implementation of the ALEPlot algorithm for ordered levels."""
    x = X[:, feature]
    local_effects = []
    for lower, upper in pairwise(levels):
        X_up, X_down = X[x == lower].copy(), X[x == upper].copy()
        X_up_moved, X_down_moved = X_up.copy(), X_down.copy()
        X_up_moved[:, feature] = upper
        X_down_moved[:, feature] = lower
        differences = np.concatenate(
            [
                predict(X_up_moved) - predict(X_up),
                predict(X_down) - predict(X_down_moved),
            ]
        )
        local_effects.append(differences.mean())
    ale = np.concatenate([[0], np.cumsum(local_effects)])
    counts = np.array([np.sum(x == level) for level in levels])
    return ale - np.sum(ale * counts) / counts.sum()


def test_categorical_matches_reference():
    rng = np.random.RandomState(0)
    X = rng.normal(size=(300, 3))
    X[:, 1] = rng.choice([0, 2, 3, 7], size=300)
    y = np.sin(X[:, 1]) + X[:, 0] * X[:, 1] + rng.normal(scale=0.1, size=300)
    reg = DecisionTreeRegressor(max_depth=5, random_state=0).fit(X, y)

    result = accumulated_local_effects(reg, X, 1, categorical_features=[1])
    assert_allclose(result.grid_values[0], [0, 2, 3, 7])
    # Recompute the ALE with the loop-based reference, on the same levels.
    expected = _ale_reference_categorical(reg.predict, X, 1, result.grid_values[0])
    assert_allclose(result.effects[0], expected)
    # The number of samples at each level.
    assert_allclose(result.bin_sizes, [np.sum(X[:, 1] == v) for v in [0, 2, 3, 7]])


@pytest.mark.parametrize("grid_resolution", [5, 50])
def test_grid_values_are_observed_inverted_cdf_quantiles(grid_resolution):
    rng = np.random.RandomState(0)
    X = np.c_[rng.exponential(size=200).round(1), rng.normal(size=200)]
    reg = LinearRegression().fit(X, X[:, 0])

    result = accumulated_local_effects(reg, X, 0, grid_resolution=grid_resolution)
    expected = np.unique(
        np.percentile(
            X[:, 0], np.linspace(0, 100, grid_resolution), method="inverted_cdf"
        )
    )
    assert_allclose(result.grid_values[0], expected)
    assert np.isin(result.grid_values[0], X[:, 0]).all()
    assert (result.bin_sizes > 0).all()
    assert result.bin_sizes.sum() == X.shape[0]


def test_missing_values_are_ignored():
    rng = np.random.RandomState(0)
    X = rng.normal(size=(200, 2))
    reg = HistGradientBoostingRegressor(max_iter=10).fit(X, X[:, 0] + X[:, 1])
    missing = rng.rand(200) < 0.2
    X_missing = X.copy()
    X_missing[missing, 0] = np.nan

    result = accumulated_local_effects(reg, X_missing, 0, grid_resolution=10)
    expected = accumulated_local_effects(reg, X[~missing], 0, grid_resolution=10)
    assert_allclose(result.grid_values[0], expected.grid_values[0])
    assert_allclose(result.effects, expected.effects)


def test_two_way_missing_values_are_ignored():
    rng = np.random.RandomState(0)
    X = rng.normal(size=(300, 3))
    reg = HistGradientBoostingRegressor(max_iter=10).fit(X, X[:, 0] * X[:, 1])
    missing_0, missing_1 = rng.rand(300) < 0.1, rng.rand(300) < 0.1
    X_missing = X.copy()
    X_missing[missing_0, 0] = np.nan
    X_missing[missing_1, 1] = np.nan
    kept = ~(missing_0 | missing_1)

    result = accumulated_local_effects(reg, X_missing, [0, 1], grid_resolution=6)
    expected = accumulated_local_effects(reg, X[kept], [0, 1], grid_resolution=6)
    assert_allclose(result.effects, expected.effects)
    assert_allclose(result.bin_sizes, expected.bin_sizes)


@pytest.mark.parametrize(
    "estimator, params, err_msg",
    [
        (
            KMeans(n_clusters=2, n_init=1),
            {"features": 0},
            "'estimator' must be a fitted regressor or classifier",
        ),
        (DecisionTreeClassifier(), {"features": 0}, "Multiclass-multioutput"),
        (LinearRegression(), {"features": [0, 1, 2]}, "single feature or a pair"),
        (LinearRegression(), {"features": [-1]}, "all features must be in"),
        (
            LinearRegression(),
            {"features": 0, "categorical_features": [True]},
            "must have shape",
        ),
        (
            LinearRegression(),
            {"features": 0, "response_method": "predict_proba"},
            "The response_method parameter is ignored for regressors",
        ),
        (LinearRegression(), {"features": 3}, "at least two distinct values"),
    ],
)
def test_errors(estimator, params, err_msg):
    rng = np.random.RandomState(0)
    X = np.c_[rng.normal(size=(20, 3)), np.ones(20)]
    # The classifier is fitted on two targets with three classes each, to make
    # it multiclass-multioutput.
    y = np.c_[rng.randint(0, 3, 20), rng.randint(0, 3, 20)]
    estimator = clone(estimator).fit(X, y if is_classifier(estimator) else X[:, 0])
    with pytest.raises(ValueError, match=err_msg):
        accumulated_local_effects(estimator, X, **params)


def test_raises_non_fitted_error():
    with pytest.raises(NotFittedError):
        accumulated_local_effects(LinearRegression(), np.ones((5, 2)), 0)


def test_all_missing_values_error():
    X = np.c_[np.arange(20.0), np.arange(20.0)]
    reg = HistGradientBoostingRegressor(max_iter=5).fit(X, X[:, 1])
    X[:, 0] = np.nan
    with pytest.raises(ValueError, match="No sample is left"):
        accumulated_local_effects(reg, X, 0)


def test_categorical_features_from_dtype():
    pd = pytest.importorskip("pandas")
    rng = np.random.RandomState(0)
    colors = rng.choice(["blue", "green", "red"], size=100)
    X = pd.DataFrame({"a": rng.normal(size=100), "color": pd.Categorical(colors)})
    y = X["a"] + (colors == "red")
    model = make_pipeline(
        make_column_transformer((OneHotEncoder(), ["color"]), remainder="passthrough"),
        LinearRegression(),
    ).fit(X, y)

    from_dtype = accumulated_local_effects(model, X, "color")
    explicit = accumulated_local_effects(
        model, X, "color", categorical_features=["color"]
    )
    assert list(from_dtype.grid_values[0]) == ["blue", "green", "red"]
    assert_allclose(from_dtype.effects, explicit.effects)


def test_categorical_features_from_dtype_numpy():
    X = np.random.RandomState(0).randint(0, 4, size=(50, 2)).astype(float)
    reg = LinearRegression().fit(X, X[:, 0])
    from_dtype = accumulated_local_effects(reg, X, 0, grid_resolution=3)
    continuous = accumulated_local_effects(
        reg, X, 0, categorical_features=None, grid_resolution=3
    )
    assert_allclose(from_dtype.effects, continuous.effects)


@pytest.mark.parametrize(
    "estimator, n_classes", [(LinearRegression(), 1), (LogisticRegression(), 3)]
)
def test_two_way_output_shape(estimator, n_classes):
    rng = np.random.RandomState(0)
    X = rng.normal(size=(200, 3))
    y = np.digitize(X[:, 0], [-0.5, 0.5]) if n_classes == 3 else X[:, 0]
    estimator.fit(X, y)
    result = accumulated_local_effects(estimator, X, [0, 2], grid_resolution=6)
    assert [g.shape for g in result.grid_values] == [(6,), (6,)]
    assert result.effects.shape == (n_classes, 6, 6)
    assert result.bin_sizes.shape == (5, 5)
    assert result.bin_sizes.sum() == X.shape[0]


def test_two_way_additive_model_has_no_interaction():
    X = np.random.RandomState(0).uniform(-1, 1, size=(300, 3))
    reg = _AdditiveRegressor().fit(X)
    result = accumulated_local_effects(reg, X, [0, 2], grid_resolution=8)
    assert_allclose(result.effects, 0, atol=1e-12)


# Reference implementation: a slow, loop-based reproduction of the second-order ALE,
# used to check the vectorised code in `accumulated_local_effects`.
def _interval_index(x, edges):
    """Index k of the interval (edges[k], edges[k + 1]], with the minimum in 0."""
    index = np.zeros(x.shape[0], dtype=int)
    for k in range(len(edges) - 1):
        in_interval = (x > edges[k]) & (x <= edges[k + 1])
        if k == 0:
            in_interval |= x == edges[0]
        index[in_interval] = k
    return index


def _ale_reference_2d(predict, X, features, edges):
    """Loop implementation of ALEPlot's second-order ALE of two numeric features.

    Empty cells get a mean local effect of 0 instead of the nearest-neighbour value
    used by ALEPlot.
    """
    j, k = features
    z1, z2 = edges
    a1, a2 = _interval_index(X[:, j], z1), _interval_index(X[:, k], z2)
    delta = np.zeros((len(z1) - 1, len(z2) - 1))
    b = np.zeros_like(delta)
    for r in range(delta.shape[0]):
        for c in range(delta.shape[1]):
            rows = (a1 == r) & (a2 == c)
            b[r, c] = rows.sum()
            if not rows.any():
                continue

            def at(u, v):
                X_moved = X[rows].copy()
                X_moved[:, j], X_moved[:, k] = u, v
                return predict(X_moved)

            delta[r, c] = np.mean(
                at(z1[r + 1], z2[c + 1])
                - at(z1[r + 1], z2[c])
                - at(z1[r], z2[c + 1])
                + at(z1[r], z2[c])
            )
    f = np.zeros((len(z1), len(z2)))
    f[1:, 1:] = delta.cumsum(axis=0).cumsum(axis=1)
    d = f[1:, :] - f[:-1, :]
    f1 = np.r_[0, np.cumsum((b * (d[:, :-1] + d[:, 1:]) / 2).sum(1) / b.sum(1))]
    d = f[:, 1:] - f[:, :-1]
    f2 = np.r_[0, np.cumsum((b * (d[:-1, :] + d[1:, :]) / 2).sum(0) / b.sum(0))]
    f = f - f1[:, None] - f2[None, :]
    f0 = np.sum(b * (f[:-1, :-1] + f[:-1, 1:] + f[1:, :-1] + f[1:, 1:]) / 4) / b.sum()
    return f - f0, b


@pytest.mark.parametrize(
    "correlated", [False, True], ids=["independent", "empty-cells"]
)
def test_two_way_matches_reference(correlated):
    rng = np.random.RandomState(0)
    X = rng.uniform(0, 1, size=(400, 3))
    if correlated:
        X[:, 1] = X[:, 0] + rng.normal(scale=0.05, size=400)
    y = X[:, 0] * X[:, 1] + np.sin(3 * X[:, 0]) + X[:, 2]
    reg = DecisionTreeRegressor(max_depth=6, random_state=0).fit(X, y)

    result = accumulated_local_effects(reg, X, [0, 1], grid_resolution=8)
    # Recompute the ALE with the loop-based reference, on the same grid.
    expected, cell_sizes = _ale_reference_2d(reg.predict, X, (0, 1), result.grid_values)
    assert_allclose(result.bin_sizes, cell_sizes)
    assert (result.bin_sizes == 0).any() == correlated
    assert_allclose(result.effects[0], expected, atol=1e-12)


def test_two_way_categorical_not_supported():
    X = np.c_[np.random.RandomState(0).randint(0, 3, size=(50, 2)), np.ones(50)]
    reg = LinearRegression().fit(X, X[:, 0])
    with pytest.raises(NotImplementedError, match="two continuous features"):
        accumulated_local_effects(reg, X, [0, 1], categorical_features=[0])
