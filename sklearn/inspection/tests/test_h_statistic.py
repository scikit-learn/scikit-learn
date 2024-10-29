# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
import pytest
from numpy.testing import assert_allclose

import sklearn
from sklearn.base import BaseEstimator, ClassifierMixin, RegressorMixin, clone
from sklearn.cluster import KMeans
from sklearn.compose import ColumnTransformer
from sklearn.datasets import make_classification, make_regression
from sklearn.ensemble import (
    HistGradientBoostingRegressor,
    RandomForestClassifier,
    RandomForestRegressor,
)
from sklearn.inspection import h_statistic
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OrdinalEncoder


class Deterministic(RegressorMixin):
    """Create estimator from deteministic prediction function"""

    def __init__(self, pred_fun):
        self._pred_fun = pred_fun
        self.is_fitted_ = True

    def fit(self, *args, **kwargs):
        # Simulate that something happens during fit()
        self.coef_ = 0
        return self

    def predict(self, X):
        return self._pred_fun(X)


@pytest.mark.parametrize(
    ["est", "is_additive"],
    (
        [LinearRegression(), True],
        [HistGradientBoostingRegressor(random_state=0, max_depth=1, max_iter=20), True],
        [RandomForestRegressor(random_state=0, max_depth=4, n_estimators=20), False],
        [HistGradientBoostingRegressor(random_state=0, max_iter=20), False],
    ),
)
@pytest.mark.parametrize("features", [None, [0, 1, 2]])
@pytest.mark.parametrize("subsample", [500, 50])
@pytest.mark.parametrize("sample_weight", [None, "ones"])
def test_h_statistic_regression(est, is_additive, features, subsample, sample_weight):
    """Checks on regression.

    Additive regression models should get statistics of 0, while non-additive
    models should get a value > 0.
    """
    N_FEAT = 4
    N = 200
    X, y = make_regression(n_samples=N, n_features=N_FEAT, random_state=0)
    w = np.ones(N) if sample_weight == "ones" else None

    model = est.fit(X, y, sample_weight=w)
    result = h_statistic(
        model,
        X,
        features=features,
        sample_weight=w,
        subsample=subsample,
        random_state=1,
    )

    m = N_FEAT if features is None else len(features)
    expected_length = m * (m - 1) / 2

    assert result.h_squared_pairwise.shape == (expected_length, 1)
    if is_additive:
        assert_allclose(result.h_squared_pairwise, 0)
    else:
        assert any(result.h_squared_pairwise > 1e-5)


def test_h_statistic_additive_classification():
    """Checks that additive classification gets statistics of 0.

    The presence of link functions (especially for GradientBoosting with
    more than 2 classes) would make further tests tricky.
    """
    N_FEAT = 4
    N_CLASSES = 4
    X, y = make_classification(
        n_samples=200,
        n_features=N_FEAT,
        n_informative=N_FEAT,
        n_redundant=0,
        n_classes=N_CLASSES,
        random_state=0,
    )
    model = RandomForestClassifier(n_estimators=20, random_state=0, max_depth=1)
    model.fit(X, y)
    result = h_statistic(model, X, random_state=1)

    expected_length = N_FEAT * (N_FEAT - 1) / 2

    assert result.h_squared_pairwise.shape == (expected_length, N_CLASSES)
    assert_allclose(result.h_squared_pairwise, 0)


def test_h_statistic_binary_classification_shape():
    """Checks that 0 class is dropped from binary classification."""
    N_FEAT = 4
    X, y = make_classification(
        n_samples=200,
        n_features=N_FEAT,
        n_informative=N_FEAT,
        n_redundant=0,
        n_classes=2,
        random_state=0,
    )
    model = RandomForestClassifier(n_estimators=20, random_state=0, max_depth=3)
    model.fit(X, y)
    result = h_statistic(model, X, random_state=1)

    expected_length = N_FEAT * (N_FEAT - 1) / 2

    assert result.h_squared_pairwise.shape == (expected_length, 1)


@pytest.mark.parametrize("sample_weight", [None, "ones"])
@pytest.mark.parametrize("subsample", [500, 50])
def test_h_statistic_equivalence_array_dataframe(subsample, sample_weight):
    """Checks that index operations give the same for numpy arrays or dataframes."""
    pd = pytest.importorskip("pandas")

    N = 200
    X, y = make_regression(n_samples=N, n_features=4, random_state=0)
    X_df = pd.DataFrame(X, columns=list("abcd"))
    w = np.ones(N) if sample_weight == "ones" else None

    # Numpy
    model_np = RandomForestRegressor(random_state=0, max_depth=4, n_estimators=20)
    model_np.fit(X, y, sample_weight=w)
    result_np = h_statistic(
        model_np, X, sample_weight=w, subsample=subsample, random_state=1
    )

    # Pandas
    model_pd = clone(model_np)
    model_pd.fit(X_df, y, sample_weight=w)
    result_pd = h_statistic(
        model_pd, X_df, sample_weight=w, subsample=subsample, random_state=1
    )

    assert_allclose(result_np.h_squared_pairwise, result_pd.h_squared_pairwise)


def test_h_statistic_reflect_interaction_constraints():
    """Checks that interaction constraints are reflected."""
    X, y = make_regression(n_samples=200, n_features=4, random_state=0)

    model = HistGradientBoostingRegressor(
        random_state=0, max_iter=20, interaction_cst=[[0, 1], [2], [3]]
    )
    model.fit(X, y)
    result = h_statistic(model, X, random_state=1)

    assert result.h_squared_pairwise[0] > 1e-5
    assert_allclose(result.h_squared_pairwise[1:], 0)


def test_h_statistic_matches_other_implementations():
    """Check against other implementations

    Check for one specific example that the result matches the one from two
    independent R packages ({iml} of Christoph Molnar, {hstats} of Michael Mayer).
    The latter also allows for sample weights.
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
    """

    X = np.array([[1, 1, 1, 2], [0, 0, 2, 1], [2, 2, 1, 4]]).T

    # Model with interaction in last feature pair
    model = Deterministic(lambda X: X[:, 0] + X[:, 1] * np.sqrt(X[:, 2]))

    assert_allclose(h_statistic(model, X=X).h_squared_pairwise[2], 0.078125)
    assert_allclose(
        h_statistic(model, X=X, sample_weight=range(4)).h_squared_pairwise[2],
        0.1010093,
    )


def test_h_statistic_matches_analytic_result():
    """Check against analytic result

    By definition of the H-squared, a model consisting of an interaction
    and without main effects would result in a value of 1.
    """

    X = np.array([[1, 1, -1, -1], [1, -1, 1, -1]]).T

    # Model consistent of pure interaction (and zero main effects)
    model = Deterministic(lambda X: X[:, 0] * X[:, 1])

    assert_allclose(h_statistic(model, X=X).h_squared_pairwise[0], 1)


@pytest.mark.parametrize("sample_weight", [None, "ones"])
@pytest.mark.parametrize("subsample", [500, 50])
def test_h_statistic_does_not_change_pandas_input(subsample, sample_weight):
    """Checks that pandas data is unchanged by the function call."""
    pd = pytest.importorskip("pandas")

    N = 200
    X, y = make_regression(n_samples=N, n_features=4, random_state=0)
    X_df = pd.DataFrame(X)
    X_df_orig = X_df.copy()
    w = np.ones(N) if sample_weight == "ones" else None

    # Pandas
    model = RandomForestRegressor(random_state=0, max_depth=4, n_estimators=20)
    model.fit(X_df, y, sample_weight=w)
    _ = h_statistic(model, X_df, sample_weight=w, subsample=subsample, random_state=1)

    pd.testing.assert_frame_equal(X_df, X_df_orig)


def test_h_statistic_works_with_mixed_type_pandas_data():
    """Mixed-type input should work even when np.unique() fails."""
    pd = pytest.importorskip("pandas")

    y = np.array([0, 1, 2, 3])
    X = np.array([y, y]).T
    X_mixed = pd.DataFrame({"x1": y, "x2": list("abcd")})

    # Model on homogeneous numpy array
    model = RandomForestRegressor(random_state=0).fit(X, y)

    # Identical model on mixed-type pandas df
    preprocessor = ColumnTransformer(
        transformers=[
            ("keep", "passthrough", ["x1"]),
            ("ordinal", OrdinalEncoder(), ["x2"]),
        ]
    )
    model_mixed = Pipeline(
        steps=[
            ("prep", preprocessor),
            ("rf", RandomForestRegressor(random_state=0)),
        ]
    )
    model_mixed.fit(X_mixed, y)

    assert_allclose(
        h_statistic(model, X).h_squared_pairwise,
        h_statistic(model_mixed, X_mixed).h_squared_pairwise,
    )


@pytest.mark.parametrize(
    "Estimator",
    (
        sklearn.tree.DecisionTreeClassifier,
        sklearn.ensemble.RandomForestClassifier,
    ),
)
def test_h_statistic_raises_for_multioutput(Estimator):
    """Make sure error is raised for multiclass-multioutput classifiers

    Slightly modified from test_partial_dependence.py
    """
    # make multiclass-multioutput dataset
    X, y = make_classification(n_classes=3, n_clusters_per_class=1, random_state=0)
    y = np.array([y, y]).T
    est = Estimator().fit(X, y)

    with pytest.raises(
        ValueError, match="Multiclass-multioutput estimators are not supported"
    ):
        h_statistic(est, X, [0, 1])


class RegressorWithoutPredict(RegressorMixin, BaseEstimator):
    def fit(self, X, y):
        # Simulate that something happens during fit()
        self.coef_ = 0
        return self

    def predict_proba(self, X):
        return np.zeros(X.shape[0])  # pragma: no cover


class ClassifierWithoutPredictProba(ClassifierMixin, BaseEstimator):
    def fit(self, X, y):
        # Simulate that we have some classes
        self.classes_ = [0, 1]
        return self

    def predict(self, X):
        return np.zeros(X.shape[0])  # pragma: no cover


@pytest.mark.parametrize(
    "estimator, err_msg",
    [
        (
            RegressorWithoutPredict(),
            "The regressor has no predict method",
        ),
        (
            ClassifierWithoutPredictProba(),
            "The classifier has no predict_proba method",
        ),
        (
            KMeans(n_clusters=2, random_state=0, n_init="auto"),
            "'estimator' must be a regressor or classifier",
        ),
    ],
)
def test_h_statistic_error(estimator, err_msg):
    y = np.array([0, 0, 1, 1])
    X = np.zeros((len(y), 2))
    estimator.fit(X, y)

    with pytest.raises(ValueError, match=err_msg):
        h_statistic(estimator, X)
