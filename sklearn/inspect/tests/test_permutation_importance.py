import pytest
import numpy as np

from sklearn.compose import ColumnTransformer
from sklearn.datasets import load_boston, load_iris
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspect import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


@pytest.mark.parametrize("to_pd", [True, False])
def test_permutation_importance_correlated_feature_regression(to_pd):
    rng = np.random.RandomState(42)

    bunch = load_boston()
    X, y = bunch.data, bunch.target
    y_with_little_noise = (
        y + rng.normal(scale=0.001, size=y.shape[0])).reshape(-1, 1)

    if to_pd:
        pd = pytest.importorskip("pandas")
        X = pd.DataFrame(X, columns=bunch.feature_names)
        X['correlated_feature'] = y_with_little_noise
    else:
        # Adds correlated feature to X
        X = np.hstack([X, y_with_little_noise])

    rf = RandomForestRegressor(n_estimators=10, random_state=rng)
    rf.fit(X, y)
    permute_scores = permutation_importance(rf, X, y, n_bootstrap=10,
                                            random_state=rng,
                                            scoring="neg_mean_absolute_error")

    assert permute_scores.shape == (X.shape[1], 10)

    permute_score_means = np.mean(permute_scores, axis=-1)

    # correlated feature has the highest importance
    assert np.all(permute_score_means[-1] > permute_score_means[:-1])


def test_permutation_importance_correlated_feature_column_transframer():
    pd = pytest.importorskip("pandas")
    rng = np.random.RandomState(42)

    bunch = load_iris()
    X, y = bunch.data, bunch.target
    y_with_little_noise = (
        y + rng.normal(scale=0.001, size=y.shape[0])).reshape(-1, 1)

    df = pd.DataFrame(X, columns=bunch.feature_names)
    df["correlated_feature"] = y_with_little_noise

    column_trans = ColumnTransformer(
        [("scale", StandardScaler(),
         ["sepal length (cm)", "sepal width (cm)",
          "petal length (cm)", "petal width (cm)"])],
        remainder='passthrough')
    model = Pipeline([("preprocessing", column_trans),
                      ("estimator",
                       LogisticRegression(multi_class='auto',
                                          solver='lbfgs'))])
    model.fit(df, y)
    permute_scores = permutation_importance(model, df, y, n_bootstrap=10,
                                            random_state=rng)
    assert permute_scores.shape == (df.shape[1], 10)

    permute_score_means = np.mean(permute_scores, axis=-1)

    # correlated feature has the highest importance
    assert np.all(permute_score_means[-1] > permute_score_means[:-1])
