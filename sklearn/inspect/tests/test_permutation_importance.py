import numpy as np

from sklearn.datasets import load_boston
from sklearn.inspect import permutation_importance
from sklearn.ensemble import RandomForestRegressor


def test_permutation_importance_correlated_feature_is_important():
    scoring = "neg_mean_absolute_error"
    rng = np.random.RandomState(42)
    X, y = load_boston(return_X_y=True)

    # Adds correlated feature to X
    y_with_little_noise = y + rng.normal(scale=0.001, size=y.shape[0])
    X = np.hstack([X, y_with_little_noise.reshape(-1, 1)])

    rf = RandomForestRegressor(n_estimators=50, random_state=rng)
    rf.fit(X, y)
    permute_scores = permutation_importance(rf, X, y, n_bootstrap=30,
                                            random_state=rng, scoring=scoring)

    assert permute_scores.shape == (X.shape[1], 30)

    permuate_score_means = np.mean(permute_scores, axis=-1)
    assert np.all(permuate_score_means[-1] > permuate_score_means[:-1])
