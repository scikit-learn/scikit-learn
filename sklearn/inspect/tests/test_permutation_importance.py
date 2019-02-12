import pytest

import numpy as np

from sklearn.datasets import load_boston
from sklearn.inspect import permutation_importance
from sklearn.ensemble import RandomForestRegressor


@pytest.mark.parametrize("columns", [
    None, [0, 2, 4, 6, 8, 10, 12, 13], [1, 3, 5, 7, 9, 11, 13]
])
@pytest.mark.parametrize("scoring", [
    None, "neg_mean_absolute_error"
])
def test_permutation_importance_correlated_feature_is_important(
        columns, scoring):
    rng = np.random.RandomState(42)
    X, y = load_boston(return_X_y=True)

    # Adds correlated feature to X
    y_with_little_noise = y + rng.normal(scale=0.001, size=y.shape[0])
    X = np.hstack([X, y_with_little_noise.reshape(-1, 1)])

    rf = RandomForestRegressor(n_estimators=50, random_state=42)
    permute_scores = permutation_importance(rf, X, y, columns=columns, cv=4,
                                            random_state=42, scoring=scoring)

    if columns is None:
        assert permute_scores.shape == (X.shape[1], 4)
    else:
        assert permute_scores.shape == (len(columns), 4)

    permuate_score_means = np.mean(permute_scores, axis=-1)
    assert np.all(permuate_score_means[-1] > permuate_score_means[:-1])
