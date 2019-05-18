import pytest
import numpy as np

from sklearn.datasets import load_boston
from sklearn.datasets import load_iris
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.utils.testing import assert_array_almost_equal


@pytest.mark.parametrize("convert_to_df", [True, False])
@pytest.mark.parametrize("load_dataset,clf", [
    (load_boston, RandomForestRegressor(n_estimators=10, random_state=42)),
    (load_iris, RandomForestClassifier(n_estimators=10, random_state=42))
])
def test_permutation_importance_correlated_feature_regression(
        convert_to_df, load_dataset, clf):
    # Make sure that feature highly correlated to the target have a higher
    # importance
    rng = np.random.RandomState(42)
    n_rounds = 10

    dataset = load_dataset()
    X, y = dataset.data, dataset.target
    y_with_little_noise = (
        y + rng.normal(scale=0.001, size=y.shape[0])).reshape(-1, 1)

    # Adds feature correlated as the last column with target
    if convert_to_df:
        pd = pytest.importorskip("pandas")
        X = pd.DataFrame(X, columns=dataset.feature_names)
        X['correlated_feature'] = y_with_little_noise
    else:
        X = np.hstack([X, y_with_little_noise])

    clf.fit(X, y)

    X_before = X.copy()
    permute_scores = permutation_importance(clf, X, y, n_rounds=n_rounds,
                                            random_state=rng)

    assert permute_scores.shape == (X.shape[1], n_rounds)

    # permutation_importance does not change X
    assert_array_almost_equal(X_before, X)

    permute_score_means = np.mean(permute_scores, axis=-1)

    # the correlated feature was added as the last column and should
    # have the highest importance
    assert np.all(permute_score_means[-1] > permute_score_means[:-1])
