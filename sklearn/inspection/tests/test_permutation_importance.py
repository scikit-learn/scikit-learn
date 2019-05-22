import pytest
import numpy as np

from sklearn.compose import ColumnTransformer
from sklearn.datasets import load_boston
from sklearn.datasets import load_iris
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.inspection import permutation_importance
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
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

    # Adds feature correlated with y as the last column
    if convert_to_df:
        pd = pytest.importorskip("pandas")
        X = pd.DataFrame(X, columns=dataset.feature_names)
        X['correlated_feature'] = y_with_little_noise
    else:
        X = np.hstack([X, y_with_little_noise])

    clf.fit(X, y)

    X_before = X.copy()
    permute_imp = permutation_importance(clf, X, y, n_rounds=n_rounds,
                                         random_state=rng)

    assert permute_imp.shape == (X.shape[1], n_rounds)

    # permutation_importance does not change X
    assert_array_almost_equal(X_before, X)

    permute_score_means = np.mean(permute_imp, axis=-1)

    # the correlated feature with y was added as the last column and should
    # have the highest importance
    assert np.all(permute_score_means[-1] > permute_score_means[:-1])


@pytest.mark.parametrize("convert_to_df", [True, False])
def test_permutation_importance_mixed_types(convert_to_df):
    rng = np.random.RandomState(42)
    n_rounds = 10

    # Last column is correlated with y
    X = np.array([[1, 2, 3, 4], ['a', 'b', 'a', 'b']]).T
    y = np.array([0, 1, 0, 1])

    if convert_to_df:
        pd = pytest.importorskip("pandas")
        X = pd.DataFrame(X, columns=['num_col', 'cat_col'])
        X['num_col'] = X['num_col'].astype(int)
        preprocess = ColumnTransformer([
            ('num', StandardScaler(), ['num_col']),
            ('cat', OneHotEncoder(), ['cat_col'])
        ])
        clf = make_pipeline(preprocess, LogisticRegression(solver='lbfgs'))
    else:
        clf = make_pipeline(OneHotEncoder(),
                            LogisticRegression(solver='lbfgs'))

    clf.fit(X, y)

    X_before = X.copy()
    permute_imp = permutation_importance(clf, X, y, n_rounds=n_rounds,
                                         random_state=rng)

    assert permute_imp.shape == (X.shape[1], n_rounds)

    # permutation_importance does not change X
    assert np.all(X_before == X)

    permute_score_means = np.mean(permute_imp, axis=-1)

    # the correlated feature with y is the last column and should
    # have the highest importance
    assert np.all(permute_score_means[-1] > permute_score_means[:-1])
