import pytest
import numpy as np

from sklearn.compose import ColumnTransformer
from sklearn.datasets import load_boston
from sklearn.datasets import load_iris
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import StandardScaler


def test_permutation_importance_correlated_feature_regression():
    # Make sure that feature highly correlated to the target have a higher
    # importance
    rng = np.random.RandomState(42)
    n_rounds = 5

    dataset = load_boston()
    X, y = dataset.data, dataset.target
    y_with_little_noise = (
        y + rng.normal(scale=0.001, size=y.shape[0])).reshape(-1, 1)

    X = np.hstack([X, y_with_little_noise])

    clf = RandomForestRegressor(n_estimators=10, random_state=42)
    clf.fit(X, y)

    permute_imp = permutation_importance(clf, X, y, n_rounds=n_rounds,
                                         random_state=rng)

    assert permute_imp.shape == (X.shape[1], n_rounds)
    permute_score_means = np.mean(permute_imp, axis=-1)

    # the correlated feature with y was added as the last column and should
    # have the highest importance
    assert np.all(permute_score_means[-1] > permute_score_means[:-1])


def test_permutation_importance_correlated_feature_regression_pandas():
    pd = pytest.importorskip("pandas")

    # Make sure that feature highly correlated to the target have a higher
    # importance
    rng = np.random.RandomState(42)
    n_rounds = 5

    dataset = load_iris()
    X, y = dataset.data, dataset.target
    y_with_little_noise = (
        y + rng.normal(scale=0.001, size=y.shape[0])).reshape(-1, 1)

    # Adds feature correlated with y as the last column
    X = pd.DataFrame(X, columns=dataset.feature_names)
    X['correlated_feature'] = y_with_little_noise

    clf = RandomForestClassifier(n_estimators=10, random_state=42)
    clf.fit(X, y)

    permute_imp = permutation_importance(clf, X, y, n_rounds=n_rounds,
                                         random_state=rng)

    assert permute_imp.shape == (X.shape[1], n_rounds)
    permute_score_means = np.mean(permute_imp, axis=-1)

    # the correlated feature with y was added as the last column and should
    # have the highest importance
    assert np.all(permute_score_means[-1] > permute_score_means[:-1])


def test_permutation_importance_mixed_types():
    rng = np.random.RandomState(42)
    n_rounds = 3

    # Last column is correlated with y
    X = np.array([[1.0, 2.0, 3.0, np.nan], [2, 1, 2, 1]]).T
    y = np.array([0, 1, 0, 1])

    clf = make_pipeline(SimpleImputer(),
                        LogisticRegression(solver='lbfgs'))
    clf.fit(X, y)
    permute_imp = permutation_importance(clf, X, y, n_rounds=n_rounds,
                                         random_state=rng)

    assert permute_imp.shape == (X.shape[1], n_rounds)
    permute_score_means = np.mean(permute_imp, axis=-1)

    # the correlated feature with y is the last column and should
    # have the highest importance
    assert np.all(permute_score_means[-1] > permute_score_means[:-1])


def test_permutation_importance_mixed_types_pandas():
    pd = pytest.importorskip("pandas")
    rng = np.random.RandomState(42)
    n_rounds = 5

    # Last column is correlated with y
    X = pd.DataFrame({'col1': [1.0, 2.0, 3.0, np.nan],
                      'col2': ['a', 'b', 'a', 'b']})
    y = np.array([0, 1, 0, 1])

    print(X)
    num_preprocess = make_pipeline(SimpleImputer(), StandardScaler())
    preprocess = ColumnTransformer([
        ('num', num_preprocess, ['col1']),
        ('cat', OneHotEncoder(), ['col2'])
    ])
    clf = make_pipeline(preprocess, LogisticRegression(solver='lbfgs'))
    clf.fit(X, y)

    permute_imp = permutation_importance(clf, X, y, n_rounds=n_rounds,
                                         random_state=rng)

    assert permute_imp.shape == (X.shape[1], n_rounds)
    permute_score_means = np.mean(permute_imp, axis=-1)
    # the correlated feature with y is the last column and should
    # have the highest importance
    assert np.all(permute_score_means[-1] > permute_score_means[:-1])
