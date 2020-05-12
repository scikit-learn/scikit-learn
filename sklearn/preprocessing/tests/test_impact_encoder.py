import pytest
import numpy as np
from numpy.testing import assert_allclose

from sklearn.preprocessing import ImpactEncoder


@pytest.mark.parametrize('categories, unknown_X, seed', [
    (np.array([0, 1, 2], dtype=int),
     np.array([4], dtype=int), 0),
    (np.array(['cat', 'dog', 'snake'], dtype=object),
     np.array(['cow'], dtype=object), 1),
])
def test_regression(categories, unknown_X, seed):
    # checks impact encoder for regression

    X_int = np.array([[0] * 20 + [1] * 30 + [2] * 40], dtype=int).T
    X_input = categories[X_int]
    n_samples = X_input.shape[0]
    rng = np.random.RandomState(seed)
    y = rng.uniform(low=-10, high=20, size=n_samples)

    # bayes smoothing
    # y = lambda * (E[y|category] - E[y])
    # where lambda = count(category) / (count(category) + m)
    # and m = n_samples / cardinality(feature)
    cat_means = np.array([np.mean(y[:20]), np.mean(y[20:50]), np.mean(y[50:])])
    cardinality = 3
    cat_counts = np.array([20, 30, 40])
    lambdas = cat_counts / (cat_counts + n_samples / cardinality)
    cat_encoded = lambdas * (cat_means - np.mean(y))

    # shuffle
    shuffled_idx = rng.permutation(n_samples)
    X_input = X_input[shuffled_idx]
    y = y[shuffled_idx]

    enc = ImpactEncoder().fit(X_input, y)

    expected_encoding = cat_encoded[X_int]
    X_trans = enc.transform(X_input)
    assert_allclose(expected_encoding, X_trans)

    # check known test data
    X_test = np.array([[2, 0, 1]], dtype=int).T
    X_input = categories[X_test]
    X_trans = enc.transform(X_input)
    expected_encoding = cat_encoded[X_test]
    assert_allclose(expected_encoding, X_trans)

    # unknown
    X_trans = enc.transform(unknown_X)
    assert_allclose(X_trans, [[0]])


@pytest.mark.parametrize("X, categories", [
    (np.array([[0] * 10 + [1] * 10 + [3]],
              dtype=int).T,  # 3 is unknown
     [[0, 1, 2]]),
    (np.array([['cat'] * 3 + ['dog'] * 4, ['snake']],
              dtype=object).T,  # snake is unknown
     [['dog', 'cat', 'cow']]),
])
def test_custom_categories(X, categories):
    # Test custom categoires with unknow categories
    rng = np.random.RandomState(42)
    y = rng.uniform(low=-10, high=20, size=X.shape[0])

    X_trans = ImpactEncoder(categories=categories).fit_transform(X, y)
    # The last element is unknown
    assert_allclose(X_trans[-1, 0], 0)
