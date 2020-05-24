import pytest
import numpy as np
from numpy.testing import assert_allclose

from sklearn.preprocessing import TargetRegressorEncoder


@pytest.mark.parametrize('categories, unknown_X', [
    (
        np.array([0, 1, 2], dtype=int),
        np.array([4], dtype=int)
    ),
    (
        np.array(['cat', 'dog', 'snake'], dtype=object),
        np.array(['cow'], dtype=object)
    ),
])
@pytest.mark.parametrize('seed', range(3))
def test_regression(categories, unknown_X, seed):
    # checks impact encoder for regression

    X_int = np.array([[0] * 20 + [1] * 30 + [2] * 40], dtype=int).T
    X_input = categories[X_int]
    n_samples = X_input.shape[0]
    rng = np.random.RandomState(seed)
    y = rng.uniform(low=-10, high=20, size=n_samples)

    # multilevel partial pooling
    y_sections = [y[:20], y[20:50], y[50:]]
    cat_means = np.array([np.mean(sect) for sect in y_sections])
    cat_variance = np.array([np.var(sect) for sect in y_sections])
    cat_counts = np.array([20, 30, 40])

    y_variance = np.var(y)
    y_mean = np.mean(y)

    # multilevel partial pooling directly
    cat_encoded = ((cat_counts * cat_means / cat_variance +
                    y_mean / y_variance))
    cat_encoded /= (cat_counts / cat_variance + 1 / y_variance)

    # shuffle
    shuffled_idx = rng.permutation(n_samples)
    X_input = X_input[shuffled_idx]
    y = y[shuffled_idx]

    enc = TargetRegressorEncoder().fit(X_input, y)

    assert len(enc.cat_encodings_) == 1
    assert enc.y_mean_ == pytest.approx(y_mean)
    assert_allclose(enc.cat_encodings_[0], cat_encoded)

    expected_encoding = np.take(cat_encoded, X_int[shuffled_idx, :])
    X_trans = enc.transform(X_input)
    assert_allclose(expected_encoding, X_trans)

    # check known test data
    X_test = np.array([[2, 0, 1]], dtype=int).T
    X_input = categories[X_test]
    X_trans = enc.transform(X_input)
    expected_encoding = cat_encoded[X_test]
    assert_allclose(expected_encoding, X_trans)

    # unknown
    X_trans = enc.transform([unknown_X])
    assert_allclose(X_trans, [[y_mean]])


@pytest.mark.parametrize('categories, unknown_X', [
    (
        np.array([0, 1, 2], dtype=int),
        np.array([4], dtype=int)
    ),
    (
        np.array(['cat', 'dog', 'snake'], dtype=object),
        np.array(['cow'], dtype=object)
    ),
])
def test_zero_variance_group(categories, unknown_X):
    # The first group is constant results in using the mean for the group
    # as the encoding
    X_int = np.array([[0] * 20 + [1] * 30 + [2] * 40], dtype=int).T
    X_input = categories[X_int]

    # The first group is a constant (10) and have no variance
    y = np.array([10] * 20 + [-4] * 15 + [9] * 15 + [-6] * 30 + [25] * 10)

    enc = TargetRegressorEncoder().fit(X_input, y)
    X_test = np.array([[0]], dtype=int).T
    X_input = categories[X_test]
    X_trans = enc.transform(X_input)
    assert_allclose(X_trans, [[10]])

    # unknown
    X_trans = enc.transform([unknown_X])
    assert_allclose(X_trans, [[y.mean()]])


@pytest.mark.parametrize('categories, unknown_X', [
    (
        np.array([0, 1, 2], dtype=int),
        np.array([4], dtype=int)
    ),
    (
        np.array(['cat', 'dog', 'snake'], dtype=object),
        np.array(['cow'], dtype=object)
    ),
])
def test_zero_variance_target(categories, unknown_X):
    # if the target has zero variance, then the mean of the target is used.
    X_int = np.array([[0] * 20 + [1] * 30 + [2] * 40], dtype=int).T
    X_input = categories[X_int]
    n_samples = X_input.shape[0]

    y = np.ones(n_samples)

    enc = TargetRegressorEncoder()
    X_trans = enc.fit_transform(X_input, y)
    expected_trans = np.full((n_samples, 1), fill_value=y.mean(), dtype=float)
    assert_allclose(X_trans, expected_trans)

    # unknown
    X_trans = enc.transform([unknown_X])
    assert_allclose(X_trans, [[y.mean()]])


@pytest.mark.parametrize("X, categories", [
    (
        np.array([[0] * 10 + [1] * 10 + [3]],
                 dtype=int).T,  # 3 is unknown
        [[0, 1, 2]]),
    (
        np.array([['cat'] * 3 + ['dog'] * 4 + ['snake']],
                 dtype=object).T,  # snake is unknown
        [['dog', 'cat', 'cow']]
    ),
])
def test_custom_categories(X, categories):
    # Test custom categoires with unknow categories
    rng = np.random.RandomState(42)
    y = rng.uniform(low=-10, high=20, size=X.shape[0])

    enc = TargetRegressorEncoder(categories=categories)
    X_trans = enc.fit_transform(X, y)

    # The last element is unknown
    assert_allclose(X_trans[-1], [y.mean()])

    assert len(enc.cat_encodings_) == 1
    # unknown category seen during fitting is mapped to the mean
    assert enc.cat_encodings_[0][-1] == pytest.approx(y.mean())


def test_multiple_categories_santiy():
    X = np.array([
        [1, 1],
        [0, 1],
        [1, 1],
        [0, 1],
        [1, 0],
        [0, 1],
        [1, 0],
        [0, 0]
    ], dtype=int)
    y = np.array([0, 1, 2, 3, 4, 5, 10, 7])
    y_mean = np.mean(y)

    # manually compute multilevel partial pooling
    feat_0_cat_0_encoding = ((4 * 4.0 / 5.0 + 4.0 / 9.5) /
                             (4 / 5.0 + 1 / 9.5))
    feat_0_cat_1_encoding = ((4 * 4.0 / 14.0 + 4.0 / 9.5) /
                             (4 / 14.0 + 1 / 9.5))

    feat_1_cat_0_encoding = ((3 * 7.0 / 6.0 + 4.0 / 9.5) /
                             (3 / 6.0 + 1 / 9.5))
    feat_1_cat_1_encoding = ((5 * 2.2 / 2.96 + 4.0 / 9.5) /
                             (5 / 2.96 + 1 / 9.5))

    expected_encoding = [
        [feat_0_cat_0_encoding, feat_0_cat_1_encoding],
        [feat_1_cat_0_encoding, feat_1_cat_1_encoding]
    ]

    enc = TargetRegressorEncoder().fit(X, y)
    assert_allclose(expected_encoding, enc.cat_encodings_)
    assert enc.y_mean_ == pytest.approx(y_mean)

    X_test = np.array([
        [0, 1],
        [1, 0],
        [2, 10],  # unknown
    ], dtype=int)

    X_trans = enc.transform(X_test)
    X_trans_expected = np.array([
       [feat_0_cat_0_encoding, feat_1_cat_1_encoding],
       [feat_0_cat_1_encoding, feat_1_cat_0_encoding],
       [y_mean, y_mean],  # unknown maps to y_mean
    ])
    assert_allclose(X_trans, X_trans_expected)


def test_multiple_categories_santiy_pandas():
    pd = pytest.importorskip('pandas')

    df = pd.DataFrame({
        'feat1': [1, 0, 1, 0, 1, 0],
        'feat2': ['dog', 'dog', 'cow', 'cow', 'snake', 'snake']
    }, columns=['feat1', 'feat2'])
    y = np.array([0, 2, 4, 8, 16, 32])
    y_mean = y.mean()
    y_var = y.var()

    # manually compute multilevel partial pooling
    feat_0_cat_0_encoding = ((3 * 14.0 / 168.0 + y_mean / y_var) /
                             (3 / 168.0 + 1 / y_var))
    feat_0_cat_1_encoding = ((3 * np.mean([0, 4, 16]) / np.var([0, 4, 16]) +
                              y_mean / y_var) /
                             (3 / np.var([0, 4, 16]) + 1 / y_var))

    feat_1_cat_dog_encoding = ((2 * 1.0 / 1.0 + y_mean / y_var) /
                               (2 / 1.0 + 1 / y_var))
    feat_1_cat_cow_encoding = ((2 * 6.0 / 4.0 + y_mean / y_var) /
                               (2 / 4.0 + 1 / y_var))
    feat_1_cat_snake_encoding = ((2 * 24.0 / 64.0 + y_mean / y_var) /
                                 (2 / 64.0 + 1 / y_var))

    expected_encoding = [
        [feat_0_cat_0_encoding, feat_0_cat_1_encoding],
        [feat_1_cat_cow_encoding, feat_1_cat_dog_encoding,
         feat_1_cat_snake_encoding]
    ]

    enc = TargetRegressorEncoder().fit(df, y)
    assert_allclose(enc.cat_encodings_[0], expected_encoding[0])
    assert_allclose(enc.cat_encodings_[1], expected_encoding[1])

    X_test = np.array([
        [1, 'dog'],
        [0, 'cow'],
        [2, 'snake'],
        [10, 'human']
    ], dtype=object)
    X_trans = enc.transform(X_test)

    expected_X_trans = np.array([
      [feat_0_cat_1_encoding, feat_1_cat_dog_encoding],
      [feat_0_cat_0_encoding, feat_1_cat_cow_encoding],
      [y_mean, feat_1_cat_snake_encoding],
      [y_mean, y_mean]
    ])
    assert_allclose(X_trans, expected_X_trans)
