import pytest
import numpy as np
from numpy.testing import assert_allclose

from sklearn.preprocessing import TargetRegressorEncoder


@pytest.mark.parametrize('categories', [
        np.array([0, 1, 2], dtype=int),
        np.array(['cat', 'dog', 'snake'], dtype=object)
])
@pytest.mark.parametrize('seed', range(3))
def test_regression(categories, seed):
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

    assert len(enc.encodings_) == 1
    assert enc.encoding_mean_ == pytest.approx(y_mean)
    assert_allclose(enc.encodings_[0], cat_encoded)

    expected_encoding = np.take(cat_encoded, X_int[shuffled_idx, :])
    X_trans = enc.transform(X_input)
    assert_allclose(expected_encoding, X_trans)

    # test on unknown category
    X_trans = enc.transform(np.array([[5]], dtype=categories.dtype))
    assert_allclose(X_trans, [[y_mean]])


@pytest.mark.parametrize('categories', [
    np.array([0, 1, 2], dtype=int),
    np.array(['cat', 'dog', 'snake'], dtype=object),
])
def test_zero_variance_category(categories):
    # When the target is constant for a given category, the category encoding
    # should correspond to the that constant target values
    X_int = np.array([[0] * 20 + [1] * 30 + [2] * 40], dtype=int).T
    X_input = categories[X_int]

    # The target of the first category are constant and have no variance
    y = np.array([10] * 20 + [-4] * 15 + [9] * 15 + [-6] * 30 + [25] * 10)

    enc = TargetRegressorEncoder().fit(X_input, y)
    X_test = np.array([[0]], dtype=int).T
    X_input = categories[X_test]
    X_trans = enc.transform(X_input)
    assert_allclose(X_trans, [[10]])


@pytest.mark.parametrize('categories', [
    np.array([0, 1, 2], dtype=int),
    np.array(['cat', 'dog', 'snake'], dtype=object),
])
def test_zero_variance_target(categories):
    # if the target has zero variance, then the mean of the target is used.
    X_int = np.array([[0] * 20 + [1] * 30 + [2] * 40], dtype=int).T
    X_input = categories[X_int]
    n_samples = X_input.shape[0]

    y = np.ones(n_samples)

    enc = TargetRegressorEncoder()
    X_trans = enc.fit_transform(X_input, y)
    expected_trans = np.full((n_samples, 1), fill_value=y.mean(), dtype=float)
    assert_allclose(X_trans, expected_trans)


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
    # Test custom categoires with unknown categories
    rng = np.random.RandomState(42)
    y = rng.uniform(low=-10, high=20, size=X.shape[0])

    enc = TargetRegressorEncoder(categories=categories)
    X_trans = enc.fit_transform(X, y)

    # The last element is unknown
    assert_allclose(X_trans[-1], [y.mean()])

    assert len(enc.encodings_) == 1
    # known category that is unseen during fit time is mapped to the mean
    assert enc.encodings_[0][-1] == pytest.approx(y.mean())


@pytest.mark.parametrize('to_pandas', [True, False])
def test_multiple_features_sanity(to_pandas):
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

    X_test = np.array([
        [0, 1],
        [1, 0],
        [2, 10],  # unknown
    ], dtype=int)

    if to_pandas:
        pd = pytest.importorskip('pandas')
        # convert second feature to a object
        X_obj = np.array(['cat', 'dog'], dtype=object)[X[:, 1]]
        X = pd.DataFrame({
            'feat0': X[:, 0], 'feat1': X_obj}, columns=['feat0', 'feat1']
        )
        X_test = pd.DataFrame({
            'feat0': X_test[:, 0], 'feat1': ['dog', 'cat', 'snake']
        })

    # manually compute multilevel partial pooling
    feat_0_cat_0_encoding = ((4 * 4. / 5. + 4. / 9.5) /
                             (4 / 5. + 1 / 9.5))
    feat_0_cat_1_encoding = ((4 * 4. / 14. + 4. / 9.5) /
                             (4 / 14. + 1 / 9.5))

    feat_1_cat_0_encoding = ((3 * 7. / 6. + 4. / 9.5) /
                             (3 / 6. + 1 / 9.5))
    feat_1_cat_1_encoding = ((5 * 2.2 / 2.96 + 4. / 9.5) /
                             (5 / 2.96 + 1 / 9.5))

    expected_encoding = [
        [feat_0_cat_0_encoding, feat_0_cat_1_encoding],
        [feat_1_cat_0_encoding, feat_1_cat_1_encoding]
    ]

    enc = TargetRegressorEncoder().fit(X, y)
    assert_allclose(expected_encoding, enc.encodings_)
    assert enc.encoding_mean_ == pytest.approx(y_mean)


    X_trans = enc.transform(X_test)
    X_trans_expected = np.array([
       [feat_0_cat_0_encoding, feat_1_cat_1_encoding],
       [feat_0_cat_1_encoding, feat_1_cat_0_encoding],
       [y_mean, y_mean],  # unknown maps to y_mean
    ])
    assert_allclose(X_trans, X_trans_expected)
