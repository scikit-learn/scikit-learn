import numpy as np
from numpy.testing import assert_allclose
from numpy.testing import assert_array_equal
import pytest

from sklearn.preprocessing import TargetRegressorEncoder
from sklearn.model_selection import KFold


@pytest.mark.parametrize(
    "categories, unknown_value",
    [
        ([np.array([0, 1, 2], dtype=np.int64)], 4),
        ([np.array([1.0, 3.0, np.nan], dtype=np.float64)], 6.0),
        ([np.array(["cat", "dog", "snake"], dtype=object)], "bear"),
        ("auto", 3),
    ],
)
@pytest.mark.parametrize("seed", range(2))
@pytest.mark.parametrize("smooth", [5.0, 10.0])
def test_regression(categories, unknown_value, seed, smooth):
    """Check regression encoding."""

    X_int = np.array([[0] * 20 + [1] * 30 + [2] * 40], dtype=np.int64).T
    n_categories = 3
    n_samples = X_int.shape[0]

    if isinstance(categories, str) and categories == "auto":
        X_input = X_int
    else:
        X_input = categories[0][X_int]

    rng = np.random.RandomState(seed)
    y = rng.uniform(low=-10, high=20, size=n_samples)

    # compute encodings for all data to validate `transform`
    y_mean = np.mean(y)
    smooth_sum = smooth * y_mean

    expected_encoding_0 = (np.sum(y[:20]) + smooth_sum) / (20.0 + smooth)
    expected_encoding_1 = (np.sum(y[20:50]) + smooth_sum) / (30.0 + smooth)
    expected_encoding_2 = (np.sum(y[50:]) + smooth_sum) / (40.0 + smooth)
    expected_encodings = np.asarray(
        [expected_encoding_0, expected_encoding_1, expected_encoding_2]
    )

    shuffled_idx = rng.permutation(n_samples)
    X_int = X_int[shuffled_idx]
    X_input = X_input[shuffled_idx]
    y = y[shuffled_idx]

    # Get encodings for cv splits to validate `fit_transform`
    expected_X_fit_transform = np.empty_like(X_int, dtype=np.float64)
    kfold = KFold(n_splits=3)
    for train_idx, test_idx in kfold.split(X_input):
        cur_encodings = np.zeros(n_categories, dtype=np.float64)
        X_, y_ = X_int[train_idx, :], y[train_idx]
        y_train_mean = np.mean(y_)
        for c in range(n_categories):
            y_subset = y_[X_[:, 0] == c]
            current_sum = np.sum(y_subset) + y_train_mean * smooth
            current_cnt = y_subset.shape[0] + smooth
            cur_encodings[c] = current_sum / current_cnt

        expected_X_fit_transform[test_idx, 0] = np.take(
            cur_encodings, indices=X_int[test_idx, 0]
        )

    target_encoder = TargetRegressorEncoder(
        smooth=smooth, categories=categories, cv=kfold
    )

    X_fit_transform = target_encoder.fit_transform(X_input, y)
    assert_allclose(X_fit_transform, expected_X_fit_transform)
    assert len(target_encoder.encodings_) == 1
    assert_allclose(target_encoder.encodings_[0], expected_encodings)
    assert target_encoder.encoding_mean_ == pytest.approx(y_mean)

    X_test = np.array([[0, 1, 2]], dtype=np.int64).T
    expected_X_test_transform = np.concatenate(
        (expected_encodings, np.asarray([y_mean]))
    ).reshape(-1, 1)

    if isinstance(categories, str) and categories == "auto":
        X_test_input = X_test
    else:
        X_test_input = categories[0][X_test]

    X_test_input = np.concatenate((X_test_input, [[unknown_value]]))

    X_test_transform = target_encoder.transform(X_test_input)
    assert_allclose(X_test_transform, expected_X_test_transform)


@pytest.mark.parametrize(
    "X, categories",
    [
        (
            np.array([[0] * 10 + [1] * 10 + [3]], dtype=np.int64).T,  # 3 is unknown
            [[0, 1, 2]],
        ),
        (
            np.array(
                [["cat"] * 10 + ["dog"] * 10 + ["snake"]], dtype=object
            ).T,  # snake is unknown
            [["dog", "cat", "cow"]],
        ),
    ],
)
def test_regression_custom_categories(X, categories):
    """custom categoires with unknown categories that are not in training data."""
    rng = np.random.RandomState(42)
    y = rng.uniform(low=-10, high=20, size=X.shape[0])
    enc = TargetRegressorEncoder(categories=categories).fit(X, y)

    # The last element is unknown and encoded as the mean
    y_mean = y.mean()
    X_trans = enc.transform(X[-1:])
    assert X_trans[0, 0] == pytest.approx(y_mean)

    assert len(enc.encodings_) == 1
    # custom category that is not in training data
    assert enc.encodings_[0][-1] == pytest.approx(y_mean)


@pytest.mark.parametrize(
    "y, msg",
    [
        ([1, 2, 0, 1], "Found input variables with inconsistent"),
        (np.asarray([[1, 2, 0], [1, 2, 3]]).T, "y should be a 1d array"),
        (["cat", "dog", "bear"], "could not convert string to float"),
    ],
)
def test_regression_errors(y, msg):
    """Check invalidate input."""
    X = np.asarray([[1, 0, 1]]).T

    enc = TargetRegressorEncoder()
    with pytest.raises(ValueError, match=msg):
        enc.fit_transform(X, y)


def test_regression_feature_names_out_set_output():
    """Check TargetRegressorEncoder works with set_output."""
    pd = pytest.importorskip("pandas")

    X_df = pd.DataFrame({"A": ["a", "b"] * 10, "B": [1, 2] * 10})
    y = [1, 2] * 10

    enc_default = TargetRegressorEncoder(cv=2, smooth=3.0).set_output(
        transform="default"
    )
    enc_pandas = TargetRegressorEncoder(cv=2, smooth=3.0).set_output(transform="pandas")

    X_default = enc_default.fit_transform(X_df, y)
    X_pandas = enc_pandas.fit_transform(X_df, y)

    assert_allclose(X_pandas.to_numpy(), X_default)
    assert_array_equal(enc_pandas.get_feature_names_out(), X_pandas.columns)


@pytest.mark.parametrize("to_pandas", [True, False])
@pytest.mark.parametrize("smooth", [1.0, 2.0])
def test_regression_multiple_features_quick(to_pandas, smooth):
    """Check regression encoder with multiple features."""
    X_int = np.array(
        [[1, 1], [0, 1], [1, 1], [2, 1], [1, 0], [0, 1], [1, 0], [0, 0]], dtype=np.int64
    )
    y = np.array([3, 5, 2, 3, 4, 5, 10, 7])
    y_mean = np.mean(y)
    categories = [[0, 1, 2], [0, 1]]

    X_test = np.array(
        [
            [0, 1],
            [3, 0],  # 3 is unknown
            [1, 10],  # 10 is unknown
        ],
        dtype=np.int64,
    )

    if to_pandas:
        pd = pytest.importorskip("pandas")
        # convert second feature to a object
        X_obj = np.array(["cat", "dog"], dtype=object)[X_int[:, 1]]
        X_input = pd.DataFrame(
            {"feat0": X_int[:, 0], "feat1": X_obj}, columns=["feat0", "feat1"]
        )
        X_test = pd.DataFrame({"feat0": X_test[:, 0], "feat1": ["dog", "cat", "snake"]})
    else:
        X_input = X_int

    cv_splits = [([0, 2, 4, 6], [1, 3, 5, 7]), ([1, 3, 5, 7], [0, 2, 4, 6])]

    # manually compute encoding for fit_transform
    expected_X_fit_transform = np.empty_like(X_int, dtype=np.float64)
    for f_idx, cats in enumerate(categories):
        for train_idx, test_idx in cv_splits:
            n_cats = len(cats)
            current_encoding = np.zeros(n_cats, dtype=np.float64)
            X_, y_ = X_int[train_idx, :], y[train_idx]
            y_train_mean = np.mean(y_)
            for c in range(n_cats):
                y_subset = y_[X_[:, f_idx] == c]
                current_sum = np.sum(y_subset) + smooth * y_train_mean
                current_cnt = y_subset.shape[0] + smooth
                current_encoding[c] = current_sum / current_cnt

            expected_X_fit_transform[test_idx, f_idx] = np.take(
                current_encoding, indices=X_int[test_idx, f_idx]
            )

    # manually compute encoding for transform
    expected_encodings = []
    for f_idx, cats in enumerate(categories):
        n_cats = len(cats)
        current_encoding = np.zeros(n_cats, dtype=np.float64)
        for c in range(n_cats):
            y_subset = y[X_int[:, f_idx] == c]
            current_sum = np.sum(y_subset) + smooth * y_mean
            current_cnt = y_subset.shape[0] + smooth
            current_encoding[c] = current_sum / current_cnt
        expected_encodings.append(current_encoding)

    expected_X_test_transform = np.array(
        [
            [expected_encodings[0][0], expected_encodings[1][1]],
            [y_mean, expected_encodings[1][0]],
            [expected_encodings[0][1], y_mean],
        ],
        dtype=np.float64,
    )

    enc = TargetRegressorEncoder(smooth=smooth, cv=cv_splits)
    X_fit_transform = enc.fit_transform(X_input, y)
    assert_allclose(X_fit_transform, expected_X_fit_transform)

    assert len(enc.encodings_) == 2
    for i in range(2):
        assert_allclose(enc.encodings_[i], expected_encodings[i])

    X_test_transform = enc.transform(X_test)
    assert_allclose(X_test_transform, expected_X_test_transform)
