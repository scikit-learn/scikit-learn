import numpy as np
from numpy.testing import assert_allclose
from numpy.testing import assert_array_equal
import pytest

from sklearn.preprocessing import TargetEncoder, LabelEncoder
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
@pytest.mark.parametrize("target_type", ["binary", "continuous"])
def test_encoding(categories, unknown_value, seed, smooth, target_type):
    """Check encoding for binary and continuous targets."""

    X_int = np.array([[0] * 20 + [1] * 30 + [2] * 40], dtype=np.int64).T
    n_categories = 3
    n_samples = X_int.shape[0]

    if isinstance(categories, str) and categories == "auto":
        X_input = X_int
    else:
        X_input = categories[0][X_int]

    rng = np.random.RandomState(seed)

    if target_type == "binary":
        y_num = rng.randint(low=0, high=2, size=n_samples)
        target_names = np.asarray(["cat", "dog"], dtype=object)
        y_input = target_names[y_num]
    else:  # target_type == continuous
        y_num = rng.uniform(low=-10, high=20, size=n_samples)
        y_input = y_num

    # compute encodings for all data to validate `transform`
    y_mean = np.mean(y_num)
    smooth_sum = smooth * y_mean

    expected_encoding_0 = (np.sum(y_num[:20]) + smooth_sum) / (20.0 + smooth)
    expected_encoding_1 = (np.sum(y_num[20:50]) + smooth_sum) / (30.0 + smooth)
    expected_encoding_2 = (np.sum(y_num[50:]) + smooth_sum) / (40.0 + smooth)
    expected_encodings = np.asarray(
        [expected_encoding_0, expected_encoding_1, expected_encoding_2]
    )

    shuffled_idx = rng.permutation(n_samples)
    X_int = X_int[shuffled_idx]
    X_input = X_input[shuffled_idx]
    y_input = y_input[shuffled_idx]
    y_num = y_num[shuffled_idx]

    # Get encodings for cv splits to validate `fit_transform`
    expected_X_fit_transform = np.empty_like(X_int, dtype=np.float64)
    kfold = KFold(n_splits=3)
    for train_idx, test_idx in kfold.split(X_input):
        cur_encodings = np.zeros(n_categories, dtype=np.float64)
        X_, y_ = X_int[train_idx, :], y_num[train_idx]
        y_train_mean = np.mean(y_)
        for c in range(n_categories):
            y_subset = y_[X_[:, 0] == c]
            current_sum = np.sum(y_subset) + y_train_mean * smooth
            current_cnt = y_subset.shape[0] + smooth
            cur_encodings[c] = current_sum / current_cnt

        expected_X_fit_transform[test_idx, 0] = cur_encodings[X_int[test_idx, 0]]

    target_encoder = TargetEncoder(smooth=smooth, categories=categories, cv=kfold)

    X_fit_transform = target_encoder.fit_transform(X_input, y_input)

    assert target_encoder.target_type_ == target_type
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
def test_custom_categories(X, categories):
    """custom categoires with unknown categories that are not in training data."""
    rng = np.random.RandomState(42)
    y = rng.uniform(low=-10, high=20, size=X.shape[0])
    enc = TargetEncoder(categories=categories).fit(X, y)

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
        (
            np.asarray([[1, 2, 0], [1, 2, 3]]).T,
            "Target type was inferred to be 'multiclass-multioutput'",
        ),
        (["cat", "dog", "bear"], "Target type was inferred to be 'multiclass'"),
    ],
)
def test_errors(y, msg):
    """Check invalidate input."""
    X = np.asarray([[1, 0, 1]]).T

    enc = TargetEncoder()
    with pytest.raises(ValueError, match=msg):
        enc.fit_transform(X, y)


def test_feature_names_out_set_output():
    """Check TargetEncoder works with set_output."""
    pd = pytest.importorskip("pandas")

    X_df = pd.DataFrame({"A": ["a", "b"] * 10, "B": [1, 2] * 10})
    y = [1, 2] * 10

    enc_default = TargetEncoder(cv=2, smooth=3.0).set_output(transform="default")
    enc_pandas = TargetEncoder(cv=2, smooth=3.0).set_output(transform="pandas")

    X_default = enc_default.fit_transform(X_df, y)
    X_pandas = enc_pandas.fit_transform(X_df, y)

    assert_allclose(X_pandas.to_numpy(), X_default)
    assert_array_equal(enc_pandas.get_feature_names_out(), X_pandas.columns)


@pytest.mark.parametrize("to_pandas", [True, False])
@pytest.mark.parametrize("smooth", [1.0, 2.0])
@pytest.mark.parametrize("target_type", ["binary-ints", "binary-str", "continuous"])
def test_multiple_features_quick(to_pandas, smooth, target_type):
    """Check target encoder with multiple features."""
    X_int = np.array(
        [[1, 1], [0, 1], [1, 1], [2, 1], [1, 0], [0, 1], [1, 0], [0, 0]], dtype=np.int64
    )
    if target_type == "binary-str":
        y_input = np.array(["a", "b", "a", "a", "b", "b", "a", "b"])
        y_num = LabelEncoder().fit_transform(y_input)
    elif target_type == "binary-ints":
        y_input = np.array([3, 4, 3, 3, 3, 4, 4, 4])
        y_num = LabelEncoder().fit_transform(y_input)
    else:
        y_input = np.array([3.0, 5.1, 2.4, 3.5, 4.1, 5.5, 10.3, 7.3])
        y_num = y_input
    y_mean = np.mean(y_num)
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
            X_, y_ = X_int[train_idx, :], y_num[train_idx]
            y_train_mean = np.mean(y_)
            for c in range(n_cats):
                y_subset = y_[X_[:, f_idx] == c]
                current_sum = np.sum(y_subset) + smooth * y_train_mean
                current_cnt = y_subset.shape[0] + smooth
                current_encoding[c] = current_sum / current_cnt

            expected_X_fit_transform[test_idx, f_idx] = current_encoding[
                X_int[test_idx, f_idx]
            ]

    # manually compute encoding for transform
    expected_encodings = []
    for f_idx, cats in enumerate(categories):
        n_cats = len(cats)
        current_encoding = np.zeros(n_cats, dtype=np.float64)
        for c in range(n_cats):
            y_subset = y_num[X_int[:, f_idx] == c]
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

    enc = TargetEncoder(smooth=smooth, cv=cv_splits)
    X_fit_transform = enc.fit_transform(X_input, y_input)
    assert_allclose(X_fit_transform, expected_X_fit_transform)

    assert len(enc.encodings_) == 2
    for i in range(2):
        assert_allclose(enc.encodings_[i], expected_encodings[i])

    X_test_transform = enc.transform(X_test)
    assert_allclose(X_test_transform, expected_X_test_transform)
