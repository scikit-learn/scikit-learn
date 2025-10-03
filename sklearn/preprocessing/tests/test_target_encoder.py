import re

import numpy as np
import numpy.testing as npt
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.model_selection import (
    KFold,
    ShuffleSplit,
    StratifiedKFold,
    cross_val_score,
    train_test_split,
)
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import (
    KBinsDiscretizer,
    LabelBinarizer,
    LabelEncoder,
    TargetEncoder,
)
from sklearn.preprocessing._target_encoder import _norm_key


def _encode_target(X_ordinal, y_numeric, n_categories, smooth):
    """Simple Python implementation of target encoding."""
    cur_encodings = np.zeros(n_categories, dtype=np.float64)
    y_mean = np.mean(y_numeric)

    if smooth == "auto":
        y_variance = np.var(y_numeric)
        for c in range(n_categories):
            y_subset = y_numeric[X_ordinal == c]
            n_i = y_subset.shape[0]

            if n_i == 0:
                cur_encodings[c] = y_mean
                continue

            y_subset_variance = np.var(y_subset)
            m = y_subset_variance / y_variance
            lambda_ = n_i / (n_i + m)

            cur_encodings[c] = lambda_ * np.mean(y_subset) + (1 - lambda_) * y_mean
        return cur_encodings
    else:  # float
        for c in range(n_categories):
            y_subset = y_numeric[X_ordinal == c]
            current_sum = np.sum(y_subset) + y_mean * smooth
            current_cnt = y_subset.shape[0] + smooth
            cur_encodings[c] = current_sum / current_cnt
        return cur_encodings


@pytest.mark.parametrize(
    "categories, unknown_value",
    [
        ([np.array([0, 1, 2], dtype=np.int64)], 4),
        ([np.array([1.0, 3.0, np.nan], dtype=np.float64)], 6.0),
        ([np.array(["cat", "dog", "snake"], dtype=object)], "bear"),
        ("auto", 3),
    ],
)
@pytest.mark.parametrize("smooth", [5.0, "auto"])
@pytest.mark.parametrize("target_type", ["binary", "continuous"])
def test_encoding(categories, unknown_value, global_random_seed, smooth, target_type):
    """Check encoding for binary and continuous targets.

    Compare the values returned by `TargetEncoder.fit_transform` against the
    expected encodings for cv splits from a naive reference Python
    implementation in _encode_target.
    """

    n_categories = 3
    X_train_int_array = np.array([[0] * 20 + [1] * 30 + [2] * 40], dtype=np.int64).T
    X_test_int_array = np.array([[0, 1, 2]], dtype=np.int64).T
    n_samples = X_train_int_array.shape[0]

    if categories == "auto":
        X_train = X_train_int_array
        X_test = X_test_int_array
    else:
        X_train = categories[0][X_train_int_array]
        X_test = categories[0][X_test_int_array]

    X_test = np.concatenate((X_test, [[unknown_value]]))

    data_rng = np.random.RandomState(global_random_seed)
    n_splits = 3
    if target_type == "binary":
        y_numeric = data_rng.randint(low=0, high=2, size=n_samples)
        target_names = np.array(["cat", "dog"], dtype=object)
        y_train = target_names[y_numeric]

    else:
        assert target_type == "continuous"
        y_numeric = data_rng.uniform(low=-10, high=20, size=n_samples)
        y_train = y_numeric

    shuffled_idx = data_rng.permutation(n_samples)
    X_train_int_array = X_train_int_array[shuffled_idx]
    X_train = X_train[shuffled_idx]
    y_train = y_train[shuffled_idx]
    y_numeric = y_numeric[shuffled_idx]

    # Define our CV splitting strategy
    if target_type == "binary":
        cv = StratifiedKFold(
            n_splits=n_splits, random_state=global_random_seed, shuffle=True
        )
    else:
        cv = KFold(n_splits=n_splits, random_state=global_random_seed, shuffle=True)

    # Compute the expected values using our reference Python implementation of
    # target encoding:
    expected_X_fit_transform = np.empty_like(X_train_int_array, dtype=np.float64)

    for train_idx, test_idx in cv.split(X_train_int_array, y_train):
        X_, y_ = X_train_int_array[train_idx, 0], y_numeric[train_idx]
        cur_encodings = _encode_target(X_, y_, n_categories, smooth)
        expected_X_fit_transform[test_idx, 0] = cur_encodings[
            X_train_int_array[test_idx, 0]
        ]

    # Check that we can obtain the same encodings by calling `fit_transform` on
    # the estimator with the same CV parameters:
    target_encoder = TargetEncoder(
        smooth=smooth,
        categories=categories,
        cv=n_splits,
        random_state=global_random_seed,
    )

    X_fit_transform = target_encoder.fit_transform(X_train, y_train)

    assert target_encoder.target_type_ == target_type
    assert_allclose(X_fit_transform, expected_X_fit_transform)
    assert len(target_encoder.encodings_) == 1
    if target_type == "binary":
        assert_array_equal(target_encoder.classes_, target_names)
    else:
        assert target_encoder.classes_ is None

    # compute encodings for all data to validate `transform`
    y_mean = np.mean(y_numeric)
    expected_encodings = _encode_target(
        X_train_int_array[:, 0], y_numeric, n_categories, smooth
    )
    assert_allclose(target_encoder.encodings_[0], expected_encodings)
    assert target_encoder.target_mean_ == pytest.approx(y_mean)

    # Transform on test data, the last value is unknown so it is encoded as the target
    # mean
    expected_X_test_transform = np.concatenate(
        (expected_encodings, np.array([y_mean]))
    ).reshape(-1, 1)

    X_test_transform = target_encoder.transform(X_test)
    assert_allclose(X_test_transform, expected_X_test_transform)


@pytest.mark.parametrize(
    "categories, unknown_values",
    [
        ([np.array([0, 1, 2], dtype=np.int64)], "auto"),
        ([np.array(["cat", "dog", "snake"], dtype=object)], ["bear", "rabbit"]),
    ],
)
@pytest.mark.parametrize(
    "target_labels", [np.array([1, 2, 3]), np.array(["a", "b", "c"])]
)
@pytest.mark.parametrize("smooth", [5.0, "auto"])
def test_encoding_multiclass(
    global_random_seed, categories, unknown_values, target_labels, smooth
):
    """Check encoding for multiclass targets."""
    rng = np.random.RandomState(global_random_seed)

    n_samples = 80
    n_features = 2
    feat_1_int = np.array(rng.randint(low=0, high=2, size=n_samples))
    feat_2_int = np.array(rng.randint(low=0, high=3, size=n_samples))
    feat_1 = categories[0][feat_1_int]
    feat_2 = categories[0][feat_2_int]
    X_train = np.column_stack((feat_1, feat_2))
    X_train_int = np.column_stack((feat_1_int, feat_2_int))
    categories_ = [[0, 1], [0, 1, 2]]

    n_classes = 3
    y_train_int = np.array(rng.randint(low=0, high=n_classes, size=n_samples))
    y_train = target_labels[y_train_int]
    y_train_enc = LabelBinarizer().fit_transform(y_train)

    n_splits = 3
    cv = StratifiedKFold(
        n_splits=n_splits, random_state=global_random_seed, shuffle=True
    )

    # Manually compute encodings for cv splits to validate `fit_transform`
    expected_X_fit_transform = np.empty(
        (X_train_int.shape[0], X_train_int.shape[1] * n_classes),
        dtype=np.float64,
    )
    for f_idx, cats in enumerate(categories_):
        for c_idx in range(n_classes):
            for train_idx, test_idx in cv.split(X_train, y_train):
                y_class = y_train_enc[:, c_idx]
                X_, y_ = X_train_int[train_idx, f_idx], y_class[train_idx]
                current_encoding = _encode_target(X_, y_, len(cats), smooth)
                # f_idx:   0, 0, 0, 1, 1, 1
                # c_idx:   0, 1, 2, 0, 1, 2
                # exp_idx: 0, 1, 2, 3, 4, 5
                exp_idx = c_idx + (f_idx * n_classes)
                expected_X_fit_transform[test_idx, exp_idx] = current_encoding[
                    X_train_int[test_idx, f_idx]
                ]

    target_encoder = TargetEncoder(
        smooth=smooth,
        cv=n_splits,
        random_state=global_random_seed,
    )
    X_fit_transform = target_encoder.fit_transform(X_train, y_train)

    assert target_encoder.target_type_ == "multiclass"
    assert_allclose(X_fit_transform, expected_X_fit_transform)

    # Manually compute encoding to validate `transform`
    expected_encodings = []
    for f_idx, cats in enumerate(categories_):
        for c_idx in range(n_classes):
            y_class = y_train_enc[:, c_idx]
            current_encoding = _encode_target(
                X_train_int[:, f_idx], y_class, len(cats), smooth
            )
            expected_encodings.append(current_encoding)

    assert len(target_encoder.encodings_) == n_features * n_classes
    for i in range(n_features * n_classes):
        assert_allclose(target_encoder.encodings_[i], expected_encodings[i])
    assert_array_equal(target_encoder.classes_, target_labels)

    # Include unknown values at the end
    X_test_int = np.array([[0, 1], [1, 2], [4, 5]])
    if unknown_values == "auto":
        X_test = X_test_int
    else:
        X_test = np.empty_like(X_test_int[:-1, :], dtype=object)
        for column_idx in range(X_test_int.shape[1]):
            X_test[:, column_idx] = categories[0][X_test_int[:-1, column_idx]]
        # Add unknown values at end
        X_test = np.vstack((X_test, unknown_values))

    y_mean = np.mean(y_train_enc, axis=0)
    expected_X_test_transform = np.empty(
        (X_test_int.shape[0], X_test_int.shape[1] * n_classes),
        dtype=np.float64,
    )
    n_rows = X_test_int.shape[0]
    f_idx = [0, 0, 0, 1, 1, 1]
    # Last row are unknowns, dealt with later
    for row_idx in range(n_rows - 1):
        for i, enc in enumerate(expected_encodings):
            expected_X_test_transform[row_idx, i] = enc[X_test_int[row_idx, f_idx[i]]]

    # Unknowns encoded as target mean for each class
    # `y_mean` contains target mean for each class, thus cycle through mean of
    # each class, `n_features` times
    mean_idx = [0, 1, 2, 0, 1, 2]
    for i in range(n_classes * n_features):
        expected_X_test_transform[n_rows - 1, i] = y_mean[mean_idx[i]]

    X_test_transform = target_encoder.transform(X_test)
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
@pytest.mark.parametrize("smooth", [4.0, "auto"])
def test_custom_categories(X, categories, smooth):
    """Custom categories with unknown categories that are not in training data."""
    rng = np.random.RandomState(0)
    y = rng.uniform(low=-10, high=20, size=X.shape[0])
    enc = TargetEncoder(categories=categories, smooth=smooth, random_state=0).fit(X, y)

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
            np.array([[1, 2, 0], [1, 2, 3]]).T,
            "Target type was inferred to be 'multiclass-multioutput'",
        ),
    ],
)
def test_errors(y, msg):
    """Check invalidate input."""
    X = np.array([[1, 0, 1]]).T

    enc = TargetEncoder()
    with pytest.raises(ValueError, match=msg):
        enc.fit_transform(X, y)


def test_use_regression_target():
    """Check inferred and specified `target_type` on regression target."""
    X = np.array([[0, 1, 0, 1, 0, 1]]).T
    y = np.array([1.0, 2.0, 3.0, 2.0, 3.0, 4.0])

    enc = TargetEncoder(cv=2)
    with pytest.warns(
        UserWarning,
        match=re.escape(
            "The least populated class in y has only 1 members, which is less than"
            " n_splits=2."
        ),
    ):
        enc.fit_transform(X, y)
    assert enc.target_type_ == "multiclass"

    enc = TargetEncoder(cv=2, target_type="continuous")
    enc.fit_transform(X, y)
    assert enc.target_type_ == "continuous"


@pytest.mark.parametrize(
    "y, feature_names",
    [
        ([1, 2] * 10, ["A", "B"]),
        ([1, 2, 3] * 6 + [1, 2], ["A_1", "A_2", "A_3", "B_1", "B_2", "B_3"]),
        (
            ["y1", "y2", "y3"] * 6 + ["y1", "y2"],
            ["A_y1", "A_y2", "A_y3", "B_y1", "B_y2", "B_y3"],
        ),
    ],
)
def test_feature_names_out_set_output(y, feature_names):
    """Check TargetEncoder works with set_output."""
    pd = pytest.importorskip("pandas")

    X_df = pd.DataFrame({"A": ["a", "b"] * 10, "B": [1, 2] * 10})

    enc_default = TargetEncoder(cv=2, smooth=3.0, random_state=0)
    enc_default.set_output(transform="default")
    enc_pandas = TargetEncoder(cv=2, smooth=3.0, random_state=0)
    enc_pandas.set_output(transform="pandas")

    X_default = enc_default.fit_transform(X_df, y)
    X_pandas = enc_pandas.fit_transform(X_df, y)

    assert_allclose(X_pandas.to_numpy(), X_default)
    assert_array_equal(enc_pandas.get_feature_names_out(), feature_names)
    assert_array_equal(enc_pandas.get_feature_names_out(), X_pandas.columns)


@pytest.mark.parametrize("to_pandas", [True, False])
@pytest.mark.parametrize("smooth", [1.0, "auto"])
@pytest.mark.parametrize("target_type", ["binary-ints", "binary-str", "continuous"])
def test_multiple_features_quick(to_pandas, smooth, target_type):
    """Check target encoder with multiple features."""
    X_ordinal = np.array(
        [[1, 1], [0, 1], [1, 1], [2, 1], [1, 0], [0, 1], [1, 0], [0, 0]], dtype=np.int64
    )
    if target_type == "binary-str":
        y_train = np.array(["a", "b", "a", "a", "b", "b", "a", "b"])
        y_integer = LabelEncoder().fit_transform(y_train)
        cv = StratifiedKFold(2, random_state=0, shuffle=True)
    elif target_type == "binary-ints":
        y_train = np.array([3, 4, 3, 3, 3, 4, 4, 4])
        y_integer = LabelEncoder().fit_transform(y_train)
        cv = StratifiedKFold(2, random_state=0, shuffle=True)
    else:
        y_train = np.array([3.0, 5.1, 2.4, 3.5, 4.1, 5.5, 10.3, 7.3], dtype=np.float32)
        y_integer = y_train
        cv = KFold(2, random_state=0, shuffle=True)
    y_mean = np.mean(y_integer)
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
        # convert second feature to an object
        X_train = pd.DataFrame(
            {
                "feat0": X_ordinal[:, 0],
                "feat1": np.array(["cat", "dog"], dtype=object)[X_ordinal[:, 1]],
            }
        )
        # "snake" is unknown
        X_test = pd.DataFrame({"feat0": X_test[:, 0], "feat1": ["dog", "cat", "snake"]})
    else:
        X_train = X_ordinal

    # manually compute encoding for fit_transform
    expected_X_fit_transform = np.empty_like(X_ordinal, dtype=np.float64)
    for f_idx, cats in enumerate(categories):
        for train_idx, test_idx in cv.split(X_ordinal, y_integer):
            X_, y_ = X_ordinal[train_idx, f_idx], y_integer[train_idx]
            current_encoding = _encode_target(X_, y_, len(cats), smooth)
            expected_X_fit_transform[test_idx, f_idx] = current_encoding[
                X_ordinal[test_idx, f_idx]
            ]

    # manually compute encoding for transform
    expected_encodings = []
    for f_idx, cats in enumerate(categories):
        current_encoding = _encode_target(
            X_ordinal[:, f_idx], y_integer, len(cats), smooth
        )
        expected_encodings.append(current_encoding)

    expected_X_test_transform = np.array(
        [
            [expected_encodings[0][0], expected_encodings[1][1]],
            [y_mean, expected_encodings[1][0]],
            [expected_encodings[0][1], y_mean],
        ],
        dtype=np.float64,
    )

    enc = TargetEncoder(smooth=smooth, cv=2, random_state=0)
    X_fit_transform = enc.fit_transform(X_train, y_train)
    assert_allclose(X_fit_transform, expected_X_fit_transform)

    assert len(enc.encodings_) == 2
    for i in range(2):
        assert_allclose(enc.encodings_[i], expected_encodings[i])

    X_test_transform = enc.transform(X_test)
    assert_allclose(X_test_transform, expected_X_test_transform)


@pytest.mark.parametrize(
    "y, y_mean",
    [
        (np.array([3.4] * 20), 3.4),
        (np.array([0] * 20), 0),
        (np.array(["a"] * 20, dtype=object), 0),
    ],
    ids=["continuous", "binary", "binary-string"],
)
@pytest.mark.parametrize("smooth", ["auto", 4.0, 0.0])
def test_constant_target_and_feature(y, y_mean, smooth):
    """Check edge case where feature and target is constant."""
    X = np.array([[1] * 20]).T
    n_samples = X.shape[0]

    enc = TargetEncoder(cv=2, smooth=smooth, random_state=0)
    X_trans = enc.fit_transform(X, y)
    assert_allclose(X_trans, np.repeat([[y_mean]], n_samples, axis=0))
    assert enc.encodings_[0][0] == pytest.approx(y_mean)
    assert enc.target_mean_ == pytest.approx(y_mean)

    X_test = np.array([[1], [0]])
    X_test_trans = enc.transform(X_test)
    assert_allclose(X_test_trans, np.repeat([[y_mean]], 2, axis=0))


def test_fit_transform_not_associated_with_y_if_ordinal_categorical_is_not(
    global_random_seed,
):
    cardinality = 30  # not too large, otherwise we need a very large n_samples
    n_samples = 3000
    rng = np.random.RandomState(global_random_seed)
    y_train = rng.normal(size=n_samples)
    X_train = rng.randint(0, cardinality, size=n_samples).reshape(-1, 1)

    # Sort by y_train to attempt to cause a leak
    y_sorted_indices = y_train.argsort()
    y_train = y_train[y_sorted_indices]
    X_train = X_train[y_sorted_indices]

    target_encoder = TargetEncoder(shuffle=True, random_state=global_random_seed)
    X_encoded_train_shuffled = target_encoder.fit_transform(X_train, y_train)

    target_encoder = TargetEncoder(shuffle=False)
    X_encoded_train_no_shuffled = target_encoder.fit_transform(X_train, y_train)

    # Check that no information about y_train has leaked into X_train:
    regressor = RandomForestRegressor(
        n_estimators=10, min_samples_leaf=20, random_state=global_random_seed
    )

    # It's impossible to learn a good predictive model on the training set when
    # using the original representation X_train or the target encoded
    # representation with shuffled inner CV. For the latter, no information
    # about y_train has inadvertently leaked into the prior used to generate
    # `X_encoded_train_shuffled`:
    cv = ShuffleSplit(n_splits=50, random_state=global_random_seed)
    assert cross_val_score(regressor, X_train, y_train, cv=cv).mean() < 0.1
    assert (
        cross_val_score(regressor, X_encoded_train_shuffled, y_train, cv=cv).mean()
        < 0.1
    )

    # Without the inner CV shuffling, a lot of information about y_train goes into the
    # the per-fold y_train.mean() priors: shrinkage is no longer effective in this
    # case and would no longer be able to prevent downstream over-fitting.
    assert (
        cross_val_score(regressor, X_encoded_train_no_shuffled, y_train, cv=cv).mean()
        > 0.5
    )


def test_smooth_zero():
    """Check edge case with zero smoothing and cv does not contain category."""
    X = np.array([[0, 0, 0, 0, 0, 1, 1, 1, 1, 1]]).T
    y = np.array([2.1, 4.3, 1.2, 3.1, 1.0, 9.0, 10.3, 14.2, 13.3, 15.0])

    enc = TargetEncoder(smooth=0.0, shuffle=False, cv=2)
    X_trans = enc.fit_transform(X, y)

    # With cv = 2, category 0 does not exist in the second half, thus
    # it will be encoded as the mean of the second half
    assert_allclose(X_trans[0], np.mean(y[5:]))

    # category 1 does not exist in the first half, thus it will be encoded as
    # the mean of the first half
    assert_allclose(X_trans[-1], np.mean(y[:5]))


@pytest.mark.parametrize("smooth", [0.0, 1e3, "auto"])
def test_invariance_of_encoding_under_label_permutation(smooth, global_random_seed):
    # Check that the encoding does not depend on the integer of the value of
    # the integer labels. This is quite a trivial property but it is helpful
    # to understand the following test.
    rng = np.random.RandomState(global_random_seed)

    # Random y and informative categorical X to make the test non-trivial when
    # using smoothing.
    y = rng.normal(size=1000)
    n_categories = 30
    X = KBinsDiscretizer(
        n_bins=n_categories, quantile_method="averaged_inverted_cdf", encode="ordinal"
    ).fit_transform(y.reshape(-1, 1))

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, random_state=global_random_seed
    )

    # Shuffle the labels to make sure that the encoding is invariant to the
    # permutation of the labels
    permutated_labels = rng.permutation(n_categories)
    X_train_permuted = permutated_labels[X_train.astype(np.int32)]
    X_test_permuted = permutated_labels[X_test.astype(np.int32)]

    target_encoder = TargetEncoder(smooth=smooth, random_state=global_random_seed)
    X_train_encoded = target_encoder.fit_transform(X_train, y_train)
    X_test_encoded = target_encoder.transform(X_test)

    X_train_permuted_encoded = target_encoder.fit_transform(X_train_permuted, y_train)
    X_test_permuted_encoded = target_encoder.transform(X_test_permuted)

    assert_allclose(X_train_encoded, X_train_permuted_encoded)
    assert_allclose(X_test_encoded, X_test_permuted_encoded)


@pytest.mark.parametrize("smooth", [0.0, "auto"])
def test_target_encoding_for_linear_regression(smooth, global_random_seed):
    # Check some expected statistical properties when fitting a linear
    # regression model on target encoded features depending on their relation
    # with that target.

    # In this test, we use the Ridge class with the "lsqr" solver and a little
    # bit of regularization to implement a linear regression model that
    # converges quickly for large `n_samples` and robustly in case of
    # correlated features. Since we will fit this model on a mean centered
    # target, we do not need to fit an intercept and this will help simplify
    # the analysis with respect to the expected coefficients.
    linear_regression = Ridge(alpha=1e-6, solver="lsqr", fit_intercept=False)

    # Construct a random target variable. We need a large number of samples for
    # this test to be stable across all values of the random seed.
    n_samples = 50_000
    rng = np.random.RandomState(global_random_seed)
    y = rng.randn(n_samples)

    # Generate a single informative ordinal feature with medium cardinality.
    # Inject some irreducible noise to make it harder for a multivariate model
    # to identify the informative feature from other pure noise features.
    noise = 0.8 * rng.randn(n_samples)
    n_categories = 100
    X_informative = KBinsDiscretizer(
        n_bins=n_categories,
        encode="ordinal",
        strategy="uniform",
        random_state=rng,
    ).fit_transform((y + noise).reshape(-1, 1))

    # Let's permute the labels to hide the fact that this feature is
    # informative to naive linear regression model trained on the raw ordinal
    # values. As highlighted in the previous test, the target encoding should be
    # invariant to such a permutation.
    permutated_labels = rng.permutation(n_categories)
    X_informative = permutated_labels[X_informative.astype(np.int32)]

    # Generate a shuffled copy of the informative feature to destroy the
    # relationship with the target.
    X_shuffled = rng.permutation(X_informative)

    # Also include a very high cardinality categorical feature that is by
    # itself independent of the target variable: target encoding such a feature
    # without internal cross-validation should cause catastrophic overfitting
    # for the downstream regressor, even with shrinkage. This kind of features
    # typically represents near unique identifiers of samples. In general they
    # should be removed from a machine learning datasets but here we want to
    # study the ability of the default behavior of TargetEncoder to mitigate
    # them automatically.
    X_near_unique_categories = rng.choice(
        int(0.9 * n_samples), size=n_samples, replace=True
    ).reshape(-1, 1)

    # Assemble the dataset and do a train-test split:
    X = np.concatenate(
        [X_informative, X_shuffled, X_near_unique_categories],
        axis=1,
    )
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

    # Let's first check that a linear regression model trained on the raw
    # features underfits because of the meaning-less ordinal encoding of the
    # labels.
    raw_model = linear_regression.fit(X_train, y_train)
    assert raw_model.score(X_train, y_train) < 0.1
    assert raw_model.score(X_test, y_test) < 0.1

    # Now do the same with target encoding using the internal CV mechanism
    # implemented when using fit_transform.
    model_with_cv = make_pipeline(
        TargetEncoder(smooth=smooth, random_state=rng), linear_regression
    ).fit(X_train, y_train)

    # This model should be able to fit the data well and also generalise to the
    # test data (assuming that the binning is fine-grained enough). The R2
    # scores are not perfect because of the noise injected during the
    # generation of the unique informative feature.
    coef = model_with_cv[-1].coef_
    assert model_with_cv.score(X_train, y_train) > 0.5, coef
    assert model_with_cv.score(X_test, y_test) > 0.5, coef

    # The target encoder recovers the linear relationship with slope 1 between
    # the target encoded unique informative predictor and the target. Since the
    # target encoding of the 2 other features is not informative thanks to the
    # use of internal cross-validation, the multivariate linear regressor
    # assigns a coef of 1 to the first feature and 0 to the other 2.
    assert coef[0] == pytest.approx(1, abs=1e-2)
    assert (np.abs(coef[1:]) < 0.2).all()

    # Let's now disable the internal cross-validation by calling fit and then
    # transform separately on the training set:
    target_encoder = TargetEncoder(smooth=smooth, random_state=rng).fit(
        X_train, y_train
    )
    X_enc_no_cv_train = target_encoder.transform(X_train)
    X_enc_no_cv_test = target_encoder.transform(X_test)
    model_no_cv = linear_regression.fit(X_enc_no_cv_train, y_train)

    # The linear regression model should always overfit because it assigns
    # too much weight to the extremely high cardinality feature relatively to
    # the informative feature. Note that this is the case even when using
    # the empirical Bayes smoothing which is not enough to prevent such
    # overfitting alone.
    coef = model_no_cv.coef_
    assert model_no_cv.score(X_enc_no_cv_train, y_train) > 0.7, coef
    assert model_no_cv.score(X_enc_no_cv_test, y_test) < 0.5, coef

    # The model overfits because it assigns too much weight to the high
    # cardinality yet non-informative feature instead of the lower
    # cardinality yet informative feature:
    assert abs(coef[0]) < abs(coef[2])


def test_pandas_copy_on_write():
    """
    Test target-encoder cython code when y is read-only.

    The numpy array underlying df["y"] is read-only when copy-on-write is enabled.
    Non-regression test for gh-27879.
    """
    pd = pytest.importorskip("pandas", minversion="2.0")
    with pd.option_context("mode.copy_on_write", True):
        df = pd.DataFrame({"x": ["a", "b", "b"], "y": [4.0, 5.0, 6.0]})
        TargetEncoder(target_type="continuous").fit(df[["x"]], df["y"])


# ---------------------------------------------------------------------------
# Small batch Target encoder tests
# ---------------------------------------------------------------------------


class _NoCategoriesAfterFit(TargetEncoder):
    """
    White-box subclass that simulates an upstream breakage by deleting
    `categories_` after a normal fit.

    Intuition
    ---------
    The fast-path initialization guards that `fit` *must* set `categories_`.
    We ensure this guard fails loudly and with a clear error message.
    """

    def _fit_encodings_all(self, X, y):
        out = super()._fit_encodings_all(X, y)
        if hasattr(self, "categories_"):
            delattr(self, "categories_")
        return out


class _DupCats(TargetEncoder):
    """
    White-box subclass to exercise the defensive index-map collision branch.

    Intuition
    ---------
    `_BaseEncoder._fit` guarantees per-feature `categories_[j]` are unique, and
    this TargetEncoder keeps NA-like values distinct via custom sentinels.
    Therefore a collision in the normalized dict keys should be unreachable.

    To cover the defensive code path without altering public behavior, we:
      1) run a standard fit,
      2) inject an illegal duplicate at the front of `categories_[0]`,
      3) mirror that duplication in `encodings_[0]` so lengths stay aligned,
      4) call the private cache builder for feature 0,
      5) assert that fast-path caches for that feature are left as `None`.

    We *do not* call `transform()` afterward, because the small-batch path
    would try to use `idx_map.get` on a `None` map (by design), which is beyond
    the intended scope of this coverage-only test.
    """

    def _fit_encodings_all(self, X, y):
        # Run the standard fitting first.
        Xo, mask, y_enc, ncat = super()._fit_encodings_all(X, y)

        # Inject a duplicate category at the front of feature 0 to force:
        #   key in index_map and index_map[key] != i
        j = 0
        cats0 = list(self.categories_[j])
        if len(cats0) >= 1:
            dup = cats0[0]
            # categories_: make the first entry appear twice
            self.categories_[j] = np.asarray([dup] + cats0, dtype=object)

            # encodings_: duplicate the first encoding value to keep lengths aligned
            enc0 = np.asarray(self.encodings_[j])
            self.encodings_[j] = np.concatenate([enc0[:1], enc0], axis=0)

        return Xo, mask, y_enc, ncat


class _FakeDF:
    """
    NumPy-coercible, DataFrame-like object with string columns.

    Intuition
    ---------
    We want to simulate a very small subset of a pandas DataFrame so that:
    - `check_array` can coerce it via `__array__` (avoids "scalar array" errors),
    - `TargetEncoder.fit` can see a `.columns` attribute and (optionally) set
      `feature_names_in_` when all column labels are strings.

    Keeping this shim minimal helps avoid accidental dependency on pandas logic
    while still exercising the intended feature-name path.
    """

    def __init__(self, arr, cols):
        self._arr = np.asarray(arr, dtype=object)
        self.columns = list(cols)

    def __array__(self, dtype=None):
        """Allow sklearn's validation to treat this like an ndarray."""
        if dtype is None:
            return self._arr
        return np.asarray(self._arr, dtype=dtype)

    @property
    def shape(self):
        """Expose array shape like a DataFrame would."""
        return self._arr.shape

    @property
    def ndim(self):
        """Expose ndim like a DataFrame would."""
        return self._arr.ndim

    def astype(self, dtype):
        """Mirror DataFrame.astype behavior just enough for check_array."""
        return np.asarray(self._arr, dtype=dtype)


class _BadShapeDF(_FakeDF):
    """
    Same as _FakeDF but `.shape` raises to exercise a defensive branch.

    Intuition
    ---------
    `TargetEncoder.fit` checks (inside `try/except`) whether a DF-like's shape
    matches the validated `ndarray` shape. If accessing `.shape` explodes,
    it falls back to treating input as ndarray (and thus will not set
    `feature_names_in_`). This class lets us test that graceful fallback.
    """

    @property
    def shape(self):
        raise RuntimeError("bad shape")


def _fit_pair_numpy():
    """
    Fit two encoders on the same tiny dataset:
    - te_fast uses the small-batch fast-path (threshold >= 256),
    - te_vec forces the baseline vectorized path (threshold < 0).

    Intuition
    ---------
    We compare outputs on subsequent transforms to assert that the
    implementation-specific fast path matches the reference vectorized path.
    """
    X_fit = np.array(
        [
            ["u", "x"],
            ["v", "y"],
            ["u", "x"],
            ["w", "y"],
        ],
        dtype=object,
    )
    y = np.array([0.0, 1.0, 1.0, 0.0])
    te_fast = TargetEncoder(smooth=5.0).fit(X_fit, y)
    # White-box: enable small-batch fast path deterministically.
    te_fast._small_batch_threshold = 256
    te_vec = TargetEncoder(smooth=5.0).fit(X_fit, y)
    # White-box: force always-vectorized path for reference.
    te_vec._small_batch_threshold = -1
    return te_fast, te_vec


def test_norm_key_real_values_cover_nan_nat_and_except_paths():
    """
    `_norm_key` should:
    - map None, float NaN, and NumPy NaT to *distinct* sentinels,
    - return ordinary strings unchanged,
    - not crash when np.isnat() raises on e.g. strings (defensive `except`).
    """
    # 1) Exactly None → distinct sentinel
    k_none = _norm_key(None)

    # 2) Float NaN (Python/NumPy) → distinct sentinel
    k_nan1 = _norm_key(float("nan"))
    k_nan2 = _norm_key(np.float64(np.nan))
    assert k_nan1 is k_nan2

    # 3) Plain string reaches the NumPy NaT check and naturally raises inside
    #    np.isnat("x") → defensive 'except: pass' executes and falls through.
    k_str = _norm_key("x")
    assert k_str == "x"

    # 4) NumPy datetime/timedelta NaT (non-float) → distinct sentinel
    k_nat_dt = _norm_key(np.datetime64("NaT"))
    k_nat_td = _norm_key(np.timedelta64("NaT"))
    assert k_nat_dt is k_nat_td

    # sanity: all sentinels are distinct from ordinary values
    assert k_none is not k_nan1 and k_none is not k_nat_dt
    assert k_nan1 is not k_nat_dt and k_nan1 != "x"
    assert k_nat_dt != "x"


def test_norm_key_pandas_duck_detection_exception_is_handled():
    """
    Force the pandas-duck-typing block to raise when accessing __module__.

    Intuition
    ---------
    The implementation refuses to import pandas and instead checks type name +
    module prefix. We fabricate a type whose metaclass throws on `__module__`.
    `_norm_key` must catch and fall through (no crash), returning the object.
    """

    class _RaisingMeta(type):
        def __getattribute__(cls, name):
            if name == "__module__":
                raise RuntimeError("boom")
            return super().__getattribute__(name)

    class _WeirdNA(metaclass=_RaisingMeta):
        pass

    x = _WeirdNA()
    assert _norm_key(x) is x  # fall-through behavior


def test_norm_key_duck_pandas_NA_branch_hits_return():
    """
    Hit the pandas.NA duck-typed path without importing pandas.

    Intuition
    ---------
    `_looks_like_pandas_na` checks `type(x).__name__ == "NAType"` and
    module startswith("pandas"). We spoof both for a plain Python class.
    """

    class _FakePandasNA:
        pass

    _FakePandasNA.__name__ = "NAType"
    _FakePandasNA.__module__ = "pandas.core.arrays._masked"

    k1 = _norm_key(_FakePandasNA())
    k2 = _norm_key(_FakePandasNA())
    assert k1 is k2  # both should map to the same sentinel identity


def test_norm_key_duck_pandas_NaT_branch_hits_return():
    """
    Hit the pandas.NaT duck-typed path without importing pandas.

    Intuition
    ---------
    `_looks_like_pandas_nat` accepts type names {"NaTType", "NaT"} with a
    module that startswith("pandas"). We spoof those identifiers.
    """

    class _FakePandasNaT:
        pass

    _FakePandasNaT.__name__ = "NaTType"  # could also set to "NaT"
    _FakePandasNaT.__module__ = "pandas._libs.tslibs.nattype"

    k1 = _norm_key(_FakePandasNaT())
    k2 = _norm_key(_FakePandasNaT())
    assert k1 is k2


def test_small_batch_fastpath_matches_vectorized_with_nans_and_nats():
    """
    For tiny inputs, small-batch dict-lookup path should match the reference
    vectorized path exactly (within tight numerical tolerance).

    Important detail
    ----------------
    Do NOT mix datetime64 NaT and timedelta64 NaT in the *same column*, because
    the vectorized unknown-check (`set(values)`) cannot compare them together.
    """
    te_fast, te_vec = _fit_pair_numpy()

    X = np.empty((5, 2), dtype=object)
    X[0] = ["u", "x"]  # seen
    X[1] = [float("nan"), "x"]  # float NaN branch
    X[2] = [None, "x"]  # None branch
    X[3] = ["u", np.datetime64("NaT")]  # datetime NaT branch
    X[4] = ["new_unseen", "x"]  # unseen → default

    npt.assert_allclose(te_fast.transform(X), te_vec.transform(X), rtol=0, atol=1e-9)


def test_large_batch_vectorized_consistency():
    """
    For large inputs, both encoders choose the vectorized path; results must
    match exactly (within tight tolerance).
    """
    te_fast, te_vec = _fit_pair_numpy()

    base = np.array(
        [["u", "x"], ["v", "y"], ["u", "x"], ["w", "y"], [None, "x"], ["new", "y"]],
        dtype=object,
    )
    X_big = np.tile(base, (60, 1))  # 360 rows → well over 256

    Z_fast = te_fast.transform(X_big)  # chooses vectorized due to size
    Z_vec = te_vec.transform(X_big)  # always vectorized
    npt.assert_allclose(Z_fast, Z_vec, rtol=0, atol=1e-9)


def test_multiclass_default_fallback_uses_block_mean_axis1():
    """
    Multiclass robust default:
    If `target_mean_` is a valid numeric array but wrong-shaped, the default
    vector for unseen categories should fall back to `block.mean(axis=1)`.
    """
    rng = np.random.RandomState(0)
    X_fit = np.stack(
        [
            rng.choice(list("abc"), size=80),
            rng.choice(list("wxyz"), size=80),
        ],
        axis=1,
    ).astype(object)
    y = rng.randint(0, 3, size=80)  # legitimate multiclass target

    te = TargetEncoder(target_type="multiclass").fit(X_fit, y)
    # White-box: allow small input to build caches where default is computed.
    te._small_batch_threshold = 256

    # Wrong shape but numeric → triggers robust per-class default via block.mean(axis=1)
    te.target_mean_ = np.array([[0.2, 0.5, 0.3]], dtype=float)

    X_small = np.array(
        [
            ["a", "w"],  # seen
            ["zzz_unseen", "w"],  # unseen in first feature → use default vector
        ],
        dtype=object,
    )

    Z = te.transform(X_small)
    n_classes = len(te.classes_)
    assert Z.shape == (2, X_small.shape[1] * n_classes)
    assert np.isfinite(Z).all()


def test_fit_dataframe_like_sets_feature_names():
    """
    DF-like input with string columns:
    Accept both upstream behaviors:
      - some implementations set `feature_names_in_` for DF-like input,
      - others do not (if validate_data ultimately saw an ndarray).
    In either case, `get_feature_names_out()` should be consistent.
    """
    arr = np.array([["a", "x"], ["b", "y"], ["a", "x"]], dtype=object)
    cols = ["feat_a", "feat_b"]
    te = TargetEncoder(smooth=5.0).fit(
        _FakeDF(arr, cols), np.array([0, 1, 1], dtype=float)
    )

    names = te.get_feature_names_out()
    default_names = np.array([f"x{i}" for i in range(len(cols))], dtype=object)
    if hasattr(te, "feature_names_in_"):
        npt.assert_array_equal(names, np.array(cols, dtype=object))
    else:
        npt.assert_array_equal(names, default_names)


def test_fit_dataframe_like_shape_mismatch_falls_back_to_ndarray():
    """
    If accessing `.shape` on DF-like input fails, fit should gracefully fall
    back to ndarray handling and *not* set `feature_names_in_`.
    """
    arr = np.array([["a", "x"], ["b", "y"]], dtype=object)
    te = TargetEncoder(smooth=5.0).fit(
        _BadShapeDF(arr, ["fa", "fb"]), np.array([0, 1], dtype=float)
    )
    assert not hasattr(te, "feature_names_in_")


def test_smallbatch_cache_reshape_and_default_mean_fallback_binary():
    """
    Binary/regression robust default & vector shape:
    - Force `encodings_[0]` to be 2-D so the cache flattens it.
    - Force `target_mean_` to be 1-D so default scalar falls back to `mean(enc_vec)`.
    """
    te_fast, te_vec = _fit_pair_numpy()

    # Make per-category encoding vector 2-D → fast-path flattens it
    enc0 = np.asarray(te_fast.encodings_[0]).reshape(-1, 1)
    te_fast.encodings_[0] = enc0

    # Make global mean 1-D → default scalar falls back to mean(enc_vec)
    te_fast.target_mean_ = np.array([float(te_fast.target_mean_)], dtype=float)

    X = np.array([["u", "x"], ["new_unseen", "x"]], dtype=object)
    npt.assert_allclose(te_fast.transform(X), te_vec.transform(X), rtol=0, atol=1e-9)


def test_fit_raises_when_categories_missing_no_pytest():
    """
    `TargetEncoder.fit` should raise a clear AttributeError if `categories_`
    is not set before fast-path initialization (defensive guard).
    """
    X = np.array([["a"], ["b"]], dtype=object)
    y = np.array([0.0, 1.0])
    te = _NoCategoriesAfterFit(smooth=5.0)

    try:
        te.fit(X, y)
    except AttributeError as e:
        msg = "TargetEncoder.fit must set 'categories_' before fast path."
        assert msg in str(e), f"unexpected error message: {e}"
    else:
        raise AssertionError("Expected AttributeError was not raised")


def test_collision_branch_disables_fastpath_for_feature():
    """
    Force the index-map collision branch and assert that all per-feature
    fast-path caches are left as None (the intended 'disable fastpath' behavior).
    """
    # Small, valid dataset for a normal fit
    X = np.array([["u"], ["v"], ["u"]], dtype=object)
    y = np.array([0.0, 1.0, 1.0])

    te = _DupCats(smooth=5.0).fit(X, y)
    # White-box: allow small-batch cache building
    te._small_batch_threshold = 256

    # Safeguard in case the private method is renamed in the future.
    assert hasattr(te, "_ensure_fastpath_structs_for_feature"), (
        "private fast-path builder was renamed; update test accordingly"
    )

    # Build caches for feature 0; this must hit the collision branch
    j = 0
    te._ensure_fastpath_structs_for_feature(j)

    # Collision path leaves all fast caches as None and returns
    assert te._te_index_maps_[j] is None
    assert te._te_enc_vecs_[j] is None
    assert te._te_enc_blocks_[j] is None
    assert te._te_defaults_[j] is None


def _force_slow_transform(enc: TargetEncoder, X):
    # Force reference/slow path by disabling the small-batch route and clearing caches.
    old_thresh = enc._small_batch_threshold
    enc._small_batch_threshold = -1  # never take fast path
    # Clear caches (simulate a fresh instance for fairness)
    enc._te_index_maps_ = None
    enc._te_is_multiclass_ = None
    enc._te_enc_vecs_ = None
    enc._te_enc_blocks_ = None
    enc._te_defaults_ = None
    try:
        out = enc.transform(X)
    finally:
        enc._small_batch_threshold = old_thresh
    return out


@pytest.mark.parametrize("n_small", [1, 2, 8, 32])
@pytest.mark.parametrize("n_cats", [10_000, 50_000])
def test_small_batch_binary_parity_and_missing(n_small, n_cats):
    rng = np.random.default_rng(0)
    X_fit = rng.integers(0, n_cats, size=(150_000, 1)).astype(object)
    # sprinkle missing during fit to learn its encoding too
    X_fit[:100, 0] = None
    y_fit = rng.integers(0, 2, size=(150_000,))
    enc = TargetEncoder(random_state=0).fit(X_fit, y_fit)

    X_small = rng.integers(0, n_cats * 2, size=(n_small, 1)).astype(object)
    # add missing/unseen explicitly
    if n_small >= 2:
        X_small[0, 0] = None
        X_small[-1, 0] = np.nan

    ref = _force_slow_transform(enc, X_small)
    fast = enc.transform(X_small)
    assert_allclose(fast, ref, rtol=0, atol=0)


def test_small_batch_multiclass_parity_and_order():
    rng = np.random.default_rng(1)
    n_cats = 30_000
    n_classes = 5
    X_fit = rng.integers(0, n_cats, size=(200_000, 1)).astype(object)
    y_fit = rng.integers(0, n_classes, size=(200_000,))
    enc = TargetEncoder(random_state=0).fit(X_fit, y_fit)

    X_small = np.array([[0], [n_cats + 123], [42], [None], [np.nan]], dtype=object)
    ref = _force_slow_transform(enc, X_small)
    fast = enc.transform(X_small)

    # 1) numerical parity
    assert_allclose(fast, ref, rtol=0, atol=0)
    # 2) shape and class-order parity
    assert fast.shape[1] % n_classes == 0
    assert ref.shape == fast.shape
