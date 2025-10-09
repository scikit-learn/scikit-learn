import warnings

import numpy as np
import pytest
from scipy import sparse

from sklearn.exceptions import NotFittedError
from sklearn.preprocessing import OneHotEncoder, OrdinalEncoder
from sklearn.utils._missing import is_scalar_nan
from sklearn.utils._testing import (
    _convert_container,
    assert_allclose,
    assert_array_equal,
)


def test_one_hot_encoder_sparse_dense():
    # check that sparse and dense will give the same results

    X = np.array([[3, 2, 1], [0, 1, 1]])
    enc_sparse = OneHotEncoder()
    enc_dense = OneHotEncoder(sparse_output=False)

    X_trans_sparse = enc_sparse.fit_transform(X)
    X_trans_dense = enc_dense.fit_transform(X)

    assert X_trans_sparse.shape == (2, 5)
    assert X_trans_dense.shape == (2, 5)

    assert sparse.issparse(X_trans_sparse)
    assert not sparse.issparse(X_trans_dense)

    # check outcome
    assert_array_equal(
        X_trans_sparse.toarray(), [[0.0, 1.0, 0.0, 1.0, 1.0], [1.0, 0.0, 1.0, 0.0, 1.0]]
    )
    assert_array_equal(X_trans_sparse.toarray(), X_trans_dense)


@pytest.mark.parametrize("handle_unknown", ["ignore", "infrequent_if_exist", "warn"])
def test_one_hot_encoder_handle_unknown(handle_unknown):
    X = np.array([[0, 2, 1], [1, 0, 3], [1, 0, 2]])
    X2 = np.array([[4, 1, 1]])

    # Test that one hot encoder raises error for unknown features
    # present during transform.
    oh = OneHotEncoder(handle_unknown="error")
    oh.fit(X)
    with pytest.raises(ValueError, match="Found unknown categories"):
        oh.transform(X2)

    # Test the ignore option, ignores unknown features (giving all 0's)
    oh = OneHotEncoder(handle_unknown=handle_unknown)
    oh.fit(X)
    X2_passed = X2.copy()
    assert_array_equal(
        oh.transform(X2_passed).toarray(),
        np.array([[0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]]),
    )
    # ensure transformed data was not modified in place
    assert_allclose(X2, X2_passed)


@pytest.mark.parametrize("handle_unknown", ["ignore", "infrequent_if_exist", "warn"])
def test_one_hot_encoder_handle_unknown_strings(handle_unknown):
    X = np.array(["11111111", "22", "333", "4444"]).reshape((-1, 1))
    X2 = np.array(["55555", "22"]).reshape((-1, 1))
    # Non Regression test for the issue #12470
    # Test the ignore option, when categories are numpy string dtype
    # particularly when the known category strings are larger
    # than the unknown category strings
    oh = OneHotEncoder(handle_unknown=handle_unknown)
    oh.fit(X)
    X2_passed = X2.copy()
    assert_array_equal(
        oh.transform(X2_passed).toarray(),
        np.array([[0.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]]),
    )
    # ensure transformed data was not modified in place
    assert_array_equal(X2, X2_passed)


@pytest.mark.parametrize("output_dtype", [np.int32, np.float32, np.float64])
@pytest.mark.parametrize("input_dtype", [np.int32, np.float32, np.float64])
def test_one_hot_encoder_dtype(input_dtype, output_dtype):
    X = np.asarray([[0, 1]], dtype=input_dtype).T
    X_expected = np.asarray([[1, 0], [0, 1]], dtype=output_dtype)

    oh = OneHotEncoder(categories="auto", dtype=output_dtype)
    assert_array_equal(oh.fit_transform(X).toarray(), X_expected)
    assert_array_equal(oh.fit(X).transform(X).toarray(), X_expected)

    oh = OneHotEncoder(categories="auto", dtype=output_dtype, sparse_output=False)
    assert_array_equal(oh.fit_transform(X), X_expected)
    assert_array_equal(oh.fit(X).transform(X), X_expected)


@pytest.mark.parametrize("output_dtype", [np.int32, np.float32, np.float64])
def test_one_hot_encoder_dtype_pandas(output_dtype):
    pd = pytest.importorskip("pandas")

    X_df = pd.DataFrame({"A": ["a", "b"], "B": [1, 2]})
    X_expected = np.array([[1, 0, 1, 0], [0, 1, 0, 1]], dtype=output_dtype)

    oh = OneHotEncoder(dtype=output_dtype)
    assert_array_equal(oh.fit_transform(X_df).toarray(), X_expected)
    assert_array_equal(oh.fit(X_df).transform(X_df).toarray(), X_expected)

    oh = OneHotEncoder(dtype=output_dtype, sparse_output=False)
    assert_array_equal(oh.fit_transform(X_df), X_expected)
    assert_array_equal(oh.fit(X_df).transform(X_df), X_expected)


def test_one_hot_encoder_feature_names():
    enc = OneHotEncoder()
    X = [
        ["Male", 1, "girl", 2, 3],
        ["Female", 41, "girl", 1, 10],
        ["Male", 51, "boy", 12, 3],
        ["Male", 91, "girl", 21, 30],
    ]

    enc.fit(X)
    feature_names = enc.get_feature_names_out()

    assert_array_equal(
        [
            "x0_Female",
            "x0_Male",
            "x1_1",
            "x1_41",
            "x1_51",
            "x1_91",
            "x2_boy",
            "x2_girl",
            "x3_1",
            "x3_2",
            "x3_12",
            "x3_21",
            "x4_3",
            "x4_10",
            "x4_30",
        ],
        feature_names,
    )

    feature_names2 = enc.get_feature_names_out(["one", "two", "three", "four", "five"])

    assert_array_equal(
        [
            "one_Female",
            "one_Male",
            "two_1",
            "two_41",
            "two_51",
            "two_91",
            "three_boy",
            "three_girl",
            "four_1",
            "four_2",
            "four_12",
            "four_21",
            "five_3",
            "five_10",
            "five_30",
        ],
        feature_names2,
    )

    with pytest.raises(ValueError, match="input_features should have length"):
        enc.get_feature_names_out(["one", "two"])


def test_one_hot_encoder_feature_names_unicode():
    enc = OneHotEncoder()
    X = np.array([["c❤t1", "dat2"]], dtype=object).T
    enc.fit(X)
    feature_names = enc.get_feature_names_out()
    assert_array_equal(["x0_c❤t1", "x0_dat2"], feature_names)
    feature_names = enc.get_feature_names_out(input_features=["n👍me"])
    assert_array_equal(["n👍me_c❤t1", "n👍me_dat2"], feature_names)


def test_one_hot_encoder_custom_feature_name_combiner():
    """Check the behaviour of `feature_name_combiner` as a callable."""

    def name_combiner(feature, category):
        return feature + "_" + repr(category)

    enc = OneHotEncoder(feature_name_combiner=name_combiner)
    X = np.array([["None", None]], dtype=object).T
    enc.fit(X)
    feature_names = enc.get_feature_names_out()
    assert_array_equal(["x0_'None'", "x0_None"], feature_names)
    feature_names = enc.get_feature_names_out(input_features=["a"])
    assert_array_equal(["a_'None'", "a_None"], feature_names)

    def wrong_combiner(feature, category):
        # we should be returning a Python string
        return 0

    enc = OneHotEncoder(feature_name_combiner=wrong_combiner).fit(X)
    err_msg = (
        "When `feature_name_combiner` is a callable, it should return a Python string."
    )
    with pytest.raises(TypeError, match=err_msg):
        enc.get_feature_names_out()


def test_one_hot_encoder_set_params():
    X = np.array([[1, 2]]).T
    oh = OneHotEncoder()
    # set params on not yet fitted object
    oh.set_params(categories=[[0, 1, 2, 3]])
    assert oh.get_params()["categories"] == [[0, 1, 2, 3]]
    assert oh.fit_transform(X).toarray().shape == (2, 4)
    # set params on already fitted object
    oh.set_params(categories=[[0, 1, 2, 3, 4]])
    assert oh.fit_transform(X).toarray().shape == (2, 5)


def check_categorical_onehot(X):
    enc = OneHotEncoder(categories="auto")
    Xtr1 = enc.fit_transform(X)

    enc = OneHotEncoder(categories="auto", sparse_output=False)
    Xtr2 = enc.fit_transform(X)

    assert_allclose(Xtr1.toarray(), Xtr2)

    assert sparse.issparse(Xtr1) and Xtr1.format == "csr"
    return Xtr1.toarray()


@pytest.mark.parametrize(
    "X",
    [
        [["def", 1, 55], ["abc", 2, 55]],
        np.array([[10, 1, 55], [5, 2, 55]]),
        np.array([["b", "A", "cat"], ["a", "B", "cat"]], dtype=object),
        np.array([["b", 1, "cat"], ["a", np.nan, "cat"]], dtype=object),
        np.array([["b", 1, "cat"], ["a", float("nan"), "cat"]], dtype=object),
        np.array([[None, 1, "cat"], ["a", 2, "cat"]], dtype=object),
        np.array([[None, 1, None], ["a", np.nan, None]], dtype=object),
        np.array([[None, 1, None], ["a", float("nan"), None]], dtype=object),
    ],
    ids=[
        "mixed",
        "numeric",
        "object",
        "mixed-nan",
        "mixed-float-nan",
        "mixed-None",
        "mixed-None-nan",
        "mixed-None-float-nan",
    ],
)
def test_one_hot_encoder(X):
    Xtr = check_categorical_onehot(np.array(X)[:, [0]])
    assert_allclose(Xtr, [[0, 1], [1, 0]])

    Xtr = check_categorical_onehot(np.array(X)[:, [0, 1]])
    assert_allclose(Xtr, [[0, 1, 1, 0], [1, 0, 0, 1]])

    Xtr = OneHotEncoder(categories="auto").fit_transform(X)
    assert_allclose(Xtr.toarray(), [[0, 1, 1, 0, 1], [1, 0, 0, 1, 1]])


@pytest.mark.parametrize("handle_unknown", ["ignore", "infrequent_if_exist", "warn"])
@pytest.mark.parametrize("sparse_", [False, True])
@pytest.mark.parametrize("drop", [None, "first"])
def test_one_hot_encoder_inverse(handle_unknown, sparse_, drop):
    X = [["abc", 2, 55], ["def", 1, 55], ["abc", 3, 55]]
    enc = OneHotEncoder(sparse_output=sparse_, drop=drop)
    X_tr = enc.fit_transform(X)
    exp = np.array(X, dtype=object)
    assert_array_equal(enc.inverse_transform(X_tr), exp)

    X = [[2, 55], [1, 55], [3, 55]]
    enc = OneHotEncoder(sparse_output=sparse_, categories="auto", drop=drop)
    X_tr = enc.fit_transform(X)
    exp = np.array(X)
    assert_array_equal(enc.inverse_transform(X_tr), exp)

    if drop is None:
        # with unknown categories
        # drop is incompatible with handle_unknown=ignore
        X = [["abc", 2, 55], ["def", 1, 55], ["abc", 3, 55]]
        enc = OneHotEncoder(
            sparse_output=sparse_,
            handle_unknown=handle_unknown,
            categories=[["abc", "def"], [1, 2], [54, 55, 56]],
        )
        X_tr = enc.fit_transform(X)
        exp = np.array(X, dtype=object)
        exp[2, 1] = None
        assert_array_equal(enc.inverse_transform(X_tr), exp)

        # with an otherwise numerical output, still object if unknown
        X = [[2, 55], [1, 55], [3, 55]]
        enc = OneHotEncoder(
            sparse_output=sparse_,
            categories=[[1, 2], [54, 56]],
            handle_unknown=handle_unknown,
        )
        X_tr = enc.fit_transform(X)
        exp = np.array(X, dtype=object)
        exp[2, 0] = None
        exp[:, 1] = None
        assert_array_equal(enc.inverse_transform(X_tr), exp)

    # incorrect shape raises
    X_tr = np.array([[0, 1, 1], [1, 0, 1]])
    msg = "Shape of the passed X data is not correct"
    with pytest.raises(ValueError, match=msg):
        enc.inverse_transform(X_tr)


@pytest.mark.parametrize("sparse_", [False, True])
@pytest.mark.parametrize(
    "X, X_trans",
    [
        ([[2, 55], [1, 55], [2, 55]], [[0, 1, 1], [0, 0, 0], [0, 1, 1]]),
        (
            [["one", "a"], ["two", "a"], ["three", "b"], ["two", "a"]],
            [[0, 0, 0, 0, 0], [0, 0, 0, 0, 1], [0, 1, 0, 0, 0]],
        ),
    ],
)
def test_one_hot_encoder_inverse_transform_raise_error_with_unknown(
    X, X_trans, sparse_
):
    """Check that `inverse_transform` raise an error with unknown samples, no
    dropped feature, and `handle_unknow="error`.
    Non-regression test for:
    https://github.com/scikit-learn/scikit-learn/issues/14934
    """
    enc = OneHotEncoder(sparse_output=sparse_).fit(X)
    msg = (
        r"Samples \\[(\\d )*\\d\\] can not be inverted when drop=None and "
        r"handle_unknown='error' because they contain all zeros"
    )

    if sparse_:
        # emulate sparse data transform by a one-hot encoder sparse.
        X_trans = _convert_container(X_trans, "sparse")
    with pytest.raises(ValueError, match=msg):
        enc.inverse_transform(X_trans)


def test_one_hot_encoder_inverse_if_binary():
    X = np.array([["Male", 1], ["Female", 3], ["Female", 2]], dtype=object)
    ohe = OneHotEncoder(drop="if_binary", sparse_output=False)
    X_tr = ohe.fit_transform(X)
    assert_array_equal(ohe.inverse_transform(X_tr), X)


@pytest.mark.parametrize("drop", ["if_binary", "first", None])
@pytest.mark.parametrize("reset_drop", ["if_binary", "first", None])
def test_one_hot_encoder_drop_reset(drop, reset_drop):
    # check that resetting drop option without refitting does not throw an error
    X = np.array([["Male", 1], ["Female", 3], ["Female", 2]], dtype=object)
    ohe = OneHotEncoder(drop=drop, sparse_output=False)
    ohe.fit(X)
    X_tr = ohe.transform(X)
    feature_names = ohe.get_feature_names_out()
    ohe.set_params(drop=reset_drop)
    assert_array_equal(ohe.inverse_transform(X_tr), X)
    assert_allclose(ohe.transform(X), X_tr)
    assert_array_equal(ohe.get_feature_names_out(), feature_names)


@pytest.mark.parametrize("method", ["fit", "fit_transform"])
@pytest.mark.parametrize("X", [[1, 2], np.array([3.0, 4.0])])
def test_X_is_not_1D(X, method):
    oh = OneHotEncoder()

    msg = "Expected 2D array, got 1D array instead"
    with pytest.raises(ValueError, match=msg):
        getattr(oh, method)(X)


@pytest.mark.parametrize("method", ["fit", "fit_transform"])
def test_X_is_not_1D_pandas(method):
    pd = pytest.importorskip("pandas")
    X = pd.Series([6, 3, 4, 6])
    oh = OneHotEncoder()

    msg = f"Expected a 2-dimensional container but got {type(X)} instead."
    with pytest.raises(ValueError, match=msg):
        getattr(oh, method)(X)


@pytest.mark.parametrize(
    "X, cat_exp, cat_dtype",
    [
        ([["abc", 55], ["def", 55]], [["abc", "def"], [55]], np.object_),
        (np.array([[1, 2], [3, 2]]), [[1, 3], [2]], np.integer),
        (
            np.array([["A", "cat"], ["B", "cat"]], dtype=object),
            [["A", "B"], ["cat"]],
            np.object_,
        ),
        (np.array([["A", "cat"], ["B", "cat"]]), [["A", "B"], ["cat"]], np.str_),
        (np.array([[1, 2], [np.nan, 2]]), [[1, np.nan], [2]], np.float64),
        (
            np.array([["A", np.nan], [None, np.nan]], dtype=object),
            [["A", None], [np.nan]],
            np.object_,
        ),
        (
            np.array([["A", float("nan")], [None, float("nan")]], dtype=object),
            [["A", None], [float("nan")]],
            np.object_,
        ),
    ],
    ids=[
        "mixed",
        "numeric",
        "object",
        "string",
        "missing-float",
        "missing-np.nan-object",
        "missing-float-nan-object",
    ],
)
def test_one_hot_encoder_categories(X, cat_exp, cat_dtype):
    # order of categories should not depend on order of samples
    for Xi in [X, X[::-1]]:
        enc = OneHotEncoder(categories="auto")
        enc.fit(Xi)
        # assert enc.categories == 'auto'
        assert isinstance(enc.categories_, list)
        for res, exp in zip(enc.categories_, cat_exp):
            res_list = res.tolist()
            if is_scalar_nan(exp[-1]):
                assert is_scalar_nan(res_list[-1])
                assert res_list[:-1] == exp[:-1]
            else:
                assert res.tolist() == exp
            assert np.issubdtype(res.dtype, cat_dtype)


@pytest.mark.parametrize("handle_unknown", ["ignore", "infrequent_if_exist", "warn"])
@pytest.mark.parametrize(
    "X, X2, cats, cat_dtype",
    [
        (
            np.array([["a", "b"]], dtype=object).T,
            np.array([["a", "d"]], dtype=object).T,
            [["a", "b", "c"]],
            np.object_,
        ),
        (
            np.array([[1, 2]], dtype="int64").T,
            np.array([[1, 4]], dtype="int64").T,
            [[1, 2, 3]],
            np.int64,
        ),
        (
            np.array([["a", "b"]], dtype=object).T,
            np.array([["a", "d"]], dtype=object).T,
            [np.array(["a", "b", "c"])],
            np.object_,
        ),
        (
            np.array([[None, "a"]], dtype=object).T,
            np.array([[None, "b"]], dtype=object).T,
            [[None, "a", "z"]],
            object,
        ),
        (
            np.array([["a", "b"]], dtype=object).T,
            np.array([["a", np.nan]], dtype=object).T,
            [["a", "b", "z"]],
            object,
        ),
        (
            np.array([["a", None]], dtype=object).T,
            np.array([["a", np.nan]], dtype=object).T,
            [["a", None, "z"]],
            object,
        ),
    ],
    ids=[
        "object",
        "numeric",
        "object-string",
        "object-string-none",
        "object-string-nan",
        "object-None-and-nan",
    ],
)
def test_one_hot_encoder_specified_categories(X, X2, cats, cat_dtype, handle_unknown):
    enc = OneHotEncoder(categories=cats)
    exp = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    assert_array_equal(enc.fit_transform(X).toarray(), exp)
    assert list(enc.categories[0]) == list(cats[0])
    assert enc.categories_[0].tolist() == list(cats[0])
    # manually specified categories should have same dtype as
    # the data when coerced from lists
    assert enc.categories_[0].dtype == cat_dtype

    # when specifying categories manually, unknown categories should already
    # raise when fitting
    enc = OneHotEncoder(categories=cats)
    with pytest.raises(ValueError, match="Found unknown categories"):
        enc.fit(X2)
    enc = OneHotEncoder(categories=cats, handle_unknown=handle_unknown)
    exp = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    assert_array_equal(enc.fit(X2).transform(X2).toarray(), exp)


def test_one_hot_encoder_unsorted_categories():
    X = np.array([["a", "b"]], dtype=object).T

    enc = OneHotEncoder(categories=[["b", "a", "c"]])
    exp = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0]])
    assert_array_equal(enc.fit(X).transform(X).toarray(), exp)
    assert_array_equal(enc.fit_transform(X).toarray(), exp)
    assert enc.categories_[0].tolist() == ["b", "a", "c"]
    assert np.issubdtype(enc.categories_[0].dtype, np.object_)

    # unsorted passed categories still raise for numerical values
    X = np.array([[1, 2]]).T
    enc = OneHotEncoder(categories=[[2, 1, 3]])
    msg = "Unsorted categories are not supported"
    with pytest.raises(ValueError, match=msg):
        enc.fit_transform(X)


@pytest.mark.parametrize("Encoder", [OneHotEncoder, OrdinalEncoder])
def test_encoder_nan_ending_specified_categories(Encoder):
    """Test encoder for specified categories that nan is at the end.

    Non-regression test for:
    https://github.com/scikit-learn/scikit-learn/issues/27088
    """
    cats = [np.array([0, np.nan, 1])]
    enc = Encoder(categories=cats)
    X = np.array([[0, 1]], dtype=object).T
    with pytest.raises(ValueError, match="Nan should be the last element"):
        enc.fit(X)


def test_one_hot_encoder_specified_categories_mixed_columns():
    # multiple columns
    X = np.array([["a", "b"], [0, 2]], dtype=object).T
    enc = OneHotEncoder(categories=[["a", "b", "c"], [0, 1, 2]])
    exp = np.array([[1.0, 0.0, 0.0, 1.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0, 0.0, 1.0]])
    assert_array_equal(enc.fit_transform(X).toarray(), exp)
    assert enc.categories_[0].tolist() == ["a", "b", "c"]
    assert np.issubdtype(enc.categories_[0].dtype, np.object_)
    assert enc.categories_[1].tolist() == [0, 1, 2]
    # integer categories but from object dtype data
    assert np.issubdtype(enc.categories_[1].dtype, np.object_)


def test_one_hot_encoder_pandas():
    pd = pytest.importorskip("pandas")

    X_df = pd.DataFrame({"A": ["a", "b"], "B": [1, 2]})

    Xtr = check_categorical_onehot(X_df)
    assert_allclose(Xtr, [[1, 0, 1, 0], [0, 1, 0, 1]])


@pytest.mark.parametrize(
    "drop, expected_names",
    [
        ("first", ["x0_c", "x2_b"]),
        ("if_binary", ["x0_c", "x1_2", "x2_b"]),
        (["c", 2, "b"], ["x0_b", "x2_a"]),
    ],
    ids=["first", "binary", "manual"],
)
def test_one_hot_encoder_feature_names_drop(drop, expected_names):
    X = [["c", 2, "a"], ["b", 2, "b"]]

    ohe = OneHotEncoder(drop=drop)
    ohe.fit(X)
    feature_names = ohe.get_feature_names_out()
    assert_array_equal(expected_names, feature_names)


def test_one_hot_encoder_drop_equals_if_binary():
    # Canonical case
    X = [[10, "yes"], [20, "no"], [30, "yes"]]
    expected = np.array(
        [[1.0, 0.0, 0.0, 1.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 1.0]]
    )
    expected_drop_idx = np.array([None, 0])

    ohe = OneHotEncoder(drop="if_binary", sparse_output=False)
    result = ohe.fit_transform(X)
    assert_array_equal(ohe.drop_idx_, expected_drop_idx)
    assert_allclose(result, expected)

    # with only one cat, the behaviour is equivalent to drop=None
    X = [["true", "a"], ["false", "a"], ["false", "a"]]
    expected = np.array([[1.0, 1.0], [0.0, 1.0], [0.0, 1.0]])
    expected_drop_idx = np.array([0, None])

    ohe = OneHotEncoder(drop="if_binary", sparse_output=False)
    result = ohe.fit_transform(X)
    assert_array_equal(ohe.drop_idx_, expected_drop_idx)
    assert_allclose(result, expected)


@pytest.mark.parametrize(
    "X",
    [
        [["abc", 2, 55], ["def", 1, 55]],
        np.array([[10, 2, 55], [20, 1, 55]]),
        np.array([["a", "B", "cat"], ["b", "A", "cat"]], dtype=object),
    ],
    ids=["mixed", "numeric", "object"],
)
def test_ordinal_encoder(X):
    enc = OrdinalEncoder()
    exp = np.array([[0, 1, 0], [1, 0, 0]], dtype="int64")
    assert_array_equal(enc.fit_transform(X), exp.astype("float64"))
    enc = OrdinalEncoder(dtype="int64")
    assert_array_equal(enc.fit_transform(X), exp)


@pytest.mark.parametrize(
    "X, X2, cats, cat_dtype",
    [
        (
            np.array([["a", "b"]], dtype=object).T,
            np.array([["a", "d"]], dtype=object).T,
            [["a", "b", "c"]],
            np.object_,
        ),
        (
            np.array([[1, 2]], dtype="int64").T,
            np.array([[1, 4]], dtype="int64").T,
            [[1, 2, 3]],
            np.int64,
        ),
        (
            np.array([["a", "b"]], dtype=object).T,
            np.array([["a", "d"]], dtype=object).T,
            [np.array(["a", "b", "c"])],
            np.object_,
        ),
    ],
    ids=["object", "numeric", "object-string-cat"],
)
def test_ordinal_encoder_specified_categories(X, X2, cats, cat_dtype):
    enc = OrdinalEncoder(categories=cats)
    exp = np.array([[0.0], [1.0]])
    assert_array_equal(enc.fit_transform(X), exp)
    assert list(enc.categories[0]) == list(cats[0])
    assert enc.categories_[0].tolist() == list(cats[0])
    # manually specified categories should have same dtype as
    # the data when coerced from lists
    assert enc.categories_[0].dtype == cat_dtype

    # when specifying categories manually, unknown categories should already
    # raise when fitting
    enc = OrdinalEncoder(categories=cats)
    with pytest.raises(ValueError, match="Found unknown categories"):
        enc.fit(X2)


def test_ordinal_encoder_inverse():
    X = [["abc", 2, 55], ["def", 1, 55]]
    enc = OrdinalEncoder()
    X_tr = enc.fit_transform(X)
    exp = np.array(X, dtype=object)
    assert_array_equal(enc.inverse_transform(X_tr), exp)

    # incorrect shape raises
    X_tr = np.array([[0, 1, 1, 2], [1, 0, 1, 0]])
    msg = "Shape of the passed X data is not correct"
    with pytest.raises(ValueError, match=msg):
        enc.inverse_transform(X_tr)


def test_ordinal_encoder_handle_unknowns_string():
    enc = OrdinalEncoder(handle_unknown="use_encoded_value", unknown_value=-2)
    X_fit = np.array([["a", "x"], ["b", "y"], ["c", "z"]], dtype=object)
    X_trans = np.array([["c", "xy"], ["bla", "y"], ["a", "x"]], dtype=object)
    enc.fit(X_fit)

    X_trans_enc = enc.transform(X_trans)
    exp = np.array([[2, -2], [-2, 1], [0, 0]], dtype="int64")
    assert_array_equal(X_trans_enc, exp)

    X_trans_inv = enc.inverse_transform(X_trans_enc)
    inv_exp = np.array([["c", None], [None, "y"], ["a", "x"]], dtype=object)
    assert_array_equal(X_trans_inv, inv_exp)


@pytest.mark.parametrize("dtype", [float, int])
def test_ordinal_encoder_handle_unknowns_numeric(dtype):
    enc = OrdinalEncoder(handle_unknown="use_encoded_value", unknown_value=-999)
    X_fit = np.array([[1, 7], [2, 8], [3, 9]], dtype=dtype)
    X_trans = np.array([[3, 12], [23, 8], [1, 7]], dtype=dtype)
    enc.fit(X_fit)

    X_trans_enc = enc.transform(X_trans)
    exp = np.array([[2, -999], [-999, 1], [0, 0]], dtype="int64")
    assert_array_equal(X_trans_enc, exp)

    X_trans_inv = enc.inverse_transform(X_trans_enc)
    inv_exp = np.array([[3, None], [None, 8], [1, 7]], dtype=object)
    assert_array_equal(X_trans_inv, inv_exp)


def test_ordinal_encoder_handle_unknowns_nan():
    # Make sure unknown_value=np.nan properly works

    enc = OrdinalEncoder(handle_unknown="use_encoded_value", unknown_value=np.nan)

    X_fit = np.array([[1], [2], [3]])
    enc.fit(X_fit)
    X_trans = enc.transform([[1], [2], [4]])
    assert_array_equal(X_trans, [[0], [1], [np.nan]])


def test_ordinal_encoder_handle_unknowns_nan_non_float_dtype():
    # Make sure an error is raised when unknown_value=np.nan and the dtype
    # isn't a float dtype
    enc = OrdinalEncoder(
        handle_unknown="use_encoded_value", unknown_value=np.nan, dtype=int
    )

    X_fit = np.array([[1], [2], [3]])
    with pytest.raises(ValueError, match="dtype parameter should be a float dtype"):
        enc.fit(X_fit)


def test_ordinal_encoder_raise_categories_shape():
    X = np.array([["Low", "Medium", "High", "Medium", "Low"]], dtype=object).T
    cats = ["Low", "Medium", "High"]
    enc = OrdinalEncoder(categories=cats)
    msg = "Shape mismatch: if categories is an array,"

    with pytest.raises(ValueError, match=msg):
        enc.fit(X)


def test_encoder_dtypes():
    # check that dtypes are preserved when determining categories
    enc = OneHotEncoder(categories="auto")
    exp = np.array([[1.0, 0.0, 1.0, 0.0], [0.0, 1.0, 0.0, 1.0]], dtype="float64")

    for X in [
        np.array([[1, 2], [3, 4]], dtype="int64"),
        np.array([[1, 2], [3, 4]], dtype="float64"),
        np.array([["a", "b"], ["c", "d"]]),  # str dtype
        np.array([[b"a", b"b"], [b"c", b"d"]]),  # bytes dtype
        np.array([[1, "a"], [3, "b"]], dtype="object"),
    ]:
        enc.fit(X)
        assert all([enc.categories_[i].dtype == X.dtype for i in range(2)])
        assert_array_equal(enc.transform(X).toarray(), exp)

    X = [[1, 2], [3, 4]]
    enc.fit(X)
    assert all([np.issubdtype(enc.categories_[i].dtype, np.integer) for i in range(2)])
    assert_array_equal(enc.transform(X).toarray(), exp)

    X = [[1, "a"], [3, "b"]]
    enc.fit(X)
    assert all([enc.categories_[i].dtype == "object" for i in range(2)])
    assert_array_equal(enc.transform(X).toarray(), exp)


def test_encoder_dtypes_pandas():
    # check dtype (similar to test_categorical_encoder_dtypes for dataframes)
    pd = pytest.importorskip("pandas")

    enc = OneHotEncoder(categories="auto")
    exp = np.array(
        [[1.0, 0.0, 1.0, 0.0, 1.0, 0.0], [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]],
        dtype="float64",
    )

    X = pd.DataFrame({"A": [1, 2], "B": [3, 4], "C": [5, 6]}, dtype="int64")
    enc.fit(X)
    assert all([enc.categories_[i].dtype == "int64" for i in range(2)])
    assert_array_equal(enc.transform(X).toarray(), exp)

    X = pd.DataFrame({"A": [1, 2], "B": ["a", "b"], "C": [3.0, 4.0]})
    expected_cat_type = ["int64", "object", "float64"]
    enc.fit(X)
    assert all([enc.categories_[i].dtype == expected_cat_type[i] for i in range(3)])
    assert_array_equal(enc.transform(X).toarray(), exp)


def test_one_hot_encoder_warning():
    enc = OneHotEncoder()
    X = [["Male", 1], ["Female", 3]]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        enc.fit_transform(X)


@pytest.mark.parametrize("drop", ["if_binary", "first"])
def test_ohe_handle_unknown_warn(drop):
    """Check handle_unknown='warn' works correctly."""

    X = [["a", 0], ["b", 2], ["b", 1]]

    ohe = OneHotEncoder(
        drop=drop,
        sparse_output=False,
        handle_unknown="warn",
        categories=[["b", "a"], [1, 2]],
    )
    ohe.fit(X)

    X_test = [["c", 1]]
    X_expected = np.array([[0, 0]])

    warn_msg = (
        r"Found unknown categories in columns \[0\] during transform. "
        r"These unknown categories will be encoded as all zeros"
    )
    with pytest.warns(UserWarning, match=warn_msg):
        X_trans = ohe.transform(X_test)
    assert_allclose(X_trans, X_expected)


@pytest.mark.parametrize("missing_value", [np.nan, None, float("nan")])
def test_one_hot_encoder_drop_manual(missing_value):
    cats_to_drop = ["def", 12, 3, 56, missing_value]
    enc = OneHotEncoder(drop=cats_to_drop)
    X = [
        ["abc", 12, 2, 55, "a"],
        ["def", 12, 1, 55, "a"],
        ["def", 12, 3, 56, missing_value],
    ]
    trans = enc.fit_transform(X).toarray()
    exp = [[1, 0, 1, 1, 1], [0, 1, 0, 1, 1], [0, 0, 0, 0, 0]]
    assert_array_equal(trans, exp)
    assert enc.drop is cats_to_drop

    dropped_cats = [
        cat[feature] for cat, feature in zip(enc.categories_, enc.drop_idx_)
    ]
    X_inv_trans = enc.inverse_transform(trans)
    X_array = np.array(X, dtype=object)

    # last value is np.nan
    if is_scalar_nan(cats_to_drop[-1]):
        assert_array_equal(dropped_cats[:-1], cats_to_drop[:-1])
        assert is_scalar_nan(dropped_cats[-1])
        assert is_scalar_nan(cats_to_drop[-1])
        # do not include the last column which includes missing values
        assert_array_equal(X_array[:, :-1], X_inv_trans[:, :-1])

        # check last column is the missing value
        assert_array_equal(X_array[-1, :-1], X_inv_trans[-1, :-1])
        assert is_scalar_nan(X_array[-1, -1])
        assert is_scalar_nan(X_inv_trans[-1, -1])
    else:
        assert_array_equal(dropped_cats, cats_to_drop)
        assert_array_equal(X_array, X_inv_trans)


@pytest.mark.parametrize("drop", [["abc", 3], ["abc", 3, 41, "a"]])
def test_invalid_drop_length(drop):
    enc = OneHotEncoder(drop=drop)
    err_msg = "`drop` should have length equal to the number"
    with pytest.raises(ValueError, match=err_msg):
        enc.fit([["abc", 2, 55], ["def", 1, 55], ["def", 3, 59]])


@pytest.mark.parametrize("density", [True, False], ids=["sparse", "dense"])
@pytest.mark.parametrize("drop", ["first", ["a", 2, "b"]], ids=["first", "manual"])
def test_categories(density, drop):
    ohe_base = OneHotEncoder(sparse_output=density)
    ohe_test = OneHotEncoder(sparse_output=density, drop=drop)
    X = [["c", 1, "a"], ["a", 2, "b"]]
    ohe_base.fit(X)
    ohe_test.fit(X)
    assert_array_equal(ohe_base.categories_, ohe_test.categories_)
    if drop == "first":
        assert_array_equal(ohe_test.drop_idx_, 0)
    else:
        for drop_cat, drop_idx, cat_list in zip(
            drop, ohe_test.drop_idx_, ohe_test.categories_
        ):
            assert cat_list[int(drop_idx)] == drop_cat
    assert isinstance(ohe_test.drop_idx_, np.ndarray)
    assert ohe_test.drop_idx_.dtype == object


@pytest.mark.parametrize("Encoder", [OneHotEncoder, OrdinalEncoder])
def test_encoders_has_categorical_tags(Encoder):
    assert Encoder().__sklearn_tags__().input_tags.categorical


@pytest.mark.parametrize(
    "kwargs",
    [
        {"max_categories": 2},
        {"min_frequency": 11},
        {"min_frequency": 0.29},
        {"max_categories": 2, "min_frequency": 6},
        {"max_categories": 4, "min_frequency": 12},
    ],
)
@pytest.mark.parametrize("categories", ["auto", [["a", "b", "c", "d"]]])
def test_ohe_infrequent_two_levels(kwargs, categories):
    """Test that different parameters for combine 'a', 'c', and 'd' into
    the infrequent category works as expected."""

    X_train = np.array([["a"] * 5 + ["b"] * 20 + ["c"] * 10 + ["d"] * 3]).T
    ohe = OneHotEncoder(
        categories=categories,
        handle_unknown="infrequent_if_exist",
        sparse_output=False,
        **kwargs,
    ).fit(X_train)
    assert_array_equal(ohe.infrequent_categories_, [["a", "c", "d"]])

    X_test = [["b"], ["a"], ["c"], ["d"], ["e"]]
    expected = np.array([[1, 0], [0, 1], [0, 1], [0, 1], [0, 1]])

    X_trans = ohe.transform(X_test)
    assert_allclose(expected, X_trans)

    expected_inv = [[col] for col in ["b"] + ["infrequent_sklearn"] * 4]
    X_inv = ohe.inverse_transform(X_trans)
    assert_array_equal(expected_inv, X_inv)

    feature_names = ohe.get_feature_names_out()
    assert_array_equal(["x0_b", "x0_infrequent_sklearn"], feature_names)


@pytest.mark.parametrize("drop", ["if_binary", "first", None])
def test_ohe_infrequent_two_levels_drop_frequent(drop):
    """Test two levels and dropping the frequent category."""

    X_train = np.array([["a"] * 5 + ["b"] * 20 + ["c"] * 10 + ["d"] * 3]).T
    ohe = OneHotEncoder(
        handle_unknown="infrequent_if_exist",
        sparse_output=False,
        max_categories=2,
        drop=drop,
    ).fit(X_train)
    assert ohe.categories_[0][ohe.drop_idx_[0]] == "b"

    X_test = np.array([["b"], ["c"]])
    X_trans = ohe.transform(X_test)
    assert_allclose([[0], [1]], X_trans)

    feature_names = ohe.get_feature_names_out()
    assert_array_equal(["x0_infrequent_sklearn"], feature_names)

    X_inverse = ohe.inverse_transform(X_trans)
    assert_array_equal([["b"], ["infrequent_sklearn"]], X_inverse)


@pytest.mark.parametrize("drop", [["a"], ["d"]])
def test_ohe_infrequent_two_levels_drop_infrequent_errors(drop):
    """Test two levels and dropping any infrequent category removes the
    whole infrequent category."""

    X_train = np.array([["a"] * 5 + ["b"] * 20 + ["c"] * 10 + ["d"] * 3]).T
    ohe = OneHotEncoder(
        handle_unknown="infrequent_if_exist",
        sparse_output=False,
        max_categories=2,
        drop=drop,
    )

    msg = f"Unable to drop category {drop[0]!r} from feature 0 because it is infrequent"
    with pytest.raises(ValueError, match=msg):
        ohe.fit(X_train)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"max_categories": 3},
        {"min_frequency": 6},
        {"min_frequency": 9},
        {"min_frequency": 0.24},
        {"min_frequency": 0.16},
        {"max_categories": 3, "min_frequency": 8},
        {"max_categories": 4, "min_frequency": 6},
    ],
)
def test_ohe_infrequent_three_levels(kwargs):
    """Test that different parameters for combing 'a', and 'd' into
    the infrequent category works as expected."""

    X_train = np.array([["a"] * 5 + ["b"] * 20 + ["c"] * 10 + ["d"] * 3]).T
    ohe = OneHotEncoder(
        handle_unknown="infrequent_if_exist", sparse_output=False, **kwargs
    ).fit(X_train)
    assert_array_equal(ohe.infrequent_categories_, [["a", "d"]])

    X_test = [["b"], ["a"], ["c"], ["d"], ["e"]]
    expected = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0], [0, 0, 1], [0, 0, 1]])

    X_trans = ohe.transform(X_test)
    assert_allclose(expected, X_trans)

    expected_inv = [
        ["b"],
        ["infrequent_sklearn"],
        ["c"],
        ["infrequent_sklearn"],
        ["infrequent_sklearn"],
    ]
    X_inv = ohe.inverse_transform(X_trans)
    assert_array_equal(expected_inv, X_inv)

    feature_names = ohe.get_feature_names_out()
    assert_array_equal(["x0_b", "x0_c", "x0_infrequent_sklearn"], feature_names)


@pytest.mark.parametrize("drop", ["first", ["b"]])
def test_ohe_infrequent_three_levels_drop_frequent(drop):
    """Test three levels and dropping the frequent category."""

    X_train = np.array([["a"] * 5 + ["b"] * 20 + ["c"] * 10 + ["d"] * 3]).T
    ohe = OneHotEncoder(
        handle_unknown="infrequent_if_exist",
        sparse_output=False,
        max_categories=3,
        drop=drop,
    ).fit(X_train)

    X_test = np.array([["b"], ["c"], ["d"]])
    assert_allclose([[0, 0], [1, 0], [0, 1]], ohe.transform(X_test))

    # Check handle_unknown="ignore"
    ohe.set_params(handle_unknown="ignore").fit(X_train)
    msg = "Found unknown categories"
    with pytest.warns(UserWarning, match=msg):
        X_trans = ohe.transform([["b"], ["e"]])

    assert_allclose([[0, 0], [0, 0]], X_trans)


@pytest.mark.parametrize("drop", [["a"], ["d"]])
def test_ohe_infrequent_three_levels_drop_infrequent_errors(drop):
    """Test three levels and dropping the infrequent category."""
    X_train = np.array([["a"] * 5 + ["b"] * 20 + ["c"] * 10 + ["d"] * 3]).T
    ohe = OneHotEncoder(
        handle_unknown="infrequent_if_exist",
        sparse_output=False,
        max_categories=3,
        drop=drop,
    )

    msg = f"Unable to drop category {drop[0]!r} from feature 0 because it is infrequent"
    with pytest.raises(ValueError, match=msg):
        ohe.fit(X_train)


def test_ohe_infrequent_handle_unknown_error():
    """Test that different parameters for combining 'a', and 'd' into
    the infrequent category works as expected."""

    X_train = np.array([["a"] * 5 + ["b"] * 20 + ["c"] * 10 + ["d"] * 3]).T
    ohe = OneHotEncoder(
        handle_unknown="error", sparse_output=False, max_categories=3
    ).fit(X_train)
    assert_array_equal(ohe.infrequent_categories_, [["a", "d"]])

    # all categories are known
    X_test = [["b"], ["a"], ["c"], ["d"]]
    expected = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0], [0, 0, 1]])

    X_trans = ohe.transform(X_test)
    assert_allclose(expected, X_trans)

    # 'bad' is not known and will error
    X_test = [["bad"]]
    msg = r"Found unknown categories \['bad'\] in column 0"
    with pytest.raises(ValueError, match=msg):
        ohe.transform(X_test)


@pytest.mark.parametrize(
    "kwargs", [{"max_categories": 3, "min_frequency": 1}, {"min_frequency": 4}]
)
def test_ohe_infrequent_two_levels_user_cats_one_frequent(kwargs):
    """'a' is the only frequent category, all other categories are infrequent."""

    X_train = np.array([["a"] * 5 + ["e"] * 30], dtype=object).T
    ohe = OneHotEncoder(
        categories=[["c", "d", "a", "b"]],
        sparse_output=False,
        handle_unknown="infrequent_if_exist",
        **kwargs,
    ).fit(X_train)

    X_test = [["a"], ["b"], ["c"], ["d"], ["e"]]
    expected = np.array([[1, 0], [0, 1], [0, 1], [0, 1], [0, 1]])

    X_trans = ohe.transform(X_test)
    assert_allclose(expected, X_trans)

    # 'a' is dropped
    drops = ["first", "if_binary", ["a"]]
    X_test = [["a"], ["c"]]
    for drop in drops:
        ohe.set_params(drop=drop).fit(X_train)
        assert_allclose([[0], [1]], ohe.transform(X_test))


def test_ohe_infrequent_two_levels_user_cats():
    """Test that the order of the categories provided by a user is respected."""
    X_train = np.array(
        [["a"] * 5 + ["b"] * 20 + ["c"] * 10 + ["d"] * 3], dtype=object
    ).T
    ohe = OneHotEncoder(
        categories=[["c", "d", "a", "b"]],
        sparse_output=False,
        handle_unknown="infrequent_if_exist",
        max_categories=2,
    ).fit(X_train)

    assert_array_equal(ohe.infrequent_categories_, [["c", "d", "a"]])

    X_test = [["b"], ["a"], ["c"], ["d"], ["e"]]
    expected = np.array([[1, 0], [0, 1], [0, 1], [0, 1], [0, 1]])

    X_trans = ohe.transform(X_test)
    assert_allclose(expected, X_trans)

    # 'infrequent' is used to denote the infrequent categories for
    # `inverse_transform`
    expected_inv = [[col] for col in ["b"] + ["infrequent_sklearn"] * 4]
    X_inv = ohe.inverse_transform(X_trans)
    assert_array_equal(expected_inv, X_inv)


def test_ohe_infrequent_three_levels_user_cats():
    """Test that the order of the categories provided by a user is respected.
    In this case 'c' is encoded as the first category and 'b' is encoded
    as the second one."""

    X_train = np.array(
        [["a"] * 5 + ["b"] * 20 + ["c"] * 10 + ["d"] * 3], dtype=object
    ).T
    ohe = OneHotEncoder(
        categories=[["c", "d", "b", "a"]],
        sparse_output=False,
        handle_unknown="infrequent_if_exist",
        max_categories=3,
    ).fit(X_train)

    assert_array_equal(ohe.infrequent_categories_, [["d", "a"]])

    X_test = [["b"], ["a"], ["c"], ["d"], ["e"]]
    expected = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0], [0, 0, 1], [0, 0, 1]])

    X_trans = ohe.transform(X_test)
    assert_allclose(expected, X_trans)

    # 'infrequent' is used to denote the infrequent categories for
    # `inverse_transform`
    expected_inv = [
        ["b"],
        ["infrequent_sklearn"],
        ["c"],
        ["infrequent_sklearn"],
        ["infrequent_sklearn"],
    ]
    X_inv = ohe.inverse_transform(X_trans)
    assert_array_equal(expected_inv, X_inv)


def test_ohe_infrequent_mixed():
    """Test infrequent categories where feature 0 has infrequent categories,
    and feature 1 does not."""

    # X[:, 0] 1 and 2 are infrequent
    # X[:, 1] nothing is infrequent
    X = np.c_[[0, 1, 3, 3, 3, 3, 2, 0, 3], [0, 0, 0, 0, 1, 1, 1, 1, 1]]

    ohe = OneHotEncoder(max_categories=3, drop="if_binary", sparse_output=False)
    ohe.fit(X)

    X_test = [[3, 0], [1, 1]]
    X_trans = ohe.transform(X_test)

    # feature 1 is binary so it drops a category 0
    assert_allclose(X_trans, [[0, 1, 0, 0], [0, 0, 1, 1]])


def test_ohe_infrequent_multiple_categories():
    """Test infrequent categories with feature matrix with 3 features."""

    X = np.c_[
        [0, 1, 3, 3, 3, 3, 2, 0, 3],
        [0, 0, 5, 1, 1, 10, 5, 5, 0],
        [1, 0, 1, 0, 1, 0, 1, 0, 1],
    ]

    ohe = OneHotEncoder(
        categories="auto", max_categories=3, handle_unknown="infrequent_if_exist"
    )
    # X[:, 0] 1 and 2 are infrequent
    # X[:, 1] 1 and 10 are infrequent
    # X[:, 2] nothing is infrequent

    X_trans = ohe.fit_transform(X).toarray()
    assert_array_equal(ohe.infrequent_categories_[0], [1, 2])
    assert_array_equal(ohe.infrequent_categories_[1], [1, 10])
    assert_array_equal(ohe.infrequent_categories_[2], None)

    # 'infrequent' is used to denote the infrequent categories
    # For the first column, 1 and 2 have the same frequency. In this case,
    # 1 will be chosen to be the feature name because is smaller lexiconically
    feature_names = ohe.get_feature_names_out()
    assert_array_equal(
        [
            "x0_0",
            "x0_3",
            "x0_infrequent_sklearn",
            "x1_0",
            "x1_5",
            "x1_infrequent_sklearn",
            "x2_0",
            "x2_1",
        ],
        feature_names,
    )

    expected = [
        [1, 0, 0, 1, 0, 0, 0, 1],
        [0, 0, 1, 1, 0, 0, 1, 0],
        [0, 1, 0, 0, 1, 0, 0, 1],
        [0, 1, 0, 0, 0, 1, 1, 0],
        [0, 1, 0, 0, 0, 1, 0, 1],
        [0, 1, 0, 0, 0, 1, 1, 0],
        [0, 0, 1, 0, 1, 0, 0, 1],
        [1, 0, 0, 0, 1, 0, 1, 0],
        [0, 1, 0, 1, 0, 0, 0, 1],
    ]

    assert_allclose(expected, X_trans)

    X_test = [[3, 1, 2], [4, 0, 3]]

    X_test_trans = ohe.transform(X_test)

    # X[:, 2] does not have an infrequent category, thus it is encoded as all
    # zeros
    expected = [[0, 1, 0, 0, 0, 1, 0, 0], [0, 0, 1, 1, 0, 0, 0, 0]]
    assert_allclose(expected, X_test_trans.toarray())

    X_inv = ohe.inverse_transform(X_test_trans)
    expected_inv = np.array(
        [[3, "infrequent_sklearn", None], ["infrequent_sklearn", 0, None]], dtype=object
    )
    assert_array_equal(expected_inv, X_inv)

    # error for unknown categories
    ohe = OneHotEncoder(
        categories="auto", max_categories=3, handle_unknown="error"
    ).fit(X)
    with pytest.raises(ValueError, match="Found unknown categories"):
        ohe.transform(X_test)

    # only infrequent or known categories
    X_test = [[1, 1, 1], [3, 10, 0]]
    X_test_trans = ohe.transform(X_test)

    expected = [[0, 0, 1, 0, 0, 1, 0, 1], [0, 1, 0, 0, 0, 1, 1, 0]]
    assert_allclose(expected, X_test_trans.toarray())

    X_inv = ohe.inverse_transform(X_test_trans)

    expected_inv = np.array(
        [["infrequent_sklearn", "infrequent_sklearn", 1], [3, "infrequent_sklearn", 0]],
        dtype=object,
    )
    assert_array_equal(expected_inv, X_inv)


def test_ohe_infrequent_multiple_categories_dtypes():
    """Test infrequent categories with a pandas dataframe with multiple dtypes."""

    pd = pytest.importorskip("pandas")
    categorical_dtype = pd.CategoricalDtype(["bird", "cat", "dog", "snake"])
    X = pd.DataFrame(
        {
            "str": ["a", "f", "c", "f", "f", "a", "c", "b", "b"],
            "int": [5, 3, 0, 10, 10, 12, 0, 3, 5],
            "categorical": pd.Series(
                ["dog"] * 4 + ["cat"] * 3 + ["snake"] + ["bird"],
                dtype=categorical_dtype,
            ),
        },
        columns=["str", "int", "categorical"],
    )

    ohe = OneHotEncoder(
        categories="auto", max_categories=3, handle_unknown="infrequent_if_exist"
    )
    # X[:, 0] 'a', 'b', 'c' have the same frequency. 'a' and 'b' will be
    # considered infrequent because they appear first when sorted

    # X[:, 1] 0, 3, 5, 10 has frequency 2 and 12 has frequency 1.
    # 0, 3, 12 will be considered infrequent because they appear first when
    # sorted.

    # X[:, 2] "snake" and "bird" or infrequent

    assert_array_equal(ohe.infrequent_categories_[0], ["a", "b"])
    assert_array_equal(ohe.infrequent_categories_[1], [0, 3, 12])
    assert_array_equal(ohe.infrequent_categories_[2], ["bird", "snake"])

    X_test = pd.DataFrame(
        {
            "str": ["a", "b", "f", "c"],
            "int": [12, 0, 10, 5],
            "categorical": pd.Series(
                ["cat"] + ["snake"] + ["bird"] + ["dog"],
                dtype=categorical_dtype,
            ),
        },
        columns=["str", "int", "categorical"],
    )
    expected_trans = [[2, 2, 0], [2, 2, 2], [1, 1, 2], [0, 0, 1]]

    X_trans = ohe.transform(X_test)
    assert_allclose(X_trans, expected_trans)


def test_ordinal_encoder_infrequent_three_levels(kwargs):
    """Test parameters for grouping 'a', and 'd' into the infrequent category."""

    X_train = np.array([["a"] * 5 + ["b"] * 20 + ["c"] * 10 + ["d"] * 3]).T
    ordinal = OrdinalEncoder(
        handle_unknown="use_encoded_value", unknown_value=-1, **kwargs
    ).fit(X_train)
    assert_array_equal(ordinal.categories_, [["a", "b", "c", "d"]])
    assert_array_equal(ordinal.infrequent_categories_, [["a", "d"]])

    X_test = [["a"], ["b"], ["c"], ["d"], ["z"]]
    expected_trans = [[2], [0], [1], [2], [-1]]

    X_trans = ordinal.transform(X_test)
    assert_allclose(X_trans, expected_trans)

    X_inverse = ordinal.inverse_transform(X_trans)
    expected_inverse = [
        ["infrequent_sklearn"],
        ["b"],
        ["c"],
        ["infrequent_sklearn"],
        [None],
    ]
    assert_array_equal(X_inverse, expected_inverse)


def test_ordinal_encoder_infrequent_three_levels_user_cats():
    """Test that the order of the categories provided by a user is respected.

    In this case 'c' is encoded as the first category and 'b' is encoded
    as the second one.
    """

    X_train = np.array(
        [["a"] * 5 + ["b"] * 20 + ["c"] * 10 + ["d"] * 3], dtype=object
    ).T
    ordinal = OrdinalEncoder(
        categories=[["c", "d", "b", "a"]],
        max_categories=3,
        handle_unknown="use_encoded_value",
        unknown_value=-1,
    ).fit(X_train)
    assert_array_equal(ordinal.categories_, [["c", "d", "b", "a"]])
    assert_array_equal(ordinal.infrequent_categories_, [["d", "a"]])

    X_test = [["a"], ["b"], ["c"], ["d"], ["z"]]
    expected_trans = [[2], [1], [0], [2], [-1]]

    X_trans = ordinal.transform(X_test)
    assert_allclose(X_trans, expected_trans)

    X_inverse = ordinal.inverse_transform(X_trans)
    expected_inverse = [
        ["infrequent_sklearn"],
        ["b"],
        ["c"],
        ["infrequent_sklearn"],
        [None],
    ]
    assert_array_equal(X_inverse, expected_inverse)


def test_ordinal_encoder_infrequent_mixed():
    """Test when feature 0 has infrequent categories and feature 1 does not."""

    X = np.column_stack(([0, 1, 3, 3, 3, 3, 2, 0, 3], [0, 0, 0, 0, 1, 1, 1, 1, 1]))

    ordinal = OrdinalEncoder(max_categories=3).fit(X)

    assert_array_equal(ordinal.infrequent_categories_[0], [1, 2])
    assert ordinal.infrequent_categories_[1] is None

    X_test = [[3, 0], [1, 1]]
    expected_trans = [[1, 0], [2, 1]]

    X_trans = ordinal.transform(X_test)
    assert_allclose(X_trans, expected_trans)

    X_inverse = ordinal.inverse_transform(X_trans)
    expected_inverse = np.array([[3, 0], ["infrequent_sklearn", 1]], dtype=object)
    assert_array_equal(X_inverse, expected_inverse)


def test_ordinal_encoder_infrequent_multiple_categories_dtypes():
    """Test infrequent categories with a pandas DataFrame with multiple dtypes."""

    pd = pytest.importorskip("pandas")
    categorical_dtype = pd.CategoricalDtype(["bird", "cat", "dog", "snake"])
    X = pd.DataFrame(
        {
            "str": ["a", "f", "c", "f", "f", "a", "c", "b", "b"],
            "int": [5, 3, 0, 10, 10, 12, 0, 3, 5],
            "categorical": pd.Series(
                ["dog"] * 4 + ["cat"] * 3 + ["snake"] + ["bird"],
                dtype=categorical_dtype,
            ),
        },
        columns=["str", "int", "categorical"],
    )

    ordinal = OrdinalEncoder(max_categories=3).fit(X)
    # X[:, 0] 'a', 'b', 'c' have the same frequency. 'a' and 'b' will be
    # considered infrequent because they appear first when sorted

    # X[:, 1] 0, 3, 5, 10 has frequency 2 and 12 has frequency 1.
    # 0, 3, 12 will be considered infrequent because they appear first when
    # sorted.

    # X[:, 2] "snake" and "bird" or infrequent

    assert_array_equal(ordinal.infrequent_categories_[0], ["a", "b"])
    assert_array_equal(ordinal.infrequent_categories_[1], [0, 3, 12])
    assert_array_equal(ordinal.infrequent_categories_[2], ["bird", "snake"])

    X_test = pd.DataFrame(
        {
            "str": ["a", "b", "f", "c"],
            "int": [12, 0, 10, 5],
            "categorical": pd.Series(
                ["cat"] + ["snake"] + ["bird"] + ["dog"],
                dtype=categorical_dtype,
            ),
        },
        columns=["str", "int", "categorical"],
    )
    expected_trans = [[2, 2, 0], [2, 2, 2], [1, 1, 2], [0, 0, 1]]

    X_trans = ordinal.transform(X_test)
    assert_allclose(X_trans, expected_trans)


def test_ordinal_encoder_infrequent_custom_mapping():
    """Check behavior of unknown_value and encoded_missing_value with infrequent."""
    X_train = np.array(
        [["a"] * 5 + ["b"] * 20 + ["c"] * 10 + ["d"] * 3 + [np.nan]], dtype=object
    ).T

    ordinal = OrdinalEncoder(
        handle_unknown="use_encoded_value",
        unknown_value=2,
        max_categories=2,
        encoded_missing_value=3,
    ).fit(X_train)
    assert_array_equal(ordinal.infrequent_categories_, [["a", "c", "d"]])

    X_test = np.array([["a"], ["b"], ["c"], ["d"], ["e"], [np.nan]], dtype=object)
    expected_trans = [[1], [0], [1], [1], [2], [3]]

    X_trans = ordinal.transform(X_test)
    assert_allclose(X_trans, expected_trans)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"max_categories": 6},
        {"min_frequency": 2},
    ],
)
def test_ordinal_encoder_all_frequent(kwargs):
    """All categories are considered frequent have same encoding as default encoder."""
    X_train = np.array(
        [["a"] * 5 + ["b"] * 20 + ["c"] * 10 + ["d"] * 3], dtype=object
    ).T

    adjusted_encoder = OrdinalEncoder(
        **kwargs, handle_unknown="use_encoded_value", unknown_value=-1
    ).fit(X_train)
    default_encoder = OrdinalEncoder(
        handle_unknown="use_encoded_value", unknown_value=-1
    ).fit(X_train)

    X_test = [["a"], ["b"], ["c"], ["d"], ["e"]]

    assert_allclose(
        adjusted_encoder.transform(X_test), default_encoder.transform(X_test)
    )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"max_categories": 1},
        {"min_frequency": 100},
    ],
)
def test_ordinal_encoder_all_infrequent(kwargs):
    """When all categories are infrequent, they are all encoded as zero."""
    X_train = np.array(
        [["a"] * 5 + ["b"] * 20 + ["c"] * 10 + ["d"] * 3], dtype=object
    ).T
    encoder = OrdinalEncoder(
        **kwargs, handle_unknown="use_encoded_value", unknown_value=-1
    ).fit(X_train)

    X_test = [["a"], ["b"], ["c"], ["d"], ["e"]]
    assert_allclose(encoder.transform(X_test), [[0], [0], [0], [0], [-1]])


def test_ordinal_encoder_missing_appears_frequent():
    """Check behavior when missing value appears frequently."""
    X = np.array(
        [[np.nan] * 20 + ["dog"] * 10 + ["cat"] * 5 + ["snake"] + ["deer"]],
        dtype=object,
    ).T
    ordinal = OrdinalEncoder(max_categories=3).fit(X)

    X_test = np.array([["snake", "cat", "dog", np.nan]], dtype=object).T
    X_trans = ordinal.transform(X_test)
    assert_allclose(X_trans, [[2], [0], [1], [np.nan]])


def test_ordinal_encoder_missing_appears_infrequent():
    """Check behavior when missing value appears infrequently."""

    # feature 0 has infrequent categories
    # feature 1 has no infrequent categories
    X = np.array(
        [
            [np.nan] + ["dog"] * 10 + ["cat"] * 5 + ["snake"] + ["deer"],
            ["red"] * 9 + ["green"] * 9,
        ],
        dtype=object,
    ).T
    ordinal = OrdinalEncoder(min_frequency=4).fit(X)

    X_test = np.array(
        [
            ["snake", "red"],
            ["deer", "green"],
            [np.nan, "green"],
            ["dog", "green"],
            ["cat", "red"],
        ],
        dtype=object,
    )
    X_trans = ordinal.transform(X_test)
    assert_allclose(X_trans, [[2, 1], [2, 0], [np.nan, 0], [1, 0], [0, 1]])


@pytest.mark.parametrize("Encoder", [OneHotEncoder, OrdinalEncoder])
def test_encoder_not_fitted(Encoder):
    """Check that we raise a `NotFittedError` by calling transform before fit with
    the encoders.

    One could expect that the passing the `categories` argument to the encoder
    would make it stateless. However, `fit` is making a couple of check, such as the
    position of `np.nan`.
    """
    X = np.array([["A"], ["B"], ["C"]], dtype=object)
    encoder = Encoder(categories=[["A", "B", "C"]])
    with pytest.raises(NotFittedError):
        encoder.transform(X)
