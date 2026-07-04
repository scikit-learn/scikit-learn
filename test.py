import numpy as np
import pytest
from scipy import sparse
import re
from sklearn.exceptions import NotFittedError
from sklearn.preprocessing import OneHotEncoder, OrdinalEncoder
from sklearn.utils.fixes import CSR_CONTAINERS
from sklearn.utils._testing import (
    _convert_container,
    assert_allclose,
    assert_array_equal,
)
def main_same():
    enc = OrdinalEncoder(
        handle_unknown="use_encoded_value", unknown_value=-2,categories='frequency',
    )

    X_fit = np.array(
        [
            ["a", "x"],
            ["a", "x"],
            ["a", "x"],
            ["a", "x"],
            ["b", "y"],
            ["b", "y"],
            ["b", "z"],
            ["b", "z"]
        ],
        dtype=object,
    )
    X_trans = np.array(
        [
            ["a", "xy"],
            ["bla", "y"],
            ["a", "x"],
            ["a", "x"],
            ["b", "x"],
            ["b", "x"],
            ["b", "x"],
            ["b", "x"]
        ],
        dtype=object,
    )
    enc.fit(X_fit)

    X_trans_enc = enc.transform(X_trans)
    print("X_trans_enc",X_trans_enc)
    exp = np.array(
        [
            [ 0.0, -2.0],
            [-2.0,  1.0],
            [ 0.0,  0.0],
            [ 0.0,  0.0],
            [ 1.0,  0.0],
            [ 1.0,  0.0],
            [ 1.0,  0.0],
            [ 1.0,  0.0]
        ], dtype="int64",
    )
    
    #assert_array_equal(X_trans_enc, exp)
    X_trans_inv = enc.inverse_transform(X_trans_enc)

def main():
    enc = OrdinalEncoder(
        handle_unknown="use_encoded_value", unknown_value=-2,categories='frequency',
    )
    #m  
    #same values
    #cover edge cases
    #handle ties correctly with weird config
    #handle error messages correctly.
    #unknown values

    X_fit = np.array(
        [
            ["b", "x"],
            ["b", "x"],
            ["b", "x"],
            ["b", "x"],
            ["a", "y"],
            ["a", "y"],
            ["a", "z"],
            ["a", "z"],
            ["g", "j"],
            ["n", "m"],
            ["z", "e"]
        ],
        dtype=object,
    )
    X_trans = np.array(
        [
            ["a", "xy"],
            ["bla", "y"],
            ["a", "x"],
            ["a", "x"],
            ["b", "x"],
            ["b", "x"],
            ["b", "x"],
            ["b", "x"]
        ],
        dtype=object,
    )
    enc.fit(X_fit)

    X_trans_enc = enc.transform(X_trans)
    print("X_trans_enc",X_trans_enc)
    exp = np.array(
        [
            [ 0.0, -2.0],
            [-2.0,  1.0],
            [ 0.0,  0.0],
            [ 0.0,  0.0],
            [ 1.0,  0.0],
            [ 1.0,  0.0],
            [ 1.0,  0.0],
            [ 1.0,  0.0]
        ], dtype="int64",
    )
    
    #assert_array_equal(X_trans_enc, exp)
    X_trans_inv = enc.inverse_transform(X_trans_enc)
    inv_exp = np.array(
        [
            ["a", None],
            [None, "y"],
            ["a", "x"],
            ["a", "x"],
            ["b", "x"],
            ["b", "x"],
            ["b" ,"x"],
            ["b" ,"x"]
        ],
        dtype=object,
    )
    print("X_trans_inv",X_trans_inv)

    assert_array_equal(X_trans_inv, inv_exp)

main()
# def test_ordinal_encoder_inverse():
#     X = [["abc", 2, 55], ["def", 1, 55]]
#     enc = OrdinalEncoder()
#     X_tr = enc.fit_transform(X)
#     exp = np.array(X, dtype=object)
#     assert_array_equal(enc.inverse_transform(X_tr), exp)

#     # incorrect shape raises
#     X_tr = np.array([[0, 1, 1, 2], [1, 0, 1, 0]])
#     msg = re.escape("Shape of the passed X data is not correct")
#     with pytest.raises(ValueError, match=msg):
#         enc.inverse_transform(X_tr)
#test_ordinal_encoder_inverse()

def test_ordinal_encoder_tie():
    # enc = OrdinalEncoder()
    # exp = np.array([[0, 1, 0], [1, 0, 0]], dtype="int64")
    # assert_array_equal(enc.fit_transform(X), exp.astype("float64"))
    # enc = OrdinalEncoder(dtype="int64")
    # assert_array_equal(enc.fit_transform(X), exp)
    enc = OrdinalEncoder(
        handle_unknown="use_encoded_value", unknown_value=-2,categories='frequency',
    )
    X_fit = np.array(
        [
            ["a", "x"],
            ["a", "x"],
            ["a", "x"],
            ["a", "x"],
            ["b", "y"],
            ["b", "y"],
            ["b", "y"],
            ["b", "y"]
        ],
        dtype=object,
    )
    X_trans = np.array(
        [   ["a", "x"],
            ["a", "x"],
            ["a", "x"],
            ["a", "x"],
            ["b", "y"],
            ["b", "y"],
            ["b", "y"],
            ["b", "y"]
        ],
        dtype=object,
    )
    enc.fit(X_fit)

    X_trans_enc = enc.transform(X_trans)
    print("X_trans_enc",X_trans_enc)
    exp = np.array(
        [
            [ 0.0,  0.0],
            [ 0.0,  0.0],
            [ 0.0,  0.0],
            [ 0.0,  0.0],
            [ 1.0,  1.0],
            [ 1.0,  1.0],
            [ 1.0,  1.0],
            [ 1.0,  1.0]
        ], dtype="int64",
    )
    
    assert_array_equal(X_trans_enc, exp)
    X_trans_inv = enc.inverse_transform(X_trans_enc)
    assert_array_equal(X_trans_inv, X_trans)

def test_ordinal_encoder_tie_inv_fit():
    # trying to break the fit by switching the order of fitting
    enc = OrdinalEncoder(
        handle_unknown="use_encoded_value", unknown_value=-2,categories='frequency',
    )
    X_fit = np.array(
        [

            ["b", "y"],
            ["b", "y"],
            ["b", "y"],
            ["b", "y"],
            ["a", "x"],
            ["a", "x"],
            ["a", "x"],
            ["a", "x"]
        ],
        dtype=object,
    )
    X_trans = np.array(
        [   ["a", "x"],
            ["a", "x"],
            ["a", "x"],
            ["a", "x"],
            ["b", "y"],
            ["b", "y"],
            ["b", "y"],
            ["b", "y"]
        ],
        dtype=object,
    )
    enc.fit(X_fit)

    X_trans_enc = enc.transform(X_trans)
    print("X_trans_enc",X_trans_enc)
    exp = np.array(
        [
            [ 0.0,  0.0],
            [ 0.0,  0.0],
            [ 0.0,  0.0],
            [ 0.0,  0.0],
            [ 1.0,  1.0],
            [ 1.0,  1.0],
            [ 1.0,  1.0],
            [ 1.0,  1.0]
        ], dtype="int64",
    )
    
    assert_array_equal(X_trans_enc, exp)
    X_trans_inv = enc.inverse_transform(X_trans_enc)
    assert_array_equal(X_trans_inv, X_trans)

def test_ordinal_encoder_none():
    # enc = OrdinalEncoder()
    # exp = np.array([[0, 1, 0], [1, 0, 0]], dtype="int64")
    # assert_array_equal(enc.fit_transform(X), exp.astype("float64"))
    # enc = OrdinalEncoder(dtype="int64")
    # assert_array_equal(enc.fit_transform(X), exp)
    enc = OrdinalEncoder(
        handle_unknown="use_encoded_value", unknown_value=-2,categories='frequency',
    )
    X_fit = np.array(
        [
            ["a", "x"],
            ["a", "x"],
            ["a", "x"],
            ["a", "x"],
            ["b", "y"],
            ["b", "y"],
            ["b", "y"],
            ["b", "y"]
        ],
        dtype=object,
    )
    X_trans = np.array(
        [   ["a", "x"],
            ["ab", "x"],
            ["a", "x"],
            ["a", "x"],
            ["b", "y"],
            ["b", "y"],
            ["b", "y"],
            ["b", "y"]
        ],
        dtype=object,
    )
    enc.fit(X_fit)

    X_trans_enc = enc.transform(X_trans)
    print("X_trans_enc",X_trans_enc)
    exp = np.array(
        [
            [ 0.0,  0.0],
            [ -2.0,  0.0],
            [ 0.0,  0.0],
            [ 0.0,  0.0],
            [ 1.0,  1.0],
            [ 1.0,  1.0],
            [ 1.0,  1.0],
            [ 1.0,  1.0]
        ], dtype="int64",
    )
    
    assert_array_equal(X_trans_enc, exp)
    X_trans_inv = enc.inverse_transform(X_trans_enc)
   #assert_array_equal(X_trans_inv, X_trans)

def test_ordinal_encoder_multy_letters():
    # enc = OrdinalEncoder()
    # exp = np.array([[0, 1, 0], [1, 0, 0]], dtype="int64")
    # assert_array_equal(enc.fit_transform(X), exp.astype("float64"))
    # enc = OrdinalEncoder(dtype="int64")
    # assert_array_equal(enc.fit_transform(X), exp)
    enc = OrdinalEncoder(
        handle_unknown="use_encoded_value", unknown_value=-2,categories='frequency',
    )
    X_fit = np.array(
        [
            ["aa", "xx"],
            ["aa", "xx"],
            ["aa", "xx"],
            ["aa", "xx"],
            ["bb", "yy"],
            ["bb", "yy"],
            ["bb", "yy"],
            ["bb", "yy"]
        ],
        dtype=object,
    )
    X_trans = np.array(
        [   ["aa", "xx"],
            ["aa", "xx"],
            ["aa", "xx"],
            ["aa", "xx"],
            ["bb", "yy"],
            ["bb", "yy"],
            ["bb", "yy"],
            ["bb", "yy"]
        ],
        dtype=object,
    )
    enc.fit(X_fit)

    X_trans_enc = enc.transform(X_trans)
    print("X_trans_enc",X_trans_enc)
    exp = np.array(
        [
            [ 0.0,  0.0],
            [ 0.0,  0.0],
            [ 0.0,  0.0],
            [ 0.0,  0.0],
            [ 1.0,  1.0],
            [ 1.0,  1.0],
            [ 1.0,  1.0],
            [ 1.0,  1.0]
        ], dtype="int64",
    )
    
    assert_array_equal(X_trans_enc, exp)
    X_trans_inv = enc.inverse_transform(X_trans_enc)
    assert_array_equal(X_trans_inv, X_trans)

def test_ordinal_encoder_numbers():
    # enc = OrdinalEncoder()
    # exp = np.array([[0, 1, 0], [1, 0, 0]], dtype="int64")
    # assert_array_equal(enc.fit_transform(X), exp.astype("float64"))
    # enc = OrdinalEncoder(dtype="int64")
    # assert_array_equal(enc.fit_transform(X), exp)
    enc = OrdinalEncoder(
        handle_unknown="use_encoded_value", unknown_value=-2,categories='frequency',
    )
    X_fit = np.array(
        [
            [2, 7],
            [2, 7],
            [2, 7],
            [1, 5],
            [1, 5],
            [1, 5],
            [1, 5]
        ],
        dtype=object,
    )
    X_trans = np.array(
        [   [2, 7],
            [2, 7],
            [2, 7],
            [2, 7],
            [1, 5],
            [1, 5],
            [1, 5],
            [1, 5],
        ],
        dtype=object,
    )
    enc.fit(X_fit)

    X_trans_enc = enc.transform(X_trans)
    print("X_trans_enc",X_trans_enc)
    exp = np.array(
        [
            [ 0.0,  0.0],
            [ 0.0,  0.0],
            [ 0.0,  0.0],
            [ 0.0,  0.0],
            [ 1.0,  1.0],
            [ 1.0,  1.0],
            [ 1.0,  1.0],
            [ 1.0,  1.0]
        ], dtype="int64",
    )
    
    assert_array_equal(X_trans_enc, exp)
    X_trans_inv = enc.inverse_transform(X_trans_enc)
    assert_array_equal(X_trans_inv, X_trans)
###########################

def test_ordinal_encoder_inverse():
    X = [["abc", 2, 55], ["def", 1, 55]]
    enc = OrdinalEncoder(categories='frequency')
    X_tr = enc.fit_transform(X)
    exp = np.array(X, dtype=object)
    assert_array_equal(enc.inverse_transform(X_tr), exp)

    # incorrect shape raises
    X_tr = np.array([[0, 1, 1, 2], [1, 0, 1, 0]])
    msg = re.escape("Shape of the passed X data is not correct")
    with pytest.raises(ValueError, match=msg):
        enc.inverse_transform(X_tr)


def test_ordinal_encoder_handle_unknowns_string():
    enc = OrdinalEncoder(handle_unknown="use_encoded_value", unknown_value=-2,categories='frequency')
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
    enc = OrdinalEncoder(handle_unknown="use_encoded_value", unknown_value=-999,categories='frequency')
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

    enc = OrdinalEncoder(handle_unknown="use_encoded_value", unknown_value=np.nan,categories='frequency')

    X_fit = np.array([[1], [2], [3]])
    enc.fit(X_fit)
    X_trans = enc.transform([[1], [2], [4]])
    assert_array_equal(X_trans, [[0], [1], [np.nan]])


def test_ordinal_encoder_handle_unknowns_nan_non_float_dtype():
    # Make sure an error is raised when unknown_value=np.nan and the dtype
    # isn't a float dtype
    enc = OrdinalEncoder(
        handle_unknown="use_encoded_value", unknown_value=np.nan, dtype=int,categories='frequency'
    )

    X_fit = np.array([[1], [2], [3]])
    with pytest.raises(ValueError, match="dtype parameter should be a float dtype"):
        enc.fit(X_fit)










@pytest.mark.parametrize(
    "X, expected_X_trans, X_test",
    [
        (
            np.array([[1.0, np.nan, 3.0]]).T,
            np.array([[0.0, np.nan, 1.0]]).T,
            np.array([[4.0]]),
        ),
        (
            np.array([[1.0, 4.0, 3.0]]).T,
            np.array([[0.0, 2.0, 1.0]]).T,
            np.array([[np.nan]]),
        ),
        (
            np.array([["c", np.nan, "b"]], dtype=object).T,
            np.array([[1.0, np.nan, 0.0]]).T,
            np.array([["d"]], dtype=object),
        ),
        (
            np.array([["c", "a", "b"]], dtype=object).T,
            np.array([[2.0, 0.0, 1.0]]).T,
            np.array([[np.nan]], dtype=object),
        ),
    ],
)
def test_ordinal_encoder_handle_missing_and_unknown(X, expected_X_trans, X_test):
    """Test the interaction between missing values and handle_unknown"""

    oe = OrdinalEncoder(handle_unknown="use_encoded_value",categories='frequency', unknown_value=-1)

    X_trans = oe.fit_transform(X)
    assert_allclose(X_trans, expected_X_trans)

    assert_allclose(oe.transform(X_test), [[-1.0]])


@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
def test_ordinal_encoder_sparse(csr_container):
    """Check that we raise proper error with sparse input in OrdinalEncoder.
    Non-regression test for:
    https://github.com/scikit-learn/scikit-learn/issues/19878
    """
    X = np.array([[3, 2, 1], [0, 1, 1]])
    X_sparse = csr_container(X)

    encoder = OrdinalEncoder(categories='frequency')

    err_msg = "Sparse data was passed, but dense data is required"
    with pytest.raises(TypeError, match=err_msg):
        encoder.fit(X_sparse)
    with pytest.raises(TypeError, match=err_msg):
        encoder.fit_transform(X_sparse)

    X_trans = encoder.fit_transform(X)
    X_trans_sparse = csr_container(X_trans)
    with pytest.raises(TypeError, match=err_msg):
        encoder.inverse_transform(X_trans_sparse)



@pytest.mark.parametrize(
    "X_train",
    [
        [["AA", "B"]],
        np.array([["AA", "B"]], dtype="O"),
        np.array([["AA", "B"]], dtype="U"),
    ],
)
@pytest.mark.parametrize(
    "X_test",
    [
        [["A", "B"]],
        np.array([["A", "B"]], dtype="O"),
        np.array([["A", "B"]], dtype="U"),
    ],
)
def test_ordinal_encoder_handle_unknown_string_dtypes(X_train, X_test):
    """Checks that `OrdinalEncoder` transforms string dtypes.
    Non-regression test for:
    https://github.com/scikit-learn/scikit-learn/issues/19872
    """
    enc = OrdinalEncoder(handle_unknown="use_encoded_value",categories='frequency', unknown_value=-9)
    enc.fit(X_train)

    X_trans = enc.transform(X_test)
    assert_allclose(X_trans, [[-9, 0]])


def test_ordinal_encoder_python_integer():
    """Check that `OrdinalEncoder` accepts Python integers that are potentially
    larger than 64 bits.
    Non-regression test for:
    https://github.com/scikit-learn/scikit-learn/issues/20721
    """
    X = np.array(
        [
            44253463435747313673,
            9867966753463435747313673,
            44253462342215747313673,
            442534634357764313673,
        ]
    ).reshape(-1, 1)
    encoder = OrdinalEncoder(categories='frequency').fit(X)
    assert_array_equal(encoder.categories_, np.sort(X, axis=0).T)
    X_trans = encoder.transform(X)
    assert_array_equal(X_trans, [[0], [3], [2], [1]])


def test_ordinal_encoder_features_names_out_pandas():
    """Check feature names out is same as the input."""
    pd = pytest.importorskip("pandas")

    names = ["b", "c", "a"]
    X = pd.DataFrame([[1, 2, 3]], columns=names)
    enc = OrdinalEncoder(categories='frequency').fit(X)

    feature_names_out = enc.get_feature_names_out()
    assert_array_equal(names, feature_names_out)


def test_ordinal_encoder_passthrough_missing_values_float_errors_dtype():
    """Test ordinal encoder with nan passthrough fails when dtype=np.int32."""

    X = np.array([[np.nan, 3.0, 1.0, 3.0]]).T
    oe = OrdinalEncoder(dtype=np.int32,categories='frequency')

    msg = (
        r"There are missing values in features \[0\]. For OrdinalEncoder "
        f"to encode missing values with dtype: {np.int32}"
    )
    with pytest.raises(ValueError, match=msg):
        oe.fit(X)

