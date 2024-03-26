import warnings
from copy import copy
from unittest import SkipTest

import numpy as np
import pytest

import sklearn
from sklearn.externals._packaging.version import parse as parse_version
from sklearn.utils import _safe_indexing, resample, shuffle
from sklearn.utils._array_api import yield_namespace_device_dtype_combinations
from sklearn.utils._indexing import (
    _determine_key_type,
    _get_column_indices,
    _safe_assign,
)
from sklearn.utils._mocking import MockDataFrame
from sklearn.utils._testing import (
    _array_api_for_tests,
    _convert_container,
    assert_allclose_dense_sparse,
    assert_array_equal,
    skip_if_array_api_compat_not_configured,
)
from sklearn.utils.fixes import CSC_CONTAINERS, CSR_CONTAINERS

# toy array
X_toy = np.arange(9).reshape((3, 3))


def test_polars_indexing():
    """Check _safe_indexing for polars as expected."""
    pl = pytest.importorskip("polars", minversion="0.18.2")
    df = pl.DataFrame(
        {"a": [1, 2, 3, 4], "b": [4, 5, 6, 8], "c": [1, 4, 1, 10]}, orient="row"
    )

    from polars.testing import assert_frame_equal

    str_keys = [["b"], ["a", "b"], ["b", "a", "c"], ["c"], ["a"]]

    for key in str_keys:
        out = _safe_indexing(df, key, axis=1)
        assert_frame_equal(df[key], out)

    bool_keys = [([True, False, True], ["a", "c"]), ([False, False, True], ["c"])]

    for bool_key, str_key in bool_keys:
        out = _safe_indexing(df, bool_key, axis=1)
        assert_frame_equal(df[:, str_key], out)

    int_keys = [([0, 1], ["a", "b"]), ([2], ["c"])]

    for int_key, str_key in int_keys:
        out = _safe_indexing(df, int_key, axis=1)
        assert_frame_equal(df[:, str_key], out)

    axis_0_keys = [[0, 1], [1, 3], [3, 2]]
    for key in axis_0_keys:
        out = _safe_indexing(df, key, axis=0)
        assert_frame_equal(df[key], out)


@pytest.mark.parametrize(
    "key, dtype",
    [
        (0, "int"),
        ("0", "str"),
        (True, "bool"),
        (np.bool_(True), "bool"),
        ([0, 1, 2], "int"),
        (["0", "1", "2"], "str"),
        ((0, 1, 2), "int"),
        (("0", "1", "2"), "str"),
        (slice(None, None), None),
        (slice(0, 2), "int"),
        (np.array([0, 1, 2], dtype=np.int32), "int"),
        (np.array([0, 1, 2], dtype=np.int64), "int"),
        (np.array([0, 1, 2], dtype=np.uint8), "int"),
        ([True, False], "bool"),
        ((True, False), "bool"),
        (np.array([True, False]), "bool"),
        ("col_0", "str"),
        (["col_0", "col_1", "col_2"], "str"),
        (("col_0", "col_1", "col_2"), "str"),
        (slice("begin", "end"), "str"),
        (np.array(["col_0", "col_1", "col_2"]), "str"),
        (np.array(["col_0", "col_1", "col_2"], dtype=object), "str"),
    ],
)
def test_determine_key_type(key, dtype):
    assert _determine_key_type(key) == dtype


def test_determine_key_type_error():
    with pytest.raises(ValueError, match="No valid specification of the"):
        _determine_key_type(1.0)


def test_determine_key_type_slice_error():
    with pytest.raises(TypeError, match="Only array-like or scalar are"):
        _determine_key_type(slice(0, 2, 1), accept_slice=False)


@skip_if_array_api_compat_not_configured
@pytest.mark.parametrize(
    "array_namespace, device, dtype_name", yield_namespace_device_dtype_combinations()
)
def test_determine_key_type_array_api(array_namespace, device, dtype_name):
    xp = _array_api_for_tests(array_namespace, device)

    with sklearn.config_context(array_api_dispatch=True):
        int_array_key = xp.asarray([1, 2, 3])
        assert _determine_key_type(int_array_key) == "int"

        bool_array_key = xp.asarray([True, False, True])
        assert _determine_key_type(bool_array_key) == "bool"

        try:
            complex_array_key = xp.asarray([1 + 1j, 2 + 2j, 3 + 3j])
        except TypeError:
            # Complex numbers are not supported by all Array API libraries.
            complex_array_key = None

        if complex_array_key is not None:
            with pytest.raises(ValueError, match="No valid specification of the"):
                _determine_key_type(complex_array_key)


@pytest.mark.parametrize(
    "array_type", ["list", "array", "sparse", "dataframe", "polars"]
)
@pytest.mark.parametrize("indices_type", ["list", "tuple", "array", "series", "slice"])
def test_safe_indexing_2d_container_axis_0(array_type, indices_type):
    indices = [1, 2]
    if indices_type == "slice" and isinstance(indices[1], int):
        indices[1] += 1
    array = _convert_container([[1, 2, 3], [4, 5, 6], [7, 8, 9]], array_type)
    indices = _convert_container(indices, indices_type)
    subset = _safe_indexing(array, indices, axis=0)
    assert_allclose_dense_sparse(
        subset, _convert_container([[4, 5, 6], [7, 8, 9]], array_type)
    )


@pytest.mark.parametrize("array_type", ["list", "array", "series", "polars_series"])
@pytest.mark.parametrize("indices_type", ["list", "tuple", "array", "series", "slice"])
def test_safe_indexing_1d_container(array_type, indices_type):
    indices = [1, 2]
    if indices_type == "slice" and isinstance(indices[1], int):
        indices[1] += 1
    array = _convert_container([1, 2, 3, 4, 5, 6, 7, 8, 9], array_type)
    indices = _convert_container(indices, indices_type)
    subset = _safe_indexing(array, indices, axis=0)
    assert_allclose_dense_sparse(subset, _convert_container([2, 3], array_type))


@pytest.mark.parametrize("array_type", ["array", "sparse", "dataframe", "polars"])
@pytest.mark.parametrize("indices_type", ["list", "tuple", "array", "series", "slice"])
@pytest.mark.parametrize("indices", [[1, 2], ["col_1", "col_2"]])
def test_safe_indexing_2d_container_axis_1(array_type, indices_type, indices):
    # validation of the indices
    # we make a copy because indices is mutable and shared between tests
    indices_converted = copy(indices)
    if indices_type == "slice" and isinstance(indices[1], int):
        indices_converted[1] += 1

    columns_name = ["col_0", "col_1", "col_2"]
    array = _convert_container(
        [[1, 2, 3], [4, 5, 6], [7, 8, 9]], array_type, columns_name
    )
    indices_converted = _convert_container(indices_converted, indices_type)

    if isinstance(indices[0], str) and array_type not in ("dataframe", "polars"):
        err_msg = (
            "Specifying the columns using strings is only supported for dataframes"
        )
        with pytest.raises(ValueError, match=err_msg):
            _safe_indexing(array, indices_converted, axis=1)
    else:
        subset = _safe_indexing(array, indices_converted, axis=1)
        assert_allclose_dense_sparse(
            subset, _convert_container([[2, 3], [5, 6], [8, 9]], array_type)
        )


@pytest.mark.parametrize("array_read_only", [True, False])
@pytest.mark.parametrize("indices_read_only", [True, False])
@pytest.mark.parametrize("array_type", ["array", "sparse", "dataframe", "polars"])
@pytest.mark.parametrize("indices_type", ["array", "series"])
@pytest.mark.parametrize(
    "axis, expected_array", [(0, [[4, 5, 6], [7, 8, 9]]), (1, [[2, 3], [5, 6], [8, 9]])]
)
def test_safe_indexing_2d_read_only_axis_1(
    array_read_only, indices_read_only, array_type, indices_type, axis, expected_array
):
    array = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    if array_read_only:
        array.setflags(write=False)
    array = _convert_container(array, array_type)
    indices = np.array([1, 2])
    if indices_read_only:
        indices.setflags(write=False)
    indices = _convert_container(indices, indices_type)
    subset = _safe_indexing(array, indices, axis=axis)
    assert_allclose_dense_sparse(subset, _convert_container(expected_array, array_type))


@pytest.mark.parametrize("array_type", ["list", "array", "series", "polars_series"])
@pytest.mark.parametrize("indices_type", ["list", "tuple", "array", "series"])
def test_safe_indexing_1d_container_mask(array_type, indices_type):
    indices = [False] + [True] * 2 + [False] * 6
    array = _convert_container([1, 2, 3, 4, 5, 6, 7, 8, 9], array_type)
    indices = _convert_container(indices, indices_type)
    subset = _safe_indexing(array, indices, axis=0)
    assert_allclose_dense_sparse(subset, _convert_container([2, 3], array_type))


@pytest.mark.parametrize("array_type", ["array", "sparse", "dataframe", "polars"])
@pytest.mark.parametrize("indices_type", ["list", "tuple", "array", "series"])
@pytest.mark.parametrize(
    "axis, expected_subset",
    [(0, [[4, 5, 6], [7, 8, 9]]), (1, [[2, 3], [5, 6], [8, 9]])],
)
def test_safe_indexing_2d_mask(array_type, indices_type, axis, expected_subset):
    columns_name = ["col_0", "col_1", "col_2"]
    array = _convert_container(
        [[1, 2, 3], [4, 5, 6], [7, 8, 9]], array_type, columns_name
    )
    indices = [False, True, True]
    indices = _convert_container(indices, indices_type)

    subset = _safe_indexing(array, indices, axis=axis)
    assert_allclose_dense_sparse(
        subset, _convert_container(expected_subset, array_type)
    )


@pytest.mark.parametrize(
    "array_type, expected_output_type",
    [
        ("list", "list"),
        ("array", "array"),
        ("sparse", "sparse"),
        ("dataframe", "series"),
        ("polars", "polars_series"),
    ],
)
def test_safe_indexing_2d_scalar_axis_0(array_type, expected_output_type):
    array = _convert_container([[1, 2, 3], [4, 5, 6], [7, 8, 9]], array_type)
    indices = 2
    subset = _safe_indexing(array, indices, axis=0)
    expected_array = _convert_container([7, 8, 9], expected_output_type)
    assert_allclose_dense_sparse(subset, expected_array)


@pytest.mark.parametrize("array_type", ["list", "array", "series", "polars_series"])
def test_safe_indexing_1d_scalar(array_type):
    array = _convert_container([1, 2, 3, 4, 5, 6, 7, 8, 9], array_type)
    indices = 2
    subset = _safe_indexing(array, indices, axis=0)
    assert subset == 3


@pytest.mark.parametrize(
    "array_type, expected_output_type",
    [
        ("array", "array"),
        ("sparse", "sparse"),
        ("dataframe", "series"),
        ("polars", "polars_series"),
    ],
)
@pytest.mark.parametrize("indices", [2, "col_2"])
def test_safe_indexing_2d_scalar_axis_1(array_type, expected_output_type, indices):
    columns_name = ["col_0", "col_1", "col_2"]
    array = _convert_container(
        [[1, 2, 3], [4, 5, 6], [7, 8, 9]], array_type, columns_name
    )

    if isinstance(indices, str) and array_type not in ("dataframe", "polars"):
        err_msg = (
            "Specifying the columns using strings is only supported for dataframes"
        )
        with pytest.raises(ValueError, match=err_msg):
            _safe_indexing(array, indices, axis=1)
    else:
        subset = _safe_indexing(array, indices, axis=1)
        expected_output = [3, 6, 9]
        if expected_output_type == "sparse":
            # sparse matrix are keeping the 2D shape
            expected_output = [[3], [6], [9]]
        expected_array = _convert_container(expected_output, expected_output_type)
        assert_allclose_dense_sparse(subset, expected_array)


@pytest.mark.parametrize("array_type", ["list", "array", "sparse"])
def test_safe_indexing_None_axis_0(array_type):
    X = _convert_container([[1, 2, 3], [4, 5, 6], [7, 8, 9]], array_type)
    X_subset = _safe_indexing(X, None, axis=0)
    assert_allclose_dense_sparse(X_subset, X)


def test_safe_indexing_pandas_no_matching_cols_error():
    pd = pytest.importorskip("pandas")
    err_msg = "No valid specification of the columns."
    X = pd.DataFrame(X_toy)
    with pytest.raises(ValueError, match=err_msg):
        _safe_indexing(X, [1.0], axis=1)


@pytest.mark.parametrize("axis", [None, 3])
def test_safe_indexing_error_axis(axis):
    with pytest.raises(ValueError, match="'axis' should be either 0"):
        _safe_indexing(X_toy, [0, 1], axis=axis)


@pytest.mark.parametrize("X_constructor", ["array", "series", "polars_series"])
def test_safe_indexing_1d_array_error(X_constructor):
    # check that we are raising an error if the array-like passed is 1D and
    # we try to index on the 2nd dimension
    X = list(range(5))
    if X_constructor == "array":
        X_constructor = np.asarray(X)
    elif X_constructor == "series":
        pd = pytest.importorskip("pandas")
        X_constructor = pd.Series(X)
    elif X_constructor == "polars_series":
        pl = pytest.importorskip("polars")
        X_constructor = pl.Series(values=X)

    err_msg = "'X' should be a 2D NumPy array, 2D sparse matrix or dataframe"
    with pytest.raises(ValueError, match=err_msg):
        _safe_indexing(X_constructor, [0, 1], axis=1)


def test_safe_indexing_container_axis_0_unsupported_type():
    indices = ["col_1", "col_2"]
    array = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    err_msg = "String indexing is not supported with 'axis=0'"
    with pytest.raises(ValueError, match=err_msg):
        _safe_indexing(array, indices, axis=0)


def test_safe_indexing_pandas_no_settingwithcopy_warning():
    # Using safe_indexing with an array-like indexer gives a copy of the
    # DataFrame -> ensure it doesn't raise a warning if modified
    pd = pytest.importorskip("pandas")

    pd_version = parse_version(pd.__version__)
    pd_base_version = parse_version(pd_version.base_version)

    if pd_base_version >= parse_version("3"):
        raise SkipTest("SettingWithCopyWarning has been removed in pandas 3.0.0.dev")

    X = pd.DataFrame({"a": [1, 2, 3], "b": [3, 4, 5]})
    subset = _safe_indexing(X, [0, 1], axis=0)
    if hasattr(pd.errors, "SettingWithCopyWarning"):
        SettingWithCopyWarning = pd.errors.SettingWithCopyWarning
    else:
        # backward compatibility for pandas < 1.5
        SettingWithCopyWarning = pd.core.common.SettingWithCopyWarning
    with warnings.catch_warnings():
        warnings.simplefilter("error", SettingWithCopyWarning)
        subset.iloc[0, 0] = 10
    # The original dataframe is unaffected by the assignment on the subset:
    assert X.iloc[0, 0] == 1


@pytest.mark.parametrize("indices", [0, [0, 1], slice(0, 2), np.array([0, 1])])
def test_safe_indexing_list_axis_1_unsupported(indices):
    """Check that we raise a ValueError when axis=1 with input as list."""
    X = [[1, 2], [4, 5], [7, 8]]
    err_msg = "axis=1 is not supported for lists"
    with pytest.raises(ValueError, match=err_msg):
        _safe_indexing(X, indices, axis=1)


@pytest.mark.parametrize("array_type", ["array", "sparse", "dataframe"])
def test_safe_assign(array_type):
    """Check that `_safe_assign` works as expected."""
    rng = np.random.RandomState(0)
    X_array = rng.randn(10, 5)

    row_indexer = [1, 2]
    values = rng.randn(len(row_indexer), X_array.shape[1])
    X = _convert_container(X_array, array_type)
    _safe_assign(X, values, row_indexer=row_indexer)

    assigned_portion = _safe_indexing(X, row_indexer, axis=0)
    assert_allclose_dense_sparse(
        assigned_portion, _convert_container(values, array_type)
    )

    column_indexer = [1, 2]
    values = rng.randn(X_array.shape[0], len(column_indexer))
    X = _convert_container(X_array, array_type)
    _safe_assign(X, values, column_indexer=column_indexer)

    assigned_portion = _safe_indexing(X, column_indexer, axis=1)
    assert_allclose_dense_sparse(
        assigned_portion, _convert_container(values, array_type)
    )

    row_indexer, column_indexer = None, None
    values = rng.randn(*X.shape)
    X = _convert_container(X_array, array_type)
    _safe_assign(X, values, column_indexer=column_indexer)

    assert_allclose_dense_sparse(X, _convert_container(values, array_type))


@pytest.mark.parametrize(
    "key, err_msg",
    [
        (10, r"all features must be in \[0, 2\]"),
        ("whatever", "A given column is not a column of the dataframe"),
        (object(), "No valid specification of the columns"),
    ],
)
def test_get_column_indices_error(key, err_msg):
    pd = pytest.importorskip("pandas")
    X_df = pd.DataFrame(X_toy, columns=["col_0", "col_1", "col_2"])

    with pytest.raises(ValueError, match=err_msg):
        _get_column_indices(X_df, key)


@pytest.mark.parametrize(
    "key", [["col1"], ["col2"], ["col1", "col2"], ["col1", "col3"], ["col2", "col3"]]
)
def test_get_column_indices_pandas_nonunique_columns_error(key):
    pd = pytest.importorskip("pandas")
    toy = np.zeros((1, 5), dtype=int)
    columns = ["col1", "col1", "col2", "col3", "col2"]
    X = pd.DataFrame(toy, columns=columns)

    err_msg = "Selected columns, {}, are not unique in dataframe".format(key)
    with pytest.raises(ValueError) as exc_info:
        _get_column_indices(X, key)
    assert str(exc_info.value) == err_msg


def test_get_column_indices_interchange():
    """Check _get_column_indices for edge cases with the interchange"""
    pd = pytest.importorskip("pandas", minversion="1.5")

    df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], columns=["a", "b", "c"])

    # Hide the fact that this is a pandas dataframe to trigger the dataframe protocol
    # code path.
    class MockDataFrame:
        def __init__(self, df):
            self._df = df

        def __getattr__(self, name):
            return getattr(self._df, name)

    df_mocked = MockDataFrame(df)

    key_results = [
        (slice(1, None), [1, 2]),
        (slice(None, 2), [0, 1]),
        (slice(1, 2), [1]),
        (["b", "c"], [1, 2]),
        (slice("a", "b"), [0, 1]),
        (slice("a", None), [0, 1, 2]),
        (slice(None, "a"), [0]),
        (["c", "a"], [2, 0]),
        ([], []),
    ]
    for key, result in key_results:
        assert _get_column_indices(df_mocked, key) == result

    msg = "A given column is not a column of the dataframe"
    with pytest.raises(ValueError, match=msg):
        _get_column_indices(df_mocked, ["not_a_column"])

    msg = "key.step must be 1 or None"
    with pytest.raises(NotImplementedError, match=msg):
        _get_column_indices(df_mocked, slice("a", None, 2))


def test_resample():
    # Border case not worth mentioning in doctests
    assert resample() is None

    # Check that invalid arguments yield ValueError
    with pytest.raises(ValueError):
        resample([0], [0, 1])
    with pytest.raises(ValueError):
        resample([0, 1], [0, 1], replace=False, n_samples=3)

    # Issue:6581, n_samples can be more when replace is True (default).
    assert len(resample([1, 2], n_samples=5)) == 5


def test_resample_stratified():
    # Make sure resample can stratify
    rng = np.random.RandomState(0)
    n_samples = 100
    p = 0.9
    X = rng.normal(size=(n_samples, 1))
    y = rng.binomial(1, p, size=n_samples)

    _, y_not_stratified = resample(X, y, n_samples=10, random_state=0, stratify=None)
    assert np.all(y_not_stratified == 1)

    _, y_stratified = resample(X, y, n_samples=10, random_state=0, stratify=y)
    assert not np.all(y_stratified == 1)
    assert np.sum(y_stratified) == 9  # all 1s, one 0


def test_resample_stratified_replace():
    # Make sure stratified resampling supports the replace parameter
    rng = np.random.RandomState(0)
    n_samples = 100
    X = rng.normal(size=(n_samples, 1))
    y = rng.randint(0, 2, size=n_samples)

    X_replace, _ = resample(
        X, y, replace=True, n_samples=50, random_state=rng, stratify=y
    )
    X_no_replace, _ = resample(
        X, y, replace=False, n_samples=50, random_state=rng, stratify=y
    )
    assert np.unique(X_replace).shape[0] < 50
    assert np.unique(X_no_replace).shape[0] == 50

    # make sure n_samples can be greater than X.shape[0] if we sample with
    # replacement
    X_replace, _ = resample(
        X, y, replace=True, n_samples=1000, random_state=rng, stratify=y
    )
    assert X_replace.shape[0] == 1000
    assert np.unique(X_replace).shape[0] == 100


def test_resample_stratify_2dy():
    # Make sure y can be 2d when stratifying
    rng = np.random.RandomState(0)
    n_samples = 100
    X = rng.normal(size=(n_samples, 1))
    y = rng.randint(0, 2, size=(n_samples, 2))
    X, y = resample(X, y, n_samples=50, random_state=rng, stratify=y)
    assert y.ndim == 2


@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
def test_resample_stratify_sparse_error(csr_container):
    # resample must be ndarray
    rng = np.random.RandomState(0)
    n_samples = 100
    X = rng.normal(size=(n_samples, 2))
    y = rng.randint(0, 2, size=n_samples)
    stratify = csr_container(y.reshape(-1, 1))
    with pytest.raises(TypeError, match="Sparse data was passed"):
        X, y = resample(X, y, n_samples=50, random_state=rng, stratify=stratify)


def test_shuffle_on_ndim_equals_three():
    def to_tuple(A):  # to make the inner arrays hashable
        return tuple(tuple(tuple(C) for C in B) for B in A)

    A = np.array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]])  # A.shape = (2,2,2)
    S = set(to_tuple(A))
    shuffle(A)  # shouldn't raise a ValueError for dim = 3
    assert set(to_tuple(A)) == S


@pytest.mark.parametrize("csc_container", CSC_CONTAINERS)
def test_shuffle_dont_convert_to_array(csc_container):
    # Check that shuffle does not try to convert to numpy arrays with float
    # dtypes can let any indexable datastructure pass-through.
    a = ["a", "b", "c"]
    b = np.array(["a", "b", "c"], dtype=object)
    c = [1, 2, 3]
    d = MockDataFrame(np.array([["a", 0], ["b", 1], ["c", 2]], dtype=object))
    e = csc_container(np.arange(6).reshape(3, 2))
    a_s, b_s, c_s, d_s, e_s = shuffle(a, b, c, d, e, random_state=0)

    assert a_s == ["c", "b", "a"]
    assert type(a_s) == list  # noqa: E721

    assert_array_equal(b_s, ["c", "b", "a"])
    assert b_s.dtype == object

    assert c_s == [3, 2, 1]
    assert type(c_s) == list  # noqa: E721

    assert_array_equal(d_s, np.array([["c", 2], ["b", 1], ["a", 0]], dtype=object))
    assert type(d_s) == MockDataFrame  # noqa: E721

    assert_array_equal(e_s.toarray(), np.array([[4, 5], [2, 3], [0, 1]]))
