"""Test the openml loader."""
import gzip
import json
import os
import re
import warnings
from functools import partial
from importlib import resources
from io import BytesIO
from urllib.error import HTTPError

import numpy as np
import scipy.sparse
import pytest

import sklearn
from sklearn import config_context
from sklearn.utils import Bunch
from sklearn.utils._testing import assert_allclose, assert_array_equal
from sklearn.utils._testing import fails_if_pypy

from sklearn.datasets import fetch_openml as fetch_openml_orig
from sklearn.datasets._openml import (
    _open_openml_url,
    _get_local_path,
    _retry_with_clean_cache,
)
from sklearn.datasets.tests.test_common import check_return_X_y


OPENML_TEST_DATA_MODULE = "sklearn.datasets.tests.data.openml"
# if True, urlopen will be monkey patched to only use local files
test_offline = True


# Do not use a cache for `fetch_openml` to avoid concurrent writing
# issues with `pytest-xdist`.
# Furthermore sklearn/datasets/tests/data/openml/ is not always consistent
# with the version on openml.org. If one were to load the dataset outside of
# the tests, it may result in data that does not represent openml.org.
fetch_openml = partial(fetch_openml_orig, data_home=None)


def _fetch_dataset_from_openml(
    data_id,
    data_name,
    data_version,
    target_column,
    expected_observations,
    expected_features,
    expected_missing,
    expected_data_dtype,
    expected_target_dtype,
    expect_sparse,
    compare_default_target,
):
    # fetches a dataset in three various ways from OpenML, using the
    # fetch_openml function, and does various checks on the validity of the
    # result. Note that this function can be mocked (by invoking
    # _monkey_patch_webbased_functions before invoking this function)
    data_by_name_id = fetch_openml(
        name=data_name, version=data_version, cache=False, as_frame=False
    )
    assert int(data_by_name_id.details["id"]) == data_id

    # Please note that cache=False is crucial, as the monkey patched files are
    # not consistent with reality
    with warnings.catch_warnings():
        # See discussion in PR #19373
        # Catching UserWarnings about multiple versions of dataset
        warnings.simplefilter("ignore", category=UserWarning)
        fetch_openml(name=data_name, cache=False, as_frame=False)
    # without specifying the version, there is no guarantee that the data id
    # will be the same

    # fetch with dataset id
    data_by_id = fetch_openml(
        data_id=data_id, cache=False, target_column=target_column, as_frame=False
    )
    assert data_by_id.details["name"] == data_name
    assert data_by_id.data.shape == (expected_observations, expected_features)
    if isinstance(target_column, str):
        # single target, so target is vector
        assert data_by_id.target.shape == (expected_observations,)
        assert data_by_id.target_names == [target_column]
    elif isinstance(target_column, list):
        # multi target, so target is array
        assert data_by_id.target.shape == (expected_observations, len(target_column))
        assert data_by_id.target_names == target_column
    assert data_by_id.data.dtype == expected_data_dtype
    assert data_by_id.target.dtype == expected_target_dtype
    assert len(data_by_id.feature_names) == expected_features
    for feature in data_by_id.feature_names:
        assert isinstance(feature, str)

    # TODO: pass in a list of expected nominal features
    for feature, categories in data_by_id.categories.items():
        feature_idx = data_by_id.feature_names.index(feature)

        # TODO: Remove when https://github.com/numpy/numpy/issues/19300 gets fixed
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                category=DeprecationWarning,
                message="elementwise comparison failed",
            )
            values = np.unique(data_by_id.data[:, feature_idx])
        values = values[np.isfinite(values)]
        assert set(values) <= set(range(len(categories)))

    if compare_default_target:
        # check whether the data by id and data by id target are equal
        data_by_id_default = fetch_openml(data_id=data_id, cache=False, as_frame=False)
        np.testing.assert_allclose(data_by_id.data, data_by_id_default.data)
        if data_by_id.target.dtype == np.float64:
            np.testing.assert_allclose(data_by_id.target, data_by_id_default.target)
        else:
            assert np.array_equal(data_by_id.target, data_by_id_default.target)

    if expect_sparse:
        assert isinstance(data_by_id.data, scipy.sparse.csr_matrix)
    else:
        assert isinstance(data_by_id.data, np.ndarray)
        # np.isnan doesn't work on CSR matrix
        assert np.count_nonzero(np.isnan(data_by_id.data)) == expected_missing

    # test return_X_y option
    fetch_func = partial(
        fetch_openml,
        data_id=data_id,
        cache=False,
        target_column=target_column,
        as_frame=False,
    )
    check_return_X_y(data_by_id, fetch_func)
    return data_by_id


class _MockHTTPResponse:
    def __init__(self, data, is_gzip):
        self.data = data
        self.is_gzip = is_gzip

    def read(self, amt=-1):
        return self.data.read(amt)

    def close(self):
        self.data.close()

    def info(self):
        if self.is_gzip:
            return {"Content-Encoding": "gzip"}
        return {}

    def __iter__(self):
        return iter(self.data)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return False


def _monkey_patch_webbased_functions(context, data_id, gzip_response):
    # monkey patches the urlopen function. Important note: Do NOT use this
    # in combination with a regular cache directory, as the files that are
    # stored as cache should not be mixed up with real openml datasets
    url_prefix_data_description = "https://openml.org/api/v1/json/data/"
    url_prefix_data_features = "https://openml.org/api/v1/json/data/features/"
    url_prefix_download_data = "https://openml.org/data/v1/"
    url_prefix_data_list = "https://openml.org/api/v1/json/data/list/"

    path_suffix = ".gz"
    read_fn = gzip.open

    data_module = OPENML_TEST_DATA_MODULE + "." + f"id_{data_id}"

    def _file_name(url, suffix):
        output = (
            re.sub(r"\W", "-", url[len("https://openml.org/") :]) + suffix + path_suffix
        )
        # Shorten the filenames to have better compatibility with windows 10
        # and filenames > 260 characters
        return (
            output.replace("-json-data-list", "-jdl")
            .replace("-json-data-features", "-jdf")
            .replace("-json-data-qualities", "-jdq")
            .replace("-json-data", "-jd")
            .replace("-data_name", "-dn")
            .replace("-download", "-dl")
            .replace("-limit", "-l")
            .replace("-data_version", "-dv")
            .replace("-status", "-s")
            .replace("-deactivated", "-dact")
            .replace("-active", "-act")
        )

    def _mock_urlopen_shared(url, has_gzip_header, expected_prefix, suffix):
        assert url.startswith(expected_prefix)

        data_file_name = _file_name(url, suffix)

        with resources.open_binary(data_module, data_file_name) as f:
            if has_gzip_header and gzip_response:
                fp = BytesIO(f.read())
                return _MockHTTPResponse(fp, True)
            else:
                decompressed_f = read_fn(f, "rb")
                fp = BytesIO(decompressed_f.read())
                return _MockHTTPResponse(fp, False)

    def _mock_urlopen_data_description(url, has_gzip_header):
        return _mock_urlopen_shared(
            url=url,
            has_gzip_header=has_gzip_header,
            expected_prefix=url_prefix_data_description,
            suffix=".json",
        )

    def _mock_urlopen_data_features(url, has_gzip_header):
        return _mock_urlopen_shared(
            url=url,
            has_gzip_header=has_gzip_header,
            expected_prefix=url_prefix_data_features,
            suffix=".json",
        )

    def _mock_urlopen_download_data(url, has_gzip_header):
        return _mock_urlopen_shared(
            url=url,
            has_gzip_header=has_gzip_header,
            expected_prefix=url_prefix_download_data,
            suffix=".arff",
        )

    def _mock_urlopen_data_list(url, has_gzip_header):
        assert url.startswith(url_prefix_data_list)

        data_file_name = _file_name(url, ".json")

        # load the file itself, to simulate a http error
        with resources.open_binary(data_module, data_file_name) as f:
            decompressed_f = read_fn(f, "rb")
            decoded_s = decompressed_f.read().decode("utf-8")
            json_data = json.loads(decoded_s)
        if "error" in json_data:
            raise HTTPError(
                url=None, code=412, msg="Simulated mock error", hdrs=None, fp=None
            )

        with resources.open_binary(data_module, data_file_name) as f:
            if has_gzip_header:
                fp = BytesIO(f.read())
                return _MockHTTPResponse(fp, True)
            else:
                decompressed_f = read_fn(f, "rb")
                fp = BytesIO(decompressed_f.read())
                return _MockHTTPResponse(fp, False)

    def _mock_urlopen(request):
        url = request.get_full_url()
        has_gzip_header = request.get_header("Accept-encoding") == "gzip"
        if url.startswith(url_prefix_data_list):
            return _mock_urlopen_data_list(url, has_gzip_header)
        elif url.startswith(url_prefix_data_features):
            return _mock_urlopen_data_features(url, has_gzip_header)
        elif url.startswith(url_prefix_download_data):
            return _mock_urlopen_download_data(url, has_gzip_header)
        elif url.startswith(url_prefix_data_description):
            return _mock_urlopen_data_description(url, has_gzip_header)
        else:
            raise ValueError("Unknown mocking URL pattern: %s" % url)

    # XXX: Global variable
    if test_offline:
        context.setattr(sklearn.datasets._openml, "urlopen", _mock_urlopen)


###############################################################################
# Test the behaviour of `fetch_openml` depending of the input parameters.

# Known failure of PyPy for OpenML. See the following issue:
# https://github.com/scikit-learn/scikit-learn/issues/18906
@fails_if_pypy
@pytest.mark.parametrize(
    "data_id, n_samples, n_features, n_targets",
    [
        (61, 150, 4, 1),  # iris
        (2, 11, 38, 1),  # anneal
        (561, 209, 7, 1),  # cpu
        (40589, 13, 72, 6),  # emotions
        (1119, 10, 14, 1),  # adult-census
        (40945, 1309, 13, 1),  # titanic
    ],
)
@pytest.mark.parametrize("parser", ["liac-arff", "pandas"])
@pytest.mark.parametrize("infer_casting", [True, False])
@pytest.mark.parametrize("gzip_response", [True, False])
def test_fetch_openml_as_frame_true(
    monkeypatch,
    data_id,
    n_samples,
    n_features,
    n_targets,
    parser,
    infer_casting,
    gzip_response,
):
    """Check the behaviour of `fetch_openml` with `as_frame=True`.

    We explicitely check for the type of the returned container and the shape
    of these containers.
    """
    pd = pytest.importorskip("pandas")

    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response=gzip_response)
    bunch = fetch_openml(
        data_id=data_id,
        as_frame=True,
        cache=False,
        parser=parser,
        infer_casting=infer_casting,
    )

    assert isinstance(bunch, Bunch)

    assert isinstance(bunch.frame, pd.DataFrame)
    assert bunch.frame.shape == (n_samples, n_features + n_targets)

    assert isinstance(bunch.data, pd.DataFrame)
    assert bunch.data.shape == (n_samples, n_features)

    if n_targets == 1:
        assert isinstance(bunch.target, pd.Series)
        assert bunch.target.shape == (n_samples,)
    else:
        assert isinstance(bunch.target, pd.DataFrame)
        assert bunch.target.shape == (n_samples, n_targets)

    assert bunch.categories is None


# Known failure of PyPy for OpenML. See the following issue:
# https://github.com/scikit-learn/scikit-learn/issues/18906
@fails_if_pypy
@pytest.mark.parametrize(
    "data_id, n_samples, n_features, n_targets",
    [
        (61, 150, 4, 1),  # iris
        (2, 11, 38, 1),  # anneal
        (561, 209, 7, 1),  # cpu
        (40589, 13, 72, 6),  # emotions
        (1119, 10, 14, 1),  # adult-census
        (40966, 7, 77, 1),  # miceprotein
    ],
)
@pytest.mark.parametrize("parser", ["liac-arff", "pandas"])
@pytest.mark.parametrize("infer_casting", [True, False])
def test_fetch_openml_as_frame_false(
    monkeypatch, data_id, n_samples, n_features, n_targets, parser, infer_casting
):
    """Check the behaviour of `fetch_openml` with `as_frame=False`.

    We explicitely check for the type of the returned container and the shape
    of these containers.
    """
    pytest.importorskip("pandas")

    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response=True)
    bunch = fetch_openml(
        data_id=data_id,
        as_frame=False,
        cache=False,
        parser=parser,
        infer_casting=infer_casting,
    )

    assert isinstance(bunch, Bunch)

    assert bunch.frame is None

    assert isinstance(bunch.data, np.ndarray)
    assert bunch.data.shape == (n_samples, n_features)

    assert isinstance(bunch.target, np.ndarray)
    if n_targets == 1:
        assert bunch.target.shape == (n_samples,)
    else:
        assert bunch.target.shape == (n_samples, n_targets)

    assert isinstance(bunch.categories, dict)


# Known failure of PyPy for OpenML. See the following issue:
# https://github.com/scikit-learn/scikit-learn/issues/18906
@fails_if_pypy
@pytest.mark.parametrize("parser", ["liac-arff", "pandas"])
@pytest.mark.parametrize("infer_casting", [True, False])
def test_fetch_openml_equivalence_array_dataframe(monkeypatch, parser, infer_casting):
    """Check the equivalence of the dataset when using `as_frame=False` and
    `as_dataframe=True`.
    """
    pytest.importorskip("pandas")

    data_id = 61
    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response=True)
    bunch_as_frame_true = fetch_openml(
        data_id=data_id,
        as_frame=True,
        cache=False,
        parser=parser,
        infer_casting=infer_casting,
    )

    bunch_as_frame_false = fetch_openml(
        data_id=data_id,
        as_frame=True,
        cache=False,
        parser=parser,
        infer_casting=infer_casting,
    )

    assert_allclose(bunch_as_frame_false.data, bunch_as_frame_true.data)
    assert_array_equal(bunch_as_frame_false.target, bunch_as_frame_true.target)


# Known failure of PyPy for OpenML. See the following issue:
# https://github.com/scikit-learn/scikit-learn/issues/18906
@fails_if_pypy
@pytest.mark.parametrize("parser", ["liac-arff", "pandas"])
@pytest.mark.parametrize("infer_casting", [True, False])
def test_fetch_openml_iris_pandas(monkeypatch, parser, infer_casting):
    """Check fetching on a numerical only dataset with string labels."""
    pd = pytest.importorskip("pandas")
    CategoricalDtype = pd.api.types.CategoricalDtype
    data_id = 61
    data_shape = (150, 4)
    target_shape = (150,)
    frame_shape = (150, 5)

    target_dtype = CategoricalDtype(
        ["Iris-setosa", "Iris-versicolor", "Iris-virginica"]
    )
    data_dtypes = [np.float64] * 4
    data_names = ["sepallength", "sepalwidth", "petallength", "petalwidth"]
    target_name = "class"

    _monkey_patch_webbased_functions(monkeypatch, data_id, True)

    bunch = fetch_openml(
        data_id=data_id,
        as_frame=True,
        cache=False,
        parser=parser,
        infer_casting=infer_casting,
    )
    data = bunch.data
    target = bunch.target
    frame = bunch.frame

    assert isinstance(data, pd.DataFrame)
    assert np.all(data.dtypes == data_dtypes)
    assert data.shape == data_shape
    assert np.all(data.columns == data_names)
    assert np.all(bunch.feature_names == data_names)
    assert bunch.target_names == [target_name]

    assert isinstance(target, pd.Series)
    assert target.dtype == target_dtype
    assert target.shape == target_shape
    assert target.name == target_name
    assert target.index.is_unique

    assert isinstance(frame, pd.DataFrame)
    assert frame.shape == frame_shape
    assert np.all(frame.dtypes == data_dtypes + [target_dtype])
    assert frame.index.is_unique


# Known failure of PyPy for OpenML. See the following issue:
# https://github.com/scikit-learn/scikit-learn/issues/18906
@fails_if_pypy
@pytest.mark.parametrize("parser", ["liac-arff", "pandas"])
@pytest.mark.parametrize("infer_casting", [True, False])
@pytest.mark.parametrize("target_column", ["petalwidth", ["petalwidth", "petallength"]])
def test_fetch_openml_forcing_targets(
    monkeypatch, parser, infer_casting, target_column
):
    """Check that we can force the target to not be the default target."""
    pd = pytest.importorskip("pandas")

    data_id = 61
    _monkey_patch_webbased_functions(monkeypatch, data_id, True)
    bunch_forcing_target = fetch_openml(
        data_id=data_id,
        as_frame=True,
        cache=False,
        target_column=target_column,
        parser=parser,
        infer_casting=infer_casting,
    )
    bunch_default = fetch_openml(
        data_id=data_id,
        as_frame=True,
        cache=False,
        parser=parser,
        infer_casting=infer_casting,
    )

    pd.testing.assert_frame_equal(bunch_forcing_target.frame, bunch_default.frame)
    if isinstance(target_column, list):
        pd.testing.assert_index_equal(
            bunch_forcing_target.target.columns, pd.Index(target_column)
        )
        assert bunch_forcing_target.data.shape == (150, 3)
    else:
        assert bunch_forcing_target.target.name == target_column
        assert bunch_forcing_target.data.shape == (150, 4)


# Known failure of PyPy for OpenML. See the following issue:
# https://github.com/scikit-learn/scikit-learn/issues/18906
@fails_if_pypy
@pytest.mark.parametrize("data_id", [61, 2, 561, 40589, 1119])
@pytest.mark.parametrize("parser", ["liac-arff", "pandas"])
@pytest.mark.parametrize("infer_casting", [True, False])
def test_fetch_openml_equivalence_frame_return_X_y(
    monkeypatch, data_id, parser, infer_casting
):
    """Check the behaviour of `return_X_y=True` when `as_frame=True`."""
    pd = pytest.importorskip("pandas")

    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response=True)
    bunch = fetch_openml(
        data_id=data_id,
        as_frame=True,
        cache=False,
        return_X_y=False,
        parser=parser,
        infer_casting=infer_casting,
    )
    X, y = fetch_openml(
        data_id=data_id,
        as_frame=True,
        cache=False,
        return_X_y=True,
        parser=parser,
        infer_casting=infer_casting,
    )

    pd.testing.assert_frame_equal(bunch.data, X)
    if isinstance(y, pd.Series):
        pd.testing.assert_series_equal(bunch.target, y)
    else:
        pd.testing.assert_frame_equal(bunch.target, y)


# Known failure of PyPy for OpenML. See the following issue:
# https://github.com/scikit-learn/scikit-learn/issues/18906
@fails_if_pypy
@pytest.mark.parametrize("data_id", [61, 561, 40589, 1119])
@pytest.mark.parametrize("parser", ["liac-arff", "pandas"])
@pytest.mark.parametrize("infer_casting", [True, False])
def test_fetch_openml_equivalence_array_return_X_y(
    monkeypatch, data_id, parser, infer_casting
):
    """Check the behaviour of `return_X_y=True` when `as_frame=False`."""
    pytest.importorskip("pandas")

    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response=True)
    bunch = fetch_openml(
        data_id=data_id,
        as_frame=False,
        cache=False,
        return_X_y=False,
        parser=parser,
        infer_casting=infer_casting,
    )
    X, y = fetch_openml(
        data_id=data_id,
        as_frame=False,
        cache=False,
        return_X_y=True,
        parser=parser,
        infer_casting=infer_casting,
    )

    assert_array_equal(bunch.data, X)
    assert_array_equal(bunch.target, y)


###############################################################################
# Test the ARFF parsing on several dataset to check if detect the correct
# types (categories, intgers, floats).


@pytest.fixture(scope="module")
def datasets_column_names():
    """Returns the columns names for each dataset."""
    return {
        61: ["sepallength", "sepalwidth", "petallength", "petalwidth", "class"],
        2: [
            "family",
            "product-type",
            "steel",
            "carbon",
            "hardness",
            "temper_rolling",
            "condition",
            "formability",
            "strength",
            "non-ageing",
            "surface-finish",
            "surface-quality",
            "enamelability",
            "bc",
            "bf",
            "bt",
            "bw%2Fme",
            "bl",
            "m",
            "chrom",
            "phos",
            "cbond",
            "marvi",
            "exptl",
            "ferro",
            "corr",
            "blue%2Fbright%2Fvarn%2Fclean",
            "lustre",
            "jurofm",
            "s",
            "p",
            "shape",
            "thick",
            "width",
            "len",
            "oil",
            "bore",
            "packing",
            "class",
        ],
        561: ["vendor", "MYCT", "MMIN", "MMAX", "CACH", "CHMIN", "CHMAX", "class"],
        40589: [
            "Mean_Acc1298_Mean_Mem40_Centroid",
            "Mean_Acc1298_Mean_Mem40_Rolloff",
            "Mean_Acc1298_Mean_Mem40_Flux",
            "Mean_Acc1298_Mean_Mem40_MFCC_0",
            "Mean_Acc1298_Mean_Mem40_MFCC_1",
            "Mean_Acc1298_Mean_Mem40_MFCC_2",
            "Mean_Acc1298_Mean_Mem40_MFCC_3",
            "Mean_Acc1298_Mean_Mem40_MFCC_4",
            "Mean_Acc1298_Mean_Mem40_MFCC_5",
            "Mean_Acc1298_Mean_Mem40_MFCC_6",
            "Mean_Acc1298_Mean_Mem40_MFCC_7",
            "Mean_Acc1298_Mean_Mem40_MFCC_8",
            "Mean_Acc1298_Mean_Mem40_MFCC_9",
            "Mean_Acc1298_Mean_Mem40_MFCC_10",
            "Mean_Acc1298_Mean_Mem40_MFCC_11",
            "Mean_Acc1298_Mean_Mem40_MFCC_12",
            "Mean_Acc1298_Std_Mem40_Centroid",
            "Mean_Acc1298_Std_Mem40_Rolloff",
            "Mean_Acc1298_Std_Mem40_Flux",
            "Mean_Acc1298_Std_Mem40_MFCC_0",
            "Mean_Acc1298_Std_Mem40_MFCC_1",
            "Mean_Acc1298_Std_Mem40_MFCC_2",
            "Mean_Acc1298_Std_Mem40_MFCC_3",
            "Mean_Acc1298_Std_Mem40_MFCC_4",
            "Mean_Acc1298_Std_Mem40_MFCC_5",
            "Mean_Acc1298_Std_Mem40_MFCC_6",
            "Mean_Acc1298_Std_Mem40_MFCC_7",
            "Mean_Acc1298_Std_Mem40_MFCC_8",
            "Mean_Acc1298_Std_Mem40_MFCC_9",
            "Mean_Acc1298_Std_Mem40_MFCC_10",
            "Mean_Acc1298_Std_Mem40_MFCC_11",
            "Mean_Acc1298_Std_Mem40_MFCC_12",
            "Std_Acc1298_Mean_Mem40_Centroid",
            "Std_Acc1298_Mean_Mem40_Rolloff",
            "Std_Acc1298_Mean_Mem40_Flux",
            "Std_Acc1298_Mean_Mem40_MFCC_0",
            "Std_Acc1298_Mean_Mem40_MFCC_1",
            "Std_Acc1298_Mean_Mem40_MFCC_2",
            "Std_Acc1298_Mean_Mem40_MFCC_3",
            "Std_Acc1298_Mean_Mem40_MFCC_4",
            "Std_Acc1298_Mean_Mem40_MFCC_5",
            "Std_Acc1298_Mean_Mem40_MFCC_6",
            "Std_Acc1298_Mean_Mem40_MFCC_7",
            "Std_Acc1298_Mean_Mem40_MFCC_8",
            "Std_Acc1298_Mean_Mem40_MFCC_9",
            "Std_Acc1298_Mean_Mem40_MFCC_10",
            "Std_Acc1298_Mean_Mem40_MFCC_11",
            "Std_Acc1298_Mean_Mem40_MFCC_12",
            "Std_Acc1298_Std_Mem40_Centroid",
            "Std_Acc1298_Std_Mem40_Rolloff",
            "Std_Acc1298_Std_Mem40_Flux",
            "Std_Acc1298_Std_Mem40_MFCC_0",
            "Std_Acc1298_Std_Mem40_MFCC_1",
            "Std_Acc1298_Std_Mem40_MFCC_2",
            "Std_Acc1298_Std_Mem40_MFCC_3",
            "Std_Acc1298_Std_Mem40_MFCC_4",
            "Std_Acc1298_Std_Mem40_MFCC_5",
            "Std_Acc1298_Std_Mem40_MFCC_6",
            "Std_Acc1298_Std_Mem40_MFCC_7",
            "Std_Acc1298_Std_Mem40_MFCC_8",
            "Std_Acc1298_Std_Mem40_MFCC_9",
            "Std_Acc1298_Std_Mem40_MFCC_10",
            "Std_Acc1298_Std_Mem40_MFCC_11",
            "Std_Acc1298_Std_Mem40_MFCC_12",
            "BH_LowPeakAmp",
            "BH_LowPeakBPM",
            "BH_HighPeakAmp",
            "BH_HighPeakBPM",
            "BH_HighLowRatio",
            "BHSUM1",
            "BHSUM2",
            "BHSUM3",
            "amazed.suprised",
            "happy.pleased",
            "relaxing.calm",
            "quiet.still",
            "sad.lonely",
            "angry.aggresive",
        ],
        1119: [
            "age",
            "workclass",
            "fnlwgt:",
            "education:",
            "education-num:",
            "marital-status:",
            "occupation:",
            "relationship:",
            "race:",
            "sex:",
            "capital-gain:",
            "capital-loss:",
            "hours-per-week:",
            "native-country:",
            "class",
        ],
        40966: [
            "DYRK1A_N",
            "ITSN1_N",
            "BDNF_N",
            "NR1_N",
            "NR2A_N",
            "pAKT_N",
            "pBRAF_N",
            "pCAMKII_N",
            "pCREB_N",
            "pELK_N",
            "pERK_N",
            "pJNK_N",
            "PKCA_N",
            "pMEK_N",
            "pNR1_N",
            "pNR2A_N",
            "pNR2B_N",
            "pPKCAB_N",
            "pRSK_N",
            "AKT_N",
            "BRAF_N",
            "CAMKII_N",
            "CREB_N",
            "ELK_N",
            "ERK_N",
            "GSK3B_N",
            "JNK_N",
            "MEK_N",
            "TRKA_N",
            "RSK_N",
            "APP_N",
            "Bcatenin_N",
            "SOD1_N",
            "MTOR_N",
            "P38_N",
            "pMTOR_N",
            "DSCR1_N",
            "AMPKA_N",
            "NR2B_N",
            "pNUMB_N",
            "RAPTOR_N",
            "TIAM1_N",
            "pP70S6_N",
            "NUMB_N",
            "P70S6_N",
            "pGSK3B_N",
            "pPKCG_N",
            "CDK5_N",
            "S6_N",
            "ADARB1_N",
            "AcetylH3K9_N",
            "RRP1_N",
            "BAX_N",
            "ARC_N",
            "ERBB4_N",
            "nNOS_N",
            "Tau_N",
            "GFAP_N",
            "GluR3_N",
            "GluR4_N",
            "IL1B_N",
            "P3525_N",
            "pCASP9_N",
            "PSD95_N",
            "SNCA_N",
            "Ubiquitin_N",
            "pGSK3B_Tyr216_N",
            "SHH_N",
            "BAD_N",
            "BCL2_N",
            "pS6_N",
            "pCFOS_N",
            "SYP_N",
            "H3AcK18_N",
            "EGR1_N",
            "H3MeK4_N",
            "CaNA_N",
            "class",
        ],
        40945: [
            "pclass",
            "survived",
            "name",
            "sex",
            "age",
            "sibsp",
            "parch",
            "ticket",
            "fare",
            "cabin",
            "embarked",
            "boat",
            "body",
            "home.dest",
        ],
    }


@pytest.fixture(scope="module")
def datasets_missing_values():
    return {
        61: {
            "sepallength": 0,
            "sepalwidth": 0,
            "petallength": 0,
            "petalwidth": 0,
            "class": 0,
        },
        2: {
            "family": 11,
            "product-type": 0,
            "steel": 0,
            "carbon": 0,
            "hardness": 0,
            "temper_rolling": 9,
            "condition": 2,
            "formability": 4,
            "strength": 0,
            "non-ageing": 10,
            "surface-finish": 11,
            "surface-quality": 0,
            "enamelability": 11,
            "bc": 11,
            "bf": 10,
            "bt": 11,
            "bw%2Fme": 8,
            "bl": 9,
            "m": 11,
            "chrom": 11,
            "phos": 11,
            "cbond": 10,
            "marvi": 11,
            "exptl": 11,
            "ferro": 11,
            "corr": 11,
            "blue%2Fbright%2Fvarn%2Fclean": 11,
            "lustre": 8,
            "jurofm": 11,
            "s": 11,
            "p": 11,
            "shape": 0,
            "thick": 0,
            "width": 0,
            "len": 0,
            "oil": 10,
            "bore": 0,
            "packing": 11,
            "class": 0,
        },
        561: {
            "vendor": 0,
            "MYCT": 0,
            "MMIN": 0,
            "MMAX": 0,
            "CACH": 0,
            "CHMIN": 0,
            "CHMAX": 0,
            "class": 0,
        },
        40589: {
            "Mean_Acc1298_Mean_Mem40_Centroid": 0,
            "Mean_Acc1298_Mean_Mem40_Rolloff": 0,
            "Mean_Acc1298_Mean_Mem40_Flux": 0,
            "Mean_Acc1298_Mean_Mem40_MFCC_0": 0,
            "Mean_Acc1298_Mean_Mem40_MFCC_1": 0,
            "Mean_Acc1298_Mean_Mem40_MFCC_2": 0,
            "Mean_Acc1298_Mean_Mem40_MFCC_3": 0,
            "Mean_Acc1298_Mean_Mem40_MFCC_4": 0,
            "Mean_Acc1298_Mean_Mem40_MFCC_5": 0,
            "Mean_Acc1298_Mean_Mem40_MFCC_6": 0,
            "Mean_Acc1298_Mean_Mem40_MFCC_7": 0,
            "Mean_Acc1298_Mean_Mem40_MFCC_8": 0,
            "Mean_Acc1298_Mean_Mem40_MFCC_9": 0,
            "Mean_Acc1298_Mean_Mem40_MFCC_10": 0,
            "Mean_Acc1298_Mean_Mem40_MFCC_11": 0,
            "Mean_Acc1298_Mean_Mem40_MFCC_12": 0,
            "Mean_Acc1298_Std_Mem40_Centroid": 0,
            "Mean_Acc1298_Std_Mem40_Rolloff": 0,
            "Mean_Acc1298_Std_Mem40_Flux": 0,
            "Mean_Acc1298_Std_Mem40_MFCC_0": 0,
            "Mean_Acc1298_Std_Mem40_MFCC_1": 0,
            "Mean_Acc1298_Std_Mem40_MFCC_2": 0,
            "Mean_Acc1298_Std_Mem40_MFCC_3": 0,
            "Mean_Acc1298_Std_Mem40_MFCC_4": 0,
            "Mean_Acc1298_Std_Mem40_MFCC_5": 0,
            "Mean_Acc1298_Std_Mem40_MFCC_6": 0,
            "Mean_Acc1298_Std_Mem40_MFCC_7": 0,
            "Mean_Acc1298_Std_Mem40_MFCC_8": 0,
            "Mean_Acc1298_Std_Mem40_MFCC_9": 0,
            "Mean_Acc1298_Std_Mem40_MFCC_10": 0,
            "Mean_Acc1298_Std_Mem40_MFCC_11": 0,
            "Mean_Acc1298_Std_Mem40_MFCC_12": 0,
            "Std_Acc1298_Mean_Mem40_Centroid": 0,
            "Std_Acc1298_Mean_Mem40_Rolloff": 0,
            "Std_Acc1298_Mean_Mem40_Flux": 0,
            "Std_Acc1298_Mean_Mem40_MFCC_0": 0,
            "Std_Acc1298_Mean_Mem40_MFCC_1": 0,
            "Std_Acc1298_Mean_Mem40_MFCC_2": 0,
            "Std_Acc1298_Mean_Mem40_MFCC_3": 0,
            "Std_Acc1298_Mean_Mem40_MFCC_4": 0,
            "Std_Acc1298_Mean_Mem40_MFCC_5": 0,
            "Std_Acc1298_Mean_Mem40_MFCC_6": 0,
            "Std_Acc1298_Mean_Mem40_MFCC_7": 0,
            "Std_Acc1298_Mean_Mem40_MFCC_8": 0,
            "Std_Acc1298_Mean_Mem40_MFCC_9": 0,
            "Std_Acc1298_Mean_Mem40_MFCC_10": 0,
            "Std_Acc1298_Mean_Mem40_MFCC_11": 0,
            "Std_Acc1298_Mean_Mem40_MFCC_12": 0,
            "Std_Acc1298_Std_Mem40_Centroid": 0,
            "Std_Acc1298_Std_Mem40_Rolloff": 0,
            "Std_Acc1298_Std_Mem40_Flux": 0,
            "Std_Acc1298_Std_Mem40_MFCC_0": 0,
            "Std_Acc1298_Std_Mem40_MFCC_1": 0,
            "Std_Acc1298_Std_Mem40_MFCC_2": 0,
            "Std_Acc1298_Std_Mem40_MFCC_3": 0,
            "Std_Acc1298_Std_Mem40_MFCC_4": 0,
            "Std_Acc1298_Std_Mem40_MFCC_5": 0,
            "Std_Acc1298_Std_Mem40_MFCC_6": 0,
            "Std_Acc1298_Std_Mem40_MFCC_7": 0,
            "Std_Acc1298_Std_Mem40_MFCC_8": 0,
            "Std_Acc1298_Std_Mem40_MFCC_9": 0,
            "Std_Acc1298_Std_Mem40_MFCC_10": 0,
            "Std_Acc1298_Std_Mem40_MFCC_11": 0,
            "Std_Acc1298_Std_Mem40_MFCC_12": 0,
            "BH_LowPeakAmp": 0,
            "BH_LowPeakBPM": 0,
            "BH_HighPeakAmp": 0,
            "BH_HighPeakBPM": 0,
            "BH_HighLowRatio": 0,
            "BHSUM1": 0,
            "BHSUM2": 0,
            "BHSUM3": 0,
            "amazed.suprised": 0,
            "happy.pleased": 0,
            "relaxing.calm": 0,
            "quiet.still": 0,
            "sad.lonely": 0,
            "angry.aggresive": 0,
        },
        1119: {
            "age": 0,
            "workclass": 0,
            "fnlwgt:": 0,
            "education:": 0,
            "education-num:": 0,
            "marital-status:": 0,
            "occupation:": 0,
            "relationship:": 0,
            "race:": 0,
            "sex:": 0,
            "capital-gain:": 0,
            "capital-loss:": 0,
            "hours-per-week:": 0,
            "native-country:": 0,
            "class": 0,
        },
        40966: {
            "DYRK1A_N": 0,
            "ITSN1_N": 0,
            "BDNF_N": 0,
            "NR1_N": 0,
            "NR2A_N": 0,
            "pAKT_N": 0,
            "pBRAF_N": 0,
            "pCAMKII_N": 0,
            "pCREB_N": 0,
            "pELK_N": 0,
            "pERK_N": 0,
            "pJNK_N": 0,
            "PKCA_N": 0,
            "pMEK_N": 0,
            "pNR1_N": 0,
            "pNR2A_N": 0,
            "pNR2B_N": 0,
            "pPKCAB_N": 0,
            "pRSK_N": 0,
            "AKT_N": 0,
            "BRAF_N": 0,
            "CAMKII_N": 0,
            "CREB_N": 0,
            "ELK_N": 0,
            "ERK_N": 0,
            "GSK3B_N": 0,
            "JNK_N": 0,
            "MEK_N": 0,
            "TRKA_N": 0,
            "RSK_N": 0,
            "APP_N": 0,
            "Bcatenin_N": 0,
            "SOD1_N": 0,
            "MTOR_N": 0,
            "P38_N": 0,
            "pMTOR_N": 0,
            "DSCR1_N": 0,
            "AMPKA_N": 0,
            "NR2B_N": 0,
            "pNUMB_N": 0,
            "RAPTOR_N": 0,
            "TIAM1_N": 0,
            "pP70S6_N": 0,
            "NUMB_N": 0,
            "P70S6_N": 0,
            "pGSK3B_N": 0,
            "pPKCG_N": 0,
            "CDK5_N": 0,
            "S6_N": 0,
            "ADARB1_N": 0,
            "AcetylH3K9_N": 0,
            "RRP1_N": 0,
            "BAX_N": 0,
            "ARC_N": 0,
            "ERBB4_N": 0,
            "nNOS_N": 0,
            "Tau_N": 0,
            "GFAP_N": 0,
            "GluR3_N": 0,
            "GluR4_N": 0,
            "IL1B_N": 0,
            "P3525_N": 0,
            "pCASP9_N": 0,
            "PSD95_N": 0,
            "SNCA_N": 0,
            "Ubiquitin_N": 0,
            "pGSK3B_Tyr216_N": 0,
            "SHH_N": 0,
            "BAD_N": 0,
            "BCL2_N": 7,
            "pS6_N": 0,
            "pCFOS_N": 0,
            "SYP_N": 0,
            "H3AcK18_N": 0,
            "EGR1_N": 0,
            "H3MeK4_N": 0,
            "CaNA_N": 0,
            "class": 0,
        },
        40945: {
            "pclass": 0,
            "survived": 0,
            "name": 0,
            "sex": 0,
            "age": 263,
            "sibsp": 0,
            "parch": 0,
            "ticket": 0,
            "fare": 1,
            "cabin": 1014,
            "embarked": 2,
            "boat": 823,
            "body": 1188,
            "home.dest": 564,
        },
    }


# Known failure of PyPy for OpenML. See the following issue:
# https://github.com/scikit-learn/scikit-learn/issues/18906
@fails_if_pypy
@pytest.mark.parametrize(
    "data_id, parser, infer_casting, expected_n_categories, expected_n_floats,"
    " expected_n_ints",
    [
        # iris dataset
        (61, "liac-arff", False, 1, 4, 0),
        (61, "liac-arff", True, 1, 4, 0),
        (61, "pandas", False, 1, 4, 0),
        (61, "pandas", True, 1, 4, 0),
        # anneal dataset
        (2, "liac-arff", False, 33, 6, 0),
        (2, "liac-arff", True, 33, 2, 4),
        (2, "pandas", False, 33, 2, 4),
        (2, "pandas", True, 33, 2, 4),
        # cpu dataset
        (561, "liac-arff", False, 1, 7, 0),
        (561, "liac-arff", True, 1, 0, 7),
        (561, "pandas", False, 1, 0, 7),
        (561, "pandas", True, 1, 0, 7),
        # emotions dataset
        (40589, "liac-arff", False, 6, 72, 0),
        (40589, "liac-arff", True, 6, 69, 3),
        (40589, "pandas", False, 6, 69, 3),
        (40589, "pandas", True, 6, 69, 3),
        # adult-census dataset
        (1119, "liac-arff", False, 9, 6, 0),
        (1119, "liac-arff", True, 9, 0, 6),
        (1119, "pandas", False, 9, 0, 6),
        (1119, "pandas", True, 9, 0, 6),
        # miceprotein
        # 1 column has only missing values with object dtype
        (40966, "liac-arff", False, 1, 76, 0),
        # with casting it will be transformed to either float or Int64
        (40966, "liac-arff", True, 1, 76, 1),
        (40966, "pandas", False, 1, 77, 0),
        (40966, "pandas", True, 1, 76, 1),
        # titanic
        (40945, "liac-arff", False, 3, 5, 0),
        (40945, "liac-arff", True, 3, 2, 4),
        (40945, "pandas", False, 3, 3, 3),
        (40945, "pandas", True, 3, 2, 4),
    ],
)
@pytest.mark.parametrize("gzip_response", [True, False])
def test_fetch_openml_types_inference(
    monkeypatch,
    data_id,
    parser,
    infer_casting,
    expected_n_categories,
    expected_n_floats,
    expected_n_ints,
    gzip_response,
    datasets_column_names,
    datasets_missing_values,
):
    """Check that `fetch_openml` infer the right number of categories, integers, and
    floats."""
    pd = pytest.importorskip("pandas")
    CategoricalDtype = pd.api.types.CategoricalDtype

    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response=gzip_response)

    bunch = fetch_openml(
        data_id=data_id,
        as_frame=True,
        cache=False,
        parser=parser,
        infer_casting=infer_casting,
    )
    frame = bunch.frame

    n_categories = len(
        [dtype for dtype in frame.dtypes if isinstance(dtype, CategoricalDtype)]
    )
    n_floats = len([dtype for dtype in frame.dtypes if dtype.kind == "f"])
    n_ints = len([dtype for dtype in frame.dtypes if dtype.kind == "i"])

    assert n_categories == expected_n_categories
    assert n_floats == expected_n_floats
    assert n_ints == expected_n_ints

    assert frame.columns.tolist() == datasets_column_names[data_id]
    assert frame.isna().sum().to_dict() == datasets_missing_values[data_id]


###############################################################################
# Test some more specific behaviour


@pytest.mark.filterwarnings("ignore:Version 1 of dataset Australian is inactive")
@pytest.mark.parametrize("parser", ["liac-arff", "pandas"])
@pytest.mark.parametrize("infer_casting", [True, False])
def test_fetch_openml_australian_pandas_error_sparse(
    monkeypatch, parser, infer_casting
):
    """Check that we raise an error if a dataset is sparse and we try to request a
    dataframe.
    """
    data_id = 292

    _monkey_patch_webbased_functions(monkeypatch, data_id, True)

    msg = "Cannot return dataframe with sparse data. Switch to `as_frame=False`."
    with pytest.raises(ValueError, match=msg):
        fetch_openml(
            data_id=data_id,
            as_frame=True,
            cache=False,
            parser=parser,
            infer_casting=infer_casting,
        )


@pytest.mark.filterwarnings("ignore:Version 1 of dataset Australian is inactive")
def test_fetch_openml_pandas_parser_error_on_sparse(monkeypatch):
    """Check that we raise an error using the pandas parser on a sparse dataset."""
    data_id = 292

    _monkey_patch_webbased_functions(monkeypatch, data_id, True)

    msg = (
        "Cannot use `parser='pandas'` with sparse ARFF file. Switch to "
        "`parser='liac-arff'`."
    )
    with pytest.raises(ValueError, match=msg):
        fetch_openml(
            data_id=data_id,
            as_frame=False,
            cache=False,
            parser="pandas",
        )


# Known failure of PyPy for OpenML. See the following issue:
# https://github.com/scikit-learn/scikit-learn/issues/18906
@fails_if_pypy
@pytest.mark.filterwarnings("ignore:Version 1 of dataset Australian is inactive")
@pytest.mark.parametrize(
    "data_id, data_type",
    [
        (61, "dataframe"),  # iris dataset version 1
        (292, "sparse"),  # Australian dataset version 1
    ],
)
def test_fetch_openml_auto_mode(monkeypatch, data_id, data_type):
    """Check the auto mode of `fetch_openml`."""
    pd = pytest.importorskip("pandas")

    _monkey_patch_webbased_functions(monkeypatch, data_id, True)
    data = fetch_openml(data_id=data_id, as_frame="auto", parser="auto", cache=False)
    klass = pd.DataFrame if data_type == "dataframe" else scipy.sparse.csr_matrix
    assert isinstance(data.data, klass)


# Known failure of PyPy for OpenML. See the following issue:
# https://github.com/scikit-learn/scikit-learn/issues/18906
@fails_if_pypy
@pytest.mark.parametrize("infer_casting", [True, False])
def test_convert_arff_data_dataframe_warning_low_memory_pandas(
    monkeypatch, infer_casting
):
    """Check that we raise a warning regarding the working memory when using
    LIAC-ARFF parser."""
    pytest.importorskip("pandas")

    data_id = 1119
    _monkey_patch_webbased_functions(monkeypatch, data_id, True)

    msg = "Could not adhere to working_memory config."
    with pytest.warns(UserWarning, match=msg):
        with config_context(working_memory=1e-6):
            fetch_openml(
                data_id=data_id,
                as_frame=True,
                cache=False,
                parser="liac-arff",
                infer_casting=infer_casting,
            )


@pytest.mark.parametrize("gzip_response", [True, False])
@pytest.mark.parametrize("parser", ["pandas", "liac-arff"])
def test_fetch_openml_iris_warn_multiple_version(monkeypatch, gzip_response, parser):
    """Check that a warning is raised when multiple versions exist and no version is
    requested."""
    data_id = 61
    data_name = "iris"

    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response)

    msg = (
        "Multiple active versions of the dataset matching the name"
        " iris exist. Versions may be fundamentally different, "
        "returning version 1."
    )
    with pytest.warns(UserWarning, match=msg):
        fetch_openml(
            name=data_name,
            as_frame=False,
            cache=False,
            parser=parser,
        )


@pytest.mark.parametrize("gzip_response", [True, False])
def test_fetch_openml_no_target(monkeypatch, gzip_response):
    """Check that we can get a dataset without target."""
    data_id = 61
    target_column = None
    expected_observations = 150
    expected_features = 5

    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response)
    data = fetch_openml(
        data_id=data_id, target_column=target_column, cache=False, as_frame=False
    )
    assert data.data.shape == (expected_observations, expected_features)
    assert data.target is None


@pytest.mark.parametrize("gzip_response", [True, False])
@pytest.mark.parametrize("parser", ["liac-arff", "pandas"])
@pytest.mark.parametrize("infer_casting", [True, False])
def test_missing_values_pandas(monkeypatch, gzip_response, parser, infer_casting):
    """check that missing values in categories are compatible with pandas
    categorical"""
    pytest.importorskip("pandas")

    data_id = 42585
    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response=gzip_response)
    penguins = fetch_openml(
        data_id=data_id,
        cache=False,
        as_frame=True,
        parser=parser,
        infer_casting=infer_casting,
    )

    cat_dtype = penguins.data.dtypes["sex"]
    # there are nans in the categorical
    assert penguins.data["sex"].isna().any()
    assert_array_equal(cat_dtype.categories, ["FEMALE", "MALE", "_"])


@pytest.mark.parametrize("gzip_response", [True, False])
@pytest.mark.parametrize(
    "dataset_params",
    [
        {"data_id": 40675},
        {"data_id": None, "name": "glass2", "version": 1},
    ],
)
def test_fetch_openml_inactive(monkeypatch, gzip_response, dataset_params):
    """Check that we raise a warning when the dataset is inactive."""
    data_id = 40675
    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response)
    msg = "Version 1 of dataset glass2 is inactive,"
    with pytest.warns(UserWarning, match=msg):
        glass2 = fetch_openml(cache=False, as_frame=False, **dataset_params)
    assert glass2.data.shape == (163, 9)
    assert glass2.details["id"] == "40675"


@pytest.mark.parametrize("gzip_response", [True, False])
@pytest.mark.parametrize(
    "data_id, params, err_type, err_msg",
    [
        (40675, {"name": "glass2"}, ValueError, "No active dataset glass2 found"),
        (
            61,
            {"data_id": 61, "target_column": ["sepalwidth", "class"]},
            ValueError,
            "Can only handle homogeneous multi-target datasets",
        ),
        (
            40945,
            {"data_id": 40945, "as_frame": False},
            ValueError,
            "STRING attributes are not supported for array representation. Try"
            " as_frame=True",
        ),
        (
            2,
            {"data_id": 2, "target_column": "family", "as_frame": True},
            ValueError,
            "Target column 'family'",
        ),
        (
            2,
            {"data_id": 2, "target_column": "family", "as_frame": False},
            ValueError,
            "Target column 'family'",
        ),
        (
            61,
            {"data_id": 61, "target_column": "undefined"},
            KeyError,
            "Could not find target_column='undefined'",
        ),
        (
            61,
            {"data_id": 61, "target_column": ["undefined", "class"]},
            KeyError,
            "Could not find target_column='undefined'",
        ),
    ],
)
def test_fetch_openml_error(
    monkeypatch, gzip_response, data_id, params, err_type, err_msg
):
    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response)
    with pytest.raises(err_type, match=err_msg):
        fetch_openml(cache=False, **params)


@pytest.mark.parametrize(
    "params, err_type, err_msg",
    [
        (
            {"data_id": -1, "name": None, "version": "version"},
            ValueError,
            "Dataset data_id=-1 and version=version passed, but you can only",
        ),
        (
            {"data_id": -1, "name": "nAmE"},
            ValueError,
            "Dataset data_id=-1 and name=name passed, but you can only",
        ),
        (
            {"data_id": -1, "name": "nAmE", "version": "version"},
            ValueError,
            "Dataset data_id=-1 and name=name passed, but you can only",
        ),
        (
            {},
            ValueError,
            "Neither name nor data_id are provided. Please provide name or data_id.",
        ),
    ],
)
def test_fetch_openml_raises_illegal_argument(params, err_type, err_msg):
    with pytest.raises(err_type, match=err_msg):
        fetch_openml(**params)


@pytest.mark.parametrize("gzip_response", [True, False])
def test_warn_ignore_attribute(monkeypatch, gzip_response):
    data_id = 40966
    expected_row_id_msg = "target_column='{}' has flag is_row_identifier."
    expected_ignore_msg = "target_column='{}' has flag is_ignore."
    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response)
    # single column test
    target_col = "MouseID"
    msg = expected_row_id_msg.format(target_col)
    with pytest.warns(UserWarning, match=msg):
        fetch_openml(
            data_id=data_id, target_column=target_col, cache=False, as_frame=False
        )
    target_col = "Genotype"
    msg = expected_ignore_msg.format(target_col)
    with pytest.warns(UserWarning, match=msg):
        fetch_openml(
            data_id=data_id, target_column=target_col, cache=False, as_frame=False
        )
    # multi column test
    target_col = "MouseID"
    msg = expected_row_id_msg.format(target_col)
    with pytest.warns(UserWarning, match=msg):
        fetch_openml(
            data_id=data_id,
            target_column=[target_col, "class"],
            cache=False,
            as_frame=False,
        )
    target_col = "Genotype"
    msg = expected_ignore_msg.format(target_col)
    with pytest.warns(UserWarning, match=msg):
        fetch_openml(
            data_id=data_id,
            target_column=[target_col, "class"],
            cache=False,
            as_frame=False,
        )


@pytest.mark.parametrize("gzip_response", [True, False])
def test_dataset_with_openml_error(monkeypatch, gzip_response):
    data_id = 1
    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response)
    msg = "OpenML registered a problem with the dataset. It might be unusable. Error:"
    with pytest.warns(UserWarning, match=msg):
        fetch_openml(data_id=data_id, cache=False, as_frame=False)


@pytest.mark.parametrize("gzip_response", [True, False])
def test_dataset_with_openml_warning(monkeypatch, gzip_response):
    data_id = 3
    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response)
    msg = "OpenML raised a warning on the dataset. It might be unusable. Warning:"
    with pytest.warns(UserWarning, match=msg):
        fetch_openml(data_id=data_id, cache=False, as_frame=False)


###############################################################################
# Test cache, retry mechanisms, checksum, etc.


@pytest.mark.parametrize("gzip_response", [True, False])
def test_open_openml_url_cache(monkeypatch, gzip_response, tmpdir):
    data_id = 61

    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response)
    openml_path = sklearn.datasets._openml._DATA_FILE.format(data_id)
    cache_directory = str(tmpdir.mkdir("scikit_learn_data"))
    # first fill the cache
    response1 = _open_openml_url(openml_path, cache_directory)
    # assert file exists
    location = _get_local_path(openml_path, cache_directory)
    assert os.path.isfile(location)
    # redownload, to utilize cache
    response2 = _open_openml_url(openml_path, cache_directory)
    assert response1.read() == response2.read()


@pytest.mark.parametrize("write_to_disk", [True, False])
def test_open_openml_url_unlinks_local_path(monkeypatch, tmpdir, write_to_disk):
    data_id = 61
    openml_path = sklearn.datasets._openml._DATA_FILE.format(data_id)
    cache_directory = str(tmpdir.mkdir("scikit_learn_data"))
    location = _get_local_path(openml_path, cache_directory)

    def _mock_urlopen(request):
        if write_to_disk:
            with open(location, "w") as f:
                f.write("")
        raise ValueError("Invalid request")

    monkeypatch.setattr(sklearn.datasets._openml, "urlopen", _mock_urlopen)

    with pytest.raises(ValueError, match="Invalid request"):
        _open_openml_url(openml_path, cache_directory)

    assert not os.path.exists(location)


def test_retry_with_clean_cache(tmpdir):
    data_id = 61
    openml_path = sklearn.datasets._openml._DATA_FILE.format(data_id)
    cache_directory = str(tmpdir.mkdir("scikit_learn_data"))
    location = _get_local_path(openml_path, cache_directory)
    os.makedirs(os.path.dirname(location))

    with open(location, "w") as f:
        f.write("")

    @_retry_with_clean_cache(openml_path, cache_directory)
    def _load_data():
        # The first call will raise an error since location exists
        if os.path.exists(location):
            raise Exception("File exist!")
        return 1

    warn_msg = "Invalid cache, redownloading file"
    with pytest.warns(RuntimeWarning, match=warn_msg):
        result = _load_data()
    assert result == 1


def test_retry_with_clean_cache_http_error(tmpdir):
    data_id = 61
    openml_path = sklearn.datasets._openml._DATA_FILE.format(data_id)
    cache_directory = str(tmpdir.mkdir("scikit_learn_data"))

    @_retry_with_clean_cache(openml_path, cache_directory)
    def _load_data():
        raise HTTPError(
            url=None, code=412, msg="Simulated mock error", hdrs=None, fp=None
        )

    error_msg = "Simulated mock error"
    with pytest.raises(HTTPError, match=error_msg):
        _load_data()


@pytest.mark.parametrize("gzip_response", [True, False])
def test_fetch_openml_cache(monkeypatch, gzip_response, tmpdir):
    def _mock_urlopen_raise(request):
        raise ValueError(
            "This mechanism intends to test correct cache"
            "handling. As such, urlopen should never be "
            "accessed. URL: %s"
            % request.get_full_url()
        )

    data_id = 61
    cache_directory = str(tmpdir.mkdir("scikit_learn_data"))
    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response)
    X_fetched, y_fetched = fetch_openml(
        data_id=data_id,
        cache=True,
        data_home=cache_directory,
        return_X_y=True,
        as_frame=False,
    )

    monkeypatch.setattr(sklearn.datasets._openml, "urlopen", _mock_urlopen_raise)

    X_cached, y_cached = fetch_openml(
        data_id=data_id,
        cache=True,
        data_home=cache_directory,
        return_X_y=True,
        as_frame=False,
    )
    np.testing.assert_array_equal(X_fetched, X_cached)
    np.testing.assert_array_equal(y_fetched, y_cached)


# Known failure of PyPy for OpenML. See the following issue:
# https://github.com/scikit-learn/scikit-learn/issues/18906
@fails_if_pypy
@pytest.mark.parametrize("as_frame", [True, False])
@pytest.mark.parametrize("parser", ["pandas", "liac-arff"])
def test_fetch_openml_verify_checksum(monkeypatch, as_frame, cache, tmpdir, parser):
    """Check that the checksum is working as expected."""
    if as_frame:
        pytest.importorskip("pandas")

    data_id = 2
    _monkey_patch_webbased_functions(monkeypatch, data_id, True)

    # create a temporary modified arff file
    original_data_module = OPENML_TEST_DATA_MODULE + "." + f"id_{data_id}"
    original_data_file_name = "data-v1-dl-1666876.arff.gz"
    corrupt_copy_path = tmpdir / "test_invalid_checksum.arff"
    with resources.open_binary(
        original_data_module, original_data_file_name
    ) as orig_file:
        orig_gzip = gzip.open(orig_file, "rb")
        data = bytearray(orig_gzip.read())
        data[len(data) - 1] = 37

    with gzip.GzipFile(corrupt_copy_path, "wb") as modified_gzip:
        modified_gzip.write(data)

    # Requests are already mocked by monkey_patch_webbased_functions.
    # We want to re-use that mock for all requests except file download,
    # hence creating a thin mock over the original mock
    mocked_openml_url = sklearn.datasets._openml.urlopen

    def swap_file_mock(request):
        url = request.get_full_url()
        if url.endswith("data/v1/download/1666876"):
            return _MockHTTPResponse(open(corrupt_copy_path, "rb"), is_gzip=True)
        else:
            return mocked_openml_url(request)

    monkeypatch.setattr(sklearn.datasets._openml, "urlopen", swap_file_mock)

    # validate failed checksum
    with pytest.raises(ValueError) as exc:
        sklearn.datasets.fetch_openml(
            data_id=data_id, cache=False, as_frame=as_frame, parser=parser
        )
    # exception message should have file-path
    assert exc.match("1666876")


###############################################################################
# Non-regressiont tests


@pytest.mark.parametrize("gzip_response", [True, False])
@pytest.mark.parametrize("parser", ["pandas", "liac-arff"])
def test_fetch_openml_with_ignored_feature(monkeypatch, gzip_response, parser):
    """Check that we can load the "zoo" dataset.
    Non-regression test for:
    https://github.com/scikit-learn/scikit-learn/issues/14340
    """
    data_id = 62
    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response)

    dataset = sklearn.datasets.fetch_openml(
        data_id=data_id, cache=False, as_frame=False, parser=parser
    )
    assert dataset is not None
    # The dataset has 17 features, including 1 ignored (animal),
    # so we assert that we don't have the ignored feature in the final Bunch
    assert dataset["data"].shape == (101, 16)
    assert "animal" not in dataset["feature_names"]
