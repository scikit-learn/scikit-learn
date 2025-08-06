import hashlib
import io
import os
import re
import shutil
import tempfile
import warnings
from functools import partial
from importlib import resources
from pathlib import Path
from pickle import dumps, loads
from unittest.mock import Mock
from urllib.error import HTTPError
from urllib.parse import urlparse

import numpy as np
import pytest

from sklearn.datasets import (
    clear_data_home,
    fetch_file,
    get_data_home,
    load_breast_cancer,
    load_diabetes,
    load_digits,
    load_files,
    load_iris,
    load_linnerud,
    load_sample_image,
    load_sample_images,
    load_wine,
)
from sklearn.datasets._base import (
    RemoteFileMetadata,
    _derive_folder_and_filename_from_url,
    _fetch_remote,
    load_csv_data,
    load_gzip_compressed_csv_data,
)
from sklearn.datasets.tests.test_common import check_as_frame
from sklearn.preprocessing import scale
from sklearn.utils import Bunch


class _DummyPath:
    """Minimal class that implements the os.PathLike interface."""

    def __init__(self, path):
        self.path = path

    def __fspath__(self):
        return self.path


def _remove_dir(path):
    if os.path.isdir(path):
        shutil.rmtree(path)


@pytest.fixture(scope="module")
def data_home(tmpdir_factory):
    tmp_file = str(tmpdir_factory.mktemp("scikit_learn_data_home_test"))
    yield tmp_file
    _remove_dir(tmp_file)


@pytest.fixture(scope="module")
def load_files_root(tmpdir_factory):
    tmp_file = str(tmpdir_factory.mktemp("scikit_learn_load_files_test"))
    yield tmp_file
    _remove_dir(tmp_file)


@pytest.fixture
def test_category_dir_1(load_files_root):
    test_category_dir1 = tempfile.mkdtemp(dir=load_files_root)
    sample_file = tempfile.NamedTemporaryFile(dir=test_category_dir1, delete=False)
    sample_file.write(b"Hello World!\n")
    sample_file.close()
    yield str(test_category_dir1)
    _remove_dir(test_category_dir1)


@pytest.fixture
def test_category_dir_2(load_files_root):
    test_category_dir2 = tempfile.mkdtemp(dir=load_files_root)
    yield str(test_category_dir2)
    _remove_dir(test_category_dir2)


@pytest.mark.thread_unsafe
@pytest.mark.parametrize("path_container", [None, Path, _DummyPath])
def test_data_home(path_container, data_home):
    # get_data_home will point to a pre-existing folder
    if path_container is not None:
        data_home = path_container(data_home)
    data_home = get_data_home(data_home=data_home)
    assert data_home == data_home
    assert os.path.exists(data_home)

    # clear_data_home will delete both the content and the folder it-self
    if path_container is not None:
        data_home = path_container(data_home)
    clear_data_home(data_home=data_home)
    assert not os.path.exists(data_home)

    # if the folder is missing it will be created again
    data_home = get_data_home(data_home=data_home)
    assert os.path.exists(data_home)


def test_default_empty_load_files(load_files_root):
    res = load_files(load_files_root)
    assert len(res.filenames) == 0
    assert len(res.target_names) == 0
    assert res.DESCR is None


def test_default_load_files(test_category_dir_1, test_category_dir_2, load_files_root):
    res = load_files(load_files_root)
    assert len(res.filenames) == 1
    assert len(res.target_names) == 2
    assert res.DESCR is None
    assert res.data == [b"Hello World!\n"]


def test_load_files_w_categories_desc_and_encoding(
    test_category_dir_1, test_category_dir_2, load_files_root
):
    category = os.path.abspath(test_category_dir_1).split(os.sep).pop()
    res = load_files(
        load_files_root, description="test", categories=[category], encoding="utf-8"
    )

    assert len(res.filenames) == 1
    assert len(res.target_names) == 1
    assert res.DESCR == "test"
    assert res.data == ["Hello World!\n"]


def test_load_files_wo_load_content(
    test_category_dir_1, test_category_dir_2, load_files_root
):
    res = load_files(load_files_root, load_content=False)
    assert len(res.filenames) == 1
    assert len(res.target_names) == 2
    assert res.DESCR is None
    assert res.get("data") is None


@pytest.mark.parametrize("allowed_extensions", ([".txt"], [".txt", ".json"]))
def test_load_files_allowed_extensions(tmp_path, allowed_extensions):
    """Check the behaviour of `allowed_extension` in `load_files`."""
    d = tmp_path / "sub"
    d.mkdir()
    files = ("file1.txt", "file2.json", "file3.json", "file4.md")
    paths = [d / f for f in files]
    for p in paths:
        p.write_bytes(b"hello")
    res = load_files(tmp_path, allowed_extensions=allowed_extensions)
    assert set([str(p) for p in paths if p.suffix in allowed_extensions]) == set(
        res.filenames
    )


@pytest.mark.parametrize(
    "filename, expected_n_samples, expected_n_features, expected_target_names",
    [
        ("wine_data.csv", 178, 13, ["class_0", "class_1", "class_2"]),
        ("iris.csv", 150, 4, ["setosa", "versicolor", "virginica"]),
        ("breast_cancer.csv", 569, 30, ["malignant", "benign"]),
    ],
)
def test_load_csv_data(
    filename, expected_n_samples, expected_n_features, expected_target_names
):
    actual_data, actual_target, actual_target_names = load_csv_data(filename)
    assert actual_data.shape[0] == expected_n_samples
    assert actual_data.shape[1] == expected_n_features
    assert actual_target.shape[0] == expected_n_samples
    np.testing.assert_array_equal(actual_target_names, expected_target_names)


def test_load_csv_data_with_descr():
    data_file_name = "iris.csv"
    descr_file_name = "iris.rst"

    res_without_descr = load_csv_data(data_file_name=data_file_name)
    res_with_descr = load_csv_data(
        data_file_name=data_file_name, descr_file_name=descr_file_name
    )
    assert len(res_with_descr) == 4
    assert len(res_without_descr) == 3

    np.testing.assert_array_equal(res_with_descr[0], res_without_descr[0])
    np.testing.assert_array_equal(res_with_descr[1], res_without_descr[1])
    np.testing.assert_array_equal(res_with_descr[2], res_without_descr[2])

    assert res_with_descr[-1].startswith(".. _iris_dataset:")


@pytest.mark.parametrize(
    "filename, kwargs, expected_shape",
    [
        ("diabetes_data_raw.csv.gz", {}, [442, 10]),
        ("diabetes_target.csv.gz", {}, [442]),
        ("digits.csv.gz", {"delimiter": ","}, [1797, 65]),
    ],
)
def test_load_gzip_compressed_csv_data(filename, kwargs, expected_shape):
    actual_data = load_gzip_compressed_csv_data(filename, **kwargs)
    assert actual_data.shape == tuple(expected_shape)


def test_load_gzip_compressed_csv_data_with_descr():
    data_file_name = "diabetes_target.csv.gz"
    descr_file_name = "diabetes.rst"

    expected_data = load_gzip_compressed_csv_data(data_file_name=data_file_name)
    actual_data, descr = load_gzip_compressed_csv_data(
        data_file_name=data_file_name,
        descr_file_name=descr_file_name,
    )

    np.testing.assert_array_equal(actual_data, expected_data)
    assert descr.startswith(".. _diabetes_dataset:")


def test_load_sample_images():
    try:
        res = load_sample_images()
        assert len(res.images) == 2
        assert len(res.filenames) == 2
        images = res.images

        # assert is china image
        assert np.all(images[0][0, 0, :] == np.array([174, 201, 231], dtype=np.uint8))
        # assert is flower image
        assert np.all(images[1][0, 0, :] == np.array([2, 19, 13], dtype=np.uint8))
        assert res.DESCR
    except ImportError:
        warnings.warn("Could not load sample images, PIL is not available.")


def test_load_sample_image():
    try:
        china = load_sample_image("china.jpg")
        assert china.dtype == "uint8"
        assert china.shape == (427, 640, 3)
    except ImportError:
        warnings.warn("Could not load sample images, PIL is not available.")


def test_load_diabetes_raw():
    """Test to check that we load a scaled version by default but that we can
    get an unscaled version when setting `scaled=False`."""
    diabetes_raw = load_diabetes(scaled=False)
    assert diabetes_raw.data.shape == (442, 10)
    assert diabetes_raw.target.size == 442
    assert len(diabetes_raw.feature_names) == 10
    assert diabetes_raw.DESCR

    diabetes_default = load_diabetes()

    np.testing.assert_allclose(
        scale(diabetes_raw.data) / (442**0.5), diabetes_default.data, atol=1e-04
    )


@pytest.mark.parametrize(
    "loader_func, data_shape, target_shape, n_target, has_descr, filenames",
    [
        (load_breast_cancer, (569, 30), (569,), 2, True, ["filename"]),
        (load_wine, (178, 13), (178,), 3, True, []),
        (load_iris, (150, 4), (150,), 3, True, ["filename"]),
        (
            load_linnerud,
            (20, 3),
            (20, 3),
            3,
            True,
            ["data_filename", "target_filename"],
        ),
        (load_diabetes, (442, 10), (442,), None, True, []),
        (load_digits, (1797, 64), (1797,), 10, True, []),
        (partial(load_digits, n_class=9), (1617, 64), (1617,), 10, True, []),
    ],
)
def test_loader(loader_func, data_shape, target_shape, n_target, has_descr, filenames):
    bunch = loader_func()

    assert isinstance(bunch, Bunch)
    assert bunch.data.shape == data_shape
    assert bunch.target.shape == target_shape
    if hasattr(bunch, "feature_names"):
        assert len(bunch.feature_names) == data_shape[1]
    if n_target is not None:
        assert len(bunch.target_names) == n_target
    if has_descr:
        assert bunch.DESCR
    if filenames:
        assert "data_module" in bunch
        assert all(
            [
                f in bunch
                and (resources.files(bunch["data_module"]) / bunch[f]).is_file()
                for f in filenames
            ]
        )


@pytest.mark.parametrize(
    "loader_func, data_dtype, target_dtype",
    [
        (load_breast_cancer, np.float64, int),
        (load_diabetes, np.float64, np.float64),
        (load_digits, np.float64, int),
        (load_iris, np.float64, int),
        (load_linnerud, np.float64, np.float64),
        (load_wine, np.float64, int),
    ],
)
def test_toy_dataset_frame_dtype(loader_func, data_dtype, target_dtype):
    default_result = loader_func()
    check_as_frame(
        default_result,
        loader_func,
        expected_data_dtype=data_dtype,
        expected_target_dtype=target_dtype,
    )


def test_loads_dumps_bunch():
    bunch = Bunch(x="x")
    bunch_from_pkl = loads(dumps(bunch))
    bunch_from_pkl.x = "y"
    assert bunch_from_pkl["x"] == bunch_from_pkl.x


def test_bunch_pickle_generated_with_0_16_and_read_with_0_17():
    bunch = Bunch(key="original")
    # This reproduces a problem when Bunch pickles have been created
    # with scikit-learn 0.16 and are read with 0.17. Basically there
    # is a surprising behaviour because reading bunch.key uses
    # bunch.__dict__ (which is non empty for 0.16 Bunch objects)
    # whereas assigning into bunch.key uses bunch.__setattr__. See
    # https://github.com/scikit-learn/scikit-learn/issues/6196 for
    # more details
    bunch.__dict__["key"] = "set from __dict__"
    bunch_from_pkl = loads(dumps(bunch))
    # After loading from pickle the __dict__ should have been ignored
    assert bunch_from_pkl.key == "original"
    assert bunch_from_pkl["key"] == "original"
    # Making sure that changing the attr does change the value
    # associated with __getitem__ as well
    bunch_from_pkl.key = "changed"
    assert bunch_from_pkl.key == "changed"
    assert bunch_from_pkl["key"] == "changed"


def test_bunch_dir():
    # check that dir (important for autocomplete) shows attributes
    data = load_iris()
    assert "data" in dir(data)


def test_load_boston_error():
    """Check that we raise the ethical warning when trying to import `load_boston`."""
    msg = "The Boston housing prices dataset has an ethical problem"
    with pytest.raises(ImportError, match=msg):
        from sklearn.datasets import load_boston  # noqa: F401

    # other non-existing function should raise the usual import error
    msg = "cannot import name 'non_existing_function' from 'sklearn.datasets'"
    with pytest.raises(ImportError, match=msg):
        from sklearn.datasets import non_existing_function  # noqa: F401


def test_fetch_remote_raise_warnings_with_invalid_url(monkeypatch):
    """Check retry mechanism in _fetch_remote."""

    url = "https://scikit-learn.org/this_file_does_not_exist.tar.gz"
    invalid_remote_file = RemoteFileMetadata("invalid_file", url, None)
    urlretrieve_mock = Mock(
        side_effect=HTTPError(
            url=url, code=404, msg="Not Found", hdrs=None, fp=io.BytesIO()
        )
    )
    monkeypatch.setattr("sklearn.datasets._base.urlretrieve", urlretrieve_mock)

    with pytest.warns(UserWarning, match="Retry downloading") as record:
        with pytest.raises(HTTPError, match="HTTP Error 404"):
            _fetch_remote(invalid_remote_file, n_retries=3, delay=0)

        assert urlretrieve_mock.call_count == 4

        for r in record:
            assert str(r.message) == f"Retry downloading from url: {url}"
        assert len(record) == 3


def test_derive_folder_and_filename_from_url():
    folder, filename = _derive_folder_and_filename_from_url(
        "https://example.com/file.tar.gz"
    )
    assert folder == "example.com"
    assert filename == "file.tar.gz"

    folder, filename = _derive_folder_and_filename_from_url(
        "https://example.com/نمونه نماینده.data"
    )
    assert folder == "example.com"
    assert filename == "نمونه-نماینده.data"

    folder, filename = _derive_folder_and_filename_from_url(
        "https://example.com/path/to-/.file.tar.gz"
    )
    assert folder == "example.com/path_to"
    assert filename == "file.tar.gz"

    folder, filename = _derive_folder_and_filename_from_url("https://example.com/")
    assert folder == "example.com"
    assert filename == "downloaded_file"

    folder, filename = _derive_folder_and_filename_from_url("https://example.com")
    assert folder == "example.com"
    assert filename == "downloaded_file"

    folder, filename = _derive_folder_and_filename_from_url(
        "https://example.com/path/@to/data.json?param=value"
    )
    assert folder == "example.com/path_to"
    assert filename == "data.json"

    folder, filename = _derive_folder_and_filename_from_url(
        "https://example.com/path/@@to._/-_.data.json.#anchor"
    )
    assert folder == "example.com/path_to"
    assert filename == "data.json"

    folder, filename = _derive_folder_and_filename_from_url(
        "https://example.com//some_file.txt"
    )
    assert folder == "example.com"
    assert filename == "some_file.txt"

    folder, filename = _derive_folder_and_filename_from_url(
        "http://example/../some_file.txt"
    )
    assert folder == "example"
    assert filename == "some_file.txt"

    folder, filename = _derive_folder_and_filename_from_url(
        "https://example.com/!.'.,/some_file.txt"
    )
    assert folder == "example.com"
    assert filename == "some_file.txt"

    folder, filename = _derive_folder_and_filename_from_url(
        "https://example.com/a/!.'.,/b/some_file.txt"
    )
    assert folder == "example.com/a_b"
    assert filename == "some_file.txt"

    folder, filename = _derive_folder_and_filename_from_url("https://example.com/!.'.,")
    assert folder == "example.com"
    assert filename == "downloaded_file"

    with pytest.raises(ValueError, match="Invalid URL"):
        _derive_folder_and_filename_from_url("https:/../")


def _mock_urlretrieve(server_side):
    def _urlretrieve_mock(url, local_path):
        server_root = Path(server_side)
        file_path = urlparse(url).path.strip("/")
        if not (server_root / file_path).exists():
            raise HTTPError(url, 404, "Not Found", None, None)
        shutil.copy(server_root / file_path, local_path)

    return Mock(side_effect=_urlretrieve_mock)


def test_fetch_file_using_data_home(monkeypatch, tmpdir):
    tmpdir = Path(tmpdir)
    server_side = tmpdir / "server_side"
    server_side.mkdir()
    data_file = server_side / "data.jsonl"
    server_data = '{"a": 1, "b": 2}\n'
    data_file.write_text(server_data, encoding="utf-8")

    server_subfolder = server_side / "subfolder"
    server_subfolder.mkdir()
    other_data_file = server_subfolder / "other_file.txt"
    other_data_file.write_text("Some important text data.", encoding="utf-8")

    data_home = tmpdir / "data_home"
    data_home.mkdir()

    urlretrieve_mock = _mock_urlretrieve(server_side)
    monkeypatch.setattr("sklearn.datasets._base.urlretrieve", urlretrieve_mock)

    monkeypatch.setattr(
        "sklearn.datasets._base.get_data_home", Mock(return_value=data_home)
    )
    fetched_file_path = fetch_file(
        "https://example.com/data.jsonl",
    )
    assert fetched_file_path == data_home / "example.com" / "data.jsonl"
    assert fetched_file_path.read_text(encoding="utf-8") == server_data

    fetched_file_path = fetch_file(
        "https://example.com/subfolder/other_file.txt",
    )
    assert (
        fetched_file_path == data_home / "example.com" / "subfolder" / "other_file.txt"
    )
    assert fetched_file_path.read_text(encoding="utf-8") == other_data_file.read_text(
        "utf-8"
    )

    expected_warning_msg = re.escape(
        "Retry downloading from url: https://example.com/subfolder/invalid.txt"
    )
    with pytest.raises(HTTPError):
        with pytest.warns(match=expected_warning_msg):
            fetch_file(
                "https://example.com/subfolder/invalid.txt",
                delay=0,
            )

    local_subfolder = data_home / "example.com" / "subfolder"
    assert sorted(local_subfolder.iterdir()) == [local_subfolder / "other_file.txt"]


def test_fetch_file_without_sha256(monkeypatch, tmpdir):
    server_side = tmpdir.mkdir("server_side")
    data_file = Path(server_side / "data.jsonl")
    server_data = '{"a": 1, "b": 2}\n'
    data_file.write_text(server_data, encoding="utf-8")

    client_side = tmpdir.mkdir("client_side")

    urlretrieve_mock = _mock_urlretrieve(server_side)
    monkeypatch.setattr("sklearn.datasets._base.urlretrieve", urlretrieve_mock)

    # The first call should trigger a download:
    fetched_file_path = fetch_file(
        "https://example.com/data.jsonl",
        folder=client_side,
    )
    assert fetched_file_path == client_side / "data.jsonl"
    assert fetched_file_path.read_text(encoding="utf-8") == server_data
    assert urlretrieve_mock.call_count == 1

    # Fetching again the same file to the same folder should do nothing:
    fetched_file_path = fetch_file(
        "https://example.com/data.jsonl",
        folder=client_side,
    )
    assert fetched_file_path == client_side / "data.jsonl"
    assert fetched_file_path.read_text(encoding="utf-8") == server_data
    assert urlretrieve_mock.call_count == 1

    # Deleting and calling again should re-download
    fetched_file_path.unlink()
    fetched_file_path = fetch_file(
        "https://example.com/data.jsonl",
        folder=client_side,
    )
    assert fetched_file_path == client_side / "data.jsonl"
    assert fetched_file_path.read_text(encoding="utf-8") == server_data
    assert urlretrieve_mock.call_count == 2


def test_fetch_file_with_sha256(monkeypatch, tmpdir):
    server_side = tmpdir.mkdir("server_side")
    data_file = Path(server_side / "data.jsonl")
    server_data = '{"a": 1, "b": 2}\n'
    data_file.write_text(server_data, encoding="utf-8")
    expected_sha256 = hashlib.sha256(data_file.read_bytes()).hexdigest()

    client_side = tmpdir.mkdir("client_side")

    urlretrieve_mock = _mock_urlretrieve(server_side)
    monkeypatch.setattr("sklearn.datasets._base.urlretrieve", urlretrieve_mock)

    # The first call should trigger a download.
    fetched_file_path = fetch_file(
        "https://example.com/data.jsonl", folder=client_side, sha256=expected_sha256
    )
    assert fetched_file_path == client_side / "data.jsonl"
    assert fetched_file_path.read_text(encoding="utf-8") == server_data
    assert urlretrieve_mock.call_count == 1

    # Fetching again the same file to the same folder should do nothing when
    # the sha256 match:
    fetched_file_path = fetch_file(
        "https://example.com/data.jsonl", folder=client_side, sha256=expected_sha256
    )
    assert fetched_file_path == client_side / "data.jsonl"
    assert fetched_file_path.read_text(encoding="utf-8") == server_data
    assert urlretrieve_mock.call_count == 1

    # Corrupting the local data should yield a warning and trigger a new download:
    fetched_file_path.write_text("corrupted contents", encoding="utf-8")
    expected_msg = (
        r"SHA256 checksum of existing local file data.jsonl "
        rf"\(.*\) differs from expected \({expected_sha256}\): "
        r"re-downloading from https://example.com/data.jsonl \."
    )
    with pytest.warns(match=expected_msg):
        fetched_file_path = fetch_file(
            "https://example.com/data.jsonl", folder=client_side, sha256=expected_sha256
        )
        assert fetched_file_path == client_side / "data.jsonl"
        assert fetched_file_path.read_text(encoding="utf-8") == server_data
        assert urlretrieve_mock.call_count == 2

    # Calling again should do nothing:
    fetched_file_path = fetch_file(
        "https://example.com/data.jsonl", folder=client_side, sha256=expected_sha256
    )
    assert fetched_file_path == client_side / "data.jsonl"
    assert fetched_file_path.read_text(encoding="utf-8") == server_data
    assert urlretrieve_mock.call_count == 2

    # Deleting the local file and calling again should redownload without warning:
    fetched_file_path.unlink()
    fetched_file_path = fetch_file(
        "https://example.com/data.jsonl", folder=client_side, sha256=expected_sha256
    )
    assert fetched_file_path == client_side / "data.jsonl"
    assert fetched_file_path.read_text(encoding="utf-8") == server_data
    assert urlretrieve_mock.call_count == 3

    # Calling without a sha256 should also work without redownloading:
    fetched_file_path = fetch_file(
        "https://example.com/data.jsonl",
        folder=client_side,
    )
    assert fetched_file_path == client_side / "data.jsonl"
    assert fetched_file_path.read_text(encoding="utf-8") == server_data
    assert urlretrieve_mock.call_count == 3

    # Calling with a wrong sha256 should raise an informative exception:
    non_matching_sha256 = "deadbabecafebeef"
    expected_warning_msg = "differs from expected"
    expected_error_msg = re.escape(
        f"The SHA256 checksum of data.jsonl ({expected_sha256}) differs from "
        f"expected ({non_matching_sha256})."
    )
    with pytest.raises(OSError, match=expected_error_msg):
        with pytest.warns(match=expected_warning_msg):
            fetch_file(
                "https://example.com/data.jsonl",
                folder=client_side,
                sha256=non_matching_sha256,
            )
