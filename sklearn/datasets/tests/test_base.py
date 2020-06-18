import os
import shutil
import tempfile
import warnings
import numpy
from pickle import loads
from pickle import dumps
from functools import partial

import pytest

import numpy as np
from sklearn.datasets import get_data_home
from sklearn.datasets import clear_data_home
from sklearn.datasets import load_files
from sklearn.datasets import load_sample_images
from sklearn.datasets import load_sample_image
from sklearn.datasets import load_digits
from sklearn.datasets import load_diabetes
from sklearn.datasets import load_linnerud
from sklearn.datasets import load_iris
from sklearn.datasets import load_breast_cancer
from sklearn.datasets import load_boston
from sklearn.datasets import load_wine
from sklearn.utils import Bunch
from sklearn.datasets.tests.test_common import check_as_frame

from sklearn.externals._pilutil import pillow_installed

from sklearn.utils import IS_PYPY


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
    sample_file = tempfile.NamedTemporaryFile(dir=test_category_dir1,
                                              delete=False)
    sample_file.write(b"Hello World!\n")
    sample_file.close()
    yield str(test_category_dir1)
    _remove_dir(test_category_dir1)


@pytest.fixture
def test_category_dir_2(load_files_root):
    test_category_dir2 = tempfile.mkdtemp(dir=load_files_root)
    yield str(test_category_dir2)
    _remove_dir(test_category_dir2)


def test_data_home(data_home):
    # get_data_home will point to a pre-existing folder
    data_home = get_data_home(data_home=data_home)
    assert data_home == data_home
    assert os.path.exists(data_home)

    # clear_data_home will delete both the content and the folder it-self
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


def test_default_load_files(test_category_dir_1, test_category_dir_2,
                            load_files_root):
    if IS_PYPY:
        pytest.xfail('[PyPy] fails due to string containing NUL characters')
    res = load_files(load_files_root)
    assert len(res.filenames) == 1
    assert len(res.target_names) == 2
    assert res.DESCR is None
    assert res.data == [b"Hello World!\n"]


def test_load_files_w_categories_desc_and_encoding(
        test_category_dir_1, test_category_dir_2, load_files_root):
    if IS_PYPY:
        pytest.xfail('[PyPy] fails due to string containing NUL characters')
    category = os.path.abspath(test_category_dir_1).split('/').pop()
    res = load_files(load_files_root, description="test",
                     categories=category, encoding="utf-8")
    assert len(res.filenames) == 1
    assert len(res.target_names) == 1
    assert res.DESCR == "test"
    assert res.data == ["Hello World!\n"]


def test_load_files_wo_load_content(
        test_category_dir_1, test_category_dir_2, load_files_root):
    res = load_files(load_files_root, load_content=False)
    assert len(res.filenames) == 1
    assert len(res.target_names) == 2
    assert res.DESCR is None
    assert res.get('data') is None


def test_load_sample_images():
    try:
        res = load_sample_images()
        assert len(res.images) == 2
        assert len(res.filenames) == 2
        images = res.images

        # assert is china image
        assert np.all(images[0][0, 0, :] ==
                      np.array([174, 201, 231], dtype=np.uint8))
        # assert is flower image
        assert np.all(images[1][0, 0, :] ==
                      np.array([2, 19, 13], dtype=np.uint8))
        assert res.DESCR
    except ImportError:
        warnings.warn("Could not load sample images, PIL is not available.")


def test_load_sample_image():
    try:
        china = load_sample_image('china.jpg')
        assert china.dtype == 'uint8'
        assert china.shape == (427, 640, 3)
    except ImportError:
        warnings.warn("Could not load sample images, PIL is not available.")


def test_load_missing_sample_image_error():
    if pillow_installed:
        with pytest.raises(AttributeError):
            load_sample_image('blop.jpg')
    else:
        warnings.warn("Could not load sample images, PIL is not available.")


@pytest.mark.parametrize(
    "loader_func, data_shape, target_shape, n_target, has_descr, filenames",
    [(load_breast_cancer, (569, 30), (569,), 2, True, ["filename"]),
     (load_wine, (178, 13), (178,), 3, True, []),
     (load_iris, (150, 4), (150,), 3, True, ["filename"]),
     (load_linnerud, (20, 3), (20, 3), 3, True,
      ["data_filename", "target_filename"]),
     (load_diabetes, (442, 10), (442,), None, True, []),
     (load_digits, (1797, 64), (1797,), 10, True, []),
     (partial(load_digits, n_class=9), (1617, 64), (1617,), 10, True, []),
     (load_boston, (506, 13), (506,), None, True, ["filename"])]
)
def test_loader(loader_func, data_shape, target_shape, n_target, has_descr,
                filenames):
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
        assert all([os.path.exists(bunch.get(f, False)) for f in filenames])


@pytest.mark.parametrize("loader_func, data_dtype, target_dtype", [
    (load_breast_cancer, np.float64, np.int64),
    (load_diabetes, np.float64, np.float64),
    (load_digits, np.float64, np.int64),
    (load_iris, np.float64, np.int64),
    (load_linnerud, np.float64, np.float64),
    (load_wine, np.float64, np.int64),
])
def test_toy_dataset_frame_dtype(loader_func, data_dtype, target_dtype):
    default_result = loader_func()
    check_as_frame(default_result, loader_func,
                   expected_data_dtype=data_dtype,
                   expected_target_dtype=target_dtype)


def test_loads_dumps_bunch():
    bunch = Bunch(x="x")
    bunch_from_pkl = loads(dumps(bunch))
    bunch_from_pkl.x = "y"
    assert bunch_from_pkl['x'] == bunch_from_pkl.x


def test_bunch_pickle_generated_with_0_16_and_read_with_0_17():
    bunch = Bunch(key='original')
    # This reproduces a problem when Bunch pickles have been created
    # with scikit-learn 0.16 and are read with 0.17. Basically there
    # is a surprising behaviour because reading bunch.key uses
    # bunch.__dict__ (which is non empty for 0.16 Bunch objects)
    # whereas assigning into bunch.key uses bunch.__setattr__. See
    # https://github.com/scikit-learn/scikit-learn/issues/6196 for
    # more details
    bunch.__dict__['key'] = 'set from __dict__'
    bunch_from_pkl = loads(dumps(bunch))
    # After loading from pickle the __dict__ should have been ignored
    assert bunch_from_pkl.key == 'original'
    assert bunch_from_pkl['key'] == 'original'
    # Making sure that changing the attr does change the value
    # associated with __getitem__ as well
    bunch_from_pkl.key = 'changed'
    assert bunch_from_pkl.key == 'changed'
    assert bunch_from_pkl['key'] == 'changed'


def test_bunch_dir():
    # check that dir (important for autocomplete) shows attributes
    data = load_iris()
    assert "data" in dir(data)
