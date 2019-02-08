import os
import shutil
import tempfile
import warnings
import numpy
from pickle import loads
from pickle import dumps
from functools import partial

import pytest

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
from sklearn.datasets import SampleImages   # noqa: F401
from sklearn.datasets import Digits         # noqa: F401
from sklearn.datasets import Diabetes       # noqa: F401
from sklearn.datasets import Linnerud       # noqa: F401
from sklearn.datasets import Iris           # noqa: F401
from sklearn.datasets import BreastCancer   # noqa: F401
from sklearn.datasets import Boston         # noqa: F401
from sklearn.datasets import Wine           # noqa: F401
from sklearn.datasets.base import Bunch
from sklearn.datasets.tests.test_common import check_return_X_y

from sklearn.externals._pilutil import pillow_installed

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises

import numpy.testing.utils as np_test_util


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
    assert_equal(data_home, data_home)
    assert os.path.exists(data_home)

    # clear_data_home will delete both the content and the folder it-self
    clear_data_home(data_home=data_home)
    assert not os.path.exists(data_home)

    # if the folder is missing it will be created again
    data_home = get_data_home(data_home=data_home)
    assert os.path.exists(data_home)


def test_default_empty_load_files(load_files_root):
    res = load_files(load_files_root)
    assert_equal(len(res.filenames), 0)
    assert_equal(len(res.target_names), 0)
    assert_equal(res.DESCR, None)


def test_default_load_files(test_category_dir_1, test_category_dir_2,
                            load_files_root):
    res = load_files(load_files_root)
    assert_equal(len(res.filenames), 1)
    assert_equal(len(res.target_names), 2)
    assert_equal(res.DESCR, None)
    assert_equal(res.data, [b"Hello World!\n"])


def test_load_files_w_categories_desc_and_encoding(
        test_category_dir_1, test_category_dir_2, load_files_root):
    category = os.path.abspath(test_category_dir_1).split('/').pop()
    res = load_files(load_files_root, description="test",
                     categories=category, encoding="utf-8")
    assert_equal(len(res.filenames), 1)
    assert_equal(len(res.target_names), 1)
    assert_equal(res.DESCR, "test")
    assert_equal(res.data, ["Hello World!\n"])


@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_load_files_wo_load_content(
        test_category_dir_1, test_category_dir_2, load_files_root):
    res = load_files(load_files_root, load_content=False)
    assert_equal(len(res.filenames), 1)
    assert_equal(len(res.target_names), 2)
    assert_equal(res.DESCR, None)
    assert_equal(res.get('data'), None)


@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_load_sample_images():
    try:
        res = load_sample_images()
        assert_equal(len(res.images), 2)
        assert_equal(len(res.filenames), 2)
        assert res.DESCR
    except ImportError:
        warnings.warn("Could not load sample images, PIL is not available.")


@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_load_digits():
    digits = load_digits()
    assert_equal(digits.data.shape, (1797, 64))
    assert_equal(numpy.unique(digits.target).size, 10)

    # test return_X_y option
    check_return_X_y(digits, partial(load_digits))


@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_load_digits_n_class_lt_10():
    digits = load_digits(9)
    assert_equal(digits.data.shape, (1617, 64))
    assert_equal(numpy.unique(digits.target).size, 9)


@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_load_sample_image():
    try:
        china = load_sample_image('china.jpg')
        assert_equal(china.dtype, 'uint8')
        assert_equal(china.shape, (427, 640, 3))
    except ImportError:
        warnings.warn("Could not load sample images, PIL is not available.")


@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_load_missing_sample_image_error():
    if pillow_installed:
        assert_raises(AttributeError, load_sample_image,
                      'blop.jpg')
    else:
        warnings.warn("Could not load sample images, PIL is not available.")


@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_load_diabetes():
    res = load_diabetes()
    assert_equal(res.data.shape, (442, 10))
    assert res.target.size, 442
    assert_equal(len(res.feature_names), 10)
    assert res.DESCR

    # test return_X_y option
    check_return_X_y(res, partial(load_diabetes))


@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_load_linnerud():
    res = load_linnerud()
    assert_equal(res.data.shape, (20, 3))
    assert_equal(res.target.shape, (20, 3))
    assert_equal(len(res.target_names), 3)
    assert res.DESCR
    assert os.path.exists(res.data_filename)
    assert os.path.exists(res.target_filename)

    # test return_X_y option
    check_return_X_y(res, partial(load_linnerud))


@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_load_iris():
    res = load_iris()
    assert_equal(res.data.shape, (150, 4))
    assert_equal(res.target.size, 150)
    assert_equal(res.target_names.size, 3)
    assert res.DESCR
    assert os.path.exists(res.filename)

    # test return_X_y option
    check_return_X_y(res, partial(load_iris))


@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_load_wine():
    res = load_wine()
    assert_equal(res.data.shape, (178, 13))
    assert_equal(res.target.size, 178)
    assert_equal(res.target_names.size, 3)
    assert res.DESCR

    # test return_X_y option
    check_return_X_y(res, partial(load_wine))


@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_load_breast_cancer():
    res = load_breast_cancer()
    assert_equal(res.data.shape, (569, 30))
    assert_equal(res.target.size, 569)
    assert_equal(res.target_names.size, 2)
    assert res.DESCR
    assert os.path.exists(res.filename)

    # test return_X_y option
    check_return_X_y(res, partial(load_breast_cancer))


@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_load_boston():
    res = load_boston()
    assert_equal(res.data.shape, (506, 13))
    assert_equal(res.target.size, 506)
    assert_equal(res.feature_names.size, 13)
    assert res.DESCR
    assert os.path.exists(res.filename)

    # test return_X_y option
    check_return_X_y(res, partial(load_boston))


def test_loads_dumps_bunch():
    bunch = Bunch(x="x")
    bunch_from_pkl = loads(dumps(bunch))
    bunch_from_pkl.x = "y"
    assert_equal(bunch_from_pkl['x'], bunch_from_pkl.x)


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
    assert_equal(bunch_from_pkl.key, 'original')
    assert_equal(bunch_from_pkl['key'], 'original')
    # Making sure that changing the attr does change the value
    # associated with __getitem__ as well
    bunch_from_pkl.key = 'changed'
    assert_equal(bunch_from_pkl.key, 'changed')
    assert_equal(bunch_from_pkl['key'], 'changed')


@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_bunch_dir():
    # check that dir (important for autocomplete) shows attributes
    data = load_iris()
    assert "data" in dir(data)


@pytest.mark.parametrize('test_loader,expected_target_dtype', [
    ('Iris', 'int'),
    ('Boston', 'float'),
    ('BreastCancer', 'int'),
    ('Digits', 'int'),
    ('Diabetes', 'int'),
    ('Linnerud', 'float'),
    ('Wine', 'int')
])
def test_dataset_loader_dtype(test_loader, expected_target_dtype):
    assert eval(test_loader)().load().target.dtype == expected_target_dtype


@pytest.mark.parametrize('test_loader', [
    'Iris',
    'Boston',
    'BreastCancer',
    'Digits',
    'Diabetes',
    'Linnerud',
    'Wine'
])
def test_dataset_loader_bunch_paths(test_loader):
    bunch = eval(test_loader)().load()
    paths = bunch.filename, bunch.data_filename, bunch.target_filename
    assert all(list(map(os.path.exists, paths))) is True


@pytest.mark.parametrize('test_loader,exp_features,exp_targets,exp_n', [
    ('Iris', 4, 3, 150),
    ('Boston', 13, 1, 506),
    ('BreastCancer', 30, 2, 569),
    ('Digits', 64, 10, 1797),
    ('Diabetes', 10, 1, 442),
    ('Linnerud', 3, 3, 20),
    ('Wine', 13, 3, 178)])
def test_dataset_loader_shape(test_loader, exp_features,
                              exp_targets, exp_n):
    bunch = eval(test_loader)().load()
    n_features, m_features = bunch.data.shape[:2]
    n_targets, m_targets = bunch.target.shape[0], bunch.target_names.size
    assert (m_features == exp_features) and \
           (m_targets == exp_targets) and \
           (n_features == n_targets == exp_n)


@pytest.mark.parametrize('test_loader', [
    'Iris',
    'Boston',
    'BreastCancer',
    'Digits',
    'Diabetes',
    'Linnerud',
    'Wine'
])
def test_dataset_loader_check_nan(test_loader):
    bunch = eval(test_loader)().load()
    data, target = bunch.data, bunch.target
    np_test_util.assert_equal(numpy.isnan(data).any(), False)
    np_test_util.assert_equal(numpy.isnan(target).any(), False)


def test_load_data_deprecated():
    from ..base import load_data
    iris_path = Iris().local_data_paths['X']
    pytest.deprecated_call(load_data, iris_path)


@pytest.mark.parametrize('test_deprecated_fun', [
    'load_wine',
    'load_iris',
    'load_digits',
    'load_diabetes',
    'load_boston',
    'load_breast_cancer',
    'load_linnerud',
    'load_sample_images'
])
def test_functional_load_deprecated(test_deprecated_fun):
    pytest.deprecated_call(eval(test_deprecated_fun))
