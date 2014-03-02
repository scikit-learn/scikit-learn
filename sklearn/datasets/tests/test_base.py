import os
import shutil
import tempfile
import warnings
import nose
import numpy

from sklearn.datasets import get_data_home
from sklearn.datasets import clear_data_home
from sklearn.datasets import load_files
from sklearn.datasets import load_sample_images
from sklearn.datasets import load_sample_image
from sklearn.datasets import load_digits
from sklearn.datasets import load_diabetes
from sklearn.datasets import load_linnerud
from sklearn.datasets import load_iris
from sklearn.datasets import load_boston

from sklearn.externals.six import b, u

from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises


DATA_HOME = tempfile.mkdtemp(prefix="scikit_learn_data_home_test_")
LOAD_FILES_ROOT = tempfile.mkdtemp(prefix="scikit_learn_load_files_test_")
TEST_CATEGORY_DIR1 = ""
TEST_CATEGORY_DIR2 = ""


def _remove_dir(path):
    if os.path.isdir(path):
        shutil.rmtree(path)


def teardown_module():
    """Test fixture (clean up) run once after all tests of this module"""
    for path in [DATA_HOME, LOAD_FILES_ROOT]:
        _remove_dir(path)


def setup_load_files():
    global TEST_CATEGORY_DIR1
    global TEST_CATEGORY_DIR2
    TEST_CATEGORY_DIR1 = tempfile.mkdtemp(dir=LOAD_FILES_ROOT)
    TEST_CATEGORY_DIR2 = tempfile.mkdtemp(dir=LOAD_FILES_ROOT)
    sample_file = tempfile.NamedTemporaryFile(dir=TEST_CATEGORY_DIR1,
                                              delete=False)
    sample_file.write(b("Hello World!\n"))
    sample_file.close()


def teardown_load_files():
    _remove_dir(TEST_CATEGORY_DIR1)
    _remove_dir(TEST_CATEGORY_DIR2)


def test_data_home():
    # get_data_home will point to a pre-existing folder
    data_home = get_data_home(data_home=DATA_HOME)
    assert_equal(data_home, DATA_HOME)
    assert_true(os.path.exists(data_home))

    # clear_data_home will delete both the content and the folder it-self
    clear_data_home(data_home=data_home)
    assert_false(os.path.exists(data_home))

    # if the folder is missing it will be created again
    data_home = get_data_home(data_home=DATA_HOME)
    assert_true(os.path.exists(data_home))


def test_default_empty_load_files():
    res = load_files(LOAD_FILES_ROOT)
    assert_equal(len(res.filenames), 0)
    assert_equal(len(res.target_names), 0)
    assert_equal(res.DESCR, None)


@nose.tools.with_setup(setup_load_files, teardown_load_files)
def test_default_load_files():
    res = load_files(LOAD_FILES_ROOT)
    assert_equal(len(res.filenames), 1)
    assert_equal(len(res.target_names), 2)
    assert_equal(res.DESCR, None)
    assert_equal(res.data, [b("Hello World!\n")])


@nose.tools.with_setup(setup_load_files, teardown_load_files)
def test_load_files_w_categories_desc_and_encoding():
    category = os.path.abspath(TEST_CATEGORY_DIR1).split('/').pop()
    res = load_files(LOAD_FILES_ROOT, description="test",
                     categories=category, encoding="utf-8")
    assert_equal(len(res.filenames), 1)
    assert_equal(len(res.target_names), 1)
    assert_equal(res.DESCR, "test")
    assert_equal(res.data, [u("Hello World!\n")])


@nose.tools.with_setup(setup_load_files, teardown_load_files)
def test_load_files_wo_load_content():
    res = load_files(LOAD_FILES_ROOT, load_content=False)
    assert_equal(len(res.filenames), 1)
    assert_equal(len(res.target_names), 2)
    assert_equal(res.DESCR, None)
    assert_equal(res.get('data'), None)


def test_load_sample_images():
    try:
        res = load_sample_images()
        assert_equal(len(res.images), 2)
        assert_equal(len(res.filenames), 2)
        assert_true(res.DESCR)
    except ImportError:
        warnings.warn("Could not load sample images, PIL is not available.")


def test_load_digits():
    digits = load_digits()
    assert_equal(digits.data.shape, (1797, 64))
    assert_equal(numpy.unique(digits.target).size, 10)


def test_load_digits_n_class_lt_10():
    digits = load_digits(9)
    assert_equal(digits.data.shape, (1617, 64))
    assert_equal(numpy.unique(digits.target).size, 9)


def test_load_sample_image():
    try:
        china = load_sample_image('china.jpg')
        assert_equal(china.dtype, 'uint8')
        assert_equal(china.shape, (427, 640, 3))
    except ImportError:
        warnings.warn("Could not load sample images, PIL is not available.")


def test_load_missing_sample_image_error():
    have_PIL = True
    try:
        try:
            from scipy.misc import imread
        except ImportError:
            from scipy.misc.pilutil import imread
    except ImportError:
        have_PIL = False
    if have_PIL:
        assert_raises(AttributeError, load_sample_image,
                      'blop.jpg')
    else:
        warnings.warn("Could not load sample images, PIL is not available.")


def test_load_diabetes():
    res = load_diabetes()
    assert_equal(res.data.shape, (442, 10))
    assert_true(res.target.size, 442)


def test_load_linnerud():
    res = load_linnerud()
    assert_equal(res.data.shape, (20, 3))
    assert_equal(res.target.shape, (20, 3))
    assert_equal(len(res.target_names), 3)
    assert_true(res.DESCR)


def test_load_iris():
    res = load_iris()
    assert_equal(res.data.shape, (150, 4))
    assert_equal(res.target.size, 150)
    assert_equal(res.target_names.size, 3)
    assert_true(res.DESCR)


def test_load_boston():
    res = load_boston()
    assert_equal(res.data.shape, (506, 13))
    assert_equal(res.target.size, 506)
    assert_equal(res.feature_names.size, 13)
    assert_true(res.DESCR)
