import os
import shutil
import tempfile
import nose

from sklearn.datasets import get_data_home
from sklearn.datasets import clear_data_home
from sklearn.datasets import load_filenames
from sklearn.datasets import load_files

from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_equal


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
    sample_file.write("Hello World!\n")
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


def test_deprecated_load_filenames():
    res = load_filenames(LOAD_FILES_ROOT)
    assert_true(res)


@nose.tools.with_setup(setup_load_files, teardown_load_files)
def test_default_load_files():
    res = load_files(LOAD_FILES_ROOT)
    assert_equal(len(res.filenames), 1)
    assert_equal(len(res.target_names), 2)    
    assert_equal(res.DESCR, None)
    assert_equal(res.data, ["Hello World!\n"])


@nose.tools.with_setup(setup_load_files, teardown_load_files)
def test_load_files_w_categories_desc_and_charset():
    category = os.path.abspath(TEST_CATEGORY_DIR1).split('/').pop()
    res = load_files(LOAD_FILES_ROOT, description="test",
                     categories=category, charset="utf-8")
    assert_equal(len(res.filenames), 1)
    assert_equal(len(res.target_names), 1)    
    assert_equal(res.DESCR, "test")
    assert_equal(res.data, ["Hello World!\n"])


@nose.tools.with_setup(setup_load_files, teardown_load_files)
def test_load_files_wo_load_content():
    res = load_files(LOAD_FILES_ROOT, load_content=False)
    assert_equal(len(res.filenames), 1)
    assert_equal(len(res.target_names), 2)    
    assert_equal(res.DESCR, None)
    assert_equal(res.get('data'), None)
