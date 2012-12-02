import os
import shutil
import tempfile

from sklearn.datasets import get_data_home
from sklearn.datasets import clear_data_home

from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_equal


DATA_HOME = tempfile.mkdtemp(prefix="scikit_learn_data_home_test_")


def teardown_module():
    """Test fixture (clean up) run once after all tests of this module"""
    if os.path.isdir(DATA_HOME):
        shutil.rmtree(DATA_HOME)


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
