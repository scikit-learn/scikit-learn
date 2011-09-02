import os
import shutil
import tempfile

from scikits.learn.datasets import get_data_home
from scikits.learn.datasets import clear_data_home

from nose.tools import assert_false
from nose.tools import assert_true
from nose.tools import assert_equals


DATA_HOME = tempfile.mkdtemp(prefix="scikit_learn_data_home_test_")


def teardown_module():
    """Test fixture (clean up) run once after all tests of this module"""
    if os.path.isdir(DATA_HOME):
        shutil.rmtree(DATA_HOME)


def test_data_home():
    # get_data_home will point to a pre-existing folder
    data_home = get_data_home(data_home=DATA_HOME)
    assert_equals(data_home, DATA_HOME)
    assert_true(os.path.exists(data_home))

    # clear_data_home will delete both the content and the folder it-self
    clear_data_home(data_home=data_home)
    assert_false(os.path.exists(data_home))

    # if the folder is missing it will be created again
    data_home = get_data_home(data_home=DATA_HOME)
    assert_true(os.path.exists(data_home))
