from sklearn import get_config, set_config
from sklearn.utils.testing import assert_equal


def test_set_config():
    assert_equal(get_config(), {'assume_finite': False})

    # Not using as a context manager affects nothing
    set_config(assume_finite=True)
    assert_equal(get_config(), {'assume_finite': False})

    with set_config(assume_finite=True):
        assert_equal(get_config(), {'assume_finite': True})
    assert_equal(get_config(), {'assume_finite': False})

    with set_config(assume_finite=True):
        with set_config(assume_finite=None):
            assert_equal(get_config(), {'assume_finite': True})

        with set_config(assume_finite=False):
            assert_equal(get_config(), {'assume_finite': False})

            with set_config(assume_finite=None):
                assert_equal(get_config(), {'assume_finite': False})

    assert_equal(get_config(), {'assume_finite': False})
