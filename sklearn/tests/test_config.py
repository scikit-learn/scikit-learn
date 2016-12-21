from sklearn import get_config, set_config, config_context
from sklearn.utils.testing import assert_equal, assert_raises


def test_config_context():
    assert_equal(get_config(), {'assume_finite': False})

    # Not using as a context manager affects nothing
    config_context(assume_finite=True)
    assert_equal(get_config(), {'assume_finite': False})

    with config_context(assume_finite=True):
        assert_equal(get_config(), {'assume_finite': True})
    assert_equal(get_config(), {'assume_finite': False})

    with config_context(assume_finite=True):
        with config_context(assume_finite=None):
            assert_equal(get_config(), {'assume_finite': True})

        assert_equal(get_config(), {'assume_finite': True})

        with config_context(assume_finite=False):
            assert_equal(get_config(), {'assume_finite': False})

            with config_context(assume_finite=None):
                assert_equal(get_config(), {'assume_finite': False})

                # global setting will not be retained outside of context that
                # did not modify this setting
                set_config(assume_finite=True)
                assert_equal(get_config(), {'assume_finite': True})

            assert_equal(get_config(), {'assume_finite': False})

        assert_equal(get_config(), {'assume_finite': True})

    assert_equal(get_config(), {'assume_finite': False})

    # No positional arguments
    assert_raises(TypeError, config_context, True)
    # No unknown arguments
    assert_raises(TypeError, config_context(do_something_else=True).__enter__)


def test_config_context_exception():
    assert_equal(get_config(), {'assume_finite': False})
    try:
        with config_context(assume_finite=True):
            assert_equal(get_config(), {'assume_finite': True})
            raise ValueError()
    except ValueError:
        pass
    assert_equal(get_config(), {'assume_finite': False})


def test_set_config():
    assert_equal(get_config(), {'assume_finite': False})
    set_config(assume_finite=None)
    assert_equal(get_config(), {'assume_finite': False})
    set_config(assume_finite=True)
    assert_equal(get_config(), {'assume_finite': True})
    set_config(assume_finite=None)
    assert_equal(get_config(), {'assume_finite': True})
    set_config(assume_finite=False)
    assert_equal(get_config(), {'assume_finite': False})

    # No unknown arguments
    assert_raises(TypeError, set_config, do_something_else=True)
