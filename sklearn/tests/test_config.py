from sklearn import get_config, set_config, config_context
from sklearn.utils.testing import assert_equal, assert_raises

dict_true = {'assume_finite': True, 'show_default_parameters': True}
dict_false = {'assume_finite': False, 'show_default_parameters': True}


def test_config_context():
    assert_equal(get_config(), dict_false)

    # Not using as a context manager affects nothing
    config_context(assume_finite=True)
    assert_equal(get_config(), dict_false)

    with config_context(assume_finite=True):
        assert_equal(get_config(), dict_true)
    assert_equal(get_config(), dict_false)

    with config_context(assume_finite=True):
        with config_context(assume_finite=None):
            assert_equal(get_config(), dict_true)

        assert_equal(get_config(), dict_true)

        with config_context(assume_finite=False):
            assert_equal(get_config(), dict_false)

            with config_context(assume_finite=None):
                assert_equal(get_config(), dict_false)

                # global setting will not be retained outside of context that
                # did not modify this setting
                set_config(assume_finite=True)
                assert_equal(get_config(), dict_true)

            assert_equal(get_config(), dict_false)

        assert_equal(get_config(), dict_true)

    assert_equal(get_config(), dict_false)

    # No positional arguments
    assert_raises(TypeError, config_context, True)
    # No unknown arguments
    assert_raises(TypeError, config_context(do_something_else=True).__enter__)


def test_config_context_exception():
    assert_equal(get_config(), dict_false)
    try:
        with config_context(assume_finite=True):
            assert_equal(get_config(), dict_true)
            raise ValueError()
    except ValueError:
        pass
    assert_equal(get_config(), dict_false)


def test_set_config():
    assert_equal(get_config(), dict_false)
    set_config(assume_finite=None)
    assert_equal(get_config(), dict_false)
    set_config(assume_finite=True)
    assert_equal(get_config(), dict_true)
    set_config(assume_finite=None)
    assert_equal(get_config(), dict_true)
    set_config(assume_finite=False)
    assert_equal(get_config(), dict_false)

    # No unknown arguments
    assert_raises(TypeError, set_config, do_something_else=True)
