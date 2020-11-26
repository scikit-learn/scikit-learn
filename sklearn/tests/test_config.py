import time

from joblib import Parallel
import joblib
import pytest

from sklearn import get_config, set_config, config_context
from sklearn.utils._testing import assert_raises
from sklearn.utils.fixes import delayed
from sklearn.utils.fixes import parse_version


def test_config_context():
    assert get_config() == {'assume_finite': False, 'working_memory': 1024,
                            'print_changed_only': True,
                            'display': 'text'}

    # Not using as a context manager affects nothing
    config_context(assume_finite=True)
    assert get_config()['assume_finite'] is False

    with config_context(assume_finite=True):
        assert get_config() == {'assume_finite': True, 'working_memory': 1024,
                                'print_changed_only': True,
                                'display': 'text'}
    assert get_config()['assume_finite'] is False

    with config_context(assume_finite=True):
        with config_context(assume_finite=None):
            assert get_config()['assume_finite'] is True

        assert get_config()['assume_finite'] is True

        with config_context(assume_finite=False):
            assert get_config()['assume_finite'] is False

            with config_context(assume_finite=None):
                assert get_config()['assume_finite'] is False

                # global setting will not be retained outside of context that
                # did not modify this setting
                set_config(assume_finite=True)
                assert get_config()['assume_finite'] is True

            assert get_config()['assume_finite'] is False

        assert get_config()['assume_finite'] is True

    assert get_config() == {'assume_finite': False, 'working_memory': 1024,
                            'print_changed_only': True,
                            'display': 'text'}

    # No positional arguments
    assert_raises(TypeError, config_context, True)
    # No unknown arguments
    assert_raises(TypeError, config_context(do_something_else=True).__enter__)


def test_config_context_exception():
    assert get_config()['assume_finite'] is False
    try:
        with config_context(assume_finite=True):
            assert get_config()['assume_finite'] is True
            raise ValueError()
    except ValueError:
        pass
    assert get_config()['assume_finite'] is False


def test_set_config():
    assert get_config()['assume_finite'] is False
    set_config(assume_finite=None)
    assert get_config()['assume_finite'] is False
    set_config(assume_finite=True)
    assert get_config()['assume_finite'] is True
    set_config(assume_finite=None)
    assert get_config()['assume_finite'] is True
    set_config(assume_finite=False)
    assert get_config()['assume_finite'] is False

    # No unknown arguments
    assert_raises(TypeError, set_config, do_something_else=True)


def set_assume_finite(assume_finite, sleep_dur):
    with config_context(assume_finite=assume_finite):
        time.sleep(sleep_dur)
        return get_config()['assume_finite']


@pytest.mark.parametrize("backend",
                         ["loky", "multiprocessing", "threading"])
def test_config_threadsafe(backend):
    """Test that the global config is threadsafe."""

    if (parse_version(joblib.__version__) < parse_version('0.12')
            and backend == 'loky'):
        pytest.skip('loky backend does not exist in joblib <0.12')  # noqa

    booleans = [False, True]
    sleep_seconds = [0.1, 0.2]

    items = Parallel(backend=backend, n_jobs=2)(
        delayed(set_assume_finite)(value, sleep_dur)
        for value, sleep_dur in zip(booleans, sleep_seconds))

    assert items == [False, True]
