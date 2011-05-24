import numpy as np

from scikits.learn.utils import check_arrays
from scikits.learn.utils import check_random_state
from scikits.learn.utils import resample
from nose.tools import assert_raises


def test_make_rng():
    """Check the check_random_state utility function behavior"""
    assert check_random_state(None) is np.random.mtrand._rand
    assert check_random_state(np.random) is np.random.mtrand._rand

    rng_42 = np.random.RandomState(42)
    assert check_random_state(42).randint(100) == rng_42.randint(100)

    rng_42 = np.random.RandomState(42)
    assert check_random_state(rng_42) is rng_42

    rng_42 = np.random.RandomState(42)
    assert check_random_state(43).randint(100) != rng_42.randint(100)

    assert_raises(ValueError, check_random_state, "some invalid seed")


def test_resample_noarg():
    """Border case not worth mentioning in doctests"""
    assert resample() is None


def test_check_arrays_value_errors():
    """Check that invalid arguments yield ValueError"""
    assert_raises(ValueError, check_arrays, [0], [0, 1])
    assert_raises(ValueError, check_arrays, 0, [0, 1])
    assert_raises(ValueError, check_arrays, [0], 0)
    assert_raises(ValueError, check_arrays, [0, 1], [0, 1], meaning_of_life=42)


def test_resample_value_errors():
    """Check that invalid arguments yield ValueError"""
    assert_raises(ValueError, resample, [0], [0, 1])
    assert_raises(ValueError, resample, [0, 1], [0, 1], n_samples=3)
    assert_raises(ValueError, resample, [0, 1], [0, 1], meaning_of_life=42)
