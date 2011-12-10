from nose.tools import assert_equal, assert_in, assert_raises, assert_true
import warnings

import numpy as np

from sklearn.utils import check_arrays
from sklearn.utils import check_random_state
from sklearn.utils import deprecated
from sklearn.utils import resample


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


def test_check_arrays_exceptions():
    """Check that invalid arguments raise appropriate exceptions"""
    assert_raises(ValueError, check_arrays, [0], [0, 1])
    assert_raises(TypeError, check_arrays, 0, [0, 1])
    assert_raises(TypeError, check_arrays, [0], 0)
    assert_raises(ValueError, check_arrays, [0, 1], [0, 1], meaning_of_life=42)
    assert_raises(ValueError, check_arrays, [0], [0], sparse_format='fake')


def test_deprecated():
    """Test whether the deprecated decorator issues appropriate warnings"""
    # Copied almost verbatim from http://docs.python.org/library/warnings.html

    # First a function...
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")

        @deprecated()
        def ham():
            return "spam"

        spam = ham()

        assert_equal(spam, "spam")     # function must remain usable

        assert_equal(len(w), 1)
        assert_true(issubclass(w[0].category, DeprecationWarning))
        assert_in("deprecated", str(w[0].message).lower())

    # ... then a class.
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")

        @deprecated("don't use this")
        class Ham(object):
            SPAM = 1

        ham = Ham()

        assert_true(hasattr(ham, "SPAM"))

        assert_equal(len(w), 1)
        assert_true(issubclass(w[0].category, DeprecationWarning))
        assert_in("deprecated", str(w[0].message).lower())


def test_resample_value_errors():
    """Check that invalid arguments yield ValueError"""
    assert_raises(ValueError, resample, [0], [0, 1])
    assert_raises(ValueError, resample, [0, 1], [0, 1], n_samples=3)
    assert_raises(ValueError, resample, [0, 1], [0, 1], meaning_of_life=42)
