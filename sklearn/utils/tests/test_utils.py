import warnings

import numpy as np
import scipy.sparse as sp
from scipy.linalg import pinv2

from nose.tools import assert_equal, assert_raises, assert_true
from numpy.testing import assert_almost_equal

from sklearn.utils import check_random_state
from sklearn.utils import deprecated
from sklearn.utils import resample
from sklearn.utils import safe_mask
from sklearn.utils.extmath import pinvh


def test_make_rng():
    """Check the check_random_state utility function behavior"""
    assert_true(check_random_state(None) is np.random.mtrand._rand)
    assert_true(check_random_state(np.random) is np.random.mtrand._rand)

    rng_42 = np.random.RandomState(42)
    assert_true(check_random_state(42).randint(100) == rng_42.randint(100))

    rng_42 = np.random.RandomState(42)
    assert_true(check_random_state(rng_42) is rng_42)

    rng_42 = np.random.RandomState(42)
    assert_true(check_random_state(43).randint(100) != rng_42.randint(100))

    assert_raises(ValueError, check_random_state, "some invalid seed")


def test_resample_noarg():
    """Border case not worth mentioning in doctests"""
    assert_true(resample() is None)


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
        assert_true("deprecated" in str(w[0].message).lower())

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
        assert_true("deprecated" in str(w[0].message).lower())


def test_resample_value_errors():
    """Check that invalid arguments yield ValueError"""
    assert_raises(ValueError, resample, [0], [0, 1])
    assert_raises(ValueError, resample, [0, 1], [0, 1], n_samples=3)
    assert_raises(ValueError, resample, [0, 1], [0, 1], meaning_of_life=42)


def test_safe_mask():
    random_state = check_random_state(0)
    X = random_state.rand(5, 4)
    X_csr = sp.csr_matrix(X)
    mask = [False, False, True, True, True]

    mask = safe_mask(X, mask)
    assert_equal(X[mask].shape[0], 3)

    mask = safe_mask(X_csr, mask)
    assert_equal(X_csr[mask].shape[0], 3)


def test_pinvh_simple_real():
    a = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 10]], dtype=np.float64)
    a = np.dot(a, a.T)
    a_pinv = pinvh(a)
    assert_almost_equal(np.dot(a, a_pinv), np.eye(3))


def test_pinvh_nonpositive():
    a = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=np.float64)
    a = np.dot(a, a.T)
    u, s, vt = np.linalg.svd(a)
    s[0] *= -1
    a = np.dot(u * s, vt)  # a is now symmetric non-positive and singular
    a_pinv = pinv2(a)
    a_pinvh = pinvh(a)
    assert_almost_equal(a_pinv, a_pinvh)


def test_pinvh_simple_complex():
    a = (np.array([[1, 2, 3], [4, 5, 6], [7, 8, 10]])
         + 1j * np.array([[10, 8, 7], [6, 5, 4], [3, 2, 1]]))
    a = np.dot(a, a.conj().T)
    a_pinv = pinvh(a)
    assert_almost_equal(np.dot(a, a_pinv), np.eye(3))
