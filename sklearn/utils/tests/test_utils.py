import warnings

import numpy as np
import scipy.sparse as sp
from scipy.linalg import pinv2

from sklearn.utils.testing import (assert_equal, assert_raises, assert_true,
                                   assert_almost_equal, assert_array_equal,
                                   SkipTest)

from sklearn.utils import check_random_state
from sklearn.utils import deprecated
from sklearn.utils import resample
from sklearn.utils import safe_mask
from sklearn.utils import column_or_1d
from sklearn.utils import safe_indexing
from sklearn.utils import shuffle
from sklearn.utils.extmath import pinvh
from sklearn.utils.mocking import MockDataFrame


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


def test_column_or_1d():
    EXAMPLES = [
        ("binary", ["spam", "egg", "spam"]),
        ("binary", [0, 1, 0, 1]),
        ("continuous", np.arange(10) / 20.),
        ("multiclass", [1, 2, 3]),
        ("multiclass", [0, 1, 2, 2, 0]),
        ("multiclass", [[1], [2], [3]]),
        ("multilabel-indicator", [[0, 1, 0], [0, 0, 1]]),
        ("multiclass-multioutput", [[1, 2, 3]]),
        ("multiclass-multioutput", [[1, 1], [2, 2], [3, 1]]),
        ("multiclass-multioutput", [[5, 1], [4, 2], [3, 1]]),
        ("multiclass-multioutput", [[1, 2, 3]]),
        ("continuous-multioutput", np.arange(30).reshape((-1, 3))),
    ]

    for y_type, y in EXAMPLES:
        if y_type in ["binary", 'multiclass', "continuous"]:
            assert_array_equal(column_or_1d(y), np.ravel(y))
        else:
            assert_raises(ValueError, column_or_1d, y)


def test_safe_indexing():
    X = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    inds = np.array([1, 2])
    X_inds = safe_indexing(X, inds)
    X_arrays = safe_indexing(np.array(X), inds)
    assert_array_equal(np.array(X_inds), X_arrays)
    assert_array_equal(np.array(X_inds), np.array(X)[inds])


def test_safe_indexing_pandas():
    try:
        import pandas as pd
    except ImportError:
        raise SkipTest("Pandas not found")
    X = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    X_df = pd.DataFrame(X)
    inds = np.array([1, 2])
    X_df_indexed = safe_indexing(X_df, inds)
    X_indexed = safe_indexing(X_df, inds)
    assert_array_equal(np.array(X_df_indexed), X_indexed)


def test_safe_indexing_mock_pandas():
    X = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    X_df = MockDataFrame(X)
    inds = np.array([1, 2])
    X_df_indexed = safe_indexing(X_df, inds)
    X_indexed = safe_indexing(X_df, inds)
    assert_array_equal(np.array(X_df_indexed), X_indexed)


def test_shuffle_on_ndim_equals_three():
    def to_tuple(A):    # to make the inner arrays hashable
        return tuple(tuple(tuple(C) for C in B) for B in A)

    A = np.array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]])  # A.shape = (2,2,2)
    S = set(to_tuple(A))
    shuffle(A)  # shouldn't raise a ValueError for dim = 3
    assert_equal(set(to_tuple(A)), S)
