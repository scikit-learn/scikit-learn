from itertools import chain, product
import warnings

import pytest
import numpy as np
import scipy.sparse as sp

from sklearn.utils.testing import (assert_equal, assert_raises, assert_true,
                                   assert_array_equal,
                                   SkipTest, assert_raises_regex,
                                   assert_warns_message, assert_no_warnings)
from sklearn.utils import check_random_state
from sklearn.utils import deprecated
from sklearn.utils import resample
from sklearn.utils import safe_mask
from sklearn.utils import column_or_1d
from sklearn.utils import safe_indexing
from sklearn.utils import shuffle
from sklearn.utils import gen_even_slices
from sklearn.utils import get_chunk_n_rows
from sklearn.utils import is_scalar_nan
from sklearn.utils.mocking import MockDataFrame
from sklearn import config_context
from pytest import approx


def test_make_rng():
    # Check the check_random_state utility function behavior
    assert_true(check_random_state(None) is np.random.mtrand._rand)
    assert_true(check_random_state(np.random) is np.random.mtrand._rand)

    rng_42 = np.random.RandomState(42)
    assert_true(check_random_state(42).randint(100) == rng_42.randint(100))

    rng_42 = np.random.RandomState(42)
    assert_true(check_random_state(rng_42) is rng_42)

    rng_42 = np.random.RandomState(42)
    assert_true(check_random_state(43).randint(100) != rng_42.randint(100))

    assert_raises(ValueError, check_random_state, "some invalid seed")


def test_deprecated():
    # Test whether the deprecated decorator issues appropriate warnings
    # Copied almost verbatim from https://docs.python.org/library/warnings.html

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


def test_resample():
    # Border case not worth mentioning in doctests
    assert_true(resample() is None)

    # Check that invalid arguments yield ValueError
    assert_raises(ValueError, resample, [0], [0, 1])
    assert_raises(ValueError, resample, [0, 1], [0, 1],
                  replace=False, n_samples=3)
    assert_raises(ValueError, resample, [0, 1], [0, 1], meaning_of_life=42)
    # Issue:6581, n_samples can be more when replace is True (default).
    assert_equal(len(resample([1, 2], n_samples=5)), 5)


def test_safe_mask():
    random_state = check_random_state(0)
    X = random_state.rand(5, 4)
    X_csr = sp.csr_matrix(X)
    mask = [False, False, True, True, True]

    mask = safe_mask(X, mask)
    assert_equal(X[mask].shape[0], 3)

    mask = safe_mask(X_csr, mask)
    assert_equal(X_csr[mask].shape[0], 3)


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
    # fun with read-only data in dataframes
    # this happens in joblib memmapping
    X.setflags(write=False)
    X_df_readonly = pd.DataFrame(X)
    inds_readonly = inds.copy()
    inds_readonly.setflags(write=False)

    for this_df, this_inds in product([X_df, X_df_readonly],
                                      [inds, inds_readonly]):
        with warnings.catch_warnings(record=True):
            X_df_indexed = safe_indexing(this_df, this_inds)

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


def test_shuffle_dont_convert_to_array():
    # Check that shuffle does not try to convert to numpy arrays with float
    # dtypes can let any indexable datastructure pass-through.
    a = ['a', 'b', 'c']
    b = np.array(['a', 'b', 'c'], dtype=object)
    c = [1, 2, 3]
    d = MockDataFrame(np.array([['a', 0],
                                ['b', 1],
                                ['c', 2]],
                      dtype=object))
    e = sp.csc_matrix(np.arange(6).reshape(3, 2))
    a_s, b_s, c_s, d_s, e_s = shuffle(a, b, c, d, e, random_state=0)

    assert_equal(a_s, ['c', 'b', 'a'])
    assert_equal(type(a_s), list)

    assert_array_equal(b_s, ['c', 'b', 'a'])
    assert_equal(b_s.dtype, object)

    assert_equal(c_s, [3, 2, 1])
    assert_equal(type(c_s), list)

    assert_array_equal(d_s, np.array([['c', 2],
                                      ['b', 1],
                                      ['a', 0]],
                                     dtype=object))
    assert_equal(type(d_s), MockDataFrame)

    assert_array_equal(e_s.toarray(), np.array([[4, 5],
                                                [2, 3],
                                                [0, 1]]))


def test_gen_even_slices():
    # check that gen_even_slices contains all samples
    some_range = range(10)
    joined_range = list(chain(*[some_range[slice] for slice in
                                gen_even_slices(10, 3)]))
    assert_array_equal(some_range, joined_range)

    # check that passing negative n_chunks raises an error
    slices = gen_even_slices(10, -1)
    assert_raises_regex(ValueError, "gen_even_slices got n_packs=-1, must be"
                        " >=1", next, slices)


@pytest.mark.parametrize(
    ('row_bytes', 'max_n_rows', 'working_memory', 'expected', 'warning'),
    [(1024, None, 1, 1024, None),
     (1024, None, 0.99999999, 1023, None),
     (1023, None, 1, 1025, None),
     (1025, None, 1, 1023, None),
     (1024, None, 2, 2048, None),
     (1024, 7, 1, 7, None),
     (1024 * 1024, None, 1, 1, None),
     (1024 * 1024 + 1, None, 1, 1,
      'Could not adhere to working_memory config. '
      'Currently 1MiB, 2MiB required.'),
     ])
def test_get_chunk_n_rows(row_bytes, max_n_rows, working_memory,
                          expected, warning):
    if warning is not None:
        def check_warning(*args, **kw):
            return assert_warns_message(UserWarning, warning, *args, **kw)
    else:
        check_warning = assert_no_warnings

    actual = check_warning(get_chunk_n_rows,
                           row_bytes=row_bytes,
                           max_n_rows=max_n_rows,
                           working_memory=working_memory)

    assert actual == expected
    assert type(actual) is type(expected)
    with config_context(working_memory=working_memory):
        actual = check_warning(get_chunk_n_rows,
                               row_bytes=row_bytes,
                               max_n_rows=max_n_rows)
        assert actual == expected
        assert type(actual) is type(expected)


@pytest.mark.parametrize("value, result", [(float("nan"), True),
                                           (np.nan, True),
                                           (np.float("nan"), True),
                                           (np.float32("nan"), True),
                                           (np.float64("nan"), True),
                                           (0, False),
                                           (0., False),
                                           (None, False),
                                           ("", False),
                                           ("nan", False),
                                           ([np.nan], False)])
def test_is_scalar_nan(value, result):
    assert is_scalar_nan(value) is result

def test_init_arpack_v0():
    v0s = []
    for i in range(100):
        v0s.append(_init_arpack_v0(1000, i))
        if i > 0:
            assert not any(np.equal(v0s[i], v0s[i-1]))

    v0 = np.concatenate(v0s)
    assert np.mean(v0) == approx(0, abs=1e-2)
    assert np.std(v0) == approx(1/np.sqrt(3), abs=1e-3)
