from itertools import chain, product
import warnings
import string
import timeit

import pytest
import numpy as np
import scipy.sparse as sp

from sklearn.utils.testing import (assert_equal, assert_raises,
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
from sklearn.utils import _message_with_time, _print_elapsed_time
from sklearn.utils import get_chunk_n_rows
from sklearn.utils import is_scalar_nan
from sklearn.utils.mocking import MockDataFrame
from sklearn import config_context


def test_make_rng():
    # Check the check_random_state utility function behavior
    assert check_random_state(None) is np.random.mtrand._rand
    assert check_random_state(np.random) is np.random.mtrand._rand

    rng_42 = np.random.RandomState(42)
    assert check_random_state(42).randint(100) == rng_42.randint(100)

    rng_42 = np.random.RandomState(42)
    assert check_random_state(rng_42) is rng_42

    rng_42 = np.random.RandomState(42)
    assert check_random_state(43).randint(100) != rng_42.randint(100)

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

        assert spam == "spam"     # function must remain usable

        assert len(w) == 1
        assert issubclass(w[0].category, DeprecationWarning)
        assert "deprecated" in str(w[0].message).lower()

    # ... then a class.
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")

        @deprecated("don't use this")
        class Ham:
            SPAM = 1

        ham = Ham()

        assert hasattr(ham, "SPAM")

        assert len(w) == 1
        assert issubclass(w[0].category, DeprecationWarning)
        assert "deprecated" in str(w[0].message).lower()


def test_resample():
    # Border case not worth mentioning in doctests
    assert resample() is None

    # Check that invalid arguments yield ValueError
    assert_raises(ValueError, resample, [0], [0, 1])
    assert_raises(ValueError, resample, [0, 1], [0, 1],
                  replace=False, n_samples=3)
    assert_raises(ValueError, resample, [0, 1], [0, 1], meaning_of_life=42)
    # Issue:6581, n_samples can be more when replace is True (default).
    assert len(resample([1, 2], n_samples=5)) == 5


def test_resample_stratified():
    # Make sure resample can stratify
    rng = np.random.RandomState(0)
    n_samples = 100
    p = .9
    X = rng.normal(size=(n_samples, 1))
    y = rng.binomial(1, p, size=n_samples)

    _, y_not_stratified = resample(X, y, n_samples=10, random_state=0,
                                   stratify=None)
    assert np.all(y_not_stratified == 1)

    _, y_stratified = resample(X, y, n_samples=10, random_state=0, stratify=y)
    assert not np.all(y_stratified == 1)
    assert np.sum(y_stratified) == 9  # all 1s, one 0


def test_resample_stratified_replace():
    # Make sure stratified resampling supports the replace parameter
    rng = np.random.RandomState(0)
    n_samples = 100
    X = rng.normal(size=(n_samples, 1))
    y = rng.randint(0, 2, size=n_samples)

    X_replace, _ = resample(X, y, replace=True, n_samples=50,
                            random_state=rng, stratify=y)
    X_no_replace, _ = resample(X, y, replace=False, n_samples=50,
                               random_state=rng, stratify=y)
    assert np.unique(X_replace).shape[0] < 50
    assert np.unique(X_no_replace).shape[0] == 50

    # make sure n_samples can be greater than X.shape[0] if we sample with
    # replacement
    X_replace, _ = resample(X, y, replace=True, n_samples=1000,
                            random_state=rng, stratify=y)
    assert X_replace.shape[0] == 1000
    assert np.unique(X_replace).shape[0] == 100


def test_resample_stratify_2dy():
    # Make sure y can be 2d when stratifying
    rng = np.random.RandomState(0)
    n_samples = 100
    X = rng.normal(size=(n_samples, 1))
    y = rng.randint(0, 2, size=(n_samples, 2))
    X, y = resample(X, y, n_samples=50, random_state=rng, stratify=y)
    assert y.ndim == 2


def test_resample_stratify_sparse_error():
    # resample must be ndarray
    rng = np.random.RandomState(0)
    n_samples = 100
    X = rng.normal(size=(n_samples, 2))
    y = rng.randint(0, 2, size=n_samples)
    stratify = sp.csr_matrix(y)
    with pytest.raises(TypeError, match='A sparse matrix was passed'):
        X, y = resample(X, y, n_samples=50, random_state=rng,
                        stratify=stratify)


def test_safe_mask():
    random_state = check_random_state(0)
    X = random_state.rand(5, 4)
    X_csr = sp.csr_matrix(X)
    mask = [False, False, True, True, True]

    mask = safe_mask(X, mask)
    assert X[mask].shape[0] == 3

    mask = safe_mask(X_csr, mask)
    assert X_csr[mask].shape[0] == 3


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
    assert set(to_tuple(A)) == S


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

    assert a_s == ['c', 'b', 'a']
    assert type(a_s) == list

    assert_array_equal(b_s, ['c', 'b', 'a'])
    assert b_s.dtype == object

    assert c_s == [3, 2, 1]
    assert type(c_s) == list

    assert_array_equal(d_s, np.array([['c', 2],
                                      ['b', 1],
                                      ['a', 0]],
                                     dtype=object))
    assert type(d_s) == MockDataFrame

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


@pytest.mark.parametrize(
    ['source', 'message', 'is_long'],
    [
        ('ABC', string.ascii_lowercase, False),
        ('ABCDEF', string.ascii_lowercase, False),
        ('ABC', string.ascii_lowercase * 3, True),
        ('ABC' * 10, string.ascii_lowercase, True),
        ('ABC', string.ascii_lowercase + u'\u1048', False),
    ])
@pytest.mark.parametrize(
    ['time', 'time_str'],
    [
        (0.2, '   0.2s'),
        (20, '  20.0s'),
        (2000, '33.3min'),
        (20000, '333.3min'),
    ])
def test_message_with_time(source, message, is_long, time, time_str):
    out = _message_with_time(source, message, time)
    if is_long:
        assert len(out) > 70
    else:
        assert len(out) == 70

    assert out.startswith('[' + source + '] ')
    out = out[len(source) + 3:]

    assert out.endswith(time_str)
    out = out[:-len(time_str)]
    assert out.endswith(', total=')
    out = out[:-len(', total=')]
    assert out.endswith(message)
    out = out[:-len(message)]
    assert out.endswith(' ')
    out = out[:-1]

    if is_long:
        assert not out
    else:
        assert list(set(out)) == ['.']


@pytest.mark.parametrize(
    ['message', 'expected'],
    [
        ('hello', _message_with_time('ABC', 'hello', 0.1) + '\n'),
        ('', _message_with_time('ABC', '', 0.1) + '\n'),
        (None, ''),
    ])
def test_print_elapsed_time(message, expected, capsys, monkeypatch):
    monkeypatch.setattr(timeit, 'default_timer', lambda: 0)
    with _print_elapsed_time('ABC', message):
        monkeypatch.setattr(timeit, 'default_timer', lambda: 0.1)
    assert capsys.readouterr().out == expected


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


def dummy_func():
    pass


def test_deprecation_joblib_api(tmpdir):
    def check_warning(*args, **kw):
        return assert_warns_message(
            DeprecationWarning, "deprecated in version 0.20.1", *args, **kw)

    # Ensure that the joblib API is deprecated in sklearn.util
    from sklearn.utils import Parallel, Memory, delayed
    from sklearn.utils import cpu_count, hash, effective_n_jobs
    check_warning(Memory, str(tmpdir))
    check_warning(hash, 1)
    check_warning(Parallel)
    check_warning(cpu_count)
    check_warning(effective_n_jobs, 1)
    check_warning(delayed, dummy_func)

    # Only parallel_backend and register_parallel_backend are not deprecated in
    # sklearn.utils
    from sklearn.utils import parallel_backend, register_parallel_backend
    assert_no_warnings(parallel_backend, 'loky', None)
    assert_no_warnings(register_parallel_backend, 'failing', None)

    # Ensure that the deprecation have no side effect in sklearn.utils._joblib
    from sklearn.utils._joblib import Parallel, Memory, delayed
    from sklearn.utils._joblib import cpu_count, hash, effective_n_jobs
    from sklearn.utils._joblib import parallel_backend
    from sklearn.utils._joblib import register_parallel_backend
    assert_no_warnings(Memory, str(tmpdir))
    assert_no_warnings(hash, 1)
    assert_no_warnings(Parallel)
    assert_no_warnings(cpu_count)
    assert_no_warnings(effective_n_jobs, 1)
    assert_no_warnings(delayed, dummy_func)
    assert_no_warnings(parallel_backend, 'loky', None)
    assert_no_warnings(register_parallel_backend, 'failing', None)

    from sklearn.utils._joblib import joblib
    del joblib.parallel.BACKENDS['failing']
