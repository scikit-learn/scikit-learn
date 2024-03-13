import string
import timeit
import warnings

import numpy as np
import pytest

from sklearn.utils import (
    _message_with_time,
    _print_elapsed_time,
    _to_object_array,
    check_random_state,
    column_or_1d,
    deprecated,
    safe_mask,
)
from sklearn.utils._array_api import yield_namespace_device_dtype_combinations
from sklearn.utils._missing import is_scalar_nan
from sklearn.utils._testing import assert_array_equal, assert_no_warnings
from sklearn.utils.fixes import CSR_CONTAINERS
from sklearn.utils.validation import _is_polars_df


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

    with pytest.raises(ValueError):
        check_random_state("some invalid seed")


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

        assert spam == "spam"  # function must remain usable

        assert len(w) == 1
        assert issubclass(w[0].category, FutureWarning)
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
        assert issubclass(w[0].category, FutureWarning)
        assert "deprecated" in str(w[0].message).lower()


@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
def test_safe_mask(csr_container):
    random_state = check_random_state(0)
    X = random_state.rand(5, 4)
    X_csr = csr_container(X)
    mask = [False, False, True, True, True]

    mask = safe_mask(X, mask)
    assert X[mask].shape[0] == 3

    mask = safe_mask(X_csr, mask)
    assert X_csr[mask].shape[0] == 3


def test_column_or_1d():
    EXAMPLES = [
        ("binary", ["spam", "egg", "spam"]),
        ("binary", [0, 1, 0, 1]),
        ("continuous", np.arange(10) / 20.0),
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
        if y_type in ["binary", "multiclass", "continuous"]:
            assert_array_equal(column_or_1d(y), np.ravel(y))
        else:
            with pytest.raises(ValueError):
                column_or_1d(y)


@pytest.mark.parametrize(
    ["source", "message", "is_long"],
    [
        ("ABC", string.ascii_lowercase, False),
        ("ABCDEF", string.ascii_lowercase, False),
        ("ABC", string.ascii_lowercase * 3, True),
        ("ABC" * 10, string.ascii_lowercase, True),
        ("ABC", string.ascii_lowercase + "\u1048", False),
    ],
)
@pytest.mark.parametrize(
    ["time", "time_str"],
    [
        (0.2, "   0.2s"),
        (20, "  20.0s"),
        (2000, "33.3min"),
        (20000, "333.3min"),
    ],
)
def test_message_with_time(source, message, is_long, time, time_str):
    out = _message_with_time(source, message, time)
    if is_long:
        assert len(out) > 70
    else:
        assert len(out) == 70

    assert out.startswith("[" + source + "] ")
    out = out[len(source) + 3 :]

    assert out.endswith(time_str)
    out = out[: -len(time_str)]
    assert out.endswith(", total=")
    out = out[: -len(", total=")]
    assert out.endswith(message)
    out = out[: -len(message)]
    assert out.endswith(" ")
    out = out[:-1]

    if is_long:
        assert not out
    else:
        assert list(set(out)) == ["."]


@pytest.mark.parametrize(
    ["message", "expected"],
    [
        ("hello", _message_with_time("ABC", "hello", 0.1) + "\n"),
        ("", _message_with_time("ABC", "", 0.1) + "\n"),
        (None, ""),
    ],
)
def test_print_elapsed_time(message, expected, capsys, monkeypatch):
    monkeypatch.setattr(timeit, "default_timer", lambda: 0)
    with _print_elapsed_time("ABC", message):
        monkeypatch.setattr(timeit, "default_timer", lambda: 0.1)
    assert capsys.readouterr().out == expected


@pytest.mark.parametrize(
    "value, result",
    [
        (float("nan"), True),
        (np.nan, True),
        (float(np.nan), True),
        (np.float32(np.nan), True),
        (np.float64(np.nan), True),
        (0, False),
        (0.0, False),
        (None, False),
        ("", False),
        ("nan", False),
        ([np.nan], False),
        (9867966753463435747313673, False),  # Python int that overflows with C type
    ],
)
def test_is_scalar_nan(value, result):
    assert is_scalar_nan(value) is result
    # make sure that we are returning a Python bool
    assert isinstance(is_scalar_nan(value), bool)


def dummy_func():
    pass


def test_deprecation_joblib_api(tmpdir):
    # Only parallel_backend and register_parallel_backend are not deprecated in
    # sklearn.utils
    from sklearn.utils import parallel_backend, register_parallel_backend

    assert_no_warnings(parallel_backend, "loky", None)
    assert_no_warnings(register_parallel_backend, "failing", None)

    from sklearn.utils._joblib import joblib

    del joblib.parallel.BACKENDS["failing"]


@pytest.mark.parametrize("sequence", [[np.array(1), np.array(2)], [[1, 2], [3, 4]]])
def test_to_object_array(sequence):
    out = _to_object_array(sequence)
    assert isinstance(out, np.ndarray)
    assert out.dtype.kind == "O"
    assert out.ndim == 1


def test__is_polars_df():
    """Check that _is_polars_df return False for non-dataframe objects."""

    class LooksLikePolars:
        def __init__(self):
            self.columns = ["a", "b"]
            self.schema = ["a", "b"]

    assert not _is_polars_df(LooksLikePolars())
