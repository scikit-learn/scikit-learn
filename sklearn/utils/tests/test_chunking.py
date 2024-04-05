import warnings
from itertools import chain

import pytest

from sklearn import config_context
from sklearn.utils._chunking import gen_even_slices, get_chunk_n_rows
from sklearn.utils._testing import assert_array_equal


def test_gen_even_slices():
    # check that gen_even_slices contains all samples
    some_range = range(10)
    joined_range = list(chain(*[some_range[slice] for slice in gen_even_slices(10, 3)]))
    assert_array_equal(some_range, joined_range)


@pytest.mark.parametrize(
    ("row_bytes", "max_n_rows", "working_memory", "expected"),
    [
        (1024, None, 1, 1024),
        (1024, None, 0.99999999, 1023),
        (1023, None, 1, 1025),
        (1025, None, 1, 1023),
        (1024, None, 2, 2048),
        (1024, 7, 1, 7),
        (1024 * 1024, None, 1, 1),
    ],
)
def test_get_chunk_n_rows(row_bytes, max_n_rows, working_memory, expected):
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        actual = get_chunk_n_rows(
            row_bytes=row_bytes,
            max_n_rows=max_n_rows,
            working_memory=working_memory,
        )

    assert actual == expected
    assert type(actual) is type(expected)
    with config_context(working_memory=working_memory):
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            actual = get_chunk_n_rows(row_bytes=row_bytes, max_n_rows=max_n_rows)
        assert actual == expected
        assert type(actual) is type(expected)


def test_get_chunk_n_rows_warns():
    """Check that warning is raised when working_memory is too low."""
    row_bytes = 1024 * 1024 + 1
    max_n_rows = None
    working_memory = 1
    expected = 1

    warn_msg = (
        "Could not adhere to working_memory config. Currently 1MiB, 2MiB required."
    )
    with pytest.warns(UserWarning, match=warn_msg):
        actual = get_chunk_n_rows(
            row_bytes=row_bytes,
            max_n_rows=max_n_rows,
            working_memory=working_memory,
        )

    assert actual == expected
    assert type(actual) is type(expected)

    with config_context(working_memory=working_memory):
        with pytest.warns(UserWarning, match=warn_msg):
            actual = get_chunk_n_rows(row_bytes=row_bytes, max_n_rows=max_n_rows)
        assert actual == expected
        assert type(actual) is type(expected)
