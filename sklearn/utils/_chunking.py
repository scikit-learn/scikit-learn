import warnings
from itertools import islice
from numbers import Integral

import numpy as np

from .._config import get_config
from ._param_validation import Interval, validate_params


def chunk_generator(gen, chunksize):
    """Chunk generator, ``gen`` into lists of length ``chunksize``. The last
    chunk may have a length less than ``chunksize``."""
    while True:
        chunk = list(islice(gen, chunksize))
        if chunk:
            yield chunk
        else:
            return


@validate_params(
    {
        "n": [Interval(Integral, 1, None, closed="left")],
        "batch_size": [Interval(Integral, 1, None, closed="left")],
        "min_batch_size": [Interval(Integral, 0, None, closed="left")],
    },
    prefer_skip_nested_validation=True,
)
def gen_batches(n, batch_size, *, min_batch_size=0):
    """Generator to create slices containing `batch_size` elements from 0 to `n`.

    The last slice may contain less than `batch_size` elements, when
    `batch_size` does not divide `n`.

    Parameters
    ----------
    n : int
        Size of the sequence.
    batch_size : int
        Number of elements in each batch.
    min_batch_size : int, default=0
        Minimum number of elements in each batch.

    Yields
    ------
    slice of `batch_size` elements

    See Also
    --------
    gen_even_slices: Generator to create n_packs slices going up to n.

    Examples
    --------
    >>> from sklearn.utils import gen_batches
    >>> list(gen_batches(7, 3))
    [slice(0, 3, None), slice(3, 6, None), slice(6, 7, None)]
    >>> list(gen_batches(6, 3))
    [slice(0, 3, None), slice(3, 6, None)]
    >>> list(gen_batches(2, 3))
    [slice(0, 2, None)]
    >>> list(gen_batches(7, 3, min_batch_size=0))
    [slice(0, 3, None), slice(3, 6, None), slice(6, 7, None)]
    >>> list(gen_batches(7, 3, min_batch_size=2))
    [slice(0, 3, None), slice(3, 7, None)]
    """
    start = 0
    for _ in range(int(n // batch_size)):
        end = start + batch_size
        if end + min_batch_size > n:
            continue
        yield slice(start, end)
        start = end
    if start < n:
        yield slice(start, n)


@validate_params(
    {
        "n": [Interval(Integral, 1, None, closed="left")],
        "n_packs": [Interval(Integral, 1, None, closed="left")],
        "n_samples": [Interval(Integral, 1, None, closed="left"), None],
    },
    prefer_skip_nested_validation=True,
)
def gen_even_slices(n, n_packs, *, n_samples=None):
    """Generator to create `n_packs` evenly spaced slices going up to `n`.

    If `n_packs` does not divide `n`, except for the first `n % n_packs`
    slices, remaining slices may contain fewer elements.

    Parameters
    ----------
    n : int
        Size of the sequence.
    n_packs : int
        Number of slices to generate.
    n_samples : int, default=None
        Number of samples. Pass `n_samples` when the slices are to be used for
        sparse matrix indexing; slicing off-the-end raises an exception, while
        it works for NumPy arrays.

    Yields
    ------
    `slice` representing a set of indices from 0 to n.

    See Also
    --------
    gen_batches: Generator to create slices containing batch_size elements
        from 0 to n.

    Examples
    --------
    >>> from sklearn.utils import gen_even_slices
    >>> list(gen_even_slices(10, 1))
    [slice(0, 10, None)]
    >>> list(gen_even_slices(10, 10))
    [slice(0, 1, None), slice(1, 2, None), ..., slice(9, 10, None)]
    >>> list(gen_even_slices(10, 5))
    [slice(0, 2, None), slice(2, 4, None), ..., slice(8, 10, None)]
    >>> list(gen_even_slices(10, 3))
    [slice(0, 4, None), slice(4, 7, None), slice(7, 10, None)]
    """
    start = 0
    for pack_num in range(n_packs):
        this_n = n // n_packs
        if pack_num < n % n_packs:
            this_n += 1
        if this_n > 0:
            end = start + this_n
            if n_samples is not None:
                end = min(n_samples, end)
            yield slice(start, end, None)
            start = end


def get_chunk_n_rows(row_bytes, *, max_n_rows=None, working_memory=None):
    """Calculate how many rows can be processed within `working_memory`.

    Parameters
    ----------
    row_bytes : int
        The expected number of bytes of memory that will be consumed
        during the processing of each row.
    max_n_rows : int, default=None
        The maximum return value.
    working_memory : int or float, default=None
        The number of rows to fit inside this number of MiB will be
        returned. When None (default), the value of
        ``sklearn.get_config()['working_memory']`` is used.

    Returns
    -------
    int
        The number of rows which can be processed within `working_memory`.

    Warns
    -----
    Issues a UserWarning if `row_bytes exceeds `working_memory` MiB.
    """

    if working_memory is None:
        working_memory = get_config()["working_memory"]

    chunk_n_rows = int(working_memory * (2**20) // row_bytes)
    if max_n_rows is not None:
        chunk_n_rows = min(chunk_n_rows, max_n_rows)
    if chunk_n_rows < 1:
        warnings.warn(
            "Could not adhere to working_memory config. "
            "Currently %.0fMiB, %.0fMiB required."
            % (working_memory, np.ceil(row_bytes * 2**-20))
        )
        chunk_n_rows = 1
    return chunk_n_rows
