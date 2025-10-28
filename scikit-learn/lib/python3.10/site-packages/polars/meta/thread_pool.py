from __future__ import annotations

import contextlib

from polars._utils.deprecation import deprecated

with contextlib.suppress(ImportError):  # Module not available when building docs
    import polars._plr as plr

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import sys

    if sys.version_info >= (3, 13):
        from warnings import deprecated
    else:
        from typing_extensions import deprecated  # noqa: TC004


def thread_pool_size() -> int:
    """
    Return the number of threads in the Polars thread pool.

    Notes
    -----
    The thread pool size can be overridden by setting the `POLARS_MAX_THREADS`
    environment variable before process start. The thread pool is not behind a
    lock, so it cannot be modified once set. A reasonable use case for this might
    be temporarily limiting the number of threads before importing Polars in a
    PySpark UDF or similar context. Otherwise, it is strongly recommended not to
    override this value as it will be set automatically by the engine.

    Examples
    --------
    >>> pl.thread_pool_size()  # doctest: +SKIP
    16
    """
    return plr.thread_pool_size()


@deprecated("`threadpool_size` was renamed; use `thread_pool_size` instead.")
def threadpool_size() -> int:
    """
    Return the number of threads in the Polars thread pool.

    .. deprecated:: 0.20.7
        This function has been renamed to :func:`thread_pool_size`.
    """
    return thread_pool_size()
