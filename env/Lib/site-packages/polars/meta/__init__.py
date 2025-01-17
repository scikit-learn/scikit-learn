"""Public functions that provide information about the Polars package or the environment it runs in."""  # noqa: W505

from polars.meta.build import build_info
from polars.meta.index_type import get_index_type
from polars.meta.thread_pool import thread_pool_size, threadpool_size
from polars.meta.versions import show_versions

__all__ = [
    "build_info",
    "get_index_type",
    "show_versions",
    "thread_pool_size",
    "threadpool_size",
]
