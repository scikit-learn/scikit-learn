"""Caching of formatted files with feature-based invalidation."""

import os
import pickle
from pathlib import Path
import tempfile
from typing import Dict, Iterable, Set, Tuple

from appdirs import user_cache_dir

from black.mode import Mode

from _black_version import version as __version__


# types
Timestamp = float
FileSize = int
CacheInfo = Tuple[Timestamp, FileSize]
Cache = Dict[str, CacheInfo]


CACHE_DIR = Path(user_cache_dir("black", version=__version__))


def read_cache(mode: Mode) -> Cache:
    """Read the cache if it exists and is well formed.

    If it is not well formed, the call to write_cache later should resolve the issue.
    """
    cache_file = get_cache_file(mode)
    if not cache_file.exists():
        return {}

    with cache_file.open("rb") as fobj:
        try:
            cache: Cache = pickle.load(fobj)
        except (pickle.UnpicklingError, ValueError):
            return {}

    return cache


def get_cache_file(mode: Mode) -> Path:
    return CACHE_DIR / f"cache.{mode.get_cache_key()}.pickle"


def get_cache_info(path: Path) -> CacheInfo:
    """Return the information used to check if a file is already formatted or not."""
    stat = path.stat()
    return stat.st_mtime, stat.st_size


def filter_cached(cache: Cache, sources: Iterable[Path]) -> Tuple[Set[Path], Set[Path]]:
    """Split an iterable of paths in `sources` into two sets.

    The first contains paths of files that modified on disk or are not in the
    cache. The other contains paths to non-modified files.
    """
    todo, done = set(), set()
    for src in sources:
        res_src = src.resolve()
        if cache.get(str(res_src)) != get_cache_info(res_src):
            todo.add(src)
        else:
            done.add(src)
    return todo, done


def write_cache(cache: Cache, sources: Iterable[Path], mode: Mode) -> None:
    """Update the cache file."""
    cache_file = get_cache_file(mode)
    try:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        new_cache = {
            **cache,
            **{str(src.resolve()): get_cache_info(src) for src in sources},
        }
        with tempfile.NamedTemporaryFile(dir=str(cache_file.parent), delete=False) as f:
            pickle.dump(new_cache, f, protocol=4)
        os.replace(f.name, cache_file)
    except OSError:
        pass
