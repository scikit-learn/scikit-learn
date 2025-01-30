# SPDX-License-Identifier: Apache-2.0
# Copyright 2018 The Meson development team

from __future__ import annotations

import itertools
import fnmatch
import concurrent.futures
from pathlib import Path

from ..compilers import lang_suffixes
from ..mesonlib import quiet_git
import typing as T

if T.TYPE_CHECKING:
    import subprocess

def parse_pattern_file(fname: Path) -> T.List[str]:
    patterns = []
    try:
        with fname.open(encoding='utf-8') as f:
            for line in f:
                pattern = line.strip()
                if pattern and not pattern.startswith('#'):
                    patterns.append(pattern)
    except FileNotFoundError:
        pass
    return patterns

def run_tool(name: str, srcdir: Path, builddir: Path, fn: T.Callable[..., subprocess.CompletedProcess], *args: T.Any) -> int:
    patterns = parse_pattern_file(srcdir / f'.{name}-include')
    globs: T.Union[T.List[T.List[Path]], T.List[T.Generator[Path, None, None]]]
    if patterns:
        globs = [srcdir.glob(p) for p in patterns]
    else:
        r, o = quiet_git(['ls-files'], srcdir)
        if r:
            globs = [[Path(srcdir, f) for f in o.splitlines()]]
        else:
            globs = [srcdir.glob('**/*')]
    patterns = parse_pattern_file(srcdir / f'.{name}-ignore')
    ignore = [str(builddir / '*')]
    ignore.extend([str(srcdir / p) for p in patterns])
    suffixes = set(lang_suffixes['c']).union(set(lang_suffixes['cpp']))
    suffixes.add('h')
    suffixes = {f'.{s}' for s in suffixes}
    futures = []
    returncode = 0
    e = concurrent.futures.ThreadPoolExecutor()
    try:
        for f in itertools.chain(*globs):
            strf = str(f)
            if f.is_dir() or f.suffix not in suffixes or \
                    any(fnmatch.fnmatch(strf, i) for i in ignore):
                continue
            futures.append(e.submit(fn, f, *args))
        concurrent.futures.wait(
            futures,
            return_when=concurrent.futures.FIRST_EXCEPTION
        )
    finally:
        # We try to prevent new subprocesses from being started by canceling
        # the futures, but this is not water-tight: some may have started
        # between the wait being interrupted or exited and the futures being
        # canceled. (A fundamental fix would probably require the ability to
        # terminate such subprocesses upon cancellation of the future.)
        for x in futures: # Python >=3.9: e.shutdown(cancel_futures=True)
            x.cancel()
        e.shutdown()
    if futures:
        returncode = max(x.result().returncode for x in futures)
    return returncode
