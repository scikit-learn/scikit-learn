# SPDX-License-Identifier: Apache-2.0
# Copyright 2018 The Meson development team

from __future__ import annotations

import asyncio.subprocess
import fnmatch
import itertools
import json
import signal
import sys
from pathlib import Path

from .. import mlog
from ..compilers import lang_suffixes
from ..mesonlib import quiet_git, join_args, determine_worker_count
from ..mtest import complete_all
import typing as T

Info = T.TypeVar("Info")

async def run_with_buffered_output(cmdlist: T.List[str]) -> int:
    """Run the command in cmdlist, buffering the output so that it is
       not mixed for multiple child processes.  Kill the child on
       cancellation."""
    quoted_cmdline = join_args(cmdlist)
    p: T.Optional[asyncio.subprocess.Process] = None
    try:
        p = await asyncio.create_subprocess_exec(*cmdlist,
                                                 stdin=asyncio.subprocess.DEVNULL,
                                                 stdout=asyncio.subprocess.PIPE,
                                                 stderr=asyncio.subprocess.STDOUT)
        stdo, _ = await p.communicate()
    except FileNotFoundError as e:
        print(mlog.blue('>>>'), quoted_cmdline, file=sys.stderr)
        print(mlog.red('not found:'), e.filename, file=sys.stderr)
        return 1
    except asyncio.CancelledError:
        if p:
            p.kill()
            await p.wait()
            return p.returncode or 1
        else:
            return 0

    if stdo:
        print(mlog.blue('>>>'), quoted_cmdline, flush=True)
        sys.stdout.buffer.write(stdo)
    return p.returncode

async def _run_workers(infos: T.Iterable[Info],
                       fn: T.Callable[[Info], T.Iterable[T.Coroutine[None, None, int]]]) -> int:
    futures: T.List[asyncio.Future[int]] = []
    semaphore = asyncio.Semaphore(determine_worker_count())

    async def run_one(worker_coro: T.Coroutine[None, None, int]) -> int:
        try:
            async with semaphore:
                return await worker_coro
        except asyncio.CancelledError as e:
            worker_coro.throw(e)
            return await worker_coro

    def sigterm_handler() -> None:
        for f in futures:
            f.cancel()

    if sys.platform != 'win32':
        loop = asyncio.get_running_loop()
        loop.add_signal_handler(signal.SIGINT, sigterm_handler)
        loop.add_signal_handler(signal.SIGTERM, sigterm_handler)

    for i in infos:
        futures.extend((asyncio.ensure_future(run_one(x)) for x in fn(i)))
    if not futures:
        return 0

    try:
        await complete_all(futures)
    except BaseException:
        for f in futures:
            f.cancel()
        raise

    return max(f.result() for f in futures if f.done() and not f.cancelled())

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

def all_clike_files(name: str, srcdir: Path, builddir: Path) -> T.Iterable[Path]:
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
    for f in itertools.chain.from_iterable(globs):
        strf = str(f)
        if f.is_dir() or f.suffix not in suffixes or \
                any(fnmatch.fnmatch(strf, i) for i in ignore):
            continue
        yield f

def run_clang_tool(name: str, srcdir: Path, builddir: Path, fn: T.Callable[..., T.Coroutine[None, None, int]], *args: T.Any) -> int:
    if sys.platform == 'win32':
        asyncio.set_event_loop_policy(asyncio.WindowsProactorEventLoopPolicy())

    def wrapper(path: Path) -> T.Iterable[T.Coroutine[None, None, int]]:
        yield fn(path, *args)
    return asyncio.run(_run_workers(all_clike_files(name, srcdir, builddir), wrapper))

def run_clang_tool_on_sources(name: str, srcdir: Path, builddir: Path, fn: T.Callable[..., T.Coroutine[None, None, int]], *args: T.Any) -> int:
    if sys.platform == 'win32':
        asyncio.set_event_loop_policy(asyncio.WindowsProactorEventLoopPolicy())

    source_files = set()
    with open('meson-info/intro-targets.json', encoding='utf-8') as fp:
        targets = json.load(fp)

        for target in targets:
            for target_source in target.get('target_sources') or []:
                for source in target_source.get('sources') or []:
                    source_files.add(Path(source))

    clike_files = set(all_clike_files(name, srcdir, builddir))
    source_files = source_files.intersection(clike_files)

    def wrapper(path: Path) -> T.Iterable[T.Coroutine[None, None, int]]:
        yield fn(path, *args)
    return asyncio.run(_run_workers(source_files, wrapper))

def run_tool_on_targets(fn: T.Callable[[T.Dict[str, T.Any]],
                                       T.Iterable[T.Coroutine[None, None, int]]]) -> int:
    if sys.platform == 'win32':
        asyncio.set_event_loop_policy(asyncio.WindowsProactorEventLoopPolicy())

    with open('meson-info/intro-targets.json', encoding='utf-8') as fp:
        targets = json.load(fp)
    return asyncio.run(_run_workers(targets, fn))
