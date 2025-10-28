# SPDX-License-Identifier: Apache-2.0
# Copyright 2024 The Meson development team

from __future__ import annotations
from collections import defaultdict
import os
import tempfile
import typing as T

from .run_tool import run_tool_on_targets, run_with_buffered_output
from .. import build, mlog
from ..mesonlib import MachineChoice, PerMachine

if T.TYPE_CHECKING:
    from ..compilers.rust import RustCompiler

class ClippyDriver:
    def __init__(self, build: build.Build, tempdir: str):
        self.tools: PerMachine[T.List[str]] = PerMachine([], [])
        self.warned: T.DefaultDict[str, bool] = defaultdict(lambda: False)
        self.tempdir = tempdir
        for machine in MachineChoice:
            compilers = build.environment.coredata.compilers[machine]
            if 'rust' in compilers:
                compiler = T.cast('RustCompiler', compilers['rust'])
                self.tools[machine] = compiler.get_rust_tool('clippy-driver', build.environment)

    def warn_missing_clippy(self, machine: str) -> None:
        if self.warned[machine]:
            return
        mlog.warning(f'clippy-driver not found for {machine} machine')
        self.warned[machine] = True

    def __call__(self, target: T.Dict[str, T.Any]) -> T.Iterable[T.Coroutine[None, None, int]]:
        for src_block in target['target_sources']:
            if 'compiler' in src_block and src_block['language'] == 'rust':
                clippy = getattr(self.tools, src_block['machine'])
                if not clippy:
                    self.warn_missing_clippy(src_block['machine'])
                    continue

                cmdlist = list(clippy)
                prev = None
                lints_cap = None
                for arg in src_block['parameters']:
                    if prev == '--cap-lints':
                        cmdlist.append(prev)
                        lints_cap = arg
                        prev = None
                    elif prev:
                        prev = None
                        continue
                    if arg in {'--emit', '--out-dir', '--cap-lints'}:
                        prev = arg
                    else:
                        cmdlist.append(arg)

                # no use in running clippy if it wouldn't print anything anyway
                if lints_cap == 'allow':
                    break

                cmdlist.extend(src_block['sources'])
                # the default for --emit is to go all the way to linking,
                # and --emit dep-info= is not enough for clippy to do
                # enough analysis, so use --emit metadata.
                cmdlist.append('--emit')
                cmdlist.append('metadata')
                cmdlist.append('--out-dir')
                cmdlist.append(self.tempdir)
                yield run_with_buffered_output(cmdlist)

def run(args: T.List[str]) -> int:
    os.chdir(args[0])
    build_data = build.load(os.getcwd())
    with tempfile.TemporaryDirectory() as d:
        return run_tool_on_targets(ClippyDriver(build_data, d))
