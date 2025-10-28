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
from ..wrap import WrapMode, wrap

if T.TYPE_CHECKING:
    from ..compilers.rust import RustCompiler

async def run_and_confirm_success(cmdlist: T.List[str], crate: str) -> int:
    returncode = await run_with_buffered_output(cmdlist)
    if returncode == 0:
        print(mlog.green('Generated'), os.path.join('doc', crate))
    return returncode

class Rustdoc:
    def __init__(self, build: build.Build, tempdir: str, subprojects: T.Set[str]) -> None:
        self.tools: PerMachine[T.List[str]] = PerMachine([], [])
        self.warned: T.DefaultDict[str, bool] = defaultdict(lambda: False)
        self.tempdir = tempdir
        self.subprojects = subprojects
        for machine in MachineChoice:
            compilers = build.environment.coredata.compilers[machine]
            if 'rust' in compilers:
                compiler = T.cast('RustCompiler', compilers['rust'])
                self.tools[machine] = compiler.get_rust_tool('rustdoc', build.environment)

    def warn_missing_rustdoc(self, machine: str) -> None:
        if self.warned[machine]:
            return
        mlog.warning(f'rustdoc not found for {machine} machine')
        self.warned[machine] = True

    def __call__(self, target: T.Dict[str, T.Any]) -> T.Iterable[T.Coroutine[None, None, int]]:
        if target['subproject'] is not None and target['subproject'] not in self.subprojects:
            return

        for src_block in target['target_sources']:
            if 'compiler' in src_block and src_block['language'] == 'rust':
                rustdoc = getattr(self.tools, src_block['machine'])
                if not rustdoc:
                    self.warn_missing_rustdoc(src_block['machine'])
                    continue

                cmdlist = list(rustdoc)
                prev = None
                crate_name = None
                is_test = False
                for arg in src_block['parameters']:
                    if prev:
                        if prev == '--crate-name':
                            cmdlist.extend((prev, arg))
                            crate_name = arg
                        prev = None
                        continue

                    if arg == '--test':
                        is_test = True
                        break
                    elif arg in {'--crate-name', '--emit', '--out-dir', '-l', '--crate-type'}:
                        prev = arg
                    elif arg != '-g' and not arg.startswith('-l'):
                        cmdlist.append(arg)

                if is_test:
                    # --test has a completely different meaning for rustc and rustdoc;
                    # when using rust.test(), only the non-test target is documented
                    continue
                if crate_name:
                    cmdlist.extend(src_block['sources'])
                    # Assume documentation is generated for the developer's use
                    cmdlist.append('--document-private-items')
                    cmdlist.append('-o')
                    cmdlist.append('doc')
                    yield run_and_confirm_success(cmdlist, crate_name)
                else:
                    print(mlog.yellow('Skipping'), target['name'], '(no crate name)')

def get_nonwrap_subprojects(build_data: build.Build) -> T.Set[str]:
    wrap_resolver = wrap.Resolver(
        build_data.environment.get_source_dir(),
        build_data.subproject_dir,
        wrap_mode=WrapMode.nodownload)
    return set(sp
               for sp in build_data.environment.coredata.initialized_subprojects
               if sp and (sp not in wrap_resolver.wraps or wrap_resolver.wraps[sp].type is None))

def run(args: T.List[str]) -> int:
    os.chdir(args[0])
    build_data = build.load(os.getcwd())
    subproject_list = get_nonwrap_subprojects(build_data)
    with tempfile.TemporaryDirectory() as d:
        return run_tool_on_targets(Rustdoc(build_data, d, subproject_list))
