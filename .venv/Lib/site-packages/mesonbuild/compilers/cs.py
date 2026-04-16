# SPDX-License-Identifier: Apache-2.0
# Copyright 2012-2017 The Meson development team

from __future__ import annotations

import os.path
import textwrap
import typing as T

from ..linkers import RSPFileSyntax

from .compilers import Compiler
from .mixins.islinker import BasicLinkerIsCompilerMixin

if T.TYPE_CHECKING:
    from ..dependencies import Dependency
    from ..environment import Environment
    from ..mesonlib import MachineChoice

cs_optimization_args: T.Dict[str, T.List[str]] = {
                        'plain': [],
                        '0': [],
                        'g': [],
                        '1': ['-optimize+'],
                        '2': ['-optimize+'],
                        '3': ['-optimize+'],
                        's': ['-optimize+'],
                        }


class CsCompiler(BasicLinkerIsCompilerMixin, Compiler):

    language = 'cs'

    def __init__(self, exelist: T.List[str], version: str, for_machine: MachineChoice,
                 env: Environment, runner: T.Optional[str] = None):
        super().__init__([], exelist, version, for_machine, env)
        self.runner = runner

    @classmethod
    def get_display_language(cls) -> str:
        return 'C sharp'

    def get_always_args(self) -> T.List[str]:
        return ['/nologo']

    def get_linker_always_args(self) -> T.List[str]:
        return ['/nologo']

    def get_output_args(self, fname: str) -> T.List[str]:
        return ['-out:' + fname]

    def get_link_args(self, fname: str) -> T.List[str]:
        return ['-r:' + fname]

    def get_werror_args(self) -> T.List[str]:
        return ['-warnaserror']

    def get_pic_args(self) -> T.List[str]:
        return []

    def get_dependency_compile_args(self, dep: Dependency) -> T.List[str]:
        # Historically we ignored all compile args.  Accept what we can, but
        # filter out -I arguments, which are in some pkg-config files and
        # aren't accepted by mcs.
        return [a for a in dep.get_compile_args() if not a.startswith('-I')]

    def compute_parameters_with_absolute_paths(self, parameter_list: T.List[str],
                                               build_dir: str) -> T.List[str]:
        for idx, i in enumerate(parameter_list):
            if i[:2] == '-L':
                parameter_list[idx] = i[:2] + os.path.normpath(os.path.join(build_dir, i[2:]))
            if i[:5] == '-lib:':
                parameter_list[idx] = i[:5] + os.path.normpath(os.path.join(build_dir, i[5:]))

        return parameter_list

    def get_pch_use_args(self, pch_dir: str, header: str) -> T.List[str]:
        return []

    def get_pch_name(self, header_name: str) -> str:
        return ''

    def _sanity_check_source_code(self) -> str:
        return textwrap.dedent('''
            public class Sanity {
                static public void Main () {
                }
            }
            ''')

    def _sanity_check_run_with_exe_wrapper(self, command: T.List[str]) -> T.List[str]:
        if self.runner:
            return [self.runner] + command
        return command

    def needs_static_linker(self) -> bool:
        return False

    def get_debug_args(self, is_debug: bool) -> T.List[str]:
        return ['-debug'] if is_debug else []

    def get_optimization_args(self, optimization_level: str) -> T.List[str]:
        return cs_optimization_args[optimization_level]


class MonoCompiler(CsCompiler):

    id = 'mono'

    def __init__(self, exelist: T.List[str], version: str, for_machine: MachineChoice,
                 env: Environment):
        super().__init__(exelist, version, for_machine, env, runner='mono')

    def rsp_file_syntax(self) -> 'RSPFileSyntax':
        return RSPFileSyntax.GCC


class VisualStudioCsCompiler(CsCompiler):

    id = 'csc'

    def get_debug_args(self, is_debug: bool) -> T.List[str]:
        if is_debug:
            return ['-debug'] if self.info.is_windows() else ['-debug:portable']
        else:
            return []

    def rsp_file_syntax(self) -> 'RSPFileSyntax':
        return RSPFileSyntax.MSVC
