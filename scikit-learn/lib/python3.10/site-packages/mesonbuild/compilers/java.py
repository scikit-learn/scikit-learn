# SPDX-License-Identifier: Apache-2.0
# Copyright 2012-2017 The Meson development team

from __future__ import annotations

import os
import os.path
import shutil
import subprocess
import textwrap
import typing as T

from ..mesonlib import EnvironmentException
from .compilers import Compiler
from .mixins.islinker import BasicLinkerIsCompilerMixin

if T.TYPE_CHECKING:
    from ..envconfig import MachineInfo
    from ..environment import Environment
    from ..mesonlib import MachineChoice


java_debug_args: T.Dict[bool, T.List[str]] = {
    False: ['-g:none'],
    True: ['-g']
}

class JavaCompiler(BasicLinkerIsCompilerMixin, Compiler):

    language = 'java'
    id = 'unknown'

    _WARNING_LEVELS: T.Dict[str, T.List[str]] = {
        '0': ['-nowarn'],
        '1': ['-Xlint:all'],
        '2': ['-Xlint:all', '-Xdoclint:all'],
        '3': ['-Xlint:all', '-Xdoclint:all'],
    }

    def __init__(self, exelist: T.List[str], version: str, for_machine: MachineChoice,
                 info: 'MachineInfo', full_version: T.Optional[str] = None):
        super().__init__([], exelist, version, for_machine, info, full_version=full_version)
        self.javarunner = 'java'

    def get_warn_args(self, level: str) -> T.List[str]:
        return self._WARNING_LEVELS[level]

    def get_werror_args(self) -> T.List[str]:
        return ['-Werror']

    def get_output_args(self, outputname: str) -> T.List[str]:
        if outputname == '':
            outputname = './'
        return ['-d', outputname, '-s', outputname]

    def get_pic_args(self) -> T.List[str]:
        return []

    def get_pch_use_args(self, pch_dir: str, header: str) -> T.List[str]:
        return []

    def get_pch_name(self, name: str) -> str:
        return ''

    def compute_parameters_with_absolute_paths(self, parameter_list: T.List[str],
                                               build_dir: str) -> T.List[str]:
        for idx, i in enumerate(parameter_list):
            if i in {'-cp', '-classpath', '-sourcepath'} and idx + 1 < len(parameter_list):
                path_list = parameter_list[idx + 1].split(os.pathsep)
                path_list = [os.path.normpath(os.path.join(build_dir, x)) for x in path_list]
                parameter_list[idx + 1] = os.pathsep.join(path_list)

        return parameter_list

    def sanity_check(self, work_dir: str, environment: 'Environment') -> None:
        src = 'SanityCheck.java'
        obj = 'SanityCheck'
        source_name = os.path.join(work_dir, src)
        with open(source_name, 'w', encoding='utf-8') as ofile:
            ofile.write(textwrap.dedent(
                '''class SanityCheck {
                  public static void main(String[] args) {
                    int i;
                  }
                }
                '''))
        pc = subprocess.Popen(self.exelist + [src], cwd=work_dir)
        pc.wait()
        if pc.returncode != 0:
            raise EnvironmentException(f'Java compiler {self.name_string()} cannot compile programs.')
        runner = shutil.which(self.javarunner)
        if runner:
            cmdlist = [runner, '-cp', '.', obj]
            self.run_sanity_check(environment, cmdlist, work_dir, use_exe_wrapper_for_cross=False)
        else:
            m = "Java Virtual Machine wasn't found, but it's needed by Meson. " \
                "Please install a JRE.\nIf you have specific needs where this " \
                "requirement doesn't make sense, please open a bug at " \
                "https://github.com/mesonbuild/meson/issues/new and tell us " \
                "all about it."
            raise EnvironmentException(m)

    def needs_static_linker(self) -> bool:
        return False

    def get_optimization_args(self, optimization_level: str) -> T.List[str]:
        return []

    def get_debug_args(self, is_debug: bool) -> T.List[str]:
        return java_debug_args[is_debug]
