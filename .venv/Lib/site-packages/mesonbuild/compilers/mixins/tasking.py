# SPDX-License-Identifier: Apache-2.0
# Copyright 2012-2023 The Meson development team
from __future__ import annotations

"""Representations specific to the TASKING embedded C/C++ compiler family."""

import os
import typing as T

from ...mesonlib import EnvironmentException
from ...options import OptionKey

if T.TYPE_CHECKING:
    from ...compilers.compilers import Compiler
else:
    # This is a bit clever, for mypy we pretend that these mixins descend from
    # Compiler, so we get all of the methods and attributes defined for us, but
    # for runtime we make them descend from object (which all classes normally
    # do). This gives us DRYer type checking, with no runtime impact
    Compiler = object

tasking_buildtype_args: T.Mapping[str, T.List[str]] = {
    'plain': [],
    'debug': [],
    'debugoptimized': [],
    'release': [],
    'minsize': [],
    'custom': []
}

tasking_optimization_args: T.Mapping[str, T.List[str]] = {
    'plain': [],
    '0': ['-O0'],
    'g': ['-O1'], # There is no debug specific level, O1 is recommended by the compiler
    '1': ['-O1'],
    '2': ['-O2'],
    '3': ['-O3'],
    's': ['-Os']
}

tasking_debug_args: T.Mapping[bool, T.List[str]] = {
    False: [],
    True: ['-g3']
}

class TaskingCompiler(Compiler):
    '''
    Functionality that is common to all TASKING family compilers.
    '''

    LINKER_PREFIX = '-Wl'

    def __init__(self) -> None:
        if not self.is_cross:
            raise EnvironmentException(f'{id} supports only cross-compilation.')

        self.base_options = {
            OptionKey(o) for o in [
                'b_lto',
                'b_staticpic',
                'b_ndebug'
            ]
        }

        default_warn_args = [] # type: T.List[str]
        self.warn_args = {'0': [],
                          '1': default_warn_args,
                          '2': default_warn_args + [],
                          '3': default_warn_args + [],
                          'everything': default_warn_args + []} # type: T.Dict[str, T.List[str]]
        # TODO: add additional compilable files so that meson can detect it
        self.can_compile_suffixes.add('asm')

    def get_pic_args(self) -> T.List[str]:
        return ['--pic']

    def get_buildtype_args(self, buildtype: str) -> T.List[str]:
        return tasking_buildtype_args[buildtype]

    def get_debug_args(self, is_debug: bool) -> T.List[str]:
        return tasking_debug_args[is_debug]

    def get_compile_only_args(self) -> T.List[str]:
        return ['-c']

    def get_dependency_gen_args(self, outtarget: str, outfile: str) -> T.List[str]:
        return [f'--dep-file={outfile}']

    def get_depfile_suffix(self) -> str:
        return 'dep'

    def get_no_stdinc_args(self) -> T.List[str]:
        return ['--no-stdinc']

    def get_werror_args(self) -> T.List[str]:
        return ['--warnings-as-errors']

    def get_no_stdlib_link_args(self) -> T.List[str]:
        return ['--no-default-libraries']

    def get_output_args(self, outputname: str) -> T.List[str]:
        return ['-o', outputname]

    def get_include_args(self, path: str, is_system: bool) -> T.List[str]:
        if path == '':
            path = '.'
        return ['-I' + path]

    def get_optimization_args(self, optimization_level: str) -> T.List[str]:
        return tasking_optimization_args[optimization_level]

    def get_no_optimization_args(self) -> T.List[str]:
        return ['-O0']

    def get_prelink_args(self, prelink_name: str, obj_list: T.List[str]) -> T.Tuple[T.List[str], T.List[str]]:
        mil_link_list = []
        obj_file_list = []
        for obj in obj_list:
            if obj.endswith('.mil'):
                mil_link_list.append(obj)
            else:
                obj_file_list.append(obj)
        obj_file_list.append(prelink_name)

        return obj_file_list, ['--mil-link', '-o', prelink_name, '-c'] + mil_link_list

    def get_prelink_append_compile_args(self) -> bool:
        return True

    def compute_parameters_with_absolute_paths(self, parameter_list: T.List[str], build_dir: str) -> T.List[str]:
        for idx, i in enumerate(parameter_list):
            if i[:2] == '-I' or i[:2] == '-L':
                parameter_list[idx] = i[:2] + os.path.normpath(os.path.join(build_dir, i[2:]))

        return parameter_list

    def get_preprocess_only_args(self) -> T.List[str]:
        return ['-E']
