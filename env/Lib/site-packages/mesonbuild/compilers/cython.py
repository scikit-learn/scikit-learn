# SPDX-License-Identifier: Apache-2.0
# Copyright © 2021 Intel Corporation
from __future__ import annotations

"""Abstraction for Cython language compilers."""

import typing as T

from .. import options
from ..mesonlib import EnvironmentException, version_compare
from .compilers import Compiler

if T.TYPE_CHECKING:
    from ..coredata import MutableKeyedOptionDictType, KeyedOptionDictType
    from ..environment import Environment


class CythonCompiler(Compiler):

    """Cython Compiler."""

    language = 'cython'
    id = 'cython'

    def needs_static_linker(self) -> bool:
        # We transpile into C, so we don't need any linker
        return False

    def get_always_args(self) -> T.List[str]:
        return ['--fast-fail']

    def get_werror_args(self) -> T.List[str]:
        return ['-Werror']

    def get_output_args(self, outputname: str) -> T.List[str]:
        return ['-o', outputname]

    def get_optimization_args(self, optimization_level: str) -> T.List[str]:
        # Cython doesn't have optimization levels itself, the underlying
        # compiler might though
        return []

    def get_dependency_gen_args(self, outtarget: str, outfile: str) -> T.List[str]:
        if version_compare(self.version, '>=0.29.33'):
            return ['-M']
        return []

    def get_depfile_suffix(self) -> str:
        return 'dep'

    def sanity_check(self, work_dir: str, environment: 'Environment') -> None:
        code = 'print("hello world")'
        with self.cached_compile(code, environment.coredata) as p:
            if p.returncode != 0:
                raise EnvironmentException(f'Cython compiler {self.id!r} cannot compile programs')

    def get_pic_args(self) -> T.List[str]:
        # We can lie here, it's fine
        return []

    def compute_parameters_with_absolute_paths(self, parameter_list: T.List[str],
                                               build_dir: str) -> T.List[str]:
        new: T.List[str] = []
        for i in parameter_list:
            new.append(i)

        return new

    def get_options(self) -> 'MutableKeyedOptionDictType':
        return self.update_options(
            super().get_options(),
            self.create_option(options.UserComboOption,
                               self.form_compileropt_key('version'),
                               'Python version to target',
                               ['2', '3'],
                               '3'),
            self.create_option(options.UserComboOption,
                               self.form_compileropt_key('language'),
                               'Output C or C++ files',
                               ['c', 'cpp'],
                               'c'),
        )

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        args: T.List[str] = []
        key = self.form_compileropt_key('version')
        version = options.get_value(key)
        args.append(f'-{version}')
        key = self.form_compileropt_key('language')
        lang = options.get_value(key)
        if lang == 'cpp':
            args.append('--cplus')
        return args
