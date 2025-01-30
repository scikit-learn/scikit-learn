# SPDX-License-Identifier: Apache-2.0
# Copyright © 2024 Intel Corporation

"""Provides mixins for Apple compilers."""

from __future__ import annotations
import typing as T

from ...mesonlib import MesonException

if T.TYPE_CHECKING:
    from ..._typing import ImmutableListProtocol
    from ...environment import Environment
    from ..compilers import Compiler
else:
    # This is a bit clever, for mypy we pretend that these mixins descend from
    # Compiler, so we get all of the methods and attributes defined for us, but
    # for runtime we make them descend from object (which all classes normally
    # do). This gives up DRYer type checking, with no runtime impact
    Compiler = object


class AppleCompilerMixin(Compiler):

    """Handle differences between Vanilla Clang and the Clang shipped with XCode."""

    __BASE_OMP_FLAGS: ImmutableListProtocol[str] = ['-Xpreprocessor', '-fopenmp']

    def openmp_flags(self, env: Environment) -> T.List[str]:
        """Flags required to compile with OpenMP on Apple.

        The Apple Clang Compiler doesn't have builtin support for OpenMP, it
        must be provided separately. As such, we need to add the -Xpreprocessor
        argument so that an external OpenMP can be found.

        :return: A list of arguments
        """
        m = env.machines[self.for_machine]
        assert m is not None, 'for mypy'
        if m.cpu_family.startswith('x86'):
            root = '/usr/local'
        else:
            root = '/opt/homebrew'
        return self.__BASE_OMP_FLAGS + [f'-I{root}/opt/libomp/include']

    def openmp_link_flags(self, env: Environment) -> T.List[str]:
        m = env.machines[self.for_machine]
        assert m is not None, 'for mypy'
        if m.cpu_family.startswith('x86'):
            root = '/usr/local'
        else:
            root = '/opt/homebrew'

        link = self.find_library('omp', env, [f'{root}/opt/libomp/lib'])
        if not link:
            raise MesonException("Couldn't find libomp")
        return self.__BASE_OMP_FLAGS + link

    def get_prelink_args(self, prelink_name: str, obj_list: T.List[str]) -> T.List[str]:
        # The objects are prelinked through the compiler, which injects -lSystem
        return ['-nostdlib', '-r', '-o', prelink_name] + obj_list
