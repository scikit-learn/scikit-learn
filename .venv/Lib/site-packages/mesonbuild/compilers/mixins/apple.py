# SPDX-License-Identifier: Apache-2.0
# Copyright © 2024-2025 Intel Corporation

"""Provides mixins for Apple compilers."""

from __future__ import annotations
import functools
import subprocess
import typing as T

from ...mesonlib import MesonException


@functools.lru_cache(maxsize=None)
def _get_libomp_prefix() -> T.Optional[str]:
    """Call `brew --prefix libomp` once and cache it. Returns None if unavailable."""
    try:
        return subprocess.run(
            ['brew', '--prefix', 'libomp'],
            capture_output=True,
            encoding='utf-8',
            check=True,
        ).stdout.strip()
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None


def _get_homebrew_libomp_root(cpu_family: str, is_cross: bool) -> str:
    """Return the libomp root, preferring dynamic detection with arch-based fallback."""
    if not is_cross:
        libomp_prefix = _get_libomp_prefix()
        if libomp_prefix is not None:
            return libomp_prefix
    # Fallback: brew not on PATH, use historical defaults based on architecture
    if cpu_family.startswith('x86'):
        return '/usr/local/opt/libomp'
    return '/opt/homebrew/opt/libomp'


if T.TYPE_CHECKING:
    from ..._typing import ImmutableListProtocol
    from ...envconfig import MachineInfo
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

    if T.TYPE_CHECKING:
        # Older versions of mypy can't figure this out
        info: MachineInfo

    def openmp_flags(self) -> T.List[str]:
        """Flags required to compile with OpenMP on Apple.

        The Apple Clang Compiler doesn't have builtin support for OpenMP, it
        must be provided separately. As such, we need to add the -Xpreprocessor
        argument so that an external OpenMP can be found.

        :return: A list of arguments
        """
        root = _get_homebrew_libomp_root(self.info.cpu_family, self.is_cross)
        return self.__BASE_OMP_FLAGS + [f'-I{root}/include']

    def openmp_link_flags(self) -> T.List[str]:
        root = _get_homebrew_libomp_root(self.info.cpu_family, self.is_cross)
        link = self.find_library('omp', [f'{root}/lib'])
        if not link:
            raise MesonException("Couldn't find libomp")
        return self.__BASE_OMP_FLAGS + link

    def get_prelink_args(self, prelink_name: str, obj_list: T.List[str]) -> T.Tuple[T.List[str], T.List[str]]:
        # The objects are prelinked through the compiler, which injects -lSystem
        return [prelink_name], ['-nostdlib', '-r', '-o', prelink_name] + obj_list


class AppleCStdsMixin(Compiler):

    """Provide version overrides for the Apple Compilers."""

    _C17_VERSION = '>=10.0.0'
    _C18_VERSION = '>=11.0.0'
    _C2X_VERSION = '>=11.0.3'
    _C23_VERSION = '>=17.0.0'
    _C2Y_VERSION = '>=17.0.0'


class AppleCPPStdsMixin(Compiler):

    """Provide version overrides for the Apple C++ Compilers."""

    _CPP23_VERSION = '>=13.0.0'
    _CPP26_VERSION = '>=16.0.0'
