# SPDX-License-Identifier: Apache-2.0
# Copyright 2012-2019 The Meson development team

from __future__ import annotations

"""Representations specific to the Microchip XC C/C++ compiler family."""

import os
import typing as T

from .gnu import GnuCStds, GnuCPPStds
from ..compilers import Compiler
from ...mesonlib import EnvironmentException, version_compare

if T.TYPE_CHECKING:
    from ...build import BuildTarget
    from ...envconfig import MachineInfo

    CompilerBase = Compiler
else:
    # This is a bit clever, for mypy we pretend that these mixins descend from
    # Compiler, so we get all of the methods and attributes defined for us, but
    # for runtime we make them descend from object (which all classes normally
    # do). This gives up DRYer type checking, with no runtime impact
    CompilerBase = object

xc16_optimization_args: T.Dict[str, T.List[str]] = {
    'plain': [],
    '0': ['-O0'],
    'g': ['-O0'],
    '1': ['-O1'],
    '2': ['-O2'],
    '3': ['-O3'],
    's': ['-Os']
}

xc16_debug_args: T.Dict[bool, T.List[str]] = {
    False: [],
    True: []
}


class Xc16Compiler(Compiler):

    id = 'xc16'

    def __init__(self) -> None:
        if not self.is_cross:
            raise EnvironmentException('xc16 supports only cross-compilation.')
        # Assembly
        self.can_compile_suffixes.add('s')
        self.can_compile_suffixes.add('sx')
        default_warn_args: T.List[str] = []
        self.warn_args = {'0': [],
                          '1': default_warn_args,
                          '2': default_warn_args + [],
                          '3': default_warn_args + [],
                          'everything': default_warn_args + []}

    def get_always_args(self) -> T.List[str]:
        return []

    def get_pic_args(self) -> T.List[str]:
        # PIC support is not enabled by default for xc16,
        # if users want to use it, they need to add the required arguments explicitly
        return []

    def get_pch_suffix(self) -> str:
        return 'pch'

    def get_pch_use_args(self, pch_dir: str, header: str) -> T.List[str]:
        return []

    def thread_flags(self) -> T.List[str]:
        return []

    def get_coverage_args(self) -> T.List[str]:
        return []

    def get_no_stdinc_args(self) -> T.List[str]:
        return ['-nostdinc']

    def get_no_stdlib_link_args(self) -> T.List[str]:
        return ['--nostdlib']

    def get_optimization_args(self, optimization_level: str) -> T.List[str]:
        return xc16_optimization_args[optimization_level]

    def get_debug_args(self, is_debug: bool) -> T.List[str]:
        return xc16_debug_args[is_debug]

    @classmethod
    def _unix_args_to_native(cls, args: T.List[str], info: MachineInfo) -> T.List[str]:
        result = []
        for i in args:
            if i.startswith('-D'):
                i = '-D' + i[2:]
            if i.startswith('-I'):
                i = '-I' + i[2:]
            if i.startswith('-Wl,-rpath='):
                continue
            elif i == '--print-search-dirs':
                continue
            elif i.startswith('-L'):
                continue
            result.append(i)
        return result

    def compute_parameters_with_absolute_paths(self, parameter_list: T.List[str], build_dir: str) -> T.List[str]:
        for idx, i in enumerate(parameter_list):
            if i[:9] == '-I':
                parameter_list[idx] = i[:9] + os.path.normpath(os.path.join(build_dir, i[9:]))

        return parameter_list


class Xc32Compiler(CompilerBase):

    """Microchip XC32 compiler mixin. GCC based with some options disabled."""

    id = 'xc32-gcc'

    gcc_version = '4.5.1'  # Defaults to GCC version used by first XC32 release (v1.00).

    _COLOR_VERSION = '>=3.0'       # XC32 version based on GCC 8.3.1+
    _WPEDANTIC_VERSION = '>=1.40'  # XC32 version based on GCC 4.8.3+
    _LTO_AUTO_VERSION = '>=5.00'   # XC32 version based on GCC 13.2.1+
    _LTO_CACHE_VERSION = '==-1'
    _USE_MOLD_VERSION = '==-1'

    def __init__(self) -> None:
        if not self.is_cross:
            raise EnvironmentException('XC32 supports only cross-compilation.')

    def get_instruction_set_args(self, instruction_set: str) -> T.Optional[T.List[str]]:
        return None

    def thread_flags(self) -> T.List[str]:
        return []

    def openmp_flags(self) -> T.List[str]:
        return Compiler.openmp_flags(self)

    def get_pic_args(self) -> T.List[str]:
        return Compiler.get_pic_args(self)

    def get_pie_args(self) -> T.List[str]:
        return Compiler.get_pie_args(self)

    def get_profile_generate_args(self) -> T.List[str]:
        return Compiler.get_profile_generate_args(self)

    def get_profile_use_args(self) -> T.List[str]:
        return Compiler.get_profile_use_args(self)

    def sanitizer_compile_args(self, target: T.Optional[BuildTarget], value: T.List[str]) -> T.List[str]:
        return []

    @classmethod
    def use_linker_args(cls, linker: str, version: str) -> T.List[str]:
        return []

    def get_coverage_args(self) -> T.List[str]:
        return []

    def get_largefile_args(self) -> T.List[str]:
        return []

    def get_prelink_args(self, prelink_name: str, obj_list: T.List[str]) -> T.Tuple[T.List[str], T.List[str]]:
        return Compiler.get_prelink_args(self, prelink_name, obj_list)

    def get_prelink_append_compile_args(self) -> bool:
        return False

    def supported_warn_args(self, warn_args_by_version: T.Dict[str, T.List[str]]) -> T.List[str]:
        result: T.List[str] = []
        for version, warn_args in warn_args_by_version.items():
            if version_compare(self.gcc_version, '>=' + version):
                result += warn_args
        return result

class Xc32CStds(GnuCStds):

    """Mixin for setting C standards based on XC32 version."""

    _C18_VERSION = '>=3.0'
    _C2X_VERSION = '>=5.00'
    _C23_VERSION = '==-1'
    _C2Y_VERSION = '==-1'

class Xc32CPPStds(GnuCPPStds):

    """Mixin for setting C++ standards based on XC32 version."""

    _CPP23_VERSION = '>=5.00'
    _CPP26_VERSION = '==-1'
