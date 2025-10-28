# SPDX-License-Identifier: Apache-2.0
# Copyright 2012-2017 The Meson development team

from __future__ import annotations

import typing as T

from ..options import OptionKey, UserStdOption

from .cpp import ALL_STDS
from .compilers import Compiler
from .mixins.apple import AppleCPPStdsMixin
from .mixins.gnu import GnuCompiler, GnuCPPStds, gnu_common_warning_args, gnu_objc_warning_args
from .mixins.clang import ClangCompiler, ClangCPPStds
from .mixins.clike import CLikeCompiler

if T.TYPE_CHECKING:
    from ..envconfig import MachineInfo
    from ..environment import Environment
    from ..linkers.linkers import DynamicLinker
    from ..mesonlib import MachineChoice
    from ..build import BuildTarget
    from ..options import MutableKeyedOptionDictType


class ObjCPPCompiler(CLikeCompiler, Compiler):

    language = 'objcpp'

    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice,
                 is_cross: bool, info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        Compiler.__init__(self, ccache, exelist, version, for_machine, info,
                          is_cross=is_cross, full_version=full_version,
                          linker=linker)
        CLikeCompiler.__init__(self)

    def form_compileropt_key(self, basename: str) -> OptionKey:
        if basename == 'std':
            return OptionKey('cpp_std', machine=self.for_machine)
        return super().form_compileropt_key(basename)

    def make_option_name(self, key: OptionKey) -> str:
        if key.name == 'std':
            return 'cpp_std'
        return super().make_option_name(key)

    @staticmethod
    def get_display_language() -> str:
        return 'Objective-C++'

    def sanity_check(self, work_dir: str, environment: 'Environment') -> None:
        code = '#import<stdio.h>\nclass MyClass;int main(void) { return 0; }\n'
        return self._sanity_check_impl(work_dir, environment, 'sanitycheckobjcpp.mm', code)

    def get_options(self) -> MutableKeyedOptionDictType:
        opts = super().get_options()
        key = self.form_compileropt_key('std')
        opts.update({
            key: UserStdOption('cpp', ALL_STDS),
        })
        return opts


class GnuObjCPPCompiler(GnuCPPStds, GnuCompiler, ObjCPPCompiler):
    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice,
                 is_cross: bool, info: 'MachineInfo',
                 defines: T.Optional[T.Dict[str, str]] = None,
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        ObjCPPCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                                info, linker=linker, full_version=full_version)
        GnuCompiler.__init__(self, defines)
        default_warn_args = ['-Wall', '-Winvalid-pch']
        self.warn_args = {'0': [],
                          '1': default_warn_args,
                          '2': default_warn_args + ['-Wextra'],
                          '3': default_warn_args + ['-Wextra', '-Wpedantic'],
                          'everything': (default_warn_args + ['-Wextra', '-Wpedantic'] +
                                         self.supported_warn_args(gnu_common_warning_args) +
                                         self.supported_warn_args(gnu_objc_warning_args))}

    def get_option_std_args(self, target: BuildTarget, env: Environment, subproject: T.Optional[str] = None) -> T.List[str]:
        args: T.List[str] = []
        key = OptionKey('cpp_std', subproject=subproject, machine=self.for_machine)
        if target:
            std = env.coredata.get_option_for_target(target, key)
        else:
            std = env.coredata.optstore.get_value_for(key)
        assert isinstance(std, str)
        if std != 'none':
            args.append('-std=' + std)
        return args

class ClangObjCPPCompiler(ClangCPPStds, ClangCompiler, ObjCPPCompiler):

    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice,
                 is_cross: bool, info: 'MachineInfo',
                 defines: T.Optional[T.Dict[str, str]] = None,
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        ObjCPPCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                                info, linker=linker, full_version=full_version)
        ClangCompiler.__init__(self, defines)
        default_warn_args = ['-Wall', '-Winvalid-pch']
        self.warn_args = {'0': [],
                          '1': default_warn_args,
                          '2': default_warn_args + ['-Wextra'],
                          '3': default_warn_args + ['-Wextra', '-Wpedantic'],
                          'everything': ['-Weverything']}

    def get_option_std_args(self, target: BuildTarget, env: Environment, subproject: T.Optional[str] = None) -> T.List[str]:
        args = []
        key = OptionKey('cpp_std', machine=self.for_machine)
        std = self.get_compileropt_value(key, env, target, subproject)
        assert isinstance(std, str)
        if std != 'none':
            args.append('-std=' + std)
        return args


class AppleClangObjCPPCompiler(AppleCPPStdsMixin, ClangObjCPPCompiler):

    """Handle the differences between Apple's clang and vanilla clang."""
