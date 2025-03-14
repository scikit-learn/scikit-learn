# SPDX-License-Identifier: Apache-2.0
# Copyright 2012-2020 The Meson development team

from __future__ import annotations

import os.path
import typing as T

from .. import options
from .. import mlog
from ..mesonlib import MesonException, version_compare
from .c_function_attributes import C_FUNC_ATTRIBUTES
from .mixins.apple import AppleCompilerMixin
from .mixins.clike import CLikeCompiler
from .mixins.ccrx import CcrxCompiler
from .mixins.xc16 import Xc16Compiler
from .mixins.compcert import CompCertCompiler
from .mixins.ti import TICompiler
from .mixins.arm import ArmCompiler, ArmclangCompiler
from .mixins.visualstudio import MSVCCompiler, ClangClCompiler
from .mixins.gnu import GnuCompiler
from .mixins.gnu import gnu_common_warning_args, gnu_c_warning_args
from .mixins.intel import IntelGnuLikeCompiler, IntelVisualStudioLikeCompiler
from .mixins.clang import ClangCompiler
from .mixins.elbrus import ElbrusCompiler
from .mixins.pgi import PGICompiler
from .mixins.emscripten import EmscriptenMixin
from .mixins.metrowerks import MetrowerksCompiler
from .mixins.metrowerks import mwccarm_instruction_set_args, mwcceppc_instruction_set_args
from .mixins.tasking import TaskingCompiler
from .compilers import (
    gnu_winlibs,
    msvc_winlibs,
    Compiler,
)

if T.TYPE_CHECKING:
    from ..coredata import MutableKeyedOptionDictType, KeyedOptionDictType
    from ..dependencies import Dependency
    from ..envconfig import MachineInfo
    from ..environment import Environment
    from ..linkers.linkers import DynamicLinker
    from ..mesonlib import MachineChoice
    from .compilers import CompileCheckMode

    CompilerMixinBase = Compiler
else:
    CompilerMixinBase = object

_ALL_STDS = ['c89', 'c9x', 'c90', 'c99', 'c1x', 'c11', 'c17', 'c18', 'c2x', 'c23']
_ALL_STDS += [f'gnu{std[1:]}' for std in _ALL_STDS]
_ALL_STDS += ['iso9899:1990', 'iso9899:199409', 'iso9899:1999', 'iso9899:2011', 'iso9899:2017', 'iso9899:2018']


class CCompiler(CLikeCompiler, Compiler):
    def attribute_check_func(self, name: str) -> str:
        try:
            return C_FUNC_ATTRIBUTES[name]
        except KeyError:
            raise MesonException(f'Unknown function attribute "{name}"')

    language = 'c'

    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice, is_cross: bool,
                 info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        # If a child ObjC or CPP class has already set it, don't set it ourselves
        Compiler.__init__(self, ccache, exelist, version, for_machine, info,
                          is_cross=is_cross, full_version=full_version, linker=linker)
        CLikeCompiler.__init__(self)

    def get_no_stdinc_args(self) -> T.List[str]:
        return ['-nostdinc']

    def sanity_check(self, work_dir: str, environment: 'Environment') -> None:
        code = 'int main(void) { int class=0; return class; }\n'
        return self._sanity_check_impl(work_dir, environment, 'sanitycheckc.c', code)

    def has_header_symbol(self, hname: str, symbol: str, prefix: str,
                          env: 'Environment', *,
                          extra_args: T.Union[None, T.List[str], T.Callable[['CompileCheckMode'], T.List[str]]] = None,
                          dependencies: T.Optional[T.List['Dependency']] = None) -> T.Tuple[bool, bool]:
        fargs = {'prefix': prefix, 'header': hname, 'symbol': symbol}
        t = '''{prefix}
        #include <{header}>
        int main(void) {{
            /* If it's not defined as a macro, try to use as a symbol */
            #ifndef {symbol}
                {symbol};
            #endif
            return 0;
        }}'''
        return self.compiles(t.format(**fargs), env, extra_args=extra_args,
                             dependencies=dependencies)

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = super().get_options()
        key = self.form_compileropt_key('std')
        opts.update({
            key: options.UserStdOption('C', _ALL_STDS),
        })
        return opts


class _ClangCStds(CompilerMixinBase):

    """Mixin class for clang based compilers for setting C standards.

    This is used by both ClangCCompiler and ClangClCompiler, as they share
    the same versions
    """

    _C17_VERSION = '>=6.0.0'
    _C18_VERSION = '>=8.0.0'
    _C2X_VERSION = '>=9.0.0'
    _C23_VERSION = '>=18.0.0'

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = super().get_options()
        stds = ['c89', 'c99', 'c11']
        # https://releases.llvm.org/6.0.0/tools/clang/docs/ReleaseNotes.html
        # https://en.wikipedia.org/wiki/Xcode#Latest_versions
        if version_compare(self.version, self._C17_VERSION):
            stds += ['c17']
        if version_compare(self.version, self._C18_VERSION):
            stds += ['c18']
        if version_compare(self.version, self._C2X_VERSION):
            stds += ['c2x']
        if version_compare(self.version, self._C23_VERSION):
            stds += ['c23']
        key = self.form_compileropt_key('std')
        std_opt = opts[key]
        assert isinstance(std_opt, options.UserStdOption), 'for mypy'
        std_opt.set_versions(stds, gnu=True)
        return opts


class ClangCCompiler(_ClangCStds, ClangCompiler, CCompiler):

    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice, is_cross: bool,
                 info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 defines: T.Optional[T.Dict[str, str]] = None,
                 full_version: T.Optional[str] = None):
        CCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross, info, linker=linker, full_version=full_version)
        ClangCompiler.__init__(self, defines)
        default_warn_args = ['-Wall', '-Winvalid-pch']
        self.warn_args = {'0': [],
                          '1': default_warn_args,
                          '2': default_warn_args + ['-Wextra'],
                          '3': default_warn_args + ['-Wextra', '-Wpedantic'],
                          'everything': ['-Weverything']}

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = super().get_options()
        if self.info.is_windows() or self.info.is_cygwin():
            self.update_options(
                opts,
                self.create_option(options.UserArrayOption,
                                   self.form_compileropt_key('winlibs'),
                                   'Standard Windows libs to link against',
                                   gnu_winlibs),
            )
        return opts

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        args = []
        key = self.form_compileropt_key('std')
        std = options.get_value(key)
        if std != 'none':
            args.append('-std=' + std)
        return args

    def get_option_link_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        if self.info.is_windows() or self.info.is_cygwin():
            # without a typedict mypy can't understand this.
            key = self.form_compileropt_key('winlibs')
            libs = options.get_value(key).copy()
            assert isinstance(libs, list)
            for l in libs:
                assert isinstance(l, str)
            return libs
        return []


class ArmLtdClangCCompiler(ClangCCompiler):

    id = 'armltdclang'


class AppleClangCCompiler(AppleCompilerMixin, ClangCCompiler):

    """Handle the differences between Apple Clang and Vanilla Clang.

    Right now this just handles the differences between the versions that new
    C standards were added.
    """

    _C17_VERSION = '>=10.0.0'
    _C18_VERSION = '>=11.0.0'
    _C2X_VERSION = '>=11.0.0'


class EmscriptenCCompiler(EmscriptenMixin, ClangCCompiler):

    id = 'emscripten'

    # Emscripten uses different version numbers than Clang; `emcc -v` will show
    # the Clang version number used as well (but `emcc --version` does not).
    # See https://github.com/pyodide/pyodide/discussions/4762 for more on
    # emcc <--> clang versions. Note that c17/c18/c2x are always available, since
    # the lowest supported Emscripten version used a new-enough Clang version.
    _C17_VERSION = '>=1.38.35'
    _C18_VERSION = '>=1.38.35'
    _C2X_VERSION = '>=1.38.35'  # 1.38.35 used Clang 9.0.0
    _C23_VERSION = '>=3.1.45'    # 3.1.45 used Clang 18.0.0

    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice, is_cross: bool,
                 info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 defines: T.Optional[T.Dict[str, str]] = None,
                 full_version: T.Optional[str] = None):
        if not is_cross:
            raise MesonException('Emscripten compiler can only be used for cross compilation.')
        if not version_compare(version, '>=1.39.19'):
            raise MesonException('Meson requires Emscripten >= 1.39.19')
        ClangCCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                                info, linker=linker,
                                defines=defines, full_version=full_version)


class ArmclangCCompiler(ArmclangCompiler, CCompiler):
    '''
    Keil armclang
    '''

    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice, is_cross: bool,
                 info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                           info, linker=linker, full_version=full_version)
        ArmclangCompiler.__init__(self)
        default_warn_args = ['-Wall', '-Winvalid-pch']
        self.warn_args = {'0': [],
                          '1': default_warn_args,
                          '2': default_warn_args + ['-Wextra'],
                          '3': default_warn_args + ['-Wextra', '-Wpedantic'],
                          'everything': ['-Weverything']}

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = CCompiler.get_options(self)
        key = self.form_compileropt_key('std')
        std_opt = opts[key]
        assert isinstance(std_opt, options.UserStdOption), 'for mypy'
        std_opt.set_versions(['c90', 'c99', 'c11'], gnu=True)
        return opts

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        args = []
        key = self.form_compileropt_key('std')
        std = options.get_value(key)
        if std != 'none':
            args.append('-std=' + std)
        return args

    def get_option_link_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        return []


class GnuCCompiler(GnuCompiler, CCompiler):

    _C18_VERSION = '>=8.0.0'
    _C2X_VERSION = '>=9.0.0'
    _C23_VERSION = '>=14.0.0'
    _INVALID_PCH_VERSION = ">=3.4.0"

    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice, is_cross: bool,
                 info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 defines: T.Optional[T.Dict[str, str]] = None,
                 full_version: T.Optional[str] = None):
        CCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross, info, linker=linker, full_version=full_version)
        GnuCompiler.__init__(self, defines)
        default_warn_args = ['-Wall']
        if version_compare(self.version, self._INVALID_PCH_VERSION):
            default_warn_args += ['-Winvalid-pch']
        self.warn_args = {'0': [],
                          '1': default_warn_args,
                          '2': default_warn_args + ['-Wextra'],
                          '3': default_warn_args + ['-Wextra', '-Wpedantic'],
                          'everything': (default_warn_args + ['-Wextra', '-Wpedantic'] +
                                         self.supported_warn_args(gnu_common_warning_args) +
                                         self.supported_warn_args(gnu_c_warning_args))}

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = CCompiler.get_options(self)
        stds = ['c89', 'c99', 'c11']
        if version_compare(self.version, self._C18_VERSION):
            stds += ['c17', 'c18']
        if version_compare(self.version, self._C2X_VERSION):
            stds += ['c2x']
        if version_compare(self.version, self._C23_VERSION):
            stds += ['c23']
        key = self.form_compileropt_key('std')
        std_opt = opts[key]
        assert isinstance(std_opt, options.UserStdOption), 'for mypy'
        std_opt.set_versions(stds, gnu=True)
        if self.info.is_windows() or self.info.is_cygwin():
            self.update_options(
                opts,
                self.create_option(options.UserArrayOption,
                                   key.evolve('c_winlibs'),
                                   'Standard Windows libs to link against',
                                   gnu_winlibs),
            )
        return opts

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        args = []
        key = self.form_compileropt_key('std')
        std = options.get_value(key)
        if std != 'none':
            args.append('-std=' + std)
        return args

    def get_option_link_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        if self.info.is_windows() or self.info.is_cygwin():
            # without a typeddict mypy can't figure this out
            key = self.form_compileropt_key('winlibs')
            libs: T.List[str] = options.get_value(key).copy()
            assert isinstance(libs, list)
            for l in libs:
                assert isinstance(l, str)
            return libs
        return []

    def get_pch_use_args(self, pch_dir: str, header: str) -> T.List[str]:
        return ['-fpch-preprocess', '-include', os.path.basename(header)]


class PGICCompiler(PGICompiler, CCompiler):
    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice, is_cross: bool,
                 info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                           info, linker=linker, full_version=full_version)
        PGICompiler.__init__(self)


class NvidiaHPC_CCompiler(PGICompiler, CCompiler):

    id = 'nvidia_hpc'

    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice, is_cross: bool,
                 info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                           info, linker=linker, full_version=full_version)
        PGICompiler.__init__(self)

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = CCompiler.get_options(self)
        cppstd_choices = ['c89', 'c90', 'c99', 'c11', 'c17', 'c18']
        std_opt = opts[self.form_compileropt_key('std')]
        assert isinstance(std_opt, options.UserStdOption), 'for mypy'
        std_opt.set_versions(cppstd_choices, gnu=True)
        return opts


class ElbrusCCompiler(ElbrusCompiler, CCompiler):
    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice, is_cross: bool,
                 info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 defines: T.Optional[T.Dict[str, str]] = None,
                 full_version: T.Optional[str] = None):
        CCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                           info, linker=linker, full_version=full_version)
        ElbrusCompiler.__init__(self)

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = CCompiler.get_options(self)
        stds = ['c89', 'c9x', 'c99', 'gnu89', 'gnu9x', 'gnu99']
        stds += ['iso9899:1990', 'iso9899:199409', 'iso9899:1999']
        if version_compare(self.version, '>=1.20.00'):
            stds += ['c11', 'gnu11']
        if version_compare(self.version, '>=1.21.00') and version_compare(self.version, '<1.22.00'):
            stds += ['c90', 'c1x', 'gnu90', 'gnu1x', 'iso9899:2011']
        if version_compare(self.version, '>=1.23.00'):
            stds += ['c90', 'c1x', 'gnu90', 'gnu1x', 'iso9899:2011']
        if version_compare(self.version, '>=1.26.00'):
            stds += ['c17', 'c18', 'iso9899:2017', 'iso9899:2018', 'gnu17', 'gnu18']
        key = self.form_compileropt_key('std')
        std_opt = opts[key]
        assert isinstance(std_opt, options.UserStdOption), 'for mypy'
        std_opt.set_versions(stds)
        return opts

    # Elbrus C compiler does not have lchmod, but there is only linker warning, not compiler error.
    # So we should explicitly fail at this case.
    def has_function(self, funcname: str, prefix: str, env: 'Environment', *,
                     extra_args: T.Optional[T.List[str]] = None,
                     dependencies: T.Optional[T.List['Dependency']] = None) -> T.Tuple[bool, bool]:
        if funcname == 'lchmod':
            return False, False
        else:
            return super().has_function(funcname, prefix, env,
                                        extra_args=extra_args,
                                        dependencies=dependencies)


class IntelCCompiler(IntelGnuLikeCompiler, CCompiler):
    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice, is_cross: bool,
                 info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                           info, linker=linker, full_version=full_version)
        IntelGnuLikeCompiler.__init__(self)
        self.lang_header = 'c-header'
        default_warn_args = ['-Wall', '-w3']
        self.warn_args = {'0': [],
                          '1': default_warn_args + ['-diag-disable:remark'],
                          '2': default_warn_args + ['-Wextra', '-diag-disable:remark'],
                          '3': default_warn_args + ['-Wextra', '-diag-disable:remark'],
                          'everything': default_warn_args + ['-Wextra']}

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = CCompiler.get_options(self)
        stds = ['c89', 'c99']
        if version_compare(self.version, '>=16.0.0'):
            stds += ['c11']
        key = self.form_compileropt_key('std')
        std_opt = opts[key]
        assert isinstance(std_opt, options.UserStdOption), 'for mypy'
        std_opt.set_versions(stds, gnu=True)
        return opts

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        args = []
        key = self.form_compileropt_key('std')
        std = options.get_value(key)
        if std != 'none':
            args.append('-std=' + std)
        return args


class IntelLLVMCCompiler(ClangCCompiler):

    id = 'intel-llvm'


class VisualStudioLikeCCompilerMixin(CompilerMixinBase):

    """Shared methods that apply to MSVC-like C compilers."""

    def get_options(self) -> MutableKeyedOptionDictType:
        return self.update_options(
            super().get_options(),
            self.create_option(
                options.UserArrayOption,
                self.form_compileropt_key('winlibs'),
                'Standard Windows libs to link against',
                msvc_winlibs,
            ),
        )

    def get_option_link_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        # need a TypeDict to make this work
        key = self.form_compileropt_key('winlibs')
        libs = options.get_value(key).copy()
        assert isinstance(libs, list)
        for l in libs:
            assert isinstance(l, str)
        return libs


class VisualStudioCCompiler(MSVCCompiler, VisualStudioLikeCCompilerMixin, CCompiler):

    _C11_VERSION = '>=19.28'
    _C17_VERSION = '>=19.28'

    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice,
                 is_cross: bool, info: 'MachineInfo', target: str,
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                           info, linker=linker,
                           full_version=full_version)
        MSVCCompiler.__init__(self, target)

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = super().get_options()
        stds = ['c89', 'c99']
        if version_compare(self.version, self._C11_VERSION):
            stds += ['c11']
        if version_compare(self.version, self._C17_VERSION):
            stds += ['c17', 'c18']
        key = self.form_compileropt_key('std')
        std_opt = opts[key]
        assert isinstance(std_opt, options.UserStdOption), 'for mypy'
        std_opt.set_versions(stds, gnu=True, gnu_deprecated=True)
        return opts

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        args = []
        key = self.form_compileropt_key('std')
        std = options.get_value(key)
        # As of MVSC 16.8, /std:c11 and /std:c17 are the only valid C standard options.
        if std == 'c11':
            args.append('/std:c11')
        elif std in {'c17', 'c18'}:
            args.append('/std:c17')
        return args


class ClangClCCompiler(_ClangCStds, ClangClCompiler, VisualStudioLikeCCompilerMixin, CCompiler):
    def __init__(self, exelist: T.List[str], version: str, for_machine: MachineChoice,
                 is_cross: bool, info: 'MachineInfo', target: str,
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CCompiler.__init__(self, [], exelist, version, for_machine, is_cross,
                           info, linker=linker,
                           full_version=full_version)
        ClangClCompiler.__init__(self, target)

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        key = self.form_compileropt_key('std')
        std = options.get_value(key)
        if std != "none":
            return [f'/clang:-std={std}']
        return []


class IntelClCCompiler(IntelVisualStudioLikeCompiler, VisualStudioLikeCCompilerMixin, CCompiler):

    """Intel "ICL" compiler abstraction."""

    def __init__(self, exelist: T.List[str], version: str, for_machine: MachineChoice,
                 is_cross: bool, info: 'MachineInfo', target: str,
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CCompiler.__init__(self, [], exelist, version, for_machine, is_cross,
                           info, linker=linker,
                           full_version=full_version)
        IntelVisualStudioLikeCompiler.__init__(self, target)

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = super().get_options()
        key = self.form_compileropt_key('std')
        # To shut up mypy.
        if isinstance(opts, dict):
            raise RuntimeError('This is a transitory issue that should not happen. Please report with full backtrace.')
        std_opt = opts.get_value_object(key)
        assert isinstance(std_opt, options.UserStdOption), 'for mypy'
        std_opt.set_versions(['c89', 'c99', 'c11'])
        return opts

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        args = []
        key = self.form_compileropt_key('std')
        std = options.get_value(key)
        if std == 'c89':
            mlog.log("ICL doesn't explicitly implement c89, setting the standard to 'none', which is close.", once=True)
        elif std != 'none':
            args.append('/Qstd:' + std)
        return args


class IntelLLVMClCCompiler(IntelClCCompiler):

    id = 'intel-llvm-cl'


class ArmCCompiler(ArmCompiler, CCompiler):
    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice,
                 is_cross: bool, info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                           info, linker=linker,
                           full_version=full_version)
        ArmCompiler.__init__(self)

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = CCompiler.get_options(self)
        key = self.form_compileropt_key('std')
        std_opt = opts[key]
        assert isinstance(std_opt, options.UserStdOption), 'for mypy'
        std_opt.set_versions(['c89', 'c99', 'c11'])
        return opts

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        args = []
        key = self.form_compileropt_key('std')
        std = options.get_value(key)
        if std != 'none':
            args.append('--' + std)
        return args


class CcrxCCompiler(CcrxCompiler, CCompiler):
    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice,
                 is_cross: bool, info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                           info, linker=linker, full_version=full_version)
        CcrxCompiler.__init__(self)

    # Override CCompiler.get_always_args
    def get_always_args(self) -> T.List[str]:
        return ['-nologo']

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = CCompiler.get_options(self)
        key = self.form_compileropt_key('std')
        std_opt = opts[key]
        assert isinstance(std_opt, options.UserStdOption), 'for mypy'
        std_opt.set_versions(['c89', 'c99'])
        return opts

    def get_no_stdinc_args(self) -> T.List[str]:
        return []

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        args = []
        key = self.form_compileropt_key('std')
        std = options.get_value(key)
        if std == 'c89':
            args.append('-lang=c')
        elif std == 'c99':
            args.append('-lang=c99')
        return args

    def get_compile_only_args(self) -> T.List[str]:
        return []

    def get_no_optimization_args(self) -> T.List[str]:
        return ['-optimize=0']

    def get_output_args(self, target: str) -> T.List[str]:
        return [f'-output=obj={target}']

    def get_werror_args(self) -> T.List[str]:
        return ['-change_message=error']

    def get_include_args(self, path: str, is_system: bool) -> T.List[str]:
        if path == '':
            path = '.'
        return ['-include=' + path]


class Xc16CCompiler(Xc16Compiler, CCompiler):
    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice,
                 is_cross: bool, info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                           info, linker=linker, full_version=full_version)
        Xc16Compiler.__init__(self)

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = CCompiler.get_options(self)
        key = self.form_compileropt_key('std')
        std_opt = opts[key]
        assert isinstance(std_opt, options.UserStdOption), 'for mypy'
        std_opt.set_versions(['c89', 'c99'], gnu=True)
        return opts

    def get_no_stdinc_args(self) -> T.List[str]:
        return []

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        args = []
        key = self.form_compileropt_key('std')
        std = options.get_value(key)
        if std != 'none':
            args.append('-ansi')
            args.append('-std=' + std)
        return args

    def get_compile_only_args(self) -> T.List[str]:
        return []

    def get_no_optimization_args(self) -> T.List[str]:
        return ['-O0']

    def get_output_args(self, target: str) -> T.List[str]:
        return [f'-o{target}']

    def get_werror_args(self) -> T.List[str]:
        return ['-change_message=error']

    def get_include_args(self, path: str, is_system: bool) -> T.List[str]:
        if path == '':
            path = '.'
        return ['-I' + path]

class CompCertCCompiler(CompCertCompiler, CCompiler):
    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice,
                 is_cross: bool, info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                           info, linker=linker, full_version=full_version)
        CompCertCompiler.__init__(self)

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = CCompiler.get_options(self)
        key = self.form_compileropt_key('std')
        std_opt = opts[key]
        assert isinstance(std_opt, options.UserStdOption), 'for mypy'
        std_opt.set_versions(['c89', 'c99'])
        return opts

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        return []

    def get_no_optimization_args(self) -> T.List[str]:
        return ['-O0']

    def get_output_args(self, target: str) -> T.List[str]:
        return [f'-o{target}']

    def get_werror_args(self) -> T.List[str]:
        return ['-Werror']

    def get_include_args(self, path: str, is_system: bool) -> T.List[str]:
        if path == '':
            path = '.'
        return ['-I' + path]

class TICCompiler(TICompiler, CCompiler):
    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice,
                 is_cross: bool, info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                           info, linker=linker, full_version=full_version)
        TICompiler.__init__(self)

    # Override CCompiler.get_always_args
    def get_always_args(self) -> T.List[str]:
        return []

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = CCompiler.get_options(self)
        key = self.form_compileropt_key('std')
        std_opt = opts[key]
        assert isinstance(std_opt, options.UserStdOption), 'for mypy'
        std_opt.set_versions(['c89', 'c99', 'c11'])
        return opts

    def get_no_stdinc_args(self) -> T.List[str]:
        return []

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        args = []
        key = self.form_compileropt_key('std')
        std = options.get_value(key)
        if std != 'none':
            args.append('--' + std)
        return args

class C2000CCompiler(TICCompiler):
    # Required for backwards compat with projects created before ti-cgt support existed
    id = 'c2000'

class C6000CCompiler(TICCompiler):
    id = 'c6000'

class MetrowerksCCompilerARM(MetrowerksCompiler, CCompiler):
    id = 'mwccarm'

    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice,
                 is_cross: bool, info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                           info, linker=linker, full_version=full_version)
        MetrowerksCompiler.__init__(self)

    def get_instruction_set_args(self, instruction_set: str) -> T.Optional[T.List[str]]:
        return mwccarm_instruction_set_args.get(instruction_set, None)

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = CCompiler.get_options(self)
        c_stds = ['c99']
        key = self.form_compileropt_key('std')
        opts[key].choices = ['none'] + c_stds
        return opts

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        args = []
        key = self.form_compileropt_key('std')
        std = options.get_value(key)
        if std != 'none':
            args.append('-lang')
            args.append(std)
        return args

class MetrowerksCCompilerEmbeddedPowerPC(MetrowerksCompiler, CCompiler):
    id = 'mwcceppc'

    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice,
                 is_cross: bool, info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                           info, linker=linker, full_version=full_version)
        MetrowerksCompiler.__init__(self)

    def get_instruction_set_args(self, instruction_set: str) -> T.Optional[T.List[str]]:
        return mwcceppc_instruction_set_args.get(instruction_set, None)

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = CCompiler.get_options(self)
        c_stds = ['c99']
        key = self.form_compileropt_key('std')
        opts[key].choices = ['none'] + c_stds
        return opts

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        args = []
        key = self.form_compileropt_key('std')
        std = options.get_value(key)
        if std != 'none':
            args.append('-lang ' + std)
        return args

class TaskingCCompiler(TaskingCompiler, CCompiler):
    id = 'tasking'

    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice,
                 is_cross: bool, info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                           info, linker=linker, full_version=full_version)
        TaskingCompiler.__init__(self)
