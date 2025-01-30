# SPDX-License-Identifier: Apache-2.0
# Copyright 2012-2017 The Meson development team

from __future__ import annotations

import copy
import functools
import os.path
import typing as T

from .. import options
from .. import mlog
from ..mesonlib import MesonException, version_compare

from .compilers import (
    gnu_winlibs,
    msvc_winlibs,
    Compiler,
    CompileCheckMode,
)
from .c_function_attributes import CXX_FUNC_ATTRIBUTES, C_FUNC_ATTRIBUTES
from .mixins.apple import AppleCompilerMixin
from .mixins.clike import CLikeCompiler
from .mixins.ccrx import CcrxCompiler
from .mixins.ti import TICompiler
from .mixins.arm import ArmCompiler, ArmclangCompiler
from .mixins.visualstudio import MSVCCompiler, ClangClCompiler
from .mixins.gnu import GnuCompiler, gnu_common_warning_args, gnu_cpp_warning_args
from .mixins.intel import IntelGnuLikeCompiler, IntelVisualStudioLikeCompiler
from .mixins.clang import ClangCompiler
from .mixins.elbrus import ElbrusCompiler
from .mixins.pgi import PGICompiler
from .mixins.emscripten import EmscriptenMixin
from .mixins.metrowerks import MetrowerksCompiler
from .mixins.metrowerks import mwccarm_instruction_set_args, mwcceppc_instruction_set_args

if T.TYPE_CHECKING:
    from ..coredata import MutableKeyedOptionDictType, KeyedOptionDictType
    from ..dependencies import Dependency
    from ..envconfig import MachineInfo
    from ..environment import Environment
    from ..linkers.linkers import DynamicLinker
    from ..mesonlib import MachineChoice
    CompilerMixinBase = CLikeCompiler
else:
    CompilerMixinBase = object

_ALL_STDS = ['c++98', 'c++0x', 'c++03', 'c++1y', 'c++1z', 'c++11', 'c++14', 'c++17', 'c++2a', 'c++20', 'c++23', 'c++26']
_ALL_STDS += [f'gnu{std[1:]}' for std in _ALL_STDS]
_ALL_STDS += ['vc++11', 'vc++14', 'vc++17', 'vc++20', 'vc++latest', 'c++latest']


def non_msvc_eh_options(eh: str, args: T.List[str]) -> None:
    if eh == 'none':
        args.append('-fno-exceptions')
    elif eh in {'s', 'c'}:
        mlog.warning(f'non-MSVC compilers do not support {eh} exception handling. '
                     'You may want to set eh to \'default\'.', fatal=False)

class CPPCompiler(CLikeCompiler, Compiler):
    def attribute_check_func(self, name: str) -> str:
        try:
            return CXX_FUNC_ATTRIBUTES.get(name, C_FUNC_ATTRIBUTES[name])
        except KeyError:
            raise MesonException(f'Unknown function attribute "{name}"')

    language = 'cpp'

    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice, is_cross: bool,
                 info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        # If a child ObjCPP class has already set it, don't set it ourselves
        Compiler.__init__(self, ccache, exelist, version, for_machine, info,
                          is_cross=is_cross, linker=linker,
                          full_version=full_version)
        CLikeCompiler.__init__(self)

    @classmethod
    def get_display_language(cls) -> str:
        return 'C++'

    def get_no_stdinc_args(self) -> T.List[str]:
        return ['-nostdinc++']

    def get_no_stdlib_link_args(self) -> T.List[str]:
        return ['-nostdlib++']

    def sanity_check(self, work_dir: str, environment: 'Environment') -> None:
        code = 'class breakCCompiler;int main(void) { return 0; }\n'
        return self._sanity_check_impl(work_dir, environment, 'sanitycheckcpp.cc', code)

    def get_compiler_check_args(self, mode: CompileCheckMode) -> T.List[str]:
        # -fpermissive allows non-conforming code to compile which is necessary
        # for many C++ checks. Particularly, the has_header_symbol check is
        # too strict without this and always fails.
        return super().get_compiler_check_args(mode) + ['-fpermissive']

    def has_header_symbol(self, hname: str, symbol: str, prefix: str,
                          env: 'Environment', *,
                          extra_args: T.Union[None, T.List[str], T.Callable[[CompileCheckMode], T.List[str]]] = None,
                          dependencies: T.Optional[T.List['Dependency']] = None) -> T.Tuple[bool, bool]:
        # Check if it's a C-like symbol
        found, cached = super().has_header_symbol(hname, symbol, prefix, env,
                                                  extra_args=extra_args,
                                                  dependencies=dependencies)
        if found:
            return True, cached
        # Check if it's a class or a template
        if extra_args is None:
            extra_args = []
        t = f'''{prefix}
        #include <{hname}>
        using {symbol};
        int main(void) {{ return 0; }}'''
        return self.compiles(t, env, extra_args=extra_args,
                             dependencies=dependencies)

    def _test_cpp_std_arg(self, cpp_std_value: str) -> bool:
        # Test whether the compiler understands a -std=XY argument
        assert cpp_std_value.startswith('-std=')

        # This test does not use has_multi_arguments() for two reasons:
        # 1. has_multi_arguments() requires an env argument, which the compiler
        #    object does not have at this point.
        # 2. even if it did have an env object, that might contain another more
        #    recent -std= argument, which might lead to a cascaded failure.
        CPP_TEST = 'int i = static_cast<int>(0);'
        with self.compile(CPP_TEST, extra_args=[cpp_std_value], mode=CompileCheckMode.COMPILE) as p:
            if p.returncode == 0:
                mlog.debug(f'Compiler accepts {cpp_std_value}:', 'YES')
                return True
            else:
                mlog.debug(f'Compiler accepts {cpp_std_value}:', 'NO')
                return False

    @functools.lru_cache()
    def _find_best_cpp_std(self, cpp_std: str) -> str:
        # The initial version mapping approach to make falling back
        # from '-std=c++14' to '-std=c++1y' was too brittle. For instance,
        # Apple's Clang uses a different versioning scheme to upstream LLVM,
        # making the whole detection logic awfully brittle. Instead, let's
        # just see if feeding GCC or Clang our '-std=' setting works, and
        # if not, try the fallback argument.
        CPP_FALLBACKS = {
            'c++11': 'c++0x',
            'gnu++11': 'gnu++0x',
            'c++14': 'c++1y',
            'gnu++14': 'gnu++1y',
            'c++17': 'c++1z',
            'gnu++17': 'gnu++1z',
            'c++20': 'c++2a',
            'gnu++20': 'gnu++2a',
            'c++23': 'c++2b',
            'gnu++23': 'gnu++2b',
            'c++26': 'c++2c',
            'gnu++26': 'gnu++2c',
        }

        # Currently, remapping is only supported for Clang, Elbrus and GCC
        assert self.id in frozenset(['clang', 'lcc', 'gcc', 'emscripten', 'armltdclang', 'intel-llvm'])

        if cpp_std not in CPP_FALLBACKS:
            # 'c++03' and 'c++98' don't have fallback types
            return '-std=' + cpp_std

        for i in (cpp_std, CPP_FALLBACKS[cpp_std]):
            cpp_std_value = '-std=' + i
            if self._test_cpp_std_arg(cpp_std_value):
                return cpp_std_value

        raise MesonException(f'C++ Compiler does not support -std={cpp_std}')

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = super().get_options()
        key = self.form_compileropt_key('std')
        opts.update({
            key: options.UserStdOption('C++', _ALL_STDS),
        })
        return opts


class _StdCPPLibMixin(CompilerMixinBase):

    """Detect whether to use libc++ or libstdc++."""

    def language_stdlib_provider(self, env: Environment) -> str:
        # https://stackoverflow.com/a/31658120
        header = 'version' if self.has_header('version', '', env)[0] else 'ciso646'
        is_libcxx = self.has_header_symbol(header, '_LIBCPP_VERSION', '', env)[0]
        lib = 'c++' if is_libcxx else 'stdc++'
        return lib

    @functools.lru_cache(None)
    def language_stdlib_only_link_flags(self, env: Environment) -> T.List[str]:
        """Detect the C++ stdlib and default search dirs

        As an optimization, this method will cache the value, to avoid building the same values over and over

        :param env: An Environment object
        :raises MesonException: If a stdlib cannot be determined
        """

        # We need to apply the search prefix here, as these link arguments may
        # be passed to a different compiler with a different set of default
        # search paths, such as when using Clang for C/C++ and gfortran for
        # fortran.
        search_dirs = [f'-L{d}' for d in self.get_compiler_dirs(env, 'libraries')]

        machine = env.machines[self.for_machine]
        assert machine is not None, 'for mypy'

        lib = self.language_stdlib_provider(env)
        if self.find_library(lib, env, []) is not None:
            return search_dirs + [f'-l{lib}']

        # TODO: maybe a bug exception?
        raise MesonException('Could not detect either libc++ or libstdc++ as your C++ stdlib implementation.')


class ClangCPPCompiler(_StdCPPLibMixin, ClangCompiler, CPPCompiler):

    _CPP23_VERSION = '>=12.0.0'
    _CPP26_VERSION = '>=17.0.0'

    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice, is_cross: bool,
                 info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 defines: T.Optional[T.Dict[str, str]] = None,
                 full_version: T.Optional[str] = None):
        CPPCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                             info, linker=linker, full_version=full_version)
        ClangCompiler.__init__(self, defines)
        default_warn_args = ['-Wall', '-Winvalid-pch']
        self.warn_args = {'0': [],
                          '1': default_warn_args,
                          '2': default_warn_args + ['-Wextra'],
                          '3': default_warn_args + ['-Wextra', '-Wpedantic'],
                          'everything': ['-Weverything']}

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = CPPCompiler.get_options(self)
        self.update_options(
            opts,
            self.create_option(options.UserComboOption,
                               self.form_compileropt_key('eh'),
                               'C++ exception handling type.',
                               ['none', 'default', 'a', 's', 'sc'],
                               'default'),
            self.create_option(options.UserBooleanOption,
                               self.form_compileropt_key('rtti'),
                               'Enable RTTI',
                               True),
            self.create_option(options.UserBooleanOption,
                               self.form_compileropt_key('debugstl'),
                               'STL debug mode',
                               False),
        )
        cppstd_choices = [
            'c++98', 'c++03', 'c++11', 'c++14', 'c++17', 'c++1z', 'c++2a', 'c++20',
        ]
        if version_compare(self.version, self._CPP23_VERSION):
            cppstd_choices.append('c++23')
        if version_compare(self.version, self._CPP26_VERSION):
            cppstd_choices.append('c++26')
        std_opt = opts[self.form_compileropt_key('std')]
        assert isinstance(std_opt, options.UserStdOption), 'for mypy'
        std_opt.set_versions(cppstd_choices, gnu=True)
        if self.info.is_windows() or self.info.is_cygwin():
            self.update_options(
                opts,
                self.create_option(options.UserArrayOption,
                                   self.form_compileropt_key('winlibs'),
                                   'Standard Win libraries to link against',
                                   gnu_winlibs),
            )
        return opts

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        args: T.List[str] = []
        key = self.form_compileropt_key('std')
        std = options.get_value(key)
        if std != 'none':
            args.append(self._find_best_cpp_std(std))

        key = self.form_compileropt_key('eh')
        non_msvc_eh_options(options.get_value(key), args)

        key = self.form_compileropt_key('debugstl')
        if options.get_value(key):
            args.append('-D_GLIBCXX_DEBUG=1')

            # We can't do _LIBCPP_DEBUG because it's unreliable unless libc++ was built with it too:
            # https://discourse.llvm.org/t/building-a-program-with-d-libcpp-debug-1-against-a-libc-that-is-not-itself-built-with-that-define/59176/3
            # Note that unlike _GLIBCXX_DEBUG, _MODE_DEBUG doesn't break ABI. It's just slow.
            if version_compare(self.version, '>=18'):
                args.append('-D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_DEBUG')

        key = self.form_compileropt_key('rtti')
        if not options.get_value(key):
            args.append('-fno-rtti')

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

    def get_assert_args(self, disable: bool, env: 'Environment') -> T.List[str]:
        if disable:
            return ['-DNDEBUG']

        # Don't inject the macro if the compiler already has it pre-defined.
        for macro in ['_GLIBCXX_ASSERTIONS', '_LIBCPP_HARDENING_MODE', '_LIBCPP_ENABLE_ASSERTIONS']:
            if self.defines.get(macro) is not None:
                return []

        if self.language_stdlib_provider(env) == 'stdc++':
            return ['-D_GLIBCXX_ASSERTIONS=1']
        else:
            if version_compare(self.version, '>=18'):
                return ['-D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_FAST']
            elif version_compare(self.version, '>=15'):
                return ['-D_LIBCPP_ENABLE_ASSERTIONS=1']

        return []


class ArmLtdClangCPPCompiler(ClangCPPCompiler):

    id = 'armltdclang'


class AppleClangCPPCompiler(AppleCompilerMixin, ClangCPPCompiler):

    _CPP23_VERSION = '>=13.0.0'
    # TODO: We don't know which XCode version will include LLVM 17 yet, so
    # use something absurd.
    _CPP26_VERSION = '>=99.0.0'


class EmscriptenCPPCompiler(EmscriptenMixin, ClangCPPCompiler):

    id = 'emscripten'

    # Emscripten uses different version numbers than Clang; `emcc -v` will show
    # the Clang version number used as well (but `emcc --version` does not).
    # See https://github.com/pyodide/pyodide/discussions/4762 for more on
    # emcc <--> clang versions. Note, although earlier versions claim to be the
    # Clang versions 12.0.0 and 17.0.0 required for these C++ standards, they
    # only accept the flags in the later versions below.
    _CPP23_VERSION = '>=2.0.10'
    _CPP26_VERSION = '>=3.1.39'

    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice, is_cross: bool,
                 info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 defines: T.Optional[T.Dict[str, str]] = None,
                 full_version: T.Optional[str] = None):
        if not is_cross:
            raise MesonException('Emscripten compiler can only be used for cross compilation.')
        if not version_compare(version, '>=1.39.19'):
            raise MesonException('Meson requires Emscripten >= 1.39.19')
        ClangCPPCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                                  info, linker=linker,
                                  defines=defines, full_version=full_version)

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        args: T.List[str] = []
        key = self.form_compileropt_key('std')
        std = options.get_value(key)
        if std != 'none':
            args.append(self._find_best_cpp_std(std))
        return args


class ArmclangCPPCompiler(ArmclangCompiler, CPPCompiler):
    '''
    Keil armclang
    '''

    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice, is_cross: bool,
                 info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CPPCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                             info, linker=linker, full_version=full_version)
        ArmclangCompiler.__init__(self)
        default_warn_args = ['-Wall', '-Winvalid-pch']
        self.warn_args = {'0': [],
                          '1': default_warn_args,
                          '2': default_warn_args + ['-Wextra'],
                          '3': default_warn_args + ['-Wextra', '-Wpedantic'],
                          'everything': ['-Weverything']}

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = CPPCompiler.get_options(self)
        key = self.form_compileropt_key('std')
        self.update_options(
            opts,
            self.create_option(options.UserComboOption,
                               key.evolve('eh'),
                               'C++ exception handling type.',
                               ['none', 'default', 'a', 's', 'sc'],
                               'default'),
        )
        std_opt = opts[key]
        assert isinstance(std_opt, options.UserStdOption), 'for mypy'
        std_opt.set_versions(['c++98', 'c++03', 'c++11', 'c++14', 'c++17'], gnu=True)
        return opts

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        args: T.List[str] = []
        key = self.form_compileropt_key('std')
        std = options.get_value(key)
        if std != 'none':
            args.append('-std=' + std)

        key = self.form_compileropt_key('eh')
        non_msvc_eh_options(options.get_value(key), args)

        return args

    def get_option_link_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        return []


class GnuCPPCompiler(_StdCPPLibMixin, GnuCompiler, CPPCompiler):
    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice, is_cross: bool,
                 info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 defines: T.Optional[T.Dict[str, str]] = None,
                 full_version: T.Optional[str] = None):
        CPPCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                             info, linker=linker, full_version=full_version)
        GnuCompiler.__init__(self, defines)
        default_warn_args = ['-Wall', '-Winvalid-pch']
        self.warn_args = {'0': [],
                          '1': default_warn_args,
                          '2': default_warn_args + ['-Wextra'],
                          '3': default_warn_args + ['-Wextra', '-Wpedantic'],
                          'everything': (default_warn_args + ['-Wextra', '-Wpedantic'] +
                                         self.supported_warn_args(gnu_common_warning_args) +
                                         self.supported_warn_args(gnu_cpp_warning_args))}

    def get_options(self) -> 'MutableKeyedOptionDictType':
        key = self.form_compileropt_key('std')
        opts = CPPCompiler.get_options(self)
        self.update_options(
            opts,
            self.create_option(options.UserComboOption,
                               self.form_compileropt_key('eh'),
                               'C++ exception handling type.',
                               ['none', 'default', 'a', 's', 'sc'],
                               'default'),
            self.create_option(options.UserBooleanOption,
                               self.form_compileropt_key('rtti'),
                               'Enable RTTI',
                               True),
            self.create_option(options.UserBooleanOption,
                               self.form_compileropt_key('debugstl'),
                               'STL debug mode',
                               False),
        )
        cppstd_choices = [
            'c++98', 'c++03', 'c++11', 'c++14', 'c++17', 'c++1z',
            'c++2a', 'c++20',
        ]
        if version_compare(self.version, '>=11.0.0'):
            cppstd_choices.append('c++23')
        if version_compare(self.version, '>=14.0.0'):
            cppstd_choices.append('c++26')
        std_opt = opts[key]
        assert isinstance(std_opt, options.UserStdOption), 'for mypy'
        std_opt.set_versions(cppstd_choices, gnu=True)
        if self.info.is_windows() or self.info.is_cygwin():
            self.update_options(
                opts,
                self.create_option(options.UserArrayOption,
                                   key.evolve('cpp_winlibs'),
                                   'Standard Win libraries to link against',
                                   gnu_winlibs),
            )
        return opts

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        args: T.List[str] = []
        stdkey = self.form_compileropt_key('std')
        ehkey = self.form_compileropt_key('eh')
        rttikey = self.form_compileropt_key('rtti')
        debugstlkey = self.form_compileropt_key('debugstl')

        std = options.get_value(stdkey)
        if std != 'none':
            args.append(self._find_best_cpp_std(std))

        non_msvc_eh_options(options.get_value(ehkey), args)

        if not options.get_value(rttikey):
            args.append('-fno-rtti')

        if options.get_value(debugstlkey):
            args.append('-D_GLIBCXX_DEBUG=1')
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

    def get_assert_args(self, disable: bool, env: 'Environment') -> T.List[str]:
        if disable:
            return ['-DNDEBUG']

        # Don't inject the macro if the compiler already has it pre-defined.
        for macro in ['_GLIBCXX_ASSERTIONS', '_LIBCPP_HARDENING_MODE', '_LIBCPP_ENABLE_ASSERTIONS']:
            if self.defines.get(macro) is not None:
                return []

        if self.language_stdlib_provider(env) == 'stdc++':
            return ['-D_GLIBCXX_ASSERTIONS=1']
        else:
            if version_compare(self.version, '>=18'):
                return ['-D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_FAST']
            elif version_compare(self.version, '>=15'):
                return ['-D_LIBCPP_ENABLE_ASSERTIONS=1']

        return []

    def get_pch_use_args(self, pch_dir: str, header: str) -> T.List[str]:
        return ['-fpch-preprocess', '-include', os.path.basename(header)]


class PGICPPCompiler(PGICompiler, CPPCompiler):
    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice, is_cross: bool,
                 info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CPPCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                             info, linker=linker, full_version=full_version)
        PGICompiler.__init__(self)


class NvidiaHPC_CPPCompiler(PGICompiler, CPPCompiler):

    id = 'nvidia_hpc'

    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice, is_cross: bool,
                 info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CPPCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                             info, linker=linker, full_version=full_version)
        PGICompiler.__init__(self)

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = CPPCompiler.get_options(self)
        cppstd_choices = [
            'c++98', 'c++03', 'c++11', 'c++14', 'c++17', 'c++20', 'c++23',
            'gnu++98', 'gnu++03', 'gnu++11', 'gnu++14', 'gnu++17', 'gnu++20'
        ]
        std_opt = opts[self.form_compileropt_key('std')]
        assert isinstance(std_opt, options.UserStdOption), 'for mypy'
        std_opt.set_versions(cppstd_choices)
        return opts


class ElbrusCPPCompiler(ElbrusCompiler, CPPCompiler):
    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice, is_cross: bool,
                 info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 defines: T.Optional[T.Dict[str, str]] = None,
                 full_version: T.Optional[str] = None):
        CPPCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                             info, linker=linker, full_version=full_version)
        ElbrusCompiler.__init__(self)

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = CPPCompiler.get_options(self)

        cpp_stds = ['c++98']
        if version_compare(self.version, '>=1.20.00'):
            cpp_stds += ['c++03', 'c++0x', 'c++11']
        if version_compare(self.version, '>=1.21.00') and version_compare(self.version, '<1.22.00'):
            cpp_stds += ['c++14', 'c++1y']
        if version_compare(self.version, '>=1.22.00'):
            cpp_stds += ['c++14']
        if version_compare(self.version, '>=1.23.00'):
            cpp_stds += ['c++1y']
        if version_compare(self.version, '>=1.24.00'):
            cpp_stds += ['c++1z', 'c++17']
        if version_compare(self.version, '>=1.25.00'):
            cpp_stds += ['c++2a']
        if version_compare(self.version, '>=1.26.00'):
            cpp_stds += ['c++20']

        key = self.form_compileropt_key('std')
        self.update_options(
            opts,
            self.create_option(options.UserComboOption,
                               self.form_compileropt_key('eh'),
                               'C++ exception handling type.',
                               ['none', 'default', 'a', 's', 'sc'],
                               'default'),
            self.create_option(options.UserBooleanOption,
                               self.form_compileropt_key('debugstl'),
                               'STL debug mode',
                               False),
        )
        std_opt = opts[key]
        assert isinstance(std_opt, options.UserStdOption), 'for mypy'
        std_opt.set_versions(cpp_stds, gnu=True)
        return opts

    # Elbrus C++ compiler does not have lchmod, but there is only linker warning, not compiler error.
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

    # Elbrus C++ compiler does not support RTTI, so don't check for it.
    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        args: T.List[str] = []
        key = self.form_compileropt_key('std')
        std = options.get_value(key)
        if std != 'none':
            args.append(self._find_best_cpp_std(std))

        key = self.form_compileropt_key('eh')
        non_msvc_eh_options(options.get_value(key), args)

        key = self.form_compileropt_key('debugstl')
        if options.get_value(key):
            args.append('-D_GLIBCXX_DEBUG=1')
        return args


class IntelCPPCompiler(IntelGnuLikeCompiler, CPPCompiler):
    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice, is_cross: bool,
                 info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CPPCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                             info, linker=linker, full_version=full_version)
        IntelGnuLikeCompiler.__init__(self)
        self.lang_header = 'c++-header'
        default_warn_args = ['-Wall', '-w3', '-Wpch-messages']
        self.warn_args = {'0': [],
                          '1': default_warn_args + ['-diag-disable:remark'],
                          '2': default_warn_args + ['-Wextra', '-diag-disable:remark'],
                          '3': default_warn_args + ['-Wextra', '-diag-disable:remark'],
                          'everything': default_warn_args + ['-Wextra']}

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = CPPCompiler.get_options(self)
        # Every Unix compiler under the sun seems to accept -std=c++03,
        # with the exception of ICC. Instead of preventing the user from
        # globally requesting C++03, we transparently remap it to C++98
        c_stds = ['c++98', 'c++03']
        g_stds = ['gnu++98', 'gnu++03']
        if version_compare(self.version, '>=15.0.0'):
            c_stds += ['c++11', 'c++14']
            g_stds += ['gnu++11']
        if version_compare(self.version, '>=16.0.0'):
            c_stds += ['c++17']
        if version_compare(self.version, '>=17.0.0'):
            g_stds += ['gnu++14']
        if version_compare(self.version, '>=19.1.0'):
            c_stds += ['c++2a']
            g_stds += ['gnu++2a']

        key = self.form_compileropt_key('std')
        self.update_options(
            opts,
            self.create_option(options.UserComboOption,
                               self.form_compileropt_key('eh'),
                               'C++ exception handling type.',
                               ['none', 'default', 'a', 's', 'sc'],
                               'default'),
            self.create_option(options.UserBooleanOption,
                               self.form_compileropt_key('rtti'),
                               'Enable RTTI',
                               True),
            self.create_option(options.UserBooleanOption,
                               self.form_compileropt_key('debugstl'),
                               'STL debug mode',
                               False),
        )
        std_opt = opts[key]
        assert isinstance(std_opt, options.UserStdOption), 'for mypy'
        std_opt.set_versions(c_stds + g_stds)
        return opts

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        args: T.List[str] = []
        key = self.form_compileropt_key('std')
        std = options.get_value(key)
        if std != 'none':
            remap_cpp03 = {
                'c++03': 'c++98',
                'gnu++03': 'gnu++98'
            }
            args.append('-std=' + remap_cpp03.get(std, std))
        if options.get_value(key.evolve('eh')) == 'none':
            args.append('-fno-exceptions')
        if not options.get_value(key.evolve('rtti')):
            args.append('-fno-rtti')
        if options.get_value(key.evolve('debugstl')):
            args.append('-D_GLIBCXX_DEBUG=1')
        return args

    def get_option_link_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        return []


class IntelLLVMCPPCompiler(ClangCPPCompiler):

    id = 'intel-llvm'


class VisualStudioLikeCPPCompilerMixin(CompilerMixinBase):

    """Mixin for C++ specific method overrides in MSVC-like compilers."""

    VC_VERSION_MAP = {
        'none': (True, None),
        'vc++11': (True, 11),
        'vc++14': (True, 14),
        'vc++17': (True, 17),
        'vc++20': (True, 20),
        'vc++latest': (True, "latest"),
        'c++11': (False, 11),
        'c++14': (False, 14),
        'c++17': (False, 17),
        'c++20': (False, 20),
        'c++latest': (False, "latest"),
    }

    def get_option_link_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        # need a typeddict for this
        key = self.form_compileropt_key('winlibs')
        return T.cast('T.List[str]', options.get_value(key)[:])

    def _get_options_impl(self, opts: 'MutableKeyedOptionDictType', cpp_stds: T.List[str]) -> 'MutableKeyedOptionDictType':
        key = self.form_compileropt_key('std')
        self.update_options(
            opts,
            self.create_option(options.UserComboOption,
                               self.form_compileropt_key('eh'),
                               'C++ exception handling type.',
                               ['none', 'default', 'a', 's', 'sc'],
                               'default'),
            self.create_option(options.UserBooleanOption,
                               self.form_compileropt_key('rtti'),
                               'Enable RTTI',
                               True),
            self.create_option(options.UserArrayOption,
                               self.form_compileropt_key('winlibs'),
                               'Windows libs to link against.',
                               msvc_winlibs),
        )
        std_opt = opts[key]
        assert isinstance(std_opt, options.UserStdOption), 'for mypy'
        std_opt.set_versions(cpp_stds)
        return opts

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        args: T.List[str] = []
        key = self.form_compileropt_key('std')

        eh = options.get_value(self.form_compileropt_key('eh'))
        if eh == 'default':
            args.append('/EHsc')
        elif eh == 'none':
            args.append('/EHs-c-')
        else:
            args.append('/EH' + eh)

        if not options.get_value(self.form_compileropt_key('rtti')):
            args.append('/GR-')

        permissive, ver = self.VC_VERSION_MAP[options.get_value(key)]

        if ver is not None:
            args.append(f'/std:c++{ver}')

        if not permissive:
            args.append('/permissive-')

        return args

    def get_compiler_check_args(self, mode: CompileCheckMode) -> T.List[str]:
        # XXX: this is a hack because so much GnuLike stuff is in the base CPPCompiler class.
        return Compiler.get_compiler_check_args(self, mode)


class CPP11AsCPP14Mixin(CompilerMixinBase):

    """Mixin class for VisualStudio and ClangCl to replace C++11 std with C++14.

    This is a limitation of Clang and MSVC that ICL doesn't share.
    """

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        # Note: there is no explicit flag for supporting C++11; we attempt to do the best we can
        # which means setting the C++ standard version to C++14, in compilers that support it
        # (i.e., after VS2015U3)
        # if one is using anything before that point, one cannot set the standard.
        key = self.form_compileropt_key('std')
        if options.get_value(key) in {'vc++11', 'c++11'}:
            mlog.warning(self.id, 'does not support C++11;',
                         'attempting best effort; setting the standard to C++14',
                         once=True, fatal=False)
            # Don't mutate anything we're going to change, we need to use
            # deepcopy since we're messing with members, and we can't simply
            # copy the members because the option proxy doesn't support it.
            options = copy.deepcopy(options)
            if options.get_value(key) == 'vc++11':
                options.set_value(key, 'vc++14')
            else:
                options.set_value(key,  'c++14')
        return super().get_option_compile_args(options)


class VisualStudioCPPCompiler(CPP11AsCPP14Mixin, VisualStudioLikeCPPCompilerMixin, MSVCCompiler, CPPCompiler):

    id = 'msvc'

    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice,
                 is_cross: bool, info: 'MachineInfo', target: str,
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CPPCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                             info, linker=linker, full_version=full_version)
        MSVCCompiler.__init__(self, target)

        # By default, MSVC has a broken __cplusplus define that pretends to be c++98:
        # https://docs.microsoft.com/en-us/cpp/build/reference/zc-cplusplus?view=msvc-160
        # Pass the flag to enable a truthful define, if possible.
        if version_compare(self.version, '>= 19.14.26428'):
            self.always_args = self.always_args + ['/Zc:__cplusplus']

    def get_options(self) -> 'MutableKeyedOptionDictType':
        cpp_stds = ['none', 'c++11', 'vc++11']
        # Visual Studio 2015 and later
        if version_compare(self.version, '>=19'):
            cpp_stds.extend(['c++14', 'c++latest', 'vc++latest'])
        # Visual Studio 2017 and later
        if version_compare(self.version, '>=19.11'):
            cpp_stds.extend(['vc++14', 'c++17', 'vc++17'])
        if version_compare(self.version, '>=19.29'):
            cpp_stds.extend(['c++20', 'vc++20'])
        return self._get_options_impl(super().get_options(), cpp_stds)

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        key = self.form_compileropt_key('std')
        if options.get_value(key) != 'none' and version_compare(self.version, '<19.00.24210'):
            mlog.warning('This version of MSVC does not support cpp_std arguments', fatal=False)
            options = copy.copy(options)
            options.set_value(key, 'none')

        args = super().get_option_compile_args(options)

        if version_compare(self.version, '<19.11'):
            try:
                i = args.index('/permissive-')
            except ValueError:
                return args
            del args[i]
        return args

class ClangClCPPCompiler(CPP11AsCPP14Mixin, VisualStudioLikeCPPCompilerMixin, ClangClCompiler, CPPCompiler):

    id = 'clang-cl'

    def __init__(self, exelist: T.List[str], version: str, for_machine: MachineChoice,
                 is_cross: bool, info: 'MachineInfo', target: str,
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CPPCompiler.__init__(self, [], exelist, version, for_machine, is_cross,
                             info, linker=linker, full_version=full_version)
        ClangClCompiler.__init__(self, target)

    def get_options(self) -> 'MutableKeyedOptionDictType':
        cpp_stds = ['none', 'c++11', 'vc++11', 'c++14', 'vc++14', 'c++17', 'vc++17', 'c++20', 'vc++20', 'c++latest']
        return self._get_options_impl(super().get_options(), cpp_stds)


class IntelClCPPCompiler(VisualStudioLikeCPPCompilerMixin, IntelVisualStudioLikeCompiler, CPPCompiler):

    def __init__(self, exelist: T.List[str], version: str, for_machine: MachineChoice,
                 is_cross: bool, info: 'MachineInfo', target: str,
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CPPCompiler.__init__(self, [], exelist, version, for_machine, is_cross,
                             info, linker=linker, full_version=full_version)
        IntelVisualStudioLikeCompiler.__init__(self, target)

    def get_options(self) -> 'MutableKeyedOptionDictType':
        # This has only been tested with version 19.0,
        cpp_stds = ['none', 'c++11', 'vc++11', 'c++14', 'vc++14', 'c++17', 'vc++17', 'c++latest']
        return self._get_options_impl(super().get_options(), cpp_stds)

    def get_compiler_check_args(self, mode: CompileCheckMode) -> T.List[str]:
        # XXX: this is a hack because so much GnuLike stuff is in the base CPPCompiler class.
        return IntelVisualStudioLikeCompiler.get_compiler_check_args(self, mode)


class IntelLLVMClCPPCompiler(IntelClCPPCompiler):

    id = 'intel-llvm-cl'


class ArmCPPCompiler(ArmCompiler, CPPCompiler):
    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice, is_cross: bool,
                 info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CPPCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                             info, linker=linker, full_version=full_version)
        ArmCompiler.__init__(self)

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = CPPCompiler.get_options(self)
        std_opt = self.form_compileropt_key('std')
        assert isinstance(std_opt, options.UserStdOption), 'for mypy'
        std_opt.set_versions(['c++03', 'c++11'])
        return opts

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        args: T.List[str] = []
        key = self.form_compileropt_key('std')
        std = options.get_value(key)
        if std == 'c++11':
            args.append('--cpp11')
        elif std == 'c++03':
            args.append('--cpp')
        return args

    def get_option_link_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        return []

    def get_compiler_check_args(self, mode: CompileCheckMode) -> T.List[str]:
        return []


class CcrxCPPCompiler(CcrxCompiler, CPPCompiler):
    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice, is_cross: bool,
                 info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CPPCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                             info, linker=linker, full_version=full_version)
        CcrxCompiler.__init__(self)

    # Override CCompiler.get_always_args
    def get_always_args(self) -> T.List[str]:
        return ['-nologo', '-lang=cpp']

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        return []

    def get_compile_only_args(self) -> T.List[str]:
        return []

    def get_output_args(self, outputname: str) -> T.List[str]:
        return [f'-output=obj={outputname}']

    def get_option_link_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        return []

    def get_compiler_check_args(self, mode: CompileCheckMode) -> T.List[str]:
        return []

class TICPPCompiler(TICompiler, CPPCompiler):
    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice, is_cross: bool,
                 info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CPPCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                             info, linker=linker, full_version=full_version)
        TICompiler.__init__(self)

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = CPPCompiler.get_options(self)
        key = self.form_compileropt_key('std')
        std_opt = opts[key]
        assert isinstance(std_opt, options.UserStdOption), 'for mypy'
        std_opt.set_versions(['c++03'])
        return opts

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        args: T.List[str] = []
        key = self.form_compileropt_key('std')
        std = options.get_value(key)
        if std != 'none':
            args.append('--' + std)
        return args

    def get_always_args(self) -> T.List[str]:
        return []

    def get_option_link_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        return []

class C2000CPPCompiler(TICPPCompiler):
    # Required for backwards compat with projects created before ti-cgt support existed
    id = 'c2000'

class C6000CPPCompiler(TICPPCompiler):
    id = 'c6000'

class MetrowerksCPPCompilerARM(MetrowerksCompiler, CPPCompiler):
    id = 'mwccarm'

    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice,
                 is_cross: bool, info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CPPCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                             info, linker=linker, full_version=full_version)
        MetrowerksCompiler.__init__(self)

    def get_instruction_set_args(self, instruction_set: str) -> T.Optional[T.List[str]]:
        return mwccarm_instruction_set_args.get(instruction_set, None)

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = CPPCompiler.get_options(self)
        key = self.form_compileropt_key('std')
        opts[key].choices = ['none']
        return opts

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        args: T.List[str] = []
        key = self.form_compileropt_key('std')
        std = options.get_value(key)
        if std != 'none':
            args.append('-lang')
            args.append(std)
        return args

class MetrowerksCPPCompilerEmbeddedPowerPC(MetrowerksCompiler, CPPCompiler):
    id = 'mwcceppc'

    def __init__(self, ccache: T.List[str], exelist: T.List[str], version: str, for_machine: MachineChoice,
                 is_cross: bool, info: 'MachineInfo',
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None):
        CPPCompiler.__init__(self, ccache, exelist, version, for_machine, is_cross,
                             info, linker=linker, full_version=full_version)
        MetrowerksCompiler.__init__(self)

    def get_instruction_set_args(self, instruction_set: str) -> T.Optional[T.List[str]]:
        return mwcceppc_instruction_set_args.get(instruction_set, None)

    def get_options(self) -> 'MutableKeyedOptionDictType':
        opts = CPPCompiler.get_options(self)
        key = self.form_compileropt_key('std')
        opts[key].choices = ['none']
        return opts

    def get_option_compile_args(self, options: 'KeyedOptionDictType') -> T.List[str]:
        args: T.List[str] = []
        key = self.form_compileropt_key('std')
        std = options.get_value(key)
        if std != 'none':
            args.append('-lang ' + std)
        return args
