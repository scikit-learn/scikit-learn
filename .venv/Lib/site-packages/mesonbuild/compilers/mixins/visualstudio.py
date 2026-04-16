# SPDX-License-Identifier: Apache-2.0
# Copyright 2019 The meson development team
# Copyright © 2023 Intel Corporation

from __future__ import annotations

"""Abstractions to simplify compilers that implement an MSVC compatible
interface.
"""

import abc
import os
import typing as T

from ... import arglist
from ... import mesonlib
from mesonbuild.compilers.compilers import CompileCheckMode
from ...options import OptionKey
from mesonbuild.linkers.linkers import ClangClDynamicLinker, MSVCDynamicLinker

if T.TYPE_CHECKING:
    from ...build import BuildTarget
    from ...environment import Environment
    from .clike import CLikeCompiler as Compiler
else:
    # This is a bit clever, for mypy we pretend that these mixins descend from
    # Compiler, so we get all of the methods and attributes defined for us, but
    # for runtime we make them descend from object (which all classes normally
    # do). This gives up DRYer type checking, with no runtime impact
    Compiler = object

vs32_instruction_set_args: T.Dict[str, T.Optional[T.List[str]]] = {
    'mmx': ['/arch:SSE'], # There does not seem to be a flag just for MMX
    'sse': ['/arch:SSE'],
    'sse2': ['/arch:SSE2'],
    'sse3': ['/arch:AVX'], # VS leaped from SSE2 directly to AVX.
    'sse41': ['/arch:AVX'],
    'sse42': ['/arch:AVX'],
    'avx': ['/arch:AVX'],
    'avx2': ['/arch:AVX2'],
    'neon': None,
}

# The 64 bit compiler defaults to /arch:avx.
vs64_instruction_set_args: T.Dict[str, T.Optional[T.List[str]]] = {
    'mmx': ['/arch:AVX'],
    'sse': ['/arch:AVX'],
    'sse2': ['/arch:AVX'],
    'sse3': ['/arch:AVX'],
    'ssse3': ['/arch:AVX'],
    'sse41': ['/arch:AVX'],
    'sse42': ['/arch:AVX'],
    'avx': ['/arch:AVX'],
    'avx2': ['/arch:AVX2'],
    'neon': None,
}

msvc_optimization_args: T.Dict[str, T.List[str]] = {
    'plain': [],
    '0': ['/Od'],
    'g': [], # No specific flag to optimize debugging, /Zi or /ZI will create debug information
    '1': ['/O1'],
    '2': ['/O2'],
    '3': ['/O2', '/Gw'],
    's': ['/O1', '/Gw'],
}


class VisualStudioLikeCompiler(Compiler, metaclass=abc.ABCMeta):

    """A common interface for all compilers implementing an MSVC-style
    interface.

    A number of compilers attempt to mimic MSVC, with varying levels of
    success, such as Clang-CL and ICL (the Intel C/C++ Compiler for Windows).
    This class implements as much common logic as possible.
    """

    std_warn_args = ['/W3']
    std_opt_args = ['/O2']
    ignore_libs = arglist.UNIXY_COMPILER_INTERNAL_LIBS + ['execinfo']
    internal_libs: T.List[str] = []

    crt_args: T.Dict[str, T.List[str]] = {
        'none': [],
        'md': ['/MD'],
        'mdd': ['/MDd'],
        'mt': ['/MT'],
        'mtd': ['/MTd'],
    }

    # /showIncludes is needed for build dependency tracking in Ninja
    # See: https://ninja-build.org/manual.html#_deps
    # Assume UTF-8 sources by default, but self.unix_args_to_native() removes it
    # if `/source-charset` is set too.
    # It is also dropped if Visual Studio 2013 or earlier is used, since it would
    # not be supported in that case.
    always_args = ['/nologo', '/showIncludes', '/utf-8']
    warn_args: T.Dict[str, T.List[str]] = {
        '0': [],
        '1': ['/W2'],
        '2': ['/W3'],
        '3': ['/W4'],
        'everything': ['/Wall'],
    }

    USED_FOR_SEPARATE_LINKING_STEP = False

    def __init__(self, target: str):
        self.base_options = {OptionKey(o) for o in ['b_pch', 'b_ndebug', 'b_vscrt']} # FIXME add lto, pgo and the like
        self.target = target
        self.is_64 = ('x64' in target) or ('x86_64' in target)
        # do some canonicalization of target machine
        if 'x86_64' in target:
            self.machine = 'x64'
        elif '86' in target:
            self.machine = 'x86'
        elif 'aarch64' in target:
            self.machine = 'arm64'
        elif 'arm' in target:
            self.machine = 'arm'
        else:
            self.machine = target
        if mesonlib.version_compare(self.version, '>=19.28.29910'): # VS 16.9.0 includes cl 19.28.29910
            self.base_options.add(OptionKey('b_sanitize'))
        assert self.linker is not None
        self.linker.machine = self.machine

    # Override CCompiler.get_always_args
    def get_always_args(self) -> T.List[str]:
        # TODO: use ImmutableListProtocol[str] here instead
        return self.always_args.copy()

    def get_pch_suffix(self) -> str:
        return 'pch'

    def get_pch_name(self, name: str) -> str:
        chopped = os.path.basename(name).split('.')[:-1]
        chopped.append(self.get_pch_suffix())
        pchname = '.'.join(chopped)
        return pchname

    def get_pch_base_name(self, header: str) -> str:
        # This needs to be implemented by inheriting classes
        raise NotImplementedError

    def get_pch_use_args(self, pch_dir: str, header: str) -> T.List[str]:
        base = self.get_pch_base_name(header)
        pchname = self.get_pch_name(header)
        return ['/FI' + base, '/Yu' + base, '/Fp' + os.path.join(pch_dir, pchname)]

    def get_preprocess_only_args(self) -> T.List[str]:
        return ['/EP']

    def get_preprocess_to_file_args(self) -> T.List[str]:
        return ['/EP', '/P']

    def get_compile_only_args(self) -> T.List[str]:
        return ['/c']

    def get_no_optimization_args(self) -> T.List[str]:
        return ['/Od', '/Oi-']

    def sanitizer_compile_args(self, target: T.Optional[BuildTarget], value: T.List[str]) -> T.List[str]:
        if not value:
            return value
        return [f'/fsanitize={",".join(value)}']

    def get_output_args(self, outputname: str) -> T.List[str]:
        if self.mode == 'PREPROCESSOR':
            return ['/Fi' + outputname]
        if outputname.endswith('.exe'):
            return ['/Fe' + outputname]
        return ['/Fo' + outputname]

    def get_debug_args(self, is_debug: bool) -> T.List[str]:
        if is_debug:
            return ['/Z7']
        else:
            return []

    def get_optimization_args(self, optimization_level: str) -> T.List[str]:
        args = msvc_optimization_args[optimization_level]
        if mesonlib.version_compare(self.version, '<18.0'):
            args = [arg for arg in args if arg != '/Gw']
        return args

    def linker_to_compiler_args(self, args: T.List[str]) -> T.List[str]:
        return ['/link'] + [arg for arg in args if arg != '/link']

    def get_pic_args(self) -> T.List[str]:
        return [] # PIC is handled by the loader on Windows

    def gen_pch_args(self, header: str, source: str, pchname: str) -> T.Tuple[str, T.List[str]]:
        objname = os.path.splitext(source)[0] + '.obj'
        return objname, ['/Yc' + header, '/Fp' + pchname, '/Fo' + objname]

    def openmp_flags(self) -> T.List[str]:
        return ['/openmp']

    def openmp_link_flags(self) -> T.List[str]:
        return []

    # FIXME, no idea what these should be.
    def thread_flags(self) -> T.List[str]:
        return []

    @classmethod
    def include_arg_to_native(cls, opt: str, path: str) -> str:
        # msvc does not have a concept of system header dirs.
        return f'/I{path}'

    @classmethod
    def unix_args_to_native(cls, args: T.List[str]) -> T.List[str]:
        result: T.List[str] = []
        prev = None
        for i in args:
            if prev:
                i = cls.include_arg_to_native(prev, i)
                prev = None
            # -mms-bitfields is specific to MinGW-GCC
            # -pthread is only valid for GCC
            elif i in {'-mms-bitfields', '-pthread'}:
                continue
            elif i.startswith('-LIBPATH:'):
                i = '/LIBPATH:' + i[9:]
            elif i.startswith('-L'):
                i = '/LIBPATH:' + i[2:]
            # Translate GNU-style -lfoo library name to the import library
            elif i.startswith('-l'):
                name = i[2:]
                if name in cls.ignore_libs:
                    # With MSVC, these are provided by the C runtime which is
                    # linked in by default
                    continue
                else:
                    i = name + '.lib'
            elif i.startswith(('-iquote=', '-isystem=', '-idirafter=')):
                opt, i = i.split('=',  1)
                i = cls.include_arg_to_native(opt, i)
            elif i in {'-iquote', '-isystem', '-idirafter'}:
                prev = i
                continue
            elif i.startswith('-iquote'):
                i = cls.include_arg_to_native('-iquote', i[7:])
            elif i.startswith('-isystem'):
                i = cls.include_arg_to_native('-isystem', i[8:])
            elif i.startswith('-idirafter'):
                i = cls.include_arg_to_native('-idirafter', i[10:])
            # cl.exe does not allow specifying both, so remove /utf-8 that we
            # added automatically in the case the user overrides it manually.
            elif (i.startswith('/source-charset:')
                    or i.startswith('/execution-charset:')
                    or i == '/validate-charset-'):
                try:
                    result.remove('/utf-8')
                except ValueError:
                    pass
            result.append(i)
        return result

    @classmethod
    def native_args_to_unix(cls, args: T.List[str]) -> T.List[str]:
        result: T.List[str] = []
        for arg in args:
            if arg.startswith(('/LIBPATH:', '-LIBPATH:')):
                result.append('-L' + arg[9:])
            elif arg.endswith(('.a', '.lib')) and not mesonlib.path_has_root(arg):
                result.append('-l' + arg)
            else:
                result.append(arg)
        return result

    def get_werror_args(self) -> T.List[str]:
        return ['/WX']

    def get_include_args(self, path: str, is_system: bool) -> T.List[str]:
        if path == '':
            path = '.'
        if is_system:
            # fixed up by unix_args_to_native() for Microsoft cl.exe
            return ['-isystem', path]
        return ['-I' + path]

    def compute_parameters_with_absolute_paths(self, parameter_list: T.List[str], build_dir: str) -> T.List[str]:
        for idx, i in enumerate(parameter_list):
            if i[:2] == '-I' or i[:2] == '/I':
                parameter_list[idx] = i[:2] + os.path.normpath(os.path.join(build_dir, i[2:]))
            elif i[:9] == '/LIBPATH:':
                parameter_list[idx] = i[:9] + os.path.normpath(os.path.join(build_dir, i[9:]))

        return parameter_list

    # Visual Studio is special. It ignores some arguments it does not
    # understand and you can't tell it to error out on those.
    # http://stackoverflow.com/questions/15259720/how-can-i-make-the-microsoft-c-compiler-treat-unknown-flags-as-errors-rather-t
    def has_arguments(self, args: T.List[str], code: str, mode: CompileCheckMode) -> T.Tuple[bool, bool]:
        warning_text = '4044' if mode == CompileCheckMode.LINK else '9002'
        with self._build_wrapper(code, extra_args=args, mode=mode) as p:
            if p.returncode != 0:
                return False, p.cached
            return not (warning_text in p.stderr or warning_text in p.stdout), p.cached

    def get_compile_debugfile_args(self, rel_obj: str, pch: bool = False) -> T.List[str]:
        pdbarr = rel_obj.split('.')[:-1]
        pdbarr += ['pdb']
        args = ['/Fd' + '.'.join(pdbarr)]
        return args

    def get_instruction_set_args(self, instruction_set: str) -> T.Optional[T.List[str]]:
        if self.is_64:
            return vs64_instruction_set_args.get(instruction_set, None)
        return vs32_instruction_set_args.get(instruction_set, None)

    def get_default_include_dirs(self) -> T.List[str]:
        if 'INCLUDE' not in os.environ:
            return []
        return os.environ['INCLUDE'].split(os.pathsep)

    def get_crt_compile_args(self, crt_val: str, env: Environment) -> T.List[str]:
        crt_val = self.get_crt_val(crt_val, env)
        return self.crt_args[crt_val]

    def has_func_attribute(self, name: str) -> T.Tuple[bool, bool]:
        # MSVC doesn't have __attribute__ like Clang and GCC do, so just return
        # false without compiling anything
        return name in {'dllimport', 'dllexport'}, False

    @staticmethod
    def get_argument_syntax() -> str:
        return 'msvc'

    def symbols_have_underscore_prefix(self) -> bool:
        '''
        Check if the compiler prefixes an underscore to global C symbols.

        This overrides the Clike method, as for MSVC checking the
        underscore prefix based on the compiler define never works,
        so do not even try.
        '''
        # Try to consult a hardcoded list of cases we know
        # absolutely have an underscore prefix
        result = self._symbols_have_underscore_prefix_list()
        if result is not None:
            return result

        # As a last resort, try search in a compiled binary
        return self._symbols_have_underscore_prefix_searchbin()

    def get_pie_args(self) -> T.List[str]:
        return []

class MSVCCompiler(VisualStudioLikeCompiler):

    """Specific to the Microsoft Compilers."""

    id = 'msvc'

    def __init__(self, target: str):
        super().__init__(target)

        self.base_options.update({OptionKey('b_lto'), OptionKey('b_lto_mode'), OptionKey('b_pgo')})

        # Visual Studio 2013 and earlier don't support the /utf-8 argument.
        # We want to remove it. We also want to make an explicit copy so we
        # don't mutate class constant state
        if mesonlib.version_compare(self.version, '<19.00') and '/utf-8' in self.always_args:
            self.always_args = [r for r in self.always_args if r != '/utf-8']

    def get_compile_debugfile_args(self, rel_obj: str, pch: bool = False) -> T.List[str]:
        args = super().get_compile_debugfile_args(rel_obj, pch)
        # When generating a PDB file with PCH, all compile commands write
        # to the same PDB file. Hence, we need to serialize the PDB
        # writes using /FS since we do parallel builds. This slows down the
        # build obviously, which is why we only do this when PCH is on.
        # This was added in Visual Studio 2013 (MSVC 18.0). Before that it was
        # always on: https://msdn.microsoft.com/en-us/library/dn502518.aspx
        if pch and mesonlib.version_compare(self.version, '>=18.0'):
            args = ['/FS'] + args
        return args

    # Override CCompiler.get_always_args
    # We want to drop '/utf-8' for Visual Studio 2013 and earlier
    def get_always_args(self) -> T.List[str]:
        return self.always_args

    def get_instruction_set_args(self, instruction_set: str) -> T.Optional[T.List[str]]:
        if self.version.split('.')[0] == '16' and instruction_set == 'avx':
            # VS documentation says that this exists and should work, but
            # it does not. The headers do not contain AVX intrinsics
            # and they cannot be called.
            return None
        return super().get_instruction_set_args(instruction_set)

    def get_pch_base_name(self, header: str) -> str:
        return os.path.basename(header)

    # MSVC requires linking to the generated object file when linking a build target
    # that uses a precompiled header
    def should_link_pch_object(self) -> bool:
        return True

    def get_lto_compile_args(self, *, target: T.Optional[BuildTarget] = None, threads: int = 0,
                             mode: str = 'default') -> T.List[str]:
        args: T.List[str] = ['/GL']
        if mode == 'thin':
            args.append('/Gy')
        return args

    def get_lto_link_args(self, *, target: T.Optional[BuildTarget] = None, threads: int = 0,
                          mode: str = 'default', thinlto_cache_dir: T.Optional[str] = None) -> T.List[str]:
        args: T.List[str] = []
        # LTO data generated by MSVC is only usable by link
        if not isinstance(self.linker, MSVCDynamicLinker):
            raise mesonlib.MesonException(f"MSVC's LTCG only works with link, not {self.linker.id}")
        if mode == 'default':
            args.append('/LTCG')
        elif mode == 'thin':
            args.append('/LTCG:INCREMENTAL')
        return args

    def get_profile_generate_args(self) -> T.List[str]:
        if not isinstance(self.linker, MSVCDynamicLinker):
            raise mesonlib.MesonException(f"MSVC's PGO only works with link, not {self.linker.id}")
        return self.linker_to_compiler_args(['/GENPROFILE'])

    def get_profile_use_args(self) -> T.List[str]:
        if not isinstance(self.linker, MSVCDynamicLinker):
            raise mesonlib.MesonException(f"MSVC's PGO only works with link, not {self.linker.id}")
        return self.linker_to_compiler_args(['/USEPROFILE'])

class ClangClCompiler(VisualStudioLikeCompiler):

    """Specific to Clang-CL."""

    id = 'clang-cl'

    @classmethod
    def include_arg_to_native(cls, opt: str, path: str) -> str:
        # clang-cl does not seem to like a syntax like -iquote=...
        # but unix_args_to_native() canonicalizes opt to not have
        # a trailing equals sign
        return f'/clang:{opt}{path}'

    def __init__(self, target: str):
        super().__init__(target)

        self.base_options.update(
            {OptionKey('b_lto_threads'), OptionKey('b_lto'), OptionKey('b_lto_mode'), OptionKey('b_thinlto_cache'),
             OptionKey('b_thinlto_cache_dir')})

        # Assembly
        self.can_compile_suffixes.add('s')
        self.can_compile_suffixes.add('sx')

    def sanitizer_compile_args(self, target: T.Optional[BuildTarget], value: T.List[str]) -> T.List[str]:
        if not value:
            return value
        args = ['/clang:-fsanitize=' + ','.join(value)]
        if 'address' in value:
            args.append('/clang:-fno-omit-frame-pointer')
        return args

    def has_arguments(self, args: T.List[str], code: str, mode: CompileCheckMode) -> T.Tuple[bool, bool]:
        if mode != CompileCheckMode.LINK:
            args = args + [
                '-Werror=unknown-argument',
                '-Werror=unknown-warning-option',
                '-Werror=unused-command-line-argument',
            ]
        return super().has_arguments(args, code, mode)

    def get_pch_base_name(self, header: str) -> str:
        return header

    @classmethod
    def use_linker_args(cls, linker: str, version: str) -> T.List[str]:
        # Clang additionally can use a linker specified as a path, unlike MSVC.
        if linker == 'lld-link':
            return ['-fuse-ld=lld-link']
        return super().use_linker_args(linker, version)

    def linker_to_compiler_args(self, args: T.List[str]) -> T.List[str]:
        # clang-cl forwards arguments span-wise with the /LINK flag
        # therefore -Wl will be received by lld-link or LINK and rejected
        return super().use_linker_args(self.linker.id, '') + super().linker_to_compiler_args([flag[4:] if flag.startswith('-Wl,') else flag for flag in args])

    def openmp_link_flags(self) -> T.List[str]:
        # see https://github.com/mesonbuild/meson/issues/5298
        libs = self.find_library('libomp', [])
        if libs is None:
            raise mesonlib.MesonBugException('Could not find libomp')
        return super().openmp_link_flags() + libs

    def get_lto_compile_args(self, *, target: T.Optional[BuildTarget] = None, threads: int = 0,
                             mode: str = 'default') -> T.List[str]:
        args: T.List[str] = []
        if mode == 'thin':
            # LTO data generated by clang-cl is only usable by lld-link
            if not isinstance(self.linker, ClangClDynamicLinker):
                raise mesonlib.MesonException(f"LLVM's ThinLTO only works with lld-link, not {self.linker.id}")
            args.append(f'-flto={mode}')
        else:
            assert mode == 'default', 'someone forgot to wire something up'
            args.extend(super().get_lto_compile_args(target=target, threads=threads))
        return args

    def get_lto_link_args(self, *, target: T.Optional[BuildTarget] = None, threads: int = 0,
                          mode: str = 'default', thinlto_cache_dir: T.Optional[str] = None) -> T.List[str]:
        args = []
        if mode == 'thin' and thinlto_cache_dir is not None:
            args.extend(self.linker.get_thinlto_cache_args(thinlto_cache_dir))
        # lld-link /threads:N has the same behaviour as -flto-jobs=N in lld
        if threads > 0:
            # clang-cl was released after clang already had LTO support, so it
            # is safe to assume that all versions of clang-cl support LTO
            args.append(f'/threads:{threads}')
        return args
