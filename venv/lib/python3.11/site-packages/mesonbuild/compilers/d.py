# SPDX-License-Identifier: Apache-2.0
# Copyright 2012-2022 The Meson development team

from __future__ import annotations

import os.path
import re
import subprocess
import typing as T

from .. import mesonlib
from ..arglist import CompilerArgs
from ..linkers import RSPFileSyntax
from ..mesonlib import (
    EnvironmentException, version_compare, is_windows
)
from ..options import OptionKey

from .compilers import (
    clike_debug_args,
    Compiler,
    CompileCheckMode,
)
from .mixins.gnu import GnuCompiler
from .mixins.gnu import gnu_common_warning_args

if T.TYPE_CHECKING:
    from . import compilers
    from ..build import DFeatures
    from ..dependencies import Dependency
    from ..envconfig import MachineInfo
    from ..environment import Environment
    from ..linkers.linkers import DynamicLinker
    from ..mesonlib import MachineChoice

    CompilerMixinBase = Compiler
else:
    CompilerMixinBase = object

d_feature_args: T.Dict[str, T.Dict[str, str]] = {
    'gcc':  {
        'unittest': '-funittest',
        'debug': '-fdebug',
        'version': '-fversion',
        'import_dir': '-J'
    },
    'llvm': {
        'unittest': '-unittest',
        'debug': '-d-debug',
        'version': '-d-version',
        'import_dir': '-J'
    },
    'dmd':  {
        'unittest': '-unittest',
        'debug': '-debug',
        'version': '-version',
        'import_dir': '-J'
    }
}

ldc_optimization_args: T.Dict[str, T.List[str]] = {
    'plain': [],
    '0': [],
    'g': [],
    '1': ['-O1'],
    '2': ['-O2', '-enable-inlining', '-Hkeep-all-bodies'],
    '3': ['-O3', '-enable-inlining', '-Hkeep-all-bodies'],
    's': ['-Oz'],
}

dmd_optimization_args: T.Dict[str, T.List[str]] = {
    'plain': [],
    '0': [],
    'g': [],
    '1': ['-O'],
    '2': ['-O', '-inline'],
    '3': ['-O', '-inline'],
    's': ['-O'],
}

gdc_optimization_args: T.Dict[str, T.List[str]] = {
    'plain': [],
    '0': ['-O0'],
    'g': ['-Og'],
    '1': ['-O1'],
    '2': ['-O2', '-finline-functions'],
    '3': ['-O3', '-finline-functions'],
    's': ['-Os'],
}


class DmdLikeCompilerMixin(CompilerMixinBase):

    """Mixin class for DMD and LDC.

    LDC has a number of DMD like arguments, and this class allows for code
    sharing between them as makes sense.
    """

    def __init__(self, dmd_frontend_version: T.Optional[str]):
        if dmd_frontend_version is None:
            self._dmd_has_depfile = False
        else:
            # -makedeps switch introduced in 2.095 frontend
            self._dmd_has_depfile = version_compare(dmd_frontend_version, ">=2.095.0")

    if T.TYPE_CHECKING:
        mscrt_args: T.Dict[str, T.List[str]] = {}

        def _get_target_arch_args(self) -> T.List[str]: ...

    LINKER_PREFIX = '-L='

    def get_output_args(self, outputname: str) -> T.List[str]:
        return ['-of=' + outputname]

    def get_linker_output_args(self, outputname: str) -> T.List[str]:
        return ['-of=' + outputname]

    def get_include_args(self, path: str, is_system: bool) -> T.List[str]:
        if path == "":
            path = "."
        return ['-I=' + path]

    def compute_parameters_with_absolute_paths(self, parameter_list: T.List[str],
                                               build_dir: str) -> T.List[str]:
        for idx, i in enumerate(parameter_list):
            if i[:3] == '-I=':
                parameter_list[idx] = i[:3] + os.path.normpath(os.path.join(build_dir, i[3:]))
            if i[:4] == '-L-L':
                parameter_list[idx] = i[:4] + os.path.normpath(os.path.join(build_dir, i[4:]))
            if i[:5] == '-L=-L':
                parameter_list[idx] = i[:5] + os.path.normpath(os.path.join(build_dir, i[5:]))
            if i[:6] == '-Wl,-L':
                parameter_list[idx] = i[:6] + os.path.normpath(os.path.join(build_dir, i[6:]))

        return parameter_list

    def get_warn_args(self, level: str) -> T.List[str]:
        return ['-wi']

    def get_werror_args(self) -> T.List[str]:
        return ['-w']

    def get_coverage_args(self) -> T.List[str]:
        return ['-cov']

    def get_coverage_link_args(self) -> T.List[str]:
        return []

    def get_preprocess_only_args(self) -> T.List[str]:
        return ['-E']

    def get_compile_only_args(self) -> T.List[str]:
        return ['-c']

    def get_depfile_suffix(self) -> str:
        return 'deps'

    def get_dependency_gen_args(self, outtarget: str, outfile: str) -> T.List[str]:
        if self._dmd_has_depfile:
            return [f'-makedeps={outfile}']
        return []

    def get_pic_args(self) -> T.List[str]:
        if self.info.is_windows():
            return []
        return ['-fPIC']

    def get_optimization_link_args(self, optimization_level: str) -> T.List[str]:
        if optimization_level != 'plain':
            return self._get_target_arch_args()
        return []

    def gen_import_library_args(self, implibname: str) -> T.List[str]:
        return self.linker.import_library_args(implibname)

    def build_rpath_args(self, env: 'Environment', build_dir: str, from_dir: str,
                         rpath_paths: T.Tuple[str, ...], build_rpath: str,
                         install_rpath: str) -> T.Tuple[T.List[str], T.Set[bytes]]:
        if self.info.is_windows():
            return ([], set())

        # GNU ld, solaris ld, and lld acting like GNU ld
        if self.linker.id.startswith('ld'):
            # The way that dmd and ldc pass rpath to gcc is different than we would
            # do directly, each argument -rpath and the value to rpath, need to be
            # split into two separate arguments both prefaced with the -L=.
            args: T.List[str] = []
            (rpath_args, rpath_dirs_to_remove) = super().build_rpath_args(
                    env, build_dir, from_dir, rpath_paths, build_rpath, install_rpath)
            for r in rpath_args:
                if ',' in r:
                    a, b = r.split(',', maxsplit=1)
                    args.append(a)
                    args.append(self.LINKER_PREFIX + b)
                else:
                    args.append(r)
            return (args, rpath_dirs_to_remove)

        return super().build_rpath_args(
            env, build_dir, from_dir, rpath_paths, build_rpath, install_rpath)

    @classmethod
    def _translate_args_to_nongnu(cls, args: T.List[str], info: MachineInfo, link_id: str) -> T.List[str]:
        # Translate common arguments to flags the LDC/DMD compilers
        # can understand.
        # The flags might have been added by pkg-config files,
        # and are therefore out of the user's control.
        dcargs: T.List[str] = []
        # whether we hit a linker argument that expect another arg
        # see the comment in the "-L" section
        link_expect_arg = False
        link_flags_with_arg = [
            '-rpath', '-rpath-link', '-soname', '-compatibility_version', '-current_version',
        ]
        for arg in args:
            # Translate OS specific arguments first.
            osargs: T.List[str] = []
            if info.is_windows():
                osargs = cls.translate_arg_to_windows(arg)
            elif info.is_darwin():
                osargs = cls._translate_arg_to_osx(arg)
            if osargs:
                dcargs.extend(osargs)
                continue

            # Translate common D arguments here.
            if arg == '-pthread':
                continue
            if arg.startswith('-fstack-protector'):
                continue
            if arg.startswith('-D') and not (arg == '-D' or arg.startswith(('-Dd', '-Df'))):
                # ignore all '-D*' flags (like '-D_THREAD_SAFE')
                # unless they are related to documentation
                continue
            if arg.startswith('-Wl,'):
                # Translate linker arguments here.
                linkargs = arg[arg.index(',') + 1:].split(',')
                for la in linkargs:
                    dcargs.append('-L=' + la.strip())
                continue
            elif arg.startswith(('-link-defaultlib', '-linker', '-link-internally', '-linkonce-templates', '-lib')):
                # these are special arguments to the LDC linker call,
                # arguments like "-link-defaultlib-shared" do *not*
                # denote a library to be linked, but change the default
                # Phobos/DRuntime linking behavior, while "-linker" sets the
                # default linker.
                dcargs.append(arg)
                continue
            elif arg.startswith('-l'):
                # translate library link flag
                dcargs.append('-L=' + arg)
                continue
            elif arg.startswith('-isystem'):
                # translate -isystem system include path
                # this flag might sometimes be added by C library Cflags via
                # pkg-config.
                # NOTE: -isystem and -I are not 100% equivalent, so this is just
                # a workaround for the most common cases.
                if arg.startswith('-isystem='):
                    dcargs.append('-I=' + arg[9:])
                else:
                    dcargs.append('-I' + arg[8:])
                continue
            elif arg.startswith('-idirafter'):
                # same as -isystem, but appends the path instead
                if arg.startswith('-idirafter='):
                    dcargs.append('-I=' + arg[11:])
                else:
                    dcargs.append('-I' + arg[10:])
                continue
            elif arg.startswith('-L'):
                # The D linker expect library search paths in the form of -L=-L/path (the '=' is optional).
                #
                # This function receives a mix of arguments already prepended
                # with -L for the D linker driver and other linker arguments.
                # The arguments starting with -L can be:
                #  - library search path (with or without a second -L)
                #     - it can come from pkg-config (a single -L)
                #     - or from the user passing linker flags (-L-L would be expected)
                #  - arguments like "-L=-rpath" that expect a second argument (also prepended with -L)
                #  - arguments like "-L=@rpath/xxx" without a second argument (on Apple platform)
                #  - arguments like "-L=/SUBSYSTEM:CONSOLE (for Windows linker)
                #
                # The logic that follows tries to detect all these cases (some may be missing)
                # in order to prepend a -L only for the library search paths with a single -L

                if arg.startswith('-L='):
                    suffix = arg[3:]
                else:
                    suffix = arg[2:]

                if link_expect_arg:
                    # flags like rpath and soname expect a path or filename respectively,
                    # we must not alter it (i.e. prefixing with -L for a lib search path)
                    dcargs.append(arg)
                    link_expect_arg = False
                    continue

                if suffix in link_flags_with_arg:
                    link_expect_arg = True

                if suffix.startswith('-') or suffix.startswith('@'):
                    # this is not search path
                    dcargs.append(arg)
                    continue

                # linker flag such as -L=/DEBUG must pass through
                if info.is_windows() and link_id == 'link' and suffix.startswith('/'):
                    dcargs.append(arg)
                    continue

                # Make sure static library files are passed properly to the linker.
                if arg.endswith('.a') or arg.endswith('.lib'):
                    if len(suffix) > 0 and not suffix.startswith('-'):
                        dcargs.append('-L=' + suffix)
                        continue

                dcargs.append('-L=' + arg)
                continue
            elif not arg.startswith('-') and arg.endswith(('.a', '.lib')):
                # ensure static libraries are passed through to the linker
                dcargs.append('-L=' + arg)
                continue
            else:
                dcargs.append(arg)

        return dcargs

    @classmethod
    def translate_arg_to_windows(cls, arg: str) -> T.List[str]:
        args: T.List[str] = []
        if arg.startswith('-Wl,'):
            # Translate linker arguments here.
            linkargs = arg[arg.index(',') + 1:].split(',')
            for la in linkargs:
                if la.startswith('--out-implib='):
                    # Import library name
                    args.append('-L=/IMPLIB:' + la[13:].strip())
        elif arg.startswith('-mscrtlib='):
            args.append(arg)
            mscrtlib = arg[10:].lower()
            if cls is LLVMDCompiler:
                # Default crt libraries for LDC2 must be excluded for other
                # selected crt options.
                if mscrtlib != 'libcmt':
                    args.append('-L=/NODEFAULTLIB:libcmt')
                    args.append('-L=/NODEFAULTLIB:libvcruntime')

                # Fixes missing definitions for printf-functions in VS2017
                if mscrtlib.startswith('msvcrt'):
                    args.append('-L=/DEFAULTLIB:legacy_stdio_definitions.lib')

        return args

    @classmethod
    def _translate_arg_to_osx(cls, arg: str) -> T.List[str]:
        args: T.List[str] = []
        if arg.startswith('-install_name'):
            args.append('-L=' + arg)
        return args

    @classmethod
    def _unix_args_to_native(cls, args: T.List[str], info: MachineInfo, link_id: str = '') -> T.List[str]:
        return cls._translate_args_to_nongnu(args, info, link_id)

    def get_debug_args(self, is_debug: bool) -> T.List[str]:
        ddebug_args = []
        if is_debug:
            ddebug_args = [d_feature_args[self.id]['debug']]

        return clike_debug_args[is_debug] + ddebug_args

    def _get_crt_args(self, crt_val: str, buildtype: str) -> T.List[str]:
        if not self.info.is_windows():
            return []
        return self.mscrt_args[self.get_crt_val(crt_val, buildtype)]

    def get_soname_args(self, env: 'Environment', prefix: str, shlib_name: str,
                        suffix: str, soversion: str,
                        darwin_versions: T.Tuple[str, str]) -> T.List[str]:
        sargs = super().get_soname_args(env, prefix, shlib_name, suffix,
                                        soversion, darwin_versions)

        # LDC and DMD actually do use a linker, but they proxy all of that with
        # their own arguments
        soargs: T.List[str] = []
        if self.linker.id.startswith('ld.'):
            for arg in sargs:
                a, b = arg.split(',', maxsplit=1)
                soargs.append(a)
                soargs.append(self.LINKER_PREFIX + b)
            return soargs
        elif self.linker.id.startswith('ld64'):
            for arg in sargs:
                if not arg.startswith(self.LINKER_PREFIX):
                    soargs.append(self.LINKER_PREFIX + arg)
                else:
                    soargs.append(arg)
            return soargs
        else:
            return sargs

    def get_allow_undefined_link_args(self) -> T.List[str]:
        args = self.linker.get_allow_undefined_args()
        if self.info.is_darwin():
            # On macOS we're passing these options to the C compiler, but
            # they're linker options and need -Wl, so clang/gcc knows what to
            # do with them. I'm assuming, but don't know for certain, that
            # ldc/dmd do some kind of mapping internally for arguments they
            # understand, but pass arguments they don't understand directly.
            args = [a.replace('-L=', '-Xcc=-Wl,') for a in args]
        return args


class DCompilerArgs(CompilerArgs):
    prepend_prefixes = ('-I', '-L')
    dedup2_prefixes = ('-I', )


class DCompiler(Compiler):
    mscrt_args = {
        'none': ['-mscrtlib='],
        'md': ['-mscrtlib=msvcrt'],
        'mdd': ['-mscrtlib=msvcrtd'],
        'mt': ['-mscrtlib=libcmt'],
        'mtd': ['-mscrtlib=libcmtd'],
    }

    language = 'd'

    def __init__(self, exelist: T.List[str], version: str, for_machine: MachineChoice,
                 info: 'MachineInfo', arch: str, *,
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None,
                 is_cross: bool = False):
        super().__init__([], exelist, version, for_machine, info, linker=linker,
                         full_version=full_version, is_cross=is_cross)
        self.arch = arch

    def sanity_check(self, work_dir: str, environment: 'Environment') -> None:
        source_name = os.path.join(work_dir, 'sanity.d')
        output_name = os.path.join(work_dir, 'dtest')
        with open(source_name, 'w', encoding='utf-8') as ofile:
            ofile.write('''void main() { }''')

        compile_cmdlist = self.exelist + self.get_output_args(output_name) + self._get_target_arch_args() + [source_name]

        # If cross-compiling, we can't run the sanity check, only compile it.
        if self.is_cross and not environment.has_exe_wrapper():
            compile_cmdlist += self.get_compile_only_args()

        pc = subprocess.Popen(compile_cmdlist, cwd=work_dir)
        pc.wait()
        if pc.returncode != 0:
            raise EnvironmentException('D compiler %s cannot compile programs.' % self.name_string())

        if self.is_cross:
            if not environment.has_exe_wrapper():
                # Can't check if the binaries run so we have to assume they do
                return
            cmdlist = environment.exe_wrapper.get_command() + [output_name]
        else:
            cmdlist = [output_name]
        if subprocess.call(cmdlist) != 0:
            raise EnvironmentException('Executables created by D compiler %s are not runnable.' % self.name_string())

    def needs_static_linker(self) -> bool:
        return True

    def get_depfile_suffix(self) -> str:
        return 'deps'

    def get_pic_args(self) -> T.List[str]:
        if self.info.is_windows():
            return []
        return ['-fPIC']

    def get_feature_args(self, kwargs: DFeatures, build_to_src: str) -> T.List[str]:
        res: T.List[str] = []
        unittest_arg = d_feature_args[self.id]['unittest']
        if not unittest_arg:
            raise EnvironmentException('D compiler %s does not support the "unittest" feature.' % self.name_string())
        if kwargs['unittest']:
            res.append(unittest_arg)

        debug_level = -1
        debug_arg = d_feature_args[self.id]['debug']
        if not debug_arg:
            raise EnvironmentException('D compiler %s does not support conditional debug identifiers.' % self.name_string())

        # Parse all debug identifiers and the largest debug level identifier
        for d in kwargs['debug']:
            if isinstance(d, int):
                debug_level = max(debug_level, d)
            elif isinstance(d, str) and d.isdigit():
                debug_level = max(debug_level, int(d))
            else:
                res.append(f'{debug_arg}={d}')

        if debug_level >= 0:
            res.append(f'{debug_arg}={debug_level}')

        version_level = -1
        version_arg = d_feature_args[self.id]['version']
        if not version_arg:
            raise EnvironmentException('D compiler %s does not support conditional version identifiers.' % self.name_string())

        # Parse all version identifiers and the largest version level identifier
        for v in kwargs['versions']:
            if isinstance(v, int):
                version_level = max(version_level, v)
            elif isinstance(v, str) and v.isdigit():
                version_level = max(version_level, int(v))
            else:
                res.append(f'{version_arg}={v}')

        if version_level >= 0:
            res.append(f'{version_arg}={version_level}')

        import_dir_arg = d_feature_args[self.id]['import_dir']
        if not import_dir_arg:
            raise EnvironmentException('D compiler %s does not support the "string import directories" feature.' % self.name_string())
        # TODO: ImportDirs.to_string_list(), but we need both the project source
        # root and project build root for that.
        for idir_obj in kwargs['import_dirs']:
            basedir = idir_obj.get_curdir()
            for idir in idir_obj.get_incdirs():
                bldtreedir = os.path.join(basedir, idir)
                # Avoid superfluous '/.' at the end of paths when d is '.'
                if idir not in ('', '.'):
                    expdir = bldtreedir
                else:
                    expdir = basedir
                srctreedir = os.path.join(build_to_src, expdir)
                res.append(f'{import_dir_arg}{srctreedir}')
                res.append(f'{import_dir_arg}{bldtreedir}')

        return res

    def get_optimization_link_args(self, optimization_level: str) -> T.List[str]:
        if optimization_level != 'plain':
            return self._get_target_arch_args()
        return []

    def compiler_args(self, args: T.Optional[T.Iterable[str]] = None) -> DCompilerArgs:
        return DCompilerArgs(self, args)

    def has_multi_arguments(self, args: T.List[str], env: 'Environment') -> T.Tuple[bool, bool]:
        return self.compiles('int i;\n', env, extra_args=args)

    def _get_target_arch_args(self) -> T.List[str]:
        # LDC2 on Windows targets to current OS architecture, but
        # it should follow the target specified by the MSVC toolchain.
        if self.info.is_windows():
            if self.is_cross:
                return [f'-mtriple={self.arch}-windows-msvc']
            elif self.arch == 'x86_64':
                return ['-m64']
            return ['-m32']
        return []

    def get_crt_compile_args(self, crt_val: str, buildtype: str) -> T.List[str]:
        return []

    def get_crt_link_args(self, crt_val: str, buildtype: str) -> T.List[str]:
        return []

    def _get_compile_extra_args(self, extra_args: T.Union[T.List[str], T.Callable[[CompileCheckMode], T.List[str]], None] = None) -> T.List[str]:
        args = self._get_target_arch_args()
        if extra_args:
            if callable(extra_args):
                extra_args = extra_args(CompileCheckMode.COMPILE)
            if isinstance(extra_args, list):
                args.extend(extra_args)
            elif isinstance(extra_args, str):
                args.append(extra_args)
        return args

    def run(self, code: 'mesonlib.FileOrString', env: 'Environment',
            extra_args: T.Union[T.List[str], T.Callable[[CompileCheckMode], T.List[str]], None] = None,
            dependencies: T.Optional[T.List['Dependency']] = None,
            run_env: T.Optional[T.Dict[str, str]] = None,
            run_cwd: T.Optional[str] = None) -> compilers.RunResult:
        extra_args = self._get_compile_extra_args(extra_args)
        return super().run(code, env, extra_args, dependencies, run_env, run_cwd)

    def sizeof(self, typename: str, prefix: str, env: 'Environment', *,
               extra_args: T.Union[None, T.List[str], T.Callable[[CompileCheckMode], T.List[str]]] = None,
               dependencies: T.Optional[T.List['Dependency']] = None) -> T.Tuple[int, bool]:
        if extra_args is None:
            extra_args = []
        t = f'''
        import std.stdio : writeln;
        {prefix}
        void main() {{
            writeln(({typename}).sizeof);
        }}
        '''
        res = self.cached_run(t, env, extra_args=extra_args,
                              dependencies=dependencies)
        if not res.compiled:
            return -1, False
        if res.returncode != 0:
            raise mesonlib.EnvironmentException('Could not run sizeof test binary.')
        return int(res.stdout), res.cached

    def alignment(self, typename: str, prefix: str, env: 'Environment', *,
                  extra_args: T.Optional[T.List[str]] = None,
                  dependencies: T.Optional[T.List['Dependency']] = None) -> T.Tuple[int, bool]:
        if extra_args is None:
            extra_args = []
        t = f'''
        import std.stdio : writeln;
        {prefix}
        void main() {{
            writeln(({typename}).alignof);
        }}
        '''
        res = self.run(t, env, extra_args=extra_args,
                       dependencies=dependencies)
        if not res.compiled:
            raise mesonlib.EnvironmentException('Could not compile alignment test.')
        if res.returncode != 0:
            raise mesonlib.EnvironmentException('Could not run alignment test binary.')
        align = int(res.stdout)
        if align == 0:
            raise mesonlib.EnvironmentException(f'Could not determine alignment of {typename}. Sorry. You might want to file a bug.')
        return align, res.cached

    def has_header(self, hname: str, prefix: str, env: 'Environment', *,
                   extra_args: T.Union[None, T.List[str], T.Callable[['CompileCheckMode'], T.List[str]]] = None,
                   dependencies: T.Optional[T.List['Dependency']] = None,
                   disable_cache: bool = False) -> T.Tuple[bool, bool]:

        extra_args = self._get_compile_extra_args(extra_args)
        code = f'''{prefix}
        import {hname};
        '''
        return self.compiles(code, env, extra_args=extra_args,
                             dependencies=dependencies, mode=CompileCheckMode.COMPILE, disable_cache=disable_cache)

class GnuDCompiler(GnuCompiler, DCompiler):

    # we mostly want DCompiler, but that gives us the Compiler.LINKER_PREFIX instead
    LINKER_PREFIX = GnuCompiler.LINKER_PREFIX
    id = 'gcc'

    def __init__(self, exelist: T.List[str], version: str, for_machine: MachineChoice,
                 info: 'MachineInfo', arch: str, *,
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None,
                 is_cross: bool = False):
        DCompiler.__init__(self, exelist, version, for_machine, info, arch,
                           linker=linker,
                           full_version=full_version, is_cross=is_cross)
        GnuCompiler.__init__(self, {})
        default_warn_args = ['-Wall', '-Wdeprecated']
        self.warn_args = {'0': [],
                          '1': default_warn_args,
                          '2': default_warn_args + ['-Wextra'],
                          '3': default_warn_args + ['-Wextra', '-Wpedantic'],
                          'everything': (default_warn_args + ['-Wextra', '-Wpedantic'] +
                                         self.supported_warn_args(gnu_common_warning_args))}

        self.base_options = {
            OptionKey(o) for o in [
             'b_colorout', 'b_sanitize', 'b_staticpic', 'b_vscrt',
             'b_coverage', 'b_pgo', 'b_ndebug']}

        self._has_color_support = version_compare(self.version, '>=4.9')
        # dependencies were implemented before, but broken - support was fixed in GCC 7.1+
        # (and some backported versions)
        self._has_deps_support = version_compare(self.version, '>=7.1')

    def get_colorout_args(self, colortype: str) -> T.List[str]:
        if self._has_color_support:
            super().get_colorout_args(colortype)
        return []

    def get_dependency_gen_args(self, outtarget: str, outfile: str) -> T.List[str]:
        if self._has_deps_support:
            return super().get_dependency_gen_args(outtarget, outfile)
        return []

    def get_warn_args(self, level: str) -> T.List[str]:
        return self.warn_args[level]

    def get_optimization_args(self, optimization_level: str) -> T.List[str]:
        return gdc_optimization_args[optimization_level]

    def compute_parameters_with_absolute_paths(self, parameter_list: T.List[str],
                                               build_dir: str) -> T.List[str]:
        for idx, i in enumerate(parameter_list):
            if i[:2] == '-I' or i[:2] == '-L':
                parameter_list[idx] = i[:2] + os.path.normpath(os.path.join(build_dir, i[2:]))

        return parameter_list

    def get_allow_undefined_link_args(self) -> T.List[str]:
        return self.linker.get_allow_undefined_args()

    def get_linker_always_args(self) -> T.List[str]:
        args = super().get_linker_always_args()
        if self.info.is_windows():
            return args
        return args + ['-shared-libphobos']

    def get_assert_args(self, disable: bool, env: 'Environment') -> T.List[str]:
        if disable:
            return ['-frelease']
        return []

# LDC uses the DMD frontend code to parse and analyse the code.
# It then uses LLVM for the binary code generation and optimizations.
# This function retrieves the dmd frontend version, which determines
# the common features between LDC and DMD.
# We need the complete version text because the match is not on first line
# of version_output
def find_ldc_dmd_frontend_version(version_output: T.Optional[str]) -> T.Optional[str]:
    if version_output is None:
        return None
    version_regex = re.search(r'DMD v(\d+\.\d+\.\d+)', version_output)
    if version_regex:
        return version_regex.group(1)
    return None

class LLVMDCompiler(DmdLikeCompilerMixin, DCompiler):

    id = 'llvm'

    def __init__(self, exelist: T.List[str], version: str, for_machine: MachineChoice,
                 info: 'MachineInfo', arch: str, *,
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None,
                 is_cross: bool = False, version_output: T.Optional[str] = None):
        DCompiler.__init__(self, exelist, version, for_machine, info, arch,
                           linker=linker,
                           full_version=full_version, is_cross=is_cross)
        DmdLikeCompilerMixin.__init__(self, dmd_frontend_version=find_ldc_dmd_frontend_version(version_output))
        self.base_options = {OptionKey(o) for o in ['b_coverage', 'b_colorout', 'b_vscrt', 'b_ndebug']}

    def get_colorout_args(self, colortype: str) -> T.List[str]:
        if colortype == 'always':
            return ['-enable-color']
        return []

    def get_warn_args(self, level: str) -> T.List[str]:
        if level in {'2', '3'}:
            return ['-wi', '-dw']
        elif level == '1':
            return ['-wi']
        return []

    def get_pic_args(self) -> T.List[str]:
        return ['-relocation-model=pic']

    def get_crt_link_args(self, crt_val: str, buildtype: str) -> T.List[str]:
        return self._get_crt_args(crt_val, buildtype)

    def unix_args_to_native(self, args: T.List[str]) -> T.List[str]:
        return self._unix_args_to_native(args, self.info, self.linker.id)

    def get_optimization_args(self, optimization_level: str) -> T.List[str]:
        if optimization_level != 'plain':
            return self._get_target_arch_args() + ldc_optimization_args[optimization_level]
        return ldc_optimization_args[optimization_level]

    @classmethod
    def use_linker_args(cls, linker: str, version: str) -> T.List[str]:
        return [f'-linker={linker}']

    def get_linker_always_args(self) -> T.List[str]:
        args = super().get_linker_always_args()
        if self.info.is_windows():
            return args
        return args + ['-link-defaultlib-shared']

    def get_assert_args(self, disable: bool, env: 'Environment') -> T.List[str]:
        if disable:
            return ['--release']
        return []

    def rsp_file_syntax(self) -> RSPFileSyntax:
        # We use `mesonlib.is_windows` here because we want to know what the
        # build machine is, not the host machine. This really means we would
        # have the Environment not the MachineInfo in the compiler.
        return RSPFileSyntax.MSVC if is_windows() else RSPFileSyntax.GCC


class DmdDCompiler(DmdLikeCompilerMixin, DCompiler):

    id = 'dmd'

    def __init__(self, exelist: T.List[str], version: str, for_machine: MachineChoice,
                 info: 'MachineInfo', arch: str, *,
                 linker: T.Optional['DynamicLinker'] = None,
                 full_version: T.Optional[str] = None,
                 is_cross: bool = False):
        DCompiler.__init__(self, exelist, version, for_machine, info, arch,
                           linker=linker,
                           full_version=full_version, is_cross=is_cross)
        DmdLikeCompilerMixin.__init__(self, version)
        self.base_options = {OptionKey(o) for o in ['b_coverage', 'b_colorout', 'b_vscrt', 'b_ndebug']}

    def get_colorout_args(self, colortype: str) -> T.List[str]:
        if colortype == 'always':
            return ['-color=on']
        return []

    def get_std_exe_link_args(self) -> T.List[str]:
        if self.info.is_windows():
            # DMD links against D runtime only when main symbol is found,
            # so these needs to be inserted when linking static D libraries.
            if self.arch == 'x86_64':
                return ['phobos64.lib']
            elif self.arch == 'x86_mscoff':
                return ['phobos32mscoff.lib']
            return ['phobos.lib']
        return []

    def get_std_shared_lib_link_args(self) -> T.List[str]:
        libname = 'libphobos2.so'
        if self.info.is_windows():
            if self.arch == 'x86_64':
                libname = 'phobos64.lib'
            elif self.arch == 'x86_mscoff':
                libname = 'phobos32mscoff.lib'
            else:
                libname = 'phobos.lib'
        return ['-shared', '-defaultlib=' + libname]

    def _get_target_arch_args(self) -> T.List[str]:
        # DMD32 and DMD64 on 64-bit Windows defaults to 32-bit (OMF).
        # Force the target to 64-bit in order to stay consistent
        # across the different platforms.
        if self.info.is_windows():
            if self.arch == 'x86_64':
                return ['-m64']
            elif self.arch == 'x86_mscoff':
                return ['-m32mscoff']
            return ['-m32']
        return []

    def get_crt_compile_args(self, crt_val: str, buildtype: str) -> T.List[str]:
        return self._get_crt_args(crt_val, buildtype)

    def unix_args_to_native(self, args: T.List[str]) -> T.List[str]:
        return self._unix_args_to_native(args, self.info, self.linker.id)

    def get_optimization_args(self, optimization_level: str) -> T.List[str]:
        if optimization_level != 'plain':
            return self._get_target_arch_args() + dmd_optimization_args[optimization_level]
        return dmd_optimization_args[optimization_level]

    def can_linker_accept_rsp(self) -> bool:
        return False

    def get_linker_always_args(self) -> T.List[str]:
        args = super().get_linker_always_args()
        if self.info.is_windows():
            return args
        return args + ['-defaultlib=phobos2', '-debuglib=phobos2']

    def get_assert_args(self, disable: bool, env: 'Environment') -> T.List[str]:
        if disable:
            return ['-release']
        return []

    def rsp_file_syntax(self) -> RSPFileSyntax:
        return RSPFileSyntax.MSVC
